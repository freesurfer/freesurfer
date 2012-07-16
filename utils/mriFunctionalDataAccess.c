/**
 * @file  mriFunctionalDataAccess.c
 * @brief fmri utils
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/07/16 01:29:22 $
 *    $Revision: 1.57 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "machine.h"
#include "matrix.h"
#include "stats.h"
#include "mriFunctionalDataAccess.h"
#include "xDebug.h"
#include "xVoxel.h"
#include "xUtilities.h"
#include "macros.h"
#include "error.h"
#include "mri2.h"
#include "mri_identify.h"

#ifndef min
#define min(x,y) (x<y?x:y)
#endif

#ifndef max
#define max(x,y) (x>y?x:y)
#endif

#define KEEP_NULL_CONDITION 1

#define FunD_ksStemHeaderSuffix  "dat"
#define FunD_ksRegisterFileName  "register.dat"


char *FunD_ksaErrorString [FunD_tErr_knNumErrorCodes] =
{
  "No error.",
  "Invalid object signature.",
  "Invalid ptr to volume (was probably NULL).",
  "Path not found.",
  "No stem provided and it couldn't be guessed.",
  "Data not found in specified path.",
  "Header file not found in specified path.",
  "Couldn't allocate volume (memory allocation failed).",
  "Unrecognized header format.",
  "Questionable header format (found expected types of values, "
  "but keywords were different).",
  "Couldn't find a recognizable data file.",
  "Couldn't allocate storage (memory allocation failed).",
  "Couldn't allocate MRI.",
  "Couldn't allocate matrix (not my fault).",
  "Couldn't load MRI volume.",
  "Couldn't read register file.",
  "Couldn't calculate deviations.",
  "Error performing an operation on a transform.",
  "Invalid parameter",
  "Invalid time resolution, number of time points must be "
  "evenly divisble by it.",
  "Invalid number of pre stim time points.",
  "Invalid functional voxel, out of bounds.",
  "Invalid condition (is it zero-based?)",
  "Invalid time point (is it zero-based?)",
  "Couldn't transform mask into local space.",
  "Error performing FDR.",
  "Invalid error code."
};

// ===================================================================== VOLUME

FunD_tErr FunD_New ( mriFunctionalDataRef*  opVolume,
                     char*                  isFileName,
                     FunD_tRegistrationType iRegistrationType,
                     char*                  isRegistrationFileName,
                     int                    inScalarSize,
                     mriVolumeRef           iAnatomicalVolume,
                     tBoolean               ibIsLeftHemisphere )
{
  mriFunctionalDataRef this    = NULL;
  FunD_tErr            eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_New( opVolume=%p, isFileName=%s, "
                       "iRegistrationType=%d, isRegistrationFileName=%s, "
                       "inScalarSize=%d, iAnatomicalVolume=%p)",
                       opVolume, isFileName, iRegistrationType,
                       isRegistrationFileName, inScalarSize,
                       iAnatomicalVolume) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != opVolume &&
                      NULL != isFileName ),
                     eResult, FunD_tErr_InvalidParameter );

  /* Allocate the storage */
  DebugNote( ("Allocating mriFunctionalData") );
  this = (mriFunctionalDataRef) malloc(sizeof(mriFunctionalData));
  DebugAssertThrowX( (NULL !=this), eResult, FunD_tErr_CouldntAllocateVolume );

  this->mSignature = FunD_kSignature;

  /* Save registration name. If we didn't get one, build one out of
     the base path of the file name. */
  DebugNote( ("Saving file name") );
  strcpy( this->msFileName, isFileName );
  if ( NULL != isRegistrationFileName )
  {
    DebugNote( ("Saving registration file name") );
    strcpy( this->msRegistrationFileName, isRegistrationFileName );
  }
  else
  {
    strcpy( this->msRegistrationFileName, "" );
  }

  /* Init values. */
  this->mSampleType = FunD_tSampleType_Nearest;
  this->mConvMethod = FunD_tConversionMethod_FFF;
  this->mNumTimePoints = -1;
  this->mNumConditions = -1;
  this->mbNullConditionPresent = FALSE;
  this->mpData = NULL;
  this->mMinValue = 10000;
  this->mMaxValue = -10000;
  this->mTimeResolution = 0;
  this->mNumPreStimTimePoints = 0;
  this->mIdxToIdxTransform = NULL;
  this->mOriginalIdxToIdxTransform = NULL;
  this->mbErrorDataPresent = FALSE;
  this->mCovMtx = NULL;
  this->mpResampledData = NULL;
  this->mClientXMin = 0;
  this->mClientYMin = 0;
  this->mClientZMin = 0;
  this->mClientXMax = -1;
  this->mClientYMax = -1;
  this->mClientZMax = -1;
  this->mbHaveClientBounds = FALSE;
  this->mbScalar = FALSE;

  this->mFrequencies = NULL;
  this->mNumBins = 0;

  /* Init transform objects */
  DebugNote( ("Creating idx to idx transform") );
  Trns_New( &(this->mIdxToIdxTransform) );

  /* Load the data.*/
  this->mpData = MRIread( this->msFileName );
  DebugAssertThrowX( (NULL != this->mpData),
                     eResult, FunD_tErr_CouldntReadMRI );

  /* Try to parse a stem header if we have a stem specified. */
  DebugNote( ("Trying stem header") );
  eResult = FunD_FindAndParseStemHeader_( this );
  if ( FunD_tErr_NoError != eResult )
  {

    DebugNote( ("Guess meta information") );
    eResult = FunD_GuessMetaInformation_( this );
  }

  /* Try to reshape if we can. */
  if ( inScalarSize > 0 )
  {
    DebugNote( ("Trying reshape with %d values", inScalarSize) );
    FunD_ReshapeIfScalar_( this, inScalarSize, NULL, ibIsLeftHemisphere );
  }

  /* If we're not scalar by now, parse the registration file */
  if ( !this->mbScalar )
  {
    DebugNote( ("Parsing registration file") );
    eResult = FunD_ParseRegistrationAndInitMatricies_( this,
              iRegistrationType,
              iAnatomicalVolume );
    DebugAssertThrow( (FunD_tErr_NoError == eResult) );
  }

  FunD_DebugPrint( this );

  /* Get the value range. */
  DebugNote( ("Getting value range") );
  MRIvalRange( this->mpData, &this->mMinValue, &this->mMaxValue );

  /* return what we got */
  *opVolume = this;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );

  /* delete data if no good */
  if ( NULL != this )
  {
    free( this );
    this = NULL;
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_Delete ( mriFunctionalDataRef* iopVolume )
{

  FunD_tErr            eResult = FunD_tErr_NoError;
  mriFunctionalDataRef this    = NULL;
  int                  column  = 0;

  DebugEnterFunction( ("FunD_Delete( iopVolume=%p )", iopVolume) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != iopVolume),eResult, FunD_tErr_InvalidParameter );

  /* Get the object. */
  this = *iopVolume;

  /* Make sure we're valid. */
  eResult = FunD_Verify ( this );
  DebugAssertThrow( (FunD_tErr_NoError == eResult) );

  /* mangle our signature */
  this->mSignature = 0x1;

  /* Delete other data. */
  if ( NULL != this->mpData )
  {
    MRIfree( &this->mpData );
  }
  if ( NULL != this->mCovMtx )
  {
    for ( column = 0;
          column < (this->mNumConditions-1) * this->mNumTimePoints;
          column++ )
    {
      free( this->mCovMtx[column] );
    }
    free( this->mCovMtx );
  }

  /* delete our transformers. */
  if ( NULL != this->mIdxToIdxTransform )
  {
    Trns_Delete( &(this->mIdxToIdxTransform) );
  }
  if ( NULL != this->mOriginalIdxToIdxTransform )
  {
    Trns_Delete( &(this->mOriginalIdxToIdxTransform) );
  }

  /* delete the structure. */
  free( this );

  *iopVolume = NULL;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_ReshapeIfScalar_ ( mriFunctionalDataRef this,
                                  int                  inNumValues,
                                  tBoolean*            obReshaped ,
                                  tBoolean             ibIsLeftHemisphere)
{

  FunD_tErr eResult        = FunD_tErr_NoError;
  MRI*      reshapedVolume = NULL;

  DebugEnterFunction( ("FunD_ReshapeIfScalar_( this=%p, inNumValues=%d, "
                       "obReshaped=%p )", this, inNumValues, obReshaped) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (inNumValues > 0), eResult, FunD_tErr_InvalidParameter );

  if (this->mpData->width*this->mpData->height*this->mpData->depth == 
      2*inNumValues)
  {
    MRI *mri_tmp, *mri_tmp2 ;
    this->mbScalar = 1 ;
    if (ibIsLeftHemisphere)
    {
      mri_tmp = MRIextractInto(this->mpData, NULL, 0, 0, 0, 
			       this->mpData->width/2, this->mpData->height, this->mpData->depth,
			       0,0,0) ;
      mri_tmp2 = MRIcopyFrames(mri_tmp, NULL, 0, inNumValues-1, 0) ;
    }
    else
    {
      mri_tmp = MRIextractInto(this->mpData, NULL, 0, 0, 0, 
			       this->mpData->width/2, this->mpData->height, this->mpData->depth,
			       inNumValues,0,0) ;
      mri_tmp2 = MRIcopyFrames(mri_tmp, NULL, inNumValues, 2*inNumValues-1, 0) ;
    }
    MRIfree(&mri_tmp) ;
    MRIfree(&this->mpData) ; this->mpData = mri_tmp2 ;
  }

  if ( this->mpData->width * this->mpData->height * this->mpData->depth ==
       inNumValues )
  {

    /* Reshape it so that everything is in the x direction. */
    DebugNote( ("Calling mri_reshape") );
    reshapedVolume =
      mri_reshape( this->mpData, inNumValues, 1, 1, this->mpData->nframes );
    DebugAssertThrowX( (NULL != reshapedVolume), eResult,
                       FunD_tErr_CouldntAllocateVolume );

    DebugNote( ("Freeing current data") );
    MRIfree( &this->mpData );
    this->mpData = reshapedVolume;
    reshapedVolume = NULL;

    /* If we're scalar, we don't use a registration or transform. */
    this->mbScalar = TRUE;

  }

  if ( NULL != obReshaped )
  {
    DebugNote( ("Returning whether it was reshaped") );
    *obReshaped = this->mbScalar = TRUE;
  }

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );

  if ( NULL != reshapedVolume )
  {
    DebugNote( ("Freeing reshaped volume") );
    MRIfree( &reshapedVolume );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_FindAndParseStemHeader_ ( mriFunctionalDataRef this )
{

  FunD_tErr eResult        = FunD_tErr_NoError;
  tBoolean  bFoundFile = FALSE;
  char      sFileName[FunD_knPathLen] = "",*ptmp;
  float     fPreStimSecs   = 0;
  FILE*     pHeader        = NULL;
  tBoolean  bGood          = FALSE;
  tBoolean  bSomethingRead = FALSE;
  char      sKeyword[256]  = "";
  int       nNumValues     = 0;
  int       nCovMtxCol     = 0;
  int       nCovMtxRow     = 0;
  int       nValue         = 0;
  int       nValuesRead    = 0;


  if(getenv("USEDNGDRAW"))
  {
    printf("FunD_FindAndParseStemHeader_ msFileName %s\n",this->msFileName);
    printf("Using dng draw %s\n",getenv("USEDNGDRAW"));
    fflush(stdout);
    this->mNumConditions = CountItemsInString(getenv("USEDNGDRAW")) + 1;
    printf("Using dng draw nc=%d\n",this->mNumConditions-1);
    fflush(stdout);
    bFoundFile = FALSE;
    eResult = FunD_tErr_NoError;
    //DebugAssertQuietThrowX( (bFoundFile), eResult, FunD_tErr_HeaderNotFound );
    //Cant have overlay and plot at the same time!
    return(FunD_tErr_NoError);
  }

  DebugEnterFunction( ("FunD_FindAndParseStemHeader_( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  // Get the name of any .dat file
  ptmp = IDstemFromName(this->msFileName);
  sprintf(sFileName,"%s.dat",ptmp);
  free(ptmp);
  /* Try to open it. */
  printf("tkmedit: Trying to open %s\n", sFileName);
  DebugNote( ("Trying to open %s", sFileName) );
  pHeader = fopen( sFileName, "r" );
  if ( NULL != pHeader )
  {
    fclose( pHeader );
    bFoundFile = TRUE;
  }
  /* Bail if no file. */
  DebugAssertQuietThrowX( (bFoundFile), eResult, FunD_tErr_HeaderNotFound );

  /* Try to open the file. */
  pHeader = fopen( sFileName, "r" );

  /* Start scanning for keywords... We use bSomethingRead to keep
     track if recognized a keyword and tried to read something. Since
     we are happy to skip unrecognized lines, we could have a totally
     random .dat file there that we would skip over and not get any
     valid meta information. So if end up going through the file and
     not recognizing anything, we want to bail. This will be the case
     is bSomethingRead is still false at the end of the reading. */
  bSomethingRead = FALSE;
  while (!feof( pHeader ) )
  {

    /* grab a keyword */
    DebugNote( ("Reading a keyword") );
    nValuesRead = fscanf( pHeader, "%s", sKeyword );
    if ( feof( pHeader ) )
    {
      break;
    }
    DebugAssertThrowX( (1 == nValuesRead && !feof( pHeader )),
                       eResult, FunD_tErr_UnrecognizedHeaderFormat );

    /* assume success */
    bGood = TRUE;

    /* look at the keyword */
    if ( strcmp( sKeyword, "TER" ) == 0 )
    {
      DebugNote( ("Reading TER") );
      nValuesRead = fscanf( pHeader, "%f", &this->mTimeResolution );
      bGood = (1 == nValuesRead);
      bSomethingRead = TRUE;
    }
    else if ( strcmp( sKeyword, "TPreStim" ) == 0 )
    {
      DebugNote( ("Reading TPreStim") );
      nValuesRead = fscanf( pHeader, "%f", &fPreStimSecs );
      bGood = (1 == nValuesRead);
      bSomethingRead = TRUE;
    }
    else if ( strcmp( sKeyword, "nCond" ) == 0 )
    {
      DebugNote( ("Reading nCond") );
      nValuesRead = fscanf( pHeader, "%d", &this->mNumConditions );
      bGood = (1 == nValuesRead);
      this->mbNullConditionPresent = TRUE;
      bSomethingRead = TRUE;
    }
    else if ( strcmp( sKeyword, "Nh" ) == 0 )
    {
      DebugNote( ("Reading Nh") );
      nValuesRead = fscanf( pHeader, "%d", &this->mNumTimePoints );
      bGood = (1 == nValuesRead);
      bSomethingRead = TRUE;
    }
    else if ( strcmp( sKeyword, "SumXtX" ) == 0 )
    {
      DebugNote( ("Reading SumXtX") );
      nNumValues = pow( this->mNumTimePoints * (this->mNumConditions-1), 2 );
      for ( nValue = 0; nValue < nNumValues; nValue++ )
      {
        fscanf( pHeader, "%*d" );
      }
      bSomethingRead = TRUE;
    }
    else if ( strcmp( sKeyword, "hCovMtx" ) == 0 )
    {
      /* Allocate the covariance matrix. It's the size of the number
         of time points times conditions (minus null condition) on both
         sides. */
      DebugNote( ("Allocating cov mtx") );
      this->mCovMtx =
        (float**) calloc (this->mNumTimePoints * (this->mNumConditions-1),
                          sizeof(float*) );
      DebugAssertThrowX( (NULL != this->mCovMtx),
                         eResult, FunD_tErr_CouldntAllocateStorage );

      /* The cov mtx doesn't have values for the null condition, so we
         start with 1 instead of 0 in our loops here, and use that-1. */
      for ( nCovMtxRow = 0;
            nCovMtxRow < this->mNumTimePoints * (this->mNumConditions-1);
            nCovMtxRow++ )
      {

        DebugNote( ("Allocating cov mtx col %d", nCovMtxCol) );
        this->mCovMtx[nCovMtxRow] =
          (float*) calloc (this->mNumTimePoints * (this->mNumConditions-1),
                           sizeof(float*) );
        DebugAssertThrowX( (NULL != this->mCovMtx),
                           eResult, FunD_tErr_CouldntAllocateStorage );

        for ( nCovMtxCol = 0;
              nCovMtxCol < this->mNumTimePoints * (this->mNumConditions-1);
              nCovMtxCol++ )
        {


          DebugNote( ("Reading cov mtx value %d, %d",
                      nCovMtxCol, nCovMtxRow) );
          nValuesRead = fscanf( pHeader, "%f",
                                &this->mCovMtx[nCovMtxRow][nCovMtxCol] );
          bGood = (1 == nValuesRead);
          DebugAssertThrowX( (bGood), eResult,
                             FunD_tErr_UnrecognizedHeaderFormat );
        }
      }
      bSomethingRead = TRUE;
    }

    DebugAssertThrowX( (bGood), eResult, FunD_tErr_UnrecognizedHeaderFormat );
  }

  /* If nothing was read, this isn't a real .dat file, so bail. */
  DebugAssertThrowX( (bSomethingRead),
                     eResult, FunD_tErr_UnrecognizedHeaderFormat );

  /* Divide the num prestim points by the time res. Do it here because
     we might not have gotten the timeres when we read TPreStim up
     there. */
  this->mNumPreStimTimePoints = fPreStimSecs / this->mTimeResolution;

  /* In this format we have error data. */
  this->mbErrorDataPresent = TRUE;


  DebugCatch;

  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  DebugPrint( ("Last parsed keyword: '%s' nValuesRead=%d\n",
               sKeyword, nValuesRead ) );
  EndDebugCatch;

  if ( NULL != pHeader )
  {
    fclose( pHeader );
  }

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GuessMetaInformation_ ( mriFunctionalDataRef this )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GuessMetaInformation_( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  /* Set default values from MRI info. Just one condition, with all
     frames being time points, and no error data.  */
  this->mTimeResolution = 1.0;
  this->mNumConditions = 1;
  this->mNumTimePoints = this->mpData->nframes;
  this->mNumPreStimTimePoints = 0;
  this->mbErrorDataPresent = FALSE;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr
FunD_ParseRegistrationAndInitMatricies_ ( mriFunctionalDataRef   this,
    FunD_tRegistrationType iType,
    mriVolumeRef    iAnatomicalVolume )
{
  FunD_tErr    eResult          = FunD_tErr_NoError;
  char         sBasePath[1024]  = "";
  char         sFileName[1024]  = "";
  struct stat  fileInfo;
  int          eStat            = 0;
  char*        pCurChar         = NULL;
  char*        pBaseEnd         = NULL;
  fMRI_REG*    theRegInfo       = NULL;
  FILE*        fRegister        = NULL;
  char         sLine[1024]      = "";
  char         sKeyWord[1024]   = "";
  tBoolean     bGood            = FALSE;
  FunD_tConversionMethod
  convMethod       = FunD_tConversionMethod_FCF;
  int          eMRI             = ERROR_NONE;

  DebugEnterFunction( ("FunD_ParseRegistrationAndInitMatricies_( this=%p )",
                       this) );

  DebugNote( ("Checking this") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking registration type") );
  DebugAssertThrowX( (iType >= 0 &&
                      iType < FunD_knNumRegistrationTypes), eResult,
                     FunD_tErr_InvalidParameter );
  DebugNote( ("Checking iAnatomicalVolume") );
  DebugAssertThrowX( (NULL != iAnatomicalVolume),
                     eResult, FunD_tErr_InvalidParameter );

  /* Switch on the registration type. We'll need to get some info out
     of the registration file if they specified one, or generate some
     default info. */
  switch ( iType )
  {

    /* These both try to load files, so we'll share some code. */
  case FunD_tRegistration_File:
  case FunD_tRegistration_Find:

    /* If we need to find a file, look for it now. */
    if ( FunD_tRegistration_Find == iType )
    {
      /* Construct a registration name based on our data location and
         try to open that. If we didn't find one, return an error.
         Start by finding the last slash in the data file name. */
      DebugNote( ("Looking for last slash") );
      xUtil_strncpy( sBasePath, this->msFileName, sizeof(sBasePath));
      pCurChar = sBasePath;
      while ( NULL != pCurChar && *pCurChar != '\0' )
      {
        if ( '/' == *pCurChar )
        {
          pBaseEnd = pCurChar;
        }
        pCurChar++;
      }

      /* Set this to null char, terminating the string. Then use this
         string and the header suffix to build the header file name. */
      *pBaseEnd = '\0';
      xUtil_snprintf( sFileName, sizeof(sFileName),
                      "%s/%s", sBasePath, FunD_ksRegisterFileName );

      /* Save the file name we found. */
      strncpy( this->msRegistrationFileName, sFileName,
               sizeof(this->msRegistrationFileName) );

    }

    /* See if the file exists and is a regular file. */
    DebugNote( ("Checking if %s exists", this->msRegistrationFileName) );
    eStat = stat( this->msRegistrationFileName, &fileInfo );
    DebugAssertThrowX( (0 == eStat),
                       eResult, FunD_tErr_CouldntReadRegisterFile );

    DebugNote( ("Checking S_ISREG") );
    DebugAssertThrowX( ( S_ISREG(fileInfo.st_mode) ),
                       eResult, FunD_tErr_CouldntReadRegisterFile );

    /* Read the registration info. */
    DebugNote( ("StatReadRegistration on %s\n",
                this->msRegistrationFileName ) );
    theRegInfo = StatReadRegistration( this->msRegistrationFileName );
    DebugAssertThrowX( (NULL != theRegInfo ),
                       eResult, FunD_tErr_CouldntReadRegisterFile );

    /* Grab the brightness scale from it. */
    this->mBrightnessScale = theRegInfo->brightness_scale;
    strcpy ( this->msSubjectName, theRegInfo->name );

    /* If the ps and thickness from our volume does not match the ps and
       thickness in the registration file, complain */
    if ( !FEQUAL(this->mpData->xsize, theRegInfo->in_plane_res) ||
         !FEQUAL(this->mpData->zsize, theRegInfo->slice_thickness) )
    {

      DebugPrint( ("\tWARNING: MRI xsize and zsize are not equal to "
                   "reginfo in_plane_res and slice_thickness:\n"
                   "\t    MRI:        xsize = %f           zsize = %f\n"
                   "\treginfo: in_plane_res = %f slice_thickness = %f\n",
                   this->mpData->xsize, this->mpData->zsize,
                   theRegInfo->in_plane_res, theRegInfo->slice_thickness) );
    }

    /* Get rid of the registration info. */
    StatFreeRegistration ( &theRegInfo );

    /* At this point we scan thru the register file looking for a line
       with a registration keyword on it. This will specify our
       conversion method. If nothing is found, tkregister will be
       used. */
    DebugNote( ("Opening this->msRegistrationFileName") );
    fRegister = fopen( this->msRegistrationFileName, "r" );
    DebugAssertThrowX( (NULL != theRegInfo ),
                       eResult, FunD_tErr_CouldntReadRegisterFile );

    /* Start with our default of FCF. */
    convMethod = FunD_tConversionMethod_FCF;

    DebugNote( ("Looking for conversion keyword") );
    while ( !feof( fRegister ))
    {

      /* read line and look for string */
      fgets( sLine, 1024, fRegister );
      if (feof(fRegister))
      {
        break;
      }
      bGood = sscanf( sLine, "%s", sKeyWord );
      if ( bGood )
      {

        if ( strcmp( sKeyWord, FunD_ksConversionMethod_FCF ) == 0 )
        {
          convMethod = FunD_tConversionMethod_FCF;
        }
        if ( strcmp( sKeyWord, FunD_ksConversionMethod_FFF ) == 0 )
        {
          convMethod = FunD_tConversionMethod_FFF;
        }
        if ( strcmp( sKeyWord, FunD_ksConversionMethod_Round ) == 0 )
        {
          convMethod = FunD_tConversionMethod_Round;
        }
      }
    }

    /* Close the register file. */
    fclose( fRegister );

    break;

  case FunD_tRegistration_Identity:

    /* Default values here. */
    this->mBrightnessScale = 1.0;
    strcpy( this->msSubjectName, "" );
    convMethod = FunD_tConversionMethod_FCF;

    break;
  default:
    break;
  }

  /* Generate our anatomical index -> functional index transform */
  DebugNote(("Calling MRImakeVox2VoxReg"));
  eMRI = MRImakeVox2VoxReg( iAnatomicalVolume->mpMriValues,
                            this->mpData,
                            (int)iType,
                            this->msRegistrationFileName,
                            &this->mIdxToIdxTransform );
  DebugAssertThrowX( (ERROR_NONE == eMRI ),
                     eResult, FunD_tErr_CouldntReadRegisterFile );

  /* Save the registration. */
  Trns_DeepClone( this->mIdxToIdxTransform,
                  &this->mOriginalIdxToIdxTransform );

  /* Set the conversion method. */
  FunD_SetConversionMethod( this, convMethod );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_RestoreRegistration ( mriFunctionalDataRef this )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_RestoreRegistration( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugAssertThrowX( (NULL != this->mOriginalIdxToIdxTransform), eResult,
                     FunD_tErr_ErrorAccessingTransform );

  /* Delete the matrix and clone the backup. */
  Trns_Delete( &this->mIdxToIdxTransform );
  Trns_DeepClone( this->mOriginalIdxToIdxTransform,
                  &this->mIdxToIdxTransform );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}



FunD_tErr FunD_SetClientCoordBounds ( mriFunctionalDataRef this,
                                      int                  inXMin,
                                      int                  inYMin,
                                      int                  inZMin,
                                      int                  inXMax,
                                      int                  inYMax,
                                      int                  inZMax )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction
  (
    ("FunD_SetClientCoordBounds( this=%p, inXMin=%d, "
     "inYMin=%d, inZMin=%d, inXMax=%d, inYMax=%d, inZMax=%d )",
     this, inXMin, inYMin, inZMin, inXMax, inYMax, inZMax)
  );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking max > min") );
  DebugAssertThrowX( (inXMax >= inXMin && inYMax >= inYMin && inZMax>=inZMin),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Set the bounds. */
  this->mClientXMin = inXMin;
  this->mClientYMin = inYMin;
  this->mClientZMin = inZMin;
  this->mClientXMax = inXMax;
  this->mClientYMax = inYMax;
  this->mClientZMax = inZMax;
  this->mbHaveClientBounds = TRUE;

  /* Generate the resampled volume. */
  FunD_ResampleData_ ( this );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SetConversionMethod ( mriFunctionalDataRef this,
                                     FunD_tConversionMethod iMethod )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_SetConversionMethod( this=%p, iMethod=%d )",
                       this, iMethod) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iMethod >= 0 && iMethod < FunD_knNumConversionMethods),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Set the method. */
  this->mConvMethod = iMethod;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetSampleType ( mriFunctionalDataRef this,
                               FunD_tSampleType* oType )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetSampleType( this=%p, "
                       "oType=%p )", this, oType ) );

  DebugNote( ("Verifying volume") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (oType != NULL),
                     eResult, FunD_tErr_InvalidParameter );

  /* return the type */
  DebugNote( ("Returning the type") );
  *oType = this->mSampleType;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SetSampleType ( mriFunctionalDataRef this,
                               FunD_tSampleType iType )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_SetSampleType( this=%p, iType=%d )",
                       this, (int)iType ) );

  DebugNote( ("Verifying volume") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (iType >= 0 && iType < FunD_knNumSampleTypes),
                     eResult, FunD_tErr_InvalidParameter );

  /* Set the sample type. */
  this->mSampleType = iType;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_ClientSpaceIsTkRegRAS ( mriFunctionalDataRef this )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  MATRIX*   identity = NULL;
  Trns_tErr eTrns    = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ClientSpaceIsTkRegRAS( this=%p )", this) );

  DebugNote( ("Verifying volume") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Set the A->RAS part of the transform to identity. */
  DebugNote( ("Checking if mIdxToIdxTransform exists") );
  DebugAssertThrowX( (NULL != this->mIdxToIdxTransform),
                     eResult, FunD_tErr_ErrorAccessingTransform );

  DebugNote( ("Creating identity matrix") );
  identity = MatrixIdentity (4, NULL);
  DebugAssertThrowX( (NULL != identity),
                     eResult, FunD_tErr_CouldntAllocateMatrix );

  DebugNote( ("Copying in transform") );
  eTrns = Trns_CopyAtoRAS (this->mIdxToIdxTransform, identity);
  DebugAssertThrowX( (Trns_tErr_NoErr == eTrns),
                     eResult, FunD_tErr_ErrorAccessingTransform );


  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  if ( NULL != identity )
  {
    MatrixFree (&identity);
  }

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_IsScalar ( mriFunctionalDataRef this,
                          tBoolean* obScalar )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_IsScalar( this=%p, obScalar=%p )",
                       this, obScalar ) );

  DebugNote( ("Verifying volume") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (obScalar != NULL),
                     eResult, FunD_tErr_InvalidParameter );

  /* return the value */
  DebugNote( ("Returning the scalar flag") );
  *obScalar = this->mbScalar;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


/* NOTE: I removed all the DebugEnter/ExitFunction and DebugNote calls
   because they make this call much slower. I wish there was a better
   way to do this. - RKT */
#define BE_SUPA_FAST

FunD_tErr FunD_GetData ( mriFunctionalDataRef this,
                         xVoxelRef            iClientVox,
                         int                  iCondition,
                         int                  iTimePoint,
                         float*               oValue )
{

  FunD_tErr eResult = FunD_tErr_NoError;
  xVoxel    funcIdx;
  float     fValue  = 0;

#ifndef BE_SUPA_FAST
  DebugEnterFunction( ("FunD_GetData( this=%p, iClientVox=%p, "
                       "iCondition=%d, iTimePoint=%d, oValue=%p)",
                       this, iClientVox, iCondition, iTimePoint, oValue) );
#endif

#ifndef BE_SUPA_FAST
  DebugNote( ("Checking parameters") );
#endif
  DebugAssertThrowX( (NULL != this && NULL != iClientVox && NULL != oValue),
                     eResult, FunD_tErr_InvalidParameter );
#ifndef BE_SUPA_FAST
  DebugNote( ("Checking condition") );
#endif
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );
#ifndef BE_SUPA_FAST
  DebugNote( ("Checking time point") );
#endif
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

#ifndef BE_SUPA_FAST
  DebugNote( ("Verifying object") );
#endif
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  if ( !this->mbScalar )
  {
#ifndef BE_SUPA_FAST
    DebugNote( ("Converting client to func idx") );
#endif
    FunD_ConvertClientToFloatFuncIdx_( this, iClientVox, &funcIdx );
  }
  else
  {
    xVoxl_Copy( &funcIdx, iClientVox );
  }

#ifndef BE_SUPA_FAST
  DebugNote( ("Verifying func idx") );
#endif
  eResult = FunD_VerifyFuncIdx_( this, &funcIdx );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

#ifndef BE_SUPA_FAST
  DebugNote( ("Getting value") );
#endif
  FunD_GetValue_ ( this, this->mpData, &funcIdx,
                   iCondition, iTimePoint, &fValue );

#ifndef BE_SUPA_FAST
  DebugNote( ("Setting return value") );
#endif
  *oValue = fValue;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

#ifndef BE_SUPA_FAST
  DebugExitFunction;
#endif

  return eResult;
}

FunD_tErr FunD_GetSampledData ( mriFunctionalDataRef this,
                                xVoxelRef            iClientVox,
                                int                  iCondition,
                                int                  iTimePoint,
                                float*               oValue )
{

  FunD_tErr eResult = FunD_tErr_NoError;
  xVoxel    funcIdx;
  float     fValue  = 0;

#ifndef BE_SUPA_FAST
  DebugEnterFunction( ("FunD_GetSampledData( this=%p, iClientVox=%p, "
                       "iCondition=%d, iTimePoint=%d, oValue=%p)",
                       this, iClientVox, iCondition, iTimePoint, oValue) );
#endif

#ifndef BE_SUPA_FAST
  DebugNote( ("Checking parameters") );
#endif
  DebugAssertThrowX( (NULL != this && NULL != iClientVox && NULL != oValue),
                     eResult, FunD_tErr_InvalidParameter );
#ifndef BE_SUPA_FAST
  DebugNote( ("Checking condition") );
#endif
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );
#ifndef BE_SUPA_FAST
  DebugNote( ("Checking time point") );
#endif
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

#ifndef BE_SUPA_FAST
  DebugNote( ("Verifying object") );
#endif
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  if ( !this->mbScalar )
  {
#ifndef BE_SUPA_FAST
    DebugNote( ("Converting client to func idx") );
#endif
    FunD_ConvertClientToFloatFuncIdx_( this, iClientVox, &funcIdx );
  }
  else
  {
    xVoxl_Copy( &funcIdx, iClientVox );
  }

#ifndef BE_SUPA_FAST
  DebugNote( ("Verifying func idx") );
#endif
  eResult = FunD_VerifyFuncIdx_( this, &funcIdx );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

#ifndef BE_SUPA_FAST
  DebugNote( ("Getting value") );
#endif
  FunD_GetSampledValue_ ( this, this->mpData, &funcIdx,
                          iCondition, iTimePoint, &fValue );

#ifndef BE_SUPA_FAST
  DebugNote( ("Setting return value") );
#endif
  *oValue = fValue;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

#ifndef BE_SUPA_FAST
  DebugExitFunction;
#endif

  return eResult;
}


FunD_tErr FunD_GetDataForAllTimePoints ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iCondition,
    float*               oaValue )
{

  FunD_tErr eResult    = FunD_tErr_NoError;
  xVoxel    funcIdx;
  int       nTimePoint = 0;
  float     fValue     = 0;

  DebugEnterFunction( ("FunD_GetDataForAllTimePoints( this=%p, iClientVox=%p, "
                       "iCondition=%d, oaValue=%p)",
                       this, iClientVox, iCondition, oaValue) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != iClientVox && NULL != oaValue),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking condition") );
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  if ( !this->mbScalar )
  {
    DebugNote( ("Converting client to func idx") );
    FunD_ConvertClientToFloatFuncIdx_( this, iClientVox, &funcIdx );
  }
  else
  {
    xVoxl_Copy( &funcIdx, iClientVox );
  }

  DebugNote( ("Verifying func idx") );
  eResult = FunD_VerifyFuncIdx_( this, &funcIdx );
  DebugAssertQuietThrow( (eResult == FunD_tErr_NoError) );


  for ( nTimePoint = 0; nTimePoint < this->mNumTimePoints; nTimePoint++ )
  {

    DebugNote( ("Getting value") );
    FunD_GetValue_ ( this, this->mpData, &funcIdx,
                     iCondition, nTimePoint, &fValue );

    DebugNote( ("Setting return value for tp %d", nTimePoint) );
    oaValue[nTimePoint] = fValue;
  }

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetDeviation ( mriFunctionalDataRef this,
                              xVoxelRef            iClientVox,
                              int                  iCondition,
                              int                  iTimePoint,
                              float*               oValue )
{

  FunD_tErr eResult     = FunD_tErr_NoError;
  xVoxel    funcIdx;
  float     fCovariance = 0;
  float     fSigma      = 0;
  int       nCovMtx     = 0;

  DebugEnterFunction( ("FunD_GetDeviation( this=%p, iCondition=%d, "
                       "iTimePoint=%d, oValue=%p)", this, iCondition,
                       iTimePoint, oValue) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oValue),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking condition") );
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking time point") );
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* If this is condition 0, just return 0, since there's not error
     data on the null condition. */
  if ( 0 == iCondition )
  {

    DebugNote( ("Setting return value") );
    *oValue = 0;

  }
  else
  {

    DebugNote( ("Getting cov mtx for cond %d tp %d\n",
                iCondition, iTimePoint) );
    nCovMtx = ((iCondition-1) * this->mNumTimePoints) + iTimePoint;
    fCovariance = this->mCovMtx[nCovMtx][nCovMtx];

    if ( !this->mbScalar )
    {
      DebugNote( ("Converting client to func idx") );
      FunD_ConvertClientToFuncIdx_( this, iClientVox, &funcIdx );
    }
    else
    {
      xVoxl_Copy( &funcIdx, iClientVox );
    }

    DebugNote( ("Getting sigma value") );
    FunD_GetSigma_ ( this, &funcIdx, iTimePoint, &fSigma );

    DebugNote( ("Setting return value") );
    *oValue = sqrt( fSigma*fSigma * fCovariance );
  }

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetDeviationForAllTimePoints ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iCondition,
    float*               oaValue )
{

  FunD_tErr eResult    = FunD_tErr_NoError;
  int       nTimePoint = 0;
  xVoxel    funcIdx;
  float     fCovariance = 0;
  float     fSigma      = 0;
  int       nCovMtx     = 0;

  DebugEnterFunction( ("FunD_GetDeviation( this=%p, iCondition=%d, "
                       "oaValue=%p)", this, iCondition, oaValue) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oaValue),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking condition") );
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  if ( !this->mbScalar )
  {
    DebugNote( ("Converting client to func idx") );
    FunD_ConvertClientToFuncIdx_( this, iClientVox, &funcIdx );
  }
  else
  {
    xVoxl_Copy( &funcIdx, iClientVox );
  }

  for ( nTimePoint = 0; nTimePoint < this->mNumTimePoints; nTimePoint++ )
  {

    /* If this is condition 0, just return 0, since there's not error
       data on the null condition. */
    if ( 0 == iCondition )
    {

      DebugNote( ("Setting return value") );
      oaValue[nTimePoint] = 0;

    }
    else
    {

      DebugNote( ("Getting sigma value") );
      FunD_GetSigma_ ( this, &funcIdx, nTimePoint, &fSigma );

      DebugNote( ("Getting cov mtx for cond %d tp %d\n",
                  iCondition, nTimePoint) );
      nCovMtx = ((iCondition-1) * this->mNumTimePoints) + nTimePoint;
      fCovariance = this->mCovMtx[nCovMtx][nCovMtx];

      DebugNote( ("Setting return value") );
      oaValue[nTimePoint] = sqrt( fSigma*fSigma * fCovariance );
    }
  }

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_Smooth ( mriFunctionalDataRef this,
                        int                  iTimePoint,
                        int                  iCondition,
                        float                iSigma )
{

  FunD_tErr eResult = FunD_tErr_NoError;
  MRI*      kernel  = NULL;

  DebugEnterFunction( ("FunD_Smooth( this=%p, iCondition=%d, "
                       "iTimePoint=%d, iSigma=%f)", this, iCondition,
                       iTimePoint, iSigma) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this),eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking condition") );
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking time point") );
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Allocating gaussian kenrel with MRIgaussian1d") );
  kernel = MRIgaussian1d( iSigma, 100 );
  DebugAssertThrowX( (NULL != kernel), eResult, FunD_tErr_CouldntAllocateMRI );

  DebugNote( ("Performing convolution with MRIconvolveGaussian") );
  MRIconvolveGaussian( this->mpData, this->mpData, kernel );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  if ( NULL != kernel )
  {
    MRIfree( &kernel );
  }

  return eResult;
}

FunD_tErr FunD_NormalizeOverAll ( mriFunctionalDataRef this )
{
  FunD_tErr eResult    = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_NormalizeOverAll( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this),eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Performing normalization with MRInormalizeFrames") );
  MRInormalizeFrames( this->mpData );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_ConvertTimePointToSecond ( mriFunctionalDataRef this,
    int                  iTimePoint,
    float*               oSecond )
{

  FunD_tErr eResult         = FunD_tErr_NoError;
  float     timeAtFirstPoint = 0;

  DebugEnterFunction( ("FunD_ConvertTimePointToSecond( this=%p, "
                       "iTimePoint=%d, oSecond=%p)", this, iTimePoint,
                       oSecond) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oSecond),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking time point") );
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Calc the time second. */
  timeAtFirstPoint = -((float)(this->mNumPreStimTimePoints) *
                       this->mTimeResolution);
  DebugNote( ("Setting return value") );
  *oSecond = timeAtFirstPoint + ((float)iTimePoint * this->mTimeResolution);

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_ConvertSecondToTimePoint ( mriFunctionalDataRef this,
    float                iSecond,
    int*                 oTimePoint )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_ConvertSecondToTimePoint( this=%p, iSecond=%f, "
                       "oTimePoint=%p)", this, iSecond, oTimePoint) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oTimePoint),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  // calc the time point
  DebugNote( ("Setting return value") );
  *oTimePoint = (iSecond/this->mTimeResolution) + this->mNumPreStimTimePoints;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SetTimeResolution ( mriFunctionalDataRef this,
                                   float                iTimeResolution )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_SetTimeResolution( this=%p, iTimeResolution=%f )",
                       this, iTimeResolution) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iTimeResolution > 0),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Set the time resolution. */
  this->mTimeResolution = iTimeResolution;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SetNumPreStimTimePoints ( mriFunctionalDataRef this,
    int                  iNumPoints )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_SetNumPreStimTimePoints( this=%p, iNumPoints=%d)",
                       this, iNumPoints) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iNumPoints >= 0 && iNumPoints < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Set the number of prestim time points. */
  this->mNumPreStimTimePoints = iNumPoints;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetSubjectName ( mriFunctionalDataRef this,
                                char*                oSubjectName )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetSubjectName( this=%p, oSubjectName=%s)",
                       this, oSubjectName) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oSubjectName),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Copying subject name") );
  strcpy( oSubjectName, this->msSubjectName );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetNumTimePoints ( mriFunctionalDataRef this,
                                  int*                 oNumTimePoints )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetNumTimePoints( this=%p, oNumTimePoints=%p)",
                       this, oNumTimePoints) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oNumTimePoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oNumTimePoints = this->mNumTimePoints;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetNumConditions ( mriFunctionalDataRef this,
                                  int*                 oNumConditions )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetNumConditions( this=%p, oNumConditions=%p)",
                       this, oNumConditions) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oNumConditions),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oNumConditions = this->mNumConditions;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_GetTimeResolution ( mriFunctionalDataRef this,
                                   float*               oTimeResolution )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_Get( this=%p, oTimeResolution=%p)",
                       this, oTimeResolution) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oTimeResolution),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oTimeResolution = this->mTimeResolution;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_GetNumPreStimTimePoints ( mriFunctionalDataRef this,
    int*                 oNumPoints )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetNumPreStimTimePoints( this=%p, oNumPoints=%p)",
                       this, oNumPoints) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oNumPoints),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oNumPoints = this->mNumPreStimTimePoints;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetValueRange ( mriFunctionalDataRef this,
                               float*               oMin,
                               float*               oMax )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_GetValueRange( this=%p, oMin=%p, oMax=%p)",
                       this, oMin, oMax) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oMin && NULL != oMax),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oMin = this->mMinValue;
  *oMax = this->mMaxValue;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_IsErrorDataPresent ( mriFunctionalDataRef this,
                                    tBoolean*            oPresent )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_IsErrorDataPresent( this=%p, oPresent=%p)",
                       this, oPresent) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oPresent),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Return the value. */
  DebugNote( ("Setting return value") );
  *oPresent = this->mbErrorDataPresent;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_GetBoundsInClientSpace ( mriFunctionalDataRef this,
                                        xVoxelRef            oBeginCorner,
                                        xVoxelRef            oEndCorner )
{

  FunD_tErr eResult   = FunD_tErr_NoError;
  xVoxel    curIdx;
  xVoxel    curClient;
  xVoxel    minClient;
  xVoxel    maxClient;

  DebugEnterFunction( ("FunD_GetBoundsInClientSpace( this=%p, "
                       "oBeginCorner=%p, oEndCorner=%p)",
                       this, oBeginCorner, oEndCorner) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oBeginCorner &&
                      NULL != oEndCorner),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugNote( ("Initializing min and max voxels") );
  xVoxl_Set( &minClient, 9999, 9999, 9999 );
  xVoxl_Set( &maxClient, -9999, -9999, -9999 );

  /* For all of our indices, convert it to client space and compare to
     the min/max. Save the most min/max.  */
  xVoxl_Set( &curIdx, 0, 0, 0 );
  while ( xVoxl_IncrementUntilLimits( &curIdx,
                                      this->mpData->width,
                                      this->mpData->height,
                                      this->mpData->depth ) )
  {

    DebugNote( ("Converting to client space") );
    FunD_ConvertFuncIdxToClient_( this, &curIdx, &curClient );

    DebugNote( ("Setting min voxel") );
    xVoxl_Set( &minClient,
               min( xVoxl_GetX(&minClient), xVoxl_GetX(&curClient) ),
               min( xVoxl_GetY(&minClient), xVoxl_GetY(&curClient) ),
               min( xVoxl_GetZ(&minClient), xVoxl_GetZ(&curClient) ) );

    DebugNote( ("Setting max voxel") );
    xVoxl_Set( &maxClient,
               max( xVoxl_GetX(&maxClient), xVoxl_GetX(&curClient) ),
               max( xVoxl_GetY(&maxClient), xVoxl_GetY(&curClient) ),
               max( xVoxl_GetZ(&maxClient), xVoxl_GetZ(&curClient) ) );
  }

  DebugNote( ("Setting min return voxel") );
  xVoxl_Copy( oBeginCorner, &minClient );
  DebugNote( ("Setting max return voxel") );
  xVoxl_Copy( oEndCorner, &maxClient );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_GetValueAtPercentile ( mriFunctionalDataRef this,
                                      float                iPercentile,
                                      float*               oValue )
{


  FunD_tErr eResult             = FunD_tErr_NoError;
  int       numValues           = 0;
  int       targetCount         = 0;
  int       sum                 = 0;
  int       nBin                = 0;
  float     value               = 0;

  DebugEnterFunction( ("FunD_GetValueAtPercentile( this=%p, iPercentile=%f, "
                       "oValue=%p)", this, iPercentile, oValue) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != oValue),
                     eResult, FunD_tErr_InvalidParameter );
  DebugNote( ("Checking 0 <= percentile <= 100  ") );
  DebugAssertThrowX( (iPercentile >= 0 && iPercentile <= 100),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* If we don't have our frequencies yet, calc them now with 100 bins. */
  if ( NULL == this->mFrequencies )
  {
    eResult = FunD_CalcFrequencies_( this, FunD_knNumFrequencyBins );
    DebugAssertThrow( (eResult == FunD_tErr_NoError) );
  }

  /* Calc the total number of values. Calc the precent of that. */
  numValues =
    this->mpData->width * this->mpData->height * this->mpData->depth *
    this->mNumConditions * this->mNumTimePoints;

  targetCount = (float)numValues * iPercentile/100.0;

  /* Start from the beginning. Add the sum of each bin until we get to
     the number we're looking for. */
  sum = 0;
  nBin = 0;
  while ( sum < targetCount && nBin < this->mNumBins )
  {
    sum += this->mFrequencies[nBin];
    nBin++;
  }

  /* Calc the number that is at the beginning of that bin. */
  value = this->mMinValue +
          ( ((this->mMaxValue - this->mMinValue + 1) * nBin) / this->mNumBins );

  /* Return it. */
  DebugNote( ("Setting return value") );
  *oValue = value;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

/* This function doesn't seem to work very well but is left here for
   reference . */
#if 0
FunD_tErr FunD_GetPercentileOfValue ( mriFunctionalDataRef this,
                                      float                inValue,
                                      float*               outPercentile )
{

  FunD_tErr theErr = FunD_tErr_NoError;
  int theTargetBin, theBin;
  int theSum, theNumValues;
  float thePercentile;

  /* If we don't have our frequencies yet, calc them now with 100 bins. */
  if ( NULL == this->mFrequencies )
  {
    theErr = FunD_CalcFrequencies_( this, FunD_knNumFrequencyBins );
    if ( FunD_tErr_NoError != theErr )
    {
      goto error;
    }
  }

  /* Check to make sure numBins is 100 here? */

  /* Calculate the bin in which it should be. This is done with
     floating math because they want a percentile even tho we call it
     a bin and use the integer version. */
  theTargetBin = (inValue - this->mMinValue) /
                 (((this->mMaxValue - this->mMinValue)+1.0) / (float)this->mNumBins);
  theTargetBin = ((inValue - this->mMinValue) * (float)this->mNumBins) /
                 (this->mMaxValue - this->mMinValue + 1.0);

  /* Walk through the frequency array until we hit the bin we got as
     an answer. Add up the counts along the way. */
  theSum = 0;
  theBin = 0;
  for ( theBin = 0; theBin <= theTargetBin; theBin++ )
  {
    theSum += this->mFrequencies[theBin];
  }

  /* Divide the counts by the total number of values to get our
     percentile. */
  theNumValues = this->mNumConditions * this->mNumTimePoints *
                 this->mNumRows * this->mNumCols * this->mNumSlices;
  thePercentile = (float)theSum / (float)theNumValues * 100.0;

  *outPercentile = thePercentile;

  goto cleanup;
error:

  if ( FunD_tErr_NoError != theErr )
  {
    DebugPrint( ("Error in FunD_GetPercentileOfValue (%d): %s\n",
                 theErr, FunD_GetErrorString(theErr) ) );
  }

cleanup:

  return theErr;
}
#endif


FunD_tErr FunD_CalcFDRThreshold ( mriFunctionalDataRef this,
                                  int                  iCondition,
                                  int                  iTimePoint,
                                  int                  iSign,
                                  float                iRate,
                                  MRI*                 iMaskVolume,
                                  float*               oThresholdMin )
{

  FunD_tErr eResult       = FunD_tErr_NoError;
  int       eMRI          = ERROR_NONE;
  int       nFrame        = 0;
  double    newMin        = 0;
  MRI*      pLocalMaskVol = NULL;
  MATRIX*   localToClient = NULL;

  DebugEnterFunction( ("FunD_CalcFDRThreshold( this=%p, iCondition=%d, "
                       "iTimePoint=%d, iSign=%d, iRate=%f, iMaskVolume=%p,"
                       "oThresholdMin=%p)", this, iCondition, iTimePoint,
                       iSign, iRate, iMaskVolume, oThresholdMin) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* We need to transform the mask into our space. Allocate an empty 1
     frame volume with the same dimesnions as our data, then transform
     the mask volume into it. To transform, we need the target->source
     transform, which is local->client, which is the inverse of our
     mIdxToIdxTransform, so create that too.*/
  if ( NULL != iMaskVolume )
  {
    DebugNote( ("Allocating local mask volume") );
    pLocalMaskVol = MRIalloc( this->mpData->width, this->mpData->height,
                              this->mpData->depth, this->mpData->type );
    DebugAssertThrowX( (NULL != pLocalMaskVol),
                       eResult, FunD_tErr_CouldntAllocateVolume );

    DebugNote( ("Copying header") );
    MRIcopyHeader( this->mpData, pLocalMaskVol );

    DebugNote( ("Creating inverse of transform (for local->client)") );
    localToClient = MatrixInverse( this->mIdxToIdxTransform->mAtoB, NULL );
    DebugAssertThrowX( (NULL != localToClient),
                       eResult, FunD_tErr_CouldntAllocateMatrix );

    DebugNote( ("Converting client mask to local mask") );
    eMRI = MRIvol2Vol( iMaskVolume, pLocalMaskVol, localToClient,
                       SAMPLE_NEAREST, 0 );
    DebugAssertThrowX( (ERROR_NONE == eMRI),
                       eResult, FunD_tErr_CouldntTransformMask );
  }

  /* Get a frame number. */
  FunD_GetDataFrameNumber( iCondition, iTimePoint, &nFrame );

  /* Do the FDR */
  DebugNote( ("Running MRIfdr2vwth") );
  eMRI = MRIfdr2vwth( this->mpData,
                      nFrame,        /* 0-based frame number of input vol */
                      iRate,         /* rate 0-1 */
                      iSign,         /* 0=all, +1=only pos, -1=only neg */
                      TRUE,          /* interp vals as -log10(p) */
                      pLocalMaskVol, /* optional mask volume */
                      &newMin,       /* voxel-wise thresh (vwth) 0-1,
                                        or if log10flag, -log10(vwth) */
                      NULL );        /* option thresholded output */
  DebugAssertThrowX( (ERROR_NONE == eMRI),
                     eResult, FunD_tErr_ErrorPerformingFDR );


  /* Return the new min. */
  DebugNote( ("Returning threshold min") );
  *oThresholdMin = newMin;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  if ( NULL != pLocalMaskVol )
  {
    DebugNote( ("Freeing local mask vol") );
    MRIfree( &pLocalMaskVol );
  }
  if ( NULL != localToClient )
  {
    DebugNote( ("Freeing local to client transform") );
    MatrixFree( &localToClient );
  }

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SaveRegistration ( mriFunctionalDataRef this )
{

  FunD_tErr eResult = FunD_tErr_NoError;
  fMRI_REG* regInfo = NULL;
  struct stat fileInfo;
  char      sFileName[256];
  char      sBackupFileName[256];
  char      sConversionMethod[256];
  FILE*     pBackupFile;
  FILE*     pFile;
  char      data;
  int       nBackup = 1;

  DebugEnterFunction( ("FunD_SaveRegistration( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* make a reg info struct */
  DebugNote( ("Allocating registration") );
  regInfo = (fMRI_REG*) malloc( sizeof( fMRI_REG ));
  DebugAssertThrowX( (NULL != regInfo),
                     eResult, FunD_tErr_CouldntAllocateStorage );

  /* Fill out the fmri registration struct. Notice we used the ps and
     st values we first got from the registration, instead of the
     values in the MRI. */
  DebugNote( ("Setting registration info") );
  regInfo->in_plane_res     = this->mpData->ps;
  regInfo->slice_thickness  = this->mpData->thick;
  regInfo->brightness_scale = this->mBrightnessScale;
  DebugNote( ("Copying subject name") );
  strcpy( regInfo->name, this->msSubjectName );

  /* Copy it into the regInfo structure. */
  DebugNote( ("Copying registration matrix into registration info") );
  regInfo->mri2fmri =
    MatrixCopy( this->mIdxToIdxTransform->mARAStoBRAS, NULL );
  DebugNote( ("Inversing registration matrix into registration info") );
  regInfo->fmri2mri = MatrixInverse( regInfo->mri2fmri, NULL );

  /* If a registration already exists... */
  if ( stat( sFileName, &fileInfo ) != 0 )
  {

    /* We need to backup the file. Generate a new backup filename */
    do
    {

      /* Keep appending an increasing number to the end of the
         original filename. Check if the newly named file exists. It
         it doesn't exist, the while loop will exit and we've found
         our new backup filename. */
      sprintf( sBackupFileName, "%s.%d", sFileName, nBackup++ );
      DebugNote( ("Opening file %s", sBackupFileName) );

    }
    while ( stat( sFileName, &fileInfo ) != 0 );

    /* Copy the registration file to backup.*/
    DebugNote( ("Opening file %s", sFileName) );
    pFile = fopen( sFileName, "r" );

    DebugNote( ("Opening backup file %s", sBackupFileName) );
    pBackupFile = fopen( sBackupFileName, "w" );

    DebugNote( ("Copying bytes to backup file") );
    while (!feof( pFile ) )
    {
      fread( &data, sizeof(data), 1, pFile );
      if (feof(pFile))
      {
        break;
      }
      fwrite( &data, sizeof(data), 1, pBackupFile );
    }

    DebugNote( ("Closing backup file") );
    fclose( pBackupFile );
    DebugNote( ("Closing file") );
    fclose( pFile );
  }

  /* Write it to disk */
  DebugNote( ("Writing registration with StatWriteRegistration") );
  StatWriteRegistration( regInfo, sFileName );

  DebugNote( ("Freeing registration with StatFreeRegistration") );
  StatFreeRegistration( &regInfo );

  /* Append the conversion type to the registration file. */
  switch ( this->mConvMethod )
  {
  case FunD_tConversionMethod_FFF:
    strcpy( sConversionMethod, FunD_ksConversionMethod_FFF );
    break;
  case FunD_tConversionMethod_Round:
    strcpy( sConversionMethod, FunD_ksConversionMethod_Round );
    break;
  case FunD_tConversionMethod_FCF:
    strcpy( sConversionMethod, FunD_ksConversionMethod_FCF );
    break;
  default:
    break;
  }
  DebugNote( ("Opening %s") );
  pFile = fopen( sFileName, "a" );
  if ( NULL != pFile )
  {
    DebugNote( ("Printing conversion method to file") );
    fprintf( pFile, "%s\n", sConversionMethod );
    DebugNote( ("Closing file") );
    fclose( pFile );
  }

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_SetRegistrationToIdentity ( mriFunctionalDataRef this )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  MATRIX*   identity = NULL;

  DebugEnterFunction( ("FunD_SetRegistrationToIdentity( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Create an identity matrix and copy it into the registration
     transform. */
  DebugNote( ("Creting identity matrix") );
  identity = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixIdentity( 4, identity );

  DebugNote( ("Copying identity to transform") );
  Trns_CopyARAStoBRAS( this->mIdxToIdxTransform, identity );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  if ( NULL != identity )
  {
    DebugNote( ("Freeing identity matrix") );
    MatrixFree( &identity );
  }

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_ApplyTransformToRegistration ( mriFunctionalDataRef this,
    MATRIX*           iTransform )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;
  MATRIX*   invTransform = NULL;

  DebugEnterFunction( ("FunD_ApplyTransformToRegistration( this=%p )"
                       "iTransform=%p", this, iTransform) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != iTransform),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* We apply the inverse of the transform to our idx to idx
     transform. */
  DebugNote( ("Getting inverse of matrix") );
  invTransform = MatrixInverse( iTransform, NULL );

  DebugNote( ("Applying matrix to transform") );
  eTransform = Trns_ApplyTransform( this->mIdxToIdxTransform, invTransform );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform ),
                     eResult, FunD_tErr_ErrorAccessingTransform );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  if ( NULL != invTransform )
  {
    DebugNote( ("Freeing inverse matrix") );
    MatrixFree( &invTransform );
  }

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_TranslateRegistration ( mriFunctionalDataRef this,
                                       float                ifDistance,
                                       tAxis                iAxis )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_TranslateRegistration( this=%p )"
                       "ifDistance=%f, iAxis=%d", this, ifDistance, iAxis) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iAxis >= 0 && iAxis < knNumAxes),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* do the inverse of this action */
  DebugNote( ("Translating transformation object") );
  eTransform = Trns_Translate( this->mIdxToIdxTransform, -ifDistance, iAxis );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform ),
                     eResult, FunD_tErr_ErrorAccessingTransform );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_RotateRegistration ( mriFunctionalDataRef this,
                                    float                ifDegrees,
                                    tAxis                iAxis,
                                    xVoxelRef            iCenterFuncRAS )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;
  float     fX         = 0;
  float     fY         = 0;
  float     fZ         = 0;

  DebugEnterFunction( ("FunD_RotateRegistration( this=%p, ifDegrees=%f, "
                       "iAxis=%d iCenterFuncRAS=%p", this, ifDegrees,
                       iAxis, iCenterFuncRAS) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this && NULL != iCenterFuncRAS),
                     eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iAxis >= 0 && iAxis < knNumAxes),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* get the center point */
  DebugNote( ("Getting center coordinates") );
  fX = xVoxl_GetFloatX( iCenterFuncRAS );
  fY = xVoxl_GetFloatY( iCenterFuncRAS );
  fZ = xVoxl_GetFloatZ( iCenterFuncRAS );

  /* first translate to the center we're rotating around */
  DebugNote( ("Translating registration by negative center") );
  FunD_TranslateRegistration( this, -fX, tAxis_X );
  FunD_TranslateRegistration( this, -fY, tAxis_Y );
  FunD_TranslateRegistration( this, -fZ, tAxis_Z );

  /* do the inverse of this action */
  DebugNote( ("Rotating transformation object by %f degrees around %d",
              -ifDegrees, iAxis) );
  eTransform = Trns_Rotate( this->mIdxToIdxTransform, -ifDegrees, iAxis );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform ),
                     eResult, FunD_tErr_ErrorAccessingTransform );

  /* translate back */
  DebugNote( ("Translating registration by center") );
  FunD_TranslateRegistration( this, fX, tAxis_X );
  FunD_TranslateRegistration( this, fY, tAxis_Y );
  FunD_TranslateRegistration( this, fZ, tAxis_Z );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_ScaleRegistration ( mriFunctionalDataRef this,
                                   float                ifFactor,
                                   tAxis                iAxis )
{

  FunD_tErr eResult  = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ScaleRegistration( this=%p, ifFactor=%f, "
                       "iAxis=%d", this, ifFactor, iAxis) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this),
                     eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iAxis >= 0 && iAxis < knNumAxes),
                     eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* do the inverse of this action */
  DebugNote( ("Scaling the transformation object") );
  eTransform = Trns_Scale( this->mIdxToIdxTransform, 1.0 / ifFactor, iAxis );
  DebugAssertThrowX( (Trns_tErr_NoErr == eTransform ),
                     eResult, FunD_tErr_ErrorAccessingTransform );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_DebugPrint ( mriFunctionalDataRef this )
{

  FunD_tErr eResult  = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_DebugPrint( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  DebugPrint( ("mriFunctionalData:\n" ) );
  DebugPrint( ("\tmsFileName=%s\n", this->msFileName) );
  DebugPrint( ("\tmsRegistrationFileName=%s\n", this->msRegistrationFileName));
  DebugPrint( ("\tmsSubjectName=%s\n", this->msSubjectName) );
  DebugPrint( ("\tmpData->ps=%f\n", this->mpData->ps) );
  DebugPrint( ("\tmpData->thick=%f\n", this->mpData->thick) );
  DebugPrint( ("\tmpData->width=%d\n", this->mpData->width) );
  DebugPrint( ("\tmpData->height=%d\n", this->mpData->height) );
  DebugPrint( ("\tmpData->depth=%d\n", this->mpData->depth) );
  DebugPrint( ("\tmConvMethod=%d\n", this->mConvMethod) );
  DebugPrint( ("\tmNumTimePoints=%d\n", this->mNumTimePoints) );
  DebugPrint( ("\tmNumConditions=%d\n", this->mNumConditions) );
  DebugPrint( ("\tmbNullConditionPresent=%d\n", this->mbNullConditionPresent));
  DebugPrint( ("\tmMaxValue=%f\n", this->mMaxValue) );
  DebugPrint( ("\tmMinValue=%f\n", this->mMinValue) );
  DebugPrint( ("\tmTimeResolution=%f\n", this->mTimeResolution) );
  DebugPrint( ("\tmNumPreStimTimePoints=%d\n", this->mNumPreStimTimePoints) );
  DebugPrint( ("\tmbErrorDataPresent=%d\n", this->mbErrorDataPresent) );
  DebugPrint( ("\tmNumBins=%d\n", this->mNumBins) );
  DebugPrint( ("\n") );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_ResampleData_ ( mriFunctionalDataRef this )
{

  FunD_tErr eResult    = FunD_tErr_NoError;
  MRI*      volume     = NULL;
  int       nTimePoint = 0;
  int       nCondition = 0;
  xVoxel    begin;
  xVoxel    end;
  xVoxel    clientIdx;
  xVoxel    funcIdx;
  xVoxel    resampleIdx;
  float     value      = 0;

  DebugEnterFunction( ("FunD_ResampleData_( this=%p )", this) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Make sure we have the bounds. */
  DebugAssertThrowX( (this->mbHaveClientBounds),
                     eResult, FunD_tErr_InvalidParameter );

  /* Make an MRI with the same dimensions as the client. Use the type
     and number of frames of our data volume. */
  DebugNote( ("Creating volume") );
  volume = MRIallocSequence( this->mClientXMax - this->mClientXMin,
                             this->mClientYMax - this->mClientYMin,
                             this->mClientZMax - this->mClientZMin,
                             this->mpData->type, this->mpData->nframes );
  DebugAssertThrowX( (NULL != volume), eResult, FunD_tErr_CouldntAllocateMRI );

  DebugNote( ("Getting bounds") );
  FunD_GetBoundsInClientSpace( this, &begin, &end );

  /* For each timepoint and condition... */
  for ( nCondition = 0; nCondition < this->mNumConditions; nCondition++ )
  {
    for ( nTimePoint = 0; nTimePoint < this->mNumTimePoints; nTimePoint++ )
    {

      /* For each voxel... */
      xVoxl_Copy( &clientIdx, &begin );
      while ( xVoxl_IncrementWithMinsUntilLimits( &clientIdx,
              xVoxl_GetX(&begin),
              xVoxl_GetY(&begin),
              xVoxl_GetX(&end),
              xVoxl_GetY(&end),
              xVoxl_GetZ(&end) ) )
      {

        DebugPrint( ("\rResampling: %d %d %d", nCondition, nTimePoint,
                     xVoxl_GetZ( &clientIdx ) ) );

        /* Get the func idx. */
        if ( !this->mbScalar )
        {
          DebugNote( ("Converting client to func idx") );
          FunD_ConvertClientToFuncIdx_( this, &clientIdx, &funcIdx );
        }
        else
        {
          xVoxl_Copy( &funcIdx, &clientIdx );
        }
        eResult = FunD_VerifyFuncIdx_( this, &funcIdx );
        if ( FunD_tErr_NoError == eResult )
        {
          FunD_GetValue_ ( this, this->mpData, &funcIdx,
                           nCondition, nTimePoint, &value );
        }
        else
        {
          continue;
        }

        /* Set the value. */
        xVoxl_Set( &resampleIdx,
                   xVoxl_GetX(&clientIdx) - this->mClientXMin,
                   xVoxl_GetY(&clientIdx) - this->mClientYMin,
                   xVoxl_GetZ(&clientIdx) - this->mClientZMin );
        FunD_SetValue_ ( this, volume, &resampleIdx,
                         nCondition, nTimePoint, value );

      }
    }
  }

  fprintf( stdout, "\rResampling: 100%% done.       \n" );

  if ( NULL != this->mpResampledData )
  {
    MRIfree( &this->mpResampledData );
  }
  this->mpResampledData = volume;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );

  if ( NULL != volume )
  {
    MRIfree( &volume );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}


FunD_tErr FunD_CalcFrequencies_ ( mriFunctionalDataRef this,
                                  int                  iNumBins )
{

  FunD_tErr eResult      = FunD_tErr_NoError;
  int*      aFrequencies = NULL;
  int       nCondition   = 0;
  int       nTimePoint   = 0;
  xVoxel    curIdx;
  int       nBin         = 0;
  float     value        = 0;

  DebugEnterFunction( ("FunD_CalcFrequencies_( this=%p, iNumBins=%d )",
                       this, iNumBins) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (iNumBins > 0), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Verifying object") );
  eResult = FunD_Verify( this );
  DebugAssertThrow( (eResult == FunD_tErr_NoError) );

  /* Allocate the storage. */
  DebugNote( ("Allocating storage for new frequencies") );
  aFrequencies = calloc( iNumBins, sizeof(int) );
  DebugAssertThrowX( (NULL != aFrequencies),
                     eResult, FunD_tErr_CouldntAllocateStorage );

  for ( nCondition = 0; nCondition < this->mNumConditions; nCondition++ )
  {
    for ( nTimePoint = 0; nTimePoint < this->mNumTimePoints; nTimePoint++ )
    {
      xVoxl_Set( &curIdx, 0, 0, 0 );
      while ( xVoxl_IncrementUntilLimits( &curIdx,
                                          this->mpData->width,
                                          this->mpData->height,
                                          this->mpData->depth ) )
      {

        FunD_GetValue_( this, this->mpData,
                        &curIdx, nCondition, nTimePoint, &value );

        /* Calculate the bin in which it should be. */
        nBin = (value - this->mMinValue) /
               (((this->mMaxValue - this->mMinValue)+1.0) / (float)iNumBins);

        /* Inc the count in that bin. */
        aFrequencies[nBin]++;
      }
    }
  }

  /* If we already have frequencies, free them. */
  if ( NULL != this->mFrequencies )
  {
    DebugNote( ("Deleting exisiting frequencies") );
    free( this->mFrequencies );
  }

  DebugNote( ("Assigning new frequencies") );
  this->mFrequencies = aFrequencies;
  this->mNumBins     = iNumBins;

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );

  if ( NULL != aFrequencies )
  {
    DebugNote( ("Deleting temp frequencies") );
    free( aFrequencies );
  }

  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

char* FunD_GetErrorString ( FunD_tErr iErr )
{

  if ( !(iErr >= 0 && iErr < FunD_tErr_knNumErrorCodes) )
  {
    iErr = FunD_tErr_InvalidErrorCode;
  }

  return (char*)(FunD_ksaErrorString[iErr]);
}

#ifndef FUND_USE_MACROS

void FunD_ConvertClientToFuncIdx_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iClientVox,
                                    xVoxelRef            oFuncIdx )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertClientToFuncIdx_( this=%p, iClientVox=%p, "
                       "oFuncIdx=%p", this, iClientVox, oFuncIdx) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow( (NULL != this && NULL != iClientVox && NULL != oFuncIdx ));

  DebugNote( ("Copying client voxel to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iClientVox );

  DebugNote( ("Converting client to func idx") );
  eTransform = Trns_ConvertAtoB( this->mIdxToIdxTransform,
                                 &this->mTmpVoxel1, &this->mTmpVoxel2 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  /* do the proper conversion method */
  switch ( this->mConvMethod )
  {
  case FunD_tConversionMethod_FFF:
    DebugNote( ("Using FFF conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  case FunD_tConversionMethod_Round:
    DebugNote( ("Using Round conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    rint(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    rint(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    rint(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  case FunD_tConversionMethod_FCF:
    DebugNote( ("Using FCF conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    ceil(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  default:
    break;
  }

  DebugCatch;
  xVoxl_Set( oFuncIdx, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}

void FunD_ConvertClientToFloatFuncIdx_ ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    xVoxelRef            oFuncIdx )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertClientToFloatFuncIdx_( this=%p, iClientVox=%p, "
                       "oFuncIdx=%p", this, iClientVox, oFuncIdx) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow( (NULL != this && NULL != iClientVox && NULL != oFuncIdx ));

  DebugNote( ("Copying client voxel to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iClientVox );

  DebugNote( ("Converting client to func idx") );
  eTransform = Trns_ConvertAtoB( this->mIdxToIdxTransform,
                                 &this->mTmpVoxel1, oFuncIdx );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  DebugCatch;
  xVoxl_Set( oFuncIdx, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}

#endif


void FunD_ConvertClientToFuncRAS_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iClientVox,
                                    xVoxelRef            oFuncRAS )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertClientToFuncRAS_( this=%p, iClientVox=%p, "
                       "oFuncRAS=%p", this, iClientVox, oFuncRAS) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow((NULL != this && NULL != iClientVox && NULL != oFuncRAS ));

  DebugNote( ("Copying client voxel to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iClientVox );

  DebugNote( ("Converting client to func idx") );
  eTransform = Trns_ConvertAtoB( this->mIdxToIdxTransform,
                                 &this->mTmpVoxel1, &this->mTmpVoxel2 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  DebugNote( ("Converting func idx to func RAS") );
  eTransform = Trns_ConvertBtoRAS( this->mIdxToIdxTransform,
                                   &this->mTmpVoxel2, &this->mTmpVoxel1 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  DebugNote( ("Copying into return voxel") );
  xVoxl_Copy( oFuncRAS, &this->mTmpVoxel1 );

  DebugCatch;
  xVoxl_Set( oFuncRAS, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}


void FunD_ConvertFuncIdxToClient_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iFuncIdx,
                                    xVoxelRef            oClientVox )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertFuncIdxToClient_( this=%p, iFuncIdx=%p, "
                       "oClientVox=%p", this, iFuncIdx, oClientVox) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow((NULL != this && NULL != iFuncIdx && NULL != oClientVox ));

  DebugNote( ("Copying func idx voxel to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iFuncIdx );

  DebugNote( ("Converting func idx to client") );
  eTransform = Trns_ConvertBtoA( this->mIdxToIdxTransform,
                                 &this->mTmpVoxel1, &this->mTmpVoxel2 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  DebugNote( ("Copying into return voxel") );
  xVoxl_Copy( oClientVox, &this->mTmpVoxel2 );

  DebugCatch;
  xVoxl_Set( oClientVox, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}

void FunD_ConvertRASToFuncIdx_ ( mriFunctionalDataRef this,
                                 xVoxelRef            iRAS,
                                 xVoxelRef            oFuncIdx )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertRASToFuncIdx_( this=%p, iRAS=%p, "
                       "oFuncIdx=%p", this, iRAS, oFuncIdx) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow((NULL != this && NULL != iRAS && NULL != oFuncIdx ));

  DebugNote( ("Copying RAS to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iRAS );

  DebugNote( ("Converting RAS to func idx") );
  eTransform = Trns_ConvertBRAStoB( this->mIdxToIdxTransform,
                                    &this->mTmpVoxel1, &this->mTmpVoxel2 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  /* do the proper conversion method */
  switch ( this->mConvMethod )
  {
  case FunD_tConversionMethod_FFF:
    DebugNote( ("Using FFF conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  case FunD_tConversionMethod_Round:
    DebugNote( ("Using Round conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    rint(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    rint(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    rint(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  case FunD_tConversionMethod_FCF:
    DebugNote( ("Using FCF conversion") );
    xVoxl_SetFloat( oFuncIdx,
                    floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),
                    ceil(xVoxl_GetFloatY(&this->mTmpVoxel2)),
                    floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );
    break;
  default:
    break;
  }

  DebugCatch;
  xVoxl_Set( oFuncIdx, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}

void FunD_ConvertFuncIdxToFuncRAS_ ( mriFunctionalDataRef this,
                                     xVoxelRef            iFuncIdx,
                                     xVoxelRef            oFuncRAS )
{

  Trns_tErr eTransform = Trns_tErr_NoErr;

  DebugEnterFunction( ("FunD_ConvertFuncIdxToFuncRAS_( this=%p, iFuncIdx=%p, "
                       "oFuncRAS=%p", this, iFuncIdx, oFuncRAS) );

  DebugNote( ("Checking parameters") );
  DebugAssertThrow((NULL != this && NULL != iFuncIdx && NULL != oFuncRAS ));

  DebugNote( ("Copying func idx to temp1") );
  xVoxl_Copy( &this->mTmpVoxel1, iFuncIdx );

  DebugNote( ("Covnerting func idx to func ras") );
  eTransform = Trns_ConvertBtoRAS( this->mIdxToIdxTransform,
                                   &this->mTmpVoxel1, &this->mTmpVoxel2 );
  DebugAssertThrow( (Trns_tErr_NoErr == eTransform) );

  DebugNote( ("Copying to return voxel") );
  xVoxl_Copy( oFuncRAS, &this->mTmpVoxel2 );

  DebugCatch;
  xVoxl_Set( oFuncRAS, -1, -1, -1 );
  EndDebugCatch;

  DebugExitFunction;
}

#ifndef FUND_USE_MACROS

void FunD_GetValue_ ( mriFunctionalDataRef this,
                      MRI*                 iData,
                      xVoxelRef            iIdx,
                      int                  inCondition,
                      int                  inTimePoint,
                      float*               oValue )
{

  int nFrame = 0;

  FunD_GetDataFrameNumber(inCondition,inTimePoint,&nFrame);

  MRIsampleVolumeFrameType( iData,
                            xVoxl_GetFloatX(iIdx),
                            xVoxl_GetFloatY(iIdx),
                            xVoxl_GetFloatZ(iIdx),
                            nFrame,
                            SAMPLE_NEAREST,
                            &this->mTmpReal );

  *oValue = this->mTmpReal;

#if 0
  switch ( iData->type )
  {
  case MRI_UCHAR:
    *oValue = MRIseq_vox( iData, xVoxl_GetX(iIdx),
                          xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame);
    break;
  case MRI_INT:
    *oValue = MRIIseq_vox( iData, xVoxl_GetX(iIdx),
                           xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame);
    break;
  case MRI_LONG:
    *oValue = MRILseq_vox( iData, xVoxl_GetX(iIdx),
                           xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame);
    break;
  case MRI_FLOAT:
    *oValue = MRIFseq_vox( iData, xVoxl_GetX(iIdx),
                           xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame);
    break;
  case MRI_SHORT:
    *oValue = MRISseq_vox( iData, xVoxl_GetX(iIdx),
                           xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame);
    break;
  default:
    *oValue = 0;
    break ;
  }
#endif
}

void FunD_GetSampledValue_ ( mriFunctionalDataRef this,
                             MRI*                 iData,
                             xVoxelRef            iIdx,
                             int                  inCondition,
                             int                  inTimePoint,
                             float*               oValue )
{

  int nFrame = 0;

  FunD_GetDataFrameNumber(inCondition,inTimePoint,&nFrame);

  MRIsampleVolumeFrameType( iData,
                            xVoxl_GetFloatX(iIdx),
                            xVoxl_GetFloatY(iIdx),
                            xVoxl_GetFloatZ(iIdx),
                            nFrame,
                            (int)this->mSampleType,
                            &this->mTmpReal );

  *oValue = this->mTmpReal;
}

void FunD_SetValue_ ( mriFunctionalDataRef this,
                      MRI*                 iData,
                      xVoxelRef            iIdx,
                      int                  inCondition,
                      int                  inTimePoint,
                      float                iValue )
{

  int nFrame = 0;

  FunD_GetDataFrameNumber(inCondition,inTimePoint,&nFrame);

  switch ( iData->type )
  {
  case MRI_UCHAR:
    MRIseq_vox( iData, xVoxl_GetX(iIdx),
                xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame) = iValue;
    break;
  case MRI_INT:
    MRIIseq_vox( iData, xVoxl_GetX(iIdx),
                 xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame) = iValue;
    break;
  case MRI_LONG:
    MRILseq_vox( iData, xVoxl_GetX(iIdx),
                 xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame) = iValue;
    break;
  case MRI_FLOAT:
    MRIFseq_vox( iData, xVoxl_GetX(iIdx),
                 xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame) = iValue;
    break;
  case MRI_SHORT:
    MRISseq_vox( iData, xVoxl_GetX(iIdx),
                 xVoxl_GetY(iIdx),xVoxl_GetZ(iIdx),nFrame) = iValue;
    break;
  default:
    break ;
  }
}

#endif /* FUND_USE_MACROS */

void FunD_GetSigma_ ( mriFunctionalDataRef this,
                      xVoxelRef            iFuncIdx,
                      int                  inTimePoint,
                      float*               oSigma )
{

  int   nFrame         = 0;


  FunD_GetSigmaFrameNumber(inTimePoint, &nFrame);

  switch ( this->mpData->type )
  {
  case MRI_UCHAR:
    *oSigma = MRIseq_vox( this->mpData, xVoxl_GetX(iFuncIdx),
                          xVoxl_GetY(iFuncIdx), xVoxl_GetZ(iFuncIdx), nFrame );
    break;
  case MRI_INT:
    *oSigma = MRIIseq_vox
              ( this->mpData, xVoxl_GetX(iFuncIdx),
                xVoxl_GetY(iFuncIdx), xVoxl_GetZ(iFuncIdx), nFrame );
    break;
  case MRI_LONG:
    *oSigma = MRILseq_vox
              ( this->mpData, xVoxl_GetX(iFuncIdx),
                xVoxl_GetY(iFuncIdx), xVoxl_GetZ(iFuncIdx), nFrame );
    break;
  case MRI_FLOAT:
    *oSigma = MRIFseq_vox
              ( this->mpData, xVoxl_GetX(iFuncIdx),
                xVoxl_GetY(iFuncIdx), xVoxl_GetZ(iFuncIdx), nFrame );
    break;
  case MRI_SHORT:
    *oSigma = MRISseq_vox
              ( this->mpData, xVoxl_GetX(iFuncIdx),
                xVoxl_GetY(iFuncIdx), xVoxl_GetZ(iFuncIdx), nFrame );
    break;
  default:
    *oSigma = 0;
    break ;
  }
}

FunD_tErr FunD_Verify ( mriFunctionalDataRef this )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  //  DebugEnterFunction( ("FunD_Verify( this=%p )", this) );

  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );
  DebugAssertThrowX( (FunD_kSignature == this->mSignature),
                     eResult, FunD_tErr_InvalidSignature );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  //  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_VerifyFuncIdx_ ( mriFunctionalDataRef this,
                                xVoxelRef            iFuncIdx )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugAssertThrowX( (NULL != this && NULL != iFuncIdx),
                     eResult, FunD_tErr_InvalidParameter );

  DebugAssertThrowX( (xVoxl_GetFloatX(iFuncIdx) >= 0 &&
                      xVoxl_GetFloatX(iFuncIdx) < this->mpData->width &&
                      xVoxl_GetFloatY(iFuncIdx) >= 0 &&
                      xVoxl_GetFloatY(iFuncIdx) < this->mpData->height &&
                      xVoxl_GetFloatZ(iFuncIdx) >= 0 &&
                      xVoxl_GetFloatZ(iFuncIdx) < this->mpData->depth),
                     eResult, FunD_tErr_InvalidFunctionalVoxel );

  DebugCatch;
  //  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  return eResult;
}

FunD_tErr FunD_VerifyTimePoint ( mriFunctionalDataRef this,
                                 int                  iTimePoint )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_VerifyTimePoint( this=%p, iTimePoint=%p )",
                       this, iTimePoint) );

  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Checking time point bounds") );
  DebugAssertThrowX( (iTimePoint >= 0 && iTimePoint < this->mNumTimePoints),
                     eResult, FunD_tErr_InvalidTimePoint );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}

FunD_tErr FunD_VerifyCondition ( mriFunctionalDataRef this,
                                 int                  iCondition )
{

  FunD_tErr eResult = FunD_tErr_NoError;

  DebugEnterFunction( ("FunD_VerifyTimePoint( this=%p, iCondition=%p )",
                       this, iCondition) );

  DebugAssertThrowX( (NULL != this), eResult, FunD_tErr_InvalidParameter );

  DebugNote( ("Checking condition bounds") );
  DebugAssertThrowX( (iCondition >= 0 && iCondition < this->mNumConditions),
                     eResult, FunD_tErr_InvalidConditionIndex );

  DebugCatch;
  DebugCatchError( eResult, FunD_tErr_NoError, FunD_GetErrorString );
  EndDebugCatch;

  DebugExitFunction;

  return eResult;
}
