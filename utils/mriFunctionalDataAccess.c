#include <stdlib.h>
#include <math.h>
#include "machine.h"
#include "matrix.h"
#include "stats.h"
#include "mriFunctionalDataAccess.h"
#include "xDebug.h"
#include "xVoxel.h"

#ifndef min
#define min(x,y) x<y?x:y
#endif

#ifndef max
#define max(x,y) x>y?x:y
#endif

#define KEEP_NULL_CONDITION 1

#define kFileName_StemHeaderSuffix                   "dat"
#define kFileName_BFileHeaderSuffix                  "hdr"
#define kFileName_BFileShortSuffix                   "bshort"
#define kFileName_BFileFloatSuffix                   "bfloat"
#define kFileName_Register                           "register.dat"


char FunD_ksaErrorString [FunD_tErr_knNumErrorCodes][256] = {
  "No error.",
  "Invalid ptr to volume (was probably NULL).",
  "Path not found.",
  "No stem provided and it couldn't be guessed.",
  "Data not found in specified path.",
  "Header file not found in specified path.",
  "Couldn't allocate volume (memory allocation failed).",
  "Unrecognized header format.",
  "Questionable header format (found expected types of values, but keywords were different).",
  "Couldn't find a recognizable data file.",
  "Couldn't allocate storage (memory allocation failed).",
  "Data has already been read.",
  "Slice file not found.",
  "Error reading slice data (unexpected EOF or expected data type mismatch).",
  "Couldn't allocate matrix (not my fault).",
  "Couldn't read register file.",
  "Couldn't calculate deviations.",
  "Error performing an operation on a transform.",
  "Invalid time resolution, number of time points must be evenly divisble by it.",
  "Invalid number of pre stim time points.",
  "Invalid functional voxel, out of bounds.",
  "Invalid condition (is it zero-indexed?)",
  "Invalid time point (is it zero-indexed?)",
  "Invalid error code."
};

// ===================================================================== VOLUME

FunD_tErr FunD_New ( mriFunctionalDataRef* outVolume,
         char*                 inPathName, 
         char*                 inStem,
         char*                 inHeaderStem,
         char*                 inRegistration ) {
  
  mriFunctionalDataRef this;
  FunD_tErr theError;

  /*
  DebugPrint "\nFunD_New: path=%s stem=%s", inPathName, inStem EndDebugPrint;
  if( inHeaderStem ) 
    DebugPrint " header=%s", inHeaderStem EndDebugPrint;
  if( inRegistration ) 
    DebugPrint " registration=%s", inRegistration EndDebugPrint;
  DebugPrint "\n" EndDebugPrint;
  */

  // check to see if the path exists.

  // check to see if the path with the volume data extension exists.

  // if we didn't get a stem...
  if ( NULL == inStem ) {

    // try and guess one.
    // theError = FunD_GuessStem ();
  }

  // allocate the structure.
  this = (mriFunctionalDataRef) malloc (sizeof(mriFunctionalData));
  if ( NULL == this ) {
    return FunD_tErr_CouldntAllocateVolume;
  }

  // save the path and stem.
  strcpy ( this->mPath, inPathName );
  strcpy ( this->mStem, inStem );

  // if we have a header stem, use that, else copy in the normal stem.
  if( inHeaderStem ) {
    strcpy ( this->mHeaderStem, inHeaderStem );
  } else {
    strcpy ( this->mHeaderStem, inStem );
  }

  // if we have a different registration path, use it, else use 
  // the normal path.
  if( inRegistration ) {
    strcpy ( this->mRegistration, inRegistration );
  } else {
    strcpy ( this->mRegistration, inPathName );
  }

  // initialize values.
  this->mData = NULL;
  this->mDeviations = NULL;
  this->mCovMtx = NULL;
  this->mSigma = 0;
  this->mIsErrorDataPresent = FALSE;
  this->mMaxValue = -100000;
  this->mMinValue = 100000;

  // init our transform objects.
  Trns_New( &(this->mIdxToIdxTransform) );
  Trns_New( &(this->mRASToIdxTransform) );
  
  // try parsing a stem header.
  theError = FunD_ParseStemHeader ( this );


  // if we couldn't find it...
  if ( FunD_tErr_HeaderNotFound == theError ) {

    // try the analyse.dat header
    theError = FunD_ParseAnalyseHeader ( this );
  }

  // still couldn't find it?
  if ( FunD_tErr_HeaderNotFound == theError ) {

    // try a bfile header.
    theError = FunD_ParseBFileHeader ( this );
  }

  // if we have any errors now, abort
  if ( FunD_tErr_NoError != theError ) {
    free ( this );
    return theError;
  }

  // find out what kind of data we have.
  theError = FunD_DetermineDataType ( this );
  if ( FunD_tErr_NoError != theError ) {
    free ( this );
    return theError;
  }

  // count the num of slices we have.
  this->mNumSlices = FunD_CountSliceFiles ( this );

  // calc our plane sizes
  FunD_CalcDataPlaneSizes ( this );

  // now parse the registration file.
  theError = FunD_ParseRegistrationAndInitMatricies ( this );
  if ( FunD_tErr_NoError != theError ) {
    free ( this );
    return theError;
  }

  // now parse the data.
  theError = FunD_ParseData ( this );
  if ( FunD_tErr_NoError != theError ) {
    free ( this );
    return theError;
  }

  // if we're using errors...
  if ( FunD_IsErrorDataPresent ( this ) ) {

    // calc the data.
    theError = FunD_CalcDeviations ( this );

    if ( FunD_tErr_NoError != theError ) {

      // if we couldn't, just set the flag to false.
      if ( FunD_tErr_CouldntCalculateDeviations == theError ) {
  this->mIsErrorDataPresent = FALSE;

      } else {

  // if it was another error, return it.
  free ( this );
  return theError;
      }
    }
  }

  // init our temp storage
  xVoxl_New( &this->mpAnaRAS );
  xVoxl_New( &this->mpFuncRAS );

  // return the structure.
  *outVolume = this;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_Delete ( mriFunctionalDataRef * ioVolume ) {

  FunD_tErr theErr;
  mriFunctionalDataRef this;

  this = *ioVolume;

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // check to see if the data is allocated. if so, delete the data.
  if ( NULL != this->mData )
    free ( this->mData );
  if ( NULL != this->mDeviations )
    free ( this->mDeviations );
  if ( NULL != this->mCovMtx )
    free ( this->mCovMtx );

  // delete our transformers.
  Trns_Delete( &(this->mIdxToIdxTransform) );
  Trns_Delete( &(this->mRASToIdxTransform) );

  // delete the structure.
  free ( this );
  
  // return a nil pointer.
  *ioVolume = NULL;
  
  return FunD_tErr_NoError;
}

FunD_tErr FunD_GuessStem ( char* inPathName, char* outStem ) {
  
  return FunD_tErr_CouldntGuessStem;
}


FunD_tErr FunD_ParseStemHeader ( mriFunctionalDataRef this ) {

  FunD_tErr eResult        = FunD_tErr_NoError;
  FILE*     pHeader        = NULL;
  char      sFileName[128] = "";
  tBoolean  bGood          = FALSE;
  char      sKeyword[128]  = "";
  int       nNumValues     = 0;
  int       nValue         = 0;
  int       nValuesRead    = 0;

  /* make the filename. if we have an alternate header stem */
  sprintf ( sFileName, "%s/%s.%s", this->mPath, this->mHeaderStem,
      kFileName_StemHeaderSuffix );
  
  /* try to open the file. */
  pHeader = fopen ( sFileName, "r" );
  if ( NULL == pHeader )
    return FunD_tErr_HeaderNotFound;

  /* start scanning for keywords... */
  while( !feof( pHeader ) ) {

    /* grab a keyword */
    nValuesRead = fscanf( pHeader, "%s", sKeyword );
    if( 1 != nValuesRead
  && !feof( pHeader ) ) {
      eResult = FunD_tErr_UnrecognizedHeaderFormat;
      goto error;
    }

    /* look at the keyword */
    if( strcmp( sKeyword, "TR" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mTimeResolution );
      bGood = (1 == nValuesRead);

    } else if( strcmp( sKeyword, "TPreStim" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mNumPreStimTimePoints );
      this->mNumPreStimTimePoints /= this->mTimeResolution;
      bGood = (1 == nValuesRead);

    } else if( strcmp( sKeyword, "nCond" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mNumConditions );
      bGood = (1 == nValuesRead);
#if KEEP_NULL_CONDITION
#else
      this->mNumConditions--;
#endif

    } else if( strcmp( sKeyword, "Nh" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mNumTimePoints );
      bGood = (1 == nValuesRead);

    } else if( strcmp( sKeyword, "Npercond" ) == 0 ) { 
#if KEEP_NULL_CONDITION
      nNumValues = this->mNumConditions;
#else
      nNumValues = this->mNumConditions + 1;
#endif
      for( nValue = 0; nValue < nNumValues; nValue++ )
  fscanf( pHeader, "%*d" );

    } else if( strcmp( sKeyword, "Rows" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mNumRows );
      bGood = (1 == nValuesRead);

    } else if( strcmp( sKeyword, "Cols" ) == 0 ) { 
      nValuesRead = fscanf( pHeader, "%d", &this->mNumCols );
      bGood = (1 == nValuesRead);

    } else if( strcmp( sKeyword, "SumXtX" ) == 0 ) { 
      nNumValues = pow( this->mNumTimePoints * (this->mNumConditions-1), 2 );
      for( nValue = 0; nValue < nNumValues; nValue++ )
  fscanf( pHeader, "%*d" );

    } else if( strcmp( sKeyword, "hCovMtx" ) == 0 ) { 
      nNumValues = pow( this->mNumTimePoints * (this->mNumConditions-1), 2 );
      this->mCovMtx = (float*) malloc ( sizeof(float) * nNumValues );
      for( nValue = 0; nValue < nNumValues && bGood; nValue++ ) {
  nValuesRead = fscanf( pHeader, "%f", &(this->mCovMtx[nValue]) );
  bGood = (1 == nValuesRead);
      }
    }

    if( !bGood ) {
      eResult = FunD_tErr_UnrecognizedHeaderFormat;
      goto error;
    }
  }

  /*
  DebugPrint "FunD_ParseStemHeader(): Successfully parsed %s\n", 
    sFileName EndDebugPrint;
  DebugPrint "FunD_ParseStemHeader(): cols = %d rows = %d num conds = %d num time points = %d\n", 
    this->mNumCols, this->mNumRows, this->mNumConditions, this->mNumTimePoints 
    EndDebugPrint;
  */

  // in this format we have error data.
  this->mIsErrorDataPresent = TRUE;

  goto cleanup;

 error:

  DebugPrint "Error %d in FunD_ParseStemHeader: %s. Last parsed keyword: %s\n", eResult, FunD_GetErrorString(eResult), sKeyword EndDebugPrint;

 cleanup:

  if( NULL != pHeader )
    fclose( pHeader );

  return eResult;
}

FunD_tErr FunD_ParseAnalyseHeader ( mriFunctionalDataRef this ) {

  return FunD_tErr_HeaderNotFound;
}

FunD_tErr FunD_ParseBFileHeader ( mriFunctionalDataRef this ) {

  FILE * theHeaderFile;
  char theFileName [128];
    
  // make the filename
  sprintf ( theFileName, "%s/%s_000.%s", this->mPath, this->mHeaderStem,
     kFileName_BFileHeaderSuffix );

  // try to open the file.
  theHeaderFile = fopen ( theFileName, "r" );
  if ( NULL == theHeaderFile )
    return FunD_tErr_HeaderNotFound;

  // we just have four numbers in this file. read them in.
  fscanf ( theHeaderFile, "%d %d %d %*f",
     &this->mNumCols, &this->mNumRows, &this->mNumTimePoints );

  // close the file, we're done.
  fclose ( theHeaderFile );

  // set the other values to signifiy that we have no conditions or
  // other time data.
  this->mNumConditions = 1;
  this->mNumPreStimTimePoints = 0;
  this->mTimeResolution = 1;
  this->mIsErrorDataPresent = FALSE;
  this->mCovMtx = NULL;
  
  /*
  DebugPrint "FunD_ParseBFileHeader(): Successfully parsed %s\n", 
    theFileName EndDebugPrint;
  DebugPrint "FunD_ParseBFileHeader(): cols = %d rows = %d num time points = %d\n", this->mNumCols, this->mNumRows, this->mNumTimePoints EndDebugPrint;
  */

  return FunD_tErr_NoError;
}

FunD_tErr FunD_ReadKeywordAndValue 
                        ( FILE * inFile, 
        char* inExpectedKeyword, char* inValueType, 
        int inNumValues, int inValueSize, char* inAssign ) {

  char theKeyword[256];
  char theValue[8];
  char isQuestionable = FALSE, isGood;
  int theValueIndex;

  // first see if we can read in a keyword.
  isGood = fscanf ( inFile, "%s", theKeyword );

  // if not, return failure.
  if ( !isGood ) 
    return FunD_tErr_UnrecognizedHeaderFormat;

  // if so, but it doesn't match, mark this as questionable.
  if ( strcmp ( theKeyword, inExpectedKeyword ) )
    isQuestionable = TRUE;

  // for each value they want us to read...
  for ( theValueIndex = 0; theValueIndex < inNumValues; theValueIndex++ ) {

    // read a value.
    isGood = fscanf ( inFile, inValueType, theValue );
    
    // if the conversion failed, return failure.
    if ( !isGood ) 
      return FunD_tErr_UnrecognizedHeaderFormat;
    
    // if they want us to assign it...
    if ( NULL != inAssign  ) {

      // copy it in. assume a contiguous array of values if more than one,
      // so advance the dest pointer the proper amount.
      memcpy ( inAssign+(theValueIndex*inValueSize), theValue, inValueSize );
    }
  }

  // if we got a questionable reading...
  if ( isQuestionable ) {

    // return warning.
    return FunD_tErr_QuestionableHeaderFormat;
  }
  
  return FunD_tErr_NoError;
}

FunD_tErr FunD_ParseRegistrationAndInitMatricies ( mriFunctionalDataRef this ) {

  char theFileName[256];
  fMRI_REG * theRegInfo;
  MATRIX* mTmp             = NULL;
  float   ps               = 0;
  float   st               = 0;
  float   slices           = 0;
  float   rows             = 0;
  float   cols             = 0;

  // read the registration info.
  sprintf ( theFileName, "%s/%s", this->mRegistration, kFileName_Register ); 
  theRegInfo = StatReadRegistration ( theFileName );
  if ( NULL == theRegInfo ) {
    DebugPrint "FunD_AllocateAndInitMatricies(): Couldn\'t read registration info from %s\n", theFileName EndDebugPrint;
    return FunD_tErr_CouldntReadRegisterFile;
  }

  // grab the info we need from it.
  this->mPixelSize = theRegInfo->in_plane_res;
  this->mSliceThickness = theRegInfo->slice_thickness;
  this->mBrightnessScale = theRegInfo->brightness_scale;
  strcpy ( this->mSubjectName, theRegInfo->name );

  // get our stats as floats
  ps     = this->mPixelSize;
  st     = this->mSliceThickness;
  slices = this->mNumSlices;
  rows   = this->mNumRows;
  cols   = this->mNumCols;
  
  // create the functional index to functional ras matrix
  mTmp = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixClear( mTmp );
  *MATRIX_RELT(mTmp,1,1) = -ps;
  *MATRIX_RELT(mTmp,2,3) = st;
  *MATRIX_RELT(mTmp,3,2) = -ps;
  *MATRIX_RELT(mTmp,1,4) = (ps*cols) / 2.0;
  *MATRIX_RELT(mTmp,2,4) = -(st*slices) / 2.0;
  *MATRIX_RELT(mTmp,3,4) = (ps*rows) / 2.0;
  *MATRIX_RELT(mTmp,4,4) = 1.0;
  Trns_CopyBtoRAS( this->mIdxToIdxTransform, mTmp );
  Trns_CopyBtoRAS( this->mRASToIdxTransform, mTmp );

  // create the anatomical index to anatomical ras matrix
  MatrixClear( mTmp );
  *MATRIX_RELT(mTmp,1,1) = -1.0;
  *MATRIX_RELT(mTmp,2,3) = 1.0;
  *MATRIX_RELT(mTmp,3,2) = -1.0;
  *MATRIX_RELT(mTmp,1,4) = 128.0;
  *MATRIX_RELT(mTmp,2,4) = -128.0;
  *MATRIX_RELT(mTmp,3,4) = 128.0;
  *MATRIX_RELT(mTmp,4,4) = 1.0;
  Trns_CopyAtoRAS( this->mIdxToIdxTransform, mTmp );

  /* the ras to idx transformer gets the identity matrix here */
  MatrixIdentity( 4, mTmp );
  Trns_CopyAtoRAS( this->mRASToIdxTransform, mTmp );

  /* a is anatomical, b is functional, mri2fmri takes us from a_ras
     to b_ras */
  Trns_CopyARAStoBRAS( this->mIdxToIdxTransform, theRegInfo->mri2fmri );
  Trns_CopyARAStoBRAS( this->mRASToIdxTransform, theRegInfo->mri2fmri );

  // get rid of the registration info.
  StatFreeRegistration ( &theRegInfo );

  MatrixFree( &mTmp );

  /*
  Trns_DebugPrint_( this->mIdxToIdxTransform );
  Trns_DebugPrint_( this->mRASToIdxTransform );
  */

  return FunD_tErr_NoError;
}

FunD_tErr FunD_SaveRegistration ( mriFunctionalDataRef this ) {

  fMRI_REG* regInfo;
  MATRIX*   mRegistration;
  char      sFileName[256];
  char      sBackupFileName[256];
  FILE*     pFile;
  int       nBackup = 1;

  /* make a reg info struct */
  regInfo = (fMRI_REG*) malloc( sizeof( fMRI_REG ));

  /* fill out the fmri registration struct */
  regInfo->in_plane_res = this->mPixelSize;
  regInfo->slice_thickness = this->mSliceThickness;
  regInfo->brightness_scale = this->mBrightnessScale;
  strcpy( regInfo->name, this->mSubjectName );

  /* allocate matrices and copy */
  Trns_GetARAStoBRAS( this->mIdxToIdxTransform, &mRegistration );
  regInfo->fmri2mri = MatrixCopy( mRegistration, NULL );
  regInfo->mri2fmri = MatrixInverse( regInfo->fmri2mri, NULL );

  /* if a registration already exists... */
  sprintf( sFileName, "%s/%s", this->mRegistration, kFileName_Register ); 
  pFile = fopen( sFileName, "r" );
  if( NULL != pFile ) {

    while( NULL != pFile ) {
      
      fclose( pFile );
      
      /* keep appending an increasing number to the end */
      sprintf( sBackupFileName, "%s.%d", sFileName, nBackup++ );
      pFile = fopen( sBackupFileName, "r" );
    }

    /* copy the registration file to backup */
  }

  /* write it to disk */
  StatWriteRegistration( regInfo, sFileName );

  StatFreeRegistration( &regInfo );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_SetRegistrationToIdentity ( mriFunctionalDataRef this ) {

  MATRIX* mID;

  mID = MatrixAlloc( 4, 4, MATRIX_REAL );
  MatrixIdentity( 4, mID );
  
  Trns_CopyARAStoBRAS( this->mIdxToIdxTransform, mID );
  Trns_CopyARAStoBRAS( this->mRASToIdxTransform, mID );

  MatrixFree( &mID );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_DetermineDataType ( mriFunctionalDataRef this ) {

  char theFileName [128];
  FILE * theFile;

  // make a short file name.
  this->mDataType = FunD_tDataType_Short;
  FunD_MakeSliceFileName ( this, 0, theFileName );
  theFile = fopen ( theFileName, "r" );
  if ( NULL != theFile ) {

    fclose ( theFile );
    return FunD_tErr_NoError;
  }

  // try a bfloat.
  this->mDataType = FunD_tDataType_Float;
  FunD_MakeSliceFileName ( this, 0, theFileName );
  theFile = fopen ( theFileName, "r" );
  if ( NULL != theFile ) {

    fclose ( theFile );
    return FunD_tErr_NoError;
  }

  // couldn't find anything
  return FunD_tErr_CouldntDetermineDataType;
}

FunD_tErr FunD_ParseData ( mriFunctionalDataRef this ) {

  int theNumDataValues, theDataSize, theNumValuesRead;
  int theNumValuesInPlane, thePlaneSize;
  int theSliceIndex, theTimePointIndex, theCondIndex, theRowIndex, theColIndex;
 xVoxelRef theFunctionalVoxel;
  FILE * theSliceFile;
  char theSliceFileName [256];
  char* theDestPtr;
  char* theValuePtr;
  float theFloatValue;
  short theShortValue;

  // getting compiler warnings because these weren't initialized...??
  theRowIndex = 0;
  theColIndex = 0;

  // make sure we haven't already read the data.
  if ( NULL != this->mData )
    return FunD_tErr_DataAlreadyRead;
 
  // size is rows * cols * slices * time points * conditions
  theNumDataValues = this->mNumRows * this->mNumCols * this->mNumSlices *
    this->mNumTimePoints * this->mNumConditions;

  // figure out how big our storage should be and then allocate it. also
  // point our value ptr to the right value variable.
  switch ( this->mDataType ) {
  case FunD_tDataType_Short:
    theDataSize = sizeof ( short );
    theValuePtr = (char*) &theShortValue;
    break;
  case FunD_tDataType_Float:
    theDataSize = sizeof ( float );
    theValuePtr = (char*) &theFloatValue;
    break;
  default:
    theDataSize = 1;
    theValuePtr = NULL;
  }
  this->mData = (char*) malloc ( theDataSize * theNumDataValues );

  // make sure it succeeded
  if ( NULL == this->mData )
    return FunD_tErr_CouldntAllocateStorage;

  // set the number of values.
  this->mNumDataValues = theNumDataValues;

  // start our dest ptr. 
  theDestPtr = this->mData;

  // create a voxel to store the location of the values we're reading in.
  xVoxl_New ( &theFunctionalVoxel );

  // calculate the num values in a plane.
  theNumValuesInPlane = this->mNumRows * this->mNumCols;
  thePlaneSize = theNumValuesInPlane * theDataSize;

  // for each slice...
  for ( theSliceIndex = 0; theSliceIndex < this->mNumSlices; theSliceIndex++ ){
    
    // create the file name and try to open it.
    theSliceFile = NULL;
    FunD_MakeSliceFileName ( this, theSliceIndex, theSliceFileName );

    // open the file.
    theSliceFile = fopen ( theSliceFileName, "r" );

    // if that doesn't work either, bail out.
    if ( NULL == theSliceFile ) {
      DebugPrint "\tCouldn't open slice file %s.\n",
  theSliceFileName EndDebugPrint;
      return FunD_tErr_SliceFileNotFound;
    }

#if KEEP_NULL_CONDITION
#else
    // if we have conditions, we want to skip the first null condition.
    if ( this->mNumConditions > 1 ) {
    
      // if we're reading in errors, we want to skip the null condition,
      // but read sigma from the second plane if this is slice 0.
      if ( theSliceIndex == 0 ) {
  
  // skip ahead one plane of data for each time point.
  if ( fseek ( theSliceFile,
         thePlaneSize * this->mNumTimePoints, SEEK_CUR ) ) {
    
    DebugPrint "\tAttempting to skip first plane before reading sigma, but couldn't read from file %s.\n",
      theSliceFileName EndDebugPrint;
    fclose ( theSliceFile );
    return FunD_tErr_ErrorReadingSliceData;
  }
  
  // now read a value
  fread ( theValuePtr, theDataSize, 1, theSliceFile );
  
  // order it and save it into sigma.
  switch ( this->mDataType ) {
  case FunD_tDataType_Short:
    theShortValue = orderShortBytes ( theShortValue ); 
    this->mSigma = theShortValue;
    break;
  case FunD_tDataType_Float:
    theFloatValue = orderFloatBytes ( theFloatValue ); 
    this->mSigma = theFloatValue;
    break;
  }

  // and skip the rest of the time points.
  if ( fseek ( theSliceFile,
         (this->mNumTimePoints * thePlaneSize) - theDataSize,
         SEEK_CUR ) ) {
    DebugPrint "\tAttempting to skip a plane after reading sigma, but couldn't read from file %s.\n", theSliceFileName EndDebugPrint;
    fclose ( theSliceFile );
    return FunD_tErr_ErrorReadingSliceData;
  }

      } else {
  
  // not slice 0, just skip two planes of time poitns.
  if ( fseek ( theSliceFile,
         this->mNumTimePoints * thePlaneSize * 2, SEEK_CUR ) ) {
    DebugPrint "\tAttempting to skip null condition planes, but couldn't read from file %s.\n", theSliceFileName EndDebugPrint;
    fclose ( theSliceFile );
    return FunD_tErr_ErrorReadingSliceData;
  }
      }

    }
#endif

    // for each condtion...
    for ( theCondIndex = 0; 
    theCondIndex < this->mNumConditions; 
    theCondIndex++) {

      for ( theTimePointIndex = 0;
      theTimePointIndex < this->mNumTimePoints;
      theTimePointIndex++ ) {

  for ( theRowIndex = 0;
        theRowIndex < this->mNumRows;
        theRowIndex++ ) {

    for ( theColIndex = 0; 
    theColIndex < this->mNumCols;
    theColIndex++ ) {
    
      // set the voxel to the current coords
      xVoxl_Set ( theFunctionalVoxel,
      theColIndex, theRowIndex, theSliceIndex );
      
      // read in a float.
      theNumValuesRead = fread ( theValuePtr, theDataSize,
               1, theSliceFile );
      
      // make sure we got it.
      if ( theNumValuesRead != 1 ) {
        DebugPrint "\tcond %d time %d row %d col %d slice %d: Trying to read a value but couldn't read from file %s.\n",
    theCondIndex, theTimePointIndex, theRowIndex, theColIndex,
    theSliceIndex, theSliceFileName EndDebugPrint;
        fclose ( theSliceFile );
        return FunD_tErr_ErrorReadingSliceData;
      }
      
      // order it and save it.
      switch ( this->mDataType ) {
      case FunD_tDataType_Short:
        theShortValue = orderShortBytes ( theShortValue ); 
        FunD_SetValue ( this, theFunctionalVoxel,
            theCondIndex, theTimePointIndex, 
            theShortValue );
        break;
      case FunD_tDataType_Float:
        theFloatValue = orderFloatBytes ( theFloatValue ); 
        FunD_SetValue ( this, theFunctionalVoxel,
            theCondIndex, theTimePointIndex,
            theFloatValue );
        break;
      }
    }
  }
      }
   
      // if we're reading errors...
      if ( FunD_IsErrorDataPresent ( this ) ) {
  
  // if we're reading in errors and this is the first slice,
  // read sigma for the error data.
  if ( theSliceIndex == 0 ) {
  
    // read a value
    fread ( theValuePtr, theDataSize, 1, theSliceFile );
    
    // order it and save it into sigma.
    switch ( this->mDataType ) {
    case FunD_tDataType_Short:
      theShortValue = orderShortBytes ( theShortValue ); 
      this->mSigma = theShortValue;
      break;
    case FunD_tDataType_Float:
      theFloatValue = orderFloatBytes ( theFloatValue ); 
      this->mSigma = theFloatValue;
      break;
  }

    // and skip the rest of the time points.
    if ( fseek ( theSliceFile,
           (this->mNumTimePoints * thePlaneSize) - theDataSize,
           SEEK_CUR ) ) {
      
      DebugPrint "\tATtempting to skip a plane after reading sigma, but couldn't read from file %s.\n",
        theSliceFileName EndDebugPrint;
      fclose ( theSliceFile );
      return FunD_tErr_ErrorReadingSliceData;
    }

  } else {
    
    // we skip one plane of data for every time point for
    // error values.
    if ( fseek ( theSliceFile, 
           thePlaneSize * this->mNumTimePoints, SEEK_CUR ) ) {
      
      DebugPrint "\tcond %d time %d row %d col %d slice %d: Trying skip errors but couldn't read from file %s.\n",
        theCondIndex, theTimePointIndex, theRowIndex, theColIndex,
        theSliceIndex, theSliceFileName EndDebugPrint;
      fclose ( theSliceFile );
      return FunD_tErr_ErrorReadingSliceData;
    }
  }
      }
    }

    // we should be at the end of the file. try to read one more char and make
    // sure we got an eof flag.
    fgetc ( theSliceFile );
    if ( !feof ( theSliceFile ) ) {

      DebugPrint "\t!!!!!!!!!! Not at end of file!!!\n" EndDebugPrint;

      // if we wern't, close the file and return an error.
      //fclose ( theSliceFile );
      //return FunD_tErr_ErrorReadingSliceData;

    }

    // close this slice.
    fclose ( theSliceFile );

  }

  xVoxl_Delete ( &theFunctionalVoxel );

  return FunD_tErr_NoError;

}

FunD_tErr FunD_CalcDeviations ( mriFunctionalDataRef this ) {

  int theErrorSize;
  int theErrorIndex, theCovMtxIndex;
  int theNumConditions, theNumTimePoints;
  int theCondIndex, theTimePointIndex;

  if ( NULL == this->mCovMtx ) {
    DebugPrint "FunD_CalcDeviations(): mCovMtx was null.\n" EndDebugPrint;
    return FunD_tErr_CouldntCalculateDeviations;
  }

  // number of conditions and number of time points.
  theNumConditions = (this->mNumConditions);
  theNumTimePoints = (this->mNumTimePoints);

  // allocate the deviations.
  theErrorSize = theNumConditions * theNumTimePoints;
  this->mDeviations = (float*) malloc (theErrorSize * sizeof(float));
  bzero( (this->mDeviations), theErrorSize * sizeof(float) );
  if ( NULL == this->mDeviations ) {
    DebugPrint "FunD_CalcDeviations(): Allocation of deviations storage of size %d failed.\n", theErrorSize EndDebugPrint;
    return FunD_tErr_CouldntAllocateStorage;
  }

  // start at condition one because we don't want any error data for
  // the null condition (0).
  for ( theCondIndex = 1; theCondIndex < theNumConditions; theCondIndex++ ) {
    for ( theTimePointIndex = 0; 
    theTimePointIndex < theNumTimePoints;
    theTimePointIndex++ ) {

      // use the full number of conditions when calcing the index here.
      theErrorIndex = (theNumTimePoints * theCondIndex) + theTimePointIndex;

      // since the cov mtx is based on non-null conditions, use the cond
      // index - 1 here to calc the cov mtx index.
      theCovMtxIndex = 
  ((theNumTimePoints * (theCondIndex-1) + theTimePointIndex) * 
   (theNumConditions-1) * theNumTimePoints) +
  (theNumTimePoints * (theCondIndex-1) + theTimePointIndex);

      (this->mDeviations)[theErrorIndex] = (float)this->mSigma *
  sqrt ( (this->mCovMtx)[theCovMtxIndex] );
    }
  }

  return FunD_tErr_NoError;
}

FunD_tErr FunD_ApplyTransformToRegistration( mriFunctionalDataRef this,
               MATRIX*              iTransform ){

  FunD_tErr eResult    = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;
  MATRIX*   invTransform = NULL;

  // make sure we're valid.
  eResult = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != eResult )
    return eResult;

  /* we're getting the matrix in ras space. */
  invTransform = MatrixInverse( iTransform, NULL );

  /* apply the transform to both transformations */
  eTransform = Trns_ApplyTransform( this->mIdxToIdxTransform, invTransform );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;
  eTransform = Trns_ApplyTransform( this->mRASToIdxTransform, invTransform );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;

  MatrixFree( &invTransform );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_TranslateRegistration ( mriFunctionalDataRef this,
               float                ifDistance,
               tAxis                iAxis ) {
  
  FunD_tErr eResult    = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;

  // make sure we're valid.
  eResult = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != eResult )
    return eResult;

  /* do the inverse of this action */
  eTransform = Trns_Translate( this->mIdxToIdxTransform, -ifDistance, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;
  eTransform = Trns_Translate( this->mRASToIdxTransform, -ifDistance, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_RotateRegistration ( mriFunctionalDataRef this,
            float                ifDegrees,
            tAxis                iAxis ) {

  FunD_tErr eResult    = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;

  // make sure we're valid.
  eResult = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != eResult )
    return eResult;

  /* do the inverse of this action */
  eTransform = Trns_Rotate( this->mIdxToIdxTransform, -ifDegrees, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;
  eTransform = Trns_Rotate( this->mRASToIdxTransform, -ifDegrees, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_ScaleRegistration ( mriFunctionalDataRef this,
           float                ifFactor,
           tAxis                iAxis ) {

  FunD_tErr eResult    = FunD_tErr_NoError;
  Trns_tErr eTransform = Trns_tErr_NoErr;

  // make sure we're valid.
  eResult = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != eResult )
    return eResult;

  /* do the inverse of this action */
  eTransform = Trns_Scale( this->mIdxToIdxTransform, 1.0 / ifFactor, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;
  eTransform = Trns_Scale( this->mRASToIdxTransform, 1.0 / ifFactor, iAxis );
  if( Trns_tErr_NoErr != eTransform )
    return FunD_tErr_ErrorAccessingTransform;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetDataAtAnaIdx ( mriFunctionalDataRef this,
         xVoxelRef inVoxel, 
         int inCondition, 
         int inTimePoint,
         float *outData ) {

  FunD_tErr theErr             = FunD_tErr_NoError;
  xVoxelRef theFunctionalVoxel = NULL;
  float     theValue           = 0;

  xVoxl_New ( &theFunctionalVoxel );

  *outData = 0;

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr ) {
    goto cleanup;
  }

  // convert to func index
  FunD_ConvertAnaIdxToFuncIdx( this, inVoxel, theFunctionalVoxel );

  // make sure all the indicies are valid.
  if ( !FunD_IsFunctionalVoxelValid(this,theFunctionalVoxel) ) {
    theErr =  FunD_tErr_InvalidFunctionalVoxel;
    goto cleanup;
  }
  
  if ( !FunD_IsConditionIndexValid(this,inCondition) ) {
    theErr = FunD_tErr_InvalidConditionIndex;
    goto cleanup;
  }
  
  if ( !FunD_IsTimePointValid(this,inTimePoint) ) {
    theErr = FunD_tErr_InvalidTimePoint;
    goto cleanup;
  }  

  // get the data at this point.
  theValue = FunD_GetValue ( this, theFunctionalVoxel, 
           inCondition, inTimePoint );

  // return it.
  *outData = theValue;
  
 cleanup:
  
  xVoxl_Delete ( &theFunctionalVoxel );

  return theErr;
}

FunD_tErr FunD_GetDataAtRAS ( mriFunctionalDataRef this,
            xVoxelRef            inVoxel, 
            int                  inCondition, 
            int                  inTimePoint,
            float*               outData ) {

  FunD_tErr theErr             = FunD_tErr_NoError;
 xVoxelRef         theFunctionalVoxel = NULL;
  float theValue;

  xVoxl_New ( &theFunctionalVoxel );

  *outData = 0;

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr ) {
    goto cleanup;
  }

  // convert to func index
  FunD_ConvertRASToFuncIdx( this, inVoxel, theFunctionalVoxel );

  // make sure all the indicies are valid.
  if ( !FunD_IsFunctionalVoxelValid(this,theFunctionalVoxel) ) {
    theErr =  FunD_tErr_InvalidFunctionalVoxel;
    goto cleanup;
  }
  
  if ( !FunD_IsConditionIndexValid(this,inCondition) ) {
    theErr = FunD_tErr_InvalidConditionIndex;
    goto cleanup;
  }
  
  if ( !FunD_IsTimePointValid(this,inTimePoint) ) {
    theErr = FunD_tErr_InvalidTimePoint;
    goto cleanup;
  }  

  // get the data at this point.
  theValue = FunD_GetValue ( this, theFunctionalVoxel, 
             inCondition, inTimePoint );

  // return it.
  *outData = theValue;
  
  cleanup:
  
  xVoxl_Delete ( &theFunctionalVoxel );

  return theErr;
}

FunD_tErr FunD_GetDataAtAnaIdxForAllTimePoints
                        ( mriFunctionalDataRef this,xVoxelRef inVoxel, int inCondition, 
        float *outData ) {

  FunD_tErr theErr;
 xVoxelRef theFunctionalVoxel;
  int theTimePoint;
  float theValue;

  xVoxl_New( &theFunctionalVoxel );

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // convert to func index
  FunD_ConvertAnaIdxToFuncIdx( this, inVoxel, theFunctionalVoxel );

  // make sure all the indicies are valid.
  if ( !FunD_IsFunctionalVoxelValid(this,theFunctionalVoxel) ) {
    theErr = FunD_tErr_InvalidFunctionalVoxel;
    goto cleanup;
  }

  if ( !FunD_IsConditionIndexValid(this,inCondition) ) {
    theErr =  FunD_tErr_InvalidConditionIndex;
    goto cleanup;
  }
  
  // for each time point...
  for ( theTimePoint = 0;
  theTimePoint < this->mNumTimePoints;
  theTimePoint++ ) {
    
    // get the data at this point.
    theValue = FunD_GetValue ( this, theFunctionalVoxel,
         inCondition, theTimePoint );

    // put it in output array.
    outData[theTimePoint] = theValue;
  }
  
      cleanup:
  
  xVoxl_Delete ( &theFunctionalVoxel );

  return theErr;
}

FunD_tErr FunD_GetDeviation ( mriFunctionalDataRef this,
               int inCondition, int inTimePoint,
               float * outValue ) {

  FunD_tErr theErr;

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *outValue = (this->mDeviations)[(inCondition*(this->mNumTimePoints)) +
         inTimePoint];

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetDeviationForAllTimePoints
                        ( mriFunctionalDataRef this, int inCondition, 
        float *outData ) {

  FunD_tErr theErr;
  int theTimePoint;
  float theValue;

   // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // make sure all the indicies are valid.
  if ( !FunD_IsConditionIndexValid(this,inCondition) ) {
    return FunD_tErr_InvalidConditionIndex;
  }
  
  // for each time point...
  for ( theTimePoint = 0;
  theTimePoint < this->mNumTimePoints;
  theTimePoint++ ) {
    
    // get the deviation at this point.
    theErr = FunD_GetDeviation ( this, inCondition, theTimePoint,
             &theValue );
    if ( FunD_tErr_NoError != theErr )
      return theErr;

    // put it in output array.
    outData[theTimePoint] = theValue;
  }
  
  return FunD_tErr_NoError;
}


FunD_tErr FunD_ConvertTimePointToSecond ( mriFunctionalDataRef this,
               int inTimePoint,
               int* outSecond ) {

  int theTimeResolution, theFirstTimePoint;
  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // calc the time second
  theTimeResolution = this->mTimeResolution;
  theFirstTimePoint = -(this->mNumPreStimTimePoints
      * theTimeResolution); 
  *outSecond = theFirstTimePoint + (inTimePoint* theTimeResolution);
 

  return FunD_tErr_NoError;
}

FunD_tErr FunD_ConvertSecondToTimePoint ( mriFunctionalDataRef this,
               int inSecond,
               int* outTimePoint ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // calc the time point
  *outTimePoint = (inSecond / this->mTimeResolution) + 
    this->mNumPreStimTimePoints;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_SetTimeResolution ( mriFunctionalDataRef this, int inTimeRes ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // make sure the value is in bounds
  if ( !FunD_IsTimeResolutionValid ( this, inTimeRes ) )
    return FunD_tErr_InvalidTimeResolution;

  this->mTimeResolution = inTimeRes;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_SetNumPreStimTimePoints ( mriFunctionalDataRef this,
              int inNumPoints ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  // make sure the value is in bounds
  if ( !FunD_IsNumPreStimTimePointsValid ( this, inNumPoints ) )
    return FunD_tErr_InvalidNumPreStimTimePoints;

  this->mNumPreStimTimePoints = inNumPoints;

  return FunD_tErr_NoError;
}



FunD_tErr FunD_GetPath ( mriFunctionalDataRef this, char* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  strcpy ( out, this->mPath );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetStem ( mriFunctionalDataRef this, char* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  strcpy ( out, this->mStem );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetSubjectName ( mriFunctionalDataRef this, char* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  strcpy ( out, this->mSubjectName );

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumCols ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumCols;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumRows ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumRows;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumSlices ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumSlices;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumTimePoints ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumTimePoints;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetFirstCondition ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = 0;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetLastCondition ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumConditions - 1;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumConditions ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumConditions;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetTimeResolution ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mTimeResolution;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumPreStimTimePoints ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumPreStimTimePoints;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetNumDataValues ( mriFunctionalDataRef this, int* out ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *out = this->mNumDataValues;

  return FunD_tErr_NoError;
}

FunD_tErr FunD_GetValueRange ( mriFunctionalDataRef this,
             float*               outMin,
             float*               outMax ) {

  FunD_tErr theErr;
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  *outMin = this->mMinValue;
  *outMax = this->mMaxValue;

  return FunD_tErr_NoError;
}


FunD_tErr FunD_DebugPrint ( mriFunctionalDataRef this ) {

  FunD_tErr theErr;

  // int theMatrixSize, i;

  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    return theErr;

  DebugPrint "Volume:\n" EndDebugPrint;
  DebugPrint "\tPath and stem: %s/%s\n", this->mPath,this->mStem EndDebugPrint;
  DebugPrint "\tSubjectName: %s\n", this->mSubjectName EndDebugPrint;
  DebugPrint "\tCols: %d\n", this->mNumCols EndDebugPrint;
  DebugPrint "\tRows: %d\n", this->mNumRows EndDebugPrint;
  DebugPrint "\tSlices: %d\n", this->mNumSlices EndDebugPrint;
  DebugPrint "\tPixelSize: %f\n", this->mPixelSize EndDebugPrint;
  DebugPrint "\tSliceThickness: %f\n", this->mSliceThickness EndDebugPrint;
  DebugPrint "\tNumTimePoints: %d\n", this->mNumTimePoints EndDebugPrint;
  DebugPrint "\tNumConditions: %d\n", this->mNumConditions EndDebugPrint;
  DebugPrint "\tTimeResoultion: %d\n", this->mTimeResolution EndDebugPrint;
  DebugPrint "\tNumPreStimTimePoints: %d\n", this->mNumPreStimTimePoints EndDebugPrint;
  DebugPrint "\tNumDataValues: %d\n", this->mNumDataValues EndDebugPrint;

  DebugPrint "\tdata plane sizes: %d %d %d %d %d\n",
    this->mDataPlaneSizes[0], this->mDataPlaneSizes[1],
    this->mDataPlaneSizes[2], this->mDataPlaneSizes[3],
    this->mDataPlaneSizes[4] EndDebugPrint;

  /*
  DebugPrint "Functional to Anatomical:\n" EndDebugPrint;
  MatrixPrint ( stderr, this->mFunctionalToAnatomicalMatrix );
  DebugPrint "Anatomical to Functional:\n" EndDebugPrint;
  MatrixPrint ( stderr, this->mAnatomicalToFunctionalMatrix );
  */

  DebugPrint "\tSigma:%2.5f\n", this->mSigma EndDebugPrint;

  /*
  if ( NULL != this->mCovMtx ) {
    DebugPrint "\tCovMtx: \n" EndDebugPrint;
    theMatrixSize = ( this->mNumTimePoints * this->mNumConditions ) *
                    ( this->mNumTimePoints * this->mNumConditions );
    
    for ( i = 0; i < theMatrixSize; i++ ) {
      if ( (i % (this->mNumTimePoints*this->mNumConditions)) == 0 ) {
  DebugPrint "\n\t\t" EndDebugPrint;
      }
      DebugPrint "%2.2f ", this->mCovMtx[i] EndDebugPrint;
    } 
  } 
  */

  DebugPrint "\n" EndDebugPrint;

  return FunD_tErr_NoError;
}

char* FunD_GetErrorString ( FunD_tErr inErr ) {

  if ( !(inErr >= 0 && inErr < FunD_tErr_knNumErrorCodes) )
    inErr = FunD_tErr_InvalidErrorCode;

  return (char*)(FunD_ksaErrorString[inErr]);
}

tBoolean FunD_IsErrorDataPresent ( mriFunctionalDataRef this ) {

  return this->mIsErrorDataPresent;
}


#if 0

/* FOR REFERENCE ====================================================== */

void FunD_ConvertFunctionalRASToFunctionalIndex_ ( mriFunctionalDataRef this,
               xVoxelRef inFuncRAS,
               xVoxelRef outFuncIndex ) {

  int              theVoxX            = 0;
  int              theVoxY            = 0;
  int              theVoxZ            = 0;
  float            x                  = 0;
  float            y                  = 0;
  float            z                  = 0;
  float            ps                 = 0;
  float            st                 = 0;
  float            slices             = 0;
  float            rows               = 0;
  float            cols               = 0;

  // convert the functional RAS to a functional index voxel, rounding.
  x      = xVoxl_GetFloatX(inFuncRAS);
  y      = xVoxl_GetFloatY(inFuncRAS);
  z      = xVoxl_GetFloatZ(inFuncRAS);
  ps     = this->mPixelSize;
  st     = this->mSliceThickness;
  slices = this->mNumSlices;
  rows   = this->mNumRows;
  cols   = this->mNumCols;
  theVoxX =   nint( (     (ps  * cols   / 2.0) - x) / ps );
  theVoxY = rows-1- nint((rows-1.0) -( (ps  * rows   / 2.0) - z ) / ps );
  theVoxZ =   nint( ( y - (-st * slices / 2.0)   )  / st );

  xVoxl_Set( outFuncIndex, theVoxX, theVoxY, theVoxZ );
}

void FunD_ConvertFunctionalIndexToFunctionalRAS_ ( mriFunctionalDataRef this,
               xVoxelRef inFuncIndex,
               xVoxelRef outFuncRAS ) {

  /* THIS FUNCTION IS NOT TESTED!!!!!!!!!!!!!!!!!!!! */
  float            theRASX            = 0;
  float            theRASY            = 0;
  float            theRASZ            = 0;
  float            x                  = 0;
  float            y                  = 0;
  float            z                  = 0;
  float            ps                 = 0;
  float            st                 = 0;
  float            slices             = 0;
  float            rows               = 0;
  float            cols               = 0;

  // convert the functional RAS to a functional index voxel, rounding.
  x      = xVoxl_GetFloatX(inFuncIndex);
  y      = xVoxl_GetFloatY(inFuncIndex);
  z      = xVoxl_GetFloatZ(inFuncIndex);
  ps     = this->mPixelSize;
  st     = this->mSliceThickness;
  slices = this->mNumSlices;
  rows   = this->mNumRows;
  cols   = this->mNumCols;
  theRASX = ((ps * cols) / 2.0) - (ps * x);
  theRASY = ((-st * slices) / 2.0) + (st * z);
  theRASZ = ((ps * rows) / 2.0) - (ps * y);

  xVoxl_SetFloat( outFuncRAS, theRASX, theRASY, theRASZ );
}
/* ==================================================================== */

#endif

static xVoxel sCoord1;
static xVoxel sCoord2;

void FunD_ConvertAnaIdxToFuncIdx ( mriFunctionalDataRef this,
            xVoxelRef inAnaIdx,
            xVoxelRef outFuncIdx ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  xVoxl_SetFloat( &sCoord1, xVoxl_ExpandFloat(inAnaIdx) );

  eTransform = Trns_ConvertAtoB( this->mIdxToIdxTransform, 
         &sCoord1, &sCoord2 );
  if( Trns_tErr_NoErr != eTransform ) {
    xVoxl_Set( outFuncIdx, -1, -1, -1 );
    return;
  }

  xVoxl_SetFloat( outFuncIdx, xVoxl_ExpandFloat(&sCoord2) );
}

void FunD_ConvertFuncIdxToAnaIdx ( mriFunctionalDataRef this,
            xVoxelRef inFuncIdx,
            xVoxelRef outAnaIdx ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  xVoxl_SetFloat( &sCoord1, xVoxl_ExpandFloat(inFuncIdx) );

  eTransform = Trns_ConvertBtoA( this->mIdxToIdxTransform, 
         &sCoord1, &sCoord2 );
  if( Trns_tErr_NoErr != eTransform ) {
    xVoxl_Set( outAnaIdx, -1, -1, -1 );
    return;
  }

  xVoxl_SetFloat( outAnaIdx, xVoxl_ExpandFloat(&sCoord2) );
}

void FunD_ConvertRASToFuncIdx ( mriFunctionalDataRef this,
         xVoxelRef inRAS,
         xVoxelRef outFuncIdx ) {

  Trns_tErr eTransform = Trns_tErr_NoErr;

  xVoxl_SetFloat( &sCoord1, xVoxl_ExpandFloat(inRAS) );

  eTransform = Trns_ConvertAtoB( this->mRASToIdxTransform, 
         &sCoord1, &sCoord2 );
  if( Trns_tErr_NoErr != eTransform ) {
    xVoxl_Set( outFuncIdx, -1, -1, -1 );
    return;
  }

  xVoxl_SetFloat( outFuncIdx, xVoxl_ExpandFloat(&sCoord2) );
}

FunD_tErr FunD_GetBoundsInAnatomical ( mriFunctionalDataRef this, 
          xVoxelRef  outMin,
          xVoxelRef  outMax ) {
  FunD_tErr theErr         = FunD_tErr_NoError;
 xVoxelRef         theVox         = NULL;
 xVoxelRef         theMin         = NULL;
 xVoxelRef         theMax         = NULL;
  int              nX             = 0;
  int              nY             = 0;
  int              nZ             = 0;

  xVoxl_New( &theVox );
  xVoxl_New( &theMin );
  xVoxl_New( &theMax );
  
  // make sure we're valid.
  theErr = FunD_AssertIsValid ( this );
  if ( FunD_tErr_NoError != theErr )
    goto error;
  
  /* set min and max to extremes */
  xVoxl_Set( theMin, 500, 500, 500 );
  xVoxl_Set( theMax, -500, -500, -500 );

  /* we want to check all combinations of the extremes in our indices */
  for( nZ = -10; nZ <= this->mNumSlices; nZ += this->mNumSlices+10 ) {
    for( nY = -10; nY <= this->mNumRows; nY += this->mNumRows+10 ) {
      for( nX = -10; nX <= this->mNumCols; nX += this->mNumCols+10 ) {

  /* set func index vox. */
  xVoxl_Set( theVox, nX, nY, nZ );
  
  /* convert to anatomical index. */
  FunD_ConvertFuncIdxToAnaIdx( this, theVox, theVox );

  /* set the lesser of each in the min voxel. */
  xVoxl_Set( theMin, 
       min( xVoxl_GetX(theMin), xVoxl_GetX(theVox) ),
       min( xVoxl_GetY(theMin), xVoxl_GetY(theVox) ),
       min( xVoxl_GetZ(theMin), xVoxl_GetZ(theVox) ) );

  /* set the greater of each in the max voxel. */
  xVoxl_Set( theMax, 
       max( xVoxl_GetX(theMax), xVoxl_GetX(theVox) ),
       max( xVoxl_GetY(theMax), xVoxl_GetY(theVox) ),
       max( xVoxl_GetZ(theMax), xVoxl_GetZ(theVox) ) );
      }
    }
  }

  /* return the voxels */
  xVoxl_Copy( outMin, theMin );
  xVoxl_Copy( outMax, theMax );

  goto cleanup;
    
 error:
  
  if( FunD_tErr_NoError != theErr ) {
    DebugPrint "Error in FunD_GetBoundsInAnatomicalRAS (%d): %s\n",
      theErr, FunD_GetErrorString(theErr) EndDebugPrint;
  }

 cleanup:

  xVoxl_Delete( &theVox );
  xVoxl_Delete( &theMin );
  xVoxl_Delete( &theMax );

  return theErr;

}


FunD_tErr FunD_AssertIsValid ( mriFunctionalDataRef this ) {

  // check for null pointer.
  if ( NULL == this )
    return FunD_tErr_InvalidPtr;

  // check for sensible values.

  // check if data is allocated.

  return FunD_tErr_NoError;
}

int FunD_CountSliceFiles ( mriFunctionalDataRef this ) {

  int theSliceNumber = 0;
  char theSliceFileName [256];
  FILE * theTestFile;
  char isDone = FALSE;

  // start looking for slice 0
  this->mNumSlices = 0;

  // while we haven't not found a slice...
  while ( !isDone ) {

    // make a name for this slice.
    FunD_MakeSliceFileName ( this, theSliceNumber, theSliceFileName );

    // try to open it.
    theTestFile = fopen ( theSliceFileName, "r" );

    // if we made it...
    if ( NULL != theTestFile ) {

      // close the file.
      fclose( theTestFile );

      // increment the slice count.
      this->mNumSlices++;
      theSliceNumber++;

    } else {

      // we're done.
      isDone = TRUE;
    }
  }

  return this->mNumSlices;
}

void FunD_CalcDataPlaneSizes ( mriFunctionalDataRef this ) {

  if ( this->mNumConditions > 1 ) {
    this->mDataPlaneSizes[4] = 
     this->mNumTimePoints * this->mNumSlices * this->mNumCols * this->mNumRows;
  } else {
    this->mDataPlaneSizes[4] = 0;
  }

  this->mDataPlaneSizes[0] = 1;
  this->mDataPlaneSizes[1] = this->mNumCols;
  this->mDataPlaneSizes[2] = this->mNumCols * this->mNumRows;
  this->mDataPlaneSizes[3] = this->mNumSlices * this->mNumCols * this->mNumRows;

}

inline 
int FunD_CoordsToIndex ( mriFunctionalDataRef this,
        xVoxelRef inFunctionalVoxel,
         int inConditionIndex, int inTimePoint ) {

  // the i,j,k coords are zero based but the condition and time point index
  // are one based, so we subtract one from them.
  return ( ( inConditionIndex * this->mDataPlaneSizes[4] ) +
     ( inTimePoint*  this->mDataPlaneSizes[3] ) +
     ( xVoxl_GetK(inFunctionalVoxel) * this->mDataPlaneSizes[2] ) +
     ( xVoxl_GetJ(inFunctionalVoxel) * this->mDataPlaneSizes[1] ) +
     ( xVoxl_GetI(inFunctionalVoxel) * this->mDataPlaneSizes[0] ) ) ;
}

inline 
void FunD_IndexToCoords ( mriFunctionalDataRef this, int inIndex,
         xVoxelRef outFunctionalVoxel,
          int* outConditionIndex, int* outTimePoint ) {
  
  int theI, theJ, theK;

  /*
    explanation: it's easier to think of this as taking a 5 digit number
    and stripping out the tenthousandsths digit, the thousandsths digit, 
    etc. you can do this with modulus and regular division. see below. the
    running totals on each division are in parenthasees.
    
    sum = (10000 * tenthousands) + (1000 * thousands) + 
          (100 * hundreds) + (10 * tens) + ones = 12345;
    ones = sum % 10000 (2345) % 1000 (345) % 100 (45) % 10 (5) / 1 (5);
    tens = sum % 10000 (2345) % 1000 (345) % 100 (45) / 10 (4);
    hundreds = sum % 10000 (2345) % 1000 (345) / 100 (3);
    thousndas = sum / 1000 (1);

    so we're doing the same kind of thing here, except instead of 10000, 1000
    and so one we're using the values we calculated earlier as average data
    index offsets, which are basically the number of values in each 'plane'
    of data.
  */
  theI = inIndex % 
    this->mDataPlaneSizes[4] %
    this->mDataPlaneSizes[3] %
    this->mDataPlaneSizes[2] %
    this->mDataPlaneSizes[1] /
    this->mDataPlaneSizes[0];
  
  theJ = inIndex % 
    this->mDataPlaneSizes[4] %
    this->mDataPlaneSizes[3] %
    this->mDataPlaneSizes[2] /
    this->mDataPlaneSizes[1];
  
  theK = inIndex % 
    this->mDataPlaneSizes[4] %
    this->mDataPlaneSizes[3] /
    this->mDataPlaneSizes[2];
  
  *outTimePoint = inIndex % 
    this->mDataPlaneSizes[4] /
    this->mDataPlaneSizes[3];
  
  *outConditionIndex = inIndex / 
    this->mDataPlaneSizes[4];

  xVoxl_Set ( outFunctionalVoxel, theI, theJ, theK );
}

void FunD_MakeSliceFileName ( mriFunctionalDataRef this, int inSliceNumber,
        char* inFileNameToSet ) {
  
                                      /* use the path we've saved, the
           stem, the slice number, and a suffix
           to build the file name that this 
           slice should be in. */
  switch ( this->mDataType ) {

  case FunD_tDataType_Short:
    FunD_MakeBShortSliceFileName ( this->mPath, this->mStem,
             inSliceNumber, inFileNameToSet );
    break;

  case FunD_tDataType_Float:
    FunD_MakeBFloatSliceFileName ( this->mPath, this->mStem,
             inSliceNumber, inFileNameToSet );
    break;
  }
}

void FunD_MakeBShortSliceFileName ( char* inPath, char* inStem, 
              int inSliceNumber, char* ioFileName ) {

  sprintf ( ioFileName, "%s/%s_%.3d.%s", 
      inPath, inStem, inSliceNumber, 
      kFileName_BFileShortSuffix );
}

void FunD_MakeBFloatSliceFileName ( char* inPath, char* inStem, 
              int inSliceNumber, char* ioFileName ) {

  sprintf ( ioFileName, "%s/%s_%.3d.%s", 
      inPath, inStem, inSliceNumber, 
      kFileName_BFileFloatSuffix );
}

float FunD_GetValue ( mriFunctionalDataRef this,
    xVoxelRef inFunctionalVoxel,
      int inConditionIndex, int inTimePoint ) {

  float theValue;
  int theIndex;
  theIndex = FunD_CoordsToIndex ( this, inFunctionalVoxel,
            inConditionIndex, inTimePoint );

  switch ( this->mDataType ) {
  case FunD_tDataType_Short:
    theValue = ((short*)this->mData)[theIndex];
    break;
  case FunD_tDataType_Float:
    theValue = ((float*)this->mData)[theIndex];
    break;
  default:
    theValue = 0;
  }
  
  /*
  DebugPrint "getting (%d, %d, %d) cond: %d time: %d index: %d is %2.2f\n",
    xVoxl_GetI(inFunctionalVoxel), 
    xVoxl_GetJ(inFunctionalVoxel),
    xVoxl_GetK(inFunctionalVoxel),
    inConditionIndex, inTimePoint, theIndex, theValue EndDebugPrint;
  */

  return theValue;

}
void FunD_SetValue ( mriFunctionalDataRef this,
          xVoxelRef inFunctionalVoxel,
           int inConditionIndex, int inTimePoint, float inValue ) {

  int theIndex;
  theIndex = FunD_CoordsToIndex ( this, inFunctionalVoxel,
            inConditionIndex, inTimePoint );

  /*
  DebugPrint "setting (%d, %d, %d) cond: %d time: %d index: %d to %2.2f\n",
    xVoxl_GetI(inFunctionalVoxel), 
    xVoxl_GetJ(inFunctionalVoxel),
    xVoxl_GetK(inFunctionalVoxel),
    inConditionIndex, inTimePoint, theIndex, inValue EndDebugPrint;
  */

  switch ( this->mDataType ) {
  case FunD_tDataType_Short:
    ((short*)this->mData)[theIndex] = inValue;
    break;
  case FunD_tDataType_Float:
    ((float*)this->mData)[theIndex] = inValue;
    break;
  }
  
  /* update our min and max */
  if( inValue < this->mMinValue )
    this->mMinValue = inValue;
  if( inValue > this->mMaxValue )
    this->mMaxValue = inValue;
}

// bounds checking
tBoolean FunD_IsTimeResolutionValid ( mriFunctionalDataRef this, int inTimeRes ) {

  if ( inTimeRes <= 0 )
    return FALSE;

  return TRUE;
}

tBoolean FunD_IsNumPreStimTimePointsValid ( mriFunctionalDataRef this, int inNumPoints ) {

  // must be above 0 and less than num time points
  if ( inNumPoints < 0 || inNumPoints >= this->mNumTimePoints )
    return FALSE;

  return TRUE;
}

tBoolean FunD_IsFunctionalVoxelValid ( mriFunctionalDataRef this,
               xVoxelRef            inVoxel ) {

  /*
  // i should be within the col bounds...
  if ( xVoxl_GetI(inVoxel) >= 0 && xVoxl_GetI(inVoxel) < this->mNumCols &&
       // j in the row bounds..
       xVoxl_GetJ(inVoxel) >= 0 && xVoxl_GetJ(inVoxel) < this->mNumRows &&
       // and k in the slice bounds.
       xVoxl_GetK(inVoxel) >= 0 && xVoxl_GetK(inVoxel) < this->mNumSlices )
  */

  // i should be within the col bounds...
  if ( xVoxl_GetFloatX(inVoxel) >= 0 && 
       xVoxl_GetFloatX(inVoxel) < this->mNumCols &&
       // j in the row bounds..
       xVoxl_GetFloatY(inVoxel) >= 0 && 
       xVoxl_GetFloatY(inVoxel) < this->mNumRows &&
       // and k in the slice bounds.
       xVoxl_GetFloatZ(inVoxel) >= 0 && 
       xVoxl_GetFloatZ(inVoxel) < this->mNumSlices )

    return TRUE;
  else 
    return FALSE;

}

tBoolean FunD_IsConditionIndexValid ( mriFunctionalDataRef this, int inConditionIndex ) {

  if ( inConditionIndex >= 0 && inConditionIndex < this->mNumConditions )
    
    return TRUE;
  else
    return FALSE;
}

tBoolean FunD_IsTimePointValid ( mriFunctionalDataRef this, int inTimePoint ) {

  if ( inTimePoint >= 0 && inTimePoint < this->mNumTimePoints )

    return TRUE;
  else
    return FALSE;
}

tBoolean FunD_IsIndexValid ( mriFunctionalDataRef this, int inIndex ) {

  if ( inIndex >= 0 && inIndex < this->mNumDataValues )

    return TRUE;
  else
    return FALSE;
}

