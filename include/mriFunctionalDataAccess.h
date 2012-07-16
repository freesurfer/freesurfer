/**
 * @file  mriFunctionalDataAccess.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/07/16 01:29:17 $
 *    $Revision: 1.32 $
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


#ifndef tkmFunctionalDataAccess_H
#define tkmFunctionalDataAccess_H

#include "xTypes.h"
#include "xVoxel.h"
#include "xList.h"
#include "mriTransform.h"
#include "mriVolume.h"

/* Enable this to turn macros on, see details below. */
#define FUND_USE_MACROS


/* error constants. */
typedef enum
{
  FunD_tErr_NoError = 0,
  FunD_tErr_InvalidSignature,
  FunD_tErr_InvalidPtr,
  FunD_tErr_PathNotFound,
  FunD_tErr_CouldntGuessStem,
  FunD_tErr_DataNotFound,
  FunD_tErr_HeaderNotFound,
  FunD_tErr_CouldntAllocateVolume,
  FunD_tErr_UnrecognizedHeaderFormat,
  FunD_tErr_QuestionableHeaderFormat,
  FunD_tErr_CouldntDetermineDataType,
  FunD_tErr_CouldntAllocateStorage,
  FunD_tErr_CouldntAllocateMRI,
  FunD_tErr_CouldntAllocateMatrix,
  FunD_tErr_CouldntReadMRI,
  FunD_tErr_CouldntReadRegisterFile,
  FunD_tErr_CouldntCalculateDeviations,
  FunD_tErr_ErrorAccessingTransform,
  FunD_tErr_InvalidParameter,
  FunD_tErr_InvalidTimeResolution,
  FunD_tErr_InvalidNumPreStimTimePoints,
  FunD_tErr_InvalidFunctionalVoxel,
  FunD_tErr_InvalidConditionIndex,
  FunD_tErr_InvalidTimePoint,
  FunD_tErr_CouldntTransformMask,
  FunD_tErr_ErrorPerformingFDR,
  FunD_tErr_InvalidErrorCode,
  FunD_tErr_knNumErrorCodes
} FunD_tErr;

typedef enum
{
  FunD_tDataType_Short,
  FunD_tDataType_Float
} FunD_tDataType;

/* This should be kept in sync with the SAMPLE_* constants in mri.h */
typedef enum
{
  FunD_tSampleType_Nearest = 0,
  FunD_tSampleType_Trilinear,
  FunD_tSampleType_Sinc,
  FunD_knNumSampleTypes
} FunD_tSampleType;

/* methods of func idx conversion. this is the rounding method used to
   get the functional index from an ras coord. */
typedef enum
{
  FunD_tConversionMethod_FFF = 0, /* floor, floor, floor */
  FunD_tConversionMethod_Round,   /* rint */
  FunD_tConversionMethod_FCF,     /* floor, ceil, floor */
  FunD_knNumConversionMethods
} FunD_tConversionMethod;

/* Methods of finding the registration. Specified in the New function. */
typedef enum
{
  FunD_tRegistration_None = -1,
  FunD_tRegistration_File = 0,
  FunD_tRegistration_Find,
  FunD_tRegistration_Identity,
  FunD_tRegistration_NoneNeeded,
  FunD_knNumRegistrationTypes
} FunD_tRegistrationType;

#define FunD_ksConversionMethod_FFF   "floor"
#define FunD_ksConversionMethod_Round "round"
#define FunD_ksConversionMethod_FCF   "tkregister"

#define FunD_kSignature 0xe871bc90

#define FunD_knPathLen 1024

#define FunD_knNumFrequencyBins 100

/* the volume struct. */
typedef struct
{

  long mSignature;

  char msFileName[FunD_knPathLen];
  char msRegistrationFileName[FunD_knPathLen];

  /* Conversion method */
  FunD_tConversionMethod mConvMethod;

  /* How to intrepret additional frames of data (see below). */
  int      mNumTimePoints;
  int      mNumConditions;
  tBoolean mbNullConditionPresent;

  /* The data. */
  MRI* mpData;

  /* If we're scalar data, we don't use transformation objects. */
  tBoolean mbScalar;

  /* Sample type to use. */
  FunD_tSampleType mSampleType;

  float mMaxValue;
  float mMinValue;

  /* Addtl meta data. */
  float mTimeResolution;             /* Num seconds between time points. */
  int   mNumPreStimTimePoints;       /* Num of time points before stimuleus */
  char  msSubjectName[1024];
  float mBrightnessScale;
  float mRegistrationPixelSize;
  float mRegistrationThickness;

  /* transformation objects. */
  mriTransformRef mIdxToIdxTransform;
  mriTransformRef mOriginalIdxToIdxTransform;

  /* When we know we have error values present, these are used */
  tBoolean mbErrorDataPresent;
  float** mCovMtx;       /* [conditions-1*timepts][conditions-1*timepts] */

  /* TODO: If the client tells us its bound size, we'll try making a
     resampled volume in client space to speed up data access (avoids
     coordinate conversion and bounds checking. */
  MRI*     mpResampledData;
  tBoolean mbHaveClientBounds;
  int      mClientXMin;
  int      mClientYMin;
  int      mClientZMin;
  int      mClientXMax;
  int      mClientYMax;
  int      mClientZMax;

  /* used to calculate quintiles */
  int* mFrequencies;
  int  mNumBins;

  /* temp storage in conversions */
  xVoxel mTmpVoxel1;
  xVoxel mTmpVoxel2;
  int    mTmpFrame;
  double   mTmpReal;

}
mriFunctionalData, *mriFunctionalDataRef;

/* A note on the format of the data in the MRI struct: functional data
   is probably of multiple frames. These frames are interpreted as
   time points and conditions.

   Additionally, the MRI data could contain error data, which is an
   additional frame of data for each condition, but only the first
   number is valid in each condition. Additionally, the MRI has a set
   of 'null data' and a sigma data at the first condition if error
   data is present.

   Here is how it is ordered:


   WITHOUT ERROR DATA

   condition
     time point
       frame of displayable data

                       time
   frame   condition   point   type
     0         0         0     data
     1         0         1     data
     2         1         0     data
     3         1         1     data

   WITH ERROR DATA

   condition 0
     time point
       null values (ignored)
       sigma values
   condition 1..n
     time point
       frame of displayable data
       std dev (ignored)

                       time
   frame   condition   point   type
     0         0               null (ignored)
     1         0               null (ignored)
     2         0               null (ignored)
     3         0         0     sigma
     4         0         1     sigma
     5         0         2     sigma
     6         1         0     data
     7         1         1     data
     8         1         2     data
     9         1               std dev (ignored)
     10        1               std dev (ignored)
     11        1               std dev (ignored)
     12        2         0     data
     13        2         1     data
     14        2         2     data
     15        2               std dev (ignored)
     16        2               std dev (ignored)
     17        2               std dev (ignored)

*/


/* Allocates and deletes volumes. Volumes can be interpreted as
   anatomically-based volumes that are registered onto a real
   anatomical volume, or as a scalar volume, with a list of indexed
   values. isFileName should be something appropriate for
   MRIread. iRegistrationType is one of the registration types, and is
   used if the volume turns out to be anatomical. isRegistrationFile
   is an optional registration file to be used if the registration
   type is FunD_tRegistration_File. inScalarSize is a number of scalar
   values; if specified (not -1), the size of the volumeis compared to
   this number and if it matches, the volume will be treated as a
   scalar value set, and registration will not be
   use. iAnatomicalVolume is used if the registration type is
   FunD_tRegistration_Identity to calculate an identity registration,
   and can be NULL. */

FunD_tErr FunD_New    ( mriFunctionalDataRef*  opVolume,
                        char*                  isFileName,
                        FunD_tRegistrationType iRegistrationType,
                        char*                  isRegistrationFile,
                        int                    inScalarSize,
                        mriVolumeRef           iAnatomicalVolume,
                        tBoolean               ibIsLeftHemi) ;

FunD_tErr FunD_Delete ( mriFunctionalDataRef*  iopVolume );

/* Parse a stem header if available. */
FunD_tErr FunD_FindAndParseStemHeader_ ( mriFunctionalDataRef this );

/* Get meta information from MRI data. */
FunD_tErr FunD_GuessMetaInformation_   ( mriFunctionalDataRef this );

/* Reads the register.dat file, allocates and initializes matricies. */
FunD_tErr FunD_ParseRegistrationAndInitMatricies_ ( mriFunctionalDataRef this,
    FunD_tRegistrationType iType,
    mriVolumeRef iAnatomicalVolume );

/* Restores registration matrix to original copy. */
FunD_tErr FunD_RestoreRegistration ( mriFunctionalDataRef this );

/* Calcs the deviations. Needs sigma and CovMtx, should be done after
   reading the header _and_ data files. Deviations are calced by
   multplying sigma by the square root of the diagonal of the CovMtx,
   then indexing in by time major (time,cond 0,0, 1,0, 2,0...) */
FunD_tErr FunD_CalcDeviations_ ( mriFunctionalDataRef this );


/* Sets the client information so that functional volume can use a
   resampled volume of the data in client space. */
FunD_tErr FunD_SetClientCoordBounds ( mriFunctionalDataRef this,
                                      int                  inXMin,
                                      int                  inYMin,
                                      int                  inZMin,
                                      int                  inXMax,
                                      int                  inYMax,
                                      int                  inZMax );

/* Sets the value conversion method for getting integer values. */
FunD_tErr FunD_SetConversionMethod ( mriFunctionalDataRef this,
                                     FunD_tConversionMethod iMethod );

/* Gets and sets the sample type for interpolating data. Only is
   used with the GetSampledData function. */
FunD_tErr FunD_GetSampleType ( mriFunctionalDataRef this,
                               FunD_tSampleType*    oType );

FunD_tErr FunD_SetSampleType ( mriFunctionalDataRef this,
                               FunD_tSampleType     iType );

/* Whether or not the volume is scalar and has been reshaped. */
FunD_tErr FunD_IsScalar ( mriFunctionalDataRef this,
                          tBoolean*            obScalar );

/* Use this if the coordinates in the client space are actually in
   native tkRegRAS space (e.g. surface). This will set the A->RAS
   matrix in the IdxToIdx transform to identity. */
FunD_tErr FunD_ClientSpaceIsTkRegRAS ( mriFunctionalDataRef this );

/* Value accessors. */
FunD_tErr FunD_GetData                 ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iTimePoint,
    int                  iCondition,
    float*               oData );
FunD_tErr FunD_GetSampledData          ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iTimePoint,
    int                  iCondition,
    float*               oData );
FunD_tErr FunD_GetDataForAllTimePoints ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iCondition,
    float*               oaData );

FunD_tErr FunD_GetDeviation                 ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iCondition,
    int                  iTimePoint,
    float*               oValue );
FunD_tErr FunD_GetDeviationForAllTimePoints ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    int                  iCondition,
    float*               oaData );


/* Smooths the data. */
FunD_tErr FunD_Smooth ( mriFunctionalDataRef this,
                        int                  iTimePoint,
                        int                  iCondition,
                        float                iSigma );

/* Normalizes the data over all time points. */
FunD_tErr FunD_NormalizeOverAll ( mriFunctionalDataRef this );


/* converts a time point index to a second based on the time
   resolution and time resolution */
FunD_tErr FunD_ConvertTimePointToSecond ( mriFunctionalDataRef this,
    int                  iTimePoint,
    float*               oSecond );
FunD_tErr FunD_ConvertSecondToTimePoint ( mriFunctionalDataRef this,
    float                iSecond,
    int*                 oTimePoint );

/* setting these values changes the way time points are converted to
   seconds. */
FunD_tErr FunD_SetTimeResolution       ( mriFunctionalDataRef this,
    float                inTimeRes );
FunD_tErr FunD_SetNumPreStimTimePoints ( mriFunctionalDataRef this,
    int                  inNumPoints );

/* accessors for volume information. */
FunD_tErr FunD_GetSubjectName          ( mriFunctionalDataRef this,
    char*                out );
FunD_tErr FunD_GetNumTimePoints        ( mriFunctionalDataRef this,
    int*                 out );
FunD_tErr FunD_GetNumConditions        ( mriFunctionalDataRef this,
    int*                 out );
FunD_tErr FunD_GetTimeResolution       ( mriFunctionalDataRef this,
    float*               out );
FunD_tErr FunD_GetNumPreStimTimePoints ( mriFunctionalDataRef this,
    int*                 out );
FunD_tErr FunD_GetValueRange           ( mriFunctionalDataRef this,
    float*               outMin,
    float*               outMax );
FunD_tErr FunD_IsErrorDataPresent      ( mriFunctionalDataRef this,
    tBoolean*            oPresent );
/* Gets bounds in client coords. Ghis is actually a bounding box in
   client space. If the functional data is oblique to the client, it
   will have coords that arn't valid in functional space, but it is a
   good place to start iterating for the in anatomical space. */
FunD_tErr FunD_GetBoundsInClientSpace ( mriFunctionalDataRef this,
                                        xVoxelRef             outBeginCorner,
                                        xVoxelRef             outEndCorner );

/* Gets the value at a quintile. Assumes quintile is within 0-100. */
FunD_tErr FunD_GetValueAtPercentile ( mriFunctionalDataRef this,
                                      float                inPercent,
                                      float*               outValue );


/* Saves the registration to file, making a backup if it already
   exists. */
FunD_tErr FunD_SaveRegistration ( mriFunctionalDataRef this );

/* Sets registration to identity matrix. */
FunD_tErr FunD_SetRegistrationToIdentity ( mriFunctionalDataRef this );

/* Applies a transformation to the registration, which can later be
   saved out as a new registration. */
FunD_tErr FunD_ApplyTransformToRegistration ( mriFunctionalDataRef this,
    MATRIX*             iTransform );

/* Geometric trasnformations on the registration. */
FunD_tErr FunD_TranslateRegistration        ( mriFunctionalDataRef this,
    float                ifDistance,
    tAxis                iAxis );
FunD_tErr FunD_RotateRegistration           ( mriFunctionalDataRef this,
    float                ifDegrees,
    tAxis                iAxis,
    xVoxelRef   iCenterFuncRAS );
FunD_tErr FunD_ScaleRegistration            ( mriFunctionalDataRef this,
    float                ifFactor,
    tAxis                iAxis );


/* Get the FDR threshold for a frame. */
FunD_tErr FunD_CalcFDRThreshold ( mriFunctionalDataRef this,
                                  int                  iCondition,
                                  int                  iTimePoint,
                                  int                  iSign,
                                  float                iRate,
                                  MRI*                 iMaskVolume,
                                  float*               oThresholdMin );

/* This function doesn't seem to work very well but is left here for
   reference . */
#if 0
FunD_tErr FunD_GetPercentileOfValue ( mriFunctionalDataRef this,
                                      float                inValue,
                                      float*               outPercentile );
#endif


/* internal functions */

/* If the dimensions of the data, multiplied, match the number of
   values here, the data is reshaped into the x direction and treated
   as a scalar volume with no transform. */
FunD_tErr FunD_ReshapeIfScalar_ ( mriFunctionalDataRef this,
                                  int                  inNumValues,
                                  tBoolean*            obReshaped,
                                  tBoolean             ibIsLeftHemisphere);

/* Resamples the data into client space for fast lookup. */
FunD_tErr FunD_ResampleData_ ( mriFunctionalDataRef this );

/* Calculate a histogram for the data, dividing it up into a number of
   bins. Allocates internal memory. */
FunD_tErr FunD_CalcFrequencies_ ( mriFunctionalDataRef this,
                                  int                  inNumBins );

/* Converts between client and functional voxel coordinates. */
/* NOTE: FunD_ConvertClientToFuncIdx_ is now a macro */
void FunD_ConvertFuncIdxToClient_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iFuncIdx,
                                    xVoxelRef            oClientVox );
void FunD_ConvertClientToFuncRAS_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iClientVox,
                                    xVoxelRef            oFuncRAS );


void FunD_GetSigma_ ( mriFunctionalDataRef this,
                      xVoxelRef            iFuncIdx,
                      int                  iTimePoint,
                      float*               oSigma );

#define FunD_GetSigmaFrameNumber(iTimePoint,oFrame) \
  if( this->mbErrorDataPresent ) { \
    *(oFrame) = this->mNumTimePoints + (iTimePoint); \
  } else { \
    *(oFrame) = 0 ; \
  }


/* Note that the functions in this section are implemented as
   functions and macros. The functions are slower but safer, and the
   macros are faster but don't make any checks. So you should test
   with the functions and then build with the macros turned on. And
   keep the macro and functions version synced in
   development. _grin_ Turn the flag on at the top of this file. */
#ifndef FUND_USE_MACROS

void FunD_GetValue_ ( mriFunctionalDataRef this,
                      MRI*                 iData,
                      xVoxelRef            iIdx,
                      int                  iCondition,
                      int                  iTimePoint,
                      float*               oValue );

void FunD_GetSampledValue_ ( mriFunctionalDataRef this,
                             MRI*                 iData,
                             xVoxelRef            iIdx,
                             int                  iCondition,
                             int                  iTimePoint,
                             float*               oValue );

void FunD_SetValue_ ( mriFunctionalDataRef this,
                      MRI*                 iData,
                      xVoxelRef            iIdx,
                      int                  iCondition,
                      int                  iTimePoint,
                      float                oValue );

/* This converts to functional index and performs the proper rounding
   on the result. */
void FunD_ConvertClientToFuncIdx_ ( mriFunctionalDataRef this,
                                    xVoxelRef            iClientVox,
                                    xVoxelRef            oFuncIdx );
/* This converts to func idx with floating info intact. */
void FunD_ConvertClientToFloatFuncIdx_ ( mriFunctionalDataRef this,
    xVoxelRef            iClientVox,
    xVoxelRef            oFuncIdx );

#define FunD_GetDataFrameNumber(iCondition,iTimePoint,oFrame) \
  if( this->mbErrorDataPresent ) { \
    *(oFrame) = ((iCondition) * 2 * this->mNumTimePoints) + (iTimePoint); \
  } else { \
    *(oFrame) = ((iCondition) * this->mNumTimePoints) + (iTimePoint); \
  }


#else /* FUND_USE_MACROS */


#define FunD_GetDataFrameNumber(iCondition,iTimePoint,oFrame) \
  if( this->mbErrorDataPresent ) { \
    *(oFrame) = ((iCondition) * 2 * this->mNumTimePoints) + (iTimePoint); \
  } else { \
    *(oFrame) = ((iCondition) * this->mNumTimePoints) + (iTimePoint); \
  }

#define FunD_GetValue_(this,iData,iIdx,inCondition,inTimePoint,oValue) \
  if( this->mbErrorDataPresent ) { \
    this->mTmpFrame = \
 ((inCondition) * 2 * this->mNumTimePoints) + (inTimePoint); \
  } else { \
    this->mTmpFrame = (inCondition * this->mNumTimePoints) + inTimePoint; \
  } \
  MRIsampleVolumeFrameType( (iData), \
       xVoxl_GetFloatX((iIdx)),\
       xVoxl_GetFloatY((iIdx)),\
       xVoxl_GetFloatZ((iIdx)),\
       this->mTmpFrame, \
       SAMPLE_NEAREST, \
       &this->mTmpReal ); \
  *(oValue) = this->mTmpReal;

#define FunD_GetSampledValue_(this,iData,iIdx,inCondition,inTimePoint,oValue) \
  if( this->mbErrorDataPresent ) { \
    this->mTmpFrame = \
((inCondition) * 2 * this->mNumTimePoints) + (inTimePoint); \
  } else { \
    this->mTmpFrame = (inCondition * this->mNumTimePoints) + inTimePoint; \
  } \
  MRIsampleVolumeFrameType( (iData), \
       xVoxl_GetFloatX((iIdx)),\
       xVoxl_GetFloatY((iIdx)),\
       xVoxl_GetFloatZ((iIdx)),\
       this->mTmpFrame, \
       (int)this->mSampleType, \
       &this->mTmpReal ); \
  *(oValue) = this->mTmpReal;


#define FunD_SetValue_(this,iData,iIdx,inCondition,inTimePoint,iValue) \
  if( this->mbErrorDataPresent ) { \
    this->mTmpFrame = \
((inCondition) * 2 * this->mNumTimePoints) + (inTimePoint); \
  } else { \
    this->mTmpFrame = (inCondition * this->mNumTimePoints) + inTimePoint; \
  } \
  switch( iData->type ) { \
    case MRI_UCHAR: \
      MRIseq_vox( iData, (int)(iIdx)->mfX, (int)(iIdx)->mfY, \
                  (int)(iIdx)->mfZ, this->mTmpFrame) = iValue; \
      break; \
    case MRI_INT: \
      MRIIseq_vox( iData, (int)(iIdx)->mfX,  (int)(iIdx)->mfY, \
                  (int)(iIdx)->mfZ, this->mTmpFrame) = iValue; \
      break; \
    case MRI_LONG: \
      MRILseq_vox( iData, (int)(iIdx)->mfX,  (int)(iIdx)->mfY, \
                  (int)(iIdx)->mfZ, this->mTmpFrame) = iValue; \
      break; \
    case MRI_FLOAT: \
      MRIFseq_vox( iData, (int)(iIdx)->mfX,  (int)(iIdx)->mfY, \
                  (int)(iIdx)->mfZ, this->mTmpFrame) = iValue; \
      break; \
    case MRI_SHORT: \
      MRISseq_vox( iData, (int)(iIdx)->mfX,  (int)(iIdx)->mfY, \
                  (int)(iIdx)->mfZ, this->mTmpFrame) = iValue; \
      break; \
    default: \
      break ; \
    }


#define FunD_ConvertClientToFuncIdx_(this,iClientVox,oFuncIdx) \
  \
  Trns_ConvertAtoB( this->mIdxToIdxTransform, \
      (iClientVox), &this->mTmpVoxel2 );\
  \
  switch( this->mConvMethod ) {\
  case FunD_tConversionMethod_FFF:\
    xVoxl_SetFloat( (oFuncIdx), \
      floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),\
      floor(xVoxl_GetFloatY(&this->mTmpVoxel2)),\
      floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );\
    break;\
  case FunD_tConversionMethod_Round:\
    xVoxl_SetFloat( (oFuncIdx), \
      rint(xVoxl_GetFloatX(&this->mTmpVoxel2)),\
      rint(xVoxl_GetFloatY(&this->mTmpVoxel2)),\
      rint(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );\
    break;\
  case FunD_tConversionMethod_FCF:\
    xVoxl_SetFloat( (oFuncIdx), \
      floor(xVoxl_GetFloatX(&this->mTmpVoxel2)),\
      ceil(xVoxl_GetFloatY(&this->mTmpVoxel2)),\
      floor(xVoxl_GetFloatZ(&this->mTmpVoxel2)) );\
    break;\
  default:\
    break;\
  }

#define FunD_ConvertClientToFloatFuncIdx_(this,iClientVox,oFuncIdx) \
  Trns_ConvertAtoB( this->mIdxToIdxTransform, \
      (iClientVox), (oFuncIdx) );

#endif

FunD_tErr FunD_DebugPrint ( mriFunctionalDataRef this );

char* FunD_GetErrorString ( FunD_tErr inErr );

FunD_tErr FunD_Verify         ( mriFunctionalDataRef this );
FunD_tErr FunD_VerifyFuncIdx_ ( mriFunctionalDataRef this,
                                xVoxelRef            iFuncIdx );
FunD_tErr FunD_VerifyTimePoint ( mriFunctionalDataRef this,
                                 int                  iTimePoint );
FunD_tErr FunD_VerifyCondition ( mriFunctionalDataRef this,
                                 int                  iCondition );

#endif
