#ifndef mriVolume_h
#define mriVolume_h

#include "volume_io.h"
#include "mri.h"
#include "xTypes.h"
#include "xVoxel.h"
#include "mriTypes.h"
#include "mriTransform.h"
#include "transform.h"

#define VOLM_USE_MACROS

typedef enum {
  Volm_tErr_Invalid = -1,
  Volm_tErr_NoErr = 0,
  Volm_tErr_InvalidSignature,
  Volm_tErr_InvalidParamater,
  Volm_tErr_InvalidIdx,
  Volm_tErr_AllocationFailed,
  Volm_tErr_CouldntReadVolume,
  Volm_tErr_CouldntCopyVolume,
  Volm_tErr_CouldntReadTransform,
  Volm_tErr_CouldntCopyTransform,
  Volm_tErr_CouldntNormalizeVolume,
  Volm_tErr_CouldntExportVolumeToCOR,
  Volm_tErr_NormValuesNotPresent,
  Volm_tErr_ScannerTransformNotPresent,
  Volm_tErr_IdxToRASTransformNotPresent,
  
  Volm_knNumErrorCodes
} Volm_tErr;

#define Volm_knErrStringLen 1024

typedef BUFTYPE Volm_tValue;
typedef Volm_tValue *Volm_tValueRef;
#define Volm_knMaxValue   255
#define Volm_knNumValues  256

#define Volm_kSignature 0xd9d6b2a9

#define Volm_kfDefaultBrightness 0.35
#define Volm_kfDefaultContrast   12.0

typedef struct {
  
  long mSignature;
  
  int mnDimensionX;
  int mnDimensionY;
  int mnDimensionZ;
  
  MRI*           mpMriValues;    /* normalized color values (0-255) */
  BUFTYPE*       mpSnapshot;      /* copy of normalized values */
  Volm_tValueRef mpMaxValues;     /* max projection. 3 planes, each dimension
             of a slice. 
             mpMaxValues[orientation][x][y] */
  
  mriTransformRef    mIdxToRASTransform;           /* idx -> ras (by MRI) */
  //  LTA*               mDisplayTransform;            /* buf -> index */
  mriTransformRef    mDisplayTransform;            /* buf -> index */
  mriTransformRef    mMNITalLtzToRealTalTransform; /* tal (z<0) -> real tal */
  mriTransformRef    mMNITalGtzToRealTalTransform; /* tal (z>0) -> real tal */
  mriTransformRef    mScannerTransform;            /* idx -> scnaner */
  
  char msSubjectName[mri_knSubjectNameLen];
  char msVolumeName[mri_knSubjectNameLen];
  char msOriginalPath[mri_knPathLen];
  
  float mfColorBrightness;                 /* threshold */
  float mfColorContrast;                   /* squash */
  xColor3f maColorTable[Volm_knNumValues]; /* color for each value */
  MATRIX   *m_resample_orig ;              /* matrix that takes 256^3->original volume */
  MATRIX   *m_resample ;                   /* matrix that takes 256^3->original volume */
  float     max_val ;                      /* max value in mpMriValues */
  float     min_val ;                      /* min value in mpMriValues */
} mriVolume, *mriVolumeRef;

typedef enum {
  Volm_tVisitComm_Continue = 0,
  Volm_tVisitComm_SkipRestOfRow,
  Volm_tVisitComm_SkipRestOfPlane,
  Volm_tVisitComm_Stop,
  Volm_knNumVisitCommands
} Volm_tVisitCommand;
typedef Volm_tVisitCommand(*Volm_tVisitFunction)(xVoxelRef  iAnaIdx,
             float      iValue,
             void*      ipData);

Volm_tErr Volm_New        ( mriVolumeRef* opVolume );
Volm_tErr Volm_Delete     ( mriVolumeRef* iopVolume );
Volm_tErr Volm_DeepClone  ( mriVolumeRef  this, 
          mriVolumeRef* opVolume );

Volm_tErr Volm_ImportData      ( mriVolumeRef this,
         char*        isSource );
Volm_tErr Volm_ExportNormToCOR ( mriVolumeRef this,
         char*        isPath );

Volm_tErr Volm_LoadDisplayTransform ( mriVolumeRef this,
              char*        isSource );
Volm_tErr Volm_UnloadDisplayTransform ( mriVolumeRef this );

void Volm_GetColorAtIdx        ( mriVolumeRef     this,
         xVoxelRef        iIdx,
         xColor3fRef      oColor );
void Volm_GetColorAtXYSlice    ( mriVolumeRef     this,
         mri_tOrientation iOrientation,
         xPoint2nRef      iPoint,
         int              inSlice,
         xColor3fRef      oColor );
void Volm_GetMaxColorAtIdx     ( mriVolumeRef     this,
         xVoxelRef        iIdx,
         mri_tOrientation iOrientation,
         xColor3fRef      oColor );
void Volm_GetMaxColorAtXYSlice ( mriVolumeRef     this,
         mri_tOrientation iOrientation,
         xPoint2nRef      iPoint,
         int              inSlice,
         xColor3fRef      oColor );

Volm_tErr Volm_GetDimensions        ( mriVolumeRef this,
              int*         onDimensionX, 
              int*         onDimensionY, 
              int*         onDimensionZ );

/* get a normal value. before calling the unsafe version, make sure
   the volume is valid, the index is in bounds, and you're passing
   a valid pointer. */
Volm_tErr Volm_GetValueAtIdx       ( mriVolumeRef this,
             xVoxelRef    iIdx,
             float*       oValue );
Volm_tErr Volm_GetValueAtIdxUnsafe ( mriVolumeRef this,
             xVoxelRef    iIdx,
             float*       oValue );
Volm_tErr Volm_SetValueAtIdx       ( mriVolumeRef this,
             xVoxelRef    iIdx,
             Volm_tValue );

/* coordinate conversion. idx stands for index and is the 0->dimension-1
   index of the volume. RAS space, aka world space, is centered on the
   center of the volume and is in mm. mni tal is mni's version of 
   talairach space. the other tal is mni tal with a small modification. */
Volm_tErr Volm_ConvertIdxToRAS     ( mriVolumeRef this,
             xVoxelRef    iIdx,
             xVoxelRef    oRAS );
Volm_tErr Volm_ConvertRASToIdx     ( mriVolumeRef this,
             xVoxelRef    iAndIdx,
             xVoxelRef    oIdx );
Volm_tErr Volm_ConvertIdxToMNITal  ( mriVolumeRef this,
             xVoxelRef    iIdx,
             xVoxelRef    oMNITal );
Volm_tErr Volm_ConvertIdxToTal     ( mriVolumeRef this,
             xVoxelRef    iIdx,
             xVoxelRef    oTal );
Volm_tErr Volm_ConvertTalToIdx     ( mriVolumeRef this,
             xVoxelRef    iTal,
             xVoxelRef    oIdx );
Volm_tErr Volm_ConvertIdxToScanner ( mriVolumeRef this,
             xVoxelRef    iIdx,
             xVoxelRef    oScanner );

Volm_tErr Volm_GetIdxToRASTransform ( mriVolumeRef     this,
              mriTransformRef* opTransform );

/* calls the parameter function for every voxel, passing the voxel, the
   value, and the pointer passed to it. */
Volm_tErr Volm_VisitAllVoxels ( mriVolumeRef        this,
        Volm_tVisitFunction iFunc,
        void*               ipData );


Volm_tErr Volm_FindMaxValues ( mriVolumeRef this );

Volm_tErr Volm_MakeColorTable ( mriVolumeRef this );

Volm_tErr Volm_SetBrightnessAndContrast ( mriVolumeRef this,
            float        ifBrightness,
            float        ifContrast );

Volm_tErr Volm_SaveToSnapshot      ( mriVolumeRef this );
Volm_tErr Volm_RestoreFromSnapshot ( mriVolumeRef this );

Volm_tErr Volm_Rotate    ( mriVolumeRef     this,
         mri_tOrientation iAxis,
         float            ifDegrees );
Volm_tErr Volm_Threshold ( mriVolumeRef     this,
         Volm_tValue      iThreshold,
         tBoolean         ibAbove,
         Volm_tValue      iNewValue );
Volm_tErr Volm_Flip      ( mriVolumeRef     this,
         mri_tOrientation iAxis );

Volm_tErr Volm_CopySubjectName ( mriVolumeRef this,
         char*        oSubjectName,
         int          inDestLen );
Volm_tErr Volm_CopyVolumeName  ( mriVolumeRef this,
         char*        oVolumeName,
         int          inDestLen );
Volm_tErr Volm_CopySourceDir   ( mriVolumeRef this,
         char*        oSourceDir,
         int          inDestLen );

Volm_tErr Volm_SetSubjectName           ( mriVolumeRef this,
            char*        isName );
Volm_tErr Volm_SetVolumeName            ( mriVolumeRef this,
            char*        isName );
Volm_tErr Volm_ExtractAndSetSubjectName ( mriVolumeRef this,
            char*        isSource );
Volm_tErr Volm_ExtractAndSetVolumeName  ( mriVolumeRef this,
            char*        isSource );

#ifndef VOLM_USE_MACROS

/* safer function versions of main accessing and setting functions */

void Volm_ConvertIdxToXYSlice_ ( xVoxelRef         iIdx,
         mri_tOrientation  iOrientation,
         xPoint2nRef       oPoint,
         int*              onSlice );
void Volm_ConvertXYSliceToIdx_ ( mri_tOrientation  iOrientation,
         xPoint2nRef       iPoint,
         int               inSlice,
         xVoxelRef         oIdx );

Volm_tValue Volm_GetNormValueAtXYSlice_ ( mriVolumeRef this,
            mri_tOrientation  iOrientation,
            xPoint2nRef       iPoint,
            int               inSlice );
Volm_tValue Volm_GetNormValueAtIdx_     ( mriVolumeRef     this,
            xVoxelRef        iIdx );
void        Volm_SetNormValueAtIdx_     ( mriVolumeRef     this,
            xVoxelRef        iIdx,
            Volm_tValue      iValue );
float       Volm_GetRawValueAtIdx_      ( mriVolumeRef     this,
            xVoxelRef        iIdx );

Volm_tValue Volm_GetMaxValueAtXYSlice_ ( mriVolumeRef this,
           mri_tOrientation iOrientation, 
           xPoint2nRef      iPoint );
void        Volm_SetMaxValueAtXYSlice_ ( mriVolumeRef this,
           mri_tOrientation iOrientation, 
           xPoint2nRef      iPoint,
           Volm_tValue      iValue );
int Volm_GetMaxValueIndex_ ( mriVolumeRef     this,
           mri_tOrientation iOrientation, 
           xPoint2nRef      iPoint );
#else /* VOLM_USE_MACROS */

/* macro version  of main accessing and setting functions */

#define Volm_ConvertIdxToXYSlice_(iIdx,iOrientation,oPoint,onSlice) \
  switch( iOrientation ) {                          \
  case mri_tOrientation_Coronal:                    \
    (oPoint)->mnX = (iIdx)->mfX;                    \
    (oPoint)->mnY = (iIdx)->mfY;                    \
    *(onSlice)    = (iIdx)->mfZ;                    \
    break;                                          \
  case mri_tOrientation_Horizontal:                 \
    (oPoint)->mnX = (iIdx)->mfX;                    \
    (oPoint)->mnY = (iIdx)->mfZ;                    \
    *(onSlice)    = (iIdx)->mfY;                    \
    break;                                          \
  case mri_tOrientation_Sagittal:                   \
    (oPoint)->mnX = (iIdx)->mfZ;                    \
    (oPoint)->mnY = (iIdx)->mfY;                    \
    *(onSlice)    = (iIdx)->mfX;                    \
    break;                                          \
  default:                                          \
    DebugPrintStack;                                \
    DebugPrint( ("Volm_ConvertIdxToXYSlice_ called with invalid " \
     "orientation %d", iOrientation) );               \
    (oPoint)->mnX = (oPoint)->mnY = *(onSlice) = 0;               \
    break;                                          \
  }


#define Volm_ConvertXYSliceToIdx_(iOrientation,iPoint,inSlice,oIdx) \
  switch( iOrientation ) {                                          \
  case mri_tOrientation_Coronal:                                    \
    (oIdx)->mfX = (iPoint)->mnX;                                    \
    (oIdx)->mfY = (iPoint)->mnY;                                    \
    (oIdx)->mfZ = (inSlice);                                        \
    break;                                                          \
  case mri_tOrientation_Horizontal:                                 \
    (oIdx)->mfX = (iPoint)->mnX;                                    \
    (oIdx)->mfY = (inSlice);                                        \
    (oIdx)->mfZ = (iPoint)->mnY;                                    \
    break;                                                          \
  case mri_tOrientation_Sagittal:                                   \
    (oIdx)->mfX = (inSlice);                                        \
    (oIdx)->mfY = (iPoint)->mnY;                                    \
    (oIdx)->mfZ = (iPoint)->mnX;                                    \
    break;                                                          \
  default:                                                          \
    DebugPrintStack;                                                \
    DebugPrint( ("Volm_ConvertXYSliceToIdx_ called with invalid "   \
     "orientation %d", iOrientation) );                 \
    xVoxl_Set( (oIdx), 0, 0, 0 );                                   \
    break;                                                          \
  }

#define Volm_GetNormValueAtXYSlice_(iOrientation,iPoint,inSlice,oValue) \
  switch( iOrientation ) {                                          \
  case mri_tOrientation_Coronal:                                    \
    *(oValue) = MRIvox( this->mpMriValues, (iPoint)->mnX,          \
            (iPoint)->mnY, (inSlice) );                 \
    break;                                                          \
  case mri_tOrientation_Horizontal:                                 \
    *(oValue) = MRIvox( this->mpMriValues, (iPoint)->mnX,          \
            (inSlice), (iPoint)->mnY );                 \
    break;                                                          \
  case mri_tOrientation_Sagittal:                                   \
    *(oValue) = MRIvox( this->mpMriValues, (inSlice),              \
            (iPoint)->mnY, (iPoint)->mnX );             \
    break;                                                          \
  default:                                                          \
    *(oValue) = 0;                                                  \
    break;                                                          \
  }


#define Volm_GetNormValueAtIdx_(this,iIdx) \
       (Volm_tValue) MRIvox( this->mpMriValues, xVoxl_GetX(iIdx), \
           xVoxl_GetY(iIdx), xVoxl_GetZ(iIdx) )
#define Volm_GetSincNormValueAtIdx_(this,iIdx,irValue) \
       MRIsincSampleVolume(this->mpMriValues,xVoxl_GetFloatX(iIdx), \
                      xVoxl_GetFloatY(iIdx),xVoxl_GetFloatZ(iIdx),2,irValue)
#define Volm_SetNormValueAtIdx_(this,iIdx,iValue) \
                     MRIvox( this->mpMriValues, xVoxl_GetX(iIdx), \
                       xVoxl_GetY(iIdx), xVoxl_GetZ(iIdx) ) = iValue
float       Volm_GetRawValueAtIdx_      ( mriVolumeRef     this,
            xVoxelRef        iIdx );

#define Volm_GetMaxValueAtXYSlice_(this,iOrientation,iPoint) \
    this->mpMaxValues[(iOrientation * this->mnDimensionX * this->mnDimensionY) \
                       + ((iPoint)->mnY * this->mnDimensionX) + (iPoint)->mnX]
#define Volm_SetMaxValueAtXYSlice_(this,iOrientation,iPoint,iValue ) \
    this->mpMaxValues[(iOrientation * this->mnDimensionX * this->mnDimensionY) \
               + ((iPoint)->mnY * this->mnDimensionX) + (iPoint)->mnX] = iValue

#define Volm_ApplyDisplayTransform_(this,iIdx,oIdx) \
    Trns_ConvertBtoA( this->mDisplayTransform, iIdx, oIdx )

#endif /* VOLM_USE_MACROS */

Volm_tErr Volm_Verify     ( mriVolumeRef this );
Volm_tErr Volm_VerifyIdx  ( mriVolumeRef this,
          xVoxelRef    iIdx );
Volm_tErr Volm_VerifyIdx_ ( mriVolumeRef this,
          xVoxelRef    iIdx );
char* Volm_GetErrorString ( Volm_tErr ieCode );


#endif
