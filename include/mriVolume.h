/**
 * @file  mriVolume.h
 * @brief declaration of mri volume structs
 *
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.44 $
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


#ifndef mriVolume_h
#define mriVolume_h

#include "mri.h"
#include "xTypes.h"
#include "xVoxel.h"
#include "mriTypes.h"
#include "mriTransform.h"
#include "transform.h"

/* Enable this to turn macros on, see details below. */
#define VOLM_USE_MACROS

typedef enum
{
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
  Volm_tErr_MRIVolumeNotPresent,
  Volm_tErr_ScannerTransformNotPresent,
  Volm_tErr_IdxToRASTransformNotPresent,
  Volm_tErr_FloodVisitCommandNotSupported,
  Volm_tErr_FsenvNotFoundWhileSearchingForTransformDest,
  Volm_tErr_MNI305NotFoundWhilePopulatingTramsformDest,
  Volm_tErr_ErrorReadingMNI305WhilePopulatingTramsformDest,
  Volm_knNumErrorCodes
} Volm_tErr;

#define Volm_knErrStringLen 1024

typedef float Volm_tValue;
typedef Volm_tValue *Volm_tValueRef;
#define Volm_knMaxValue   255
#define Volm_knNumValues  256
#define Volm_knNumColorTableEntries  256

#define Volm_kSignature 0xd9d6b2a9

#define Volm_kfDefaultBrightness 0.35
#define Volm_kfDefaultContrast   12.0

#define Volm_knMaxFloodIteration 20000

/* Function definition and return codes for a visit callback. The user
   supplies their own visit function and it is called in the flood and
   VisitAll algorithms. */
typedef enum
{
  Volm_tVisitComm_Continue = 0,
  Volm_tVisitComm_SkipRestOfRow,
  Volm_tVisitComm_SkipRestOfPlane,
  Volm_tVisitComm_Stop,
  Volm_knNumVisitCommands
} Volm_tVisitCommand;
typedef Volm_tVisitCommand(*Volm_tVisitFunction)(xVoxelRef  iMRIIdx,
    float      iValue,
    void*      ipData);

/* Parameters for a generic flood algorithm. */
typedef enum
{
  Volm_tValueComparator_Invalid = -1,
  Volm_tValueComparator_LTE = 0,
  Volm_tValueComparator_EQ,
  Volm_tValueComparator_GTE,
  Volm_knNumValueComparators,
} Volm_tValueComparator;
/* Optional comparator function for the flood function. Return true
   for a successful flood, i.e. it should continue on this voxel,and
   false if not. */
typedef tBoolean(*Volm_tFloodComparatorFunction)(xVoxelRef  iMRIIdx,
    float      iValue,
    void*      ipData);

/* This should be kept in sync with the SAMPLE_* constants in mri.h */
typedef enum
{
  Volm_tSampleType_Nearest = 0,
  Volm_tSampleType_Trilinear,
  Volm_tSampleType_Sinc,
  Volm_knNumSampleTypes
} Volm_tSampleType;

typedef enum
{
  Volm_tResampleMethod_RAS = 0,
  Volm_tResampleMethod_Slice,
  Volm_knNumResampleMethods
} Volm_tResampleMethod;

typedef struct
{
  xVoxel           mSourceIdx;       /* The starting voxel */
  float            mfSourceValue;    /* The value at the starting voxel */
  float            mfFuzziness;      /* The fuzziness for the comparator */
  Volm_tValueComparator mComparatorType;          /* Compare type */
  Volm_tFloodComparatorFunction mComparatorFunc;  /* Compare function */
  void*            mComparatorFuncData;           /* Compare function data */
  float            mfMaxDistance;    /* Max distance of flood (in voxels) */
  tBoolean         mb3D;             /* Should the fill be in 3D? */
  mri_tOrientation mOrientation;     /* If not, which plane? */

  Volm_tVisitFunction  mpFunction;     /* User function, called for each */
  void*                mpFunctionData; /* visited. */

}
Volm_tFloodParams;

typedef struct
{

  long mSignature;

  int mnDimensionX;
  int mnDimensionY;
  int mnDimensionZ;
  int mnDimensionFrame;

  MRI*     mpMriValues;     /* MRI struct */
  float*   mpSnapshot;      /* copy of values */
  float**   mpMaxValuesX;    /* max projection in x plane, [y][z] */
  float**   mpMaxValuesY;    /* max projection in y plane, [x][z] */
  float**   mpMaxValuesZ;    /* max projection in z plane, [x][y] */


  Volm_tSampleType mSampleType;  /* How to sample the volume */
  mriTransformRef    mMRIIdxToAnaIdxTransform; /* MRI -> ana (via scnr RAS) */
  mriTransformRef    mMNITalLtzToRealTalTransform; /* tal (z<0) -> real tal */
  mriTransformRef    mMNITalGtzToRealTalTransform; /* tal (z>0) -> real tal */
  mriTransformRef    mScannerTransform;            /* idx -> scnaner */
  TRANSFORM* mDisplayTransform;

  char msSubjectName[mri_knSubjectNameLen];
  char msVolumeName[mri_knSubjectNameLen];
  char msOriginalPath[mri_knPathLen];

  float    mfColorBrightness;               /* threshold */
  float    mfColorContrast;                 /* squash */
  float    mfColorMin;                      /* min color value for col table */
  float    mfColorMax;                      /* max color value for col table */
  xColor3f mafColorTable[Volm_knNumColorTableEntries];
  xColor3n manColorTable[Volm_knNumColorTableEntries];
  MATRIX   *m_resample_orig;                /* 256^3->original volume */
  MATRIX   *m_resample;                     /* 256^3->original volume */
  MATRIX   *m_resample_inv;                 /* original volume->256^3 */
  float    mfMinValue;                      /* max value in mpMriValues */
  float    mfMaxValue;                      /* min value in mpMriValues */

  /* There are two versions of the resample matrices; one that takes
     into effect the RAS transform and one that doesn't. Normally you
     would want to use the RAS transform, but the other is available
     for viewing the volume in its slice orientation. They both have
     pixel size component, however. */
  Volm_tResampleMethod mResampleMethod;
  MATRIX*              mResampleToRAS;
  MATRIX*              mResampleToRASOrig;
  MATRIX*              mResampleToRASInverse;
  MATRIX*              mResampleToSlice;
  MATRIX*              mResampleToSliceOrig;
  MATRIX*              mResampleToSliceInverse;


  VECTOR* mpTmpScreenIdx;     /* Used as tmp variables in macros. */
  VECTOR* mpTmpMRIIdx;
  xVoxel  mTmpVoxel;
  xVoxel  mTmpVoxel2;
  double  mTmpReal;

}
mriVolume, *mriVolumeRef;

Volm_tErr Volm_New        ( mriVolumeRef* opVolume );
Volm_tErr Volm_Delete     ( mriVolumeRef* iopVolume );
Volm_tErr Volm_DeepClone  ( mriVolumeRef  this,
                            mriVolumeRef* opVolume );

Volm_tErr Volm_CreateFromVolume( mriVolumeRef this,
                                 mriVolumeRef iVolume );
Volm_tErr Volm_ImportData      ( mriVolumeRef this,
                                 char*        isSource );
Volm_tErr Volm_ExportNormToCOR ( mriVolumeRef this,
                                 char*        isPath );

Volm_tErr Volm_Save ( mriVolumeRef iVolume,
                      char*        isFileName,
                      tBoolean     ibSaveFileName );

Volm_tErr Volm_LoadDisplayTransform ( mriVolumeRef this,
                                      char*        isSource );
Volm_tErr Volm_UnloadDisplayTransform ( mriVolumeRef this );


/* Copies the MRI's VOL_GEOM field. */
Volm_tErr Volm_CopyGeometryInformation ( mriVolumeRef this,
    VOL_GEOM*    ioVolumeGeometry );

/* Use the GetColor functions when you just need to display a color on
   the screen. It gets a sampled value and passes it through the color
   table to return an intensity color. */
void Volm_GetIntColorAtIdx        ( mriVolumeRef     this,
                                    xVoxelRef        iIdx,
                                    xColor3nRef      oColor );
void Volm_GetMaxIntColorAtIdx     ( mriVolumeRef     this,
                                    xVoxelRef        iIdx,
                                    mri_tOrientation iOrientation,
                                    xColor3nRef      oColor );

Volm_tErr Volm_GetDimensions        ( mriVolumeRef this,
                                      int*         onDimensionX,
                                      int*         onDimensionY,
                                      int*         onDimensionZ );
Volm_tErr Volm_GetNumberOfFrames    ( mriVolumeRef this,
                                      int*         onDimensionFrames );
Volm_tErr Volm_GetType              ( mriVolumeRef this,
                                      int*         onType );
Volm_tErr Volm_GetValueMinMax       ( mriVolumeRef this,
                                      float*       ofMin,
                                      float*       ofMax );

Volm_tErr Volm_GetResampleMethod    ( mriVolumeRef          this,
                                      Volm_tResampleMethod* oMethod );
Volm_tErr Volm_SetResampleMethod    ( mriVolumeRef          this,
                                      Volm_tResampleMethod  iMethod );

Volm_tErr Volm_GetSampleType  ( mriVolumeRef       this,
                                Volm_tSampleType*  oType );
Volm_tErr Volm_SetSampleType  ( mriVolumeRef       this,
                                Volm_tSampleType   iType );

/* Use the GetValue functions when you want the real value of the
   voxel. Before calling the unsafe version, make sure the volume is
   valid, the index is in bounds, and you're passing a valid
   pointer. */
Volm_tErr Volm_GetValueAtIdx       ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     float*       oValue );
Volm_tErr Volm_GetValueAtIdxUnsafe ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     float*       oValue );
Volm_tErr Volm_SetValueAtIdx       ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     float        iValue );

Volm_tErr Volm_GetValueAtIdxFrame       ( mriVolumeRef this,
    xVoxelRef    iIdx,
    int          iFrame,
    float*       oValue );
Volm_tErr Volm_GetValueAtIdxFrameUnsafe ( mriVolumeRef this,
    xVoxelRef    iIdx,
    int          iFrame,
    float*       oValue );
Volm_tErr Volm_SetValueAtIdxFrame       ( mriVolumeRef this,
    xVoxelRef    iIdx,
    int          iFrame,
    float        iValue );

/* These work on MRI indices instead of screen indices, for high-res
   volumes. */
Volm_tErr Volm_GetValueAtMRIIdx       ( mriVolumeRef this,
                                        xVoxelRef    iMRIIdx,
                                        float*       oValue );
Volm_tErr Volm_SetValueAtMRIIdx       ( mriVolumeRef this,
                                        xVoxelRef    iMRIIdx,
                                        float        iValue );

/* Returns whether or not we have a talairach transform. */
Volm_tErr Volm_HasTalTransform ( mriVolumeRef this,
                                 tBoolean*    obHasTransform );

/* coordinate conversion. idx stands for index and is the 0->dimension-1
   index of the volume. RAS space, aka world space, is centered on the
   center of the volume and is in mm. mni tal is mni's version of
   talairach space. the other tal is mni tal with a small modification. */
Volm_tErr Volm_ConvertIdxToRAS     ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     xVoxelRef    oRAS );
Volm_tErr Volm_ConvertRASToIdx     ( mriVolumeRef this,
                                     xVoxelRef    iRAS,
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
Volm_tErr Volm_ConvertMRIIdxToScanner ( mriVolumeRef this,
                                        xVoxelRef    iMRIIdx,
                                        xVoxelRef    oScanner );
Volm_tErr Volm_ConvertIdxToMRIIdx  ( mriVolumeRef this,
                                     xVoxelRef    iIdx,
                                     xVoxelRef    oMRIIdx );
Volm_tErr Volm_ConvertMRIIdxToIdx  ( mriVolumeRef this,
                                     xVoxelRef    iMRIIdx,
                                     xVoxelRef    oIdx );
Volm_tErr Volm_ConvertMRIIdxToRAS  ( mriVolumeRef this,
                                     xVoxelRef    iMRIIdx,
                                     xVoxelRef    oRAS );
Volm_tErr Volm_ConvertRASToMRIIdx  ( mriVolumeRef this,
                                     xVoxelRef    iRAS,
                                     xVoxelRef    oMRIIdx );

/* Use these when converting an RAS coming from the surface to an MRI
   index. This is different than normal RAS because surface RAS always
   has c_ras = 0, where volume RAS doesn't.  */
Volm_tErr Volm_ConvertMRIIdxToSurfaceRAS  ( mriVolumeRef this,
    xVoxelRef    iMRIIdx,
    xVoxelRef    oSurfaceRAS );
Volm_tErr Volm_ConvertSurfaceRASToMRIIdx  ( mriVolumeRef this,
    xVoxelRef    iSurfaceRAS,
    xVoxelRef    oMRIIdx );

Volm_tErr Volm_GetMRIIdxToAnaIdxTransform ( mriVolumeRef     this,
    mriTransformRef* opTransform );

/* Generic flood algorithm. Starts at a user-supplied voxel and floods
   outwards, in 3D or inplane, and for every valid voxel, calls the
   user-supplied visitation function. User needs to fill out
   Volm_tFloodParams struct. */
Volm_tErr Volm_Flood         ( mriVolumeRef        this,
                               Volm_tFloodParams*  iParams );
Volm_tErr Volm_FloodIterate_ ( mriVolumeRef        this,
                               Volm_tFloodParams*  iParams,
                               xVoxelRef           iIdx,
                               tBoolean*           visited );

/* calls the parameter function for every voxel, passing the MRI index, the
   value, and the pointer passed to it. */
Volm_tErr Volm_VisitAllVoxels ( mriVolumeRef        this,
                                Volm_tVisitFunction iFunc,
                                void*               ipData );


Volm_tErr Volm_FindMaxValues ( mriVolumeRef this );

Volm_tErr Volm_MakeColorTable ( mriVolumeRef this );

Volm_tErr Volm_SetBrightnessAndContrast ( mriVolumeRef this,
    float        ifBrightness,
    float        ifContrast );
Volm_tErr Volm_SetColorMinMax           ( mriVolumeRef this,
    float        ifMin,
    float        ifMax );
Volm_tErr Volm_GetBrightnessAndContrast ( mriVolumeRef this,
    float*       ofBrightness,
    float*       ofContrast );

Volm_tErr Volm_SaveToSnapshot      ( mriVolumeRef this );
Volm_tErr Volm_RestoreFromSnapshot ( mriVolumeRef this );

Volm_tErr Volm_Rotate    ( mriVolumeRef     this,
                           mri_tOrientation iAxis,
                           float            ifDegrees );
Volm_tErr Volm_Threshold ( mriVolumeRef     this,
                           float            iThreshold,
                           tBoolean         ibAbove,
                           float            iNewValue );
Volm_tErr Volm_Flip      ( mriVolumeRef     this,
                           mri_tOrientation iAxis );
Volm_tErr Volm_SetAllValues ( mriVolumeRef     this,
                              float            iValue );

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

/* Sets the mri voxel size and fov so that the smallest dimension is 1. */
Volm_tErr Volm_SetMinVoxelSizeToOne ( mriVolumeRef this );

/* Conforms the volume to 256^3 1mm isotropic. */
Volm_tErr Volm_Conform ( mriVolumeRef this );

/* Sets up internal data structures from an MRI volume and starts
   using the given MRI as its voxel values. */
Volm_tErr Volm_SetFromMRI_ ( mriVolumeRef this,
                             MRI*        iMRI );

/* Calculates this->mMRIIdxToAnaIdxTransform baesd on resample method. */
Volm_tErr Volm_CalculateMRIIdxToAnaIdx_ ( mriVolumeRef this );

/* For the max values in an orientation. */
Volm_tErr Volm_GetMaxValueAtMRIIdx_ ( mriVolumeRef     this,
                                      xVoxelRef        iMRIIdx,
                                      mri_tOrientation iOrientation,
                                      float*           oValue );

/* Note that the functions in this section are implemented as
   functions and macros. The functions are slower but safer, and the
   macros are faster but don't make any checks. So you should test
   with the functions and then build with the macros turned on. And
   keep the macro and functions version synced in
   development. _grin_ Turn the flag on at the top of this file. */
#ifndef VOLM_USE_MACROS

/* Note that all the public functions that refer to 'idx' refer to
   screen idx, or the index in 256^3 space. The local MRI idx might
   not be 256^3, so this conversion takes the idx from 256^3 space to
   the local MRI space. This is done using the m_resample matrix. All
   the access functions of course need to use the local MRI idx. */
void Volm_ConvertScreenIdxToMRIIdx_ ( mriVolumeRef this,
                                      xVoxelRef    iScreenIdx,
                                      xVoxelRef    oMRIIdx );
void Volm_ConvertMRIIdxToScreenIdx_ ( mriVolumeRef this,
                                      xVoxelRef    iMRIIdx,
                                      xVoxelRef    oScreenIdx );

/* safer function versions of main accessing and setting functions */
void Volm_GetValueAtIdx_             ( mriVolumeRef      this,
                                       xVoxelRef         iIdx,
                                       float*            oValue );
void Volm_GetValueAtMRIIdx_          ( mriVolumeRef      this,
                                       xVoxelRef         iMRIIdx,
                                       float*            oValue );
void Volm_GetValueAtIdxFrame_        ( mriVolumeRef      this,
                                       xVoxelRef         iIdx,
                                       int               iFrame,
                                       float*            oValue );
void Volm_GetSampledValueAtIdx_      ( mriVolumeRef      this,
                                       xVoxelRef         iIdx,
                                       float*            oValue);
void Volm_GetSampledValueAtIdxFrame_ ( mriVolumeRef      this,
                                       xVoxelRef         iIdx,
                                       int               iFrame,
                                       float*            oValue );


void Volm_SetValueAtIdx_         ( mriVolumeRef      this,
                                   xVoxelRef         iIdx,
                                   float             iValue );
void Volm_SetValueAtMRIIdx_         ( mriVolumeRef      this,
                                      xVoxelRef         iIdx,
                                      float             iValue );
void Volm_SetValueAtIdxFrame_    ( mriVolumeRef      this,
                                   xVoxelRef         iIdx,
                                   int               iFrame,
                                   float             iValue );

#else /* VOLM_USE_MACROS */

/* macro version  of main accessing and setting functions */

#define Volm_ConvertScreenIdxToMRIIdx_(this,iScreenIdx,oMRIIdx) \
  *MATRIX_RELT(this->mpTmpScreenIdx,1,1) = (iScreenIdx)->mfX; \
  *MATRIX_RELT(this->mpTmpScreenIdx,2,1) = (iScreenIdx)->mfY; \
  *MATRIX_RELT(this->mpTmpScreenIdx,3,1) = (iScreenIdx)->mfZ; \
  *MATRIX_RELT(this->mpTmpScreenIdx,4,1) = 1.0 ; \
  \
 MatrixMultiply( this->m_resample, this->mpTmpScreenIdx, this->mpTmpMRIIdx ); \
  \
  (oMRIIdx)->mfX = *MATRIX_RELT(this->mpTmpMRIIdx,1,1); \
  (oMRIIdx)->mfY = *MATRIX_RELT(this->mpTmpMRIIdx,1,2); \
  (oMRIIdx)->mfZ = *MATRIX_RELT(this->mpTmpMRIIdx,1,3); \
  \
  if( floor((oMRIIdx)->mfX+0.5) < 0 || \
floor((oMRIIdx)->mfX+0.5) >= this->mnDimensionX || \
      floor((oMRIIdx)->mfY+0.5) < 0 || \
floor((oMRIIdx)->mfY+0.5) >= this->mnDimensionY || \
      floor((oMRIIdx)->mfZ+0.5) < 0 || \
floor((oMRIIdx)->mfZ+0.5) >= this->mnDimensionZ ) { \
    xVoxl_Set( oMRIIdx, 0, 0, 0 ); \
  }

#define Volm_ConvertMRIIdxToScreenIdx_(this,iMRIIdx,oScreenIdx) \
  *MATRIX_RELT(this->mpTmpMRIIdx,1,1) = (iMRIIdx)->mfX; \
  *MATRIX_RELT(this->mpTmpMRIIdx,2,1) = (iMRIIdx)->mfY; \
  *MATRIX_RELT(this->mpTmpMRIIdx,3,1) = (iMRIIdx)->mfZ; \
  *MATRIX_RELT(this->mpTmpMRIIdx,4,1) = 1.0 ; \
 \
  MatrixMultiply( this->m_resample_inv,  \
    this->mpTmpMRIIdx, this->mpTmpScreenIdx ); \
 \
  (oScreenIdx)->mfX = *MATRIX_RELT(this->mpTmpScreenIdx,1,1); \
  (oScreenIdx)->mfY = *MATRIX_RELT(this->mpTmpScreenIdx,1,2); \
  (oScreenIdx)->mfZ = *MATRIX_RELT(this->mpTmpScreenIdx,1,3); \
 \
  if( (floor((oScreenIdx)->mfX+0.5) < 0 || \
floor((oScreenIdx)->mfX+0.5) >= 256 || \
       floor((oScreenIdx)->mfY+0.5) < 0 || \
floor((oScreenIdx)->mfY+0.5) >= 256 || \
       floor((oScreenIdx)->mfZ+0.5) < 0 || \
floor((oScreenIdx)->mfZ+0.5) >= 256) ) { \
 \
    xVoxl_Set( oScreenIdx, 0, 0, 0 ); \
  }



#define Volm_GetValueAtIdx_(this,iIdx,oValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  switch( this->mpMriValues->type ) { \
    case MRI_UCHAR: \
      *oValue =  \
 MRIvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
  (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ); \
      break; \
    case MRI_INT: \
      *oValue =  \
 MRIIvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
   (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ); \
      break; \
    case MRI_LONG: \
      *oValue =  \
 MRILvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
   (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ); \
      break; \
    case MRI_FLOAT: \
      *oValue =  \
 MRIFvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
   (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ); \
      break; \
    case MRI_SHORT: \
      *oValue =  \
 MRISvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
   (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ); \
      break; \
    default: \
      *oValue = 0; \
      break ; \
    }

#define Volm_GetValueAtMRIIdx_(this,iMRIIdx,oValue) \
  switch( this->mpMriValues->type ) { \
    case MRI_UCHAR: \
      *oValue =  \
 MRIvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
  (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ); \
      break; \
    case MRI_INT: \
      *oValue =  \
 MRIIvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
   (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ); \
      break; \
    case MRI_LONG: \
      *oValue =  \
 MRILvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
   (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ); \
      break; \
    case MRI_FLOAT: \
      *oValue =  \
 MRIFvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
   (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ); \
      break; \
    case MRI_SHORT: \
      *oValue =  \
 MRISvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
   (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ); \
      break; \
    default: \
      *oValue = 0; \
      break ; \
    }

#define Volm_GetSampledValueAtIdx_(this,iIdx,oValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  MRIsampleVolumeType( this->mpMriValues, \
         xVoxl_GetFloatX(&this->mTmpVoxel), \
         xVoxl_GetFloatY(&this->mTmpVoxel), \
         xVoxl_GetFloatZ(&this->mTmpVoxel), \
         &this->mTmpReal, \
         (int)this->mSampleType ); \
 \
  *oValue = (float)this->mTmpReal;



#define Volm_GetSampledValueAtIdxFrame_(this,iIdx,iFrame,oValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  MRIsampleVolumeFrameType( this->mpMriValues, \
       xVoxl_GetFloatX(&this->mTmpVoxel), \
       xVoxl_GetFloatY(&this->mTmpVoxel), \
       xVoxl_GetFloatZ(&this->mTmpVoxel), \
       iFrame, \
       (int)this->mSampleType, \
       &value ); \
 \
  *oValue = (float)this->mTmpReal;

#define Volm_SetValueAtIdx_(this,iIdx,iValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  switch (this->mpMriValues->type) \
    { \
    default: \
      break ; \
    case MRI_UCHAR: \
      MRIvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
       (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ) =  \
 (BUFTYPE) iValue; \
      break ; \
    case MRI_SHORT: \
      MRISvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
        (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ) =  \
 (short) iValue; \
      break ; \
    case MRI_FLOAT: \
      MRIFvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
        (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ) =  \
 (float) iValue; \
      break ; \
    case MRI_LONG: \
      MRILvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
        (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ) =  \
 (long) iValue; \
      break ; \
    case MRI_INT: \
      MRIIvox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
        (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5) ) =  \
 (int) iValue; \
      break ; \
    }

#define Volm_SetValueAtMRIIdx_(this,iMRIIdx,iValue) \
  switch (this->mpMriValues->type) \
    { \
    default: \
      break ; \
    case MRI_UCHAR: \
      MRIvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
       (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ) =  \
 (BUFTYPE) iValue; \
      break ; \
    case MRI_SHORT: \
      MRISvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
        (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ) =  \
 (short) iValue; \
      break ; \
    case MRI_FLOAT: \
      MRIFvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
        (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ) =  \
 (float) iValue; \
      break ; \
    case MRI_LONG: \
      MRILvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
        (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ) =  \
 (long) iValue; \
      break ; \
    case MRI_INT: \
      MRIIvox( this->mpMriValues, (int)floor((iMRIIdx)->mfX+0.5),  \
        (int)floor((iMRIIdx)->mfY+0.5), (int)floor((iMRIIdx)->mfZ+0.5) ) =  \
 (int) iValue; \
      break ; \
    }

#define Volm_GetValueAtIdxFrame_(this,iIdx,iFrame,oValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  switch( this->mpMriValues->type ) { \
    case MRI_UCHAR: \
      *oValue =  \
 MRIseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
         (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                    iFrame ); \
      break; \
    case MRI_INT: \
      *oValue =  \
 MRIIseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
       (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                     iFrame ); \
      break; \
    case MRI_LONG: \
      *oValue =  \
 MRILseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
       (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                     iFrame ); \
      break; \
    case MRI_FLOAT: \
      *oValue =  \
 MRIFseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
       (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                     iFrame ); \
      break; \
    case MRI_SHORT: \
      *oValue =  \
 MRISseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
       (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                     iFrame ); \
      break; \
    default: \
      *oValue = 0; \
      break ; \
    }

#define Volm_SetValueAtIdxFrame_(this,iIdx,iFrame,iValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  switch (this->mpMriValues->type) \
    { \
    default: \
      break ; \
    case MRI_UCHAR: \
      MRIseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
           (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                  iFrame ) =  \
 (BUFTYPE) iValue; \
      break ; \
    case MRI_SHORT: \
      MRISseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
            (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                   iFrame ) =  \
 (short) iValue; \
      break ; \
    case MRI_FLOAT: \
      MRIFseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
            (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                   iFrame ) =  \
 (float) iValue; \
      break ; \
    case MRI_LONG: \
      MRILseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
            (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                   iFrame ) =  \
 (long) iValue; \
      break ; \
    case MRI_INT: \
      MRIIseq_vox( this->mpMriValues, (int)floor((this->mTmpVoxel).mfX+0.5),  \
            (int)floor((this->mTmpVoxel).mfY+0.5), \
(int)floor((this->mTmpVoxel).mfZ+0.5), \
                   iFrame ) =  \
 (int) iValue; \
      break ; \
    }

#define Volm_GetSincValueAtIdx_(this,iIdx,irValue) \
  Volm_ConvertScreenIdxToMRIIdx_( this, iIdx, &this->mTmpVoxel ); \
 \
  MRIsincSampleVolume( this->mpMriValues,  \
         (this->mTmpVoxel).mfX, \
         (this->mTmpVoxel).mfY, \
         (this->mTmpVoxel).mfZ, \
         2, irValue);

#endif /* VOLM_USE_MACROS */

/* SOMEBODY changed this code so that the
   Volm_ConvertScreenIdxToMRIIdx_ function automatically sets
   out-of-bounds voxels to 0,0,0, so these functions are pretty much
   worthless because voxels always end up being valid. */
Volm_tErr Volm_Verify     ( mriVolumeRef this );
Volm_tErr Volm_VerifyIdx  ( mriVolumeRef this,
                            xVoxelRef    iIdx );
Volm_tErr Volm_VerifyIdx_ ( mriVolumeRef this,
                            xVoxelRef    iIdx );
Volm_tErr Volm_VerifyMRIIdx_ ( mriVolumeRef this,
                               xVoxelRef    iIdx );

/* So here is a function that really does return an error when an
   index is out of bounds. This can be used to actually verify indices
   so you don't end up drawing the value at 0,0,0 when you're outside
   the bounds of the volume. -RKT */
Volm_tErr Volm_VerifyIdxInMRIBounds ( mriVolumeRef this,
                                      xVoxelRef    iIdx );

Volm_tErr Volm_VerifyFrame_ ( mriVolumeRef this,
                              int          iFrame );
char* Volm_GetErrorString ( Volm_tErr ieCode );


#endif
