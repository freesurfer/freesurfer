#ifndef mriROI_h
#define mriROI_h

#include "xVoxel.h"

typedef enum {

  ROI_tErr_NoErr = 0,
  ROI_tErr_InvalidObject,
  ROI_tErr_InvalidParameter,
  ROI_tErr_InvalidSignature,
  ROI_tErr_AllocationFailed,
  ROI_tErr_InvalidErrorCode,
  ROI_knNumErrorCodes

} ROI_tErr;


typedef enum {

  ROI_tIterationPlane_X = 0,
  ROI_tIterationPlane_Y,
  ROI_tIterationPlane_Z,
  ROI_tIterationPlane_All,
  ROI_knNumIterationPlanes
} ROI_tIterationPlane;

/* voxel value pair */
typedef struct {
  xVoxel mVoxel;
  float  mfValue;
} ROI_tROIValue, *ROI_tROIValueRef;


typedef struct {

  long mSignature;

  /* list of ROI values */
  xListRef   mROIValues;

  /* name of the ROI */
  char       msName[256];

} mriROI, *mriROIRef;


ROI_tErr ROI_New    ( mriROIRef* oList,
      char*                isName );
ROI_tErr ROI_Delete ( mriROIRef* opList );

ROI_tErr ROI_Add ( mriROIRef* opList,
         xVoxelRef            iVoxel,
         float                ifValue );

/* iterate thru points */
ROI_tErr ROI_ResetIterator ( mriROIRef            this,
           ROI_tIterationPlane  iPlane,
           float                ifPlaneNumber,
           float                ifPlaneRange );
ROI_tErr ROI_NextVoxel     ( mriROIRef         this,
           ROI_tROIValueRef* opPoint );




#endif


