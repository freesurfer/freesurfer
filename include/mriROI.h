/**
 * @file  mriROI.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef mriROI_h
#define mriROI_h

#include "xTypes.h"
#include "xDebug.h"
#include "xVoxel.h"
#include "xList.h"

typedef enum
{

  ROI_tErr_NoErr = 0,
  ROI_tErr_InvalidObject,
  ROI_tErr_InvalidParameter,
  ROI_tErr_InvalidSignature,
  ROI_tErr_AllocationFailed,
  ROI_tErr_ErrorAccessingList,
  ROI_tErr_LastVoxel,
  ROI_tErr_InvalidErrorCode,
  ROI_knNumErrorCodes

} ROI_tErr;


typedef enum
{

  ROI_tIterationPlane_X = 0,
  ROI_tIterationPlane_Y,
  ROI_tIterationPlane_Z,
  ROI_tIterationPlane_All,
  ROI_knNumIterationPlanes
} ROI_tIterationPlane;

/* voxel value pair */
typedef struct
{
  xVoxel  mVoxel;
  float*  mafValues;
}
ROI_tROIVoxel, *ROI_tROIVoxelRef;

#define ROI_kSignature 0x01928374

typedef struct
{

  long mSignature;

  /* list of ROI values */
  xListRef   mVoxels;

  /* index and name of the ROI */
  int        mnIndex;
  char       msName[256];

  /* number of values per voxel and their labels */
  int        mnNumValuesPerVoxel;
  char**     masValueLabels;

  /* iterator state */
  ROI_tIterationPlane mIterPlane;
  float               mfIterPlaneNumber;
  float               mfIterPlaneRange;
  ROI_tROIVoxelRef    mIterVoxel;

}
mriROI, *mriROIRef;


ROI_tErr ROI_New    ( mriROIRef* opROI,
                      int        inIndex,
                      char*      isName,
                      int        inNumValuesPerVoxel,
                      char**     iasValueLabels );
ROI_tErr ROI_Delete ( mriROIRef* iopROI );

ROI_tErr ROI_AddVoxel ( mriROIRef this,
                        xVoxelRef  iVoxel,
                        float*     iafValues );

ROI_tErr ROI_GetIndex ( mriROIRef this,
                        int*      onIndex );
ROI_tErr ROI_GetName  ( mriROIRef this,
                        char*     osName );

ROI_tErr ROI_DebugPrint ( mriROIRef this );

ROI_tErr ROI_IsVoxelInROI ( mriROIRef this,
                            xVoxelRef iVoxel,
                            tBoolean* oIsInROI );

/* iterate thru points */
ROI_tErr ROI_ResetIterator ( mriROIRef            this,
                             ROI_tIterationPlane  iPlane,
                             float                ifPlaneNumber,
                             float                ifPlaneRange );
ROI_tErr ROI_NextVoxel     ( mriROIRef         this,
                             ROI_tROIVoxelRef* opPoint );

char* ROI_GetErrorString ( ROI_tErr ieCode );
ROI_tErr ROI_Verify ( mriROIRef this );

#endif


