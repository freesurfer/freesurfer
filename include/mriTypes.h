/**
 * @file  mriTypes.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/02 19:25:19 $
 *    $Revision: 1.8 $
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


#ifndef mriTypes_h
#define mriTypes_h

#define mri_knSubjectNameLen 256
#define mri_knPathLen 1024

typedef enum
{
  mri_tOrientation_None = -1,
  mri_tOrientation_Coronal = 0,
  mri_tOrientation_Horizontal,
  mri_tOrientation_Sagittal,
  mri_knNumOrientations
} mri_tOrientation;

typedef enum
{
  mri_tCoordSpace_None = -1,
  mri_tCoordSpace_VolumeIdx = 0,
  mri_tCoordSpace_SurfaceRAS,
  mri_tCoordSpace_RAS,
  mri_tCoordSpace_Talairach,
  mri_knNumCoordSpaces
} mri_tCoordSpace;

#endif
