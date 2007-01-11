/**
 * @file  mriTypes.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/01/11 20:15:14 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2002-2007, CorTechs Labs, Inc. (La Jolla, CA) and
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef mriTypes_h
#define mriTypes_h

#define mri_knSubjectNameLen 256
#define mri_knPathLen 1024

typedef enum {
  mri_tOrientation_None = -1,
  mri_tOrientation_Coronal = 0,
  mri_tOrientation_Horizontal,
  mri_tOrientation_Sagittal,
  mri_knNumOrientations
} mri_tOrientation;

typedef enum {
  mri_tCoordSpace_None = -1,
  mri_tCoordSpace_VolumeIdx = 0,
  mri_tCoordSpace_SurfaceRAS,
  mri_tCoordSpace_RAS,
  mri_tCoordSpace_Talairach,
  mri_knNumCoordSpaces
} mri_tCoordSpace;

#endif
