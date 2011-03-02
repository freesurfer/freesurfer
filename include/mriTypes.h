/**
 * @file  mriTypes.h
 * @brief general purpose utils
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.9 $
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
