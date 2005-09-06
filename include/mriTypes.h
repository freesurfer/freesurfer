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
  mri_knNumCoordSpaces,
	mri_tCoordSpace_SurfaceRAS
} mri_tCoordSpace;

#endif
