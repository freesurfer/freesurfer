#ifndef MRI_CIRCULARS_H
#define MRI_CIRCULARS_H

// These functions are declared here, separately from mri.h and mrisurf.h,
// because they require headers files which create circular dependencies in
// some files when compiled.
// This file should be included last in a header list where these functions
// are needed.

#include "image.h" // IMAGE
IMAGE *MRItoImage(MRI *mri, 
                  IMAGE *image, 
                  int slice) ;
MRI *ImageToMRI(IMAGE *image);
IMAGE *MRItoImageView(MRI *mri, 
                      IMAGE *image, 
                      int slice, 
                      int view, 
                      int frame) ;

#include "mrisurf.h" // MRIS
int MRIwriteAnyFormat(MRI *mri, 
                      char *fileid, 
                      char *fmt,
                      int mriframe, 
                      MRIS *surf);

#include "label.h" // LABEL
HISTOGRAM  *MRIhistogramLabelStruct(MRI *mri, 
                                    int nbins, 
                                    HISTOGRAM *histo,
                                    LABEL *label) ;

#include "transform.h" // VOL_GEOM
int MRIcopyVolGeomToMRI(MRI *mri, VOL_GEOM *vg) ;
int MRIcopyVolGeomFromMRI(MRI *mri, VOL_GEOM *vg) ;

#endif // MRI_CIRCULARS_H
