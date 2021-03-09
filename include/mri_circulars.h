/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
                      const char *fileid,
                      const char *fmt,
                      int mriframe,
                      MRIS *surf);

#include "label.h" // LABEL
HISTOGRAM  *MRIhistogramLabelStruct(MRI *mri,
                                    int nbins,
                                    HISTOGRAM *histo,
                                    LABEL *label) ;

#include "transform.h" // VOL_GEOM
int MRIcopyVolGeomToMRI( MRI *mri, const VOL_GEOM *vg ) ;
int MRIcopyVolGeomFromMRI( const MRI *mri, VOL_GEOM *vg ) ;

#endif // MRI_CIRCULARS_H
