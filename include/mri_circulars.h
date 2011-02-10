/**
 * @file  mri_circulars.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/02/10 19:17:32 $
 *    $Revision: 1.4 $
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
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
