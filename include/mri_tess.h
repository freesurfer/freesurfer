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


#include "mri.h"
#include "mrisurf.h"

#ifndef TESS_INCLUDED
#define TESS_INCLUDED

MRIS* MRIScreateSurfaceFromVolume(MRI *mri,int label,int connectivity);
MRIS** MRIScreateSurfacesFromVolume(MRI *mri,int number_of_labels, int* labelvalue,int connectivity);

#endif
