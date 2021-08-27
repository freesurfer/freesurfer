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
#include "mri_topology.h"

#ifndef MRICURV_INCLUDED
#define MRICURV_INCLUDED

float MRIcurvature(MRI *mri,int i,int j,int k,int label,int connectivity);
float Nbhcurvature(Nbh *nbh,int connectivity);


#endif





