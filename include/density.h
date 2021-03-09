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


#ifndef DENSITY_H
#define DENSITY_H

#include <stdio.h>
#include "image.h"
#include "const.h"
#include "mri.h"
#include "matrix.h"

typedef struct
{
  float   min_val1 ;
  float   min_val2 ;
  float   max_val1 ;
  float   max_val2 ;
  int     *valid1 ;
  int     *valid2 ;
  char    fname1[STRLEN] ;
  char    fname2[STRLEN] ;
  float   sigma ;
  int     dof ;
  float   min_p ;
  IMAGE   *Ipdf ;
}
DENSITY ;

DENSITY *DensityHistogramEstimate(MRI *mri1, MRI *mri2, int nbins, float sigma, int *valid1, int *valid2) ;
int      DensityWrite(DENSITY *pdf, char *fname) ;
double   DensityLogLikelihood(DENSITY *pdf, float val1, float val2) ;
DENSITY *DensityRead(char *fname) ;
MRI      *DensityLikelihoodImage(MRI *mri1, MRI *mri2, MRI *mri_ll, MATRIX *m, DENSITY *pdf, MRI *mri_seg, int inverse);

#define val2bin1(pdf,val1)   ((val1-pdf->min_val1) * (pdf->Ipdf->rows /  (pdf->max_val1 - pdf->min_val1)))
#define val2bin2(pdf,val2)   ((val2-pdf->min_val2) * (pdf->Ipdf->rows /  (pdf->max_val2 - pdf->min_val2)))

#endif
