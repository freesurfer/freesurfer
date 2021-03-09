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


#ifndef SIG_H
#define SIG_H

double fdr2vwth(double *p, int np, double fdr);
double vwth2fdr(double *p, int np, double vwth);
int doublecompar(const void *v1, const void *v2);
double sigt(double t, int df) ;
float sigchisq(double chisq, int df) ;

#endif
