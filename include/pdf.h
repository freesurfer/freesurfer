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

#ifndef PDF_INC
#define PDF_INC

unsigned long PDFtodSeed(void);
double PDFgaussian(void);
double PDFerlang(int order);
int PDFloadCDF(char *fname, double **xcdf, double **cdf, int *ncdf);
double PDFsampleCDF(double *xcdf, double *cdf, int ncdf);
int PDFsearchOrderedTable(double u, double *y, int ny);

#endif //#ifndef PDF_INC
