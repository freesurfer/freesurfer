/**
 * @file  mgh_matrix.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.6 $
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


#ifndef MGH_MATRIX_H
#define MGH_MATRIX_H

#define TINY 1.0e-20

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

FLOATTYPE* MGH_vector(int n);
int* MGH_ivector(int n);
FLOATTYPE ** MGH_matrix(int n,int m);
void print_matrix(FLOATTYPE **a,int m,int n);
void read_matrix(FILE *fptr,FLOATTYPE **a,int m,int n);
void print_vector(FLOATTYPE *v,int n);
void row_vector(FLOATTYPE **a,FLOATTYPE *v,int i,int n);
void vector_to_matrix(FLOATTYPE *v,FLOATTYPE **a,int m,int n);
void scale_matrix(FLOATTYPE **a,FLOATTYPE s,int n,int m);
void normalize_matrix(FLOATTYPE **a,int n,int m);
void matrix_copy(FLOATTYPE **a,FLOATTYPE **b,int n,int m);
void matrix_copy2(FLOATTYPE **a,FLOATTYPE **b,int n,int m,int sno,int smo,
                  int tno,int tmo);
void matrix_transpose(FLOATTYPE **a,FLOATTYPE **at,int n,int m);
void matrix_add(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c, int n,int m);
void matrix_multiply(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m);
void matrix_multiply2(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m,
                      int l);
void matrix_angles(FLOATTYPE **a,FLOATTYPE **b,FLOATTYPE **c,int n,int m);
void vector_subtract(FLOATTYPE *a,FLOATTYPE *b,FLOATTYPE *c,int n);
void vector_add(FLOATTYPE *a,FLOATTYPE *b,FLOATTYPE *c,FLOATTYPE fa,
                FLOATTYPE fb,int n);
void vector_multiply(FLOATTYPE **a,FLOATTYPE *b,FLOATTYPE *c,int n,int m);
void derivative_matrix(FLOATTYPE **D,int dy,int dx);
void MGH_identity_matrix(FLOATTYPE **I,int n);
void regularization_matrix(FLOATTYPE **R,int n);
void covariance_matrix(FLOATTYPE **R,int n);
void MGH_nrerror(char *s);
void mgh_ludcmp(FLOATTYPE **a,int n,int *indx,FLOATTYPE *d);
void mgh_lubksb(FLOATTYPE **a,int n,int *indx,FLOATTYPE *b);
void inverse(FLOATTYPE **a,FLOATTYPE **y,int n);
FLOATTYPE determinant(FLOATTYPE **a,int n);
FLOATTYPE MGH_svd(FLOATTYPE **A,FLOATTYPE **V,FLOATTYPE *z,int m,int n);
void mgh_svdcmp(FLOATTYPE **a,FLOATTYPE *w,FLOATTYPE **v,int m,int n);


#endif
