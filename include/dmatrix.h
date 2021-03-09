/**
 * @brief matrix manipulation utils using doubles instead of floats
 *
 */
/*
 * Original Author: Doug Greve
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


#ifndef DMATRIX_H
#define DMATRIX_H

#ifdef X
#undef X
#endif

#include "matrix.h"

#define DMATRIX_TYPE double

// Same as in matrix.h but using double instead of float
typedef struct
{
  short   type ;
  int     rows ;
  int     cols ;
  DMATRIX_TYPE **rptr;    /* pointer to an array of rows */
  DMATRIX_TYPE *data;     /* pointer to base of data */
}
DMATRIX, DVECTOR ;

int DMatrixTest(void);
int DMatrixCheckDims(const DMATRIX *m1, const DMATRIX *m2, const int checktype, FILE *fp, const char *str);
DMATRIX *DMatrixAlloc(const int rows, const int cols, const int type);
int DMatrixFree(DMATRIX **pmat);
int DMatrixPrintFmt(FILE *fp, const char *fmt, DMATRIX *mat);
int DMatrixPrint(FILE *fp, DMATRIX *mat);
double DMatrixCompareFMatrix(MATRIX *mf, DMATRIX *md);
DMATRIX *DMatrixCopyFMatrix(MATRIX *mf, DMATRIX *md);
DMATRIX *DMatrixMultiply(const DMATRIX *m1, const DMATRIX *m2, DMATRIX *m3);
DMATRIX *DMatrixCopy(DMATRIX *msrc, DMATRIX *mcopy);
DMATRIX *DMatrixTranspose(DMATRIX *mIn, DMATRIX *mOut);
DMATRIX *DMatrixAddMul(DMATRIX *m1, DMATRIX *m2, double v1, double v2, DMATRIX *mOut);
DMATRIX *DMatrixAdd(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout);
DMATRIX *DMatrixSubtract(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout);
DMATRIX *DMatrixScalarMul(DMATRIX *m, double v, DMATRIX *mout);
double DVectorDot(const DVECTOR *v1, const DVECTOR *v2);
DMATRIX *DMatrixConstVal(const double val, const int rows, const int cols, DMATRIX *X);
DMATRIX *DMatrixZero(const int rows, const int cols, DMATRIX *X);
double DVectorLen(const DVECTOR *v);
DVECTOR *DVectorCrossProduct(const DVECTOR *v1, const DVECTOR *v2, DVECTOR *vdst);
DMATRIX_TYPE DMatrixMaxAbs(DMATRIX *M);
DMATRIX *DMatrixDRand48(int rows, int cols, DMATRIX *m);

#endif
