#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

typedef struct
{
  int   type ;
  int   rows ;
  int   cols ;
  float **rptr;    /* pointer to an array of rows */
  float *data;     /* pointer to base of data */
} MATRIX ;

typedef struct
{
  float  real ;
  float  imag ;
} COMPLEX_FLOAT, *CPTR ;

#define MATRIX_CELT(m,r,c)      (((COMPLEX_FLOAT **)m->rptr)[r]+c)
#define MATRIX_RELT(m,r,c)      (m->rptr[r]+c)
#define MATRIX_ELT(m,r,c)       (m->type == MATRIX_REAL ? \
                                  *MATRIX_RELT(m,r,c) : \
                                  *MATRIX_CELT(m,r,c))
#define MATRIX_PTR(m,r,c)       (m->type == MATRIX_REAL ? \
                                  MATRIX_RELT(m,r,c) : \
                                  (float *)MATRIX_CELT(m,r,c))

#define MATRIX_CELT_REAL(m,r,c)  (MATRIX_CELT(m,r,c)->real)
#define MATRIX_CELT_IMAG(m,r,c)  (MATRIX_CELT(m,r,c)->imag)

#define MATRIX_REAL        1
#define MATRIX_COMPLEX     2

MATRIX  *MatrixInverse(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixAlloc(int rows, int cols, int type) ;
int     MatrixFree(MATRIX **pmat) ;
MATRIX  *MatrixMultiply(MATRIX *m1, MATRIX *m2, MATRIX *m3) ;
MATRIX  *MatrixCopy(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixReadTxt(char *fname, MATRIX *mat) ;
MATRIX  *MatrixRead(char *fname) ;
int     MatrixWrite(MATRIX *mIn, char *fname, char *name) ;
MATRIX  *MatrixIdentity(int n, MATRIX *mI) ;
int     MatrixPrint(FILE *fp, MATRIX *mat) ;
MATRIX  *MatrixTranspose(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixAdd(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixSubtract(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixScalarMul(MATRIX *mIn, float val, MATRIX *mOut) ;
MATRIX  *MatrixClear(MATRIX *mat) ;
MATRIX  *MatrixSquareElts(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixDiag(MATRIX *mDiag, MATRIX *mOut) ;
MATRIX  *MatrixCopyRegion(MATRIX *mSrc, MATRIX *mDst, int start_row, 
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
MATRIX  *MatrixCopyRealRegion(MATRIX *mSrc, MATRIX *mDst,int start_row,
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
MATRIX  *MatrixCopyImagRegion(MATRIX *mSrc, MATRIX *mDst, int start_row,
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
MATRIX *MatrixRealToComplex(MATRIX *mReal, MATRIX *mImag, MATRIX *mOut);
float  MatrixDeterminant(MATRIX *m) ;
MATRIX *MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_dst) ;
MATRIX *MatrixSVD(MATRIX *mA, float *z, MATRIX *mV) ;

#define VectorAlloc(n, type)       MatrixAlloc(n, 1, type)
#define VectorFree(pm)             MatrixFree(pm)

#define X_ROTATION   0
#define Y_ROTATION   1
#define Z_ROTATION   2

MATRIX *MatrixAllocRotation(int n, float angle, int which) ;
#define MatrixClone(mat)   MatrixCopy(mat, NULL)

#endif



