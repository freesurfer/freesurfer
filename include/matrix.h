#ifndef MATRIX_H
#define MATRIX_H

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

extern MATRIX  *MatrixInverse(MATRIX *mIn, MATRIX *mOut) ;
extern MATRIX  *MatrixAlloc(int rows, int cols, int type) ;
extern int     MatrixFree(MATRIX **pmat) ;
extern MATRIX  *MatrixMultiply(MATRIX *m1, MATRIX *m2, MATRIX *m3) ;
extern MATRIX  *MatrixCopy(MATRIX *mIn, MATRIX *mOut) ;
extern MATRIX  *MatrixRead(char *fname, MATRIX *mOut) ;
extern MATRIX  *MatrixIdentity(int n, MATRIX *mI) ;
extern void    MatrixPrint(FILE *fp, MATRIX *mat) ;
extern MATRIX  *MatrixTranspose(MATRIX *mIn, MATRIX *mOut) ;
extern MATRIX  *MatrixAdd(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
extern MATRIX  *MatrixSubtract(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
extern MATRIX  *MatrixScalarMul(MATRIX *mIn, float val, MATRIX *mOut) ;
extern MATRIX  *MatrixClear(MATRIX *mat) ;
extern MATRIX  *MatrixSquareElts(MATRIX *mIn, MATRIX *mOut) ;
extern MATRIX  *MatrixDiag(MATRIX *mDiag, MATRIX *mOut) ;
extern MATRIX  *MatrixCopyRegion(MATRIX *mSrc, MATRIX *mDst, int start_row, 
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
extern MATRIX  *MatrixCopyRealRegion(MATRIX *mSrc, MATRIX *mDst,int start_row,
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
extern MATRIX  *MatrixCopyImagRegion(MATRIX *mSrc, MATRIX *mDst, int start_row,
                                 int start_col, int rows, int cols, 
                                 int dest_row, int dest_col) ;
extern MATRIX *MatrixRealToComplex(MATRIX *mReal, MATRIX *mImag, MATRIX *mOut);

#define MatrixClone(mat)   MatrixCopy(mat, NULL)

#endif



