#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>

/* matrices and vectors are the same data type, vector will only 
   have one column (by default) or one row.
 */
typedef struct
{
  int   type ;
  int   rows ;
  int   cols ;
  float **rptr;    /* pointer to an array of rows */
  float *data;     /* pointer to base of data */
} MATRIX, VECTOR ;

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
int     MatrixPrintTranspose(FILE *fp, MATRIX *mat) ;
MATRIX  *MatrixTranspose(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixAdd(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixSubtract(MATRIX *m1, MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixScalarMul(MATRIX *mIn, float val, MATRIX *mOut) ;
MATRIX  *MatrixClear(MATRIX *mat) ;
MATRIX  *MatrixSquareElts(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixDiag(MATRIX *mDiag, MATRIX *mOut) ;
MATRIX  *MatrixMakeDiagonal(MATRIX *mSrc, MATRIX *mDst) ;
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
MATRIX *MatrixRegularize(MATRIX *mIn, MATRIX *mOut) ;
int    MatrixSingular(MATRIX *m) ;

/* determinants and eigenvectors */
float  MatrixDeterminant(MATRIX *m) ;
MATRIX *MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_dst) ;
MATRIX *MatrixSVD(MATRIX *mA, VECTOR *v_z, MATRIX *mV) ;
MATRIX *MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse) ;
float  MatrixConditionNumber(MATRIX *m) ;

/* statistical stuff */
MATRIX *MatrixCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans) ;
MATRIX *MatrixUpdateCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans);
MATRIX *MatrixUpdateMeans(MATRIX *mInputs, VECTOR *mMeans, VECTOR *mNobs) ;
MATRIX *MatrixFinalMeans(VECTOR *mMeans, VECTOR *mNobs) ;
MATRIX *MatrixFinalCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mNobs);

/* misc. I/O functions */
int    MatrixAsciiWriteInto(FILE *fp, MATRIX *m) ;
MATRIX *MatrixAsciiReadFrom(FILE *fp, MATRIX *m) ;

#define VectorAlloc(n, type)       MatrixAlloc(n, 1, type)
#define RVectorAlloc(n, type)      MatrixAlloc(1, n, type)
#define VectorFree(pm)             MatrixFree(pm)
#define VectorAdd(v1, v2, v3)      MatrixAdd(v1, v2, v3)
#define VectorSubtract(v1,v2,v3)   MatrixSubtract(v1, v2, v3)
#define VectorScalarMul(v1,val,v2) MatrixScalarMul(v1, val, v2)
#define VectorCopy(v1, v2)         MatrixCopy(v1, v2)
#define VectorClear(v)             MatrixClear(v)
#define VectorTranspose(vsrc,vdst) MatrixTranspose(vsrc, vdst)
#define VectorAsciiWriteInto       MatrixAsciiWriteInto
#define VectorAsciiReadFrom        MatrixAsciiReadFrom

#define VECTOR_ELT(v,i)            ((v)->rptr[i][1])
#define RVECTOR_ELT(v,i)            ((v)->rptr[1][i])
#define VECTOR3_LOAD(v,x,y,z)    (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z) ;
#define VECTOR_LOAD   VECTOR3_LOAD

float  VectorLen(VECTOR *v) ;
float  VectorDot(VECTOR *v1, VECTOR *v2) ;
float  VectorNormalizedDot(VECTOR *v1, VECTOR *v2) ;
float  VectorDistance(VECTOR *v1, VECTOR *v2) ;
VECTOR *MatrixColumn(MATRIX *m, VECTOR *v, int col) ;
MATRIX *VectorOuterProduct(VECTOR *v1, VECTOR *v2, MATRIX *m) ;
VECTOR *VectorCrossProduct(VECTOR *v1, VECTOR *v2, VECTOR *vdst) ;
float  VectorTripleProduct(VECTOR *v1, VECTOR *v2, VECTOR *v3) ;
VECTOR *VectorNormalize(VECTOR *vin, VECTOR *vout) ;


/* these are macro versions that work on 3-d vectors */
#define V3_X(v)      (VECTOR_ELT(v,1))
#define V3_Y(v)      (VECTOR_ELT(v,2))
#define V3_Z(v)      (VECTOR_ELT(v,3))
#define V3_LEN(v)    (sqrt(V3_X(v)*V3_X(v)+V3_Y(v)*V3_Y(v)+V3_Z(v)*V3_Z(v)))
#define V3_LEN_SQ(v) (V3_X(v)*V3_X(v)+V3_Y(v)*V3_Y(v)+V3_Z(v)*V3_Z(v))
#define V3_CROSS_PRODUCT(va,vb,vc) \
                 V3_X(vc) = V3_Y(va)*V3_Z(vb)- V3_Z(va)*V3_Y(vb),  \
                 V3_Y(vc) = V3_Z(va)*V3_X(vb)- V3_X(va)*V3_Z(vb),  \
                 V3_Z(vc) = V3_X(va)*V3_Y(vb)- V3_Y(va)*V3_X(vb) ; 
#define V3_TRIPLE_PRODUCT(va,vb,vc) \
                (V3_X(vc) * (V3_Y(va)*V3_Z(vb)- V3_Z(va)*V3_Y(vb)) + \
                 V3_Y(vc) * (V3_Z(va)*V3_X(vb)- V3_X(va)*V3_Z(vb)) + \
                 V3_Z(vc) * (V3_X(va)*V3_Y(vb)- V3_Y(va)*V3_X(vb)))
#define V3_ADD(va,vb,vc) \
                 V3_X(vc) = V3_X(va)+V3_X(vb), \
                 V3_Y(vc) = V3_Y(va)+V3_Y(vb), \
                 V3_Z(vc) = V3_Z(va)+V3_Z(vb) ;
#define V3_DOT(va,vb) (V3_X(va)*V3_X(vb)+V3_Y(va)*V3_Y(vb)+V3_Z(va)*V3_Z(vb))
#define V3_SCALAR_MUL(va,s,vb)  (V3_X(vb)=V3_X(va)*s,\
                                 V3_Y(vb)=V3_Y(va)*s,\
                                 V3_Z(vb)=V3_Z(va)*s)

#define V3_NORMALIZE(va,vb)  { float len = (V3_LEN(va)) ; \
                                  if (FZERO(len)) len = 1.0f ; \
                                  else len = 1.0f / len ; \
                                  V3_SCALAR_MUL(va,len,vb) ; }
#define V3_CLEAR(v)          (V3_X(v) = 0, V3_Y(v) = 0, V3_Z(v) = 0)

#define X_ROTATION   0
#define Y_ROTATION   1
#define Z_ROTATION   2

MATRIX *MatrixAllocRotation(int n, float angle, int which) ;
#define MatrixClone(mat)   MatrixCopy(mat, NULL)
#define VectorClone        MatrixClone



#endif



