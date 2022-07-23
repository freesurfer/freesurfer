/**
 * @brief matrix manipulation utils
 *
 */
/*
 * Original Author: Bruce Fischl
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


#ifndef MATRIX_H
#define MATRIX_H

#ifdef X
#undef X
#endif

#include "base.h"

/* matrices and vectors are the same data type, vector will only
   have one column (by default) or one row.
*/
typedef struct
{
  short   type ;
  char    inBuf;	// The matrix is in a stack buffer (see below)
  int     rows ;
  int     cols ;
  float **rptr;    /* pointer to an array of rows */
  float *data;     /* pointer to base of data */
  FILE *mmapfile;
}
MATRIX, VECTOR ;


typedef struct MatrixBuffer {
  MATRIX  matrix;
  float * rptr[5];
  float   data[5*5+2];
} MatrixBuffer;		// Upto 4x4 are so common and so small they should be stack allocated
			// so MatrixFree will just NULL the pointer

typedef struct		// This case is so important it should be optimized 
{
  float x,y,z;
} XYZ;


typedef struct
{
  float  real ;
  float  imag ;
}
COMPLEX_FLOAT, *CPTR ;

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
#define MATRIX_SYM    0
#define MATRIX_UPPER  1
#define MATRIX_LOWER  2

int     MatrixIsZero(MATRIX *m) ;
int     MatrixIsIdentity(MATRIX *m) ;
MATRIX  *MatrixReshape(MATRIX *m_src, MATRIX *m_dst, int rows, int cols) ;
int     MatrixCheck(MATRIX *m) ;
int     MatrixQRdecomposition(const MATRIX *iMatrix, MATRIX *oQ, MATRIX *oR);
MATRIX  *MatrixInverse( const MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv) ;
MATRIX  *MatrixSVDPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv) ;
MATRIX  *MatrixRightPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv) ;
#define MatrixLeftPseudoInverse MatrixPseudoInverse

MATRIX  *MatrixAlloc_wkr ( const int rows, const int cols, const int type,
    	                    const char* callSiteFile, int callSiteLine);
MATRIX  *MatrixAlloc2_wkr( const int rows, const int cols, const int type, MatrixBuffer* buf, 
    	    	    	    const char* callSiteFile, int callSiteLine);
#define MatrixAlloc(ROWS,COLS,TYPE)      MatrixAlloc_wkr ((ROWS),(COLS),(TYPE),       __FILE__,__LINE__)
#define MatrixAlloc2(ROWS,COLS,TYPE,BUF) MatrixAlloc2_wkr((ROWS),(COLS),(TYPE),(BUF), __FILE__,__LINE__)

int     MatrixFree(MATRIX **pmat) ;
MATRIX  *MatrixMultiplyD( const MATRIX *m1, const MATRIX *m2, MATRIX *m3); // use this one
MATRIX  *MatrixMultiply_wkr( const MATRIX *m1, const MATRIX *m2, MATRIX *m3, 
    	    	    	    const char* callSiteFile, int callSiteLine) ;
#define MatrixMultiply(M1,M2,M3) MatrixMultiply_wkr((M1),(M2),(M3),__FILE__,__LINE__)
    	// (r1 x c1) * (r2 x c2) = (r1 x c2)
	// c1 must equal r2

MATRIX *MatrixMultiplyElts(MATRIX *m1, MATRIX *m2, MATRIX *m12); // like matlab m1.*m2
MATRIX *MatrixDivideElts(MATRIX *num, MATRIX *den, MATRIX *quotient); // like matlab num./den
MATRIX *MatrixReplicate(MATRIX *mIn, int nr, int nc, MATRIX *mOut); // like matlab repmat()
MATRIX  *MatrixCopy( const MATRIX *mIn, MATRIX *mOut );
int     MatrixWriteTxt(const char *fname, MATRIX *mat) ;
MATRIX  *MatrixReadTxt(const char *fname, MATRIX *mat) ;
MATRIX  *MatrixRead(const char *fname) ;
int     MatrixWrite(MATRIX *mIn,const char *fname, const char *name) ;
MATRIX  *MatrixIdentity(int n, MATRIX *mI) ;
int     MatrixPrint(FILE *fp, const MATRIX *mat) ;
int     MatrixPrintFmt(FILE *fp,const char *fmt, MATRIX *mat);
int     MatrixPrintOneLine(FILE *fp, MATRIX *mat) ;
int     MatrixPrintTranspose(FILE *fp, MATRIX *mat) ;
int     MatrixPrintWithString(FILE *fp, MATRIX *m, const char *Pre, const char *Post);
MATRIX  *MatrixTranspose(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixAdd( const MATRIX *m1, const MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixSubtract( const MATRIX *m1, const MATRIX *m2, MATRIX *mOut) ;
MATRIX  *MatrixScalarMul( const MATRIX *mIn, const float val, MATRIX *mOut) ;
MATRIX  *MatrixScalarAdd( const MATRIX *mIn, const float val, MATRIX *mOut) ;
MATRIX  *VectorZeroMean(const MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixClear(MATRIX *mat) ;
MATRIX  *MatrixSquareElts(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixSignedSquareElts(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixSqrtElts(MATRIX *mIn, MATRIX *mOut) ;
MATRIX  *MatrixDiag(MATRIX *mDiag, MATRIX *mOut) ;
MATRIX  *MatrixMakeDiagonal(MATRIX *mSrc, MATRIX *mDst) ;
MATRIX  *MatrixCopyRegion( const MATRIX *mSrc, MATRIX *mDst,
			   const int start_row,
			   const int startol,
			   const int rows,
			   const int cols,
			   const int dest_row,
			   const int dest_col );
MATRIX  *MatrixCopyRealRegion( const MATRIX *mSrc, MATRIX *mDst,
			       const int start_row,
			       const int start_col,
			       const int rows,
			       const int cols,
			       const int dest_row,
			       const int dest_col  );
MATRIX  *MatrixCopyImagRegion( const MATRIX *mSrc, MATRIX *mDst,
			       const int start_row,
			       const int start_col,
			       const int rows,
			       const int cols,
			       const int dest_row,
			       const int dest_col );
MATRIX  *MatrixSetRegion(MATRIX *mSrc, MATRIX *mDst, int start_row,
                         int start_col, int rows, int cols, float val);
MATRIX *MatrixRealToComplex(MATRIX *mReal, MATRIX *mImag, MATRIX *mOut);
MATRIX *MatrixRegularize(MATRIX *mIn, MATRIX *mOut) ;
int    MatrixSingular(MATRIX *m) ;
MATRIX *MatrixToeplitz(VECTOR *v, MATRIX *T, int Type);

/* determinants and eigenvectors */
float  MatrixDeterminant(MATRIX *m) ;
MATRIX *MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_dst) ;

MATRIX *MatrixSVD(MATRIX *mA, VECTOR *v_z, MATRIX *mV) ;
MATRIX *MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse) ;
float  MatrixNSConditionNumber(MATRIX *m);
float  MatrixConditionNumber(MATRIX *m) ;
float  MatrixSVDEigenValues(MATRIX *m, float *evalues) ;
MATRIX *MatrixFactorSqrSVD(MATRIX *M, int Invert, MATRIX *D);

/* statistical stuff */
MATRIX *MatrixCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans) ;
MATRIX *MatrixUpdateCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans);
MATRIX *MatrixUpdateMeans(MATRIX *mInputs, VECTOR *mMeans, VECTOR *mNobs) ;
MATRIX *MatrixFinalMeans(VECTOR *mMeans, VECTOR *mNobs) ;
MATRIX *MatrixFinalCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mNobs);

/* misc. I/O functions */
int    MatrixAsciiWriteInto(FILE *fp, MATRIX *m) ;
MATRIX *MatrixAsciiReadFrom(FILE *fp, MATRIX *m) ;
int    MatrixAsciiWrite(const char *fname, MATRIX *m) ;
MATRIX *MatrixAsciiRead(const char *fname, MATRIX *m) ;
MATRIX *MatrixAsciiReadRaw(const char *fname, MATRIX *m) ;

int    MatrixWriteInto(FILE *fp, MATRIX *m) ;
MATRIX *MatrixReadFrom(FILE *fp, MATRIX *m) ;

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
#define V3_LOAD       VECTOR3_LOAD

#include <math.h>
#include "macros.h"

#define XYZ_LOAD(v,x,y,z)             do { XYZ* xyz = &v; xyz.x=x, xyz.y=y, xyz.z=z; } while 0
void XYZ_NORMALIZED_LOAD(XYZ* xyz, float* xyz_length, float x, float y, float z);

float XYZApproxAngle(
    XYZ const * normalizedXYZ, float x2, float y2, float z2);
    
float XYZApproxAngle_knownLength(
    XYZ const * normalizedXYZ, 
    float x2, float y2, float z2, float length2);

double Vector3Angle(VECTOR *v1, VECTOR *v2) ;
float  VectorLen( const VECTOR *v ) ;
float  VectorAngle( const VECTOR *v1, const VECTOR *v2) ;
float  VectorDot( const VECTOR *v1, const VECTOR *v2) ;
float  VectorNormalizedDot(VECTOR *v1, VECTOR *v2) ;
float  VectorDistance(VECTOR *v1, VECTOR *v2) ;
double MatrixMahalanobisDistance(VECTOR *v_mean, MATRIX *m_cov, VECTOR *v);
double MatrixTransformDistance(MATRIX *m1, MATRIX *m2, double radius);
VECTOR *MatrixColumn(MATRIX *m, VECTOR *v, int col) ;
MATRIX *VectorOuterProduct(VECTOR *v1, VECTOR *v2, MATRIX *m) ;
VECTOR *VectorCrossProduct( const VECTOR *v1, const VECTOR *v2, VECTOR *vdst) ;
VECTOR *VectorCrossProductD(const VECTOR *v1, const VECTOR *v2, VECTOR *vdst);
float  VectorTripleProduct( const VECTOR *v1,
                            const VECTOR *v2,
                            const VECTOR *v3) ;
VECTOR *VectorNormalize( const VECTOR *vin, VECTOR *vout) ;
MATRIX *MatrixNormalizeCol(MATRIX *m, MATRIX *mcnorm, MATRIX *scale);
MATRIX *MatrixNormalizeColScale(MATRIX *m, MATRIX *scale);

/* these are macro versions that work on 3-d vectors */
#define RV3_X(v)     (RVECTOR_ELT(v,1))
#define RV3_Y(v)     (RVECTOR_ELT(v,2))
#define RV3_Z(v)     (RVECTOR_ELT(v,3))
#define V3_X(v)      (VECTOR_ELT(v,1))
#define V3_Y(v)      (VECTOR_ELT(v,2))
#define V3_Z(v)      (VECTOR_ELT(v,3))
#define V3_LEN(v)    (sqrt(V3_X(v)*V3_X(v)+V3_Y(v)*V3_Y(v)+V3_Z(v)*V3_Z(v)))
#define V3_LEN_IS_ZERO(v)    (DZERO(V3_LEN_SQ(v)))
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
#define V3_NORM(v) \
         { double len ; \
             len = V3_LEN(v) ; \
             if (FZERO(len)) len = 1; \
             V3_X(v)/=len ;\
             V3_Y(v)/=len ;\
             V3_Z(v)/=len ;}

#define V3_SUBTRACT(va,vb,vc) \
                 V3_X(vc) = V3_X(va)-V3_X(vb), \
                 V3_Y(vc) = V3_Y(va)-V3_Y(vb), \
                 V3_Z(vc) = V3_Z(va)-V3_Z(vb) ;
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
MATRIX *MatrixReallocRotation(int n, float angle, int which, MATRIX *m) ;
MATRIX *MatrixAllocTranslation(int n, double *trans) ;
#define MatrixClone(mat)   MatrixZero(0, 0, MatrixCopy(mat, NULL))
#define VectorClone        MatrixClone

float MatrixTrace(MATRIX *M);
MATRIX *MatrixVertCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat);
MATRIX *MatrixHorCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat);

MATRIX *MatrixConstVal(float val, int rows, int cols, MATRIX *X);
MATRIX *MatrixZero(int rows, int cols, MATRIX *X);
double  MatrixSumElts(MATRIX *m) ;
MATRIX *MatrixSum(MATRIX *m, int dim, MATRIX *msum);
MATRIX *MatrixSumSquare(MATRIX *m, int dim, MATRIX *msumsq);
MATRIX *MatrixDRand48(int rows, int cols, MATRIX *m);
MATRIX *MatrixDRand48ZeroMean(int rows, int cols, MATRIX *m) ;
MATRIX *MatrixSimilarityTransform(MATRIX *m_src, MATRIX *m_mul, MATRIX *m_dst);

double VectorSum(const MATRIX *v);
double VectorMean(const MATRIX *v);
double VectorVar(MATRIX *v, double *pMean);
double VectorStdDev(MATRIX *v, double *pMean);
double VectorRange(MATRIX *v, double *pVmin, double *pVmax);

MATRIX *GaussianMatrix(int len, float std, int norm, MATRIX *G);
MATRIX *GaussianVector(int len, float mean, float std, int norm, MATRIX *g);
MATRIX *MatrixReorderRows(MATRIX *X, int *NewRowOrder, MATRIX *XRO);
int MatrixRandPermRows(MATRIX *X);
int MatrixColsAreNotOrthog(MATRIX *X);
#define VectorSSE(v1, v2)  MatrixSSE(v1, v2)
#define VectorRMS(v1, v2)  MatrixRMS(v1, v2)
double MatrixSSE(MATRIX *m1, MATRIX *m2) ;
double MatrixRMS(MATRIX *m1, MATRIX *m2) ;
int MatrixOrthonormalizeTransform(MATRIX *m_L) ;
int MatrixToRigidParameters(MATRIX *m, double *pxr, double *pyr, double *pzr,
                            double *pxt, double *pyt, double *pzt);
MATRIX *MatrixFromRigidParameters(MATRIX *m, double xr, double yr, double zr,
                                  double xt, double yt, double zt);

int MatrixCheckFinite(MATRIX *m);
double MatrixRowDotProduct(MATRIX *m, int row, VECTOR *v) ;
MATRIX *MatrixKron(MATRIX *m1, MATRIX *m2, MATRIX *k);
MATRIX *MatrixDemean(MATRIX *M, MATRIX *Mdm);
MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude);
MATRIX *MatrixCumTrapZ(MATRIX *y, MATRIX *t, MATRIX *yz);

MATRIX *ANOVAOmnibus(int nLevels);
MATRIX *ANOVASelectionVector(int nLevels, int Level);
MATRIX *ANOVASummingVector(int nLevels);
MATRIX *ANOVAContrast(int *FLevels, int nFactors, int *FactorList, int nFactorList);
MATRIX *MatrixMtM(MATRIX *m, MATRIX *mout);
MATRIX *MatrixAtB(MATRIX *A, MATRIX *B, MATRIX *mout);
MATRIX *MatrixSkew(MATRIX *y, MATRIX *s);
MATRIX *MatrixKurtosis(MATRIX *y, MATRIX *k);
double MatrixMaxAbsDiff(MATRIX *m1, MATRIX *m2, double dthresh);
MATRIX *MatrixColNullSpace(MATRIX *M, int *err);
MATRIX *MatrixResidualForming(MATRIX *X, MATRIX *R);
MATRIX *MatrixGlmFit(MATRIX *y, MATRIX *X, double *pRVar, MATRIX *beta);
MATRIX *MatrixACF2Kernel(MATRIX *acf, MATRIX *k, MATRIX *F);

#endif
