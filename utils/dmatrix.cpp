/**
 * @brief Matrix utilities using doubles instead of floats
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

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "diag.h"
#include "error.h"
#include "matrix.h"
#include "dmatrix.h"

/*!
  \fn int DMatrixTest(void)
  \brief Checks DMATRIX operations against their float Matrix counterparts. These
  are only negative tests (ie, expect it not to fail).
*/
int DMatrixTest(void)
{
  DMATRIX *m1d,*m2d,*md;
  MATRIX *m1f,*m2f,*mf;
  int err=0;
  double v,vf,vd;

  printf("DMatrix: testing constant value ...");
  mf = MatrixConstVal(14,7,2,NULL);
  md = DMatrixConstVal(14,7,2,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix add %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  m1f = MatrixDRand48(14, 7, NULL);
  m2f = MatrixDRand48(14, 7, NULL);
  m1d = DMatrixCopyFMatrix(m1f,NULL);
  m2d = DMatrixCopyFMatrix(m2f,NULL);

  printf("DMatrix: testing matrix add ...");
  mf = MatrixAdd(m1f,m2f,NULL);
  md = DMatrixAdd(m1d,m2d,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix add %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  printf("DMatrix: testing matrix subtract ...");
  mf = MatrixSubtract(m1f,m2f,NULL);
  md = DMatrixSubtract(m1d,m2d,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix subtract %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  printf("DMatrix: testing matrix scalar multiplication by 2.3 ...");
  mf = MatrixScalarMul(m1f,2.3,NULL);
  md = DMatrixScalarMul(m1d,2.3,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix scalar mul %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  // Need new matrices for next tests
  MatrixFree(&m1f);
  MatrixFree(&m2f);
  DMatrixFree(&m1d);
  DMatrixFree(&m2d);
  m1f = MatrixDRand48(7, 14, NULL);
  m2f = MatrixDRand48(14, 7, NULL);
  m1d = DMatrixCopyFMatrix(m1f,NULL);
  m2d = DMatrixCopyFMatrix(m2f,NULL);

  printf("DMatrix: testing matrix multiply ...");
  mf = MatrixMultiplyD(m1f,m2f,NULL);
  md = DMatrixMultiply(m1d,m2d,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix multiply %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  DMatrixFree(&md);
  MatrixFree(&mf);

  printf("DMatrix: testing matrix transpose ...");
  mf = MatrixTranspose(m1f,mf);
  md = DMatrixTranspose(m1d,md);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix transpose %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  DMatrixFree(&md);
  MatrixFree(&mf);

  // Need new matrices for next tests
  MatrixFree(&m1f);
  MatrixFree(&m2f);
  DMatrixFree(&m1d);
  DMatrixFree(&m2d);
  m1f = MatrixDRand48(3, 1, NULL);
  m2f = MatrixDRand48(3, 1, NULL);
  m1d = DMatrixCopyFMatrix(m1f,NULL);
  m2d = DMatrixCopyFMatrix(m2f,NULL);

  printf("DMatrix: testing vector cross ...");
  mf = VectorCrossProduct(m1f,m2f,NULL);
  md = DVectorCrossProduct(m1d,m2d,NULL);
  if(md==NULL) return(1);
  v = DMatrixCompareFMatrix(mf,md);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed vector cross %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = [");
    MatrixPrint(stdout,mf);
    printf("];\n");
    printf("resd = [");
    DMatrixPrintFmt(stdout,"%8.5lf",md);
    printf("];\n");
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  printf("DMatrix: testing matrix dot ...");
  vf = VectorDot(m1f,m2f);
  vd = DVectorDot(m1d,m2d);
  if(vd<0) return(1);
  v = fabs(vf-vd);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed matrix transpose %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = %lf\n",vf);
    printf("resd = %lf\n",vd);
  }
  else printf("passed\n");
  MatrixFree(&mf);
  DMatrixFree(&md);

  printf("DMatrix: testing vector length ...");
  vf = VectorLen(m1f);
  vd = DVectorLen(m1d);
  if(vd<0) return(1);
  v = fabs(vf-vd);
  if(v > .0001){
    err++;
    printf("\nDMatrix: failed vector length %g\n",v);
    printf("M1 = [");
    MatrixPrint(stdout,m1f);
    printf("];\n");
    printf("M2 = [");
    MatrixPrint(stdout,m2f);
    printf("];\n");
    printf("resf = %lf\n",vf);
    printf("resd = %lf\n",vd);
  }
  else printf("passed\n");

  // done
  MatrixFree(&m1f);
  MatrixFree(&m2f);
  DMatrixFree(&m1d);
  DMatrixFree(&m2d);

  if(err) printf("DMatrix: failed with %d errors \n",err);
  else printf("DMatrix: passed all tests\n");

  return(err);
}

/*!
  \fn int DMatrixCheckDims(const DMATRIX *m1, const DMATRIX *m2, const int checktype, FILE *fp, const char *str)
  \brief Checks the dimensions of two matrices, returns 0 if no error. Returns non-zero
    if either matrix is NULL.
  \param m1 first matrix
  \param m2 second matrix
  \param checktype 
     1 only check rows
     2 only check columns
     3 check both rows and colums
     4 check cols of m1 against rows of m2
  \param fp  Print to the given FILE stream
  \param str String to print as part of the output
*/
int DMatrixCheckDims(const DMATRIX *m1, const DMATRIX *m2, const int checktype, FILE *fp, const char *str)
{
  if(m1 == NULL){
    fprintf(fp,"ERROR: %s matrix 1 is NULL\n",str);
    fprintf(fp,"break %s:%d\n", __FILE__, __LINE__);
    return(1);
  }
  if(m2 == NULL){
    fprintf(fp,"ERROR: %s matrix 2 is NULL\n",str);
    fprintf(fp,"break %s:%d\n", __FILE__, __LINE__);
    return(1);
  }
  if((checktype == 1 || checktype == 3) && m1->rows != m2->rows){
    fprintf(fp,"ERROR: %s rows mismatch %d, %d\n",str,m1->rows,m2->rows);
    fprintf(fp,"break %s:%d\n", __FILE__, __LINE__);
    return(1);
  }
  if((checktype == 2 || checktype == 3) && m1->cols != m2->cols){
    fprintf(fp,"ERROR: %s cols mismatch %d, %d\n",str,m1->cols,m2->cols);
    fprintf(fp,"break %s:%d\n", __FILE__, __LINE__);
    return(1);
  }
  if(checktype == 4 && m1->cols != m2->rows){
    fprintf(fp,"ERROR: %s cols of m1 mismatch rows of m2 %d, %d\n",str,m1->cols,m2->rows);
    fprintf(fp,"break %s:%d\n", __FILE__, __LINE__);
    return(1);
  }
  return(0);
}

/*!
  \fn DMATRIX *DMatrixAlloc(const int rows, const int cols, const int type)
  \brief Allocates the DMATRIX
*/
DMATRIX *DMatrixAlloc(const int rows, const int cols, const int type)
{
  DMATRIX *mat;
  int row, nelts;

  if(type == MATRIX_COMPLEX) {
    printf("ERROR: DMatrixAlloc(): complex types not supported\n");
    return(NULL);
  }
  if(type != MATRIX_REAL) {
    printf("ERROR: DMatrixAlloc(): types %d not recognized\n",MATRIX_REAL);
    return(NULL);
  }

  mat = (DMATRIX *)calloc(1, sizeof(DMATRIX));
  if (!mat) ErrorExit(ERROR_NO_MEMORY, "DMatrixAlloc(%d, %d, %d): could not allocate mat", rows, cols, type);

  mat->rows  = rows;
  mat->cols  = cols;
  mat->type  = type;
  nelts = rows * cols;

  /*Allocate a single array the size of the matrix, then initialize
    row pointers to the proper locations. */
  /*Because NRC is one-based, we must leave room for a few unused (but
    valid) addresses before the start of the actual data so that
    mat->rptr[0][0] is a valid address even though it wont be used.
  */
  mat->data = (DMATRIX_TYPE *)calloc(nelts + 2, sizeof(DMATRIX_TYPE));
  if(!mat->data){
    printf("ERROR: DMatrixAlloc(%d, %d): allocation failed\n", rows, cols);
    return(NULL);
  }
  mat->data += 2;
  /*DMATRIX type requires 1-based stuff. The full data array is zero
    based, point the first row to the zeroth element, and so on.*/
  mat->rptr = (DMATRIX_TYPE **)calloc(rows + 1, sizeof(DMATRIX_TYPE *));
  if (!mat->rptr) {
    free(mat->data);
    free(mat);
    ErrorExit(ERROR_NO_MEMORY, "DMatrixAlloc(%d, %d): could not allocate rptr", rows, cols);
  }
  for (row = 1; row <= rows; row++) {
    switch (type) {
      case MATRIX_REAL:
        mat->rptr[row] = mat->data + (row - 1) * cols - 1;
        break;
    }
  }
  return (mat);
}
/*!
  \fn int DMatrixFree(DMATRIX **pmat)
  \brief Free the DMATRIX
*/
int DMatrixFree(DMATRIX **pmat)
{
  DMATRIX *mat;
  if (!pmat) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "DMatrixFree: NULL pmat POINTER!\n"));
  mat = *pmat;
  *pmat = NULL;
  if(!mat) return (0);
  mat->data -= 2; /* silly numerical recipes in C requires 1-based stuff */
  free(mat->data);
  free(mat->rptr);
  free(mat);
  return (0);
}

/*!
  \fn DMATRIX *DMatrixMultiply(const DMATRIX *m1, const DMATRIX *m2, DMATRIX *m3)
  \brief Multiplies two matrices. The accumulation is done with long double.
*/
DMATRIX *DMatrixMultiply(const DMATRIX *m1, const DMATRIX *m2, DMATRIX *m3)
{
  int col, row, i, rows, cols, m1_cols, err;
  DMATRIX_TYPE *r3;
  DMATRIX_TYPE *r1, *r2;
  long double val;

  err = DMatrixCheckDims(m1, m2, 4, stdout, "DMatrixMultiply(): ");
  if(err) return(NULL);

  if(!m3) {
    m3 = DMatrixAlloc(m1->rows, m2->cols, m1->type);
    if(!m3) return (NULL);
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols)) {
    printf("DMatrixMultiply(): m1/m2 dim mismatch\n break %s:%d\n", __FILE__, __LINE__);
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "DMatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows,m1->cols,m2->rows, m2->cols, m3->rows, m3->cols));
  }
  if((m1->type != MATRIX_REAL) || (m2->type != MATRIX_REAL) || (m3->type != MATRIX_REAL)) {
    printf("DMatrixMultiply(): all matrices must be REAL \n break %s:%d\n", __FILE__, __LINE__);
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "DMatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows,m1->cols,m2->rows, m2->cols, m3->rows, m3->cols));
  }

  cols = m3->cols;
  rows = m3->rows;
  m1_cols = m1->cols;
  for (row = 1; row <= rows; row++) {
    r3 = &m3->rptr[row][1];
    for (col = 1; col <= cols; col++) {
      val = 0.0;
      r1 = &m1->rptr[row][1];
      r2 = &m2->rptr[1][col];
      for (i = 1; i <= m1_cols; i++, r2 += cols) val += ((DMATRIX_TYPE)(*r1++) * (*r2));
      *r3++ = val;
    }
  }
  return (m3);
}

/*!
  \fn int DMatrixPrintFmt(FILE *fp, const char *fmt, DMATRIX *mat)
  \brief Print a matrix with the given format string (eg, "8.5lf")
*/
int DMatrixPrintFmt(FILE *fp, const char *fmt, DMATRIX *mat)
{
  int row, col, rows, cols;

  if (fp == NULL) {
    fp = stdout;
    ErrorPrintf(ERROR_BADPARM, "DMatrixPrint: fp = NULL!");
  }
  if(mat == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "DMatrixPrint: mat = NULL!"));

  rows = mat->rows;
  cols = mat->cols;
  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, fmt, mat->rptr[row][col]);
          break;
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: unknown type %d\n", mat->type));
      }
      if (col < cols) fprintf(fp, "  ");
    }
    fprintf(fp, "\n");  // DO NOT PRINT the semi-colon!!!
  }
  fflush(fp);
  return (NO_ERROR);
}
/*!
  \fn int DMatrixPrint(FILE *fp, DMATRIX *mat)
  \brief Print a matrix with format string "8.5lf"
*/
int DMatrixPrint(FILE *fp, DMATRIX *mat)
{
  return(DMatrixPrintFmt(fp, "%8.5lf", mat));
}

/*!
  \fn double DMatrixCompareFMatrix(MATRIX *mf, DMATRIX *md)
  \brief Compare a float matrix to a DMATRIX, return
  the max abs diff. Returns -1 on error.
*/
double DMatrixCompareFMatrix(MATRIX *mf, DMATRIX *md)
{
  int r,c;
  double d,dmax;
  if(mf == NULL){
    printf("ERROR: DMatrixCompareFMatrix(): mf is NULL\n");
    return(-1);
  }
  if(md == NULL){
    printf("ERROR: DMatrixCompareFMatrix(): md is NULL\n");
    return(-1);
  }
  if(mf->rows != md->rows){
    printf("ERROR: DMatrixCompareFMatrix(): row mismatch %d %d\n",mf->rows,md->rows);
    return(-1);
  }
  if(mf->cols != md->cols){
    printf("ERROR: DMatrixCompareFMatrix(): col mismatch %d %d\n",mf->cols,md->cols);
    return(-1);
  }
  dmax = 0;
  for(r=1; r <= mf->rows; r++){
    for(c=1; c <= mf->cols; c++){
      d = fabs(mf->rptr[r][c] - md->rptr[r][c]);
      if(dmax < d) dmax = d;
    }
  }
  return(dmax);
}

/*!
  \fn DMATRIX *DMatrixCopyFMatrix(MATRIX *mf, DMATRIX *md)
  \brief Copy a float matrix into a DMATRIX.
*/
DMATRIX *DMatrixCopyFMatrix(MATRIX *mf, DMATRIX *md)
{
  int r,c;
  if(mf == NULL){
    printf("ERROR: DMatrixCopyFMatrix(): mf is NULL\n");
    return(NULL);
  }
  if(md==NULL)
    md = DMatrixAlloc(mf->rows,mf->cols,MATRIX_REAL);
  if(mf->rows != md->rows){
    printf("ERROR: DMatrixCopyFMatrix(): row mismatch %d %d\n",mf->rows,md->rows);
    return(NULL);
  }
  if(mf->cols != md->cols){
    printf("ERROR: DMatrixCopyFMatrix(): col mismatch %d %d\n",mf->cols,md->cols);
    return(NULL);
  }
  for(r=1; r <= mf->rows; r++){
    for(c=1; c <= mf->cols; c++){
      md->rptr[r][c] = mf->rptr[r][c];
    }
  }
  return(md);
}

/*!
  \fn DMATRIX *DMatrixCopy(DMATRIX *msrc, DMATRIX *mcopy)
  \brief Copy one matrix into another
*/
DMATRIX *DMatrixCopy(DMATRIX *msrc, DMATRIX *mcopy)
{
  int r,c,err;
  if(msrc == NULL){
    printf("ERROR: DMatrixCopy(): msrc is NULL\n");
    return(NULL);
  }
  if(mcopy==NULL)
    mcopy = DMatrixAlloc(msrc->rows,msrc->cols,MATRIX_REAL);
  err = DMatrixCheckDims(msrc, mcopy, 3, stdout, "DMatrixCopy(): ");
  if(err) return(NULL);
  for(r=1; r <= msrc->rows; r++){
    for(c=1; c <= msrc->cols; c++){
      mcopy->rptr[r][c] = msrc->rptr[r][c];
    }
  }
  return(mcopy);
}

/*!
  \fn DMATRIX *DMatrixTranspose(DMATRIX *mIn, DMATRIX *mOut)
  \brief Compute transpose. Cannot be done in-place
*/
DMATRIX *DMatrixTranspose(DMATRIX *mIn, DMATRIX *mOut)
{
  int row, col;

  if(mIn == mOut){
    printf("DMatrixTranspose(): ERROR: cannot be done in-place\n");
    return(NULL);
  }
  if (!mOut) {
    mOut = DMatrixAlloc(mIn->cols, mIn->rows, mIn->type);
    if(!mOut) return (NULL);
  }
  if(mIn->rows != mOut->cols){
    printf("DMatrixTranspose(): ERROR: row/col mismatch %d %d\n",mIn->rows,mOut->cols);
    return(NULL);
  }
  if(mIn->cols != mOut->rows){
    printf("DMatrixTranspose(): ERROR: col/row mismatch %d %d\n",mIn->cols,mOut->rows);
    return(NULL);
  }
  for(row = 1; row <= mIn->rows; row++) {
    for(col = 1; col <= mIn->cols; col++) {
      mOut->rptr[col][row] = mIn->rptr[row][col];
    }
  }
  return(mOut);
}

/*!
  \fn DMATRIX *DMatrixAddMul(DMATRIX *m1, DMATRIX *m2, double v1, double v2, DMATRIX *mOut)
  \brief Computes v1*m1+v2*m2 so it can be used for adding, subtracting, and scaling.
  If a matrix is NULL or its value is 0, it is ignored.
*/
DMATRIX *DMatrixAddMul(DMATRIX *m1, DMATRIX *m2, double v1, double v2, DMATRIX *mOut)
{
  int row, err;

  if(m1 && m2){
    err = DMatrixCheckDims(m1, m2, 3, stdout, "DMatrixAddMul(): ");
    if(err) return(NULL);
  }
  if (!mOut) {
    mOut = DMatrixAlloc(m1->rows, m1->cols, m1->type);
    if(!mOut) return (NULL);
  }
  if(m1){
    err = DMatrixCheckDims(m1, mOut, 3, stdout, "DMatrixAddMul(): ");
    if(err) return(NULL);
  }
  else if(m2){
    err = DMatrixCheckDims(m2, mOut, 3, stdout, "DMatrixAddMul(): ");
    if(err) return(NULL);
  }

  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for(row = 1; row <= m1->rows; row++) {
    int col;
    double val;
    for(col = 1; col <= m1->cols; col++) {
      val = 0;
      if(m1){
	if(v1 != 0) {
	  if(v1 != 1) val += (v1*m1->rptr[row][col]);
	  else        val += (   m1->rptr[row][col]);
	}
      }
      if(m2){
	if(v2 != 0) {
	  if(v2 != 1) val += (v2*m2->rptr[row][col]);
	  else        val += (   m2->rptr[row][col]);
	}
      }
      mOut->rptr[row][col] = val;
    }
  }
  return(mOut);
}

/*! 
  \fn DMATRIX *DMatrixScalarMul(DMATRIX *m, double v, DMATRIX *mout)
  \brief Scale a matrix by the given value
*/
DMATRIX *DMatrixScalarMul(DMATRIX *m, double v, DMATRIX *mout)
{
  // check dim not needed, done in DMatrixAddMul
  return(DMatrixAddMul(m, NULL, v, 0, mout));
}

/*! 
  \fn DMATRIX *DMatrixSubtract(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout)
  \brief Subtract two matrices element-by-element
*/
DMATRIX *DMatrixSubtract(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout)
{
  // check dim not needed, done in DMatrixAddMul
  return(DMatrixAddMul(m1, m2, 1, -1, mout));
}

/*! 
  \fn DMATRIX *DMatrixAdd(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout)
  \brief Add two matrices element-by-element
*/
DMATRIX *DMatrixAdd(DMATRIX *m1, DMATRIX *m2, DMATRIX *mout)
{
  // check dim not needed, done in DMatrixAddMul
  return(DMatrixAddMul(m1, m2, 1, 1, mout));
}

/*! 
  \fn double DVectorDot(const DVECTOR *v1, const DVECTOR *v2)
  \brief Compute the dot product of two vectors
*/
double DVectorDot(const DVECTOR *v1, const DVECTOR *v2)
{
  int err, i;
  long double dot;
  err = DMatrixCheckDims(v1, v2, 3, stdout, "DDVectorDot(): ");
  if(err) return(-1);
  dot = 0.0;
  for (i = 1; i <= v1->rows; i++) dot += ((long double)v1->rptr[i][1] * v2->rptr[i][1]);
  return (dot);
}

/*! 
  \fn MATRIX *DMatrixConstVal(double val, int rows, int cols, MATRIX *X)
  \brief sets all the elements to the given value.  If X is NULL, then
  a matrix rows-by-cols is alloced. If X is non-NULL, then rows and
  cols are ignored.
*/
DMATRIX *DMatrixConstVal(const double val, const int rows, const int cols, DMATRIX *X)
{
  int r;
  if(val ==0) {
    X = DMatrixZero(rows, cols, X);
    return(X);
  }
  if(X == NULL) X = DMatrixAlloc(rows, cols, MATRIX_REAL);
  #ifdef HAVE_OPENMP
  #pragma omp parallel for 
  #endif
  for (r = 1; r <= X->rows; r++) {
    int c;
    for (c = 1; c <= X->cols; c++) {
      X->rptr[r][c] = val;
    }
  }
  return (X);
}

/*! 
  \fn MATRIX *DMatrixZero(double val, int rows, int cols, MATRIX *X)
  \brief sets all the elements to 0.  If X is NULL, then a matrix
  rows-by-cols is alloced. If X is non-NULL, then rows and cols are
  ignored. This is different than DMatrixConstVal(0) in that it
  it uses memset(), which should be much faster. Now DMatrixConstVal(0)
  will actually call DMatrixZero().
*/
DMATRIX *DMatrixZero(const int rows, const int cols, DMATRIX *X)
{
  if(X == NULL) {
    X = DMatrixAlloc(rows, cols, MATRIX_REAL);
    return(X);
  }
  memset(X->data,0,X->rows*X->cols*sizeof(DMATRIX_TYPE));
  return(X);
}
/*! 
  \fn double DVectorLen(const DVECTOR *v)
  \brief Computes the lenght of the vector;
*/
double DVectorLen(const DVECTOR *v)
{
  int i;
  long double len, vi;
  if(v==NULL){
    printf("ERROR: DVectorLen() vector is NULL\n");
    return(-1);
  }
  len = 0.0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction (+:len)
  #endif
  for (i = 1; i <= v->rows; i++) {
    vi = v->rptr[i][1];
    len += vi * vi;
  }
  len = sqrt(len);
  return (len);
}

/*!
  \fn DVECTOR *DVectorCrossProduct(const DVECTOR *v1, const DVECTOR *v2, DVECTOR *vdst)
  Calculates the cross product of two vectors and returns the result
  in vdst, allocating it if necessary. Same as VectorCrossProduct() but uses
  long doubles internally instead of floats.
*/
DVECTOR *DVectorCrossProduct(const DVECTOR *v1, const DVECTOR *v2, DVECTOR *vdst)
{
  long double x1, x2, y1, y2, z1, z2;

  if(v1 == NULL){
    printf("ERROR: DVectorCrossProduct(): v1 is NULL\n");
    return(NULL);
  }
  if(v2 == NULL){
    printf("ERROR: DVectorCrossProduct(): v2 is NULL\n");
    return(NULL);
  }
  if(v1->rows != 3 || v1->cols != 1){
    printf("ERROR: DVectorCrossProduct(): v1 is not 3x1 %d %d\n",v1->rows,v1->cols);
    return(NULL);
  }
  if(v2->rows != 3 || v2->cols != 1){
    printf("ERROR: DVectorCrossProduct(): v2 is not 3x1 %d %d\n",v2->rows,v2->cols);
    return(NULL);
  }
  if(!vdst) vdst = DMatrixAlloc(3,1,MATRIX_REAL);
  if(vdst->rows != 3 || vdst->cols != 1){
    printf("ERROR: DVectorCrossProduct(): vdst is not 3x1 %d %d\n",vdst->rows,vdst->cols);
    return(NULL);
  }

  x1 = v1->rptr[1][1];
  y1 = v1->rptr[2][1];
  z1 = v1->rptr[3][1];
  x2 = v2->rptr[1][1];
  y2 = v2->rptr[2][1];
  z2 = v2->rptr[3][1];
  vdst->rptr[1][1] = y1 * z2 - z1 * y2;
  vdst->rptr[2][1] = z1 * x2 - x1 * z2;
  vdst->rptr[3][1] = x1 * y2 - y1 * x2;

  return (vdst);
}

/*!
  \fn DMATRIX_TYPE DMatrixMaxAbs(DMATRIX *M)
  Computes the maximum of the absolute value of the given matrix
*/
DMATRIX_TYPE DMatrixMaxAbs(DMATRIX *M)
{
  DMATRIX_TYPE vmax=0;
  int r, c;
  if(M == NULL){
    printf("ERROR: DMatrixMaxAbs(): matrix is NULL\n");
    return(-1);
  }
  for (r = 1; r <= M->rows; r++) {
    for (c = 1; c <= M->cols; c++) {
      if(vmax < fabs(M->rptr[r][c])) {
	vmax = fabs(M->rptr[r][c]);
      }
    }
  }
  return(vmax);
}
/*----------------------------------------------------------------*/
DMATRIX *DMatrixDRand48(int rows, int cols, DMATRIX *m)
{
  int r, c;
  if (m == NULL) m = DMatrixAlloc(rows, cols, MATRIX_REAL);
  /* if m != NULL rows and cols are ignored */
  for (r = 1; r <= m->rows; r++) {
    for (c = 1; c <= m->cols; c++) {
      m->rptr[r][c] = drand48();
    }
  }
  return (m);
}
