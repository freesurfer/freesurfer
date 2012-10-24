/**
 * @file  matrix.c
 * @brief Matrix utilities
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/10/24 13:46:17 $
 *    $Revision: 1.129 $
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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>
#ifdef _POSIX_MAPPED_FILES
#include <sys/mman.h>
#endif
#ifdef Darwin
#include <float.h>  // defines FLT_MIN
#else
#ifdef SunOS
#include <limits.h> // defines FLT_MIN
#else
#ifdef Windows_NT
#include <float.h> // defines FLT_MIN
#else
#include <values.h> // defines FLT_MIN
#endif
#endif
#endif

#include "matrix.h"
#include "proto.h"
#include "error.h"
#include "utils.h"
#include "diag.h"
#include "macros.h"
#include "numerics.h"
#include "evschutils.h"

// private functions
MATRIX *MatrixCalculateEigenSystemHelper( MATRIX *m,
    float *evalues,
    MATRIX *m_evectors,
    int isSymmetric );

/**
 * Returns true if the matrix is symmetric (should be square too).
 */
int MatrixIsSymmetric( MATRIX *matrix );

MATRIX *
MatrixCopy( const MATRIX *mIn, MATRIX *mOut )
{
  int row, rows, cols, col ;

  if (mIn == NULL)
  {
    if (mOut)
    {
      MatrixFree(&mOut);
      mOut = NULL;
    }
    return(NULL);
  }
  rows = mIn->rows ;
  cols = mIn->cols ;

  if (!mOut)
    mOut = MatrixAlloc(rows, cols, mIn->type) ;

  if (!mOut)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixCopy: couldn't allocate mOut") ;

#if 1
  for (col = 1 ; col <= cols ; col++)
    for (row = 1 ; row <= rows ; row++)
      *MATRIX_RELT(mOut, row, col) = *MATRIX_RELT(mIn, row, col) ;
#else
  for (row = 1 ; row <= rows ; row++)
    memmove((char *)(mOut->rptr[row]), (char *)mIn->rptr[row],
           (cols+1)*sizeof(float)) ;
#endif

  return(mOut) ;
}


MATRIX *
MatrixInverse( const MATRIX *mIn, MATRIX *mOut)
{
  float  **a, **y;
  int    isError, i, j, rows, cols, alloced = 0 ;
  MATRIX *mTmp ;

  if (!mIn)
    ErrorExit(ERROR_BADPARM,
                 "MatrixInverse: NULL input matrix!\n") ;

  if (mIn->rows != mIn->cols)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixInverse: matrix (%d x %d) is not square\n",
                 mIn->rows, mIn->cols)) ;

  rows = mIn->rows ;
  cols = mIn->cols ;

  if (!mOut)
  {
    alloced = 1 ;
    mOut = MatrixAlloc(rows, cols, mIn->type) ;
  }

  /* allocate temp matrix so as not to destory contents of mIn */
  if (mIn->type == MATRIX_COMPLEX)
  {
    MATRIX *mQuad, *mInv, *mReal, *mImag ;

    mTmp = MatrixAlloc(2*rows, 2*cols, MATRIX_REAL) ;

    /*
      form square matrix of the form

      A  -C
      C  A

      where A and C are the real and imaginary components of the input
      matrix respectively.
    */
    MatrixCopyRealRegion(mIn, mTmp, 1, 1, rows, cols, 1, 1) ;
    MatrixCopyRealRegion(mIn, mTmp, 1, 1, rows, cols, rows+1, cols+1) ;
    mQuad = MatrixAlloc(rows,cols,MATRIX_REAL) ;
    MatrixCopyImagRegion(mIn, mQuad, 1, 1, rows, cols, 1, 1) ;
    MatrixScalarMul(mQuad, -1.0f, mQuad) ;
    MatrixCopyRegion(mQuad, mTmp, 1, 1, rows, cols, 1, cols+1) ;
    MatrixCopyImagRegion(mIn, mTmp, 1, 1, rows, cols, cols+1, 1) ;

#if 0
    DiagPrintf(DIAG_MATRIX, "\nextraction\n") ;
    MatrixPrint(stdout, mTmp) ;
#endif

    mReal = MatrixAlloc(rows, cols, MATRIX_REAL) ;
    mImag = MatrixAlloc(rows, cols, MATRIX_REAL) ;
    mInv = MatrixInverse(mTmp, NULL) ;

    MatrixCopyRegion(mInv, mReal, 1, 1, rows, cols, 1, 1) ;
    MatrixCopyRegion(mInv, mImag, rows+1, 1, rows, cols, 1, 1) ;
    MatrixRealToComplex(mReal, mImag, mOut) ;

#if 0
    DiagPrintf(DIAG_MATRIX, "\ninverse of extraction\n") ;
    MatrixPrint(stderr, mInv) ;
    DiagPrintf(DIAG_MATRIX, "\ninverse:\n") ;
    MatrixPrint(stderr, mOut) ;
    DiagPrintf(DIAG_MATRIX, "real + imag = \n") ;
    MatrixPrint(stderr, mReal) ;
    MatrixPrint(stderr, mImag) ;
#endif
    MatrixFree(&mQuad) ;
    MatrixFree(&mReal) ;
    MatrixFree(&mImag) ;
  }
  else
  {
    mTmp = MatrixCopy(mIn, NULL) ;

    a = mTmp->rptr ;
    y = mOut->rptr ;

    isError = OpenLUMatrixInverse(mTmp, mOut);

    if ( isError < 0 )
    {
      MatrixFree(&mTmp);
      if (alloced)
      {
        MatrixFree(&mOut);
      }
      return(NULL) ;
    }
  }

  MatrixFree(&mTmp) ;

  for (j = 1 ; j <= rows ; j++)
  {
    for (i = 1 ; i <= rows ; i++)
      switch (mOut->type)
      {
      case MATRIX_REAL:
        if (!finite(*MATRIX_RELT(mOut, i, j)))
        {
          if (alloced)
            MatrixFree(&mOut) ;
          return(NULL) ;   /* was singular */
        }
        break ;
      case MATRIX_COMPLEX:
        if (!finite(MATRIX_CELT_REAL(mOut, i, j)) ||
            !finite(MATRIX_CELT_IMAG(mOut, i, j)))
        {
          if (alloced)
            MatrixFree(&mOut) ;
          return(NULL) ;   /* was singular */
        }
        break ;
      }
  }

  return(mOut) ;
}


MATRIX *
MatrixAlloc( const int rows, const int cols, const int type)
{
  MATRIX *mat ;
  int    row, nelts;
#ifdef _POSIX_MAPPED_FILES
  int    i;
  float  f;
#endif

  mat = (MATRIX *)calloc(1, sizeof(MATRIX)) ;
  if (!mat)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixAlloc(%d, %d, %d): could not allocate mat",
              rows, cols, type) ;

  mat->rows = rows ;
  mat->cols = cols ;
  mat->type = type ;

  /*
    allocate a single array the size of the matrix, then initialize
    row pointers to the proper locations.
  */

  nelts = rows*cols ;
  if (type == MATRIX_COMPLEX)
    nelts *= 2 ;

  /*
    because NRC is one-based, we must leave room for a few unused
    (but valid) addresses before the start of the actual data so
    that mat->rptr[0][0] is a valid address even though it wont
    be used.
  */
  mat->data = (float *)calloc(nelts+2, sizeof(float)) ;

  mat->mmapfile = NULL;

#ifdef _POSIX_MAPPED_FILES
  if (!mat->data) /* First try to allocate a mmap'd tmpfile */
  {
    printf("MatrixAlloc(%d, %d): Using mmap'd tmpfile\n",
           rows, cols) ;
    if ((mat->mmapfile = tmpfile()))
    {
      /* This maintains identical behavior with calloc */
      f = 0;
      for (i = 0; i < nelts+2; ++i)
      {
        if (!(fwrite(&f, sizeof(float), 1, mat->mmapfile)))
        {
          printf("MatrixAlloc(%d, %d): fwrite failed", rows, cols) ;
          exit(1) ;
        }
      }
      /* This seems to matter with some implementations of mmap */
      fseek(mat->mmapfile, 0, 0) ;
      /* lseek(fileno(mat->mapfile), (nelts+2) * sizeof(float), 0) ;*/
      fflush(mat->mmapfile) ;

      mat->data = (float *)mmap(0, (nelts+2) * sizeof(float),
                                PROT_READ | PROT_WRITE, MAP_SHARED,
                                fileno(mat->mmapfile), 0) ;

      if (mat->data == MAP_FAILED)
      {
        mat->data = 0 ;
      }
    }
  }
#endif

  if (!mat->data) /* we _still_ couldn't get it! */
  {
    fprintf(stderr, "MatrixAlloc(%d, %d): allocation failed\n",
            rows, cols) ;
    exit(1) ;
  }
  mat->data += 2 ;

  /*
     silly numerical recipes in C requires 1-based stuff. The full
     data array is zero based, point the first row to the zeroth
     element, and so on.
  */
  mat->rptr = (float **)calloc(rows+1, sizeof(float *)) ;
  if (!mat->rptr)
  {
    free(mat->data) ;
    free(mat) ;
    ErrorExit
    (ERROR_NO_MEMORY,
     "MatrixAlloc(%d, %d): could not allocate rptr",
     rows, cols) ;
  }
  for (row = 1 ; row <= rows ; row++)
  {
    switch (type)
    {
    case MATRIX_REAL:
      mat->rptr[row] = mat->data + (row-1)*cols - 1 ;
      break ;
    case MATRIX_COMPLEX:
      mat->rptr[row] = (float *)(((CPTR)mat->data) +
                                 (row-1)*cols - 1) ;
      break ;
    default:
      ErrorReturn(NULL,
                  (ERROR_BADPARM, "MatrixAlloc: unknown type %d\n",type)) ;
    }
  }
  return(mat) ;
}


int
MatrixFree(MATRIX **pmat)
{
  MATRIX *mat ;

  if (!pmat)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixFree: NULL pmat POINTER!\n"));

  mat = *pmat ;
  *pmat = NULL;

  if (!mat)
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixFree: NULL mat POINTER!\n"));

  /* silly numerical recipes in C requires 1-based stuff */
  mat->data -= 2 ;
  if (mat->mmapfile)
  {
    int nelts ;

    nelts = mat->rows*mat->cols ;
    if (mat->type == MATRIX_COMPLEX)
      nelts *= 2 ;

#ifdef _POSIX_MAPPED_FILES
    munmap((void *) mat->data, (nelts+2) * sizeof(float)) ;
#endif
    fclose(mat->mmapfile) ;
  }
  else
  {
    free(mat->data) ;
  }

  free(mat->rptr) ;
  free(mat) ;

  return(0) ;
}


/*!
  \fn MATRIX *MatrixMultiplyD( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
  \brief Multiplies two matrices. The accumulation is done with double,
   which is more accurate than MatrixMultiplyD() which uses float.
*/
MATRIX *MatrixMultiplyD( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
{
  int   col, row, i, rows, cols, m1_cols ;
  float *r3 ;
  register float *r1, *r2 ;
  register double val;
  MATRIX   *m_tmp1 = NULL, *m_tmp2 = NULL ;

  if (!m1)
    ErrorExit(ERROR_BADPARM,
              "MatrixMultiply: m1 is null!\n") ;
  if (!m2)
    ErrorExit(ERROR_BADPARM,
              "MatrixMultiply: m2 is null!\n") ;

  if (m1->cols != m2->rows)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiply: m1 cols %d does not match m2 rows %d\n",
                 m1->cols, m2->rows)) ;

  if (!m3)
  {
    /* twitzel also did something here */
    if ((m1->type == MATRIX_COMPLEX) || (m2->type == MATRIX_COMPLEX))
      m3 = MatrixAlloc(m1->rows, m2->cols, MATRIX_COMPLEX);
    else
      m3 = MatrixAlloc(m1->rows, m2->cols, m1->type) ;
    if (!m3)
      return(NULL) ;
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows, m1->cols, m2->rows, m2->cols, m3->rows, m3->cols)) ;

  if (m3 == m2)
  {
    m_tmp1 = MatrixCopy(m2, NULL) ;
    m2 = m_tmp1 ;
  }
  if (m3 == m1)
  {
    m_tmp2 = MatrixCopy(m1, NULL) ;
    m1 = m_tmp2 ;
  }
  /*  MatrixClear(m3) ;*/
  cols = m3->cols ;
  rows = m3->rows ;
  m1_cols = m1->cols ;

  /* twitzel modified here */
  if((m1->type == MATRIX_REAL) && (m2->type == MATRIX_REAL)) {
    for(row = 1 ; row <= rows ; row++)  {
      r3 = &m3->rptr[row][1] ;
      for (col = 1 ; col <= cols ; col++){
        val = 0.0 ;
        r1 = &m1->rptr[row][1] ;
        r2 = &m2->rptr[1][col] ;
        for (i = 1 ; i <= m1_cols ; i++, r2 += cols)
          val += *r1++ * *r2 ;
        *r3++ = val ;
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_COMPLEX))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          double a, b, c, d ;  /* a + ib and c + id */

          a = MATRIX_CELT_REAL(m1,row,i) ;
          b = MATRIX_CELT_IMAG(m1,row,i) ;
          c = MATRIX_CELT_REAL(m2,i,col) ;
          d = MATRIX_CELT_IMAG(m2,i,col) ;
          MATRIX_CELT_REAL(m3,row,col) += a*c - b*d ;
          MATRIX_CELT_IMAG(m3,row,col) += a*d + b*c ;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_COMPLEX))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          double a, c, d ;  /* a + ib and c + id and b=0 here*/

          a = *MATRIX_RELT(m1,row,i);
          c = MATRIX_CELT_REAL(m2,i,col);
          d = MATRIX_CELT_IMAG(m2,i,col);
          MATRIX_CELT_REAL(m3,row,col) += a*c;
          MATRIX_CELT_IMAG(m3,row,col) += a*d;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_REAL))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          double a, b, c ;  /* a + ib and c + id and d=0 here*/

          a = MATRIX_CELT_REAL(m1,row,i);
          b = MATRIX_CELT_IMAG(m1,row,i);
          c = *MATRIX_RELT(m2,i,col);
          MATRIX_CELT_REAL(m3,row,col) += a*c;
          MATRIX_CELT_IMAG(m3,row,col) += b*c;
        }
      }
    }
  }
  if (m_tmp1)
    MatrixFree(&m_tmp1) ;
  if (m_tmp2)
    MatrixFree(&m_tmp2) ;
  return(m3) ;
}

/*!
  \fn MATRIX *MatrixMultiply( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
  \brief Multiplies two matrices. The accumulation is done with float. 
   Consider using MatrixMultiplyD() which uses double.
*/
MATRIX *MatrixMultiply( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
{
  int   col, row, i, rows, cols, m1_cols ;
  float *r3 ;
  register float val, *r1, *r2 ;
  MATRIX   *m_tmp1 = NULL, *m_tmp2 = NULL ;

  if (!m1)
    ErrorExit(ERROR_BADPARM,
              "MatrixMultiply: m1 is null!\n") ;
  if (!m2)
    ErrorExit(ERROR_BADPARM,
              "MatrixMultiply: m2 is null!\n") ;

  if (m1->cols != m2->rows)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiply: m1 cols %d does not match m2 rows %d\n",
                 m1->cols, m2->rows)) ;

  if (!m3)
  {
    /* twitzel also did something here */
    if ((m1->type == MATRIX_COMPLEX) || (m2->type == MATRIX_COMPLEX))
      m3 = MatrixAlloc(m1->rows, m2->cols, MATRIX_COMPLEX);
    else
      m3 = MatrixAlloc(m1->rows, m2->cols, m1->type) ;
    if (!m3)
      return(NULL) ;
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows, m1->cols, m2->rows, m2->cols, m3->rows, m3->cols)) ;

  if (m3 == m2)
  {
    m_tmp1 = MatrixCopy(m2, NULL) ;
    m2 = m_tmp1 ;
  }
  if (m3 == m1)
  {
    m_tmp2 = MatrixCopy(m1, NULL) ;
    m1 = m_tmp2 ;
  }
  /*  MatrixClear(m3) ;*/
  cols = m3->cols ;
  rows = m3->rows ;
  m1_cols = m1->cols ;

  /* twitzel modified here */
  if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_REAL))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      r3 = &m3->rptr[row][1] ;
      for (col = 1 ; col <= cols ; col++)
      {
        val = 0.0 ;
        r1 = &m1->rptr[row][1] ;
        r2 = &m2->rptr[1][col] ;
        for (i = 1 ; i <= m1_cols ; i++, r2 += cols)
        {
#if 0
          m3->rptr[row][col] +=
            m1->rptr[row][i] * m2->rptr[i][col] ;
#else
          val += *r1++ * *r2 ;
#endif
        }
        *r3++ = val ;
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_COMPLEX))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          float a, b, c, d ;  /* a + ib and c + id */

          a = MATRIX_CELT_REAL(m1,row,i) ;
          b = MATRIX_CELT_IMAG(m1,row,i) ;
          c = MATRIX_CELT_REAL(m2,i,col) ;
          d = MATRIX_CELT_IMAG(m2,i,col) ;
          MATRIX_CELT_REAL(m3,row,col) += a*c - b*d ;
          MATRIX_CELT_IMAG(m3,row,col) += a*d + b*c ;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_COMPLEX))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          float a, c, d ;  /* a + ib and c + id and b=0 here*/

          a = *MATRIX_RELT(m1,row,i);
          c = MATRIX_CELT_REAL(m2,i,col);
          d = MATRIX_CELT_IMAG(m2,i,col);
          MATRIX_CELT_REAL(m3,row,col) += a*c;
          MATRIX_CELT_IMAG(m3,row,col) += a*d;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_REAL))
  {
    for (row = 1 ; row <= rows ; row++)
    {
      for (col = 1 ; col <= cols ; col++)
      {
        for (i = 1 ; i <= m1->cols ; i++)
        {
          float a, b, c ;  /* a + ib and c + id and d=0 here*/

          a = MATRIX_CELT_REAL(m1,row,i);
          b = MATRIX_CELT_IMAG(m1,row,i);
          c = *MATRIX_RELT(m2,i,col);
          MATRIX_CELT_REAL(m3,row,col) += a*c;
          MATRIX_CELT_IMAG(m3,row,col) += b*c;
        }
      }
    }
  }
  if (m_tmp1)
    MatrixFree(&m_tmp1) ;
  if (m_tmp2)
    MatrixFree(&m_tmp2) ;
  return(m3) ;
}


int MatrixPrint(FILE *fp, const MATRIX *mat)
{
  int  row, col, rows, cols ;

  if (fp == NULL)
  {
    fp = stdout ;
    ErrorPrintf(ERROR_BADPARM, "MatrixPrint: fp = NULL!") ;
  }
  if (mat == NULL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "MatrixPrint: mat = NULL!")) ;

  rows = mat->rows ;
  cols = mat->cols ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      switch (mat->type)
      {
      case MATRIX_REAL:
        fprintf(fp, "% 2.3f", mat->rptr[row][col]) ;
        break ;
      case MATRIX_COMPLEX:
        fprintf(fp, "% 2.3f + % 2.3f i",
                MATRIX_CELT_REAL(mat,row,col),
                MATRIX_CELT_IMAG(mat, row, col)) ;
        break ;
      default:
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM,
                     "MatrixPrint: unknown type %d\n",mat->type)) ;
      }
#if 0
      if (col < cols)
        fprintf(fp, " | ") ;
#else
      if (col < cols)
        fprintf(fp, "  ") ;
#endif
    }
    fprintf(fp, ";\n") ;
  }
  fflush(stdout);
  return(NO_ERROR) ;
}


int MatrixPrintFmt(FILE *fp, const char *fmt, MATRIX *mat)
{
  int  row, col, rows, cols ;

  if(fp == NULL){
    fp = stdout ;
    ErrorPrintf(ERROR_BADPARM, "MatrixPrint: fp = NULL!") ;
  }
  if(mat == NULL)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM, "MatrixPrint: mat = NULL!")) ;

  rows = mat->rows ;
  cols = mat->cols ;

  for (row = 1 ; row <= rows ; row++)  {
    for (col = 1 ; col <= cols ; col++)    {
      switch (mat->type)      {
      case MATRIX_REAL:
        fprintf(fp, fmt, mat->rptr[row][col]) ;
        break ;
      case MATRIX_COMPLEX:
        fprintf(fp, fmt, MATRIX_CELT_REAL(mat,row,col));
	fprintf(fp, " + ");
        fprintf(fp, fmt, MATRIX_CELT_IMAG(mat, row, col));
	fprintf(fp, " i");
        break ;
      default:
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM,
                     "MatrixPrint: unknown type %d\n",mat->type)) ;
      }
      if (col < cols) fprintf(fp, "  ") ;
    }
    fprintf(fp, "\n") ; // DO NOT PRINT the semi-colon!!!
  }
  fflush(fp);
  return(NO_ERROR) ;
}


int
MatrixPrintOneLine(FILE *fp, MATRIX *mat)
{
  int  row, col, rows, cols ;

  rows = mat->rows ;
  cols = mat->cols ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      switch (mat->type)
      {
      case MATRIX_REAL:
        fprintf(fp, "% 2.3f", mat->rptr[row][col]) ;
        break ;
      case MATRIX_COMPLEX:
        fprintf(fp, "% 2.3f + % 2.3f i",
                MATRIX_CELT_REAL(mat,row,col),
                MATRIX_CELT_IMAG(mat, row, col)) ;
        break ;
      default:
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM,
                     "MatrixPrint: unknown type %d\n",mat->type)) ;
      }
#if 0
      if (col < cols)
        fprintf(fp, " | ") ;
#else
      if (col < cols)
        fprintf(fp, " ") ;
#endif
    }
    fprintf(fp, ";   ") ;
  }
  return(NO_ERROR) ;
}


int
MatrixPrintTranspose(FILE *fp, MATRIX *mat)
{
  int  row, col, rows, cols ;

  rows = mat->rows ;
  cols = mat->cols ;

  for (col = 1 ; col <= cols ; col++)
  {
    for (row = 1 ; row <= rows ; row++)
    {
      switch (mat->type)
      {
      case MATRIX_REAL:
        fprintf(fp, "% 2.3f", mat->rptr[row][col]) ;
        break ;
      case MATRIX_COMPLEX:
        fprintf(fp, "% 2.3f + % 2.3f i",
                MATRIX_CELT_REAL(mat,row,col),
                MATRIX_CELT_IMAG(mat, row, col)) ;
        break ;
      default:
        ErrorReturn(ERROR_BADPARM,
                    (ERROR_BADPARM,
                     "MatrixPrint: unknown type %d\n",mat->type)) ;
      }
#if 0
      if (row < rows)
        fprintf(fp, " | ") ;
#else
      if (row < rows)
        fprintf(fp, "  ") ;
#endif
    }
    fprintf(fp, "\n") ;
  }
  fflush(stdout);
  return(NO_ERROR) ;
}


MATRIX *
MatrixReadTxt(const char *fname, MATRIX *mat)
{
  FILE   *fp ;
  int     rows, cols, row, col, nlinemax, nread;
  char    line[1000] ;
  //char    *cp;
  nlinemax = 999;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NO_FILE, "MatrixRead(%s) - file open failed\n", fname));

  // Read in the first line, including the newline
  if (!fgets(line, nlinemax, fp))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixRead: could not read 1st line from %s\n", fname));

  //for(cols = 0,cp = strtok(line, " \t,") ; cp ; cp = strtok(NULL, " \t,"))
  //  cols++ ;

  cols = ItemsInString(line);
  if (cols < 1)
  {
    // Return quitely in case try to use another format.
    return(NULL);
  }

  // Count the number of lines, start at row=1 because
  // a line was read above.
  for (rows = 1 ; fgets(line, nlinemax, fp) ; rows++)
  {}

  //printf("MatrixReadTxt: %d rows and %d cols\n",rows,cols);
  mat = MatrixAlloc(rows, cols, MATRIX_REAL) ;

  rewind(fp) ;
  for (row=1; row <= rows; row++)
  {
    for (col=1; col <= cols; col++)
    {
      nread = fscanf(fp, "%f", &mat->rptr[row][col]);
      if (nread != 1)
      {
        MatrixFree(&mat) ;
        ErrorReturn(NULL,(ERROR_BADPARM,
                          "MatrixReadTxT: could not scan value [%d][%d]\n",
                          row, col));
      }
    }
  }
  fclose(fp) ;
  return(mat) ;

#if 0
  for (row = 1 ; fgets(line, 399, fp) ; row++)
  {
    for (col = 1, cp = strtok(line, " \t,") ;
         cp ;
         cp = strtok(NULL, " \t,"),
         col++)
    {
      if (sscanf(cp, "%f", &mat->rptr[row][col]) != 1)
      {
        MatrixFree(&mat) ;
        ErrorReturn
        (NULL,
         (ERROR_BADPARM,
          "MatrixRead: could not scan value [%d][%d]\n", row, col));
      }
    }
  }

  fclose(fp) ;
  return(mat) ;
#endif
}


#include "matfile.h"
#ifdef const
#undef const
#endif
MATRIX *
MatrixRead(const char *fname)
{
  return(MatlabRead(fname)) ;
}


int
MatrixWrite(MATRIX *mat, const char *fname,const char *name)
{
  if (!name)
    return(MatlabWrite(mat,fname,fname)) ;  /* name of matrix in .mat file */
  return(MatlabWrite(mat, fname, name)) ;
}


MATRIX *
MatrixIdentity(int n, MATRIX *mat)
{
  int     i ;

  if (!mat)
    mat = MatrixAlloc(n, n, MATRIX_REAL) ;
  else
    MatrixClear(mat) ;

  for (i = 1 ; i <= n ; i++)
    mat->rptr[i][i] = 1 ;

  return(mat) ;
}


MATRIX  *
MatrixTranspose(MATRIX *mIn, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->cols, mIn->rows, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  int insitu = (mIn == mOut);
  if (insitu) mIn = MatrixCopy(mOut,NULL);

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
      mOut->rptr[col][row] = mIn->rptr[row][col] ;
  }
  if (insitu) MatrixFree(&mIn);
  return(mOut) ;
}


MATRIX  *
MatrixAdd( const MATRIX *m1, const MATRIX *m2, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  rows = m1->rows ;
  cols = m1->cols ;

  if ((rows != m2->rows) || (cols != m2->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixAdd: incompatable matrices %d x %d + %d x %d\n",
                 rows, cols, m2->rows, m2->cols)) ;

  if (!mOut)
  {
    mOut = MatrixAlloc(rows, cols, m1->type) ;
    if (!mOut)
      return(NULL) ;
  }

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixAdd: incompatable matrices %d x %d = %d x %d\n",
                 rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      mOut->rptr[row][col] = m1->rptr[row][col] + m2->rptr[row][col] ;
    }
  }

  return(mOut) ;
}


MATRIX  *
MatrixSubtract( const MATRIX *m1, const MATRIX *m2, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  rows = m1->rows ;
  cols = m1->cols ;

  if ((rows != m2->rows) || (cols != m2->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSubtract: incompatable matrices %d x %d - %d x %d\n",
                 rows, cols, m2->rows, m2->cols)) ;

  if (!mOut)
  {
    mOut = MatrixAlloc(rows, cols, m1->type) ;
    if (!mOut)
      return(NULL) ;
  }

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSubtract: incompatable matrices %d x %d = %d x %d\n",
                 rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
      mOut->rptr[row][col] = m1->rptr[row][col] - m2->rptr[row][col] ;
  }

  return(mOut) ;
}


MATRIX  *
MatrixScalarMul( const MATRIX *mIn, const float val, MATRIX *mOut)
{
  int  row, col, rows, cols ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixScalarMul: incompatable matrices %d x %d != %d x %d\n",
                 rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
      mOut->rptr[row][col] = mIn->rptr[row][col] * val ;
  }
  return(mOut) ;
}


MATRIX  *
MatrixClear(MATRIX *mat)
{
  int    rows, row, cols ;

  rows = mat->rows ;
  cols = mat->cols ;
  for (row = 1 ; row <= mat->rows ; row++)
    memset((char *)mat->rptr[row], 0, (cols+1) * sizeof(float)) ;

  return(mat) ;
}


MATRIX *
MatrixSquareElts(MATRIX *mIn, MATRIX *mOut)
{
  int    row, col, rows, cols ;
  float  val ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
      rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      val = mIn->rptr[row][col] ;
      mOut->rptr[row][col] = val * val ;
    }
  }
  return(mOut) ;
}


MATRIX *
MatrixSqrtElts(MATRIX *mIn, MATRIX *mOut)
{
  int    row, col, rows, cols ;
  float  val ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
      rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      val = mIn->rptr[row][col] ;
      mOut->rptr[row][col] = sqrt(val) ;
    }
  }
  return(mOut) ;
}


MATRIX *
MatrixSignedSquareElts(MATRIX *mIn, MATRIX *mOut)
{
  int    row, col, rows, cols ;
  float  val ;

  if (!mOut)
  {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type) ;
    if (!mOut)
      return(NULL) ;
  }

  rows = mIn->rows ;
  cols = mIn->cols ;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
      rows, cols, mOut->rows, mOut->cols)) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      val = mIn->rptr[row][col] ;
      mOut->rptr[row][col] = val * val * (val < 0) ? -1 : 1 ;
    }
  }
  return(mOut) ;
}


MATRIX *
MatrixMakeDiagonal(MATRIX *mSrc, MATRIX *mDst)
{
  int  row, rows, col, cols ;

  if (!mDst)
    mDst = MatrixClone(mSrc) ;

  rows = mSrc->rows ;
  cols = mSrc->cols ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      if (row == col)
        mDst->rptr[row][col] = mSrc->rptr[row][col] ;
      else
        mDst->rptr[row][col] = 0.0f ;
    }
  }
  return(mDst) ;
}


/*
  mDiag is a column vector.
*/
MATRIX  *
MatrixDiag(MATRIX *mDiag, MATRIX *mOut)
{
  int  row, rows, col, cols, nout ;

  rows = mDiag->rows ;
  cols = mDiag->cols ;
  nout = MAX(rows,cols);

  if (!mOut)
  {
    mOut = MatrixAlloc(nout, nout, mDiag->type) ;
    if (!mOut) return(NULL) ;
  }
  else  MatrixClear(mOut) ;

  if ((nout != mOut->rows) || (nout != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixDiag: incompatable matrices %d x %d != %d x %d\n",
                 nout, nout, mOut->rows, mOut->cols)) ;


  if (rows != 1)
  {
    // column vector
    for (row = 1 ; row <= rows ; row++)
      mOut->rptr[row][row] = mDiag->rptr[row][1] ;
  }
  else
  {
    // row vector
    for (col = 1 ; col <= cols ; col++)
      mOut->rptr[col][col] = mDiag->rptr[1][col] ;
  }

  return(mOut) ;
}


MATRIX *
MatrixCopyRegion( const MATRIX *mSrc, MATRIX *mDst,
		  const int start_row,
		  const int start_col,
		  const int rows,
		  const int cols,
		  const int dest_row,
		  const int dest_col)
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col ;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixCopyRegion: bad destination (%d,%d)\n",
                 dest_row, dest_col)) ;

  srows = mSrc->rows ;
  scols = mSrc->cols ;
  end_row = start_row + rows - 1 ;
  end_col = start_col + cols - 1 ;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) ||
      (start_col > scols) || (end_row > srows) || (end_col > scols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyRegion: bad source region (%d,%d) --> (%d,%d)\n",
                 start_row, start_col, end_row, end_col)) ;

  if (!mDst)
    mDst = MatrixAlloc(rows, cols, mSrc->type) ;

  drows = mDst->rows ;
  dcols = mDst->cols ;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixCopyRegion: destination matrix not large enough (%dx%d)\n",
      rows, cols)) ;

  for (drow = dest_row, srow = start_row ; srow <= end_row ; srow++, drow++)
  {
    for (dcol = dest_col, scol = start_col ;
         scol <= end_col ;
         scol++, dcol++)
    {
      switch (mDst->type)
      {
      case MATRIX_REAL:
        *MATRIX_RELT(mDst, drow, dcol) = *MATRIX_RELT(mSrc, srow, scol) ;
        break ;
      case MATRIX_COMPLEX:
        *MATRIX_CELT(mDst, drow, dcol) = *MATRIX_CELT(mSrc,srow,scol);
        break ;
      }
    }
  }
  return(mDst) ;
}


MATRIX *
MatrixCopyRealRegion( const MATRIX *mSrc, MATRIX *mDst,
		      const int start_row,
		      const int start_col,
		      const int rows,
		      const int cols,
		      const int dest_row,
		      const int dest_col)
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col ;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyRealRegion: bad destination (%d,%d)\n",
                 dest_row, dest_col)) ;

  srows = mSrc->rows ;
  scols = mSrc->cols ;
  end_row = start_row + rows - 1 ;
  end_col = start_col + cols - 1 ;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) ||
      (start_col > scols) || (end_row > srows) || (end_col > scols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixCopyRealRegion: bad source region (%d,%d) --> (%d,%d)\n",
      start_row, start_col, end_row, end_col)) ;

  if (!mDst)
    mDst = MatrixAlloc(rows, cols, mSrc->type) ;

  drows = mDst->rows ;
  dcols = mDst->cols ;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixCopyRealRegion: destination matrix not large enough (%dx%d)\n",
      rows, cols)) ;

  for (drow = dest_row, srow = start_row ; srow <= end_row ; srow++, drow++)
  {
    for (dcol = dest_col, scol = start_col ;
         scol <= end_col ;
         scol++, dcol++)
    {
      switch (mDst->type)
      {
      case MATRIX_REAL:
        *MATRIX_RELT(mDst, drow, dcol) =
          MATRIX_CELT_REAL(mSrc, srow, scol) ;
        break ;
      case MATRIX_COMPLEX:
        MATRIX_CELT_IMAG(mDst, drow, dcol) =
          MATRIX_CELT_IMAG(mSrc,srow,scol);
        break ;
      }
    }
  }
  return(mDst) ;
}


MATRIX *
MatrixCopyImagRegion( const MATRIX *mSrc, MATRIX *mDst,
		      const int start_row,
		      const int start_col,
		      const int rows,
		      const int cols,
		      const int dest_row,
		      const int dest_col )
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col ;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyImagRegion: bad destination (%d,%d)\n",
                 dest_row, dest_col)) ;

  srows = mSrc->rows ;
  scols = mSrc->cols ;
  end_row = start_row + rows - 1 ;
  end_col = start_col + cols - 1 ;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) ||
      (start_col > scols) || (end_row > srows) || (end_col > scols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixCopyImagRegion: bad source region (%d,%d) --> (%d,%d)\n",
      start_row, start_col, end_row, end_col)) ;

  if (!mDst)
    mDst = MatrixAlloc(rows, cols, mSrc->type) ;

  drows = mDst->rows ;
  dcols = mDst->cols ;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn
    (NULL,
     (ERROR_BADPARM,
      "MatrixCopyImagRegion: destination matrix not large enough (%dx%d)\n",
      rows, cols)) ;

  for (drow = dest_row, srow = start_row ; srow <= end_row ; srow++, drow++)
  {
    for (dcol = dest_col, scol = start_col ;
         scol <= end_col ;
         scol++, dcol++)
    {
      switch (mDst->type)
      {
      case MATRIX_REAL:
        *MATRIX_RELT(mDst, drow, dcol) =
          MATRIX_CELT_IMAG(mSrc, srow, scol) ;
        break ;
      case MATRIX_COMPLEX:
        MATRIX_CELT_IMAG(mDst, drow, dcol) =
          MATRIX_CELT_IMAG(mSrc,srow,scol);
        break ;
      }
    }
  }
  return(mDst) ;
}


/*--------------------------------------------------------------------*/
MATRIX  *MatrixSetRegion(MATRIX *mSrc, MATRIX *mDst, int start_row,
                         int start_col, int rows, int cols, float val)
{
  int r, c, end_row, end_col;

  end_row = start_row + rows - 1;
  end_col = start_col + cols - 1;

  //printf("rows: %d %d (%d)\n",start_row,end_row,rows);
  //printf("cols: %d %d (%d)\n",start_col,end_col,cols);

  if (start_row > mSrc->rows)
  {
    printf("ERROR: start row > nrows\n");
    return(NULL);
  }
  if (end_row > mSrc->rows)
  {
    printf("ERROR: end row > nrows\n");
    return(NULL);
  }
  if (start_col > mSrc->cols)
  {
    printf("ERROR: start col > ncols\n");
    return(NULL);
  }
  if (end_col > mSrc->cols)
  {
    printf("ERROR: end col > ncols\n");
    return(NULL);
  }

  if (mDst != NULL)
  {
    if (mDst->rows != mSrc->rows || mDst->cols != mSrc->cols)
    {
      printf("ERROR: dimension mismatch\n");
      return(NULL);
    }
  }
  mDst = MatrixCopy(mSrc,mDst);

  for (r = start_row ; r <= end_row; r++)
    for (c = start_col ; c <= end_col; c++)
      mDst->rptr[r][c] = val;

  return(mDst);
}


/*--------------------------------------------------------------------*/
MATRIX *
MatrixRealToComplex(MATRIX *mReal, MATRIX *mImag, MATRIX *mOut)
{
  int rows, cols, row, col ;

  rows = mReal->rows ;
  cols = mReal->cols ;
  if (!mOut)
    mOut = MatrixAlloc(mReal->rows, mReal->cols, MATRIX_COMPLEX) ;

  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      MATRIX_CELT_REAL(mOut,row,col) = *MATRIX_RELT(mReal, row, col) ;
      MATRIX_CELT_IMAG(mOut,row,col) = *MATRIX_RELT(mImag, row, col) ;
    }
  }

  return(mOut) ;
}


float
MatrixDeterminant(MATRIX *mIn)
{
  return OpenMatrixDeterminant( mIn );
}

/* return the eigenvalues and eigenvectors of the symmetric matrix m in
   evalues and m_dst (columns are vectors) respectively. Note that
   the eigenvalues and eigenvectors are sorted in descending order of
   eigenvalue size.
*/


static int compare_evalues(const void *l1, const void *l2)  ;

typedef struct
{
  int   eno ;
  float evalue ;
}
EIGEN_VALUE, EVALUE ;

static int
compare_evalues(const void *l1, const void *l2)
{
  EVALUE *e1, *e2 ;

  e1 = (EVALUE *)l1 ;
  e2 = (EVALUE *)l2 ;
  return(fabs(e1->evalue) < fabs(e2->evalue) ? 1 : -1) ;
}


MATRIX *
MatrixCalculateEigenSystemHelper( MATRIX *m, float *evalues,
                                  MATRIX *m_evectors, int isSymmetric )
{
  int     col, i, nevalues, row ;
  EVALUE  *eigen_values ;
  MATRIX  *mTmp ;

  // sanity-check: input must be n-by-n
  if (m->rows != m->cols) return NULL;

  nevalues = m->rows ;
  eigen_values = (EVALUE *)calloc((UINT)nevalues, sizeof(EIGEN_VALUE));
  if (!m_evectors)
    m_evectors = MatrixAlloc(m->rows, m->cols, MATRIX_REAL) ;

  mTmp = MatrixAlloc(m->rows, m->cols, MATRIX_REAL) ;

  if (isSymmetric)
  {
    if (OpenEigenSystem(m->data, m->rows, evalues, mTmp->data) != NO_ERROR)
      return(NULL) ;
  }
  else
  {
    if (OpenNonSymmetricEigenSystem(m->data, m->rows, evalues, mTmp->data)
        != NO_ERROR)
      return(NULL) ;
  }

  /*
     sort eigenvalues in order of decreasing absolute value. The
     EIGEN_VALUE structure is needed to record final order of eigenvalue so
     we can also sort the eigenvectors.
  */
  for (i = 0 ; i < nevalues ; i++)
  {
    eigen_values[i].eno = i ;
    eigen_values[i].evalue = evalues[i] ;
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues) ;
  for (i = 0 ; i < nevalues ; i++)
    evalues[i] = eigen_values[i].evalue ;

  /* now sort the eigenvectors */
  for (col = 0 ; col < mTmp->cols ; col++)
  {
    for (row = 1 ; row <= mTmp->rows ; row++)
      m_evectors->rptr[row][col+1] =
        mTmp->rptr[row][eigen_values[col].eno+1] ;
  }

  free(eigen_values) ;
  MatrixFree(&mTmp) ;
  return(m_evectors) ;
}


int
MatrixIsSymmetric( MATRIX *matrix )
{
  int row;
  int col;
  int isSymmetric = 1;

  if (matrix->rows != matrix->cols) return 0; // non-square, so not symmetric

  for ( row=1; row<=matrix->rows; row++ )
  {
    for ( col=1; col<=matrix->cols; col++ )
    {
      if ( matrix->rptr[ row ][ col ] != matrix->rptr[ col ][ row ] )
      {
        isSymmetric = 0;
      }
    }
  }

  return isSymmetric;
}


MATRIX *
MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_evectors)
{
  int isSymmetric = MatrixIsSymmetric( m );

  return MatrixCalculateEigenSystemHelper( m, evalues, m_evectors,
         isSymmetric );
}


#if 0
/* m is rows, and n is cols */
static void svd(float **A, float **V, float *z, int m, int n)
{
  int i,j,k,count;
  float c,s,p,q,r,v,toll=0.1,pqrtoll=1;

  // identity_matrix(V,n):
  for (i = 1 ; i <= n ; i++) for (j=1 ; j <= n ; j++)
      V[i][j] = (i==j) ? 1.0 : 0.0 ;

  for (count=n*(n-1)/2;count>0;)
  {
    count=n*(n-1)/2;
    for (j=1;j<=n;j++)
    {
      for (k=j+1;k<=n;k++)
      {
        p=q=r=0;
        for (i=1;i<=m;i++)
        {
          p += A[i][j]*A[i][k];
          q += A[i][j]*A[i][j];
          r += A[i][k]*A[i][k];
        }
        if ((p==0) && (q==0) && (r==0)) break;
        pqrtoll=p*p/(q*r);
        if (pqrtoll==0) break;
        //printf("pqrtoll=%f,p=%f,q=%f,r=%f\n",pqrtoll,p,q,r);
        if ((q*r==0)||(pqrtoll<toll)) count--;
        if (q<r)
        {
          c=0;
          s=1;
        }
        else
        {
          q = q-r;
          v = (float)sqrt((double)(4.0f*p*p+q*q));
          c = (float)sqrt((double)((v+q)/(2.0f*v)));
          s = p/(v*c);
        }
        for (i=1;i<=m;i++)
        {
          r = A[i][j];
          A[i][j] = r*c+A[i][k]*s;
          A[i][k] = -r*s+A[i][k]*c;
        }
        for (i=1;i<=n;i++)
        {
          r = V[i][j];
          V[i][j] = r*c+V[i][k]*s;
          V[i][k] = -r*s+V[i][k]*c;
        }
      }
      if ((p==0) && (q==0) && (r==0)) break;
      if (pqrtoll==0) break;
    }
    if ((p==0) && (q==0) && (r==0)) break;
    if (pqrtoll==0) break;
    //printf("count=%d\n",count);
  }
  for (j=1;j<=n;j++)
  {
    q = 0;
    for (i=1;i<=m;i++) q += A[i][j]*A[i][j];
    q = sqrt(q);
    z[j-1] = q;

    for (i=1;i<=m;i++) A[i][j] /= q;

  }
}
#endif


/* z is an m->cols dimensional vector, mV is an cols x cols dimensional
   identity matrix. Note that all the matrices are (as usual) one based.
*/
/* This looks to be the same as [u s v] = svd(mA), where
   v = mV, and s = diag(v_z) */
MATRIX *MatrixSVD(MATRIX *mA, VECTOR *v_z, MATRIX *mV)
{
  //mV = MatrixIdentity(mA->rows, mV) ;
  //svd(mA->rptr, mV->rptr, v_z->data, mA->rows, mA->cols) ;

  if(mV == NULL) mV = MatrixAlloc(mA->rows,mA->rows,MATRIX_REAL);
  OpenSvdcmp ( mA, v_z, mV ) ;

  return ( mV ) ;
}


float
MatrixSVDEigenValues(MATRIX *m, float *evalues)
{
  float cond ;
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V ;
  int     row, rows, cols, nevalues, i ;
  float   wmax, wmin, wi ;
  EVALUE  *eigen_values ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;
  nevalues = m->rows ;
  memset(evalues, 0, nevalues*sizeof(evalues[0])) ;

  /* calculate condition # of matrix */
  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR)
    return(Gerror) ;

  eigen_values = (EVALUE *)calloc((UINT)nevalues, sizeof(EIGEN_VALUE));
  for (i = 0 ; i < nevalues ; i++)
  {
    eigen_values[i].eno = i ;
    eigen_values[i].evalue = RVECTOR_ELT(v_w, i+1) ;
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues) ;
  for (i = 0 ; i < nevalues ; i++)
    evalues[i] = eigen_values[i].evalue ;

  wmax = 0.0f ;
  wmin = wmax = RVECTOR_ELT(v_w,1) ;
  for (row = 2 ; row <= rows ; row++)
  {
    wi = fabs(RVECTOR_ELT(v_w,row)) ;
    if (wi > wmax)
      wmax = wi ;
    if (wi < wmin)
      wmin = wi ;
  }

  if (FZERO(wmin))
    cond = 1e8 ;  /* something big */
  else
    cond = wmax / wmin ;

  free(eigen_values) ;
  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  return(cond) ;
}


/*
  use SVD to find the inverse of m, usually useful if m is
  ill-conditioned or singular.
*/
#define TOO_SMALL   1e-4

MATRIX *
MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse)
{
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V, *m_w, *m_Ut, *m_tmp ;
  int     row, rows, cols ;
  float   wmax, wmin ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;
  m_w = MatrixAlloc(cols, cols, MATRIX_REAL) ;

  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR)
  {
    MatrixFree(&m_U) ;
    VectorFree(&v_w) ;
    MatrixFree(&m_V) ;
    MatrixFree(&m_w) ;
    return(NULL) ;
  }

  wmax = 0.0f ;
  for (row = 1 ; row <= rows ; row++)
    if (fabs(RVECTOR_ELT(v_w,row)) > wmax)
      wmax = fabs(RVECTOR_ELT(v_w,row)) ;
  wmin = TOO_SMALL * wmax ;
  for (row = 1 ; row <= rows ; row++)
  {
    if (fabs(RVECTOR_ELT(v_w, row)) < wmin)
      m_w->rptr[row][row] = 0.0f ;
    else
      m_w->rptr[row][row] = 1.0f / RVECTOR_ELT(v_w,row) ;
  }

  m_Ut = MatrixTranspose(m_U, NULL) ;
  m_tmp = MatrixMultiply(m_w, m_Ut, NULL) ;
  m_inverse = MatrixMultiply(m_V, m_tmp, m_inverse) ;

  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  MatrixFree(&m_w) ;

  MatrixFree(&m_Ut) ;
  MatrixFree(&m_tmp) ;
  return(m_inverse) ;
}


MATRIX *
MatrixAllocTranslation(int n, double *trans)
{
  MATRIX *mat ;

  mat = MatrixAlloc(n, n, MATRIX_REAL) ;
  return(mat) ;
}


MATRIX *
MatrixReallocRotation(int n, float angle, int which, MATRIX *m)
{
  float  s, c ;

  if (!m)
    m = MatrixIdentity(n, NULL) ;

  c = cos(angle) ;
  s = sin(angle) ;
  switch (which)
  {
  case X_ROTATION:
    m->rptr[2][2] = c ;
    m->rptr[2][3] = s ;
    m->rptr[3][2] = -s ;
    m->rptr[3][3] = c ;
    break ;
  case Y_ROTATION:
    m->rptr[1][1] = c ;
    m->rptr[1][3] = -s ;
    m->rptr[3][1] = s ;
    m->rptr[3][3] = c ;
    break ;
  case Z_ROTATION:
    m->rptr[1][1] = c ;
    m->rptr[1][2] = s ;
    m->rptr[2][1] = -s ;
    m->rptr[2][2] = c ;
    break ;
  }

  return(m) ;
}


MATRIX *
MatrixAllocRotation(int n, float angle, int which)
{
  return(MatrixReallocRotation(n, angle, which, NULL)) ;
}


MATRIX *
MatrixCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans)
{
  int    ninputs, nvars, input, var, var2 ;
  float  *means, covariance, obs1, obs2 ;

  ninputs = mInputs->rows ;   /* number of observations */
  nvars = mInputs->cols ;   /* number of variables */

  if (!ninputs)
  {
    if (!mCov)
      mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;
    return(mCov) ;
  }

  if (!nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixCovariance: zero size input")) ;

  if (!mCov)
    mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;

  if (!mCov)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixCovariance: could not allocate %d x %d covariance matrix",
              nvars, nvars) ;

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCovariance: incorrect covariance matrix dimensions"));

  if (!mMeans)
    means = (float *)calloc(nvars, sizeof(float)) ;
  else
    means = mMeans->data ;

  for (var = 0 ; var < nvars ; var++)
  {
    means[var] = 0.0f ;
    for (input = 0 ; input < ninputs ; input++)
      means[var] += mInputs->rptr[input+1][var+1] ;
    means[var] /= ninputs ;
  }

  for (var = 0 ; var < nvars ; var++)
  {
    for (var2 = 0 ; var2 < nvars ; var2++)
    {
      covariance = 0.0f ;
      for (input = 0 ; input < ninputs ; input++)
      {
        obs1 = mInputs->rptr[input+1][var+1] ;
        obs2 = mInputs->rptr[input+1][var2+1] ;
        covariance += (obs1-means[var]) * (obs2 - means[var2]) ;
      }
      covariance /= ninputs ;
      mCov->rptr[var+1][var2+1]=
        mCov->rptr[var2+1][var+1] = covariance ;
    }
  }

  if (!mMeans)
    free(means) ;
  return(mCov) ;
}


/*
  update the values of a covariance matrix. Note that MatrixFinalCovariance
  must be called with the means and total # of inputs before the covariance
  matrix is valid.
*/
MATRIX *
MatrixUpdateCovariance(MATRIX *mInputs, MATRIX *mCov, MATRIX *mMeans)
{
  int    ninputs, nvars, input, var, var2 ;
  float  covariance, obs1, obs2, mean1, mean2 ;

  ninputs = mInputs->rows ;   /* number of observations */
  nvars = mInputs->cols ;   /* number of variables */

  if (!mMeans || mMeans->rows != nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixUpdateCovariance: wrong mean vector size (%d x %d)",
                 mMeans->rows, mMeans->cols)) ;

  if (!ninputs)
  {
    if (!mCov)
      mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;
    return(mCov) ;
  }

  if (!nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixUpdateCovariance: zero size input")) ;

  if (!mCov)
    mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;

  if (!mCov)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixUpdateCovariance: could not allocate %d x %d covariance "
              "matrix",
              nvars, nvars) ;

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixUpdateCovariance: incorrect covariance matrix "
                 "dimensions"));

  for (var = 0 ; var < nvars ; var++)
  {
    mean1 = mMeans->rptr[var+1][1] ;
    for (var2 = 0 ; var2 < nvars ; var2++)
    {
      mean2 = mMeans->rptr[var2+1][1] ;
      covariance = 0.0f ;
      for (input = 0 ; input < ninputs ; input++)
      {
        obs1 = mInputs->rptr[input+1][var+1] ;
        obs2 = mInputs->rptr[input+1][var2+1] ;
        covariance += (obs1-mean1) * (obs2 - mean2) ;
      }

      mCov->rptr[var+1][var2+1]=
        mCov->rptr[var2+1][var+1] = covariance ;
    }
  }

  return(mCov) ;
}


/*
  update the means based on a new set of observation vectors. Notet that
  the user must keep track of the total # observations
*/
MATRIX *
MatrixUpdateMeans(MATRIX *mInputs, MATRIX *mMeans, VECTOR *mNobs)
{
  int    ninputs, nvars, input, var ;
  float  *means ;

  ninputs = mInputs->rows ;   /* number of observations */
  nvars = mInputs->cols ;     /* number of variables */
  if (!mMeans)
    mMeans = VectorAlloc(nvars, MATRIX_REAL) ;
  means = mMeans->data ;
  for (var = 0 ; var < nvars ; var++)
  {
    mNobs->rptr[var+1][1] += ninputs ;
    for (input = 0 ; input < ninputs ; input++)
      means[var] += mInputs->rptr[input+1][var+1] ;
  }
  return(mMeans) ;
}


/*
  compute the final mean vector previously generated by
  MatrixUpdateMeans. The user must supply a vector containing
  the # of observations of each variable.
*/
MATRIX *
MatrixFinalMeans(VECTOR *mMeans, VECTOR *mNobs)
{
  int    nvars, var, n ;
  float  *means ;

  nvars = mMeans->rows ;     /* number of variables */
  means = mMeans->data ;
  for (var = 0 ; var < nvars ; var++)
  {
    n = mNobs->rptr[var+1][1] ;
    if (!n)
      means[var] = 0.0f ;
    else
      means[var] /= (float)n ;
  }
  return(mMeans) ;
}


/*
  compute the final covariance matrix previously generated by
  MatrixUpdateCovariance. The user must supply a vector containing
  the # of observations per variable.
*/
MATRIX *
MatrixFinalCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mNobs)
{
  int    ninputs, nvars, var, var2 ;
  float  covariance ;

  ninputs = mInputs->rows ;   /* number of observations */
  nvars = mInputs->cols ;   /* number of variables */

  if (!ninputs)
  {
    if (!mCov)
      mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;
    return(mCov) ;
  }

  if (!nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixCovariance: zero size input")) ;

  if (!mCov)
    mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL) ;

  if (!mCov)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixCovariance: could not allocate %d x %d covariance matrix",
              nvars, nvars) ;

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCovariance: incorrect covariance matrix dimensions"));

  for (var = 0 ; var < nvars ; var++)
  {
    ninputs = mCov->rptr[var+1][1] ;
    for (var2 = 0 ; var2 < nvars ; var2++)
    {
      covariance = mCov->rptr[var+1][var2+1] ;
      if (ninputs)
        covariance /= ninputs ;
      else
        covariance = 0.0f ;
      mCov->rptr[var+1][var2+1] =
        mCov->rptr[var2+1][var+1] = covariance ;
    }
  }

  return(mCov) ;
}


int
MatrixAsciiWrite(const char *fname, MATRIX *m)
{
  FILE  *fp ;
  int   ret ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE,
                (ERROR_NO_FILE, "MatrixAsciiWrite: could not open file %s",
                 fname)) ;
  ret = MatrixAsciiWriteInto(fp, m) ;
  fclose(fp) ;
  return(ret) ;
}


MATRIX *
MatrixAsciiRead(const char *fname, MATRIX *m)
{
  FILE  *fp ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NO_FILE, "MatrixAsciiRead: could not open file %s",
                 fname)) ;
  m = MatrixAsciiReadFrom(fp, m) ;
  fclose(fp) ;
  return(m) ;
}


int
MatrixAsciiWriteInto(FILE *fp, MATRIX *m)
{
  int row, col ;

  fprintf(fp, "%d %d %d\n", m->type, m->rows, m->cols) ;
  for (row = 1 ; row <= m->rows ; row++)
  {
    for (col = 1 ; col <= m->cols ; col++)
    {
      if (m->type == MATRIX_COMPLEX)
        fprintf(fp, "%+f %+f   ",
                MATRIX_CELT_REAL(m,row,col), MATRIX_CELT_IMAG(m,row,col));
      else
        fprintf(fp, "%+f  ", m->rptr[row][col]) ;
    }
    fprintf(fp, "\n") ;
  }
  return(NO_ERROR) ;
}


MATRIX *
MatrixAsciiReadFrom(FILE *fp, MATRIX *m)
{
  int    type, rows, cols, row, col ;
  char   *cp, line[200] ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scanf parms"));
  if (sscanf(cp, "%d %d %d\n", &type, &rows, &cols) != 3)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scanf parms"));

  if (!m)
  {
    m = MatrixAlloc(rows, cols, type) ;
    if (!m)
      ErrorReturn(NULL,
                  (ERROR_BADFILE,
                   "MatrixAsciiReadFrom: could not allocate matrix")) ;
  }
  else
  {
    if (m->rows != rows || m->cols != cols || m->type != type)
      ErrorReturn
      (m,
       (ERROR_BADFILE,
        "MatrixAsciiReadFrom: specified matrix does not match file"));
  }
  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
    {
      if (m->type == MATRIX_COMPLEX)
      {
        if (fscanf(fp, "%f %f   ", &MATRIX_CELT_REAL(m,row,col),
                   &MATRIX_CELT_IMAG(m,row,col)) != 2)
          ErrorReturn
          (NULL,
           (ERROR_BADFILE,
            "MatrixAsciiReadFrom: could not scan element (%d, %d)",
            row, col)) ;
      }
      else if (fscanf(fp, "%f  ", &m->rptr[row][col]) != 1)
        ErrorReturn
        (NULL,
         (ERROR_BADFILE,
          "MatrixAsciiReadFrom: could not scan element (%d, %d)",
          row, col)) ;
    }
    fscanf(fp, "\n") ;
  }

  return(m) ;
}


/*
  calculate and return the Euclidean norm of the vector v.
*/
float
VectorLen( const VECTOR *v )
{
  int   i ;
  float len, vi ;

  for (len = 0.0f, i = 1 ; i <= v->rows ; i++)
  {
    vi = v->rptr[i][1] ;
    len += vi*vi ;
  }
  len = sqrt(len) ;
  return(len) ;
}


/*  compute the dot product of 2 vectors */
float
VectorDot( const VECTOR *v1, const VECTOR *v2 )
{
  int   i ;
  float dot ;

  for (dot = 0.0f, i = 1 ; i <= v1->rows ; i++)
    dot += v1->rptr[i][1]*v2->rptr[i][1] ;
  ;

  return(dot) ;
}


/*  compute the dot product of 2 vectors */
float
VectorNormalizedDot(VECTOR *v1, VECTOR *v2)
{
  float   dot, l1, l2 ;

  l1 = VectorLen(v1) ;
  l2 = VectorLen(v2) ;
  dot = VectorDot(v1, v2) ;
  if (FZERO(l1) || FZERO(l2))
    return(0.0f) ;

  return(dot / (l1 * l2)) ;
}


/*
  extract a column of the matrix m and return it in the
  vector v.
*/
VECTOR *
MatrixColumn(MATRIX *m, VECTOR *v, int col)
{
  int row ;

  if (!v)
    v = VectorAlloc(m->rows, MATRIX_REAL) ;

  for (row = 1 ; row <= m->rows ; row++)
    VECTOR_ELT(v,row) = m->rptr[row][col] ;
  return(v) ;
}


/* calcuate and return the Euclidean distance between two vectors */
float
VectorDistance(VECTOR *v1, VECTOR *v2)
{
  int   row ;
  float dist, d ;

  for (dist = 0.0f,row = 1 ; row <= v1->rows ; row++)
  {
    d = VECTOR_ELT(v1,row) - VECTOR_ELT(v2,row) ;
    dist += (d*d) ;
  }

  return(dist) ;
}


/*
  compute the outer product of the vectors v1 and v2.
*/
MATRIX *
VectorOuterProduct(VECTOR *v1, VECTOR *v2, MATRIX *m)
{
  int   row, col, rows, cols ;
  float r ;

  rows = v1->rows ;
  cols = v2->rows ;

  if (rows != cols)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "VectorOuterProduct: v1->rows %d != v2->rows %d",
                       rows, cols)) ;

  if (!m)
    m = MatrixAlloc(rows, cols, MATRIX_REAL) ;

  for (row = 1 ; row <= rows ; row++)
  {
    r = VECTOR_ELT(v1,row) ;
    for (col = 1 ; col <= cols ; col++)
      m->rptr[row][col] = r * VECTOR_ELT(v2,col) ;
  }

  return(m) ;
}


/*
  Add a small random diagonal matrix to mIn to make it
  non-singular.
*/
#define SMALL 1e-4
MATRIX *
MatrixRegularize(MATRIX *mIn, MATRIX *mOut)
{
  int   rows, cols, row ;
  float ran_num ;

  rows = mIn->rows ;
  cols = mIn->cols ;
  if (!mOut)
    mOut = MatrixAlloc(rows, cols, mIn->type) ;

  if (mIn != mOut)
    MatrixCopy(mIn, mOut) ;

#if 0
  ran_num = (SMALL*randomNumber(0.0, 1.0)+SMALL) ;
#else
  ran_num = SMALL ;
#endif

  for (row = 1 ; row <= rows ; row++)
    mOut->rptr[row][row] += ran_num ;

  return(mOut) ;
}


/* see if a matrix is singular */
int
MatrixSingular(MATRIX *m)
{
#if 1
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V ;
  int     row, rows, cols ;
  float   wmax, wmin, wi ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;

  /* calculate condition # of matrix */
  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR)
    return(Gerror) ;

  wmax = 0.0f ;
  wmin = wmax = RVECTOR_ELT(v_w,1) ;
  for (row = 2 ; row <= rows ; row++)
  {
    wi = fabs(RVECTOR_ELT(v_w,row)) ;
    if (wi > wmax)
      wmax = wi ;
    if (wi < wmin)
      wmin = wi ;
  }

  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  return(FZERO(wmin) ? 1 : wmin < wmax * TOO_SMALL) ;
#else
  float det ;

  det = MatrixDeterminant(m) ;
  return(FZERO(det)) ;
#endif
}


/*
  calcluate the condition # of a matrix using svd
*/
float MatrixConditionNumber(MATRIX *m)
{
  float cond ;
  VECTOR  *v_w ;
  MATRIX  *m_U, *m_V ;
  int     row, rows, cols ;
  float   wmax, wmin, wi ;

  cols = m->cols ;
  rows = m->rows ;
  m_U = MatrixCopy(m, NULL) ;
  v_w = RVectorAlloc(cols, MATRIX_REAL) ;
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;

  /* calculate condition # of matrix */
  OpenSvdcmp(m_U, v_w, m_V) ;
  wmax = 0.0f ;
  wmin = wmax = RVECTOR_ELT(v_w,1) ;
  for (row = 2 ; row <= rows ; row++)
  {
    wi = fabs(RVECTOR_ELT(v_w,row)) ;
    if (wi > wmax)
      wmax = wi ;
    if (wi < wmin)
      wmin = wi ;
  }

  if (FZERO(wmin))
    cond = 1e8 ;   /* something big */
  else
    cond = wmax / wmin ;

  MatrixFree(&m_U) ;
  VectorFree(&v_w) ;
  MatrixFree(&m_V) ;
  return(cond) ;
}


/*-----------------------------------------------------------
  MatrixNSConditionNumber() - condition of a non-square matrix.
  Works for square matrices as well.
  -----------------------------------------------------------*/
float MatrixNSConditionNumber(MATRIX *m)
{
  float cond ;
  MATRIX *mt, *p;

  if (m->rows == m->cols)
    return(MatrixConditionNumber(m));

  mt = MatrixTranspose(m,NULL);

  if (m->rows > m->cols)
    p = MatrixMultiply(mt,m,NULL);
  else
    p = MatrixMultiply(m,mt,NULL);

  cond = MatrixConditionNumber(p);

  MatrixFree(&mt);
  MatrixFree(&p);

  return(cond);
}


/*
  calculate the cross product of two vectors and return the result
  in vdst, allocating it if necessary.
*/
VECTOR *
VectorCrossProduct( const VECTOR *v1, const VECTOR *v2, VECTOR *vdst)
{
  float   x1, x2, y1, y2, z1, z2 ;

  if (v1->rows != 3 && v1->cols != 1)
    ErrorReturn(NULL,(ERROR_BADPARM, "VectorCrossProduct: must be 3-vectors"));

  if (!vdst)
    vdst = VectorClone(v1) ;

  x1 = VECTOR_ELT(v1,1) ;
  y1 = VECTOR_ELT(v1,2) ;
  z1 = VECTOR_ELT(v1,3) ;
  x2 = VECTOR_ELT(v2,1) ;
  y2 = VECTOR_ELT(v2,2) ;
  z2 = VECTOR_ELT(v2,3) ;

  VECTOR_ELT(vdst,1) = y1*z2 - z1*y2 ;
  VECTOR_ELT(vdst,2) = z1*x2 - x1*z2 ;
  VECTOR_ELT(vdst,3) = x1*y2 - y1*x2 ;
  return(vdst) ;
}


/*
  compute the triple scalar product v1 x v2 . v3
*/
float
VectorTripleProduct( const VECTOR *v1, const VECTOR *v2, const VECTOR *v3)
{
  float   x1, x2, y1, y2, z1, z2, x3, y3, z3, total ;

  if (v1->rows != 3 && v1->cols != 1)
    ErrorReturn(0.0f,(ERROR_BADPARM, "VectorCrossProduct: must be 3-vectors"));

  x1 = VECTOR_ELT(v1,1) ;
  y1 = VECTOR_ELT(v1,2) ;
  z1 = VECTOR_ELT(v1,3) ;
  x2 = VECTOR_ELT(v2,1) ;
  y2 = VECTOR_ELT(v2,2) ;
  z2 = VECTOR_ELT(v2,3) ;
  x3 = VECTOR_ELT(v3,1) ;
  y3 = VECTOR_ELT(v3,2) ;
  z3 = VECTOR_ELT(v3,3) ;

  total =  x3 * (y1*z2 - z1*y2) ;
  total += y3 * (z1*x2 - x1*z2) ;
  total += z3 * (x1*y2 - y1*x2) ;
  return(total) ;
}


VECTOR *
VectorNormalize( const VECTOR *vin, VECTOR *vout)
{
  float  len ;
  int    row, col, rows, cols ;

  if (!vout)
    vout = VectorClone(vin) ;

  len = VectorLen(vin) ;
  if (FZERO(len))
    len = 1.0f ;   /* doesn't matter - all elements are 0 */

  rows = vin->rows ;
  cols = vin->cols ;
  for (row = 1 ; row <= rows ; row++)
  {
    for (col = 1 ; col <= cols ; col++)
      vout->rptr[row][col] = vin->rptr[row][col] / len ;
  }
  return(vout) ;
}


float
VectorAngle( const VECTOR *v1, const VECTOR *v2 )
{
  float  angle, l1, l2, dot, norm ;

  l1 = VectorLen(v1) ;
  l2 = VectorLen(v2) ;
  norm = fabs(l1*l2) ;
  if (FZERO(norm))
    return(0.0f) ;
  dot = VectorDot(v1, v2) ;
  if (dot > norm)
    angle = acos(1.0) ;
  else
    angle = acos(dot / norm) ;
  return(angle) ;
}


double
Vector3Angle(VECTOR *v1, VECTOR *v2)
{
  double  angle, l1, l2, dot, norm, x, y, z ;

  x = V3_X(v1) ;
  y = V3_Y(v1) ;
  z = V3_Z(v1) ;
  l1 = sqrt(x*x+y*y+z*z) ;
  x = V3_X(v2) ;
  y = V3_Y(v2) ;
  z = V3_Z(v2) ;
  l2 = sqrt(x*x+y*y+z*z) ;
  norm = l1*l2 ;
  if (FZERO(norm))
    return(0.0f) ;
  dot = V3_DOT(v1, v2) ;
  if (fabs(dot) > fabs(norm))
    norm = fabs(dot) ;
  if (dot > norm)
    angle = acos(1.0) ;
  else
    angle = acos(dot / norm) ;
  return(angle) ;
}


int
MatrixWriteTxt(const char *fname, MATRIX *mat)
{
  FILE   *fp ;
  int     row, col ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NO_FILE,
                (ERROR_NO_FILE, "MatrixWriteTxt(%s) - file open failed\n",
                 fname));

  for (row = 1 ; row <= mat->rows ; row++)
  {
    for (col = 1 ; col <= mat->cols ; col++)
      fprintf(fp, "%+4.5f ", mat->rptr[row][col]) ;
    fprintf(fp, "\n") ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}


MATRIX *
MatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  MATRIX  *mT, *mTm, *mTm_inv ;

  if (m->rows < m->cols)
  {
    MATRIX  *mT = MatrixTranspose(m, NULL) ;
    m_pseudo_inv = MatrixPseudoInverse(mT, m_pseudo_inv) ;
    MatrixFree(&mT) ;
    mT = m_pseudo_inv ;
    m_pseudo_inv = MatrixTranspose(mT, NULL) ;
    MatrixFree(&mT) ;
    return(m_pseudo_inv) ;
  }
  mT = MatrixTranspose(m, NULL) ;

  /* build (mT m)-1 mT */
  mTm = MatrixMultiply(mT, m, NULL) ;
  mTm_inv = MatrixInverse(mTm, NULL) ;
  if (!mTm_inv)
  {
    mTm_inv = MatrixSVDInverse(mTm, NULL) ;

    if (!mTm_inv)
    {
      MatrixFree(&mT) ;
      MatrixFree(&mTm) ;
      return(NULL) ;
    }
  }
  m_pseudo_inv = MatrixMultiply(mTm_inv, mT, m_pseudo_inv) ;

  MatrixFree(&mT) ;
  MatrixFree(&mTm) ;
  MatrixFree(&mTm_inv) ;
  return(m_pseudo_inv) ;
}


MATRIX *
MatrixRightPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  MATRIX  *mT, *mmT, *mmT_inv ;

  /* build mT (m mT)-1 */
  mT = MatrixTranspose(m, NULL) ;
  mmT = MatrixMultiply(m, mT, NULL) ;
  mmT_inv = MatrixInverse(mmT, NULL) ;
  if (!mmT_inv)
  {
    MatrixFree(&mT) ;
    MatrixFree(&mmT) ;
    return(NULL) ;
  }
  m_pseudo_inv = MatrixMultiply(mT, mmT_inv, m_pseudo_inv) ;

  MatrixFree(&mT) ;
  MatrixFree(&mmT) ;
  MatrixFree(&mmT_inv) ;
  return(m_pseudo_inv) ;
}


int
MatrixCheck(MATRIX *m)
{
  int  rows,  cols, r, c ;

  rows = m->rows ;
  cols = m->cols ;
  for (r = 1 ; r <= rows ; r++)
  {
    for (c = 1 ; c <= cols ; c++)
    {
      if (!finite(*MATRIX_RELT(m, r, c)))
        return(ERROR_BADPARM) ;
    }
  }
  return(NO_ERROR) ;
}


MATRIX  *
MatrixReshape(MATRIX *m_src, MATRIX *m_dst, int rows, int cols)
{
  int   r1, c1, r2, c2 ;

  if (m_dst)
  {
    rows = m_dst->rows ;
    cols = m_dst->cols ;
  }

#if 0
  if (rows*cols > m_src->rows*m_src->cols)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MatrixReshape: (%d,%d) -> (%d,%d), lengths must be"
                       " equal", m_src->rows, m_src->cols, rows,cols)) ;
#endif
  if (!m_dst)
    m_dst = MatrixAlloc(rows, cols, m_src->type) ;

  for (r2 = c2 = r1 = 1 ; r1 <= m_src->rows ; r1++)
  {
    for (c1 = 1 ; c1 <= m_src->cols ; c1++)
    {
      *MATRIX_RELT(m_dst, r2, c2) = *MATRIX_RELT(m_src, r1, c1) ;
      if (++c2 > m_dst->cols)
      {
        c2 = 1 ;
        r2++ ;
      }
      if (r2 > rows)  /* only extract first rowsxcols elements */
        break ;
    }
    if (r2 > rows)
      break ;
  }

  return(m_dst) ;
}


/*------------------------------------------------------------*/
float MatrixTrace(MATRIX *M)
{
  int n, nmax;
  float trace;

  if (M->rows >= M->cols) nmax = M->cols;
  else                   nmax = M->rows;

  trace = 0.0;
  for (n=1; n <= nmax; n++) trace += M->rptr[n][n];

  return(trace);
}


/*------------------------------------------------------------------
  MatrixHorCat() - horizontally concatenate matrices m1 and m2, store
  the result in mcat. If mcat is NULL, a new matrix is allocated. A
  pointer to mcat (or the new matrix) is returned.  If m1 is NULL, m2
  is copied into mcat. If m2 is NULL, m1 is copied into mcat.
  -------------------------------------------------------------------*/
MATRIX *MatrixHorCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat)
{
  int r,c1,c2,c;

  if (m1==NULL && m2==NULL)
  {
    printf("ERROR: MatrixHorCat: both m1 and m2 are NULL\n");
    return(NULL);
  }

  if (m1==NULL)
  {
    mcat = MatrixCopy(m2,mcat);
    return(mcat);
  }

  if (m2==NULL)
  {
    mcat = MatrixCopy(m1,mcat);
    return(mcat);
  }

  if (m1->rows != m2->rows)
  {
    printf("ERROR: MatrixHorCat: rows of m1 (%d) not equal to m2 (%d)\n",
           m1->rows, m2->rows);
    return(NULL);
  }

  if (mcat == NULL)
    mcat = MatrixAlloc(m1->rows,m1->cols+m2->cols,MATRIX_REAL);
  else
  {
    if (mcat->rows != m1->rows)
    {
      printf("ERROR: MatrixHorCat: rows of m1 (%d) not equal to mcat (%d)\n",
             m1->rows, mcat->rows);
      return(NULL);
    }
    if (mcat->cols != (m1->cols+m2->cols))
    {
      printf
      ("ERROR: MatrixHorCat: cols of mcat (%d) do not equal to m1+m2 (%d)\n",
       mcat->cols, m1->cols+m2->cols);
      return(NULL);
    }
  }

  /* Fill in the horizontally concatenated matrix */
  for (r=0; r < mcat->rows; r++)
  {
    c = 0;
    for (c1=0; c1 < m1->cols; c1++)
    {
      mcat->rptr[r+1][c+1] = m1->rptr[r+1][c1+1];
      c++;
    }
    for (c2=0; c2 < m2->cols; c2++)
    {
      mcat->rptr[r+1][c+1] = m2->rptr[r+1][c2+1] ;
      c++;
    }
  }

  return(mcat);
}


/*------------------------------------------------------------------
  MatrixVertCat() - vertically concatenate matrices m1 and m2, store
  the result in mcat. If mcat is NULL, a new matrix is allocated. A
  pointer to mcat (or the new matrix) is returned.  If m1 is NULL, m2
  is copied into mcat. If m2 is NULL, m1 is copied into mcat.
  -------------------------------------------------------------------*/
MATRIX *MatrixVertCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat)
{
  int r,r1,r2,c;

  if (m1==NULL && m2==NULL)
  {
    printf("ERROR: MatrixVertCat: both m1 and m2 are NULL\n");
    return(NULL);
  }

  if (m1==NULL)
  {
    mcat = MatrixCopy(m2,mcat);
    return(mcat);
  }

  if (m2==NULL)
  {
    mcat = MatrixCopy(m1,mcat);
    return(mcat);
  }

  if (m1->cols != m2->cols)
  {
    printf("ERROR: MatrixVertCat: cols of m1 (%d) not equal to m2 (%d)\n",
           m1->cols, m2->cols);
    return(NULL);
  }

  if (mcat == NULL)
    mcat = MatrixAlloc(m1->rows+m2->rows,m1->cols,MATRIX_REAL);
  else
  {
    if (mcat->cols != m1->cols)
    {
      printf("ERROR: MatrixVertCat: cols of m1 (%d) not equal to mcat (%d)\n",
             m1->cols, mcat->cols);
      return(NULL);
    }
    if (mcat->rows != (m1->rows+m2->rows))
    {
      printf
      ("ERROR: MatrixVertCat: rows of mcat (%d) "
       "do not equal to m1+m2 (%d)\n",
       mcat->rows, m1->rows+m2->rows);
      return(NULL);
    }
  }

  /* Fill in the vertically concatenated matrix */
  for (c=0; c < mcat->cols; c++)
  {
    r = 0;
    for (r1=0; r1 < m1->rows; r1++)
    {
      mcat->rptr[r+1][c+1] = m1->rptr[r1+1][c+1] ;
      r++;
    }
    for (r2=0; r2 < m2->rows; r2++)
    {
      mcat->rptr[r+1][c+1] = m2->rptr[r2+1][c+1] ;
      r++;
    }

  }

  return(mcat);
}


/*-------------------------------------------------------------
  MatrixConstVal() - sets all the elements to the given value.
  If X is NULL, then a matrix rows-by-cols is alloced. If X
  is non-NULL, then rows and cols are ignored.
  -------------------------------------------------------------*/
MATRIX *MatrixConstVal(float val, int rows, int cols, MATRIX *X)
{
  int r,c;

  if (X==NULL) X = MatrixAlloc(rows,cols,MATRIX_REAL);

  for (r=1; r <= X->rows; r++)
  {
    for (c=1; c <= X->cols; c++)
    {
      X->rptr[r][c] = val;
    }
  }

  return(X);
}


/*--------------------------------------------------------------------
  MatrixZero() - sets all the elements to zero.  If X is NULL, then a
  matrix rows-by-cols is alloced. If X is non-NULL, then rows and cols
  are ignored.
  ------------------------------------------------------------------*/
MATRIX *MatrixZero(int rows, int cols, MATRIX *X)
{
  X = MatrixConstVal(0,rows,cols,X);
  return(X);
}


/*----------------------------------------------------------------*/
MATRIX *MatrixSum(MATRIX *m, int dim, MATRIX *msum)
{
  int outrows, outcols;
  int r,c;

  if (dim==1)
  { /* sum over the rows */
    outrows = 1;
    outcols = m->cols;
  }
  else
  {       /* sum over the cols */
    outrows = m->rows;
    outcols = 1;
  }

  if (msum == NULL) msum = MatrixZero(outrows,outcols,NULL);
  else
  {
    if (msum->rows != outrows || msum->cols != outcols)
    {
      printf("ERROR: MatrixSum: dimension mismatch\n");
      return(NULL);
    }
  }

  if ( (dim ==1 && m->rows > 1) || (dim == 2 && m->cols > 1) )
  {
    for (r=1; r <= m->rows; r++)
    {
      for (c=1; c <= m->cols; c++)
      {
        if (dim==1) msum->rptr[1][c] += m->rptr[r][c];
        else       msum->rptr[r][1] += m->rptr[r][c];
      }
    }
  }
  else
  { /* Just copy vector to output */
    if (dim==1)
      for (c=1; c <= m->cols; c++) msum->rptr[1][c] = m->rptr[1][c];
    else
      for (r=1; r <= m->rows; r++) msum->rptr[r][1] = m->rptr[r][1];
  }

  return(msum);
}
/*----------------------------------------------------------------*/
double MatrixSumElts(MATRIX *m)
{
  int    r,c;
  double msum ;

  for (msum = 0.0, r=1; r <= m->rows; r++)
    for (c=1; c <= m->cols; c++)
      msum += m->rptr[r][c];

  return(msum);
}


/*----------------------------------------------------------------*/
double VectorRange(MATRIX *v, double *pVmin, double *pVmax)
{
  double min,max, val;
  int r,c;

  min = v->rptr[1][1];
  max = v->rptr[1][1];
  for (r=1; r <= v->rows; r++)
  {
    for (c=1; c <= v->cols; c++)
    {
      val = v->rptr[r][c];
      if (min > val) min = val;
      if (max < val) max = val;
    }
  }

  if (pVmin != NULL) *pVmin = min;
  if (pVmax != NULL) *pVmax = max;

  return(max-min);
}


/*----------------------------------------------------------------*/
double VectorSum(MATRIX *v)
{
  double sum;
  int r,c;

  sum = 0.0;
  for (r=1; r <= v->rows; r++)
    for (c=1; c <= v->cols; c++)
      sum += (v->rptr[r][c]);
  return(sum);
}


/*----------------------------------------------------------------*/
double VectorMean(MATRIX *v)
{
  double sum,mean;

  sum = VectorSum(v);
  mean = sum/(v->rows * v->cols);

  return(mean);
}


/*----------------------------------------------------------------*/
double VectorVar(MATRIX *v, double *pMean)
{
  double f, mean, sum2, var;
  int r,c;

  mean = VectorMean(v);
  if (pMean != NULL) *pMean = mean;

  sum2 = 0.0;
  for (r=1; r <= v->rows; r++)
  {
    for (c=1; c <= v->cols; c++)
    {
      f = v->rptr[r][c] - mean;
      sum2 += (f*f);
    }
  }

  var = sum2/(v->rows * v->cols - 1);

  return(var);
}


/*----------------------------------------------------------------*/
double VectorStdDev(MATRIX *v, double *pMean)
{
  double  var;
  var = VectorVar(v,pMean);
  return(sqrt(var));
}


/*----------------------------------------------------------------*/
MATRIX *MatrixDRand48(int rows, int cols, MATRIX *m)
{
  int r,c;

  if (m==NULL) m = MatrixAlloc(rows,cols,MATRIX_REAL);
  /* if m != NULL rows and cols are ignored */

  for (r=1; r <= m->rows; r++)
  {
    for (c=1; c <= m->cols; c++)
    {
      m->rptr[r][c] = drand48();
    }
  }

  return(m);
}


/*--------------------------------------------------------------------
  MatrixFactorSqrSVD() - factors a square matrix M into D such that M
  = D'D. This is like a cholesky factorization, except that D will be
  non-causal. M must be positive semidefinite (ie, all the eigenvalues
  must be greater than zero). If the Invert flag is non-zero, then M =
  inv(D'D), in which case M must be positive definite. D will be the
  same size as M. If D is NULL, it will be allocated.

  Internally, intermediate matrices are held static and reallocated
  only if the dimensions of M change from call to call. This makes the
  code more efficient.

  Return: D
  -----------------------------------------------------------------*/
MATRIX *MatrixFactorSqrSVD(MATRIX *M, int Invert, MATRIX *D)
{
  static int matalloc = 1;
  static VECTOR *S = NULL;
  static MATRIX *U = NULL, *V = NULL, *Vt = NULL;
  static int Mrows = -1, Mcols = -1;
  float s;
  int r,c;

  if (M->rows != M->cols)
  {
    printf("ERROR: MatrixFactorSqrSVD: matrix is not square\n");
    return(NULL);
  }

  if (D == NULL)
    D = MatrixAlloc(M->rows,M->cols,MATRIX_REAL);
  else
  {
    if (D->rows != M->rows || D->cols != M->cols)
    {
      printf("ERROR: MatrixFactorSqrSVD: dimension mismatch\n");
      return(NULL);
    }
  }

  if (!matalloc)
  {
    /* Check that M is the same size as last time */
    if (Mrows != M->rows || Mcols != M->cols)
    {
      MatrixFree(&U);
      MatrixFree(&S);
      MatrixFree(&V);
      MatrixFree(&Vt);
      matalloc = 1;
    }
  }
  if (matalloc)
  {
    /* Allocate intermediate matrices */
    U  = MatrixAlloc(M->rows,M->cols,MATRIX_REAL);
    S  = MatrixAlloc(M->rows,1,MATRIX_REAL);
    V  = MatrixAlloc(M->rows,M->cols,MATRIX_REAL);
    Vt = MatrixAlloc(M->rows,M->cols,MATRIX_REAL);
  }

  /* Make a local copy of M because SVD alters it */
  MatrixCopy(M,U);

  /* Compute SVD of matrix [u s v'] = svd(M)*/
  if (MatrixSVD(U, S, V) == NULL)
  {
    printf("ERROR: MatrixFactorSqrSVD: cound not SVD\n");
    return(NULL);
  }

  /* Check for pos/posdef. Invert if necessary.*/
  for (r=1; r <= S->rows; r++)
  {
    s = S->rptr[r][1];
    if (!Invert && s < 0 )
    {
      printf("ERROR: MatrixFactorSqrSVD: matrix is not positive.\n");
      return(NULL);
    }
    if (Invert && s <= 0 )
    {
      printf("ERROR: MatrixFactorSqrSVD: matrix is not positive definite.\n");
      return(NULL);
    }
    if (Invert) s = 1.0/s;
    S->rptr[r][1] = sqrt(s);
  }

  MatrixTranspose(V,Vt);

  /* D = U*diag(S)*V' */
  /* U = U*diag(S): Multiply each column of U by S. */
  for (c=1; c <= U->cols; c++)
    for (r=1; r <= U->rows; r++)
      U->rptr[r][c] *= S->rptr[c][1];

  /* D = U*Vt = U*diag(S)*V'*/
  MatrixMultiply(U,Vt,D);

  matalloc = 0;
  Mrows = M->rows;
  Mcols = M->cols;
  return(D);
}


/*-------------------------------------------------------------
  MatrixToeplitz() - creates a toeplitz matrix T (ie, a matrix
  where all the elements on a diagonal have the same value.
  The value on the ith diagonal equals the ith value of v.
  T will be square in the number of rows of v. There are three
  subtypes of matrices possible: (1) Lower: only diagonals
  below and including the main, (2) Upper: only diagonals
  above and including the main, and (3) Symetric: both upper
  and lower.
  ---------------------------------------------------------*/
MATRIX *MatrixToeplitz(VECTOR *v, MATRIX *T, int Type)
{
  int r,c;

  if (Type != MATRIX_SYM && Type != MATRIX_UPPER && Type != MATRIX_LOWER)
  {
    printf("ERROR: Type = %d, unrecognized\n",Type);
    return(NULL);
  }

  if (T==NULL)
  {
    T = MatrixAlloc(v->rows,v->rows,MATRIX_REAL);
    if (T==NULL) return(NULL);
  }
  else
  {
    if (T->rows != v->rows || T->cols != v->rows)
    {
      printf("ERROR: dimension mismatch\n");
      return(NULL);
    }
  }

  for (r = 1; r <= T->rows; r++)
  {
    for (c = 1; c <= T->cols; c++)
    {
      if (r==c)
      {
        T->rptr[r][c] = v->rptr[1][1];
        continue;
      }
      if ( (r > c && Type == MATRIX_UPPER) ||
           (r < c && Type == MATRIX_LOWER))
      {
        T->rptr[r][c] = 0.0;
        continue;
      }
      T->rptr[r][c] = v->rptr[abs(r-c)+1][1];
    }
  }

  return(T);
}


/*!
\fn MATRIX *MatrixNormalizeColScale(MATRIX *m, MATRIX *scale)
\brief Computes the scaling used for MatrixNormalizeCol()
*/
MATRIX *MatrixNormalizeColScale(MATRIX *m, MATRIX *scale)
{
  int r,c;
  double sum2, v;

  if(scale == NULL) scale = MatrixAlloc(1,m->cols,MATRIX_REAL);
  for (c=1; c <= m->cols; c++) {
    sum2 = 0.0;
    for (r=1; r <= m->rows; r++) {
      v = m->rptr[r][c];
      sum2 += (v*v);
    }
    scale->rptr[1][c] = sqrt(sum2);
  }
  //printf("scale -----------------------------\n");
  //MatrixPrint(stdout,scale);
  //printf("------------------------------------\n");
  return(scale);
}

/*!
\fn MATRIX *MatrixNormalizeCol(MATRIX *m, MATRIX *mcnorm)
\brief Rescales m so that sum(col^2)=1
*/
MATRIX *MatrixNormalizeCol(MATRIX *m, MATRIX *mcnorm, MATRIX *scale)
{
  int r,c,FreeScale=1;
  double v;

  if (mcnorm == NULL)
  {
    mcnorm = MatrixAlloc(m->rows,m->cols,MATRIX_REAL);
    if (mcnorm == NULL)
    {
      printf("ERROR: MatrixNormalizeCol: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    if (mcnorm->rows != m->rows || mcnorm->cols != m->cols)
    {
      printf("ERROR: MatrixNormalizeCol: dimension mismatch\n");
      return(NULL);
    }
  }

  if(scale) FreeScale = 0;
  scale = MatrixNormalizeColScale(m,scale);
  for (c=1; c <= m->cols; c++)
  {
    v = scale->rptr[1][c];
    if (v != 0)
      for (r=1; r <= m->rows; r++)
        mcnorm->rptr[r][c] = (m->rptr[r][c])/v;
    else
      for (r=1; r <= m->rows; r++)
        mcnorm->rptr[r][c] = 0.0;

  }
  if(FreeScale) MatrixFree(&scale);

  //printf("m ----------------------------\n");
  //MatrixPrint(stdout,m);
  //printf("mcnorm ----------------------------\n");
  //MatrixPrint(stdout,mcnorm);

  return(mcnorm);
}


MATRIX *
MatrixSimilarityTransform(MATRIX *m_src, MATRIX *m_mul, MATRIX *m_dst)
{
  MATRIX *m_mul_T, *m_tmp ;

  m_mul_T = MatrixTranspose(m_mul, NULL) ;
  m_tmp = MatrixMultiply(m_src, m_mul_T, NULL) ;
  m_dst = MatrixMultiply(m_mul, m_tmp, m_dst) ;
  MatrixFree(&m_mul_T) ;
  MatrixFree(&m_tmp) ;
  return(m_dst) ;
}


/*---------------------------------------------------------------
  GaussianVector() - creates a vector with len rows. The gaussian
  curve is centered at the meanth row and has std. The mean can
  be non-integer. If norm == 1, then the sum is adjusted to be 1.
  ---------------------------------------------------------------*/
MATRIX *GaussianVector(int len, float mean, float std, int norm,
                       MATRIX *g)
{
  int n;
  float v,sum,var,f;

  if (g == NULL) g = MatrixAlloc(len,1,MATRIX_REAL);
  else
  {
    if (g->rows != len)
    {
      printf("ERROR: GaussianVector: dimension mismatch\n");
      return(NULL);
    }
  }

  var = std*std;
  f = 2*M_PI*std;
  sum = 0;
  for (n=0; n < len; n++)
  {
    v = exp( -(n-mean)*(n-mean)/(2*var) )/f;
    if (norm) sum += v;
    g->rptr[n+1][1] = v;
  }

  if (norm)
  {
    for (n=0; n < len; n++)
    {
      g->rptr[n+1][1] /= sum;
    }
  }

  return(g);
}


/*---------------------------------------------------------------
  GaussianMatrix() - creates a gaussian convolution matrix. Each row
  is a gaussian waveform with centered at the diagonal with standard
  deviation std.  If norm == 1, then the sum of each row is adjusted
  to be 1. The matrix will be len-by-len.
  ---------------------------------------------------------------*/
MATRIX *GaussianMatrix(int len, float std, int norm, MATRIX *G)
{
  int r, c;
  float d,v,sum,var,f;

  if (G == NULL) G = MatrixAlloc(len,len,MATRIX_REAL);
  else
  {
    if (G->rows != len || G->cols != len)
    {
      printf("ERROR: GaussianMatrix: dimension mismatch\n");
      return(NULL);
    }
  }

  var = std*std;
  f = sqrt(2*M_PI)*std;

  // Adjust scale so that sum at center line is 1
  if (norm)
  {
    sum = 0;
    for (c=0; c < len; c++)
    {
      d = c-len/2;
      v = exp( -(d*d)/(2*var) )/f;
      sum += v;
    }
    f /= sum;
  }

  for (r=0; r<len; r++)
  {
    for (c=0; c < len; c++)
    {
      d = c-r;
      v = exp( -(d*d)/(2*var) )/f;
      G->rptr[r+1][c+1] = v;
    }
  }

  //printf("std = %g\n",std);
  //MatrixWriteTxt("g.txt",G);
  //exit(1);

  return(G);
}


MATRIX *
MatrixSVDPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  int rows = m->rows ;
  int cols = m->cols ;
  if (rows < cols)
  {
    MATRIX  *mT = MatrixTranspose(m, NULL) ;
    m_pseudo_inv = MatrixSVDPseudoInverse(mT, m_pseudo_inv) ;
    MatrixFree(&mT) ;
    mT = m_pseudo_inv ;
    m_pseudo_inv = MatrixTranspose(mT, NULL) ;
  }
  else
  {
    int r,c ;
    MATRIX  *m_U, *m_V, *m_Ur, *m_Vr, *m_Sr, *m_tmp, *m_S ;
    VECTOR  *v_S ;

    m_U = MatrixCopy(m, NULL) ;
    m_V = MatrixAlloc(cols, cols, MATRIX_REAL) ;
    v_S = VectorAlloc(cols, MATRIX_REAL) ;

    OpenSvdcmp(m_U, v_S, m_V) ;

    for (r = 1 ; r <= v_S->rows ; r++)
      if (VECTOR_ELT(v_S,r)/VECTOR_ELT(v_S,1) < 1e-4)
        break ;
    r-- ; // previous one was last non-zero
    m_tmp = MatrixCopyRegion(m_U, NULL, 1, 1, m_U->rows, r, 1, 1) ;
    MatrixFree(&m_U);
    m_Ur = MatrixTranspose(m_tmp, NULL) ;
    MatrixFree(&m_tmp) ;
    m_S = MatrixDiag(v_S, NULL) ;
    VectorFree(&v_S) ;
    m_Sr = MatrixCopyRegion(m_S, NULL, 1, 1, r, r, 1, 1) ;
    MatrixFree(&m_S) ;
    for (c = 1 ; c <= m_Sr->rows ; c++)
      *MATRIX_RELT(m_Sr, c, c) = 1 / *MATRIX_RELT(m_Sr, c, c) ;
    m_Vr = MatrixCopyRegion(m_V, NULL, 1, 1, m_V->rows, r, 1, 1) ;
    MatrixFree(&m_V);
    m_tmp = MatrixMultiply(m_Vr, m_Sr, NULL) ;
    MatrixFree(&m_Sr) ;
    MatrixFree(&m_Vr) ;
    m_pseudo_inv = MatrixMultiply(m_tmp, m_Ur, NULL) ;
    MatrixFree(&m_tmp) ;
    MatrixFree(&m_Ur) ;
  }

  return(m_pseudo_inv) ;
}


/*-----------------------------------------------------------------
  MatrixReorderRows() - reorders the rows of a matrix to be that
  indicated by NewRowOrder. The first row of the output matrix
  will be the row of the input matrix indicated by NewRowOrder[0].
  Eg, if NewRowOrder[0]=3, then the 1st row of XRO will be the
  3rd row of X. The rows in NewRowOrder are 1-based. See also
  RandPermList() in evschutils.c. CANNOT be done in-place.
  -----------------------------------------------------------------*/
MATRIX *MatrixReorderRows(MATRIX *X, int *NewRowOrder, MATRIX *XRO)
{
  int r, c;
  if (XRO == NULL) XRO = MatrixAlloc(X->rows,X->cols,MATRIX_REAL);
  for (r=1; r <= X->rows; r++)
  {
    for (c=1; c <= X->cols; c++)
    {
      XRO->rptr[r][c] = X->rptr[NewRowOrder[r-1]][c];
    }
  }
  return(XRO);
}


/*-----------------------------------------------------------------
  MatrixRandPermRows() - randomly reorders the rows of the input
  matrix.
  -----------------------------------------------------------------*/
int MatrixRandPermRows(MATRIX *X)
{
  int *NewRowOrder, r;
  MATRIX *X0;

  NewRowOrder = RandPerm(X->rows,NULL);
  for (r=0; r < X->rows; r++) NewRowOrder[r]++; // Make one-based
  X0 = MatrixCopy(X,NULL);
  MatrixReorderRows(X0, NewRowOrder, X);
  MatrixFree(&X0);
  return(0);
}


/*--------------------------------------------------------------------
  MatrixColsAreNotOrthog() - returns 1 if matrix columns are not
  orthogonal. Computes X'*X and examines the off diagonals. If any
  are greater than 2*FLT_MIN, then returns 1.
  --------------------------------------------------------------------*/
int MatrixColsAreNotOrthog(MATRIX *X)
{
  MATRIX *Xt, *XtX;
  int r,c;

  Xt = MatrixTranspose(X,NULL);
  XtX = MatrixMultiply(Xt,X,NULL);
  MatrixFree(&Xt);
  for (r=1; r <= XtX->rows; r++)
  {
    for (c=r+1; c <= XtX->cols; c++)
    {
      // only upper triangle
      if (XtX->rptr[r][c] > 2*FLT_MIN)
      {
        MatrixFree(&XtX);
        return(1);
      }
    }
  }
  MatrixFree(&XtX);
  return(0);
}


double
MatrixMahalanobisDistance(VECTOR *v_mean, MATRIX *m_inv_cov, VECTOR *v)
{
  VECTOR *v_dif, *v_tmp ;
  double dist ;

  v_dif = VectorSubtract(v_mean, v,NULL) ;
  v_tmp = MatrixMultiply(m_inv_cov, v_dif, NULL) ;
  dist = VectorDot(v_dif, v_tmp) ;

  VectorFree(&v_dif) ;
  VectorFree(&v_tmp) ;
  return(dist) ;
}

/*!
\fn double MatrixTransformDistance(MATRIX *m1, MATRIX *m2, double radius)
\brief RMS distance between two affine transforms, 
the center should be in the middle of image, and radius
should include the head (100 in RAS coord)
(see Jenkinson 1999 RMS deviation - tech report
    www.fmrib.ox.ac.uk/analysis/techrep )
\param m1     4x4 affine transformation
\param m2     4x4 affine transformation (may be NULL)
\param radius of the ball to be considered
*/
double MatrixTransformDistance(MATRIX *m1, MATRIX *m2, double radius)
{

   MATRIX* drigid = MatrixCopy(m1,NULL);
   if (m2) drigid = MatrixSubtract(drigid,m2,drigid);
   else //subtract identity
   {
      MATRIX *id = MatrixIdentity(4,NULL);
      drigid = MatrixSubtract(drigid,id,drigid);
      MatrixFree(&id);
   }

   // assert we have 4x4 affine transform
   //double EPS = 0.000001;
   //assert(drigid->rows ==4 && drigid->cols == 4);
   //assert(fabs(drigid->rptr[4][1]) < EPS);
   //assert(fabs(drigid->rptr[4][2]) < EPS);
   //assert(fabs(drigid->rptr[4][3]) < EPS);

   // translation norm quadrat:
   double tdq = 0;
   int i;
   for (i=1; i <= 3; i++)
   {
      tdq += drigid->rptr[i][4] * drigid->rptr[i][4];
      drigid->rptr[i][4] = 0.0; // set last row and column to zero
      drigid->rptr[4][i] = 0.0;
   }
   drigid->rptr[4][4] = 0.0;
   
   MATRIX* dt = MatrixTranspose(drigid, NULL);
   drigid = MatrixMultiply(dt,drigid,drigid);
   MatrixFree(&dt);
   
   // Trace of A^t A (only first 3x3 submatrix)
   double tr = 0.0;
   for (i=1; i <= 3; i++)
   {
      tr += drigid->rptr[i][i];
   }
   
   MatrixFree(&drigid);
   
   return sqrt((1.0/5.0) * radius*radius* tr + tdq);

}


/* for 3d vector macros */
#include "tritri.h"
int
MatrixOrthonormalizeTransform(MATRIX *m_L)
{
  double  dot, c1[3], c2[3], c3[3], len ;
  int     i ;

  for (i = 0 ; i < 3 ; i++)
  {
    c1[i] = *MATRIX_RELT(m_L, i+1, 1) ;
    c2[i] = *MATRIX_RELT(m_L, i+1, 2) ;
    c3[i] = *MATRIX_RELT(m_L, i+1, 3) ;
  }

  /* make 1st column vector unit length */
  len = VLEN(c1) ;
  if (FZERO(len))
    len = 1.0f ;
  SCALAR_MUL(c1, 1.0/len, c1) ;

  /* project out component of 2nd vector in direction of 1st column vector */
  dot = DOT(c1, c2) ;
  for (i = 0 ; i < 3 ; i++)
    c2[i] -= dot * c1[i] ;

  /* make 2nd column vector unit length */
  len = VLEN(c2) ;
  if (FZERO(len))
    len = 1.0f ;
  SCALAR_MUL(c2, 1.0/len, c2) ;

  /* project out component of 3rd vector in direction of 1st column vector */
  dot = DOT(c1, c3) ;
  for (i = 0 ; i < 3 ; i++)
    c3[i] -= dot * c1[i] ;

  /* project out component of 3rd vector in direction of 2nd column vector */
  dot = DOT(c2, c3) ;
  for (i = 0 ; i < 3 ; i++)
    c3[i] -= dot * c2[i] ;


  /* make 3rd column vector unit length */
  len = VLEN(c3) ;
  if (FZERO(len))
    len = 1.0f ;
  SCALAR_MUL(c3, 1.0/len, c3) ;

  for (i = 0 ; i < 3 ; i++)
  {
    *MATRIX_RELT(m_L, i+1, 1) = c1[i] ;
    *MATRIX_RELT(m_L, i+1, 2) = c2[i] ;
    *MATRIX_RELT(m_L, i+1, 3) = c3[i] ;
  }
#if 0
  /* remove translation component */
  for (i = 1 ; i <= 3 ; i++)
    *MATRIX_RELT(m_L, i, 4) = 0 ;
#endif

  return(NO_ERROR) ;
}

MATRIX *
MatrixAsciiReadRaw(const char *fname, MATRIX *m)
{
  FILE   *fp ;
  int    rows, cols, row, col ;
  char   line[10*STRLEN], *cp ;

  fp = fopen(fname, "r") ;
  if (fp == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "MatrixAsciiReadRaw: could not open %s\n", fname)) ;

  cp = fgetl(line, 10*STRLEN-1, fp) ;
  for (cols = rows = 0 ; cp != NULL ; rows++)
  {
    if (rows == 0)  // count cols in first row
    {
      cp = strtok(line, " ") ;
      for (cols = 0 ; cp != NULL ; cols++)
      {
        cp = strtok(NULL, " ") ;
      }
    }
    cp = fgetl(line, 10*STRLEN-1, fp) ;
  }
  m = MatrixAlloc(rows, cols, MATRIX_REAL) ;

  rewind(fp) ;
  cp = fgetl(line, 10*STRLEN-1, fp) ;
  for (row = 1 ; row <= rows ; row++)
  {
    cp = strtok(line, " ") ;
    for (col = 1 ; col <= cols ; col++)
    {
      sscanf(cp, "%f", MATRIX_RELT(m, row, col)) ;
      cp = strtok(NULL, " ") ;
    }
    cp = fgetl(line, 10*STRLEN-1, fp) ;
  }

  fclose(fp) ;
  return(m) ;
}

int
MatrixToRigidParameters(MATRIX *m, double *pxr, double *pyr, double *pzr,
                        double *pxt, double *pyt, double *pzt)
{
  // M = Mx * My * Mz
  *pxr = atan2(*MATRIX_RELT(m, 2, 3), *MATRIX_RELT(m, 3, 3)) ;
  *pyr = asin(-*MATRIX_RELT(m, 1, 3)) ;
  *pzr = atan2(*MATRIX_RELT(m, 1, 2), *MATRIX_RELT(m, 1, 1)) ;
  *pxt = *MATRIX_RELT(m, 1, 4) ;
  *pyt = *MATRIX_RELT(m, 2, 4) ;
  *pzt = *MATRIX_RELT(m, 3, 4) ;
  return(NO_ERROR) ;
}
MATRIX *
MatrixFromRigidParameters(MATRIX *m, double xr, double yr, double zr,
                          double xt, double yt, double zt)
{
  if (m == NULL)
    m = MatrixAlloc(4, 4, MATRIX_REAL) ;

  *MATRIX_RELT(m, 1, 1) = cos(yr) * cos(zr) ;
  *MATRIX_RELT(m, 1, 2) = cos(yr) * sin(zr) ;
  *MATRIX_RELT(m, 1, 3) = -sin(yr) ;

  *MATRIX_RELT(m, 2, 1) = sin(xr)*sin(yr)*cos(zr)-cos(xr)*sin(zr) ;
  *MATRIX_RELT(m, 2, 2) = sin(xr)*sin(yr)*sin(zr)+cos(xr)*cos(zr) ;
  *MATRIX_RELT(m, 2, 3) = sin(xr)*cos(yr) ;

  *MATRIX_RELT(m, 3, 1) = cos(xr)*sin(yr)*cos(zr)+sin(xr)*sin(zr) ;
  *MATRIX_RELT(m, 3, 2) = cos(xr)*sin(yr)*sin(zr) - sin(xr)*cos(zr) ;
  *MATRIX_RELT(m, 3, 3) = cos(xr)*cos(yr) ;

  *MATRIX_RELT(m, 1, 4) = xt ;
  *MATRIX_RELT(m, 2, 4) = yt ;
  *MATRIX_RELT(m, 3, 4) = zt ;
  *MATRIX_RELT(m, 4, 4) = 1.0 ;


  return(m) ;
}

int
MatrixCheckFinite(MATRIX *m)
{
  int r, c, retval = NO_ERROR ;

  for (r = 1 ; r < m->rows ; r++)
    for (c = 1 ; c < m->cols ; c++)
    {
      if (!finite(*MATRIX_RELT(m, r, c)))
      {
        DiagBreak() ;
        retval = ERROR_BADPARM ;
      }
    }
  return(retval) ;
}

/*!
  \fn MATRIX *MatrixKron(MATRIX *m1, MATRIX *m2, MATRIX *k)
  \brief Kronecker tensor product.
*/
MATRIX *MatrixKron(MATRIX *m1, MATRIX *m2, MATRIX *k)
{
  int rows, cols;
  int r1, c1, r2, c2, r, c;
  double v1, v2;
  
  rows = m1->rows * m2->rows;
  cols = m1->cols * m2->cols;

  if(!k){
    k = MatrixAlloc(rows,cols,MATRIX_REAL);
    if(!k){
      printf("ERROR: MatrixKron: could not alloc %d %d\n",rows,cols);
      return(NULL);
    }
  }
  else {
    if(k->rows != rows || k->cols != cols){
      printf("ERROR: MatrixKron: dimension mismatch %d %d vs %d %d\n",
	     rows,cols,k->rows,k->cols);
      return(NULL);
    }
  }

  for(r1 = 1; r1 <= m1->rows; r1++){
    for(c1 = 1; c1 <= m1->cols; c1++){
      v1 = m1->rptr[r1][c1];
      for(r2 = 1; r2 <= m2->rows; r2++){
	for(c2 = 1; c2 <= m2->cols; c2++){
	  v2 = m2->rptr[r2][c2];
	  r = (r1-1)*m2->rows + r2;
	  c = (c1-1)*m2->cols + c2;
	  //printf("%d %d   %d %d   %d %d\n",r1,c1,r2,c2,r,c);
	  k->rptr[r][c] = v1*v2;
	}
      }
    }
  }
  return(k);
}

/*!
  \fn double MatrixRowDotProduct(MATRIX *m, int row, VECTOR *v)
  \brief dot product of a vector with the row of a matrix
*/
double
MatrixRowDotProduct(MATRIX *m, int row, VECTOR *v)
{
  double dot ;
  int    col ;

  for (dot = 0.0, col = 1 ; col <= m->cols ; col++)
    dot += (*MATRIX_RELT(m, row, col) * VECTOR_ELT(v, col)) ;
  return(dot) ;
}

/*!
  \fn MATRIX *MatrixDemean(MATRIX *M, MATRIX *Mdm)
  \brief Removes the mean from each column
*/
MATRIX *MatrixDemean(MATRIX *M, MATRIX *Mdm)
{
  int r,c;
  double vsum,vmean;

  if(Mdm == NULL) {
    Mdm = MatrixAlloc(M->rows,M->cols,MATRIX_REAL);
    if(Mdm == NULL) return(NULL);
  }
  if(Mdm->rows != M->rows || Mdm->cols != M->cols){
    printf("MatrixDemean: dimension mismatch\n");
    return(NULL);    
  }
  for(c=1; c <= M->cols; c++){
    vsum = 0.0;
    for(r=1; r <= M->rows; r++)  vsum += M->rptr[r][c];
    vmean = vsum/M->rows;
    for(r=1; r <= M->rows; r++) M->rptr[r][c] -= vmean;
  }    
  return(Mdm);
}


/*!
  \fn MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude)
  \brief Creates a new matrix by excluding the given set of rows.
  \param Src - source matrix.
*/
MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude)
{
  MATRIX *Trg=NULL;
  int q, n, skip, m, c, nframesNew;

  nframesNew = Src->rows - nExclude;
  
  Trg = MatrixAlloc(nframesNew,Src->cols,MATRIX_REAL);
  q = 0;
  for(n=0; n < Src->rows; n++){
    skip = 0;
    for(m=0; m < nExclude; m++) if(n == ExcludeFrames[m]) skip = 1;
    if(skip) continue;
    for(c=0; c < Src->cols; c++) {
      //printf("%d %d %d %d\n",n,q,c,Src->cols);
      Trg->rptr[q+1][c+1] = Src->rptr[n+1][c+1];
    }
    q++;
  }
  return(Trg);
}
/*!
  \fn int MatrixIsIdentity(MATRIX *m)
  \brief returns 1 if matrix m is the identity matrix, 0 otherwise
*/
int
MatrixIsIdentity(MATRIX *m)
{
  int  r, c ;

  for (r = 1 ; r < m->rows ; r++)
    for (c = 1 ; c < m->cols ; c++)
    {
      if (r == c)
      {
	if (FEQUAL(*MATRIX_RELT(m, r, c), 1) == 0)
	  return(0) ;
      }
      else
	if (FEQUAL(*MATRIX_RELT(m, r, c), 0.0) == 0)
	  return(0) ;
    }
  return(1) ;
}


/*!
  \fn MATRIX *MatrixCumTrapZ(MATRIX *y, MATRIX *t, MATRIX *yz)
  \brief Computes trapezoidal integration (like matlab cumtrapz)
*/
MATRIX *MatrixCumTrapZ(MATRIX *y, MATRIX *t, MATRIX *yz)
{
  if(yz==NULL) yz = MatrixAlloc(y->rows,y->cols,MATRIX_REAL);

  int c, f;
  double v, vprev, vsum, dt;

  for(c=0; c < y->cols; c++){
    yz->rptr[1][c+1] = 0;
    vsum = 0;
    vprev = y->rptr[1][c+1];
    for(f=1; f < y->rows; f++){
      dt = t->rptr[f+1][1] - t->rptr[f][1];
      v = y->rptr[f+1][c+1];
      vsum += (dt*((v+vprev)/2));
      yz->rptr[f+1][c+1] = vsum;
      vprev = v;
    }
  }
  return(yz);
}
