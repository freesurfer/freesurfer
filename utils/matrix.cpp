/**
 * @brief Matrix utilities
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

#include <ctype.h>
#include <math.h>
#include <cmath>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#ifdef _POSIX_MAPPED_FILES
#include <sys/mman.h>
#endif
#ifdef Darwin
#include <float.h>  // defines FLT_MIN
#else
#ifdef SunOS
#include <limits.h>  // defines FLT_MIN
#else
#ifdef Windows_NT
#include <float.h>  // defines FLT_MIN
#else
#include <values.h>  // defines FLT_MIN
#endif
#endif
#endif

#include "diag.h"
#include "error.h"
#include "evschutils.h"
#include "fio.h"
#include "macros.h"
#include "matrix.h"
#include "numerics.h"
#include "proto.h"
#include "utils.h"

#include "romp_support.h"

// private functions
MATRIX *MatrixCalculateEigenSystemHelper(MATRIX *m, float *evalues, MATRIX *m_evectors, int isSymmetric);




typedef struct MatrixStats {
    const char* file;
    int         line;
    int     	rowsCols;
    int         count;
} MatrixStats;

#define matrixStatsSizeLog2 13
#define matrixStatsSize     (1<<matrixStatsSizeLog2)
#define matrixStatsSizeMask (matrixStatsSize-1)

static MatrixStats* matrixStats = NULL;

static int stats_compare(const void* lhs_ptr, const void* rhs_ptr) {
   int lhs = *(int*)lhs_ptr;
   int rhs = *(int*)rhs_ptr;
   size_t lhsPriority = matrixStats[lhs].count;
   size_t rhsPriority = matrixStats[rhs].count;
   if (lhsPriority < rhsPriority) return +1;    // ascending order
   if (lhsPriority > rhsPriority) return -1;    // ascending order
   return 0;
}

static void matrixStatsExitHandler(void) {
    if (!matrixStats) return;

    size_t count = 0;
    size_t i;
    for (i = 0; i < matrixStatsSize; i++) {
        MatrixStats* m = &matrixStats[i];
        if (!m->line) continue;
        count++;
    }

    int* indexs = (int*)malloc(count*sizeof(int));
    count = 0;
    for (i = 0; i < matrixStatsSize; i++) {
        MatrixStats* m = &matrixStats[i];
        if (!m->line) continue;
        indexs[count++] = i;
    }

    qsort(indexs, count, sizeof(int), stats_compare);
       
    fprintf(stdout, "MatrixStats\n   file, line, rowsCols, count\n");
    for (i = 0; i < count; i++) {
        MatrixStats* m = &matrixStats[indexs[i]];
        fprintf(stdout, "%s, %d, %d, %d\n",
            m->file, m->line, m->rowsCols, m->count);
    }
    fprintf(stdout, "MatrixStats\n   file, line, rowsCols, count\n");
}

static void noteMatrixAlloced(const char* file, int line, int rows, int cols) {

    if (1) return;
    
    if (!matrixStats) {
        matrixStats = (MatrixStats*)calloc(matrixStatsSize, sizeof(MatrixStats));
        atexit(matrixStatsExitHandler);
    }
    
    int rowsCols = rows*1000000+cols;
    size_t stabs = 0;
    size_t hash = (((size_t)line ^ (size_t)file)*75321) & matrixStatsSizeMask;
    MatrixStats* m = &matrixStats[hash];
    while (m->line != line || m->file != file || m->rowsCols != rowsCols ) {
        if (m->line == 0) break; // not in chain
        hash = (hash*327 + 1) & matrixStatsSizeMask;
        m = &matrixStats[hash];
        if (++stabs > 1000) *(int*)-1 = 0;  // table too full 
    }
    
    if (m->line == 0) {     // There is a slight chance this might find the same empty cell as another thread - who cares?
        m->line = line; m->file = file; m->rowsCols = rowsCols;
    }
    
    m->count++;
}


/**
 * Returns true if the matrix is symmetric (should be square too).
 */
int MatrixIsSymmetric(MATRIX *matrix);

MATRIX *MatrixCopy(const MATRIX *mIn, MATRIX *mOut)
{
  int row, rows, cols, col;

  if (mIn == NULL) {
    if (mOut) {
      MatrixFree(&mOut);
      mOut = NULL;
    }
    return (NULL);
  }
  rows = mIn->rows;
  cols = mIn->cols;

  if (!mOut) mOut = MatrixAlloc(rows, cols, mIn->type);

  if (!mOut) ErrorExit(ERROR_NO_MEMORY, "MatrixCopy: couldn't allocate mOut");

#if 1
  for (col = 1; col <= cols; col++)
    for (row = 1; row <= rows; row++) *MATRIX_RELT(mOut, row, col) = *MATRIX_RELT(mIn, row, col);
#else
  for (row = 1; row <= rows; row++)
    memmove((char *)(mOut->rptr[row]), (char *)mIn->rptr[row], (cols + 1) * sizeof(float));
#endif

  return (mOut);
}

int MatrixQRdecomposition(const MATRIX *iMatrix, MATRIX *oQ, MATRIX *oR)
{
  int r = OpenQRdecomposition(iMatrix,oQ,oR);
  return(r);
}

MATRIX *MatrixInverse(const MATRIX *mIn, MATRIX *mOut)
{
  // float **a, **y;
  int isError, i, j, rows, cols, alloced = 0;
  MATRIX *mTmp;

  if (!mIn) {
    ErrorExit(ERROR_BADPARM, "MatrixInverse: NULL input matrix!\n");
  }
  
  if (mIn->rows != mIn->cols)
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixInverse: matrix (%d x %d) is not square\n", mIn->rows, mIn->cols));

  rows = mIn->rows;
  cols = mIn->cols;

  if (!mOut) {
    alloced = 1;
    mOut = MatrixAlloc(rows, cols, mIn->type);
  }

  /* allocate temp matrix so as not to destory contents of mIn */
  if (mIn->type == MATRIX_COMPLEX) {
    MATRIX *mQuad, *mInv, *mReal, *mImag;

    mTmp = MatrixAlloc(2 * rows, 2 * cols, MATRIX_REAL);

    /*
      form square matrix of the form

      A  -C
      C  A

      where A and C are the real and imaginary components of the input
      matrix respectively.
    */
    MatrixCopyRealRegion(mIn, mTmp, 1, 1, rows, cols, 1, 1);
    MatrixCopyRealRegion(mIn, mTmp, 1, 1, rows, cols, rows + 1, cols + 1);
    mQuad = MatrixAlloc(rows, cols, MATRIX_REAL);
    MatrixCopyImagRegion(mIn, mQuad, 1, 1, rows, cols, 1, 1);
    MatrixScalarMul(mQuad, -1.0f, mQuad);
    MatrixCopyRegion(mQuad, mTmp, 1, 1, rows, cols, 1, cols + 1);
    MatrixCopyImagRegion(mIn, mTmp, 1, 1, rows, cols, cols + 1, 1);

#if 0
    DiagPrintf(DIAG_MATRIX, "\nextraction\n") ;
    MatrixPrint(stdout, mTmp) ;
#endif

    mReal = MatrixAlloc(rows, cols, MATRIX_REAL);
    mImag = MatrixAlloc(rows, cols, MATRIX_REAL);
    mInv = MatrixInverse(mTmp, NULL);

    MatrixCopyRegion(mInv, mReal, 1, 1, rows, cols, 1, 1);
    MatrixCopyRegion(mInv, mImag, rows + 1, 1, rows, cols, 1, 1);
    MatrixRealToComplex(mReal, mImag, mOut);

#if 0
    DiagPrintf(DIAG_MATRIX, "\ninverse of extraction\n") ;
    MatrixPrint(stderr, mInv) ;
    DiagPrintf(DIAG_MATRIX, "\ninverse:\n") ;
    MatrixPrint(stderr, mOut) ;
    DiagPrintf(DIAG_MATRIX, "real + imag = \n") ;
    MatrixPrint(stderr, mReal) ;
    MatrixPrint(stderr, mImag) ;
#endif
    MatrixFree(&mQuad);
    MatrixFree(&mReal);
    MatrixFree(&mImag);
  }
  else {
    mTmp = MatrixCopy(mIn, NULL);

    // a = mTmp->rptr;
    // y = mOut->rptr;

    isError = OpenLUMatrixInverse(mTmp, mOut);

    if (isError < 0) {
      MatrixFree(&mTmp);
      if (alloced) {
        MatrixFree(&mOut);
      }
      return (NULL);
    }
  }

  MatrixFree(&mTmp);

  for (j = 1; j <= rows; j++) {
    for (i = 1; i <= rows; i++) switch (mOut->type) {
        case MATRIX_REAL:
          if (!std::isfinite(*MATRIX_RELT(mOut, i, j))) {
            if (alloced) MatrixFree(&mOut);
            return (NULL); /* was singular */
          }
          break;
        case MATRIX_COMPLEX:
          if (!std::isfinite(MATRIX_CELT_REAL(mOut, i, j)) || !std::isfinite(MATRIX_CELT_IMAG(mOut, i, j))) {
            if (alloced) MatrixFree(&mOut);
            return (NULL); /* was singular */
          }
          break;
      }
  }

  return (mOut);
}

static MATRIX *MatrixAlloc_old(const int rows, const int cols, const int type)
{
  MATRIX *mat;
  int row, nelts;
#ifdef _POSIX_MAPPED_FILES
  int i;
  float f;
#endif

  mat = (MATRIX *)calloc(1, sizeof(MATRIX));
  if (!mat) ErrorExit(ERROR_NO_MEMORY, "MatrixAlloc(%d, %d, %d): could not allocate mat", rows, cols, type);

  mat->rows  = rows;
  mat->cols  = cols;
  mat->inBuf = false;
  mat->type  = type;

  /*
    allocate a single array the size of the matrix, then initialize
    row pointers to the proper locations.
  */

  nelts = rows * cols;
  if (type == MATRIX_COMPLEX) nelts *= 2;

  /*
    because NRC is one-based, we must leave room for a few unused
    (but valid) addresses before the start of the actual data so
    that mat->rptr[0][0] is a valid address even though it wont
    be used.
  */
  mat->data = (float *)calloc(nelts + 2, sizeof(float));

  mat->mmapfile = NULL;

#ifdef _POSIX_MAPPED_FILES
  if (!mat->data) /* First try to allocate a mmap'd tmpfile */
  {
    printf("MatrixAlloc(%d, %d): Using mmap'd tmpfile\n", rows, cols);
    if ((mat->mmapfile = tmpfile())) {
      /* This maintains identical behavior with calloc */
      f = 0;
      for (i = 0; i < nelts + 2; ++i) {
        if (!(fwrite(&f, sizeof(float), 1, mat->mmapfile))) {
          printf("MatrixAlloc(%d, %d): fwrite failed", rows, cols);
          exit(1);
        }
      }
      /* This seems to matter with some implementations of mmap */
      fseek(mat->mmapfile, 0, 0);
      /* lseek(fileno(mat->mapfile), (nelts+2) * sizeof(float), 0) ;*/
      fflush(mat->mmapfile);

      mat->data =
          (float *)mmap(0, (nelts + 2) * sizeof(float), PROT_READ | PROT_WRITE, MAP_SHARED, fileno(mat->mmapfile), 0);

      if (mat->data == MAP_FAILED) {
        mat->data = 0;
      }
    }
  }
#endif

  if (!mat->data) /* we _still_ couldn't get it! */
  {
    fprintf(stderr, "MatrixAlloc(%d, %d): allocation failed\n", rows, cols);
    exit(1);
  }
  mat->data += 2;

  /*
     silly numerical recipes in C requires 1-based stuff. The full
     data array is zero based, point the first row to the zeroth
     element, and so on.
  */
  mat->rptr = (float **)calloc(rows + 1, sizeof(float *));
  if (!mat->rptr) {
    free(mat->data);
    free(mat);
    ErrorExit(ERROR_NO_MEMORY, "MatrixAlloc(%d, %d): could not allocate rptr", rows, cols);
  }
  for (row = 1; row <= rows; row++) {
    switch (type) {
      case MATRIX_REAL:
        mat->rptr[row] = mat->data + (row - 1) * cols - 1;
        break;
      case MATRIX_COMPLEX:
        mat->rptr[row] = (float *)(((CPTR)mat->data) + (row - 1) * cols - 1);
        break;
      default:
        ErrorReturn(NULL, (ERROR_BADPARM, "MatrixAlloc: unknown type %d\n", type));
    }
  }
  return (mat);
}

static MATRIX *MatrixAlloc_new(
    const int rows, 
    const int cols, 
    const int type,
    MatrixBuffer* const buf)
{

  // Calculate the storage size, to get in one allocation
  // (1) to reduce the calls to malloc by 3x, since it is an important part of mris_fix_topology
  // (2) to prepare to have a stack buffer that the MATRIX can be in, to reduce the calls to malloc for 2x2 3x3 matrices
  //
  size_t size_needed = sizeof(MATRIX);
  size_needed = (size_needed + 63) & ~63;   // round up
  
  // Decide on the number of ptrs
  //
  size_t const rptr_offset = size_needed;
  size_needed += (rows + 1) * sizeof(float *);
  size_needed = (size_needed + 63) & ~63;   // round up

  // Decide on the number of data elements
  // NRC is one-based, so leave room for a few unused (but valid) addresses before the start of the actual data so
  // that mat->rptr[0][0] is a valid address even though it wont be used.
  //
  size_t const data_offset = size_needed;
  int const nelts = ((rows * cols) + 2) * ((type == MATRIX_COMPLEX) ? 2 : 1);
  size_needed += nelts*sizeof(float);
  size_needed = (size_needed + 63) & ~63;   // round up

  // Try to get the matrix and the rptrs with out without the data
  //
  MATRIX* mat  = NULL;
  float*  data = NULL; 
  {
    static long count, limit = 128, bufsSupplied, bufsUsed;
    count++;
    if (buf) bufsSupplied++;
    
    void* memptr;
    if (buf && size_needed <= sizeof(*buf)) {
      bufsUsed++;
      memptr = &buf->matrix;
      mat    = &buf->matrix;
      data   = (float*) ((char*)memptr + data_offset);
    } else if (!posix_memalign(&memptr, 64, size_needed)) {
      mat    = (MATRIX*)memptr;
      data   = (float*) ((char*)memptr + data_offset);
    } else if (!posix_memalign(&memptr, 64, data_offset)) {
      mat    = (MATRIX*)memptr;
    }
    
    if (0 && (count >= limit)) {
      limit *= 2;
      fprintf(stdout, "%s:%d MatrixAlloc_new stats  count:%g bufsSupplied:%g bufsUsed:%g\n", __FILE__, __LINE__, 
      	(float)count, (float)bufsSupplied, (float)bufsUsed);
      fprintf(stdout, "      MatrixAlloc_new size_needed:%ld sizeof(buf):%ld rows:%d cols:%d\n",
      	size_needed, sizeof(*buf), rows, cols);
    }
  }
  
  FILE* mmapfile = NULL;
  
#ifdef _POSIX_MAPPED_FILES
  if (!data) { // Try to allocate a mmap'd tmpfile

    printf("MatrixAlloc(%d, %d): Using mmap'd tmpfile\n", rows, cols);
    if ((mmapfile = tmpfile())) {

      // This seems to matter with some implementations of mmap
      //
      fseek(mmapfile, 0, 0);
      fflush(mmapfile);

      data = (float *)mmap(0, (nelts + 2) * sizeof(float), PROT_READ | PROT_WRITE, MAP_SHARED, fileno(mat->mmapfile), 0);

      if (data == MAP_FAILED) {
        data = NULL;
      }
    }
  }
#endif

  if (!mat || !data) {
    fprintf(stderr, "MatrixAlloc(%d, %d): allocation failed\n", rows, cols);
    exit(1);
  }

  bzero(mat,  data_offset);
  bzero(data, nelts*sizeof(float));
  
  mat->rows  = rows;
  mat->cols  = cols;
  mat->inBuf = (mat == &buf->matrix);
  mat->type  = type;

  mat->data  = data;
  mat->data  += 2;

  mat->mmapfile = mmapfile;

  // silly numerical recipes in C requires 1-based stuff. The full
  // data array is zero based, point the first row to the zeroth
  // element, and so on.
  //
  float** rptr = (float**)((char*)mat + rptr_offset);
  mat->rptr = rptr;

  int row;
  for (row = 1; row <= rows; row++) {
    switch (type) {
      case MATRIX_REAL:
        rptr[row] = mat->data + (row - 1) * cols - 1;
        break;
      case MATRIX_COMPLEX:
        rptr[row] = (float *)(((CPTR)mat->data) + (row - 1) * cols - 1);
        break;
      default:
        ErrorReturn(NULL, (ERROR_BADPARM, "MatrixAlloc: unknown type %d\n", type));
    }
  }

  return (mat);
}


static int use_new_MatricAlloc() {
    static int once, result;
    if (!once) {
        once++;
        result = !getenv("FREESURFER_MatrixAlloc_old");
    }
    return result;
}

MATRIX *MatrixAlloc_wkr( const int rows, const int cols, const int type,
    	    	     const char* callSiteFile, int callSiteLine) {
		     
    noteMatrixAlloced(callSiteFile, callSiteLine, rows, cols);
    
    if (use_new_MatricAlloc()) return MatrixAlloc_new(rows, cols, type, NULL);
    else         	       return MatrixAlloc_old(rows, cols, type);
}


MATRIX *MatrixAlloc2_wkr(const int rows, const int cols, const int type, MatrixBuffer* buf,
    	    	     const char* callSiteFile, int callSiteLine) {

    noteMatrixAlloced(callSiteFile, callSiteLine, rows, cols);

    if (use_new_MatricAlloc()) return MatrixAlloc_new(rows, cols, type, buf);
    else    	    	       return MatrixAlloc_old(rows, cols, type);
}


int MatrixFree(MATRIX **pmat)
{
  MATRIX *mat;

  if (!pmat) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixFree: NULL pmat POINTER!\n"));

  mat = *pmat;
  *pmat = NULL;

  if (!mat || mat->inBuf) return (0);

  /* silly numerical recipes in C requires 1-based stuff */
  mat->data -= 2;
  if (mat->mmapfile) {
    int nelts;

    nelts = mat->rows * mat->cols;
    if (mat->type == MATRIX_COMPLEX) nelts *= 2;

#ifdef _POSIX_MAPPED_FILES
    munmap((void *)mat->data, (nelts + 2) * sizeof(float));
#endif
    fclose(mat->mmapfile);
  }
  else {
    if (!use_new_MatricAlloc()) free(mat->data);
  }

  if (!use_new_MatricAlloc()) free(mat->rptr);
  if (!mat->inBuf) free(mat);

  return (0);
}

/*!
  \fn MATRIX *MatrixMultiplyD( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
  \brief Multiplies two matrices. The accumulation is done with double,
   which is more accurate than MatrixMultiply() which uses float.
*/
MATRIX *MatrixMultiplyD(const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
{
  int col, row, i, rows, cols, m1_cols;
  float *r3;
  float *r1, *r2;
  double val;
  MATRIX *m_tmp1 = NULL, *m_tmp2 = NULL;
  char tmpstr[1000];

  if (!m1) {
    sprintf(tmpstr, "MatrixMultiplyD(): m1 is null\n break %s:%d\n", __FILE__, __LINE__);
    ErrorExit(ERROR_BADPARM, tmpstr);
  }
  if (!m2) {
    sprintf(tmpstr, "MatrixMultiplyD(): m2 is null\n break %s:%d\n", __FILE__, __LINE__);
    ErrorExit(ERROR_BADPARM, tmpstr);
  }
  if (m1->cols != m2->rows) {
    sprintf(tmpstr,
            "MatrixMultiplyD(): m1 cols %d does not match m2 rows %d\n break %s:%d\n",
            m1->cols,
            m2->rows,
            __FILE__,
            __LINE__);
    ErrorReturn(NULL, (ERROR_BADPARM, "%s", tmpstr));
  }

  if (!m3) {
    /* twitzel also did something here */
    if ((m1->type == MATRIX_COMPLEX) || (m2->type == MATRIX_COMPLEX))
      m3 = MatrixAlloc(m1->rows, m2->cols, MATRIX_COMPLEX);
    else
      m3 = MatrixAlloc(m1->rows, m2->cols, m1->type);
    if (!m3) return (NULL);
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols)) {
    printf("MatrixMultiplyD(): m1/m2 dim mismatch\n break %s:%d\n", __FILE__, __LINE__);
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiplyD: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows,
                 m1->cols,
                 m2->rows,
                 m2->cols,
                 m3->rows,
                 m3->cols));
  }

  if (m3 == m2) {
    m_tmp1 = MatrixCopy(m2, NULL);
    m2 = m_tmp1;
  }
  if (m3 == m1) {
    m_tmp2 = MatrixCopy(m1, NULL);
    m1 = m_tmp2;
  }
  /*  MatrixClear(m3) ;*/
  cols = m3->cols;
  rows = m3->rows;
  m1_cols = m1->cols;

  /* twitzel modified here */
  if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_REAL)) {
    for (row = 1; row <= rows; row++) {
      r3 = &m3->rptr[row][1];
      for (col = 1; col <= cols; col++) {
        val = 0.0;
        r1 = &m1->rptr[row][1];
        r2 = &m2->rptr[1][col];
        for (i = 1; i <= m1_cols; i++, r2 += cols) val += (double)(*r1++) * (*r2);
        *r3++ = val;
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_COMPLEX)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          double a, b, c, d; /* a + ib and c + id */

          a = MATRIX_CELT_REAL(m1, row, i);
          b = MATRIX_CELT_IMAG(m1, row, i);
          c = MATRIX_CELT_REAL(m2, i, col);
          d = MATRIX_CELT_IMAG(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c - b * d;
          MATRIX_CELT_IMAG(m3, row, col) += a * d + b * c;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_COMPLEX)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          double a, c, d; /* a + ib and c + id and b=0 here*/

          a = *MATRIX_RELT(m1, row, i);
          c = MATRIX_CELT_REAL(m2, i, col);
          d = MATRIX_CELT_IMAG(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c;
          MATRIX_CELT_IMAG(m3, row, col) += a * d;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_REAL)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          double a, b, c; /* a + ib and c + id and d=0 here*/

          a = MATRIX_CELT_REAL(m1, row, i);
          b = MATRIX_CELT_IMAG(m1, row, i);
          c = *MATRIX_RELT(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c;
          MATRIX_CELT_IMAG(m3, row, col) += b * c;
        }
      }
    }
  }
  if (m_tmp1) MatrixFree(&m_tmp1);
  if (m_tmp2) MatrixFree(&m_tmp2);
  return (m3);
}

/*!
  \fn MATRIX *MatrixMultiply( const MATRIX *m1, const MATRIX *m2, MATRIX *m3)
  \brief Multiplies two matrices. The accumulation is done with float.
   Consider using MatrixMultiplyD() which uses double.
*/
MATRIX *MatrixMultiply_wkr(const MATRIX *m1, const MATRIX *m2, MATRIX *m3, const char* callSiteFile, int callSiteLine)
{
  int col, row, i, rows, cols, m1_cols;
  float *r3;
  float val, *r1, *r2;
  MATRIX *m_tmp1 = NULL, *m_tmp2 = NULL;

  if (!m1) ErrorExit(ERROR_BADPARM, "MatrixMultiply: m1 is null!\n");
  if (!m2) ErrorExit(ERROR_BADPARM, "MatrixMultiply: m2 is null!\n");

  if (m1->cols != m2->rows) {
    printf("MatrixMultiply(): m1/m2 dim mismatch\n break %s:%d\n", __FILE__, __LINE__);
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixMultiply: m1 cols %d does not match m2 rows %d\n", m1->cols, m2->rows));
  }

  if (!m3) {
    noteMatrixAlloced(callSiteFile, callSiteLine, -m1->rows, -m2->cols);
  
    /* twitzel also did something here */
    if ((m1->type == MATRIX_COMPLEX) || (m2->type == MATRIX_COMPLEX))
      m3 = MatrixAlloc(m1->rows, m2->cols, MATRIX_COMPLEX);
    else
      m3 = MatrixAlloc(m1->rows, m2->cols, m1->type);
    if (!m3) return (NULL);
  }
  else if ((m3->rows != m1->rows) || (m3->cols != m2->cols)) {
    printf("MatrixMultiply(): m3 dim mismatch\n break %s:%d\n", __FILE__, __LINE__);
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixMultiply: (%d x %d) * (%d x %d) != (%d x %d)\n",
                 m1->rows,
                 m1->cols,
                 m2->rows,
                 m2->cols,
                 m3->rows,
                 m3->cols));
  }

  if (m3 == m2) {
    m_tmp1 = MatrixCopy(m2, NULL);
    m2 = m_tmp1;
  }
  if (m3 == m1) {
    m_tmp2 = MatrixCopy(m1, NULL);
    m1 = m_tmp2;
  }
  /*  MatrixClear(m3) ;*/
  cols = m3->cols;
  rows = m3->rows;
  m1_cols = m1->cols;

  /* twitzel modified here */
  if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_REAL)) {
    for (row = 1; row <= rows; row++) {
      r3 = &m3->rptr[row][1];
      for (col = 1; col <= cols; col++) {
        val = 0.0;
        r1 = &m1->rptr[row][1];
        r2 = &m2->rptr[1][col];
        for (i = 1; i <= m1_cols; i++, r2 += cols) {
#if 0
          m3->rptr[row][col] +=
            m1->rptr[row][i] * m2->rptr[i][col] ;
#else
          val += *r1++ * *r2;
#endif
        }
        *r3++ = val;
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_COMPLEX)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          float a, b, c, d; /* a + ib and c + id */

          a = MATRIX_CELT_REAL(m1, row, i);
          b = MATRIX_CELT_IMAG(m1, row, i);
          c = MATRIX_CELT_REAL(m2, i, col);
          d = MATRIX_CELT_IMAG(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c - b * d;
          MATRIX_CELT_IMAG(m3, row, col) += a * d + b * c;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_REAL) && (m2->type == MATRIX_COMPLEX)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          float a, c, d; /* a + ib and c + id and b=0 here*/

          a = *MATRIX_RELT(m1, row, i);
          c = MATRIX_CELT_REAL(m2, i, col);
          d = MATRIX_CELT_IMAG(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c;
          MATRIX_CELT_IMAG(m3, row, col) += a * d;
        }
      }
    }
  }
  else if ((m1->type == MATRIX_COMPLEX) && (m2->type == MATRIX_REAL)) {
    for (row = 1; row <= rows; row++) {
      for (col = 1; col <= cols; col++) {
        for (i = 1; i <= m1->cols; i++) {
          float a, b, c; /* a + ib and c + id and d=0 here*/

          a = MATRIX_CELT_REAL(m1, row, i);
          b = MATRIX_CELT_IMAG(m1, row, i);
          c = *MATRIX_RELT(m2, i, col);
          MATRIX_CELT_REAL(m3, row, col) += a * c;
          MATRIX_CELT_IMAG(m3, row, col) += b * c;
        }
      }
    }
  }
  if (m_tmp1) MatrixFree(&m_tmp1);
  if (m_tmp2) MatrixFree(&m_tmp2);
  return (m3);
}

int MatrixPrint(FILE *fp, const MATRIX *mat)
{
  int row, col, rows, cols;

  if (fp == NULL) {
    fp = stdout;
    ErrorPrintf(ERROR_BADPARM, "MatrixPrint: fp = NULL!");
  }
  if (mat == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: mat = NULL!"));

  rows = mat->rows;
  cols = mat->cols;

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, "% 2.5f", mat->rptr[row][col]);
          break;
        case MATRIX_COMPLEX:
          fprintf(fp, "% 2.3f + % 2.3f i", MATRIX_CELT_REAL(mat, row, col), MATRIX_CELT_IMAG(mat, row, col));
          break;
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: unknown type %d\n", mat->type));
      }
#if 0
      if (col < cols)
        fprintf(fp, " | ") ;
#else
      if (col < cols) fprintf(fp, "  ");
#endif
    }
    fprintf(fp, ";\n");
  }
  fflush(stdout);
  return (NO_ERROR);
}

int MatrixPrintWithString(FILE *fp, MATRIX *m, const char *Pre, const char *Post)
{
  int err;
  fprintf(fp, "%s", Pre);
  err = MatrixPrint(fp, m);
  fprintf(fp, "%s", Post);
  fflush(fp);
  return (err);
}

int MatrixPrintFmt(FILE *fp, const char *fmt, MATRIX *mat)
{
  int row, col, rows, cols;

  if (fp == NULL) {
    fp = stdout;
    ErrorPrintf(ERROR_BADPARM, "MatrixPrint: fp = NULL!");
  }
  if (mat == NULL) ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: mat = NULL!"));

  rows = mat->rows;
  cols = mat->cols;

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, fmt, mat->rptr[row][col]);
          break;
        case MATRIX_COMPLEX:
          fprintf(fp, fmt, MATRIX_CELT_REAL(mat, row, col));
          fprintf(fp, " + ");
          fprintf(fp, fmt, MATRIX_CELT_IMAG(mat, row, col));
          fprintf(fp, " i");
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

int MatrixPrintOneLine(FILE *fp, MATRIX *mat)
{
  int row, col, rows, cols;

  rows = mat->rows;
  cols = mat->cols;

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, "% 2.3f", mat->rptr[row][col]);
          break;
        case MATRIX_COMPLEX:
          fprintf(fp, "% 2.3f + % 2.3f i", MATRIX_CELT_REAL(mat, row, col), MATRIX_CELT_IMAG(mat, row, col));
          break;
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: unknown type %d\n", mat->type));
      }
#if 0
      if (col < cols)
        fprintf(fp, " | ") ;
#else
      if (col < cols) fprintf(fp, " ");
#endif
    }
    fprintf(fp, ";   ");
  }
  return (NO_ERROR);
}

int MatrixPrintTranspose(FILE *fp, MATRIX *mat)
{
  int row, col, rows, cols;

  rows = mat->rows;
  cols = mat->cols;

  for (col = 1; col <= cols; col++) {
    for (row = 1; row <= rows; row++) {
      switch (mat->type) {
        case MATRIX_REAL:
          fprintf(fp, "% 2.3f", mat->rptr[row][col]);
          break;
        case MATRIX_COMPLEX:
          fprintf(fp, "% 2.3f + % 2.3f i", MATRIX_CELT_REAL(mat, row, col), MATRIX_CELT_IMAG(mat, row, col));
          break;
        default:
          ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "MatrixPrint: unknown type %d\n", mat->type));
      }
#if 0
      if (row < rows)
        fprintf(fp, " | ") ;
#else
      if (row < rows) fprintf(fp, "  ");
#endif
    }
    fprintf(fp, "\n");
  }
  fflush(stdout);
  return (NO_ERROR);
}

MATRIX *MatrixReadTxt(const char *fname, MATRIX *mat)
{
  FILE *fp;
  int rows, cols, row, col, nlinemax, nread;
  char line[1000];
  // char    *cp;
  nlinemax = 999;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "MatrixReadTxt(%s) - file open failed\n", fname));

  // Read in the first line, including the newline
  if (!fgets(line, nlinemax, fp))
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixReadTxt: could not read 1st line from %s\n", fname));

  // for(cols = 0,cp = strtok(line, " \t,") ; cp ; cp = strtok(NULL, " \t,"))
  //  cols++ ;

  cols = ItemsInString(line);
  if (cols < 1) {
    // Return quitely in case try to use another format.
    return (NULL);
  }

  // Count the number of lines, start at row=1 because
  // a line was read above.
  for (rows = 1; fgets(line, nlinemax, fp); rows++) {
  }

  // printf("MatrixReadTxt: %d rows and %d cols\n",rows,cols);
  mat = MatrixAlloc(rows, cols, MATRIX_REAL);

  rewind(fp);
  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      nread = fscanf(fp, "%f", &mat->rptr[row][col]);
      if (nread != 1) {
        MatrixFree(&mat);
        ErrorReturn(NULL, (ERROR_BADPARM, "MatrixReadTxT: could not scan value [%d][%d]\n", row, col));
      }
    }
  }
  fclose(fp);
  return (mat);

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
          "MatrixReadTxt: could not scan value [%d][%d]\n", row, col));
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
MATRIX *MatrixRead(const char *fname) { return (MatlabRead(fname)); }

int MatrixWrite(MATRIX *mat, const char *fname, const char *name)
{
  if (!name) return (MatlabWrite(mat, fname, fname)); /* name of matrix in .mat file */
  return (MatlabWrite(mat, fname, name));
}

MATRIX *MatrixIdentity(int n, MATRIX *mat)
{
  int i;

  if (!mat)
    mat = MatrixAlloc(n, n, MATRIX_REAL);
  else
    MatrixClear(mat);

  for (i = 1; i <= n; i++) mat->rptr[i][i] = 1;

  return (mat);
}

MATRIX *MatrixTranspose(MATRIX *mIn, MATRIX *mOut)
{
  int row, col, rows, cols;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->cols, mIn->rows, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  int insitu = (mIn == mOut);
  if (insitu) mIn = MatrixCopy(mOut, NULL);

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) mOut->rptr[col][row] = mIn->rptr[row][col];
  }
  if (insitu) MatrixFree(&mIn);
  return (mOut);
}

MATRIX *MatrixAdd(const MATRIX *m1, const MATRIX *m2, MATRIX *mOut)
{
  int row, col, rows, cols;

  rows = m1->rows;
  cols = m1->cols;

  if ((rows != m2->rows) || (cols != m2->cols))
    ErrorReturn(
        NULL, (ERROR_BADPARM, "MatrixAdd: incompatable matrices %d x %d + %d x %d\n", rows, cols, m2->rows, m2->cols));

  if (!mOut) {
    mOut = MatrixAlloc(rows, cols, m1->type);
    if (!mOut) return (NULL);
  }

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(
        NULL,
        (ERROR_BADPARM, "MatrixAdd: incompatable matrices %d x %d = %d x %d\n", rows, cols, mOut->rows, mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      mOut->rptr[row][col] = m1->rptr[row][col] + m2->rptr[row][col];
    }
  }

  return (mOut);
}

MATRIX *MatrixSubtract(const MATRIX *m1, const MATRIX *m2, MATRIX *mOut)
{
  int row, col, rows, cols;

  rows = m1->rows;
  cols = m1->cols;

  if ((rows != m2->rows) || (cols != m2->cols))
    ErrorReturn(
        NULL,
        (ERROR_BADPARM, "MatrixSubtract: incompatable matrices %d x %d - %d x %d\n", rows, cols, m2->rows, m2->cols));

  if (!mOut) {
    mOut = MatrixAlloc(rows, cols, m1->type);
    if (!mOut) return (NULL);
  }

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSubtract: incompatable matrices %d x %d = %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) mOut->rptr[row][col] = m1->rptr[row][col] - m2->rptr[row][col];
  }

  return (mOut);
}

MATRIX *MatrixScalarMul(const MATRIX *mIn, const float val, MATRIX *mOut)
{
  int row, col, rows, cols;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixScalarMul: incompatable matrices %d x %d != %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) mOut->rptr[row][col] = mIn->rptr[row][col] * val;
  }
  return (mOut);
}
MATRIX *VectorZeroMean(const MATRIX *mIn, MATRIX *mOut)
{
  double mean;

  mean = VectorMean(mIn);
  return (MatrixScalarAdd(mIn, -mean, mOut));
}

MATRIX *MatrixScalarAdd(const MATRIX *mIn, const float val, MATRIX *mOut)
{
  int row, col, rows, cols;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixScalarAdd: incompatable matrices %d x %d != %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) mOut->rptr[row][col] = mIn->rptr[row][col] + val;
  }
  return (mOut);
}

MATRIX *MatrixClear(MATRIX *mat)
{
  int row, cols;
  // int rows;
  // rows = mat->rows;
  cols = mat->cols;
  for (row = 1; row <= mat->rows; row++) memset((char *)mat->rptr[row], 0, (cols + 1) * sizeof(float));

  return (mat);
}

MATRIX *MatrixSquareElts(MATRIX *mIn, MATRIX *mOut)
{
  int row, col, rows, cols;
  float val;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      val = mIn->rptr[row][col];
      mOut->rptr[row][col] = val * val;
    }
  }
  return (mOut);
}

MATRIX *MatrixSqrtElts(MATRIX *mIn, MATRIX *mOut)
{
  int row, col, rows, cols;
  float val;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      val = mIn->rptr[row][col];
      mOut->rptr[row][col] = sqrt(val);
    }
  }
  return (mOut);
}

MATRIX *MatrixSignedSquareElts(MATRIX *mIn, MATRIX *mOut)
{
  int row, col, rows, cols;
  float val;

  if (!mOut) {
    mOut = MatrixAlloc(mIn->rows, mIn->cols, mIn->type);
    if (!mOut) return (NULL);
  }

  rows = mIn->rows;
  cols = mIn->cols;

  if ((rows != mOut->rows) || (cols != mOut->cols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixSquareElts: incompatable matrices %d x %d != %d x %d\n",
                 rows,
                 cols,
                 mOut->rows,
                 mOut->cols));

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      val = mIn->rptr[row][col];
      mOut->rptr[row][col] = val * val * ((val < 0) ? -1 : 1);
    }
  }
  return (mOut);
}

MATRIX *MatrixMakeDiagonal(MATRIX *mSrc, MATRIX *mDst)
{
  int row, rows, col, cols;

  if (!mDst) mDst = MatrixClone(mSrc);

  rows = mSrc->rows;
  cols = mSrc->cols;

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      if (row == col)
        mDst->rptr[row][col] = mSrc->rptr[row][col];
      else
        mDst->rptr[row][col] = 0.0f;
    }
  }
  return (mDst);
}

/*
  mDiag is a column vector.
*/
MATRIX *MatrixDiag(MATRIX *mDiag, MATRIX *mOut)
{
  int row, rows, col, cols, nout;

  rows = mDiag->rows;
  cols = mDiag->cols;
  nout = MAX(rows, cols);

  if (!mOut) {
    mOut = MatrixAlloc(nout, nout, mDiag->type);
    if (!mOut) return (NULL);
  }
  else
    MatrixClear(mOut);

  if ((nout != mOut->rows) || (nout != mOut->cols))
    ErrorReturn(
        NULL,
        (ERROR_BADPARM, "MatrixDiag: incompatable matrices %d x %d != %d x %d\n", nout, nout, mOut->rows, mOut->cols));

  if (rows != 1) {
    // column vector
    for (row = 1; row <= rows; row++) mOut->rptr[row][row] = mDiag->rptr[row][1];
  }
  else {
    // row vector
    for (col = 1; col <= cols; col++) mOut->rptr[col][col] = mDiag->rptr[1][col];
  }

  return (mOut);
}

MATRIX *MatrixCopyRegion(const MATRIX *mSrc,
                         MATRIX *mDst,
                         const int start_row,
                         const int start_col,
                         const int rows,
                         const int cols,
                         const int dest_row,
                         const int dest_col)
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCopyRegion: bad destination (%d,%d)\n", dest_row, dest_col));

  srows = mSrc->rows;
  scols = mSrc->cols;
  end_row = start_row + rows - 1;
  end_col = start_col + cols - 1;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) || (start_col > scols) || (end_row > srows) ||
      (end_col > scols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyRegion: bad source region (%d,%d) --> (%d,%d)\n",
                 start_row,
                 start_col,
                 end_row,
                 end_col));

  if (!mDst) mDst = MatrixAlloc(rows, cols, mSrc->type);

  drows = mDst->rows;
  dcols = mDst->cols;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCopyRegion: destination matrix not large enough (%dx%d)\n", rows, cols));

  for (drow = dest_row, srow = start_row; srow <= end_row; srow++, drow++) {
    for (dcol = dest_col, scol = start_col; scol <= end_col; scol++, dcol++) {
      switch (mDst->type) {
        case MATRIX_REAL:
          *MATRIX_RELT(mDst, drow, dcol) = *MATRIX_RELT(mSrc, srow, scol);
          break;
        case MATRIX_COMPLEX:
          *MATRIX_CELT(mDst, drow, dcol) = *MATRIX_CELT(mSrc, srow, scol);
          break;
      }
    }
  }
  return (mDst);
}

MATRIX *MatrixCopyRealRegion(const MATRIX *mSrc,
                             MATRIX *mDst,
                             const int start_row,
                             const int start_col,
                             const int rows,
                             const int cols,
                             const int dest_row,
                             const int dest_col)
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCopyRealRegion: bad destination (%d,%d)\n", dest_row, dest_col));

  srows = mSrc->rows;
  scols = mSrc->cols;
  end_row = start_row + rows - 1;
  end_col = start_col + cols - 1;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) || (start_col > scols) || (end_row > srows) ||
      (end_col > scols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyRealRegion: bad source region (%d,%d) --> (%d,%d)\n",
                 start_row,
                 start_col,
                 end_row,
                 end_col));

  if (!mDst) mDst = MatrixAlloc(rows, cols, mSrc->type);

  drows = mDst->rows;
  dcols = mDst->cols;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixCopyRealRegion: destination matrix not large enough (%dx%d)\n", rows, cols));

  for (drow = dest_row, srow = start_row; srow <= end_row; srow++, drow++) {
    for (dcol = dest_col, scol = start_col; scol <= end_col; scol++, dcol++) {
      switch (mDst->type) {
        case MATRIX_REAL:
          *MATRIX_RELT(mDst, drow, dcol) = MATRIX_CELT_REAL(mSrc, srow, scol);
          break;
        case MATRIX_COMPLEX:
          MATRIX_CELT_IMAG(mDst, drow, dcol) = MATRIX_CELT_IMAG(mSrc, srow, scol);
          break;
      }
    }
  }
  return (mDst);
}

MATRIX *MatrixCopyImagRegion(const MATRIX *mSrc,
                             MATRIX *mDst,
                             const int start_row,
                             const int start_col,
                             const int rows,
                             const int cols,
                             const int dest_row,
                             const int dest_col)
{
  int srow, scol, drow, dcol, srows, scols, drows, dcols, end_row, end_col;

  if ((dest_col < 1) || (dest_row < 1))
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCopyImagRegion: bad destination (%d,%d)\n", dest_row, dest_col));

  srows = mSrc->rows;
  scols = mSrc->cols;
  end_row = start_row + rows - 1;
  end_col = start_col + cols - 1;
  if ((start_row < 1) || (start_row > srows) || (start_col < 1) || (start_col > scols) || (end_row > srows) ||
      (end_col > scols))
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixCopyImagRegion: bad source region (%d,%d) --> (%d,%d)\n",
                 start_row,
                 start_col,
                 end_row,
                 end_col));

  if (!mDst) mDst = MatrixAlloc(rows, cols, mSrc->type);

  drows = mDst->rows;
  dcols = mDst->cols;
  if ((rows > drows) || (cols > dcols))
    ErrorReturn(NULL,
                (ERROR_BADPARM, "MatrixCopyImagRegion: destination matrix not large enough (%dx%d)\n", rows, cols));

  for (drow = dest_row, srow = start_row; srow <= end_row; srow++, drow++) {
    for (dcol = dest_col, scol = start_col; scol <= end_col; scol++, dcol++) {
      switch (mDst->type) {
        case MATRIX_REAL:
          *MATRIX_RELT(mDst, drow, dcol) = MATRIX_CELT_IMAG(mSrc, srow, scol);
          break;
        case MATRIX_COMPLEX:
          MATRIX_CELT_IMAG(mDst, drow, dcol) = MATRIX_CELT_IMAG(mSrc, srow, scol);
          break;
      }
    }
  }
  return (mDst);
}

/*--------------------------------------------------------------------*/
MATRIX *MatrixSetRegion(MATRIX *mSrc, MATRIX *mDst, int start_row, int start_col, int rows, int cols, float val)
{
  int r, c, end_row, end_col;

  end_row = start_row + rows - 1;
  end_col = start_col + cols - 1;

  // printf("rows: %d %d (%d)\n",start_row,end_row,rows);
  // printf("cols: %d %d (%d)\n",start_col,end_col,cols);

  if (start_row > mSrc->rows) {
    printf("ERROR: start row > nrows\n");
    return (NULL);
  }
  if (end_row > mSrc->rows) {
    printf("ERROR: end row > nrows\n");
    return (NULL);
  }
  if (start_col > mSrc->cols) {
    printf("ERROR: start col > ncols\n");
    return (NULL);
  }
  if (end_col > mSrc->cols) {
    printf("ERROR: end col > ncols\n");
    return (NULL);
  }

  if (mDst != NULL) {
    if (mDst->rows != mSrc->rows || mDst->cols != mSrc->cols) {
      printf("ERROR: dimension mismatch\n");
      return (NULL);
    }
  }
  mDst = MatrixCopy(mSrc, mDst);

  for (r = start_row; r <= end_row; r++)
    for (c = start_col; c <= end_col; c++) mDst->rptr[r][c] = val;

  return (mDst);
}

/*--------------------------------------------------------------------*/
MATRIX *MatrixRealToComplex(MATRIX *mReal, MATRIX *mImag, MATRIX *mOut)
{
  int rows, cols, row, col;

  rows = mReal->rows;
  cols = mReal->cols;
  if (!mOut) mOut = MatrixAlloc(mReal->rows, mReal->cols, MATRIX_COMPLEX);

  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      MATRIX_CELT_REAL(mOut, row, col) = *MATRIX_RELT(mReal, row, col);
      MATRIX_CELT_IMAG(mOut, row, col) = *MATRIX_RELT(mImag, row, col);
    }
  }

  return (mOut);
}

float MatrixDeterminant(MATRIX *mIn) { return OpenMatrixDeterminant(mIn); }

/* return the eigenvalues and eigenvectors of the symmetric matrix m in
   evalues and m_dst (columns are vectors) respectively. Note that
   the eigenvalues and eigenvectors are sorted in descending order of
   eigenvalue size.
*/

static int compare_evalues(const void *l1, const void *l2);

typedef struct
{
  int eno;
  float evalue;
} EIGEN_VALUE, EVALUE;

static int compare_evalues(const void *l1, const void *l2)
{
  EVALUE *e1, *e2;

  e1 = (EVALUE *)l1;
  e2 = (EVALUE *)l2;
  return (fabs(e1->evalue) < fabs(e2->evalue) ? 1 : -1);
}

MATRIX *MatrixCalculateEigenSystemHelper(MATRIX *m, float *evalues, MATRIX *m_evectors, int isSymmetric)
{
  int col, i, nevalues, row;
  EVALUE *eigen_values;
  MATRIX *mTmp;

  // sanity-check: input must be n-by-n
  if (m->rows != m->cols) return NULL;

  nevalues = m->rows;
  eigen_values = (EVALUE *)calloc((unsigned int)nevalues, sizeof(EIGEN_VALUE));
  if (!m_evectors) m_evectors = MatrixAlloc(m->rows, m->cols, MATRIX_REAL);

  mTmp = MatrixAlloc(m->rows, m->cols, MATRIX_REAL);

  if (isSymmetric) {
    if (OpenEigenSystem(m->data, m->rows, evalues, mTmp->data) != NO_ERROR) return (NULL);
  }
  else {
    if (OpenNonSymmetricEigenSystem(m->data, m->rows, evalues, mTmp->data) != NO_ERROR) return (NULL);
  }

  /*
     sort eigenvalues in order of decreasing absolute value. The
     EIGEN_VALUE structure is needed to record final order of eigenvalue so
     we can also sort the eigenvectors.
  */
  for (i = 0; i < nevalues; i++) {
    eigen_values[i].eno = i;
    eigen_values[i].evalue = evalues[i];
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues);
  for (i = 0; i < nevalues; i++) evalues[i] = eigen_values[i].evalue;

  /* now sort the eigenvectors */
  for (col = 0; col < mTmp->cols; col++) {
    for (row = 1; row <= mTmp->rows; row++) m_evectors->rptr[row][col + 1] = mTmp->rptr[row][eigen_values[col].eno + 1];
  }

  free(eigen_values);
  MatrixFree(&mTmp);
  return (m_evectors);
}

int MatrixIsSymmetric(MATRIX *matrix)
{
  int row;
  int col;
  int isSymmetric = 1;

  if (matrix->rows != matrix->cols) return 0;  // non-square, so not symmetric

  for (row = 1; row <= matrix->rows; row++) {
    for (col = 1; col <= matrix->cols; col++) {
      if (matrix->rptr[row][col] != matrix->rptr[col][row]) {
        isSymmetric = 0;
      }
    }
  }

  return isSymmetric;
}

MATRIX *MatrixEigenSystem(MATRIX *m, float *evalues, MATRIX *m_evectors)
{
  int isSymmetric = MatrixIsSymmetric(m);

  return MatrixCalculateEigenSystemHelper(m, evalues, m_evectors, isSymmetric);
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
  // mV = MatrixIdentity(mA->rows, mV) ;
  // svd(mA->rptr, mV->rptr, v_z->data, mA->rows, mA->cols) ;

  if (mV == NULL) mV = MatrixAlloc(mA->rows, mA->rows, MATRIX_REAL);
  OpenSvdcmp(mA, v_z, mV);

  return (mV);
}

float MatrixSVDEigenValues(MATRIX *m, float *evalues)
{
  float cond;
  VECTOR *v_w;
  MATRIX *m_U, *m_V;
  int row, rows, cols, nevalues, i;
  float wmax, wmin, wi;
  EVALUE *eigen_values;

  cols = m->cols;
  rows = m->rows;
  m_U = MatrixCopy(m, NULL);
  v_w = RVectorAlloc(cols, MATRIX_REAL);
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL);
  nevalues = m->rows;
  memset(evalues, 0, nevalues * sizeof(evalues[0]));

  /* calculate condition # of matrix */
  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR) return (Gerror);

  eigen_values = (EVALUE *)calloc((unsigned int)nevalues, sizeof(EIGEN_VALUE));
  for (i = 0; i < nevalues; i++) {
    eigen_values[i].eno = i;
    eigen_values[i].evalue = RVECTOR_ELT(v_w, i + 1);
  }
  qsort((char *)eigen_values, nevalues, sizeof(EVALUE), compare_evalues);
  for (i = 0; i < nevalues; i++) evalues[i] = eigen_values[i].evalue;

  wmax = 0.0f;
  wmin = wmax = RVECTOR_ELT(v_w, 1);
  for (row = 2; row <= rows; row++) {
    wi = fabs(RVECTOR_ELT(v_w, row));
    if (wi > wmax) wmax = wi;
    if (wi < wmin) wmin = wi;
  }

  if (FZERO(wmin))
    cond = 1e8; /* something big */
  else
    cond = wmax / wmin;

  free(eigen_values);
  MatrixFree(&m_U);
  VectorFree(&v_w);
  MatrixFree(&m_V);
  return (cond);
}

/*
  use SVD to find the inverse of m, usually useful if m is
  ill-conditioned or singular.
*/
#define TOO_SMALL 1e-4

MATRIX *MatrixSVDInverse(MATRIX *m, MATRIX *m_inverse)
{
  // This function accounted for a lot of the matrix allocation and free, and they were almost always 1x3 or 3x3
  // so it has been rewritten to use the buffered matrix support
  //
  if (MatrixIsZero(m)) return (NULL);

  int const rows = m->rows;
  int const cols = m->cols;

  MatrixBuffer m_U_buffer, m_V_buffer;
  MATRIX      *m_U,       *m_V;

  m_U = MatrixAlloc2(rows, cols, MATRIX_REAL, &m_U_buffer);
  m_U = MatrixCopy(m, m_U);

  m_V = MatrixAlloc2(cols, cols, MATRIX_REAL, &m_V_buffer);

  MatrixBuffer v_w_buffer; 
  MATRIX* v_w = MatrixAlloc2(1, cols, MATRIX_REAL, &v_w_buffer);

  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR) {
    MatrixFree(&m_U);
    MatrixFree(&m_V);
    MatrixFree(&v_w);
    return (NULL);
  }

  MatrixBuffer m_w_buffer;
  MATRIX *m_w = MatrixAlloc2(cols, cols, MATRIX_REAL, &m_w_buffer);
 
  float wmax = 0.0f;
  int row;
  for (row = 1; row <= rows; row++)
    if (fabs(RVECTOR_ELT(v_w, row)) > wmax) wmax = fabs(RVECTOR_ELT(v_w, row));

  float wmin = TOO_SMALL * wmax;
  for (row = 1; row <= rows; row++) {
    if (fabs(RVECTOR_ELT(v_w, row)) < wmin)
      m_w->rptr[row][row] = 0.0f;
    else
      m_w->rptr[row][row] = 1.0f / RVECTOR_ELT(v_w, row);
  }

  // This is the WRONG way to multiply by a transpose!
  // It is faster to not transpose and perform a faster multiply
  //
  // m_w  is (cols x cols)
  // m_U  is (rows x cols)
  // m_Ut is (cols x rows)
  // (cols x cols) * (cols x rows) is (cols, rows) which is what m_tmp needs to be
  //
  MatrixBuffer m_Ut_buffer, m_tmp_buffer;
  
  MATRIX *m_Ut  = MatrixAlloc2(cols, cols, MATRIX_REAL, &m_Ut_buffer);
  m_Ut = MatrixTranspose(m_U, m_Ut);

  MATRIX *m_tmp = MatrixAlloc2(cols, rows, MATRIX_REAL, &m_tmp_buffer);
  m_tmp         = MatrixMultiply (m_w, m_Ut,  m_tmp);

  // m_V   is (cols x cols)
  // m_tmp is (cols x rows)
  // m_inverse is (cols x rows) but it is returned so can't be buffered here
  //
  m_inverse = MatrixMultiply (m_V, m_tmp, m_inverse);

  MatrixFree(&m_U);
  MatrixFree(&v_w);
  MatrixFree(&m_V);
  MatrixFree(&m_w);

  MatrixFree(&m_Ut);
  MatrixFree(&m_tmp);
  
  return (m_inverse);
}

MATRIX *MatrixAllocTranslation(int n, double *trans)
{
  MATRIX *mat;
  int i;

  mat = MatrixIdentity(n, NULL);
  for (i = 0; i < n - 1; i++) *MATRIX_RELT(mat, i + 1, 4) = trans[i];
  return (mat);
}

MATRIX *MatrixReallocRotation(int n, float angle, int which, MATRIX *m)
{
  float s, c;

  if (!m) m = MatrixIdentity(n, NULL);

  c = cos(angle);
  s = sin(angle);
  switch (which) {
    case X_ROTATION:
      m->rptr[2][2] = c;
      m->rptr[2][3] = s;
      m->rptr[3][2] = -s;
      m->rptr[3][3] = c;
      break;
    case Y_ROTATION:
      m->rptr[1][1] = c;
      m->rptr[1][3] = -s;
      m->rptr[3][1] = s;
      m->rptr[3][3] = c;
      break;
    case Z_ROTATION:
      m->rptr[1][1] = c;
      m->rptr[1][2] = s;
      m->rptr[2][1] = -s;
      m->rptr[2][2] = c;
      break;
  }

  return (m);
}

MATRIX *MatrixAllocRotation(int n, float angle, int which) { return (MatrixReallocRotation(n, angle, which, NULL)); }

MATRIX *MatrixCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mMeans)
{
  int ninputs, nvars, input, var, var2;
  float *means, covariance, obs1, obs2;

  ninputs = mInputs->rows; /* number of observations */
  nvars = mInputs->cols;   /* number of variables */

  if (!ninputs) {
    if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);
    return (mCov);
  }

  if (!nvars) ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCovariance: zero size input"));

  if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);

  if (!mCov) ErrorExit(ERROR_NO_MEMORY, "MatrixCovariance: could not allocate %d x %d covariance matrix", nvars, nvars);

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCovariance: incorrect covariance matrix dimensions"));

  if (!mMeans)
    means = (float *)calloc(nvars, sizeof(float));
  else
    means = mMeans->data;

  for (var = 0; var < nvars; var++) {
    means[var] = 0.0f;
    for (input = 0; input < ninputs; input++) means[var] += mInputs->rptr[input + 1][var + 1];
    means[var] /= ninputs;
  }

  for (var = 0; var < nvars; var++) {
    for (var2 = 0; var2 < nvars; var2++) {
      covariance = 0.0f;
      for (input = 0; input < ninputs; input++) {
        obs1 = mInputs->rptr[input + 1][var + 1];
        obs2 = mInputs->rptr[input + 1][var2 + 1];
        covariance += (obs1 - means[var]) * (obs2 - means[var2]);
      }
      covariance /= ninputs;
      mCov->rptr[var + 1][var2 + 1] = mCov->rptr[var2 + 1][var + 1] = covariance;
    }
  }

  if (!mMeans) free(means);
  return (mCov);
}

/*
  update the values of a covariance matrix. Note that MatrixFinalCovariance
  must be called with the means and total # of inputs before the covariance
  matrix is valid.
*/
MATRIX *MatrixUpdateCovariance(MATRIX *mInputs, MATRIX *mCov, MATRIX *mMeans)
{
  int ninputs, nvars, input, var, var2;
  float covariance, obs1, obs2, mean1, mean2;

  ninputs = mInputs->rows; /* number of observations */
  nvars = mInputs->cols;   /* number of variables */

  if (!mMeans || mMeans->rows != nvars)
    ErrorReturn(
        NULL, (ERROR_BADPARM, "MatrixUpdateCovariance: wrong mean vector size (%d x %d)", mMeans->rows, mMeans->cols));

  if (!ninputs) {
    if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);
    return (mCov);
  }

  if (!nvars) ErrorReturn(NULL, (ERROR_BADPARM, "MatrixUpdateCovariance: zero size input"));

  if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);

  if (!mCov)
    ErrorExit(ERROR_NO_MEMORY,
              "MatrixUpdateCovariance: could not allocate %d x %d covariance "
              "matrix",
              nvars,
              nvars);

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL,
                (ERROR_BADPARM,
                 "MatrixUpdateCovariance: incorrect covariance matrix "
                 "dimensions"));

  for (var = 0; var < nvars; var++) {
    mean1 = mMeans->rptr[var + 1][1];
    for (var2 = 0; var2 < nvars; var2++) {
      mean2 = mMeans->rptr[var2 + 1][1];
      covariance = 0.0f;
      for (input = 0; input < ninputs; input++) {
        obs1 = mInputs->rptr[input + 1][var + 1];
        obs2 = mInputs->rptr[input + 1][var2 + 1];
        covariance += (obs1 - mean1) * (obs2 - mean2);
      }

      mCov->rptr[var + 1][var2 + 1] = mCov->rptr[var2 + 1][var + 1] = covariance;
    }
  }

  return (mCov);
}

/*
  update the means based on a new set of observation vectors. Notet that
  the user must keep track of the total # observations
*/
MATRIX *MatrixUpdateMeans(MATRIX *mInputs, MATRIX *mMeans, VECTOR *mNobs)
{
  int ninputs, nvars, input, var;
  float *means;

  ninputs = mInputs->rows; /* number of observations */
  nvars = mInputs->cols;   /* number of variables */
  if (!mMeans) mMeans = VectorAlloc(nvars, MATRIX_REAL);
  means = mMeans->data;
  for (var = 0; var < nvars; var++) {
    mNobs->rptr[var + 1][1] += ninputs;
    for (input = 0; input < ninputs; input++) means[var] += mInputs->rptr[input + 1][var + 1];
  }
  return (mMeans);
}

/*
  compute the final mean vector previously generated by
  MatrixUpdateMeans. The user must supply a vector containing
  the # of observations of each variable.
*/
MATRIX *MatrixFinalMeans(VECTOR *mMeans, VECTOR *mNobs)
{
  int nvars, var, n;
  float *means;

  nvars = mMeans->rows; /* number of variables */
  means = mMeans->data;
  for (var = 0; var < nvars; var++) {
    n = mNobs->rptr[var + 1][1];
    if (!n)
      means[var] = 0.0f;
    else
      means[var] /= (float)n;
  }
  return (mMeans);
}

/*
  compute the final covariance matrix previously generated by
  MatrixUpdateCovariance. The user must supply a vector containing
  the # of observations per variable.
*/
MATRIX *MatrixFinalCovariance(MATRIX *mInputs, MATRIX *mCov, VECTOR *mNobs)
{
  int ninputs, nvars, var, var2;
  float covariance;

  ninputs = mInputs->rows; /* number of observations */
  nvars = mInputs->cols;   /* number of variables */

  if (!ninputs) {
    if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);
    return (mCov);
  }

  if (!nvars) ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCovariance: zero size input"));

  if (!mCov) mCov = MatrixAlloc(nvars, nvars, MATRIX_REAL);

  if (!mCov) ErrorExit(ERROR_NO_MEMORY, "MatrixCovariance: could not allocate %d x %d covariance matrix", nvars, nvars);

  if (mCov->cols != nvars || mCov->rows != nvars)
    ErrorReturn(NULL, (ERROR_BADPARM, "MatrixCovariance: incorrect covariance matrix dimensions"));

  for (var = 0; var < nvars; var++) {
    ninputs = mCov->rptr[var + 1][1];
    for (var2 = 0; var2 < nvars; var2++) {
      covariance = mCov->rptr[var + 1][var2 + 1];
      if (ninputs)
        covariance /= ninputs;
      else
        covariance = 0.0f;
      mCov->rptr[var + 1][var2 + 1] = mCov->rptr[var2 + 1][var + 1] = covariance;
    }
  }

  return (mCov);
}

int MatrixAsciiWrite(const char *fname, MATRIX *m)
{
  FILE *fp;
  int ret;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MatrixAsciiWrite: could not open file %s", fname));
  ret = MatrixAsciiWriteInto(fp, m);
  fclose(fp);
  return (ret);
}

MATRIX *MatrixAsciiRead(const char *fname, MATRIX *m)
{
  FILE *fp;

  fp = fopen(fname, "r");
  if (!fp) ErrorReturn(NULL, (ERROR_NO_FILE, "MatrixAsciiRead: could not open file %s", fname));
  m = MatrixAsciiReadFrom(fp, m);
  fclose(fp);
  return (m);
}

int MatrixAsciiWriteInto(FILE *fp, MATRIX *m)
{
  int row, col;

  fprintf(fp, "%d %d %d\n", m->type, m->rows, m->cols);
  for (row = 1; row <= m->rows; row++) {
    for (col = 1; col <= m->cols; col++) {
      if (m->type == MATRIX_COMPLEX)
        fprintf(fp, "%+f %+f   ", MATRIX_CELT_REAL(m, row, col), MATRIX_CELT_IMAG(m, row, col));
      else
        fprintf(fp, "%+f  ", m->rptr[row][col]);
    }
    fprintf(fp, "\n");
  }
  return (NO_ERROR);
}

MATRIX *MatrixAsciiReadFrom(FILE *fp, MATRIX *m)
{
  int type, rows, cols, row, col;
  char *cp, line[200];

  cp = fgetl(line, 199, fp);
  if (!cp) ErrorReturn(NULL, (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scanf parms"));
  if (sscanf(cp, "%d %d %d\n", &type, &rows, &cols) != 3)
    ErrorReturn(NULL, (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scanf parms"));

  if (!m) {
    m = MatrixAlloc(rows, cols, type);
    if (!m) ErrorReturn(NULL, (ERROR_BADFILE, "MatrixAsciiReadFrom: could not allocate matrix"));
  }
  else {
    if (m->rows != rows || m->cols != cols || m->type != type)
      ErrorReturn(m, (ERROR_BADFILE, "MatrixAsciiReadFrom: specified matrix does not match file"));
  }
  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) {
      if (m->type == MATRIX_COMPLEX) {
        if (fscanf(fp, "%f %f   ", &MATRIX_CELT_REAL(m, row, col), &MATRIX_CELT_IMAG(m, row, col)) != 2)
          ErrorReturn(NULL, (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scan element (%d, %d)", row, col));
      }
      else if (fscanf(fp, "%f  ", &m->rptr[row][col]) != 1)
        ErrorReturn(NULL, (ERROR_BADFILE, "MatrixAsciiReadFrom: could not scan element (%d, %d)", row, col));
    }
    // a non zero return value for fscanf would be odd, since the scan does not match any values
    // but we need to handle the return value, so this seems to be the sensible way to 'contract this'
    // even if ferror is never reached
    if(fscanf(fp, "\n") != 0 && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "MatrixAsciiReadFrom: A read error occured while scanning for newline");
    }
  }

  return (m);
}

int MatrixWriteInto(FILE *fp, MATRIX *m)
{
  int row, col;

  fwriteInt(m->type, fp);
  fwriteInt(m->rows, fp);
  fwriteInt(m->cols, fp);

  for (row = 1; row <= m->rows; row++) {
    for (col = 1; col <= m->cols; col++) {
      if (m->type == MATRIX_COMPLEX) {
        fwriteDouble(MATRIX_CELT_REAL(m, row, col), fp);
        fwriteDouble(MATRIX_CELT_IMAG(m, row, col), fp);
      }
      else
        fwriteDouble(m->rptr[row][col], fp);
    }
  }
  return (NO_ERROR);
}

MATRIX *MatrixReadFrom(FILE *fp, MATRIX *m)
{
  int row, col, rows, cols, type;

  type = freadInt(fp);
  rows = freadInt(fp);
  cols = freadInt(fp);

  if (!m) {
    m = MatrixAlloc(rows, cols, type);
    if (!m) ErrorReturn(NULL, (ERROR_BADFILE, "MatrixReadFrom: could not allocate matrix"));
  }
  else {
    if (m->rows != rows || m->cols != cols || m->type != type)
      ErrorReturn(m, (ERROR_BADFILE, "MatrixReadFrom: specified matrix does not match file"));
  }

  for (row = 1; row <= m->rows; row++) {
    for (col = 1; col <= m->cols; col++) {
      if (m->type == MATRIX_COMPLEX) {
        MATRIX_CELT_REAL(m, row, col) = freadDouble(fp);
        MATRIX_CELT_IMAG(m, row, col) = freadDouble(fp);
      }
      else
        m->rptr[row][col] = freadDouble(fp);
    }
  }

  return (m);
}

/*
  calculate and return the Euclidean norm of the vector v.
*/
float VectorLen(const VECTOR *v)
{
  int i;
  float len, vi;

  for (len = 0.0f, i = 1; i <= v->rows; i++) {
    vi = v->rptr[i][1];
    len += vi * vi;
  }
  len = sqrt(len);
  return (len);
}

/*  compute the dot product of 2 vectors */
float VectorDot(const VECTOR *v1, const VECTOR *v2)
{
  int i;
  float dot;

  for (dot = 0.0f, i = 1; i <= v1->rows; i++) dot += v1->rptr[i][1] * v2->rptr[i][1];
  ;

  return (dot);
}

/*  compute the dot product of 2 vectors */
float VectorNormalizedDot(VECTOR *v1, VECTOR *v2)
{
  float dot, l1, l2;

  l1 = VectorLen(v1);
  l2 = VectorLen(v2);
  dot = VectorDot(v1, v2);
  if (FZERO(l1) || FZERO(l2)) return (0.0f);

  return (dot / (l1 * l2));
}

/*
  extract a column of the matrix m and return it in the
  vector v.
*/
VECTOR *MatrixColumn(MATRIX *m, VECTOR *v, int col)
{
  int row;

  if (!v) v = VectorAlloc(m->rows, MATRIX_REAL);

  for (row = 1; row <= m->rows; row++) VECTOR_ELT(v, row) = m->rptr[row][col];
  return (v);
}

/* calcuate and return the Euclidean distance between two vectors */
float VectorDistance(VECTOR *v1, VECTOR *v2)
{
  int row;
  float dist, d;

  for (dist = 0.0f, row = 1; row <= v1->rows; row++) {
    d = VECTOR_ELT(v1, row) - VECTOR_ELT(v2, row);
    dist += (d * d);
  }

  return (dist);
}

/*
  compute the outer product of the vectors v1 and v2.
*/
MATRIX *VectorOuterProduct(VECTOR *v1, VECTOR *v2, MATRIX *m)
{
  int row, col, rows, cols;
  float r;

  rows = v1->rows;
  cols = v2->rows;

  if (rows != cols) ErrorReturn(NULL, (ERROR_BADPARM, "VectorOuterProduct: v1->rows %d != v2->rows %d", rows, cols));

  if (!m) m = MatrixAlloc(rows, cols, MATRIX_REAL);

  for (row = 1; row <= rows; row++) {
    r = VECTOR_ELT(v1, row);
    for (col = 1; col <= cols; col++) m->rptr[row][col] = r * VECTOR_ELT(v2, col);
  }

  return (m);
}

/*
  Add a small random diagonal matrix to mIn to make it
  non-singular.
*/
#define SMALL 1e-4
MATRIX *MatrixRegularize(MATRIX *mIn, MATRIX *mOut)
{
  int rows, cols, row;
  float ran_num;

  rows = mIn->rows;
  cols = mIn->cols;
  if (!mOut) mOut = MatrixAlloc(rows, cols, mIn->type);

  if (mIn != mOut) MatrixCopy(mIn, mOut);

#if 0
  ran_num = (SMALL*randomNumber(0.0, 1.0)+SMALL) ;
#else
  ran_num = SMALL;
#endif

  for (row = 1; row <= rows; row++) mOut->rptr[row][row] += ran_num;

  return (mOut);
}

/* see if a matrix is singular */
int MatrixSingular(MATRIX *m)
{
#if 1
  VECTOR *v_w;
  MATRIX *m_U, *m_V;
  int row, rows, cols;
  float wmax, wmin, wi;

  if (MatrixIsZero(m)) return (1);
  cols = m->cols;
  rows = m->rows;
  m_U = MatrixCopy(m, NULL);
  v_w = RVectorAlloc(cols, MATRIX_REAL);
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL);

  /* calculate condition # of matrix */
  if (OpenSvdcmp(m_U, v_w, m_V) != NO_ERROR) return (Gerror);

  wmax = 0.0f;
  wmin = wmax = RVECTOR_ELT(v_w, 1);
  for (row = 2; row <= rows; row++) {
    wi = fabs(RVECTOR_ELT(v_w, row));
    if (wi > wmax) wmax = wi;
    if (wi < wmin) wmin = wi;
  }

  MatrixFree(&m_U);
  VectorFree(&v_w);
  MatrixFree(&m_V);
  return (FZERO(wmin) ? 1 : wmin < wmax * TOO_SMALL);
#else
  float det;

  det = MatrixDeterminant(m);
  return (FZERO(det));
#endif
}

int MatrixIsZero(MATRIX *m)
{
  int row, col;
  double val;

  for (row = 1; row <= m->rows; row++)
    for (col = 1; col <= m->cols; col++) {
      val = fabs(*MATRIX_RELT(m, row, col));
      if (val > 1e-11) {
        if (val < 1e-10) DiagBreak();
        return (0);
      }
    }

  return (1);
}

/*
  calcluate the condition # of a matrix using svd
*/
float MatrixConditionNumber(MATRIX *m)
{
  float cond;
  VECTOR *v_w;
  MATRIX *m_U, *m_V;
  int row, rows, cols;
  float wmax, wmin, wi;

  if (MatrixIsZero(m)) return (1e10);
  cols = m->cols;
  rows = m->rows;
  m_U = MatrixCopy(m, NULL);
  v_w = RVectorAlloc(cols, MATRIX_REAL);
  m_V = MatrixAlloc(cols, cols, MATRIX_REAL);

  /* calculate condition # of matrix */
  OpenSvdcmp(m_U, v_w, m_V);
  wmax = 0.0f;
  wmin = wmax = RVECTOR_ELT(v_w, 1);
  for (row = 2; row <= rows; row++) {
    wi = fabs(RVECTOR_ELT(v_w, row));
    if (wi > wmax) wmax = wi;
    if (wi < wmin) wmin = wi;
  }

  if (FZERO(wmin))
    cond = 1e8; /* something big */
  else
    cond = wmax / wmin;

  MatrixFree(&m_U);
  VectorFree(&v_w);
  MatrixFree(&m_V);
  return (cond);
}

/*-----------------------------------------------------------
  MatrixNSConditionNumber() - condition of a non-square matrix.
  Works for square matrices as well.
  -----------------------------------------------------------*/
float MatrixNSConditionNumber(MATRIX *m)
{
  float cond;
  MATRIX *mt, *p;

  if (m->rows == m->cols) return (MatrixConditionNumber(m));

  mt = MatrixTranspose(m, NULL);

  if (m->rows > m->cols)
    p = MatrixMultiply(mt, m, NULL);
  else
    p = MatrixMultiply(m, mt, NULL);

  cond = MatrixConditionNumber(p);

  MatrixFree(&mt);
  MatrixFree(&p);

  return (cond);
}

/*
  calculate the cross product of two vectors and return the result
  in vdst, allocating it if necessary.
*/
VECTOR *VectorCrossProduct(const VECTOR *v1, const VECTOR *v2, VECTOR *vdst)
{
  float x1, x2, y1, y2, z1, z2;

  if (v1->rows != 3 && v1->cols != 1) ErrorReturn(NULL, (ERROR_BADPARM, "VectorCrossProduct: must be 3-vectors"));

  if (!vdst) vdst = VectorClone(v1);

  x1 = VECTOR_ELT(v1, 1);
  y1 = VECTOR_ELT(v1, 2);
  z1 = VECTOR_ELT(v1, 3);
  x2 = VECTOR_ELT(v2, 1);
  y2 = VECTOR_ELT(v2, 2);
  z2 = VECTOR_ELT(v2, 3);

  VECTOR_ELT(vdst, 1) = y1 * z2 - z1 * y2;
  VECTOR_ELT(vdst, 2) = z1 * x2 - x1 * z2;
  VECTOR_ELT(vdst, 3) = x1 * y2 - y1 * x2;
  return (vdst);
}

/*!
  \fn VECTOR *VectorCrossProductD(const VECTOR *v1, const VECTOR *v2, VECTOR *vdst)
  Calculates the cross product of two vectors and returns the result
  in vdst, allocating it if necessary. Same as VectorCrossProduct() but uses
  long doubles internally instead of floats.
*/
VECTOR *VectorCrossProductD(const VECTOR *v1, const VECTOR *v2, VECTOR *vdst)
{
  long double x1, x2, y1, y2, z1, z2;

  if (v1->rows != 3 && v1->cols != 1) return(NULL);

  if (!vdst) vdst = VectorClone(v1);

  x1 = VECTOR_ELT(v1, 1);
  y1 = VECTOR_ELT(v1, 2);
  z1 = VECTOR_ELT(v1, 3);
  x2 = VECTOR_ELT(v2, 1);
  y2 = VECTOR_ELT(v2, 2);
  z2 = VECTOR_ELT(v2, 3);

  VECTOR_ELT(vdst, 1) = y1 * z2 - z1 * y2;
  VECTOR_ELT(vdst, 2) = z1 * x2 - x1 * z2;
  VECTOR_ELT(vdst, 3) = x1 * y2 - y1 * x2;
  return (vdst);
}


/*
  compute the triple scalar product v1 x v2 . v3
*/
float VectorTripleProduct(const VECTOR *v1, const VECTOR *v2, const VECTOR *v3)
{
  float x1, x2, y1, y2, z1, z2, x3, y3, z3, total;

  if (v1->rows != 3 && v1->cols != 1) ErrorReturn(0.0f, (ERROR_BADPARM, "VectorCrossProduct: must be 3-vectors"));

  x1 = VECTOR_ELT(v1, 1);
  y1 = VECTOR_ELT(v1, 2);
  z1 = VECTOR_ELT(v1, 3);
  x2 = VECTOR_ELT(v2, 1);
  y2 = VECTOR_ELT(v2, 2);
  z2 = VECTOR_ELT(v2, 3);
  x3 = VECTOR_ELT(v3, 1);
  y3 = VECTOR_ELT(v3, 2);
  z3 = VECTOR_ELT(v3, 3);

  total = x3 * (y1 * z2 - z1 * y2);
  total += y3 * (z1 * x2 - x1 * z2);
  total += z3 * (x1 * y2 - y1 * x2);
  return (total);
}

VECTOR *VectorNormalize(const VECTOR *vin, VECTOR *vout)
{
  float len;
  int row, col, rows, cols;

  if (!vout) vout = VectorClone(vin);

  len = VectorLen(vin);
  if (FZERO(len)) len = 1.0f; /* doesn't matter - all elements are 0 */

  rows = vin->rows;
  cols = vin->cols;
  for (row = 1; row <= rows; row++) {
    for (col = 1; col <= cols; col++) vout->rptr[row][col] = vin->rptr[row][col] / len;
  }
  return (vout);
}

float VectorAngle(const VECTOR *v1, const VECTOR *v2)
{
  float angle, l1, l2, dot, norm;

  l1 = VectorLen(v1);
  l2 = VectorLen(v2);
  norm = fabs(l1 * l2);
  if (FZERO(norm)) return (0.0f);
  dot = VectorDot(v1, v2);
  if (dot > norm)
    angle = acos(1.0);
  else
    angle = acos(dot / norm);
  return (angle);
}

double Vector3Angle(VECTOR *v1, VECTOR *v2)
{
  double angle, l1, l2, dot, norm, x, y, z;

  x = V3_X(v1);
  y = V3_Y(v1);
  z = V3_Z(v1);
  l1 = sqrt(x * x + y * y + z * z);
  x = V3_X(v2);
  y = V3_Y(v2);
  z = V3_Z(v2);
  l2 = sqrt(x * x + y * y + z * z);
  norm = l1 * l2;
  if (FZERO(norm)) return (0.0f);
  dot = V3_DOT(v1, v2);
  if (fabs(dot) > fabs(norm)) norm = fabs(dot);
  if (dot > norm)
    angle = acos(1.0);
  else
    angle = acos(dot / norm);
  return (angle);
}

void XYZ_NORMALIZED_LOAD(XYZ *xyz, float *xyz_length, float x, float y, float z)
{
  float len = *xyz_length = sqrt((double)x * (double)x + (double)y * (double)y + (double)z * (double)z);

  if (len == 0.0f) {
    xyz->x = 1.0f;
    xyz->y = 0.0f;
    xyz->z = 0.0f;
  }
  else {
    float len_inv = 1.0f / len;
    xyz->x = x * len_inv;
    xyz->y = y * len_inv;
    xyz->z = z * len_inv;
  }
}


float XYZApproxAngle(XYZ const *normalizedXYZ, float x2, float y2, float z2)
{
  float norm =
      (float)sqrt((double)x2 * (double)x2 + (double)y2 * (double)y2 + (double)z2 * (double)z2);

  if (FZERO(norm)) return (0.0f);

  return XYZApproxAngle_knownLength(normalizedXYZ, x2, y2, z2, norm);
}


float XYZApproxAngle_knownLength(
    XYZ const * normalizedXYZ, 
    float x2, float y2, float z2, float norm)
{
  double x1 = normalizedXYZ->x;
  double y1 = normalizedXYZ->y;
  double z1 = normalizedXYZ->z;


  float dot = x1 * x2 + y1 * y2 + z1 * z2;

  if (fabs(dot) > norm) norm = fabs(dot);

  double acosInput = dot / norm;

  // disabling this warning because the code below fails when -fopenmp is not used during compile. 
  // This is an ongoing issue in gcc 4.9.4 - clarsen
  #if GCC_VERSION > 40408
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wunknown-pragmas"
  #endif

  // left the following here for debugging or other investigation
  if (0)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
  {
#define HL 1024
    static long histogram[HL];
    static long population;
    static long populationLimit = 128;
    int i = (int)((acosInput + 1.0f) / 2.0f * (HL - 1));
    histogram[i]++;
    population++;
    if (population == populationLimit) {
      populationLimit *= 2;
      int j, c;
      printf(
          "XYZApproxAngle histogram after %ld and after incrementing %d for acosInput:%g\n", population, i, acosInput);
      c = 0;
      for (j = 0; j < HL; j++) {
        if (!histogram[j]) continue;
        printf("  %g:%ld", (float)j / (float)(HL - 1) * 2.0f - 1.0f, histogram[j]);
        if (++c == 10) {
          c = 0;
          printf("\n");
        }
      }
      if (c > 0) printf("\n");
    }
#undef HL
  }

  // cos(x) = 1 - x**2/2 + ...
  // so when acosInput is near 1, which is almost always is,  we get
  //
  // acosInput ~= 1 - x**2/2
  // x**2/2 = 1 - acosInput
  // x      = sqrt(2.0 * (1 - acosInput))
  //
  float angle;

  if (acosInput < 0.99f)
    angle = acos(acosInput);
  else {
    angle = sqrt(2.0f * (1.0f - acosInput));
    // left the following here for debugging or other investigation
    if (0)
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
    {
      static int count;
      if (count++ < 100) printf("acos approx inp:%g approx:%g correct:%g\n", acosInput, angle, acos(acosInput));
    }
  }

  #if GCC_VERSION > 40408
  #pragma GCC diagnostic pop
  #endif
  
  return (angle);
}

int MatrixWriteTxt(const char *fname, MATRIX *mat)
{
  FILE *fp;
  int row, col;

  fp = fopen(fname, "w");
  if (!fp) ErrorReturn(ERROR_NO_FILE, (ERROR_NO_FILE, "MatrixWriteTxt(%s) - file open failed\n", fname));

  for (row = 1; row <= mat->rows; row++) {
    for (col = 1; col <= mat->cols; col++) fprintf(fp, "%+4.5f ", mat->rptr[row][col]);
    fprintf(fp, "\n");
  }

  fclose(fp);
  return (NO_ERROR);
}

MATRIX *MatrixPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  MATRIX *mT, *mTm, *mTm_inv;

  if (m->rows < m->cols) {
    MATRIX *mT = MatrixTranspose(m, NULL);
    m_pseudo_inv = MatrixPseudoInverse(mT, m_pseudo_inv);
    MatrixFree(&mT);
    mT = m_pseudo_inv;
    m_pseudo_inv = MatrixTranspose(mT, NULL);
    MatrixFree(&mT);
    return (m_pseudo_inv);
  }
  mT = MatrixTranspose(m, NULL);

  /* build (mT m)-1 mT */
  mTm = MatrixMultiply(mT, m, NULL);
  mTm_inv = MatrixInverse(mTm, NULL);
  if (!mTm_inv) {
    mTm_inv = MatrixSVDInverse(mTm, NULL);

    if (!mTm_inv) {
      MatrixFree(&mT);
      MatrixFree(&mTm);
      return (NULL);
    }
  }
  m_pseudo_inv = MatrixMultiply(mTm_inv, mT, m_pseudo_inv);

  MatrixFree(&mT);
  MatrixFree(&mTm);
  MatrixFree(&mTm_inv);
  return (m_pseudo_inv);
}

MATRIX *MatrixRightPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  MATRIX *mT, *mmT, *mmT_inv;

  /* build mT (m mT)-1 */
  mT = MatrixTranspose(m, NULL);
  mmT = MatrixMultiply(m, mT, NULL);
  mmT_inv = MatrixInverse(mmT, NULL);
  if (!mmT_inv) {
    MatrixFree(&mT);
    MatrixFree(&mmT);
    return (NULL);
  }
  m_pseudo_inv = MatrixMultiply(mT, mmT_inv, m_pseudo_inv);

  MatrixFree(&mT);
  MatrixFree(&mmT);
  MatrixFree(&mmT_inv);
  return (m_pseudo_inv);
}

int MatrixCheck(MATRIX *m)
{
  int rows, cols, r, c;

  rows = m->rows;
  cols = m->cols;
  for (r = 1; r <= rows; r++) {
    for (c = 1; c <= cols; c++) {
      if (!std::isfinite(*MATRIX_RELT(m, r, c))) return (ERROR_BADPARM);
    }
  }
  return (NO_ERROR);
}

MATRIX *MatrixReshape(MATRIX *m_src, MATRIX *m_dst, int rows, int cols)
{
  int r1, c1, r2, c2;

  if (m_dst) {
    rows = m_dst->rows;
    cols = m_dst->cols;
  }

#if 0
  if (rows*cols > m_src->rows*m_src->cols)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MatrixReshape: (%d,%d) -> (%d,%d), lengths must be"
                       " equal", m_src->rows, m_src->cols, rows,cols)) ;
#endif
  if (!m_dst) m_dst = MatrixAlloc(rows, cols, m_src->type);

  for (r2 = c2 = r1 = 1; r1 <= m_src->rows; r1++) {
    for (c1 = 1; c1 <= m_src->cols; c1++) {
      *MATRIX_RELT(m_dst, r2, c2) = *MATRIX_RELT(m_src, r1, c1);
      if (++c2 > m_dst->cols) {
        c2 = 1;
        r2++;
      }
      if (r2 > rows) /* only extract first rowsxcols elements */
        break;
    }
    if (r2 > rows) break;
  }

  return (m_dst);
}

/*------------------------------------------------------------*/
float MatrixTrace(MATRIX *M)
{
  int n, nmax;
  float trace;

  if (M->rows >= M->cols)
    nmax = M->cols;
  else
    nmax = M->rows;

  trace = 0.0;
  for (n = 1; n <= nmax; n++) trace += M->rptr[n][n];

  return (trace);
}

/*------------------------------------------------------------------
  MatrixHorCat() - horizontally concatenate matrices m1 and m2, store
  the result in mcat. If mcat is NULL, a new matrix is allocated. A
  pointer to mcat (or the new matrix) is returned.  If m1 is NULL, m2
  is copied into mcat. If m2 is NULL, m1 is copied into mcat.
  -------------------------------------------------------------------*/
MATRIX *MatrixHorCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat)
{
  int r, c1, c2, c;

  if (m1 == NULL && m2 == NULL) {
    printf("ERROR: MatrixHorCat: both m1 and m2 are NULL\n");
    return (NULL);
  }

  if (m1 == NULL) {
    mcat = MatrixCopy(m2, mcat);
    return (mcat);
  }

  if (m2 == NULL) {
    mcat = MatrixCopy(m1, mcat);
    return (mcat);
  }

  if (m1->rows != m2->rows) {
    printf("ERROR: MatrixHorCat: rows of m1 (%d) not equal to m2 (%d)\n", m1->rows, m2->rows);
    return (NULL);
  }

  if (mcat == NULL)
    mcat = MatrixAlloc(m1->rows, m1->cols + m2->cols, MATRIX_REAL);
  else {
    if (mcat->rows != m1->rows) {
      printf("ERROR: MatrixHorCat: rows of m1 (%d) not equal to mcat (%d)\n", m1->rows, mcat->rows);
      return (NULL);
    }
    if (mcat->cols != (m1->cols + m2->cols)) {
      printf("ERROR: MatrixHorCat: cols of mcat (%d) do not equal to m1+m2 (%d)\n", mcat->cols, m1->cols + m2->cols);
      return (NULL);
    }
  }

  /* Fill in the horizontally concatenated matrix */
  for (r = 0; r < mcat->rows; r++) {
    c = 0;
    for (c1 = 0; c1 < m1->cols; c1++) {
      mcat->rptr[r + 1][c + 1] = m1->rptr[r + 1][c1 + 1];
      c++;
    }
    for (c2 = 0; c2 < m2->cols; c2++) {
      mcat->rptr[r + 1][c + 1] = m2->rptr[r + 1][c2 + 1];
      c++;
    }
  }

  return (mcat);
}

/*------------------------------------------------------------------
  MatrixVertCat() - vertically concatenate matrices m1 and m2, store
  the result in mcat. If mcat is NULL, a new matrix is allocated. A
  pointer to mcat (or the new matrix) is returned.  If m1 is NULL, m2
  is copied into mcat. If m2 is NULL, m1 is copied into mcat.
  -------------------------------------------------------------------*/
MATRIX *MatrixVertCat(MATRIX *m1, MATRIX *m2, MATRIX *mcat)
{
  int r, r1, r2, c;

  if (m1 == NULL && m2 == NULL) {
    printf("ERROR: MatrixVertCat: both m1 and m2 are NULL\n");
    return (NULL);
  }

  if (m1 == NULL) {
    mcat = MatrixCopy(m2, mcat);
    return (mcat);
  }

  if (m2 == NULL) {
    mcat = MatrixCopy(m1, mcat);
    return (mcat);
  }

  if (m1->cols != m2->cols) {
    printf("ERROR: MatrixVertCat: cols of m1 (%d) not equal to m2 (%d)\n", m1->cols, m2->cols);
    return (NULL);
  }

  if (mcat == NULL)
    mcat = MatrixAlloc(m1->rows + m2->rows, m1->cols, MATRIX_REAL);
  else {
    if (mcat->cols != m1->cols) {
      printf("ERROR: MatrixVertCat: cols of m1 (%d) not equal to mcat (%d)\n", m1->cols, mcat->cols);
      return (NULL);
    }
    if (mcat->rows != (m1->rows + m2->rows)) {
      printf(
          "ERROR: MatrixVertCat: rows of mcat (%d) "
          "do not equal to m1+m2 (%d)\n",
          mcat->rows,
          m1->rows + m2->rows);
      return (NULL);
    }
  }

  /* Fill in the vertically concatenated matrix */
  for (c = 0; c < mcat->cols; c++) {
    r = 0;
    for (r1 = 0; r1 < m1->rows; r1++) {
      mcat->rptr[r + 1][c + 1] = m1->rptr[r1 + 1][c + 1];
      r++;
    }
    for (r2 = 0; r2 < m2->rows; r2++) {
      mcat->rptr[r + 1][c + 1] = m2->rptr[r2 + 1][c + 1];
      r++;
    }
  }

  return (mcat);
}

/*-------------------------------------------------------------
  MatrixConstVal() - sets all the elements to the given value.
  If X is NULL, then a matrix rows-by-cols is alloced. If X
  is non-NULL, then rows and cols are ignored.
  -------------------------------------------------------------*/
MATRIX *MatrixConstVal(float val, int rows, int cols, MATRIX *X)
{
  int r, c;

  if (X == NULL) X = MatrixAlloc(rows, cols, MATRIX_REAL);

  for (r = 1; r <= X->rows; r++) {
    for (c = 1; c <= X->cols; c++) {
      X->rptr[r][c] = val;
    }
  }

  return (X);
}

/*--------------------------------------------------------------------
  MatrixZero() - sets all the elements to zero.  If X is NULL, then a
  matrix rows-by-cols is alloced. If X is non-NULL, then rows and cols
  are ignored.
  ------------------------------------------------------------------*/
MATRIX *MatrixZero(int rows, int cols, MATRIX *X)
{
  X = MatrixConstVal(0, rows, cols, X);
  return (X);
}

/*-------------------------------------------------------------------*/
/*!
  \fn MATRIX *MatrixSum(MATRIX *m, int dim, MATRIX *msum)
  \brief Computes the sum given matrix
 */
MATRIX *MatrixSum(MATRIX *m, int dim, MATRIX *msum)
{
  int outrows, outcols;
  int r, c;

  if (dim == 1) { /* sum over the rows */
    outrows = 1;
    outcols = m->cols;
  }
  else { /* sum over the cols */
    outrows = m->rows;
    outcols = 1;
  }

  if (msum == NULL)
    msum = MatrixZero(outrows, outcols, NULL);
  else {
    if (msum->rows != outrows || msum->cols != outcols) {
      printf("ERROR: MatrixSum: dimension mismatch\n");
      return (NULL);
    }
    msum = MatrixZero(outrows, outcols, NULL);
  }

  if ((dim == 1 && m->rows > 1) || (dim == 2 && m->cols > 1)) {
    for (r = 1; r <= m->rows; r++) {
      for (c = 1; c <= m->cols; c++) {
        if (dim == 1)
          msum->rptr[1][c] += m->rptr[r][c];
        else
          msum->rptr[r][1] += m->rptr[r][c];
      }
    }
  }
  else { /* Just copy vector to output */
    if (dim == 1)
      for (c = 1; c <= m->cols; c++) msum->rptr[1][c] = m->rptr[1][c];
    else
      for (r = 1; r <= m->rows; r++) msum->rptr[r][1] = m->rptr[r][1];
  }

  return (msum);
}
/*-------------------------------------------------------------------*/
/*!
  \fn MATRIX *MatrixSumSquare(MATRIX *m, int dim, MATRIX *msumsq)
  \brief Computes the sum of the squares of the given matrix
 */
MATRIX *MatrixSumSquare(MATRIX *m, int dim, MATRIX *msumsq)
{
  int outrows, outcols;
  int r, c;

  if (dim == 1) { /* sum over the rows */
    outrows = 1;
    outcols = m->cols;
  }
  else { /* sum over the cols */
    outrows = m->rows;
    outcols = 1;
  }

  if (msumsq == NULL)
    msumsq = MatrixZero(outrows, outcols, NULL);
  else {
    if (msumsq->rows != outrows || msumsq->cols != outcols) {
      printf("ERROR: MatrixSum: dimension mismatch\n");
      return (NULL);
    }
    msumsq = MatrixZero(outrows, outcols, msumsq);
  }

  if ((dim == 1 && m->rows > 1) || (dim == 2 && m->cols > 1)) {
    for (r = 1; r <= m->rows; r++) {
      for (c = 1; c <= m->cols; c++) {
        if (dim == 1)
          msumsq->rptr[1][c] += (m->rptr[r][c] * m->rptr[r][c]);
        else
          msumsq->rptr[r][1] += (m->rptr[r][c] * m->rptr[r][c]);
      }
    }
  }
  else { /* Just copy vector to output */
    if (dim == 1)
      for (c = 1; c <= m->cols; c++) msumsq->rptr[1][c] = (m->rptr[1][c] * m->rptr[1][c]);
    else
      for (r = 1; r <= m->rows; r++) msumsq->rptr[r][1] = (m->rptr[r][1] * m->rptr[r][1]);
  }

  return (msumsq);
}
/*----------------------------------------------------------------*/
double MatrixSumElts(MATRIX *m)
{
  int r, c;
  double msum;

  for (msum = 0.0, r = 1; r <= m->rows; r++)
    for (c = 1; c <= m->cols; c++) msum += m->rptr[r][c];

  return (msum);
}

/*----------------------------------------------------------------*/
double VectorRange(MATRIX *v, double *pVmin, double *pVmax)
{
  double min, max, val;
  int r, c;

  min = v->rptr[1][1];
  max = v->rptr[1][1];
  for (r = 1; r <= v->rows; r++) {
    for (c = 1; c <= v->cols; c++) {
      val = v->rptr[r][c];
      if (min > val) min = val;
      if (max < val) max = val;
    }
  }

  if (pVmin != NULL) *pVmin = min;
  if (pVmax != NULL) *pVmax = max;

  return (max - min);
}

/*----------------------------------------------------------------*/
double VectorSum(const MATRIX *v)
{
  double sum;
  int r, c;

  sum = 0.0;
  for (r = 1; r <= v->rows; r++)
    for (c = 1; c <= v->cols; c++) sum += (v->rptr[r][c]);
  return (sum);
}

/*----------------------------------------------------------------*/
double VectorMean(const MATRIX *v)
{
  double sum, mean;

  sum = VectorSum(v);
  mean = sum / (v->rows * v->cols);

  return (mean);
}

/*----------------------------------------------------------------*/
double VectorVar(MATRIX *v, double *pMean)
{
  double f, mean, sum2, var;
  int r, c;

  mean = VectorMean(v);
  if (pMean != NULL) *pMean = mean;

  sum2 = 0.0;
  for (r = 1; r <= v->rows; r++) {
    for (c = 1; c <= v->cols; c++) {
      f = v->rptr[r][c] - mean;
      sum2 += (f * f);
    }
  }

  var = sum2 / (v->rows * v->cols - 1);

  return (var);
}

/*----------------------------------------------------------------*/
double VectorStdDev(MATRIX *v, double *pMean)
{
  double var;
  var = VectorVar(v, pMean);
  return (sqrt(var));
}

/*----------------------------------------------------------------*/
MATRIX *MatrixDRand48(int rows, int cols, MATRIX *m)
{
  int r, c;

  if (m == NULL) m = MatrixAlloc(rows, cols, MATRIX_REAL);
  /* if m != NULL rows and cols are ignored */

  for (r = 1; r <= m->rows; r++) {
    for (c = 1; c <= m->cols; c++) {
      m->rptr[r][c] = drand48();
    }
  }

  return (m);
}

/*----------------------------------------------------------------*/
MATRIX *MatrixDRand48ZeroMean(int rows, int cols, MATRIX *m)
{
  int r, c;

  if (m == NULL) m = MatrixAlloc(rows, cols, MATRIX_REAL);
  /* if m != NULL rows and cols are ignored */

  for (r = 1; r <= m->rows; r++) {
    for (c = 1; c <= m->cols; c++) {
      m->rptr[r][c] = 2 * drand48() - 1.0;
    }
  }

  return (m);
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
  int r, c;

  if (M->rows != M->cols) {
    printf("ERROR: MatrixFactorSqrSVD: matrix is not square\n");
    return (NULL);
  }

  if (D == NULL)
    D = MatrixAlloc(M->rows, M->cols, MATRIX_REAL);
  else {
    if (D->rows != M->rows || D->cols != M->cols) {
      printf("ERROR: MatrixFactorSqrSVD: dimension mismatch\n");
      return (NULL);
    }
  }

  if (!matalloc) {
    /* Check that M is the same size as last time */
    if (Mrows != M->rows || Mcols != M->cols) {
      MatrixFree(&U);
      MatrixFree(&S);
      MatrixFree(&V);
      MatrixFree(&Vt);
      matalloc = 1;
    }
  }
  if (matalloc) {
    /* Allocate intermediate matrices */
    U = MatrixAlloc(M->rows, M->cols, MATRIX_REAL);
    S = MatrixAlloc(M->rows, 1, MATRIX_REAL);
    V = MatrixAlloc(M->rows, M->cols, MATRIX_REAL);
    Vt = MatrixAlloc(M->rows, M->cols, MATRIX_REAL);
  }

  /* Make a local copy of M because SVD alters it */
  MatrixCopy(M, U);

  /* Compute SVD of matrix [u s v'] = svd(M)*/
  if (MatrixSVD(U, S, V) == NULL) {
    printf("ERROR: MatrixFactorSqrSVD: cound not SVD\n");
    return (NULL);
  }

  /* Check for pos/posdef. Invert if necessary.*/
  for (r = 1; r <= S->rows; r++) {
    s = S->rptr[r][1];
    if (!Invert && s < 0) {
      printf("ERROR: MatrixFactorSqrSVD: matrix is not positive.\n");
      return (NULL);
    }
    if (Invert && s <= 0) {
      printf("ERROR: MatrixFactorSqrSVD: matrix is not positive definite.\n");
      return (NULL);
    }
    if (Invert) s = 1.0 / s;
    S->rptr[r][1] = sqrt(s);
  }

  MatrixTranspose(V, Vt);

  /* D = U*diag(S)*V' */
  /* U = U*diag(S): Multiply each column of U by S. */
  for (c = 1; c <= U->cols; c++)
    for (r = 1; r <= U->rows; r++) U->rptr[r][c] *= S->rptr[c][1];

  /* D = U*Vt = U*diag(S)*V'*/
  MatrixMultiply(U, Vt, D);

  matalloc = 0;
  Mrows = M->rows;
  Mcols = M->cols;
  return (D);
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
  int r, c;

  if (Type != MATRIX_SYM && Type != MATRIX_UPPER && Type != MATRIX_LOWER) {
    printf("ERROR: Type = %d, unrecognized\n", Type);
    return (NULL);
  }

  if (T == NULL) {
    T = MatrixAlloc(v->rows, v->rows, MATRIX_REAL);
    if (T == NULL) return (NULL);
  }
  else {
    if (T->rows != v->rows || T->cols != v->rows) {
      printf("ERROR: dimension mismatch\n");
      return (NULL);
    }
  }

  for (r = 1; r <= T->rows; r++) {
    for (c = 1; c <= T->cols; c++) {
      if (r == c) {
        T->rptr[r][c] = v->rptr[1][1];
        continue;
      }
      if ((r > c && Type == MATRIX_UPPER) || (r < c && Type == MATRIX_LOWER)) {
        T->rptr[r][c] = 0.0;
        continue;
      }
      T->rptr[r][c] = v->rptr[abs(r - c) + 1][1];
    }
  }

  return (T);
}

/*!
\fn MATRIX *MatrixNormalizeColScale(MATRIX *m, MATRIX *scale)
\brief Computes the scaling used for MatrixNormalizeCol()
*/
MATRIX *MatrixNormalizeColScale(MATRIX *m, MATRIX *scale)
{
  int r, c;
  double sum2, v;

  if (scale == NULL) scale = MatrixAlloc(1, m->cols, MATRIX_REAL);
  for (c = 1; c <= m->cols; c++) {
    sum2 = 0.0;
    for (r = 1; r <= m->rows; r++) {
      v = m->rptr[r][c];
      sum2 += (v * v);
    }
    scale->rptr[1][c] = sqrt(sum2);
  }
  // printf("scale -----------------------------\n");
  // MatrixPrint(stdout,scale);
  // printf("------------------------------------\n");
  return (scale);
}

/*!
\fn MATRIX *MatrixNormalizeCol(MATRIX *m, MATRIX *mcnorm)
\brief Rescales m so that sum(col^2)=1
*/
MATRIX *MatrixNormalizeCol(MATRIX *m, MATRIX *mcnorm, MATRIX *scale)
{
  int r, c, FreeScale = 1;
  double v;

  if (mcnorm == NULL) {
    mcnorm = MatrixAlloc(m->rows, m->cols, MATRIX_REAL);
    if (mcnorm == NULL) {
      printf("ERROR: MatrixNormalizeCol: could not alloc\n");
      return (NULL);
    }
  }
  else {
    if (mcnorm->rows != m->rows || mcnorm->cols != m->cols) {
      printf("ERROR: MatrixNormalizeCol: dimension mismatch\n");
      return (NULL);
    }
  }

  if (scale) FreeScale = 0;
  scale = MatrixNormalizeColScale(m, scale);
  for (c = 1; c <= m->cols; c++) {
    v = scale->rptr[1][c];
    if (v != 0)
      for (r = 1; r <= m->rows; r++) mcnorm->rptr[r][c] = (m->rptr[r][c]) / v;
    else
      for (r = 1; r <= m->rows; r++) mcnorm->rptr[r][c] = 0.0;
  }
  if (FreeScale) MatrixFree(&scale);

  // printf("m ----------------------------\n");
  // MatrixPrint(stdout,m);
  // printf("mcnorm ----------------------------\n");
  // MatrixPrint(stdout,mcnorm);

  return (mcnorm);
}

MATRIX *MatrixSimilarityTransform(MATRIX *m_src, MATRIX *m_mul, MATRIX *m_dst)
{
  MATRIX *m_mul_T, *m_tmp;

  m_mul_T = MatrixTranspose(m_mul, NULL);
  m_tmp = MatrixMultiply(m_src, m_mul_T, NULL);
  m_dst = MatrixMultiply(m_mul, m_tmp, m_dst);
  MatrixFree(&m_mul_T);
  MatrixFree(&m_tmp);
  return (m_dst);
}

/*---------------------------------------------------------------
  GaussianVector() - creates a vector with len rows. The gaussian
  curve is centered at the meanth row and has std. The mean can
  be non-integer. If norm == 1, then the sum is adjusted to be 1.
  ---------------------------------------------------------------*/
MATRIX *GaussianVector(int len, float mean, float std, int norm, MATRIX *g)
{
  int n;
  float v, sum, var, f;

  if (g == NULL)
    g = MatrixAlloc(len, 1, MATRIX_REAL);
  else {
    if (g->rows != len) {
      printf("ERROR: GaussianVector: dimension mismatch\n");
      return (NULL);
    }
  }

  var = std * std;
  f = 2 * M_PI * std;
  sum = 0;
  for (n = 0; n < len; n++) {
    v = exp(-(n - mean) * (n - mean) / (2 * var)) / f;
    if (norm) sum += v;
    g->rptr[n + 1][1] = v;
  }

  if (norm) {
    for (n = 0; n < len; n++) {
      g->rptr[n + 1][1] /= sum;
    }
  }

  return (g);
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
  float d, v, sum, var, f;

  if (G == NULL)
    G = MatrixAlloc(len, len, MATRIX_REAL);
  else {
    if (G->rows != len || G->cols != len) {
      printf("ERROR: GaussianMatrix: dimension mismatch\n");
      return (NULL);
    }
  }

  var = std * std;
  f = sqrt(2 * M_PI) * std;

  // Adjust scale so that sum at center line is 1
  if (norm) {
    sum = 0;
    for (c = 0; c < len; c++) {
      d = c - len / 2;
      v = exp(-(d * d) / (2 * var)) / f;
      sum += v;
    }
    f /= sum;
  }

  for (r = 0; r < len; r++) {
    for (c = 0; c < len; c++) {
      d = c - r;
      v = exp(-(d * d) / (2 * var)) / f;
      G->rptr[r + 1][c + 1] = v;
    }
  }

  // printf("std = %g\n",std);
  // MatrixWriteTxt("g.txt",G);
  // exit(1);

  return (G);
}

MATRIX *MatrixSVDPseudoInverse(MATRIX *m, MATRIX *m_pseudo_inv)
{
  int rows = m->rows;
  int cols = m->cols;
  if (rows < cols) {
    MATRIX *mT = MatrixTranspose(m, NULL);
    m_pseudo_inv = MatrixSVDPseudoInverse(mT, m_pseudo_inv);
    MatrixFree(&mT);
    mT = m_pseudo_inv;
    m_pseudo_inv = MatrixTranspose(mT, NULL);
  }
  else {
    int r, c;
    MATRIX *m_U, *m_V, *m_Ur, *m_Vr, *m_Sr, *m_tmp, *m_S;
    VECTOR *v_S;

    m_U = MatrixCopy(m, NULL);
    m_V = MatrixAlloc(cols, cols, MATRIX_REAL);
    v_S = VectorAlloc(cols, MATRIX_REAL);

    if (MatrixIsZero(m)) return (NULL);
    OpenSvdcmp(m_U, v_S, m_V);

    for (r = 1; r <= v_S->rows; r++)
      if (VECTOR_ELT(v_S, r) / VECTOR_ELT(v_S, 1) < 1e-4) break;
    r--;  // previous one was last non-zero
    m_tmp = MatrixCopyRegion(m_U, NULL, 1, 1, m_U->rows, r, 1, 1);
    MatrixFree(&m_U);
    m_Ur = MatrixTranspose(m_tmp, NULL);
    MatrixFree(&m_tmp);
    m_S = MatrixDiag(v_S, NULL);
    VectorFree(&v_S);
    m_Sr = MatrixCopyRegion(m_S, NULL, 1, 1, r, r, 1, 1);
    MatrixFree(&m_S);
    for (c = 1; c <= m_Sr->rows; c++) *MATRIX_RELT(m_Sr, c, c) = 1 / *MATRIX_RELT(m_Sr, c, c);
    m_Vr = MatrixCopyRegion(m_V, NULL, 1, 1, m_V->rows, r, 1, 1);
    MatrixFree(&m_V);
    m_tmp = MatrixMultiply(m_Vr, m_Sr, NULL);
    MatrixFree(&m_Sr);
    MatrixFree(&m_Vr);
    m_pseudo_inv = MatrixMultiply(m_tmp, m_Ur, NULL);
    MatrixFree(&m_tmp);
    MatrixFree(&m_Ur);
  }

  return (m_pseudo_inv);
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
  if (XRO == NULL) XRO = MatrixAlloc(X->rows, X->cols, MATRIX_REAL);
  for (r = 1; r <= X->rows; r++) {
    for (c = 1; c <= X->cols; c++) {
      XRO->rptr[r][c] = X->rptr[NewRowOrder[r - 1]][c];
    }
  }
  return (XRO);
}

/*-----------------------------------------------------------------
  MatrixRandPermRows() - randomly reorders the rows of the input
  matrix.
  -----------------------------------------------------------------*/
int MatrixRandPermRows(MATRIX *X)
{
  int *NewRowOrder, r;
  MATRIX *X0;

  NewRowOrder = RandPerm(X->rows, NULL);
  for (r = 0; r < X->rows; r++) NewRowOrder[r]++;  // Make one-based
  X0 = MatrixCopy(X, NULL);
  MatrixReorderRows(X0, NewRowOrder, X);
  MatrixFree(&X0);
  free(NewRowOrder);
  return (0);
}

/*--------------------------------------------------------------------
  MatrixColsAreNotOrthog() - returns 1 if matrix columns are not
  orthogonal. Computes X'*X and examines the off diagonals. If any
  are greater than 2*FLT_MIN, then returns 1.
  --------------------------------------------------------------------*/
int MatrixColsAreNotOrthog(MATRIX *X)
{
  MATRIX *Xt, *XtX;
  int r, c;

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiply(Xt, X, NULL);
  MatrixFree(&Xt);
  for (r = 1; r <= XtX->rows; r++) {
    for (c = r + 1; c <= XtX->cols; c++) {
      // only upper triangle
      if (XtX->rptr[r][c] > 2 * FLT_MIN) {
        MatrixFree(&XtX);
        return (1);
      }
    }
  }
  MatrixFree(&XtX);
  return (0);
}

double MatrixMahalanobisDistance(VECTOR *v_mean, MATRIX *m_inv_cov, VECTOR *v)
{
  VECTOR *v_dif, *v_tmp;
  double dist;

  v_dif = VectorSubtract(v_mean, v, NULL);
  v_tmp = MatrixMultiply(m_inv_cov, v_dif, NULL);
  dist = VectorDot(v_dif, v_tmp);

  VectorFree(&v_dif);
  VectorFree(&v_tmp);
  return (dist);
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
  MATRIX *drigid = MatrixCopy(m1, NULL);
  if (m2)
    drigid = MatrixSubtract(drigid, m2, drigid);
  else  // subtract identity
  {
    MATRIX *id = MatrixIdentity(4, NULL);
    drigid = MatrixSubtract(drigid, id, drigid);
    MatrixFree(&id);
  }

  // assert we have 4x4 affine transform
  // double EPS = 0.000001;
  // assert(drigid->rows ==4 && drigid->cols == 4);
  // assert(fabs(drigid->rptr[4][1]) < EPS);
  // assert(fabs(drigid->rptr[4][2]) < EPS);
  // assert(fabs(drigid->rptr[4][3]) < EPS);

  // translation norm quadrat:
  double tdq = 0;
  int i;
  for (i = 1; i <= 3; i++) {
    tdq += drigid->rptr[i][4] * drigid->rptr[i][4];
    drigid->rptr[i][4] = 0.0;  // set last row and column to zero
    drigid->rptr[4][i] = 0.0;
  }
  drigid->rptr[4][4] = 0.0;

  MATRIX *dt = MatrixTranspose(drigid, NULL);
  drigid = MatrixMultiply(dt, drigid, drigid);
  MatrixFree(&dt);

  // Trace of A^t A (only first 3x3 submatrix)
  double tr = 0.0;
  for (i = 1; i <= 3; i++) {
    tr += drigid->rptr[i][i];
  }

  MatrixFree(&drigid);

  return sqrt((1.0 / 5.0) * radius * radius * tr + tdq);
}

/* for 3d vector macros */
#include "tritri.h"
int MatrixOrthonormalizeTransform(MATRIX *m_L)
{
  double dot, c1[3], c2[3], c3[3], len;
  int i;

  for (i = 0; i < 3; i++) {
    c1[i] = *MATRIX_RELT(m_L, i + 1, 1);
    c2[i] = *MATRIX_RELT(m_L, i + 1, 2);
    c3[i] = *MATRIX_RELT(m_L, i + 1, 3);
  }

  /* make 1st column vector unit length */
  len = VLEN(c1);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c1, 1.0 / len, c1);

  /* project out component of 2nd vector in direction of 1st column vector */
  dot = DOT(c1, c2);
  for (i = 0; i < 3; i++) c2[i] -= dot * c1[i];

  /* make 2nd column vector unit length */
  len = VLEN(c2);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c2, 1.0 / len, c2);

  /* project out component of 3rd vector in direction of 1st column vector */
  dot = DOT(c1, c3);
  for (i = 0; i < 3; i++) c3[i] -= dot * c1[i];

  /* project out component of 3rd vector in direction of 2nd column vector */
  dot = DOT(c2, c3);
  for (i = 0; i < 3; i++) c3[i] -= dot * c2[i];

  /* make 3rd column vector unit length */
  len = VLEN(c3);
  if (FZERO(len)) len = 1.0f;
  SCALAR_MUL(c3, 1.0 / len, c3);

  for (i = 0; i < 3; i++) {
    *MATRIX_RELT(m_L, i + 1, 1) = c1[i];
    *MATRIX_RELT(m_L, i + 1, 2) = c2[i];
    *MATRIX_RELT(m_L, i + 1, 3) = c3[i];
  }
#if 0
  /* remove translation component */
  for (i = 1 ; i <= 3 ; i++)
    *MATRIX_RELT(m_L, i, 4) = 0 ;
#endif

  return (NO_ERROR);
}

MATRIX *MatrixAsciiReadRaw(const char *fname, MATRIX *m)
{
  FILE *fp;
  int rows, cols, row, col;
  char line[10 * STRLEN], *cp;

  fp = fopen(fname, "r");
  if (fp == NULL) ErrorReturn(NULL, (ERROR_NOFILE, "MatrixAsciiReadRaw: could not open %s\n", fname));

  cp = fgetl(line, 10 * STRLEN - 1, fp);
  for (cols = rows = 0; cp != NULL; rows++) {
    if (rows == 0)  // count cols in first row
    {
      cp = strtok(line, " ");
      for (cols = 0; cp != NULL; cols++) {
        cp = strtok(NULL, " ");
      }
    }
    cp = fgetl(line, 10 * STRLEN - 1, fp);
  }
  m = MatrixAlloc(rows, cols, MATRIX_REAL);

  rewind(fp);
  cp = fgetl(line, 10 * STRLEN - 1, fp);
  for (row = 1; row <= rows; row++) {
    cp = strtok(line, " ");
    for (col = 1; col <= cols; col++) {
      sscanf(cp, "%f", MATRIX_RELT(m, row, col));
      cp = strtok(NULL, " ");
    }
    cp = fgetl(line, 10 * STRLEN - 1, fp);
  }

  fclose(fp);
  return (m);
}

int MatrixToRigidParameters(MATRIX *m, double *pxr, double *pyr, double *pzr, double *pxt, double *pyt, double *pzt)
{
  // M = Mx * My * Mz
  *pxr = atan2(*MATRIX_RELT(m, 2, 3), *MATRIX_RELT(m, 3, 3));
  *pyr = asin(-*MATRIX_RELT(m, 1, 3));
  *pzr = atan2(*MATRIX_RELT(m, 1, 2), *MATRIX_RELT(m, 1, 1));
  *pxt = *MATRIX_RELT(m, 1, 4);
  *pyt = *MATRIX_RELT(m, 2, 4);
  *pzt = *MATRIX_RELT(m, 3, 4);
  return (NO_ERROR);
}
MATRIX *MatrixFromRigidParameters(MATRIX *m, double xr, double yr, double zr, double xt, double yt, double zt)
{
  if (m == NULL) m = MatrixAlloc(4, 4, MATRIX_REAL);

  *MATRIX_RELT(m, 1, 1) = cos(yr) * cos(zr);
  *MATRIX_RELT(m, 1, 2) = cos(yr) * sin(zr);
  *MATRIX_RELT(m, 1, 3) = -sin(yr);

  *MATRIX_RELT(m, 2, 1) = sin(xr) * sin(yr) * cos(zr) - cos(xr) * sin(zr);
  *MATRIX_RELT(m, 2, 2) = sin(xr) * sin(yr) * sin(zr) + cos(xr) * cos(zr);
  *MATRIX_RELT(m, 2, 3) = sin(xr) * cos(yr);

  *MATRIX_RELT(m, 3, 1) = cos(xr) * sin(yr) * cos(zr) + sin(xr) * sin(zr);
  *MATRIX_RELT(m, 3, 2) = cos(xr) * sin(yr) * sin(zr) - sin(xr) * cos(zr);
  *MATRIX_RELT(m, 3, 3) = cos(xr) * cos(yr);

  *MATRIX_RELT(m, 1, 4) = xt;
  *MATRIX_RELT(m, 2, 4) = yt;
  *MATRIX_RELT(m, 3, 4) = zt;
  *MATRIX_RELT(m, 4, 4) = 1.0;

  return (m);
}

int MatrixCheckFinite(MATRIX *m)
{
  int r, c, retval = NO_ERROR;

  for (r = 1; r < m->rows; r++)
    for (c = 1; c < m->cols; c++) {
      if (!std::isfinite(*MATRIX_RELT(m, r, c))) {
        DiagBreak();
        retval = ERROR_BADPARM;
      }
    }
  return (retval);
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

  if (!k) {
    k = MatrixAlloc(rows, cols, MATRIX_REAL);
    if (!k) {
      printf("ERROR: MatrixKron: could not alloc %d %d\n", rows, cols);
      return (NULL);
    }
  }
  else {
    if (k->rows != rows || k->cols != cols) {
      printf("ERROR: MatrixKron: dimension mismatch %d %d vs %d %d\n", rows, cols, k->rows, k->cols);
      return (NULL);
    }
  }

  for (r1 = 1; r1 <= m1->rows; r1++) {
    for (c1 = 1; c1 <= m1->cols; c1++) {
      v1 = m1->rptr[r1][c1];
      for (r2 = 1; r2 <= m2->rows; r2++) {
        for (c2 = 1; c2 <= m2->cols; c2++) {
          v2 = m2->rptr[r2][c2];
          r = (r1 - 1) * m2->rows + r2;
          c = (c1 - 1) * m2->cols + c2;
          // printf("%d %d   %d %d   %d %d\n",r1,c1,r2,c2,r,c);
          k->rptr[r][c] = v1 * v2;
        }
      }
    }
  }
  return (k);
}

/*!
  \fn double MatrixRowDotProduct(MATRIX *m, int row, VECTOR *v)
  \brief dot product of a vector with the row of a matrix
*/
double MatrixRowDotProduct(MATRIX *m, int row, VECTOR *v)
{
  double dot;
  int col;

  for (dot = 0.0, col = 1; col <= m->cols; col++) dot += (*MATRIX_RELT(m, row, col) * VECTOR_ELT(v, col));
  return (dot);
}

/*!
  \fn MATRIX *MatrixDemean(MATRIX *M, MATRIX *Mdm)
  \brief Removes the mean from each column
*/
MATRIX *MatrixDemean(MATRIX *M, MATRIX *Mdm)
{
  int r, c;
  double vsum, vmean;

  if (Mdm == NULL) {
    Mdm = MatrixAlloc(M->rows, M->cols, MATRIX_REAL);
    if (Mdm == NULL) return (NULL);
  }
  if (Mdm->rows != M->rows || Mdm->cols != M->cols) {
    printf("MatrixDemean: dimension mismatch\n");
    return (NULL);
  }
  for (c = 1; c <= M->cols; c++) {
    vsum = 0.0;
    for (r = 1; r <= M->rows; r++) vsum += M->rptr[r][c];
    vmean = vsum / M->rows;
    for (r = 1; r <= M->rows; r++) M->rptr[r][c] -= vmean;
  }
  return (Mdm);
}

/*!
  \fn MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude)
  \brief Creates a new matrix by excluding the given set of rows.
  \param Src - source matrix.
*/
MATRIX *MatrixExcludeFrames(MATRIX *Src, int *ExcludeFrames, int nExclude)
{
  MATRIX *Trg = NULL;
  int q, n, skip, m, c, nframesNew;

  nframesNew = Src->rows - nExclude;

  Trg = MatrixAlloc(nframesNew, Src->cols, MATRIX_REAL);
  q = 0;
  for (n = 0; n < Src->rows; n++) {
    skip = 0;
    for (m = 0; m < nExclude; m++)
      if (n == ExcludeFrames[m]) skip = 1;
    if (skip) continue;
    for (c = 0; c < Src->cols; c++) {
      // printf("%d %d %d %d\n",n,q,c,Src->cols);
      Trg->rptr[q + 1][c + 1] = Src->rptr[n + 1][c + 1];
    }
    q++;
  }
  return (Trg);
}
/*!
  \fn int MatrixIsIdentity(MATRIX *m)
  \brief returns 1 if matrix m is the identity matrix, 0 otherwise
*/
int MatrixIsIdentity(MATRIX *m)
{
  int r, c;

  for (r = 1; r < m->rows; r++)
    for (c = 1; c < m->cols; c++) {
      if (r == c) {
        if (FEQUAL(*MATRIX_RELT(m, r, c), 1) == 0) return (0);
      }
      else if (FEQUAL(*MATRIX_RELT(m, r, c), 0.0) == 0)
        return (0);
    }
  return (1);
}

/*!
  \fn MATRIX *MatrixCumTrapZ(MATRIX *y, MATRIX *t, MATRIX *yz)
  \brief Computes trapezoidal integration (like matlab cumtrapz)
*/
MATRIX *MatrixCumTrapZ(MATRIX *y, MATRIX *t, MATRIX *yz)
{
  if (yz == NULL) yz = MatrixAlloc(y->rows, y->cols, MATRIX_REAL);

  int c, f;
  double v, vprev, vsum, dt;

  for (c = 0; c < y->cols; c++) {
    yz->rptr[1][c + 1] = 0;
    vsum = 0;
    vprev = y->rptr[1][c + 1];
    for (f = 1; f < y->rows; f++) {
      dt = t->rptr[f + 1][1] - t->rptr[f][1];
      v = y->rptr[f + 1][c + 1];
      vsum += (dt * ((v + vprev) / 2));
      yz->rptr[f + 1][c + 1] = vsum;
      vprev = v;
    }
  }
  return (yz);
}
//---------------------------------------------------------
/*!
  \fn MATRIX *ANOVAContrast(int *FLevels, int nFactors, int *FactorList, int nFactorList)
  \brief Computes a contrast matrix to test an effect or interaction in an ANOVA.
   FLevels is an array of length nFactors with the number of levels for each factor.
     Number of levels must be >= 2 or else it is not a factor.
   All factors must be discrete factors (ie, not continuous)
   FactorList is an array of length nFactorList with the Factors to test.
   Eg, FLevels = [2 2], FactorList = [1] tests for the main effect of Factor 1.
   Eg, FLevels = [2 3], FactorList = [1 2] tests for the interaction between Factors 1 and 2
   The Factors in the FactorList should be 1-based.
   The regressors should have the following order
      (eg, if factors are gender, handedness, and diagnosis)
   F1L1-F2L1-F3L1    M-L-N
   F1L1-F2L1-F3L2    M-L-AD
   F1L1-F2L2-F3L1    M-R-N
   F1L1-F2L2-F3L2    M-R-AD
   F1L2-F2L1-F3L1    F-L-N
   F1L2-F2L1-F3L2    F-L-AD
   F1L2-F2L2-F3L1    F-R-N
   F1L2-F2L2-F3L2    F-R-AD

   This should be general enough to handle any number of Factors with any number
   of levels per factor. The number of levels does not need to be the same across
   factors.

   Woodward, J. A., Bonett, D. G., & Brecht, M-L. (1990). Introduction
   to linear models and experimental design. San Diego, CA: Harcourt
   Brace Jovanovich.

*/
MATRIX *ANOVAContrast(int *FLevels, int nFactors, int *FactorList, int nFactorList)
{
  int n, nthFactor, InList;
  MATRIX *M, *Mf, *K;

  for (nthFactor = 0; nthFactor < nFactors; nthFactor++) {
    if (FLevels[nthFactor] < 2) {
      printf("ERROR: ANOVAContrast: Factor %d has only %d levels.\n", nthFactor, FactorList[nthFactor]);
      printf("       Must have at least 2 levels\n");
      return (NULL);
    }
  }
  for (nthFactor = 0; nthFactor < nFactorList; nthFactor++) {
    if (FactorList[nthFactor] > nFactors) {
      printf("ERROR: ANOVAContrast: %d > %d\n", FactorList[nthFactor], nFactors);
      return (NULL);
    }
    if (FactorList[nthFactor] < 1) {
      printf("ERROR: ANOVAContrast: Factor %d = %d < 1\n", nthFactor, FactorList[nthFactor]);
      return (NULL);
    }
  }

  M = MatrixConstVal(1, 1, 1, NULL);
  for (nthFactor = 0; nthFactor < nFactors; nthFactor++) {
    InList = 0;
    for (n = 0; n < nFactorList; n++) {
      if (nthFactor + 1 == FactorList[n]) {
        InList = 1;
        break;
      }
    }
    if (InList)
      Mf = ANOVAOmnibus(FLevels[nthFactor]);
    else
      Mf = ANOVASummingVector(FLevels[nthFactor]);
    // printf("Mf %d --------------------\n",InList);
    // MatrixPrint(stdout,Mf);
    K = MatrixKron(M, Mf, NULL);
    MatrixFree(&Mf);
    MatrixFree(&M);
    M = K;
  }
  return (M);
}
//---------------------------------------------------------
/*!
  \fn MATRIX *ANOVASummingVector(int nLevels)
  \brief Summing vector for main effect of a factor. nLevels = number
  of levels in the factor.
*/
MATRIX *ANOVASummingVector(int nLevels)
{
  MATRIX *S;
  S = MatrixConstVal(1, 1, nLevels, NULL);
  return (S);
}
//---------------------------------------------------------
/*!
  \fn MATRIX *ANOVASelectionVector(int nLevels, int Level)
  \brief Selects a given level of a given factor. nLevels = number of
  levels in the factor, Level is the one-based level number to select.
  Used to create a Simple Main Effect ANOVA matrix.
*/
MATRIX *ANOVASelectionVector(int nLevels, int Level)
{
  MATRIX *S;
  S = MatrixConstVal(0, 1, nLevels, NULL);
  S->rptr[1][Level] = 1;
  return (S);
}
//---------------------------------------------------------
MATRIX *ANOVAOmnibus(int nLevels)
{
  MATRIX *O;
  int r;
  O = MatrixAlloc(nLevels - 1, nLevels, MATRIX_REAL);
  for (r = 0; r < nLevels - 1; r++) {
    O->rptr[r + 1][r + 1] = 1;
    O->rptr[r + 1][nLevels] = -1;
  }
  return (O);
}
double MatrixSSE(MATRIX *m1, MATRIX *m2)
{
  int r, c;
  double sse, error;

  for (sse = 0.0, r = 1; r <= m1->rows; r++)
    for (c = 1; c <= m1->cols; c++) {
      error = *MATRIX_RELT(m1, r, c) - *MATRIX_RELT(m2, r, c);
      sse += SQR(error);
    }
  return (sse);
}

double MatrixRMS(MATRIX *m1, MATRIX *m2)
{
  double sse;

  sse = MatrixSSE(m1, m2);
  return (sqrt(sse / (m1->rows * m1->cols)));
}
/*
  \fn MATRIX *MatrixMtM(MATRIX *m, MATRIX *mout)
  \brief Efficiently computes M'*M. There are several optimizations:
  (1) exploits symmetry, (2) exploits sparsity, and (3) creates a LUT
  (which requires some resident storage). The LUT allows for better
  load balancing when parallelizing with Open MP (makes it almost 2x
  faster).
 */
MATRIX *MatrixMtM(MATRIX *m, MATRIX *mout)
{
  int c1, c2, n, rows, cols;
  static int *c1list = NULL, *c2list = NULL, ntot = 0;

  if (mout == NULL) mout = MatrixAlloc(m->cols, m->cols, MATRIX_REAL);
  if (mout->rows != m->cols) {
    printf("ERROR: MatrixMtM() mout cols (%d) != m cols (%d)\n", mout->cols, m->cols);
    return (NULL);
  }
  if (mout->cols != m->cols) {
    printf("ERROR: MatrixMtM() mout cols (%d) != m cols (%d)\n", mout->cols, m->cols);
    return (NULL);
  }

  rows = m->rows;
  cols = m->cols;

  if (ntot == 0) {
    // create a lookup table that maps n to c1 and c2
    if (ntot != ((cols * cols) + cols) / 2) {
      if (c1list) free(c1list);
      if (c2list) free(c2list);
    }
    ntot = ((cols * cols) + cols) / 2;
    if (Gdiag_no > 0) printf("MatrixMtM: Alloc rows=%d cols=%d ntot=%d\n", rows, cols, ntot);
    c1list = (int *)calloc(sizeof(int), ntot);
    c2list = (int *)calloc(sizeof(int), ntot);
    n = 0;
    for (c1 = 1; c1 <= cols; c1++) {
      for (c2 = c1; c2 <= cols; c2++) {
        c1list[n] = c1;
        c2list[n] = c2;
        n++;
      }
    }
  }

/* Loop over the number of distinct elements in the symetric matrix. Using
   the LUT created above is better for load balancing.  */
#ifdef HAVE_OPENMP
  #pragma omp parallel for shared(c1list, c2list, ntot, cols, rows, mout, m)
#endif
  for (n = 0; n < ntot; n++) {
    
    double v, v1, v2;
    int c1, c2, r;
    c1 = c1list[n];
    c2 = c2list[n];
    v = 0;
    for (r = 1; r <= rows; r++) {
      v1 = m->rptr[r][c1];
      if (v1 == 0) continue;
      v2 = m->rptr[r][c2];
      if (v2 == 0) continue;
      v += v1 * v2;
      continue;
    }  // row
    mout->rptr[c1][c2] = v;
    mout->rptr[c2][c1] = v;
  }

  if (0) {
    // This is a built-in test. The difference should be 0.
    MATRIX *mt, *mout2;
    double d, dmax;
    int c, r;
    mt = MatrixTranspose(m, NULL);
    mout2 = MatrixMultiplyD(mt, m, NULL);
    dmax = 0;
    for (r = 1; r <= mout->rows; r++) {
      for (c = 1; c <= mout->cols; c++) {
        d = mout->rptr[r][c] - mout2->rptr[r][c];
        if (dmax < fabs(d)) {
          dmax = fabs(d);
          if (dmax > .01) printf("%5d %5d %lf %g %g\n", r, c, d, mout->rptr[r][c], mout2->rptr[r][c]);
        }
      }
    }
    printf("MatrixMtM: test MAR %g\n", dmax);
    MatrixFree(&mt);
    MatrixFree(&mout2);
  }

  return (mout);
}

/*
  \fn MATRIX *MatrixSkew(MATRIX *y, MATRIX *s)
  \brief Computes the skew (3rd moment) of the values in each column
  of y. If y has multiple columns, a skew is computed for each one.xx
*/
MATRIX *MatrixSkew(MATRIX *y, MATRIX *s)
{
  int c, r;
  double mn, m2, m3, g1, delta, n, adj;

  if (s == NULL) s = MatrixAlloc(1, y->cols, MATRIX_REAL);

  n = y->rows;
  adj = sqrt(n * (n - 1)) / (n - 2);

  for (c = 0; c < y->cols; c++) {
    mn = 0;
    for (r = 0; r < y->rows; r++) mn += y->rptr[r + 1][c + 1];
    mn /= y->rows;
    m2 = 0;
    m3 = 0;
    for (r = 0; r < y->rows; r++) {
      delta = y->rptr[r + 1][c + 1] - mn;
      m2 += pow(delta, 2.0);  // sum of squares
      m3 += pow(delta, 3.0);  // sum of cubes
    }
    m2 /= n;
    m3 /= n;
    if (m2 != 0)
      g1 = adj * m3 / pow(m2, 1.5);
    else
      g1 = 0;
    // printf("skew: %2d mn=%g, m2=%g, m3=%g g1=%g\n",c,mn,m2,m3,g1);
    s->rptr[1][c + 1] = g1;
  }
  return (s);
}
/*
  \fn MATRIX *MatrixKurtosis(MATRIX *y, MATRIX *k)
  \brief Computes the 'unbiased' kurtosis (4th moment) of the values
  in each column of y. If y has multiple columns, a kurtosis is
  computed for each one.x
*/
MATRIX *MatrixKurtosis(MATRIX *y, MATRIX *k)
{
  int c, r;
  double mn, m4 = 0, m2 = 0, g2, delta, b1, b2, n;

  if (k == NULL) k = MatrixAlloc(1, y->cols, MATRIX_REAL);

  n = y->rows;
  b1 = (n + 1) * (n - 1) / ((n - 2) * (n - 3));
  b2 = ((n - 1) * (n - 1)) / ((n - 2) * (n - 3));
  // printf("kurt: n=%d, b1=%g b2=%g\n",(int)n,b1,b2);

  for (c = 0; c < y->cols; c++) {
    mn = 0;
    for (r = 0; r < y->rows; r++) mn += y->rptr[r + 1][c + 1];
    mn /= y->rows;
    m2 = 0;
    m4 = 0;
    for (r = 0; r < y->rows; r++) {
      delta = y->rptr[r + 1][c + 1] - mn;
      m2 += pow(delta, 2.0);  // sum of squares
      m4 += pow(delta, 4.0);  // sum of quads
    }
    m4 *= y->rows;
    // Formula below usually has a +3, but this is left off so that k has 0 mean
    if (m2 != 0)
      g2 = b1 * (m4 / (m2 * m2)) - 3 * b2;
    else
      g2 = 0;
    // printf("kurt: %2d mn=%g, m2=%g, m4=%g, g2=%g\n",c,mn,m2,m4,g2);
    k->rptr[1][c + 1] = g2;
  }
  return (k);
}

/*
  \fn MATRIX *MatrixAtB(MATRIX *A, MATRIX *B, MATRIX *mout)
  \brief Computes A'*B without computing or allocating A'
  explicitly. This can be helpful whan A is a large matrix.
  Accumlates using double. OpenMP capable.
 */
MATRIX *MatrixAtB(MATRIX *A, MATRIX *B, MATRIX *mout)
{
  int colA;

  if (A->rows != B->rows) {
    printf("ERROR: MatrixAtB(): dim mismatch: %d %d\n", A->rows, B->rows);
    return (NULL);
  }
  if (mout == NULL) {
    mout = MatrixAlloc(A->cols, B->cols, MATRIX_REAL);
    if (mout == NULL) {
      printf("ERROR: MatrixAtB(): could not alloc %d %d\n", A->cols, B->cols);
      return (NULL);
    }
  }

#ifdef HAVE_OPENMP
  #pragma omp parallel for 
#endif
  for (colA = 0; colA < A->cols; colA++) {
    int row, colB;
    double sum;
    for (colB = 0; colB < B->cols; colB++) {
      sum = 0;
      for (row = 0; row < A->rows; row++) sum += (double)A->rptr[row + 1][colA + 1] * B->rptr[row + 1][colB + 1];
      mout->rptr[colA + 1][colB + 1] = sum;
    }
  }
  
  if (0) {
    // In this test, dmax should be 0 because MatrixMultiplyD() is used
    MATRIX *At = MatrixTranspose(A, NULL);
    MATRIX *mout2 = MatrixMultiplyD(At, B, NULL);
    int r, c;
    double d, dmax;
    dmax = 0;
    for (r = 1; r <= mout->rows; r++) {
      for (c = 1; c <= mout->cols; c++) {
        d = mout->rptr[r][c] - mout2->rptr[r][c];
        if (dmax < fabs(d)) {
          dmax = fabs(d);
          if (dmax > .01) printf("%5d %5d %lf %g %g\n", r, c, d, mout->rptr[r][c], mout2->rptr[r][c]);
        }
      }
    }
    printf("MatrixAtB: test MAR %g\n", dmax);
    MatrixFree(&At);
    MatrixFree(&mout2);
  }

  return (mout);
}
/*
  \fn double MatrixMaxAbsDiff(MATRIX *m1, MATRIX *m2, double dthresh)
  \brief Computes the maximum element by element absolute difference
  between two matrices. If dthresh is > 0, then prints more info about
  each element whose diff is > dthresh
*/
double MatrixMaxAbsDiff(MATRIX *m1, MATRIX *m2, double dthresh)
{
  int r, c;
  double d, dmax;
  if (m1->rows != m2->rows || m1->cols != m2->cols) {
    printf("ERROR: MatrixMaxAbsDiff(): dim mismatch\n");
    return (-1);
  }
  dmax = 0;
  for (r = 1; r <= m1->rows; r++) {
    for (c = 1; c <= m1->cols; c++) {
      d = (double)m1->rptr[r][c] - m2->rptr[r][c];
      if (dthresh > 0 && fabs(d) > dthresh) printf("%5d %5d %lf %g %g\n", r, c, d, m1->rptr[r][c], m2->rptr[r][c]);
      if (dmax < fabs(d)) dmax = fabs(d);
    }
  }
  return (dmax);
}

/*
  \fn MATRIX *MatrixColNullSpace(MATRIX *M, int *err)
  \brief Computes the column null space of the given matrix using
  SVD. By definition, N'*M = 0 (or close to 0). If there is no null
  space then returns a null matrix. The number of rows of M must be >=
  number of columns or else a null matrix is returned and
  err=1. Otherwise err=0.
 */
MATRIX *MatrixColNullSpace(MATRIX *M, int *err)
{
  MATRIX *u, *s, *v, *M2;
  int dim, r, c, dimmax, ndim;
  double thresh = 0;
  MATRIX *N;

  *err = 0;

  if (M->rows < M->cols) {
    // Not sure why this does not work, but all singular values come out 0
    printf("ERROR: MatrixColNullSpace(): rows (%d) must be >= cols (%d)\n", M->rows, M->cols);
    *err = 1;  // can't just return a null matrix
    return (NULL);
  }

  dimmax = MAX(M->rows, M->cols);

  if (M->rows > M->cols) {
    // Create a square matrix by padding extra colums with 0
    M2 = MatrixAlloc(M->rows, M->rows, MATRIX_REAL);
    for (r = 1; r <= M->rows; r++) {
      for (c = 1; c <= M->cols; c++) {
        M2->rptr[r][c] = M->rptr[r][c];
      }
    }
  }
  else
    M2 = MatrixCopy(M, NULL);

  // Compute SVD M2 = u*s*v'
  u = MatrixCopy(M2, NULL);  // It's done in-place so make a copy
  s = RVectorAlloc(M2->cols, MATRIX_REAL);
  v = MatrixAlloc(M2->cols, M2->cols, MATRIX_REAL);
  OpenSvdcmp(u, s, v);

  // Determine dimension
  if (fabs(s->rptr[1][1]) > .00000001) {
    // printf("s1 = %e\n",s->rptr[1][1]);
    thresh = dimmax * s->rptr[1][1] * .0000001;  // not sure
    for (dim = 1; dim <= s->cols; dim++)
      if (s->rptr[1][dim] < thresh) break;
    dim--;
  }
  else
    dim = 0;

  ndim = M->rows - dim;  //  dim of the null space
  if (Gdiag_no > 0) printf("MatrixNullSpace(): dim = %d, ndim = %d, %e\n", dim, ndim, thresh);

  // Definition of null space: N'*M = 0
  N = NULL;
  if (dim != dimmax) {
    N = MatrixAlloc(M->rows, ndim, MATRIX_REAL);
    for (r = 1; r <= M->rows; r++) {
      for (c = 1; c <= ndim; c++) {
        N->rptr[r][c] = u->rptr[r][c + dim];
      }
    }
  }

  if (0) {
    // Test that N'*M is close to 0
    MATRIX *Nt, *P;
    double pmax = 0;
    Nt = MatrixTranspose(N, NULL);
    P = MatrixMultiplyD(Nt, M, NULL);
    for (r = 1; r <= Nt->rows; r++) {
      for (c = 1; c <= M->cols; c++) {
        if (fabs(P->rptr[r][c]) > pmax) pmax = fabs(P->rptr[r][c]);
        if (fabs(P->rptr[r][c]) > 10e-6) printf("TEST: MatrixNullSpace(): %d %d %e\n", r, c, P->rptr[r][c]);
      }
    }
    printf("TEST: MatrixNullSpace(): pmax %le\n", pmax);
    MatrixFree(&Nt);
    MatrixFree(&P);
  }

  // MatrixWriteTxt("m.mtx",M2);
  // MatrixWriteTxt("u.mtx",u);
  // MatrixWriteTxt("s.mtx",s);
  // MatrixWriteTxt("v.mtx",v);
  // MatrixWriteTxt("n.mtx",N);

  MatrixFree(&u);
  MatrixFree(&s);
  MatrixFree(&v);
  MatrixFree(&M2);

  return (N);
}

/*
  \fn MATRIX *MatrixResidualForming(MATRIX *X, MATRIX *R)
  \brief Computes the residual forming matrix
    R = I - X*inv(X'*X)*X';
 */
MATRIX *MatrixResidualForming(MATRIX *X, MATRIX *R)
{
  MATRIX *Xt, *XtX, *iXtX, *I, *XiXtX, *XiXtXXt;

  Xt = MatrixTranspose(X, NULL);
  XtX = MatrixMultiplyD(Xt, X, NULL);
  iXtX = MatrixInverse(XtX, NULL);
  if (iXtX == NULL) {
    printf("ERROR: MatrixResidualForming(): X is not invertable\n");
    MatrixFree(&Xt);
    MatrixFree(&XtX);
    return (NULL);
  }
  I = MatrixIdentity(X->rows, NULL);
  XiXtX = MatrixMultiplyD(X, iXtX, NULL);
  XiXtXXt = MatrixMultiplyD(XiXtX, Xt, NULL);
  R = MatrixSubtract(I, XiXtXXt, R);

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&I);
  MatrixFree(&XiXtX);
  MatrixFree(&XiXtXXt);

  return (R);
}

/*!
  MATRIX *MatrixMultiplyElts(MATRIX *m1, MATRIX *m2, MATRIX *m12)
  Multiply each element in each matrix (same as m1.*m2 in matlab)
 */
MATRIX *MatrixMultiplyElts(MATRIX *m1, MATRIX *m2, MATRIX *m12)
{
  int c,r;
  if(m1->rows != m2->rows){
    printf("ERROR: MatrixMultiplyElts(): rows not equal %d %d\n",m1->rows,m2->rows);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }
  if(m1->cols != m2->cols){
    printf("ERROR: MatrixMultiplyElts(): cols not equal %d %d\n",m1->cols,m2->cols);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }
  if(m12 == NULL) 
    m12 = MatrixAlloc(m1->rows,m1->cols,MATRIX_REAL);
  if(m12->rows != m2->rows){
    printf("ERROR: MatrixMultiplyElts(): m12 rows not equal %d %d\n",m12->rows,m2->rows);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }
  if(m12->cols != m2->cols){
    printf("ERROR: MatrixMultiplyElts(): m12 cols not equal %d %d\n",m12->cols,m2->cols);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }

  for(c=0; c < m1->cols; c++){
    for(r=0; r < m1->rows; r++){
      m12->rptr[r+1][c+1] = (double)m1->rptr[r+1][c+1] * (double)m2->rptr[r+1][c+1];
    }
  }

  return(m12);
}

/*!
  \fn MATRIX *MatrixElementDivide(MATRIX *num, MATRIX *den, MATRIX *quotient)
  \brief Element-wise matrix division. q = n/(d+FLT_EPSILON). Only works on MATRIX_REAL.
  \parameter num - numerator
  \parameter den - denominator
  \parameter quotient - result
*/
MATRIX *MatrixDivideElts(MATRIX *num, MATRIX *den, MATRIX *quotient)
{
  int r,c;

  if(num->rows != den->rows || num->cols != den->cols){
    printf("ERROR: MatrixtDivideElts(): dim mismatch\n");
    printf("%s:%d\n",__FILE__,__LINE__);
    return(NULL);
  }
  if(quotient==NULL)
    quotient = MatrixAlloc(num->rows,num->cols,MATRIX_REAL);

  for(r=0; r < num->rows; r++){
    for(c=0; c < num->cols; c++){
      quotient->rptr[r+1][c+1] = num->rptr[r+1][c+1]/(den->rptr[r+1][c+1] + FLT_EPSILON);
    }
  }
  return(quotient);
}

/*!
  MATRIX *MatrixReplicate(MATRIX *mIn, int nr, int nc, MATRIX *mOut)
  Replicate the input matrix nr times in the row direction and nc times 
  in the col direction (same as repmat(mIn,[nr nc]) in matlab)
 */
MATRIX *MatrixReplicate(MATRIX *mIn, int nr, int nc, MATRIX *mOut)
{
  int c,r,cc,rr,cout,rout;
  if(mOut == NULL) mOut = MatrixAlloc(nr*mIn->rows,nc*mIn->cols,MATRIX_REAL);    
  if(mOut->rows != nr*mIn->rows){
    printf("ERROR: MatrixReplicate(): mOut rows not equal %d %d\n",mOut->rows,nr*mIn->rows);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }
  if(mOut->cols != nc*mIn->cols){
    printf("ERROR: MatrixReplicate(): mOut cols not equal %d %d\n",mOut->cols,nc*mIn->cols);
    printf("  break %s:%d\n", __FILE__, __LINE__);
    return(NULL);
  }

  for(cc=0; cc < nc; cc++){
    for(rr=0; rr < nr; rr++){{
	for(c=0; c < mIn->cols; c++){
	  for(r=0; r < mIn->rows; r++){
	    cout = cc*mIn->cols + c;
	    rout = rr*mIn->rows + r;
	    mOut->rptr[rout+1][cout+1] = mIn->rptr[r+1][c+1];
	  }
	}
      }
    }
  }

  return(mOut);
}

/*!
  \fn MATRIX *MatrixGlmFit(MATRIX *y, MATRIX *X, double *pRVar, MATRIX *beta)
  \brief Solves the GLM
*/
MATRIX *MatrixGlmFit(MATRIX *y, MATRIX *X, double *pRVar, MATRIX *beta)
{
  MATRIX *Xt, *XtX, *iXtX, *Xty, *yhat, *res;
  double mres,rvar;

  Xt   = MatrixTranspose(X,NULL);
  XtX  = MatrixMultiplyD(Xt,X,NULL);
  iXtX = MatrixInverse(XtX,NULL);
  Xty  = MatrixMultiplyD(Xt,y,NULL);
  beta = MatrixMultiplyD(iXtX,Xty,beta);
  yhat = MatrixMultiplyD(X,beta,NULL);
  res  = MatrixSubtract(y, yhat, NULL);
  rvar = VectorVar(res,&mres);
  rvar = rvar*(X->rows-1)/(X->rows-X->cols);
  *pRVar = rvar;

  MatrixFree(&Xt);
  MatrixFree(&XtX);
  MatrixFree(&iXtX);
  MatrixFree(&Xty);
  MatrixFree(&yhat);
  MatrixFree(&res);
  return(beta);
}

/*!
\fn MATRIX *MatrixACF2Kernel(MATRIX *acf, MATRIX *k)
\brief Computes a FIR kernel from an autocorrelation function (ACF)
that, when applied to data, will induce the given ACF. The ACF should
be an Nx1 matrix where acf[1]=1, acf[2]=fistlag, etc. The returned Nx1
kernel can then be convolved with the data to realize the ACF. The
method coputes the SVD of the toeplitzized ACF, then F =
U*sqrt(S2)*U', and the kernel is extracted from F. sum(kernel) = 1;
*/
MATRIX *MatrixACF2Kernel(MATRIX *acf, MATRIX *k, MATRIX *F)
{
  // Error checking
  if(acf == NULL){
    printf("ERROR: MatrixACF2Kernel() acf is NULL\n");
    return(NULL);
  }
  if(acf->rows < 2){
    printf("ERROR: MatrixACF2Kernel() acf rows < 2\n");
    return(NULL);
  }
  if(k==NULL)  k = MatrixAlloc(acf->rows,1,MATRIX_REAL);
  if(acf->rows != k->rows){
    printf("ERROR: MatrixACF2Kernel() acf rows = %d, krows = %d\n",acf->rows,k->rows);
    return(NULL);
  }
  if(k->cols != 1){
    printf("ERROR: MatrixACF2Kernel() kernel cols = %d, shoudl be 1\n",k->cols);
    return(NULL);
  }
  if(F){
    if(acf->rows != F->rows){
      printf("ERROR: MatrixACF2Kernel() acf rows = %d, Frows = %d\n",acf->rows,F->rows);
      return(NULL);
    }
    if(acf->rows != F->cols){
      printf("ERROR: MatrixACF2Kernel() acf rows = %d, Fcols = %d\n",acf->rows,F->cols);
      return(NULL);
    }
  }

  // Create a symmetric toeplitz matrix from the ACF
  MATRIX *T = MatrixToeplitz(acf, NULL, MATRIX_SYM);

  // Compute the SVD of the toeplitz
  MATRIX *S2 = RVectorAlloc(acf->rows, MATRIX_REAL);
  MATRIX *U  = MatrixCopy(T, NULL);      // It's done in-place so make a copy
  MATRIX *V  = MatrixSVD(U, S2, NULL);  // T = U*S2*V'; U=V because sym
  if(V==NULL){
    MatrixFree(&T);
    MatrixFree(&U);
    MatrixFree(&S2);
    printf("ERROR: MatrixACF2Kernel() SVD failed\n");
    return(NULL);
  }

  // Below computes the kernel based on F=U*sqrt(S2)*U' where F is the
  // NxN filter that will impose the ACF on the data. The full matrix
  // implementation (below) is clearer, but this implementation is
  // much faster
  int nmid = round(acf->rows/2.0);
  int m,nhits=0;
  double ksum=0;
  for(m=1; m <= acf->rows; m++){
    if(S2->rptr[1][m] <= 0) continue;
    nhits ++;
    double s = sqrt(S2->rptr[1][m]);
    double v0 = U->rptr[nmid][m];
    for(int n=nmid; n <= acf->rows; n++){
      double v = U->rptr[n][m] * v0 * s;
      k->rptr[n-nmid+1][1] += v;
      ksum += v;
    }
  }
  if(nhits == 0){
    MatrixFree(&T);
    MatrixFree(&U);
    MatrixFree(&S2);
    MatrixFree(&V);
    printf("ERROR: MatrixACF2Kernel() no values > 0\n");
    return(NULL);
  }
  // Normalize so that sum(k)=1; calling app may have to do this again, eg,
  // if using in 2D or 3D context.
  for(m=1; m <= acf->rows; m++) k->rptr[m][1] /= ksum;

  //printf("acf = [");  MatrixPrint(stdout,acf); printf("];\n");
  //printf("T = [");  MatrixPrint(stdout,T); printf("];\n");
  //printf("U = [");  MatrixPrint(stdout,U); printf("];\n");
  //printf("S2 = [");  MatrixPrint(stdout,S2); printf("];\n");
  //printf("V = [");  MatrixPrint(stdout,V); printf("];\n");
  //printf("k = [");  MatrixPrint(stdout,k); printf("];\n");

  if(F){
    // This is an alternative way to compute the kernel using full
    // matrix operations. It is clearer than above but more
    // computationally expensive
    for(m=0; m < acf->rows; m++){
      if(S2->rptr[1][m+1] <= 0) S2->rptr[1][m+1]=0;
      else S2->rptr[1][m+1] = sqrt(S2->rptr[1][m+1]);
    }
    MATRIX *S = MatrixDiag(S2, NULL);
    MATRIX *Vt = MatrixTranspose(V, NULL);
    
    // This is a filter matrix which could be applied to 1D data 
    // F = U*S*V' = U*S*U';
    F = MatrixMultiplyD(U,S,F);
    F = MatrixMultiplyD(F,Vt,F);
    //printf("F = [");  MatrixPrint(stdout,F); printf("];\n");

    #if 0
    // Extract kernel from F. Should yield same result as above
    int n;
    ksum=0;
    for(n=nmid; n < acf->rows; n++){
      k->rptr[n-nmid+1][1] = F->rptr[n][nmid];
      ksum += k->rptr[n-nmid+1][1];
    }
    for(n=nmid; n < acf->rows; n++) k->rptr[n-nmid+1][1] /= ksum;
    printf("k2 = [");  MatrixPrint(stdout,k); printf("];\n");
    MatrixFree(&k2);
    #endif
  }

  MatrixFree(&T);
  MatrixFree(&U);
  MatrixFree(&S2);
  MatrixFree(&V);

  return(k);
}
