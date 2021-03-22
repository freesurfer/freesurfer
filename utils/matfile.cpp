/**
 * @brief read matlab files
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "diag.h"
#include "error.h"
#include "matrix.h"
#include "mghendian.h"
#include "proto.h"
#include "zlib.h"

#include "matfile.h"

static double **matAlloc(int rows, int ncols);
static void matFree(double **matrix, int nrows, int ncols);
static int readMatFile(FILE *fp, MATFILE *mf, double **real_matrix, double **imag_matrix);
static int znzreadMatFile(char *unbuff, MATFILE *mf, double **real_matrix, double **imag_matrix);
static void swapBytes(MATFILE *mf);
static const char *MatProgname = "matfile";

int MatFileWrite(const char *fname, float *data, int rows, int cols, const char *name);
char *MatReadHeader(FILE *fp, MATFILE *mf, long32 *compressed);
char *znzMatReadHeader(FILE *fp, MATFILE *mf, char **data);
MLFC *ReadMatlabFileContents(const char *fname);
int MLFCfree(MLFC **ppmlfc);

#if (BYTE_ORDER == LITTLE_ENDIAN)
#define DIFFERENT_ENDIAN(mf) ((mf->version == 5) ? (mf->endian != MATFILE_PC5) : (mf->type != MATFILE_PC))
#elif (BYTE_ORDER == BIG_ENDIAN)
#define DIFFERENT_ENDIAN(mf) ((mf->version == 5) ? (mf->endian == MATFILE_PC5) : (mf->type == MATFILE_PC))
#else
#error 'BYTE_ORDER not set'
#endif

static int (*mat_printf)(const char *szFormat, ...) = NULL;
int Matlab_Install_printf(int (*new_printf)(const char *szFormat, ...))
{
  mat_printf = new_printf;
  return (0);
}

int MatlabWrite(MATRIX *mat, const char *fname, const char *name)
{
  if (!name) return (MatFileWrite(fname, mat->data, mat->rows, mat->cols, "matrix"));

  return (MatFileWrite(fname, mat->data, mat->rows, mat->cols, name));
}

MATRIX *MatlabRead(const char *fname)
{
  MATRIX *mat;
  MATFILE mf;
  FILE *fp;
  char *name;
  double **real_matrix, **imag_matrix;
  // int file_type;
  int nrows, ncols, row, col;
  float *fptr = NULL;

  if (Gdiag_no > 0) printf("Using MatlabRead2()\n");
  return (MatlabRead2(fname));
  // Should never get here

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: opening file %s\n", fname);

  fp = fopen(fname, "rb");
  if (!fp) return (NULL);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: reading header\n");

  name = MatReadHeader0(fp, &mf);
  if (name == NULL) {
    ErrorPrintf(ERROR_BADFILE, "MatlabRead: readHeader returned NULL\n");
    return (NULL);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: allocating real matrix (r=%d,c=%d)\n", (int)mf.mrows, (int)mf.ncols);

  real_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: done real mtx alloc\n");

  if (mf.imagf)
    imag_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);
  else
    imag_matrix = NULL;

  // file_type =
  readMatFile(fp, &mf, real_matrix, imag_matrix);
  nrows = (int)mf.mrows;
  ncols = (int)mf.ncols;

  mat = MatrixAlloc(nrows, ncols, mf.imagf ? MATRIX_COMPLEX : MATRIX_REAL);

  /* matrix row pointer is 1 based in both row and column */
  for (row = 1; row <= nrows; row++) {
#if 0
    fptr = MATRIX_PTR(mat,row,1) ;
#else
    if (mf.imagf) {
      fptr = (float *)MATRIX_CELT(mat, row, 1);
    }
    else
      fptr = MATRIX_RELT(mat, row, 1);
#endif

    for (col = 1; col <= ncols; col++) {
      *fptr++ = (float)real_matrix[row - 1][col - 1];
      if (mf.imagf) *fptr++ = (float)imag_matrix[row - 1][col - 1];
    }
  }

  matFree(real_matrix, nrows, ncols);
  if (mf.imagf) matFree(imag_matrix, nrows, ncols);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: done\n");

  return (mat);
}

MATFILE *MatFileRead(const char *fname, int type)
{
  MATFILE *mf;
  FILE *fp;
  // char *name;
  double **real_matrix, **imag_matrix;
  char bval;
  // int file_type;
  int nrows, ncols, row, col;
  char *cptr = NULL;
  float *fptr = NULL;

  fp = fopen(fname, "rb");
  if (!fp) return (NULL);

  mf = (MATFILE *)calloc(1, sizeof(MATFILE));
  // name =
  MatReadHeader0(fp, mf);

  real_matrix = matAlloc((int)mf->mrows, (int)mf->ncols);
  if (mf->imagf)
    imag_matrix = matAlloc((int)mf->mrows, (int)mf->ncols);
  else
    imag_matrix = NULL;

  // file_type =
  readMatFile(fp, mf, real_matrix, imag_matrix);

  nrows = (int)mf->mrows;
  ncols = (int)mf->ncols;
  /*    type = MAT_BYTE ;*/
  switch (type) {
    case MAT_BYTE:
      mf->data = (char *)calloc(nrows * ncols, sizeof(char));
      cptr = mf->data;
      break;
    case MAT_FLOAT:
      mf->data = (char *)calloc(nrows * ncols, sizeof(float));
      fptr = (float *)mf->data;
      break;
    default:
      ErrorPrintf(ERROR_BADFILE, "MatFileRead: unsupported data format %d\n", type);
      break;
  }

  for (row = 0; row < nrows; row++) {
    for (col = 0; col < ncols; col++) {
      switch (type) {
        case MAT_BYTE:
          bval = (char)real_matrix[row][col];
          *cptr++ = bval;
          break;
        case MAT_FLOAT:
          *fptr++ = (float)real_matrix[row][col];
          break;
      }
    }
  }

  matFree(real_matrix, nrows, ncols);
  if (mf->imagf) matFree(imag_matrix, nrows, ncols);

  return (mf);
}

int MatFileWrite(const char *fname, float *data, int rows, int cols, const char *name)
{
  int row, col, nitems, mtype;
  float *fptr;
  FILE *fp;
  double dval;
  MATFILE mf;

  mf.namlen = (long)strlen(name) + 1;
  mf.mrows = (long)rows;
  mf.ncols = (long)cols;
  mf.imagf = 0L;
#if (BYTE_ORDER == LITTLE_ENDIAN)
  mtype = MATFILE_PC;
#elif (BYTE_ORDER == BIG_ENDIAN)
  mtype = MATFILE_SPARC;
#else
#error 'BYTE_ORDER not set'
#endif
  mf.type = mtype + MATFILE_DOUBLE + MATFILE_FULL;

  fp = fopen(fname, "wb");
  if (!fp) {
    /*ErrorPrintf(ERROR_BADFILE,
      "MatFileWrite(%s): could not open file\n",
      fname) ;
      perror(NULL) ;*/
    return (-1);
  }

  nitems = fwrite(&mf, 1, sizeof(MATHD), fp);
  if (nitems != sizeof(MATHD)) {
    fclose(fp);
    /*ErrorPrintf(ERROR_BADFILE,
      "MatFileWrite(%s): could not write header\n",
      fname) ;
      perror(NULL) ;*/
    return (-2);
  }

  nitems = fwrite(name, sizeof(char), (int)mf.namlen, fp);
  if (nitems != (int)mf.namlen) {
    fclose(fp);
    /*ErrorPrintf(ERROR_BADFILE, "MatFileWrite(%s): could not write name\n",
      fname) ;
      perror(NULL) ;*/
    return (-3);
  }

  /* write out in column-major format */
  for (col = 0; col < cols; col++) {
    for (row = 0; row < rows; row++) {
      fptr = data + row * cols + col;
      dval = (double)(*fptr);
      nitems = fwrite(&dval, 1, sizeof(double), fp);
      if (nitems != sizeof(double)) {
        fclose(fp);
        /*ErrorPrintf(ERROR_BADFILE,
          "MatFileWrite(%s): could not write (%d,%d)\n",
          fname, row, col) ;
          perror(NULL) ;*/
        return (-2);
      }
    }
  }

  fclose(fp);
  return (0);
}

static void matFree(double **matrix, int rows, int cols)
{
  int i;
  i = cols; /* prevents warning */

  for (i = 0; i < rows; i++) free(matrix[i]);

  free(matrix);
}

static double **matAlloc(int rows, int cols)
{
  double **matrix;
  int i;

  matrix = (double **)calloc(rows, sizeof(double *));
  if (!matrix) {
    ErrorPrintf(ERROR_BADFILE, "could not allocate %d x %d matrix\n", rows, cols);
    /*exit(3) ;*/
    return (NULL);
  }

  for (i = 0; i < rows; i++) {
    matrix[i] = (double *)calloc(cols, sizeof(double));
    if (!matrix[i]) {
      ErrorPrintf(ERROR_BADFILE, "could not allocate %d x %d matrix\n", rows, cols);
      /*exit(3) ;*/
      return (NULL);
    }
  }

  return (matrix);
}

static int readMatFile(FILE *fp, MATFILE *mf, double **real_matrix, double **imag_matrix)
{
  int row, col, type, nitems;
  short sval;
  long lval;
  double dval;
  char bval;
  float fval;

  if (mf->version == 5) {
    type = (int)mf->type;
    switch (type) {
      case 9:
        type = MAT_DOUBLE;
        break;
      case 7:
        type = MAT_FLOAT;
        break;
      case 6: /* long signed */
      case 5: /* long unsigned */
        type = MAT_INT;
        break;
      case 4: /* short signed */
      case 3: /* short unsigned */
        type = MAT_SHORT;
        break;
      case 2: /* bytes signed */
      case 1: /* bytes unsigned */
        type = MAT_BYTE;
        break;
      default:
        ErrorPrintf(
            ERROR_BADFILE, "unsupported matlab format %d (%s)\n", type, type == MAT_FLOAT ? "float" : "unknown");
        break;
    }
  }
  else {                                                  /* if (mf.version == 4) */
    type = (int)mf->type - ((int)mf->type / 1000) * 1000; /* remove 1000s */
    type = type - type / 100 * 100;                       /* remove 100s from type */
    type = type / 10 * 10;                                /* 10s digit specifies data type */
    switch (type) {
      case 0:
        type = MAT_DOUBLE;
        break;
      case 10:
        type = MAT_FLOAT;
        break;
      case 20:
        type = MAT_INT;
        break;
      case 30: /* signed */
      case 40: /* unsigned */
        type = MAT_SHORT;
        break;
      case 50: /* bytes */
        type = MAT_BYTE;
        break;
      default:
        ErrorPrintf(
            ERROR_BADFILE, "unsupported matlab format %d (%s)\n", type, type == MAT_FLOAT ? "float" : "unknown");
        break;
    }
  }
  /* data is stored column by column, real first then imag (if any) */
  for (col = 0; col < mf->ncols; col++) {
    for (row = 0; row < mf->mrows; row++) {
      switch (type) {
        case MAT_BYTE:
          nitems = fread(&bval, 1, sizeof(char), fp);
          if (nitems != sizeof(char)) {
            ErrorPrintf(ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file\n", MatProgname, row, col);
            /*exit(4)*/
            return (-1);
          }
          real_matrix[row][col] = (double)bval;
          break;
        case MAT_SHORT:
          nitems = fread(&sval, 1, sizeof(short), fp);
          if (nitems != sizeof(char)) {
            ErrorPrintf(ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file\n", MatProgname, row, col);
            /*exit(4)*/
            return (-1);
          }
          if (DIFFERENT_ENDIAN(mf)) sval = swapShort(sval);
          real_matrix[row][col] = (double)sval;
          break;
        case MAT_INT: /* 32 bit integer */
          nitems = fread(&lval, 1, sizeof(long), fp);
          if (nitems != sizeof(long)) {
            ErrorPrintf(ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file\n", MatProgname, row, col);
            /*exit(4)*/
            return (-1);
          }

          if (DIFFERENT_ENDIAN(mf)) lval = swapLong32(lval);
          real_matrix[row][col] = (double)lval;
          break;
        case MAT_FLOAT:
          nitems = fread(&fval, 1, sizeof(float), fp);
          if (nitems != sizeof(float)) {
            ErrorPrintf(ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file\n", MatProgname, row, col);
            /*exit(4)*/
            return (-1);
          }
          if (DIFFERENT_ENDIAN(mf)) fval = (float)swapLong32((long32)fval);
          real_matrix[row][col] = (double)fval;
          break;
        case MAT_DOUBLE:
          nitems = fread(&dval, 1, sizeof(double), fp);
          if (nitems != sizeof(double)) {
            ErrorPrintf(ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file\n", MatProgname, row, col);
            /*exit(4)*/
            return (-1);
          }
          if (DIFFERENT_ENDIAN(mf)) dval = swapDouble(dval);
          real_matrix[row][col] = dval;
          break;
        default:
          break;
      }
    }
  }

  if (imag_matrix)
    for (col = 0; col < mf->ncols; col++) {
      for (row = 0; row < mf->mrows; row++) {
        switch (type) {
          case MAT_DOUBLE:
            nitems = fread(&dval, 1, sizeof(double), fp);
            if (nitems != sizeof(double)) {
              ErrorPrintf(
                  ERROR_BADFILE, "%s: could not read val[%d, %d] from .mat file - imag\n", MatProgname, row, col);
              /*exit(4) ;*/
              return (-1);
            }
            break;
          default:
            break;
        }
        if (DIFFERENT_ENDIAN(mf)) dval = swapDouble(dval);
        imag_matrix[row][col] = dval;
      }
    }

  return (type);
}

/*---------------------------------------------------------*/
char *MatReadHeader0(FILE *fp, MATFILE *mf)
{
  int nitems;
  char *name;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_VERBOSE, "MatReadHeader: fp=%lx, mf=%lx\n", fp, mf);

  nitems = fread(mf, 1, sizeof(MATHD), fp);
  if (nitems != sizeof(MATHD)) {
    ErrorPrintf(ERROR_BADFILE, "%s: only read %d bytes of header\n", MatProgname, nitems);
    /*exit(1) ;*/
    return (NULL);
  }
  DiagPrintf(DIAG_VERBOSE, "type = %ld\n", mf->type);
  if (DIFFERENT_ENDIAN(mf)) {
    DiagPrintf(DIAG_VERBOSE, "mat file generated with different endian\n");
    swapBytes(mf);
  }
  DiagPrintf(DIAG_VERBOSE, "after swap, type = %ld\n", mf->type);

  DiagPrintf(DIAG_VERBOSE, "MatReadHeader: nitems = %d, namelen=%d\n", nitems, (int)mf->namlen + 1);

  name = (char *)calloc((int)mf->namlen + 1, sizeof(char));
  if (!name) ErrorExit(ERROR_NO_MEMORY, "matReadHeader: couldn't allocate %d bytes for name", mf->namlen + 1);

  nitems = fread(name, sizeof(char), (int)mf->namlen, fp);
  if (nitems != mf->namlen) {
    ErrorPrintf(ERROR_BADFILE, "%s: only read %d bytes of name (%ld specified)\n", MatProgname, nitems, mf->namlen);
    /*exit(1) ;*/
    return (NULL);
  }

#if 1
  DiagPrintf(DIAG_VERBOSE,
             "MATFILE: %ld x %ld, type %ld, imagf %ld, name '%s'\n",
             mf->mrows,
             mf->ncols,
             mf->type,
             mf->imagf,
             name);
#endif

  return (name);
}

static void swapBytes(MATFILE *mf)
{
  mf->type = swapLong32(mf->type);
  mf->mrows = swapLong32(mf->mrows);
  mf->ncols = swapLong32(mf->ncols);
  mf->namlen = swapLong32(mf->namlen);
}

/*---------------------------------------------------------------
  ReadMatlabFileVariable() - this function will read all the
  variables from a matlab file and return the matrix associated
  with the given variable name.
  -------------------------------------------------------------*/
MATRIX *ReadMatlabFileVariable(const char *fname, const char *varname)
{
  MLFC *mlfc;
  int n;
  MATRIX *M;

  mlfc = ReadMatlabFileContents(fname);
  if (mlfc == NULL) return (NULL);

  M = NULL;
  for (n = 0; n < mlfc->nvars; n++) {
    if (strcmp(varname, mlfc->varname[n]) == 0) {
      M = MatrixCopy(mlfc->varmtx[n], NULL);
      break;
    }
  }

  MLFCfree(&mlfc);

  if (M == NULL) printf("ERROR: varname %s not found in %s\n", varname, fname);

  return (M);
}

/*---------------------------------------------------------------
  ReadMatlabFileContents() - this function will read all the
  variables from a matlab file, including the variable name
  and data matrix.
  -------------------------------------------------------------*/
MLFC *ReadMatlabFileContents(const char *fname)
{
  MLFC *mlfc;
  MATRIX *mat;
  MATFILE mf;
  FILE *fp;
  char *name, c;
  double **real_matrix, **imag_matrix;
  // int file_type;
  int nrows, ncols, row, col, len;
  float *fptr = NULL;

  memset(&mf, 0, sizeof(MATFILE));

  fp = fopen(fname, "rb");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n", fname);
    return (NULL);
  }

  mlfc = (MLFC *)calloc(sizeof(MLFC), 1);
  len = strlen(fname);
  mlfc->mfile = (char *)calloc(sizeof(char), len + 1);
  memmove(mlfc->mfile, fname, len);

  while (1) {
    c = fgetc(fp);
    if (c == (char)EOF)
      break;
    else
      ungetc(c, fp);

    name = MatReadHeader0(fp, &mf);
    if (name == NULL) break;

    mlfc->varname[mlfc->nvars] = name;

    real_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);
    if (mf.imagf)
      imag_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);
    else
      imag_matrix = NULL;

    // file_type =
    readMatFile(fp, &mf, real_matrix, imag_matrix);
    nrows = (int)mf.mrows;
    ncols = (int)mf.ncols;

    mat = MatrixAlloc(nrows, ncols, mf.imagf ? MATRIX_COMPLEX : MATRIX_REAL);

    /* matrix row pointer is 1 based in both row and column */
    for (row = 1; row <= nrows; row++) {
      if (mf.imagf)
        fptr = (float *)MATRIX_CELT(mat, row, 1);
      else
        fptr = MATRIX_RELT(mat, row, 1);

      for (col = 1; col <= ncols; col++) {
        *fptr++ = (float)real_matrix[row - 1][col - 1];
        if (mf.imagf) *fptr++ = (float)imag_matrix[row - 1][col - 1];
      }
    }

    matFree(real_matrix, nrows, ncols);
    if (mf.imagf) matFree(imag_matrix, nrows, ncols);

    mlfc->varmtx[mlfc->nvars] = mat;
    mlfc->nvars++;
  }

  fclose(fp);

  if (mlfc->nvars == 0) {
    printf("ERROR: did not find any variable in %s\n", fname);
    free(mlfc->mfile);
    free(mlfc);
    return (NULL);
  }

  return (mlfc);
}

/*-------------------------------------------------------*/
int MLFCprint(FILE *fp, MLFC *mlfc)
{
  int n;

  fprintf(fp, "mfile: %s\n", mlfc->mfile);
  fprintf(fp, "number of variables: %d\n", mlfc->nvars);

  for (n = 0; n < mlfc->nvars; n++) {
    fprintf(fp, "n = %d ---------------------\n", n);
    fprintf(fp, "%s\n", mlfc->varname[n]);
    MatrixPrint(fp, mlfc->varmtx[n]);
  }

  return (0);
}

/*-------------------------------------------------------*/
int MLFCfree(MLFC **ppmlfc)
{
  MLFC *mlfc;
  int n;

  mlfc = *ppmlfc;
  free(mlfc->mfile);

  for (n = 0; n < mlfc->nvars; n++) {
    free(mlfc->varname[n]);
    MatrixFree(&mlfc->varmtx[n]);
  }

  free(*ppmlfc);
  return (0);
}

/*--------------------------------------------------------------------
  znzreadMatFile() - loads data from a buffer into an array. There is
  no actual decompression that happens here, but this function is run
  as part of reading in data from a file with compressed elements.
  Handles swapping of data. The arrays must have been allocated. If
  the caller wants to read in an imaginary part, then imag_matrix
  should not be NULL (but does not check that there is an imaginary
  part to read in).
  --------------------------------------------------------------*/
static int znzreadMatFile(char *unbuff, MATFILE *mf, double **real_matrix, double **imag_matrix)
{
  int row, col, type, offset = 0;
  short sval = 1;
  long lval = 1;
  double dval = 1;
  char bval = '1';
  float fval = 1;

  type = (int)mf->type;
  switch (type) {
    case 9:
      type = MAT_DOUBLE;
      break;
    case 7:
      type = MAT_FLOAT;
      break;
    case 6: /* long signed */
    case 5: /* long unsigned */
      type = MAT_INT;
      break;
    case 4: /* short signed */
    case 3: /* short unsigned */
      type = MAT_SHORT;
      break;
    case 2: /* bytes signed */
    case 1: /* bytes unsigned */
      type = MAT_BYTE;
      break;
    default:
      ErrorPrintf(ERROR_BADFILE, "unsupported matlab format %d (%s)\n", type, type == MAT_FLOAT ? "float" : "unknown");
      break;
  }

  // unbuff has header and data in it. The header is 56 bytes long,
  // so need to skip past it. May not always be 56 (if elment name
  // is small or normal). Could be a bug.
  for (col = 0; col < mf->ncols; col++) {
    for (row = 0; row < mf->mrows; row++) {
      switch (type) {
        case MAT_BYTE:
          offset = 56 + col * (int)(mf->mrows) * sizeof(char) + row * sizeof(char);
          memmove(&bval, &unbuff[offset], sizeof(char));
          real_matrix[row][col] = (double)bval;
          break;
        case MAT_SHORT:
          offset = 56 + col * (int)(mf->mrows) * sizeof(short) + row * sizeof(short);
          memmove(&sval, &unbuff[offset], sizeof(short));
          if (DIFFERENT_ENDIAN(mf)) sval = swapShort(sval);
          real_matrix[row][col] = (double)sval;
          break;
        case MAT_INT: /* 32 bit integer */
          offset = 56 + col * (int)(mf->mrows) * sizeof(long) + row * sizeof(long);
          memmove(&lval, &unbuff[offset], sizeof(long));
          if (DIFFERENT_ENDIAN(mf)) lval = swapLong32(lval);
          real_matrix[row][col] = (double)lval;
          break;
        case MAT_FLOAT:
          offset = 56 + col * (int)(mf->mrows) * sizeof(float) + row * sizeof(float);
          memmove(&fval, &unbuff[offset], sizeof(float));
          if (DIFFERENT_ENDIAN(mf)) fval = (float)swapLong32((long32)fval);
          real_matrix[row][col] = (double)fval;
          break;
        case MAT_DOUBLE:
          offset = 56 + col * (int)(mf->mrows) * sizeof(double) + row * sizeof(double);
          memmove(&dval, &unbuff[offset], sizeof(double));
          if (DIFFERENT_ENDIAN(mf)) dval = swapDouble(dval);
          real_matrix[row][col] = dval;
          break;
        default:
          break;
      }
    }
  }
  if (imag_matrix) {
    for (col = 0; col < mf->ncols; col++) {
      for (row = 0; row < mf->mrows; row++) {
        switch (type) {
          case MAT_DOUBLE:
            offset = 56 + col * (int)(mf->mrows) * sizeof(double) + row * sizeof(double);
            offset += (mf->mrows * mf->ncols * sizeof(double));  // get past real part
            memmove(&dval, &unbuff[offset], sizeof(double));
            break;
          default:
            break;
        }
        if (DIFFERENT_ENDIAN(mf)) dval = swapDouble(dval);
        imag_matrix[row][col] = dval;
      }
    }
  }
  return (0);
}
/*--------------------------------------------------------------*/
MATRIX *MatlabRead2(const char *fname)
{
  MATRIX *mat;
  MATFILE mf;
  FILE *znzfp = NULL;
  FILE *fp = NULL;
  char *name;
  double **real_matrix, **imag_matrix;
  // int file_type;
  int nrows, ncols, row, col;
  float *fptr = NULL;
  long32 compressed = 0;
  char *unbuff = NULL;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: opening file %s\n", fname);

  fp = fopen(fname, "rb");
  if (!fp) return (NULL);

  // Get the name of first data element and the MATFILE struct and
  // whether it is compressed or not. If compressed name and mf
  // will be meaningless.
  name = MatReadHeader(fp, &mf, &compressed);

  if (compressed) {
    // Read in valid name and mf struct. Decompresses the data
    // and returns it in unbuff.
    znzfp = fopen(fname, "rb");
    if (!znzfp) return (NULL);
    name = znzMatReadHeader(znzfp, &mf, &unbuff);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: reading header\n");

  if (name == NULL) {
    ErrorPrintf(ERROR_BADFILE, "MatlabRead: readHeader returned NULL\n");
    return (NULL);
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: allocating real matrix (r=%d,c=%d)\n", (int)mf.mrows, (int)mf.ncols);

  real_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: done real mtx alloc\n");

  if (mf.imagf)
    imag_matrix = matAlloc((int)mf.mrows, (int)mf.ncols);
  else
    imag_matrix = NULL;

  if (compressed) {
    // file_type =
    znzreadMatFile(unbuff, &mf, real_matrix, imag_matrix);
  }
  else {
    // file_type =
    readMatFile(fp, &mf, real_matrix, imag_matrix);
  }

  nrows = (int)mf.mrows;
  ncols = (int)mf.ncols;

  mat = MatrixAlloc(nrows, ncols, mf.imagf ? MATRIX_COMPLEX : MATRIX_REAL);

  /* matrix row pointer is 1 based in both row and column */
  for (row = 1; row <= nrows; row++) {
#if 0
    fptr = MATRIX_PTR(mat,row,1) ;
#else
    if (mf.imagf) {
      fptr = (float *)MATRIX_CELT(mat, row, 1);
    }
    else
      fptr = MATRIX_RELT(mat, row, 1);
#endif

    for (col = 1; col <= ncols; col++) {
      *fptr++ = (float)real_matrix[row - 1][col - 1];
      if (mf.imagf) *fptr++ = (float)imag_matrix[row - 1][col - 1];
    }
  }

  matFree(real_matrix, nrows, ncols);
  if (mf.imagf) matFree(imag_matrix, nrows, ncols);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_SHOW, "MatlabRead: done\n");

  return (mat);
}
/*-------------------------------------------------------------------
  MatReadHeader() - Reads header for next data element. mf is header
  for a data element.  compressed indicates whether the file is
  compressed or not. Versions 4 and 5. Does not actually read in data,
  just fills the header for the data element. Returns the name of the
  element.  If it has elements that are compressed, then it does not
  read anything. Returns NULL if something fails or if compressed.
  -------------------------------------------------------------------*/
char *MatReadHeader(FILE *fp, MATFILE *mf, long32 *compressed)
{
  int nitems, padding;
  char *name;
  char a, b, c, d, m;
  short namlen_temp;
  long32 dt;
  char ctmp[4];
  long32 tmp;
  unsigned int fourbytes = 4;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) DiagPrintf(DIAG_VERBOSE, "MatReadHeader: fp=%lx, mf=%lx\n", fp, mf);

  name = NULL;  // shut up the compiler

  m = 1;
  m = m << 4;  // shifts first bit (the 1) to the 5th bit

  /* Test the version of MAT-file*/
  // First 4 chars indicate the version. If any are 0, it's version 4,
  // otherwise version 5.
  a = fgetc(fp);
  b = fgetc(fp);
  c = fgetc(fp);
  d = fgetc(fp);
  if ((a == (char)EOF) || (b == (char)EOF) || (c == (char)EOF) || (d == (char)EOF)) {
    ErrorPrintf(ERROR_BADFILE, "%s: could not detect the version of the Matfile\n", MatProgname);
    /*exit(1) ;*/
    return (NULL);
  }
  if (!a | !b | !c | !d)
    mf->version = 4;
  else
    mf->version = 5;
  ungetc(d, fp);
  ungetc(c, fp);
  ungetc(b, fp);
  ungetc(a, fp);

  if (mf->version == 5) {
    fseek(fp, 126, SEEK_SET);

    mf->endian = fgetc(fp);

    fseek(fp, 1, SEEK_CUR);
    if (fread(&dt, fourbytes, 1, fp) != 1 && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
    }
    if (DIFFERENT_ENDIAN(mf)) dt = swapLong32(dt);
    // dt is:
    //   14 if matrix
    //   15 if the element is compressed (cannot tell what type data is)
    //   Can assume other values if not a matrix
    if (dt != 14 && dt != 15) {
      printf("WARNING: matlab elment type is %d, which is not a matrix.\n", dt);
      fflush(stdout);
    }
    if (dt == 15)
      *compressed = 1;
    else
      *compressed = 0;

    if (*compressed) {
      // data in file are compressed, use a different function
      fclose(fp);
      return (NULL);
    }

    // Not compressed
    fseek(fp, 12, SEEK_CUR);
    if (fread(&tmp, fourbytes, 1, fp) != 1 && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
    }
    if (DIFFERENT_ENDIAN(mf)) tmp = swapLong32(tmp);
    memmove(&ctmp, &tmp, 4);  // tmp is long, ctmp is char

    mf->imagf = (long)(m & ctmp[3]);  // 0=real, 1=imag
    fseek(fp, 12, SEEK_CUR);
    if (fread(&(mf->mrows), fourbytes, 1, fp) != 1 && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
    }
    if (fread(&(mf->ncols), fourbytes, 1, fp) != 1 && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
    }
    if (DIFFERENT_ENDIAN(mf)) {
      mf->mrows = swapLong32(mf->mrows);
      mf->ncols = swapLong32(mf->ncols);
    }

    // Matlab stores the name of variable in two different ways,
    // depending on the length of the name.
    // Normal - for long names (> 4 chars)
    // Small - for < 4 chars
    // This is determined by the next bytes 3 and 4
    // (or next 1 and 2 if diff endian).
    if (DIFFERENT_ENDIAN(mf)) {
      c = fgetc(fp);
      d = fgetc(fp);
      ungetc(d, fp);
      ungetc(c, fp);
    }
    else {
      a = fgetc(fp);
      b = fgetc(fp);
      c = fgetc(fp);
      d = fgetc(fp);
      ungetc(d, fp);
      ungetc(c, fp);
      ungetc(b, fp);
      ungetc(a, fp);
    }
    if (c == 0 && d == 0) { /*normal element*/
      fseek(fp, 4, SEEK_CUR);
      if (fread(&(mf->namlen), 1, fourbytes, fp) != fourbytes && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
      }
      if (DIFFERENT_ENDIAN(mf)) mf->namlen = swapLong32(mf->namlen);
      name = (char *)calloc((int)mf->namlen, sizeof(char));
      if (fread(name, (int)mf->namlen, sizeof(char), fp) != sizeof(char) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
      }
      // If name does not fill up an 8 byte segment, read past the filler
      padding = 8 - ((int)mf->namlen - (int)mf->namlen / 8 * 8);
      fseek(fp, padding, SEEK_CUR);
    }
    else { /*small data element*/
      if (fread(&(mf->namlen), 1, fourbytes, fp) != fourbytes && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
      }
      if (DIFFERENT_ENDIAN(mf)) mf->namlen = swapLong32(mf->namlen);
      memmove(&namlen_temp, &mf->namlen, sizeof(short));
      mf->namlen = (long)namlen_temp;
      name = (char *)calloc((int)mf->namlen, sizeof(char));
      if (fread(name, (int)mf->namlen, sizeof(char), fp) != sizeof(char) && ferror(fp)) {
        ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
      }
      // If name does not fill up an 4 byte segment, read past the filler
      padding = 4 - (int)mf->namlen;
      fseek(fp, padding, SEEK_CUR);
    }

    // Read in precision type (int, float, etc)
    if (fread(&(mf->type), 1, fourbytes, fp) != fourbytes && ferror(fp)) {
      ErrorPrintf(ERROR_BADFILE, "%s: encountered error while reading Matfile\n", MatProgname);
    }
    if (DIFFERENT_ENDIAN(mf)) mf->type = swapLong32(mf->type);

    // Jump past the next 4 bytes
    fseek(fp, 4, SEEK_CUR);

#if 1
    DiagPrintf(DIAG_VERBOSE, "MATFILE5: %ld x %ld, imag %d, name '%s'\n", mf->mrows, mf->ncols, mf->imagf, name);

#endif

    return (name);
    // end version 5
  }

  // Only gets here if Version 4

  nitems = fread(mf, 1, sizeof(MATHD), fp);
  if (nitems != sizeof(MATHD)) {
    ErrorPrintf(ERROR_BADFILE, "%s: only read %d bytes of header\n", MatProgname, nitems);
    /*exit(1) ;*/
    return (NULL);
  }

  DiagPrintf(DIAG_VERBOSE, "type = %ld\n", mf->type);
  if (DIFFERENT_ENDIAN(mf)) {
    DiagPrintf(DIAG_VERBOSE, "mat file generated with different endian\n");
    swapBytes(mf);
  }
  DiagPrintf(DIAG_VERBOSE, "after swap, type = %ld\n", mf->type);

  DiagPrintf(DIAG_VERBOSE, "MatReadHeader: nitems = %d, namelen=%d\n", nitems, (int)mf->namlen + 1);

  name = (char *)calloc((int)mf->namlen + 1, sizeof(char));
  if (!name) ErrorExit(ERROR_NO_MEMORY, "matReadHeader: couldn't allocate %d bytes for name", mf->namlen + 1);

  nitems = fread(name, sizeof(char), (int)mf->namlen, fp);
  if (nitems != mf->namlen) {
    ErrorPrintf(ERROR_BADFILE, "%s: only read %d bytes of name (%ld specified)\n", MatProgname, nitems, mf->namlen);
    /*exit(1) ;*/
    return (NULL);
  }

#if 1
  DiagPrintf(DIAG_VERBOSE, "MATFILE5: %ld x %ld, imag %d, name '%s'\n", mf->mrows, mf->ncols, mf->imagf, name);
#endif
  return (name);
}

/*----------------------------------------------------------------
  znzMatReadHeader() - used to read element header and data buffer for
  a compressed element. Does not actually interpret the data as a
  matrix, just uncompresses it into a buffer.
  ----------------------------------------------------------------*/
char *znzMatReadHeader(FILE *fp, MATFILE *mf, char **data)
{
  char *name = NULL;
  char m, c;
  char endian;
  long32 dt, size;
  char *buff = NULL;
  long *ptr_size_unbuff = NULL;
  long size_unbuff;
  ptr_size_unbuff = &size_unbuff;
  long32 tmp;
  char test_data[4];
  int padding;
  short namlen_temp;
  int isitok;
  char *unbuff;

  m = 1;
  m = m << 4;  // shifts first bit (the 1) to the 5th bit

  fseek(fp, 126, SEEK_SET);
  endian = fgetc(fp);
  mf->endian = endian;
  fseek(fp, 1, SEEK_CUR);

  // dt indicates type of data
  //   14 if matrix
  //   15 if the element is compressed
  // Must be 15 to get here
  if (fread(&dt, sizeof(long32), 1, fp) != 1 && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "Encountered error while reading Matfile\n");
  }
  if (DIFFERENT_ENDIAN(mf)) dt = swapLong32(dt);

  // size = size(elementheader) + size(compresseddata)
  if (fread(&size, sizeof(long32), 1, fp) != 1 && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "Encountered error while reading Matfile\n");
  }
  if (DIFFERENT_ENDIAN(mf)) size = swapLong32(size);

  // buff holds the compressed data
  buff = (char *)malloc(sizeof(char) * size);
  if (fread(buff, 1, sizeof(char) * size, fp) != 1 && ferror(fp)) {
    ErrorPrintf(ERROR_BADFILE, "Encountered error while reading Matfile\n");
  }

  // dont know the size of the data after it has been uncompressed,
  // so just assume that it is no more than a factor of 100 larger
  size_unbuff = 100 * size;
  unbuff = (char *)calloc(size_unbuff, sizeof(char));
  isitok = uncompress((Bytef *)unbuff, (uLongf *)ptr_size_unbuff, (const Bytef *)buff, size);
  if (isitok != Z_OK) {
    printf("ERROR: in uncompressing\n");
    if (isitok == Z_MEM_ERROR) ErrorPrintf(DIAG_VERBOSE, "   not enough memory.\n");
    if (isitok == Z_BUF_ERROR) ErrorPrintf(DIAG_VERBOSE, "    not enough memory in buffer\n");
  }

  // Now unbuff has some of the header and all of the data for this
  // element.  To extract header info, need to address into unbuff the
  // same way that we addressed into the file.

  // 18 is whether there is an imaginary part or not
  memmove(&c, &unbuff[18], sizeof(char));
  mf->imagf = (long)(m & c);

  memmove(&(mf->mrows), &unbuff[32], sizeof(long));
  memmove(&(mf->ncols), &unbuff[36], sizeof(long));
  if (DIFFERENT_ENDIAN(mf)) {
    mf->mrows = swapLong32(mf->mrows);
    mf->ncols = swapLong32(mf->ncols);
  }

  // Get the name of the varible. See MatReadHeader() for docs.
  memmove(&tmp, &unbuff[40], sizeof(long32));
  if (DIFFERENT_ENDIAN(mf)) tmp = swapLong32(tmp);
  memmove(&test_data, &tmp, 4);
  if (test_data[0] == 0 && test_data[1] == 0) {
    /*normal data*/
    memmove(&(mf->namlen), &unbuff[44], sizeof(long));
    name = (char *)calloc((int)mf->namlen, sizeof(char));
    memmove(name, &unbuff[48], (int)mf->namlen * sizeof(char));
    padding = 48 + ((((int)mf->namlen - 1) / 8) + 1);
    memmove(&(mf->type), &unbuff[padding], sizeof(long));
    if (DIFFERENT_ENDIAN(mf)) mf->type = swapLong32(mf->type);
  }
  else {
    /*small data*/
    memmove(&(mf->namlen), &unbuff[40], sizeof(long));
    if (DIFFERENT_ENDIAN(mf)) mf->namlen = swapLong32(mf->namlen);
    memmove(&namlen_temp, &mf->namlen, sizeof(short));
    mf->namlen = (long)namlen_temp;
    name = (char *)calloc((int)mf->namlen, sizeof(char));
    memmove(name, &unbuff[44], (int)mf->namlen * sizeof(char));
    memmove(&(mf->type), &unbuff[48], sizeof(long));
    if (DIFFERENT_ENDIAN(mf)) mf->type = swapLong32(mf->type);
  }
  *data = unbuff;
#if 1
  ErrorPrintf(DIAG_VERBOSE,
              "MATFILE5: %ld x %ld, type %ld, imagf %d, name '%s'\n",
              mf->mrows,
              mf->ncols,
              mf->type,
              mf->imagf,
              name);

#endif
  return (name);
}
