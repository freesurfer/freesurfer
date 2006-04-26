//
// matfile.c
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/04/26 23:03:00 $
// Revision       : $Revision: 1.19 $
//
////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matfile.h"
#include "matrix.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mghendian.h"

static double **matAlloc(int rows, int ncols) ;
static void   matFree(double **matrix, int nrows, int ncols) ;
static int    readMatFile(FILE *fp, MATFILE *mf, double **real_matrix,
                          double **imag_matrix) ;
static void   swapBytes(MATFILE *mf) ;
static char *MatProgname = "matfile" ;

#if (BYTE_ORDER==LITTLE_ENDIAN)
// #ifdef Linux
#define DIFFERENT_ENDIAN(mf)   (mf->type != MATFILE_PC)
#else
#define DIFFERENT_ENDIAN(mf)   (mf->type == MATFILE_PC)
#endif


static int (*mat_printf)(const char *szFormat, ...)=NULL;
int Matlab_Install_printf( int (*new_printf)(const char *szFormat, ...) )
{
  mat_printf = new_printf;
  return(0);
}


int
MatlabWrite(MATRIX *mat, const char *fname, char *name)
{
  if (!name)
    name = "matrix" ;

  return(MatFileWrite(fname, mat->data, mat->rows, mat->cols, name)) ;
}


MATRIX *
MatlabRead(const char *fname)
{
  MATRIX   *mat ;
  MATFILE  mf ;
  FILE     *fp ;
  char     *name ;
  double   **real_matrix, **imag_matrix ;
  int      file_type, nrows, ncols, row, col ;
  float    *fptr = NULL ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: opening file %s\n",fname);

  fp = fopen(fname, "rb") ;
  if (!fp)
    return(NULL) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: reading header\n");

  name = MatReadHeader(fp, &mf) ;
  if(name==NULL)
    {
      ErrorPrintf(ERROR_BADFILE, "MatlabRead: readHeader returned NULL\n");
      return(NULL);
    }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: allocating real matrix (r=%d,c=%d)\n",
               (int)mf.mrows,(int)mf.ncols);

  real_matrix = matAlloc((int)mf.mrows, (int)mf.ncols) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: done real mtx alloc\n");

  if (mf.imagf)
    imag_matrix = matAlloc((int)mf.mrows, (int)mf.ncols) ;
  else
    imag_matrix = NULL ;

  file_type = readMatFile(fp, &mf, real_matrix, imag_matrix) ;
  nrows = (int)mf.mrows ;
  ncols = (int)mf.ncols ;

  mat = MatrixAlloc(nrows, ncols,
                    mf.imagf ? MATRIX_COMPLEX : MATRIX_REAL) ;

  /* matrix row pointer is 1 based in both row and column */
  for (row = 1 ; row <= nrows ; row++)
    {
#if 0
      fptr = MATRIX_PTR(mat,row,1) ;
#else
      if (mf.imagf)
        {
          fptr = (float *)MATRIX_CELT(mat,row,1) ;
        }
      else
        fptr = MATRIX_RELT(mat,row,1) ;
#endif

      for (col = 1 ; col <= ncols ; col++)
        {
          *fptr++ = (float)real_matrix[row-1][col-1] ;
          if (mf.imagf)
            *fptr++ = (float)imag_matrix[row-1][col-1] ;
        }
    }

  matFree(real_matrix, nrows, ncols) ;
  if (mf.imagf)
    matFree(imag_matrix, nrows, ncols) ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_SHOW, "MatlabRead: done\n");

  return(mat) ;
}


MATFILE *
MatFileRead(const char *fname, int type)
{
  MATFILE  *mf ;
  FILE     *fp ;
  char     *name ;
  double   **real_matrix, **imag_matrix ;
  char     bval ;
  int      file_type, nrows, ncols, row, col ;
  char     *cptr = NULL ;
  float    *fptr = NULL ;

  fp = fopen(fname, "rb") ;
  if (!fp)
    return(NULL) ;

  mf = (MATFILE *)calloc(1, sizeof(MATFILE)) ;
  name = MatReadHeader(fp, mf) ;

  real_matrix = matAlloc((int)mf->mrows, (int)mf->ncols) ;
  if (mf->imagf)
    imag_matrix = matAlloc((int)mf->mrows, (int)mf->ncols) ;
  else
    imag_matrix = NULL ;

  file_type = readMatFile(fp, mf, real_matrix, imag_matrix) ;

  nrows = (int)mf->mrows ;
  ncols = (int)mf->ncols ;
  /*    type = MAT_BYTE ;*/
  switch (type)
    {
    case MAT_BYTE:
      mf->data = (char *)calloc(nrows*ncols, sizeof(char)) ;
      cptr = mf->data ;
      break ;
    case MAT_FLOAT:
      mf->data = (char *)calloc(nrows*ncols, sizeof(float)) ;
      fptr = (float *)mf->data ;
      break ;
    default:
      ErrorPrintf(ERROR_BADFILE, "MatFileRead: unsupported data format %d\n",
                  type) ;
      break ;
    }

  for (row = 0 ; row < nrows ; row++)
    {
      for (col = 0 ; col < ncols ; col++)
        {
          switch (type)
            {
            case MAT_BYTE:
              bval = (char)real_matrix[row][col] ;
              *cptr++ = bval ;
              break ;
            case MAT_FLOAT:
              *fptr++ = (float)real_matrix[row][col] ;
              break ;
            }
        }
    }

  matFree(real_matrix, nrows, ncols) ;
  if (mf->imagf)
    matFree(imag_matrix, nrows, ncols) ;

  return(mf) ;
}


int
MatFileWrite(const char *fname, float *data, int rows, int cols, char *name)
{
  int     row, col, nitems, mtype ;
  float   *fptr ;
  FILE    *fp ;
  double  dval ;
  MATFILE mf ;

  mf.namlen = (long)strlen(name)+1 ;
  mf.mrows = (long)rows ;
  mf.ncols = (long)cols ;
  mf.imagf = 0L ;
#ifdef Linux
  mtype = MATFILE_PC ;
#else
  mtype = MATFILE_SPARC ;
#endif
  mf.type = mtype + MATFILE_DOUBLE + MATFILE_FULL ;

  fp = fopen(fname, "wb") ;
  if (!fp)
    {
      /*ErrorPrintf(ERROR_BADFILE,
        "MatFileWrite(%s): could not open file\n",
        fname) ;
        perror(NULL) ;*/
      return(-1) ;
    }

  nitems = fwrite(&mf, 1, sizeof(MATHD), fp) ;
  if (nitems != sizeof(MATHD))
    {
      fclose(fp) ;
      /*ErrorPrintf(ERROR_BADFILE,
        "MatFileWrite(%s): could not write header\n",
        fname) ;
        perror(NULL) ;*/
      return(-2) ;
    }

  nitems = fwrite(name, sizeof(char), (int)mf.namlen, fp) ;
  if (nitems != (int)mf.namlen)
    {
      fclose(fp) ;
      /*ErrorPrintf(ERROR_BADFILE, "MatFileWrite(%s): could not write name\n",
        fname) ;
        perror(NULL) ;*/
      return(-3) ;
    }

  /* write out in column-major format */
  for (col = 0 ; col < cols ; col++)
    {
      for (row = 0 ; row < rows ; row++)
        {
          fptr = data + row * cols + col ;
          dval = (double)(*fptr) ;
          nitems = fwrite(&dval, 1, sizeof(double), fp) ;
          if (nitems != sizeof(double))
            {
              fclose(fp) ;
              /*ErrorPrintf(ERROR_BADFILE,
                "MatFileWrite(%s): could not write (%d,%d)\n",
                fname, row, col) ;
                perror(NULL) ;*/
              return(-2) ;
            }
        }
    }

  fclose(fp) ;
  return(0) ;
}


static void
matFree(double **matrix, int rows, int cols)
{
  int    i ; i=cols; /* prevents warning */

  for (i = 0 ; i < rows ; i++)
    free(matrix[i]) ;

  free(matrix) ;
}


static double **
matAlloc(int rows, int cols)
{
  double **matrix ;
  int    i ;

  matrix = (double **)calloc(rows, sizeof(double *)) ;
  if (!matrix)
    {
      ErrorPrintf
        (ERROR_BADFILE, "could not allocate %d x %d matrix\n", rows, cols) ;
      /*exit(3) ;*/
      return(NULL);
    }

  for (i = 0 ; i < rows ; i++)
    {
      matrix[i] = (double *)calloc(cols, sizeof(double)) ;
      if (!matrix[i])
        {
          ErrorPrintf
            (ERROR_BADFILE, "could not allocate %d x %d matrix\n", rows, cols);
          /*exit(3) ;*/
          return(NULL);
        }
    }

  return(matrix) ;
}


static int
readMatFile(FILE *fp, MATFILE *mf, double **real_matrix, double **imag_matrix)
{
  int     row, col, type, nitems ;
  short   sval ;
  long    lval ;
  double  dval ;
  char    bval ;
  float   fval ;

  type = (int)mf->type - ((int)mf->type / 1000) * 1000 ;   /* remove 1000s */
  type = type - type / 100 * 100 ;               /* remove 100s from type */
  type = type / 10 * 10 ;       /* 10s digit specifies data type */
  switch (type)
    {
    case 0:
      type = MAT_DOUBLE ;
      break ;
    case 10:
      type = MAT_FLOAT ;
      break ;
    case 20:
      type = MAT_INT ;
      break ;
    case 30:  /* signed */
    case 40:  /* unsigned */
      type = MAT_SHORT ;
      break ;
    case 50:     /* bytes */
      type = MAT_BYTE ;
      break ;
    default:
      ErrorPrintf(ERROR_BADFILE, "unsupported matlab format %d (%s)\n",
                  type, type == MAT_FLOAT ? "float" : "unknown") ;
      break ;
    }

  /* data is stored column by column, real first then imag (if any) */
  for (col = 0 ; col < mf->ncols ; col++)
    {
      for (row = 0 ; row < mf->mrows ; row++)
        {
          switch(type)
            {
            case MAT_BYTE:
              nitems = fread(&bval, 1, sizeof(char), fp) ;
              if (nitems != sizeof(char))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4)*/return(-1) ;
                }
              real_matrix[row][col] = (double)bval ;
              break ;
            case MAT_SHORT:
              nitems = fread(&sval, 1, sizeof(short), fp) ;
              if (nitems != sizeof(char))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4)*/return(-1) ;
                }
              if (DIFFERENT_ENDIAN(mf))
                sval = swapShort(sval) ;
              real_matrix[row][col] = (double)sval ;
              break ;
            case MAT_INT:   /* 32 bit integer */
              nitems = fread(&lval, 1, sizeof(long), fp) ;
              if (nitems != sizeof(long))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4)*/return(-1) ;
                }

              if (DIFFERENT_ENDIAN(mf))
                lval = swapLong32(lval) ;
              real_matrix[row][col] = (double)lval ;
              break ;
            case MAT_FLOAT:
              nitems = fread(&fval, 1, sizeof(float), fp) ;
              if (nitems != sizeof(float))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4)*/return(-1) ;
                }
              if (DIFFERENT_ENDIAN(mf))
                fval = (float)swapLong32((long32)fval) ;
              real_matrix[row][col] = (double)fval ;
              break ;
            case MAT_DOUBLE:
              nitems = fread(&dval, 1, sizeof(double), fp) ;
              if (nitems != sizeof(double))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4)*/return(-1) ;
                }
              if (DIFFERENT_ENDIAN(mf))
                dval = swapDouble(dval) ;
              real_matrix[row][col] = dval ;
              break ;
            default:
              break ;
            }
        }
    }

  if (imag_matrix) for (col = 0 ; col < mf->ncols ; col++)
    {
      for (row = 0 ; row < mf->mrows ; row++)
        {
          switch(type)
            {
            case MAT_DOUBLE:
              nitems = fread(&dval, 1, sizeof(double), fp) ;
              if (nitems != sizeof(double))
                {
                  ErrorPrintf
                    (ERROR_BADFILE,
                     "%s: could not read val[%d, %d] from .mat file\n",
                     MatProgname, row, col) ;
                  /*exit(4) ;*/
                  return(-1);
                }
              break ;
            default:
              break ;
            }
          if (DIFFERENT_ENDIAN(mf))
            dval = swapDouble(dval) ;
          imag_matrix[row][col] = dval ;
        }
    }

  return(type) ;
}


char *
MatReadHeader(FILE *fp, MATFILE *mf)
{
  int   nitems ;
  char  *name ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    DiagPrintf(DIAG_VERBOSE, "MatReadHeader: fp=%lx, mf=%lx\n",fp,mf);

  nitems = fread(mf, 1, sizeof(MATHD), fp) ;
  if (nitems != sizeof(MATHD))
    {
      ErrorPrintf(ERROR_BADFILE, "%s: only read %d bytes of header\n",
                  MatProgname, nitems) ;
      /*exit(1) ;*/
      return(NULL);
    }
  DiagPrintf(DIAG_VERBOSE, "type = %ld\n", mf->type) ;
  if (DIFFERENT_ENDIAN(mf))
    {
      DiagPrintf(DIAG_VERBOSE, "mat file generated with different endian\n") ;
      swapBytes(mf) ;
    }
  DiagPrintf(DIAG_VERBOSE, "after swap, type = %ld\n", mf->type) ;

  DiagPrintf(DIAG_VERBOSE, "MatReadHeader: nitems = %d, namelen=%d\n",
             nitems,(int)mf->namlen+1);

  name = (char *)calloc((int)mf->namlen+1, sizeof(char)) ;
  if (!name)
    ErrorExit(ERROR_NO_MEMORY,
              "matReadHeader: couldn't allocate %d bytes for name",
              mf->namlen+1) ;

  nitems = fread(name, sizeof(char), (int)mf->namlen, fp) ;
  if (nitems != mf->namlen)
    {
      ErrorPrintf(ERROR_BADFILE,
                  "%s: only read %d bytes of name (%ld specified)\n",
                  MatProgname, nitems, mf->namlen) ;
      /*exit(1) ;*/
      return(NULL);
    }

#if 1
  DiagPrintf(DIAG_VERBOSE,
             "MATFILE: %ld x %ld, type %ld, imagf %ld, name '%s'\n",
             mf->mrows, mf->ncols, mf->type, mf->imagf, name) ;
#endif

  return(name) ;
}


static void
swapBytes(MATFILE *mf)
{
  mf->type = swapLong32(mf->type) ;
  mf->mrows = swapLong32(mf->mrows) ;
  mf->ncols = swapLong32(mf->ncols) ;
  mf->namlen = swapLong32(mf->namlen) ;
}


/*---------------------------------------------------------------
  ReadMatlabFileVariable() - this function will read all the
  variables from a matlab file and return the matrix associated
  with the given variable name.
  -------------------------------------------------------------*/
MATRIX *ReadMatlabFileVariable(char *fname, char *varname)
{
  MLFC *mlfc;
  int n;
  MATRIX *M;

  mlfc = ReadMatlabFileContents(fname);
  if(mlfc == NULL) return(NULL);

  M = NULL;
  for(n=0; n < mlfc->nvars; n++){
    if(strcmp(varname,mlfc->varname[n]) == 0){
      M = MatrixCopy(mlfc->varmtx[n],NULL);
      break;
    }
  }

  MLFCfree(&mlfc);

  if(M==NULL)
    printf("ERROR: varname %s not found in %s\n",varname,fname);

  return(M);
}


/*---------------------------------------------------------------
  ReadMatlabFileContents() - this function will read all the
  variables from a matlab file, including the variable name
  and data matrix.
  -------------------------------------------------------------*/
MLFC *ReadMatlabFileContents(const char *fname)
{
  MLFC     *mlfc;
  MATRIX   *mat ;
  MATFILE  mf ;
  FILE     *fp ;
  char     *name, c;
  double   **real_matrix, **imag_matrix ;
  int      file_type, nrows, ncols, row, col, len ;
  float    *fptr = NULL ;

  fp = fopen(fname, "rb") ;
  if (fp == NULL){
    printf("ERROR: could not open %s\n",fname);
    return(NULL) ;
  }

  mlfc = (MLFC *) calloc(sizeof(MLFC),1);
  len = strlen(fname);
  mlfc->mfile = (char *) calloc(sizeof(char),len+1);
  memcpy(mlfc->mfile,fname,len);

  while(1){

    c = fgetc(fp);
    if(c==EOF) break;
    else ungetc(c,fp);

    name = MatReadHeader(fp, &mf) ;
    if(name==NULL) break;

    mlfc->varname[mlfc->nvars] = name;

    real_matrix = matAlloc((int)mf.mrows, (int)mf.ncols) ;
    if (mf.imagf)  imag_matrix = matAlloc((int)mf.mrows, (int)mf.ncols) ;
    else           imag_matrix = NULL ;

    file_type = readMatFile(fp, &mf, real_matrix, imag_matrix) ;
    nrows = (int)mf.mrows ;
    ncols = (int)mf.ncols ;

    mat = MatrixAlloc(nrows, ncols, mf.imagf ? MATRIX_COMPLEX : MATRIX_REAL);

    /* matrix row pointer is 1 based in both row and column */
    for (row = 1 ; row <= nrows ; row++)
      {
        if (mf.imagf) fptr = (float *)MATRIX_CELT(mat,row,1) ;
        else         fptr = MATRIX_RELT(mat,row,1) ;

        for (col = 1 ; col <= ncols ; col++){
          *fptr++ = (float)real_matrix[row-1][col-1] ;
          if (mf.imagf)
            *fptr++ = (float)imag_matrix[row-1][col-1] ;
        }
      }

    matFree(real_matrix, nrows, ncols) ;
    if (mf.imagf) matFree(imag_matrix, nrows, ncols) ;

    mlfc->varmtx[mlfc->nvars] = mat;
    mlfc->nvars++;
  }

  fclose(fp);

  if(mlfc->nvars == 0){
    printf("ERROR: did not find any variable in %s\n",fname);
    free(mlfc->mfile);
    free(mlfc);
    return(NULL);
  }

  return(mlfc) ;
}


/*-------------------------------------------------------*/
int MLFCprint(FILE *fp, MLFC *mlfc)
{
  int n;

  fprintf(fp,"mfile: %s\n",mlfc->mfile);
  fprintf(fp,"number of variables: %d\n",mlfc->nvars);

  for(n=0; n < mlfc->nvars; n++){
    fprintf(fp,"n = %d ---------------------\n",n);
    fprintf(fp,"%s\n",mlfc->varname[n]);
    MatrixPrint(fp,mlfc->varmtx[n]);
  }

  return(0);
}


/*-------------------------------------------------------*/
int MLFCfree(MLFC **ppmlfc){
  MLFC *mlfc;
  int n;

  mlfc = *ppmlfc;
  free(mlfc->mfile);

  for(n=0; n < mlfc->nvars; n++){
    free(mlfc->varname[n]);
    MatrixFree(&mlfc->varmtx[n]);
  }

  free(*ppmlfc);
  return(0) ;
}
