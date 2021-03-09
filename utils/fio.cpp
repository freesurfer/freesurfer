/**
 * @brief special file read/write routines
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

#include "fio.h"
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "bfileio.h"
#include "error.h"
#include "machine.h"
#include "mghendian.h"
#include "proto.h"
#include "utils.h"  // strcpyalloc
#include "diag.h"

#define FIO_NPUSHES_MAX 100
int fio_npushes = -1;
char fio_dirstack[FIO_NPUSHES_MAX][1000];

FILE *MGHopen_file(const char *fname, const char *rwmode)
{
  FILE *f1;

  if ((f1 = fopen(fname, rwmode)) == NULL) {
    printf("Can't open %s\n", fname);
    exit(1);
  }
  return f1;
}

int putf(float f, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  f = swapFloat(f);
#endif
  return (fwrite(&f, 4, 1, fp));
}

float getf(FILE *fp)
{
  float f;

  if (fread(&f, 4, 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "getf: fread failed");
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  f = swapFloat(f);
#endif
  return f;
}

int fread1(int *v, FILE *fp)
{
  unsigned char c;
  int ret;

  ret = fread(&c, 1, 1, fp);
  *v = c;
  return (ret);
}

int fread2(int *v, FILE *fp)
{
  short s;
  int ret;

  ret = fread(&s, 2, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  *v = s;
  return (ret);
}

int fread3(int *v, FILE *fp)
{
  unsigned int i = 0;
  int ret;

  ret = fread(&i, 3, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = (unsigned int)swapInt(i);
#endif
  *v = ((i >> 8) & 0xffffff);
  return (ret);
}

int fread4(float *v, FILE *fp)
{
  float f;
  int ret;

  ret = fread(&f, 4, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  f = swapFloat(f);
#endif
  *v = f;
  return (ret);
}

int fwrite1(int v, FILE *fp)
{
  unsigned char c = (unsigned char)v;

  return (fwrite(&c, 1, 1, fp));
}

int fwrite2(int v, FILE *fp)
{
  short s;

  if (v > 0x7fff) /* don't let it overflow */
    v = 0x7fff;
  else if (v < -0x7fff)
    v = -0x7fff;
  s = (short)v;
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  return (fwrite(&s, 2, 1, fp));
}

int fwrite3(int v, FILE *fp)
{
  unsigned int i = (unsigned int)(v << 8);

#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = (unsigned int)swapInt(i);
#endif
  return (fwrite(&i, 3, 1, fp));
}

int fwrite4(int v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapInt(v);
#endif
  return (fwrite(&v, 4, 1, fp));
}

int fwriteShort(short s, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  return (fwrite(&s, sizeof(short), 1, fp));
}
double freadDouble(FILE *fp)
{
  double d;
  int ret;

  ret = fread(&d, sizeof(double), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  d = swapDouble(d);
#endif
  if (ret != 1) ErrorPrintf(ERROR_BADFILE, "freadDouble: fread failed");
  return (d);
}

int freadInt(FILE *fp)
{
  int i, err;

  if ((err = fread(&i, sizeof(int), 1, fp)) < 0 && Gdiag_no >= 0) {
      ErrorPrintf(ERROR_BADFILE, "freadInt: fread failed %d", err);
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapInt(i);
#endif
  return (i);
}
long long freadLong(FILE *fp)
{
  long long i;

  if (fread(&i, sizeof(long long), 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "freadLong: fread failed");
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapLong64(i);
#endif
  return (i);
}

short freadShort(FILE *fp)
{
  short s;

  if (fread(&s, sizeof(short), 1, fp) != 1) {
    ErrorPrintf(ERROR_BADFILE, "freadShort: fread failed");
  }

#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif

  return (s);
}

/*******************************************************/
/* read routines which can be used for fread           */
/* the usage is                                        */
/*                                                     */
/* while (fread..Ex(., fp))                            */
/*   dosomething();                                    */
/*******************************************************/
int freadFloatEx(float *pf, FILE *fp)
{
  int ret;
  ret = fread(pf, sizeof(float), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pf = swapFloat(*pf);
#endif
  return ret;
}

int freadDoubleEx(double *pd, FILE *fp)
{
  int ret;
  ret = fread(pd, sizeof(double), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pd = swapDouble(*pd);
#endif
  return ret;
}

int freadIntEx(int *pi, FILE *fp)
{
  int nread;
  nread = fread(pi, sizeof(int), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pi = swapInt(*pi); /* swapInt(int i) */
#endif
  return (nread);
}

int freadShortEx(short *ps, FILE *fp)
{
  int nread;
  nread = fread(ps, sizeof(short), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *ps = swapShort(*ps);
#endif
  return (nread);
}

/******************************************************/
int fwriteInt(int v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapInt(v);
#endif
  return (fwrite(&v, sizeof(int), 1, fp));
}

int fwriteLong(long long v, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapLong64(v);
#endif
  return (fwrite(&v, sizeof(long long), 1, fp));
}

/*----------------------------------------*/
float freadFloat(FILE *fp)
{
  char buf[4];
  float f;
  int ret;

  ret = fread(buf, 4, 1, fp);
  // ret = fread(&f,4,1,fp); // old way
  if (ret != 1) ErrorPrintf(ERROR_BADFILE, "freadFloat: fread failed");
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf, 1);
// f = swapFloat(f);  // old way
#endif
  // error: dereferencing type-punned pointer will break strict-aliasing rules:
  //  f = *((float*)buf);
  memcpy(&f, &buf, sizeof(float));
  return (f);
}
/*----------------------------------------*/
int fwriteFloat(float f, FILE *fp)
{
  int ret;
  char buf[4];
  memmove(buf, &f, 4);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf, 1);
// f = swapFloat(f);  // old way
#endif
  ret = fwrite(buf, sizeof(float), 1, fp);
  // ret = fwrite(&f,sizeof(float),1,fp);  // old way
  return (ret);
}
/*----------------------------------------*/
int fwriteDouble(double d, FILE *fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  d = swapDouble(d);
#endif
  return (fwrite(&d, sizeof(double), 1, fp));
}

/*------ znzlib support ------------*/
/* Note: an mgz file has a variable number of fields that get written at the
  end of the file. The reader keeps reading until it gets an EOF at which
  time it knows that it is done. This EOF may be interpreted as some kind
  of error, but it is not in this case.
*/
int znzread1(int *v, znzFile fp)
{
  unsigned char c;
  int ret;

  ret = znzread(&c, 1, 1, fp);
  *v = c;
  return (ret);
}

int znzread2(int *v, znzFile fp)
{
  short s;
  int ret;

  ret = znzread(&s, 2, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  *v = s;
  return (ret);
}

int znzread3(int *v, znzFile fp)
{
  unsigned int i = 0;
  int ret;

  ret = znzread(&i, 3, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = (unsigned int)swapInt(i);
#endif
  *v = ((i >> 8) & 0xffffff);
  return (ret);
}

int znzread4(float *v, znzFile fp)
{
  float f;
  int ret;

  ret = znzread(&f, 4, 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  f = swapFloat(f);
#endif
  *v = f;
  return (ret);
}

int znzwrite1(int v, znzFile fp)
{
  unsigned char c = (unsigned char)v;

  return (znzwrite(&c, 1, 1, fp));
}

int znzwrite2(int v, znzFile fp)
{
  short s;

  if (v > 0x7fff) /* don't let it overflow */
    v = 0x7fff;
  else if (v < -0x7fff)
    v = -0x7fff;
  s = (short)v;
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  return (znzwrite(&s, 2, 1, fp));
}

int znzwrite3(int v, znzFile fp)
{
  unsigned int i = (unsigned int)(v << 8);

#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = (unsigned int)swapInt(i);
#endif
  return (znzwrite(&i, 3, 1, fp));
}

int znzwrite4(int v, znzFile fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapInt(v);
#endif
  return (znzwrite(&v, 4, 1, fp));
}

int znzwriteShort(short s, znzFile fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  return (znzwrite(&s, sizeof(short), 1, fp));
}

double znzreadDouble(znzFile fp)
{
  double d;
  int ret;

  ret = znzread(&d, sizeof(double), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  d = swapDouble(d);
#endif
  if(ret != 1 && Gdiag_no >= 0) {
    // see note above about reading mgz files
    ErrorPrintf(ERROR_BADFILE, "znzreadDouble: znzread failed");
  }
  return (d);
}

int znzreadInt(znzFile fp)
{
  int i, ret;
  ret = znzread(&i, sizeof(int), 1, fp);
  if (ret != 1 && Gdiag_no >= 0){
    // see note above about reading mgz files
    ErrorPrintf(ERROR_BADFILE, "znzreadInt: znzread failed");
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapInt(i);
#endif
  return (i);
}

long long znzreadLong(znzFile fp)
{
  long long i;
  int ret;
  ret = znzread(&i, sizeof(long long), 1, fp);
  if (ret != 1 && Gdiag_no >= 0) {
    // see note above about reading mgz files
    ErrorPrintf(ERROR_BADFILE, "znzreadLong: znzread failed");
  }
#if (BYTE_ORDER == LITTLE_ENDIAN)
  i = swapLong64(i);
#endif
  return (i);
}

short znzreadShort(znzFile fp)
{
  int nread;
  short s;

  nread = znzread(&s, sizeof(short), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  s = swapShort(s);
#endif
  if(nread != 1 && Gdiag_no >= 0) {
    // see note above about reading mgz files
    ErrorPrintf(ERROR_BADFILE, "znzreadShort: znzread failed");
  }
  return (s);
}

/*******************************************************/
/* read routines which can be used for fread           */
/* the usage is                                        */
/*                                                     */
/* while (znzread..Ex(., fp))                            */
/*   dosomething();                                    */
/*******************************************************/
/*
  Note: znzreadFloatEx() is not compatible with TAGwrite()
  because TAGwrite() does not do a byte order swap. Use
  znzTAGreadFloat() instead of znzreadFloatEx().
 */
int znzreadFloatEx(float *pf, znzFile fp)
{
  int ret;
  ret = znzread(pf, sizeof(float), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pf = swapFloat(*pf);
#endif
  return ret;
}

int znzreadDoubleEx(double *pd, znzFile fp)
{
  int ret;
  ret = znzread(pd, sizeof(double), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pd = swapDouble(*pd);
#endif
  return ret;
}

int znzreadIntEx(int *pi, znzFile fp)
{
  int nread;
  nread = znzread(pi, sizeof(int), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *pi = swapInt(*pi); /* swapInt(int i) */
#endif
  return (nread);
}

int znzreadShortEx(short *ps, znzFile fp)
{
  int nread;
  nread = znzread(ps, sizeof(short), 1, fp);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  *ps = swapShort(*ps);
#endif
  return (nread);
}

/******************************************************/
int znzwriteInt(int v, znzFile fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapInt(v);
#endif
  return (znzwrite(&v, sizeof(int), 1, fp));
}

int znzwriteLong(long long v, znzFile fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapLong64(v);
#endif
  return (znzwrite(&v, sizeof(long long), 1, fp));
}

/*----------------------------------------*/
float znzreadFloat(znzFile fp)
{
  char buf[4];
  float f;
  int ret;

  ret = znzread(buf, 4, 1, fp);
  // ret = fread(&f,4,1,fp); // old way
  if (ret != 1) ErrorPrintf(ERROR_BADFILE, "freadFloat: fread failed");
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf, 1);
// f = swapFloat(f);  // old way
#endif
  // error: dereferencing type-punned pointer will break strict-aliasing rules:
  //  f = *((float*)buf);
  memcpy(&f, &buf, sizeof(float));
  return (f);
}
/*----------------------------------------*/
int znzwriteFloat(float f, znzFile fp)
{
  int ret;
  char buf[4];
  memmove(buf, &f, 4);
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf, 1);
// f = swapFloat(f);  // old way
#endif
  ret = znzwrite(buf, sizeof(float), 1, fp);
  // ret = fwrite(&f,sizeof(float),1,fp);  // old way
  return (ret);
}
/*----------------------------------------*/
int znzwriteDouble(double d, znzFile fp)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  d = swapDouble(d);
#endif
  return (znzwrite(&d, sizeof(double), 1, fp));
}

/*------------------------------------------------------
  fio_dirname() - function to replicate the functionality
  of the unix dirname.
  Author: Douglas Greve, 9/10/2001
  ------------------------------------------------------*/
char *fio_dirname(const char *pathname)
{
  int l, n;
  char *dirname;

  if (pathname == NULL) return (NULL);

  char *pname = strcpyalloc(pathname);
  l = strlen(pname);

  /* strip off leading forward slashes */
  while (l > 0 && pname[l - 1] == '/') {
    pname[l - 1] = '\0';
    l = strlen(pname);
  }

  if (l < 2) {
    /* pname is / or . or single character */
    free(pname);
    dirname = (char *)calloc(2, sizeof(char));
    if (l == 0 || pname[0] == '/')
      dirname[0] = '/';
    else
      dirname[0] = '.';
    return (dirname);
  }

  /* Start at the end of the path name and step back
     until a forward slash is found */
  for (n = l; n >= 0; n--)
    if (pname[n] == '/') break;

  if (n < 0) {
    /* no forward slash found */
    dirname = (char *)calloc(2, sizeof(char));
    dirname[0] = '.';
    free(pname);
    return (dirname);
  }

  if (n == 0) {
    /* first forward slash is the first character */
    dirname = (char *)calloc(2, sizeof(char));
    dirname[0] = '/';
    free(pname);
    return (dirname);
  }

  dirname = (char *)calloc(n + 1, sizeof(char));
  memmove(dirname, pname, n);
  free(pname);
  return (dirname);
}
/*------------------------------------------------------
  fio_basename() - function to replicate the functionality
  of the unix basename.
  Author: Douglas Greve, 9/10/2001
  ------------------------------------------------------*/
char *fio_basename(const char *pathname, const char *ext)
{
  int l, n, lext;
  char *basename, *tmp;

  if (pathname == NULL) return (NULL);

  l = strlen(pathname);
  tmp = strcpyalloc(pathname);  // keep a copy

  /* strip off the extension if it matches ext */
  if (ext != NULL) {
    lext = strlen(ext);
    if (lext < l) {
      if (strcmp(ext, &(tmp[l - lext])) == 0) {
        memset(&(tmp[l - lext]), '\0', lext + 1);
        l = strlen(tmp);
      }
    }
  }

  /* strip off leading forward slashes */
  while (l > 0 && tmp[l - 1] == '/') {
    tmp[l - 1] = '\0';
    l = strlen(tmp);
  }

  if (l < 2) {
    /* basename is / or . or single character */
    basename = (char *)calloc(2, sizeof(char));
    if (l == 0)
      basename[0] = '/';
    else
      basename[0] = tmp[0];
    free(tmp);
    return (basename);
  }

  /* Start at the end of the path name and step back
     until a forward slash is found */
  for (n = l; n >= 0; n--)
    if (tmp[n] == '/') break;

  basename = (char *)calloc(l - n, sizeof(char));
  memmove(basename, &(tmp[n + 1]), l - n);

  free(tmp);

  return (basename);
}
/*--------------------------------------------------------------
  fio_extension() - returns the extension of the given filename.
  Author: Douglas Greve, 1/30/2002
  -------------------------------------------------------------*/
char *fio_extension(const char *pathname)
{
  int lpathname, n, lext;
  char *ext;

  if (pathname == NULL) return (NULL);

  lpathname = strlen(pathname);

  lext = 0;
  n = lpathname - 1;
  while (n >= 0 && pathname[n] != '.') {
    n--;
    lext++;
  }

  /* A dot was not found, return NULL */
  if (n < 0) return (NULL);

  /* A dot was not found at the end of the file name */
  if (lext == 0) return (NULL);

  ext = (char *)calloc(sizeof(char), lext + 1);
  memmove(ext, &(pathname[n + 1]), lext);

  return (ext);
}
/* -----------------------------------------------------
  fio_DirIsWritable(const char *dirname, int fname) -- tests
  whether the given directory is writable by creating
  and deleting a junk file there. If fname != 0, then
  dirname is treated as path to a filename. It will
  return 0 if the directory does not exist.
  ----------------------------------------------------- */
int fio_DirIsWritable(const char *dirname, int fname)
{
  FILE *fp;
  char tmpstr[2000];

  if (fname != 0)
    sprintf(tmpstr, "%s.junk.54_-_sdfklj", dirname);
  else
    sprintf(tmpstr, "%s/.junk.54_-_sdfklj", dirname);

  fp = fopen(tmpstr, "w");
  if (fp == NULL) return (0);

  fclose(fp);
  unlink(tmpstr);

  return (1);
}
/*-----------------------------------------------------
  fio_FileExistsReadable() - file exists and is readable
  -----------------------------------------------------*/
int fio_FileExistsReadable(const char *fname)
{
  FILE *fp;

  fp = fopen(fname, "r");
  if (fp != NULL) {
    fclose(fp);
    return (1);
  }
  return (0);
}
/*-----------------------------------------------------
  fio_IsDirectory(fname) - fname exists and is a directory
  -----------------------------------------------------*/
int fio_IsDirectory(const char *fname)
{
  FILE *fp;
  struct stat buf;
  int err;

  fp = fopen(fname, "r");
  if (fp == NULL) return (0);
  fclose(fp);
  err = stat(fname, &buf);
  if (err != 0) return (0);
  return (S_ISDIR(buf.st_mode));
}
/*------------------------------------------------------------
  fio_NLines() - get the number of lines. The line length
  should not exceed 4000 characters.
  ------------------------------------------------------------*/
int fio_NLines(const char *fname)
{
  FILE *fp;
  int nrows;
  char tmpstring[4001];

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", fname);
    return (-1);
  }

  nrows = 0;
  while (fgets(tmpstring, 4000, fp) != NULL) nrows++;
  fclose(fp);

  return (nrows);
}

/*------------------------------------------------------------------------*/
int fio_pushd(const char *dir)
{
  extern int fio_npushes;
  extern char fio_dirstack[FIO_NPUSHES_MAX][1000];
  int err;

  fio_npushes++;
  if (fio_npushes == FIO_NPUSHES_MAX) {
    printf("ERROR: fio_pushd: maximum number of pushes reached\n");
    return (1);
  }
  if (!getcwd(fio_dirstack[fio_npushes], 1000)) {
    printf("ERROR: getcwd: no path returned\n");
  }
  err = chdir(dir);
  if (err) {
    printf("ERROR: fio_pushd: %s\n", dir);
    fio_npushes--;
    return (1);
  }

  // printf("fio_pushd: %d %s\n",fio_npushes+1,dir);

  return (0);
}
/*------------------------------------------------------------------------*/
int fio_popd(void)
{
  extern int fio_npushes;
  extern char fio_dirstack[FIO_NPUSHES_MAX][1000];
  int err;

  if (fio_npushes == -1) {
    printf("ERROR: fio_popd: dir stack is empty\n");
    return (1);
  }

  err = chdir(fio_dirstack[fio_npushes]);
  if (err) {
    printf("ERROR: fio_popd: %d %s\n", fio_npushes, fio_dirstack[fio_npushes]);
    return (1);
  }
  fio_npushes--;
  // printf("fio_popd: %d %s\n",fio_npushes+1,fio_dirstack[fio_npushes+1]);

  return (0);
}

/*--------------------------------------------------------------------
  fio_fullpath() - gets full path to a file. The file must exist. Works
  by pushing into the file dir, getting the cwd, appending the file
  basename to the cwd to get the full path, then popping the stack.
  -------------------------------------------------------------------*/
std::string fio_fullpath(const char *fname)
{
  static char cwd[1000];
  char *dirname, *basename;
  std::string fullpath;
  int err;

  basename = fio_basename(fname, NULL);
  dirname = fio_dirname(fname);

  err = fio_pushd(dirname);
  if (err) {
    free(dirname);
    free(basename);
    return (NULL);
  }
  if (!getcwd(cwd, 1000)) {
    printf("ERROR: getcwd: no path returned\n");
  }
  fio_popd();

  fullpath = std::string(cwd) + '/' + std::string(basename);

  free(dirname);
  free(basename);

  return fullpath;
}

// Replicates mkdir -p
int fio_mkdirp(const char *path, mode_t mode)
{
  int nthseg, err;
  std::string seg, path2;
  size_t n;
  
  const size_t l = strlen(path);

  n = 0;
  nthseg = 0;
  while (n < l) {
    std::string seg;
    while (n < l && path[n] != '/') {
      seg.push_back(path[n]);
      n++;
    }
    if (nthseg == 0 && path[0] != '/') {
      path2 = seg;
    } else {
      path2 = path2 + '/' + seg;
    }
    err = mkdir(path2.c_str(), mode);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n", path2.c_str());
      perror(NULL);
      return (err);
    }
    while (n < l && path[n] == '/') {
      n++;
    }
    nthseg++;
  }

  return (0);
}

/*
  \fn int fio_FileHasCarriageReturn(char *fname)
  \brief Checks whether the file has a carriage return which probably means
  that it was created under windows. Returns 0 if no CR is found or 1 if
  a CR is found (or -1 if it could not open the file).
  The CRs can be replaced with a new line with
    cat file | sed 's/\r/\n/g' > newfile
 */
int fio_FileHasCarriageReturn(const char *fname)
{
  FILE *fp;
  char c;
  int n;

  fp = fopen(fname, "r");
  if (fp == NULL) {
    printf("ERROR: fio_FileHasCarriageReturn(): cannot open %s\n", fname);
    return (-1);
  }
  n = 0;
  while (1) {
    c = fgetc(fp);
    if (c == EOF) break;
    if (c == '\r') {
      n = 1;
      break;
    }
  }
  fclose(fp);
  return (n);
}
