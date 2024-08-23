/*
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

#include <string>

#ifndef FIO_H
#define FIO_H

#include "znzlib.h"

FILE  *MGHopen_file(const char *fname,const char *rwmode) ;
int   putf(float f, FILE *fp) ;
float getf(FILE *fp) ;

int   fread1(int *v, FILE *fp) ;
int   fread2(int *v, FILE *fp) ;
int   fread3(int *v, FILE *fp) ;
int   fread4(float *v, FILE *fp) ;
double freadDouble(FILE *fp) ;
float freadFloat(FILE *fp) ;
int   freadInt(FILE *fp) ;
long  long freadLong(FILE *fp) ;
short freadShort(FILE *fp) ;

/* return 1 if succeed, return 0 if fail */
int freadDoubleEx(double *pd, FILE *fp) ;
int freadFloatEx(float *pf, FILE *fp) ;
int freadIntEx(int *pi, FILE *fp) ;
int freadShortEx(short *ps, FILE *fp) ;

int   fwriteDouble(double d, FILE *fp) ;
int   fwriteFloat(float f, FILE *fp) ;
int   fwriteShort(short s, FILE *fp) ;
int   fwriteInt(int v, FILE *fp) ;
int   fwriteLong(long long v, FILE *fp) ;
int   fwrite1(int v,FILE *fp) ;
int   fwrite2(int v, FILE *fp) ;
int   fwrite3(int v, FILE *fp) ;
int   fwrite4(int v, FILE *fp) ;

/* znzlib support routines */
int   znzread1(int *v, znzFile fp) ;
int   znzread2(int *v, znzFile fp) ;
int   znzread3(int *v, znzFile fp) ;
int   znzread4(float *v, znzFile fp) ;
double 	znzreadDouble  (znzFile fp) ;
float   znzreadFloat	  (znzFile fp) ;
int     znzreadInt     (znzFile fp) ;
long long znzreadLong  (znzFile fp) ;
short   znzreadShort   (znzFile fp) ;

/* return 1 if succeed, return 0 if fail */
int znzreadDoubleEx  (double *pd,  znzFile fp) ;
int znzreadFloatEx   (float *pf,   znzFile fp) ;
int znzreadIntEx     (int *pi,     znzFile fp) ;
int znzreadShortEx   (short *ps,   znzFile fp) ;

int znzwriteDouble (double d,  znzFile fp) ;
int znzwriteFloat  (float f,   znzFile fp) ;
int znzwriteShort  (short s,   znzFile fp) ;
int znzwriteUShort (unsigned short s,   znzFile fp) ;
int znzwriteInt    (int v,     znzFile fp) ;
int znzwriteLong   (long long v, znzFile fp) ;
int znzwrite1      (int v,     znzFile fp) ;
int znzwrite2      (int v,     znzFile fp) ;
int znzwrite3      (int v,     znzFile fp) ;
int znzwrite4      (int v,     znzFile fp) ;

char *fio_basename(const char *pathname,const char *ext);
char *fio_dirname(const char *pathname);
char *fio_extension(const char *pathname);
int fio_DirIsWritable(const char *dirname, int fname);
int fio_FileExistsReadable(const char *fname);
int fio_IsDirectory(const char *fname);
int fio_NLines(const char *fname);

int fio_pushd(const char *dir);
int fio_popd(void);
std::string fio_fullpath(const char *fname);
int fio_mkdirp(const char *path, mode_t mode);
int fio_FileHasCarriageReturn(const char *fname);
int makelocallink(char *src, char *link, int del);

// Code to read in a tab separated value (TSV) file. The format is
// assumed to be that the first line has a list of strings (header)
// that describes each column. After that each row has a number of
// items equal to the number of header strings.
#include<vector>
class TSV {
public:
std::vector<std::string> headers;
std::vector<std::vector<double>> data;
int debug = 0;
int ncols(void){return(headers.size());}
int nrows(void){return(data.size());}
int read(char* fname);
int print(FILE *fp, std::string fmt);
int write(std::string fname, std::string fmt);
private:
};

//#define fwriteLong(l, fp)   fwrite4((int)l, fp)

#endif
