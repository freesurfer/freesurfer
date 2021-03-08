/**
 * @brief well....utils!
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef UTILS_H
#define UTILS_H


#include <string>
#include "base.h"

#define MATLAB_FILE   0
#define HIPS_FILE     1
#define LIST_FILE     2
#define UNKNOWN_FILE  3
#define TEXT_FILE     4

double randomNumber(double low, double hi) ;
int    setRandomSeed(long seed) ;
long getRandomSeed(void);
long getRandomCalls(void);

double normAngle(double angle) ;
float  deltaAngle(float angle1, float angle2) ;
double calcDeltaPhi(double phi1, double phi2) ;

#if 1
double latan2(double y, double x) ;
#else
#define latan2(y,x)  atan2(y,x)
#endif

float fastApproxAtan2f(float y, float x);

float  angleDistance(float theta1, float theta2) ;
int    QuadEqual(double a1, double a2) ;
void   fComplementCode(double *pdIn, double *pdOut, int iLen) ;
#ifndef _HOME_
char *fgetl(char *s, int n, FILE *fp) ;
#endif

  //#ifndef isfinite 
  //#define isfinite(x) (finite(x))
  //#endif

int  IntSqrt(int n) ;

char *StrRemoveSpaces(char *str) ;
char *StrUpper(char *str) ;
char *StrLower(char *str) ;
char *StrSkipNumber(char *str) ;
char *StrReplace(const char *src, char *dst, char csrc, int cdst) ;

char *FileNameOnly(const char *str, char *fname) ;
char *FileNameFromWildcard(const char *inStr, char *outStr) ;
int  FileExists(const char *fname) ;
int  FileType(const char *fname) ;
int  FileNumber(const char *fname) ;
int  FileNumberOfEntries(const char *fname) ;
char *FileName(char *full_name) ;
char *FileFullName(char *full_name) ;
char *FileTmpName(const char *basename) ;
void FileRename(const char *inName,const char *outName) ;
char *FileNameAbsolute(const char *fname, char *absFname) ;
char *FileNamePath(const char *fname, char *pathName) ;
char *FileNameRemoveExtension(const char *in_fname, char *out_fname) ;
char *FileNameExtension(const char *fname, char *ext) ;
char *AppendString(char *src, char *app);

bool stringEndsWith(const std::string& value, const std::string& ending);

bool directoryExists(std::string const &directory);
bool directoryIsWritable(std::string const &directory);
std::string randomString(int length);
std::string getTempFile(std::string const &suffix = "");
std::string getEnvironVar(std::string const &key);

int devIsinf(float value);
int devIsnan(float value);
int devFinite(float value);

int GetVmPeak(void);
int GetVmSize(void);
int getMemoryUsed(void); // return total virtual memory used by Progname 
                     // in Kbytes. works only under Linux /proc system
void printMemoryUsed(void); // print function of the above.
char *strcpyalloc(const char *str);
int  ItemsInString(const char *str);
char *deblank(const char *str);
char *str_toupper(char *str);
double sum2stddev(double xsum, double xsum2, int nx);
int compare_ints(const void *v1,const void *v2);
int compare_floats(const void *v1,const void *v2)  ;
int nunqiue_int_list(int *idlist, int nlist);

int *unqiue_int_list(int *idlist, int nlist, int *nunique);
int most_frequent_int_list(int *idlist, int nlist, int *nmax);

/* Necessary when Intel C/C++ compiler is used... */
void __ltoq(void);
void __qtol(void);

char *GetNthItemFromString(const char *str, int nth) ;
int CountItemsInString(const char *str) ;

/* Routines for Robust Gaussians (Median and MAD) */
float kth_smallest(float a[], int n, int k);
float quick_select(float a[], int n, int k);
float median(float a[],int n);
float mad(float a[], int n);

/* define nint as a function now */
int nint( double f );
int nint2( double f ); // slightly diff implmentation, see C code

/* Outputs the help files (found in utils/fsPrintHelp.c) */
int outputHelpXml(const unsigned char *text, unsigned int size);

/* Set progress callback */
extern void (*progress_callback)(int progress);
extern int global_progress_range[2];
void SetProgressCallback(void (*callback)(int), int start, int end);
void exec_progress_callback(int slice, int total_slices, int frame, int total_frames);

int  *compute_permutation(int num, int *vec)  ;
int *GetMemUsage(int *u);
int PrintMemUsage(FILE *fp);
int PrintRUsage(int who, const char *pre, FILE *fp);
int WriteRUsage(int who, const char *pre, char *fname);
double *DListStats(double *dlist, int nlist, double *stats);

#endif
