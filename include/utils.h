/**
 * @file  utils.h
 * @brief well....utils!
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/08/09 22:12:21 $
 *    $Revision: 1.27.2.3 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

#include "config.h"  // defines what we HAVE_...

#define MATLAB_FILE   0
#define HIPS_FILE     1
#define LIST_FILE     2
#define UNKNOWN_FILE  3

double randomNumber(double low, double hi) ;
int    setRandomSeed(long seed) ;
double normAngle(double angle) ;
float  deltaAngle(float angle1, float angle2) ;
double calcDeltaPhi(double phi1, double phi2) ;
#if 1
double latan2(double y, double x) ;
#else
#define latan2(y,x)  atan2(y,x)
#endif
float  angleDistance(float theta1, float theta2) ;
int    QuadEqual(double a1, double a2) ;
void   fComplementCode(double *pdIn, double *pdOut, int iLen) ;
#ifndef _HOME_
char *fgetl(char *s, int n, FILE *fp) ;
#endif

int  IntSqrt(int n) ;

char *StrRemoveSpaces(char *str) ;
char *StrUpper(char *str) ;
char *StrLower(char *str) ;
char *StrSkipNumber(char *str) ;
char *StrReplace(char *src, char *dst, char csrc, int cdst) ;

char *FileNameOnly(char *str, char *fname) ;
char *FileNameFromWildcard(char *inStr, char *outStr) ;
int  FileExists(char *fname) ;
int  FileType(char *fname) ;
int  FileNumber(char *fname) ;
int  FileNumberOfEntries(char *fname) ;
char *FileName(char *full_name) ;
char *FileFullName(char *full_name) ;
char *FileTmpName(char *base) ;
char *FileTmpName(char *basename) ;
void FileRename(char *inName, char *outName) ;
char *FileNameAbsolute(char *fname, char *absFname) ;
char *FileNamePath(char *fname, char *pathName) ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;
char *FileNameExtension(char *fname, char *ext) ;
char *AppendString(char *src, char *app);

int devIsinf(float value);
int devIsnan(float value);
int devFinite(float value);

int getMemoryUsed(void); // return total virtual memory used by Progname 
                     // in Kbytes. works only under Linux /proc system
void printMemoryUsed(void); // print function of the above.
char *strcpyalloc(char *str);
int  ItemsInString(char *str);
char *deblank(char *str);
char *str_toupper(char *str);
double sum2stddev(double xsum, double xsum2, int nx);
int compare_ints(const void *v1,const void *v2);
int nunqiue_int_list(int *idlist, int nlist);
int *unqiue_int_list(int *idlist, int nlist, int *nunique);

/* Necessary when Intel C/C++ compiler is used... */
void __ltoq(void);
void __qtol(void);

#ifndef S_ISDIR
#define S_ISDIR(m)      (((m) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(x) (((x) & S_IFMT) == S_IFREG)
#endif

#ifdef WIN32
#undef HAVE_STRNCASECMP
#undef HAVE_STRCASECMP
#endif

#ifndef HAVE_STRNCASECMP 
#define strncasecmp _strnicmp 
#endif
#ifndef HAVE_STRCASECMP 
#define strcasecmp _stricmp 
#endif

#ifdef WIN32
#define log2(x) log(x)/log(2)
#define exp2(x) pow(2.0,x)
#define random rand
#define srandom srand
#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))
#define snprintf _snprintf
double rint(double x);
// ascii to unsigned long long
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone 
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};
int gettimeofday(struct timeval *tv, struct timezone *tz);
double finite (double d);

#ifdef __cplusplus
extern "C" {
#endif


extern char* optarg;
extern int optind;
extern int opterr;
extern int optopt;

int getopt(int argc, char** argv, char* optstr);


#ifdef __cplusplus
}
#endif
#define MORE_MSVC_ODDITIES
#endif

#endif

