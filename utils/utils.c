/**
 * @file  utils.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:40 $
 *    $Revision: 1.58 $
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


/*
  @(#)utils.c  1.12
  10/16/95
*/
/*------------------------------------------------------------------------
  File Name: utils.c

  Author: Bruce Fischl

  Created: Jan. 1994

  Description: miscellaneous utility functions

  // Warning: Do not edit the following three lines.  CVS maintains them.
  // Revision Author: $Author: nicks $
  // Revision Date  : $Date: 2006/12/29 01:49:40 $
  // Revision       : $Revision: 1.58 $

  ------------------------------------------------------------------------*/

#ifdef _HOME_
#define FZERO(d)  (fabs(d) < 0.0000001)
#endif

/*------------------------------------------------------------------------
  HEADERS
  ------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/param.h>
#include <unistd.h>

#include <time.h> /* msvc (dng) */

#include "const.h"
#include "utils.h"
#include "proto.h"
#include "error.h"
#include "image.h"
#include "macros.h"
#include "mghendian.h"
#include "numerics.h"

/*------------------------------------------------------------------------
  CONSTANTS
  ------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
  STRUCTURES
  ------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
  GLOBAL DATA
  ------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
  STATIC DATA
  ------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
  STATIC PROTOTYPES
  ------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
  FUNCTIONS
  ------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
static long idum = 0L ;
int
setRandomSeed(long seed)
{
  idum = seed ;

  // also seed the 'standard' random number generators: rand() and random()
  srand(seed);
  srandom(seed);

  // seed vnl_random thingy
  OpenRan1(&idum);

  return(NO_ERROR) ;
}

double
randomNumber(double low, double hi)
{
  double val, range ;

  if (low > hi)
  {
    val = low ;
    low = hi ;
    hi = val ;
  }

  if (idum == 0L)     /* change seed from run to run */
    idum = -1L * (long)(abs((int)time(NULL))) ;

  range = hi - low ;
  val = OpenRan1(&idum) * range + low ;

  if ((val < low) || (val > hi))
    ErrorPrintf(ERROR_BADPARM, "randomNumber(%2.1f, %2.1f) - %2.1f\n",
                (float)low, (float)hi, (float)val) ;

  return(val) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double
normAngle(double angle)
{
  while (angle > PI)
    angle -= 2.0 * PI ;

  while (angle < -PI)
    angle += 2.0 * PI ;

  return(angle) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double
calcDeltaPhi(double phi1, double phi2)
{
  double delta_phi ;

  if (phi1 < 0.0)
    phi1 += 2.0 * PI ;
  else if (phi1 > 2.0 * PI)
    phi1 -= 2.0 * PI ;

  if (phi2 < 0.0)
    phi2 += 2.0 * PI ;
  else if (phi2 > 2.0 * PI)
    phi2 -= 2.0 * PI ;

  delta_phi = (phi1 - phi2) ;


  while (delta_phi > PI)
    delta_phi -= 2.0 * PI ;

  while (delta_phi < -PI)
    delta_phi += 2.0 * PI ;

  return(delta_phi) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double
latan2(double y, double x)
{
  int oerr ;
  double val ;

  oerr = errno ;
  if (FZERO(x) && FZERO(y))
    val = 0.0 ;
  else
    val = atan2(y, x) ;
  if (val < -PI)
    val += 2.0 * PI ;
  if (val > PI)
    val -= 2.0 * PI ;

  if (oerr != errno)
    ErrorPrintf(ERROR_BADPARM,
                "error %d, y %f, x %f\n", errno, (float)y, (float)x) ;
  return(val) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
int
QuadEqual(double a1, double a2)
{
  a1 = normAngle(a1) ;
  a2 = normAngle(a2) ;
  if (fabs(a1 - a2) < RADIANS(90.0))
    return(1) ;
  else
    return(0) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
void
fComplementCode(double *pdIn, double *pdOut, int iLen)
{
  int     i, iFullLen ;
  double  d ;

  iFullLen = iLen * 2 ;
  for (i = 0 ; i < iLen ; i++, pdOut++)
  {
    d = *pdIn++ ;
    /*    *pdOut = d ;*/
    pdOut[iLen] = 1.0 - d ;
  }
}

#ifdef Darwin
void srand48(long seed) ;
void
srand48(long seed)
{
  setRandomSeed(seed) ;
}
#endif

#ifdef Darwin_not_used
double drand48(void) ;
/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double
drand48(void)
{
  int  r ;

  r = rand() ;
  return((double)r / (double)RAND_MAX) ;
}
#endif
#ifndef _HOME_
/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
fgetl(char *s, int n, FILE *fp)
{
  char *cp, *cp2 ;
  int  len ;

  do
  {
    cp = fgets(s, n, fp) ;
    if (!cp)
      return(NULL) ;

    while (isspace((int)*cp))
      cp++ ;

  }
  while (((*cp) == '#') || ((*cp) == '\n') || ((*cp) == 0)) ;

  for (cp2 = cp ; *cp2 ; cp2++)
    if (*cp2 == '#')
      *cp2 = 0 ;

  len = strlen(cp) ;
  if (cp[len-1] == '\n')  /* strip newline */
    cp[len-1] = 0 ;
  return(cp) ;
}
#endif

/*
 * IntSqrt(n) -- integer square root
 *
 * Algorithm due to Isaac Newton, takes log(n)/2 iterations.
 */
int
IntSqrt(int n)
{
  register int approx, prev;

  if (n == 0)
    return 0;

  if (n < 4)
    return 1;

  approx = n/2; /* 1st approx */

  do
  {
    /*
     * f(x) = x**2 - n
     *
     * nextx = x - f(x)/f'(x)
     *   = x - (x**2-n)/(2*x) x**2 may overflow
     *   = x - (x - n/x)/2  but this will *not*
     */

    prev = approx;
    approx = prev - (prev - n / prev)/2;

  }
  while (approx != prev);

  return(approx) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
StrUpper(char *str)
{
  char *cp ;

  for (cp = str ; *cp ; cp++)
    *cp = (char)toupper(*cp) ;

  return(str) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
StrLower(char *str)
{
  char *cp ;

  for (cp = str ; *cp ; cp++)
    *cp = (char)tolower(*cp) ;

  return(str) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  replace all occurences of the character 'csrc' in the string
  'src' with the character 'cdst' in the string dst.

  Return Values:

  ------------------------------------------------------------------------*/
char *
StrReplace(char *src, char *dst, char csrc, int cdst)
{
  char *cp_src, *cp_dst ;

  for (cp_src = src, cp_dst = dst ; *cp_src ; cp_src++, cp_dst++)
  {
    if (*cp_src == csrc)
      *cp_dst = cdst ;
    else
      *cp_dst = *cp_src ;
  }

  return(dst) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract just the file name (no path) from a string.

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
FileNameOnly(char *full_name, char *fname)
{
  char *slash, *number, *at ;

  slash = strrchr(full_name, '/') ;

  if (fname)
  {
    if (!slash)
      strcpy(fname, full_name) ;
    else
      strcpy(fname, slash+1) ;
  }
  else   /* process it in place */
  {
    fname = full_name ;
    *slash = 0 ;
  }

  number = strrchr(fname, '#') ;
  if (number)
    *number = 0 ;
  at = strrchr(fname, '@') ;
  if (at)
    *at = 0 ;

  return(fname) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine whether a given file exists or not

  Return Values:
  1 if the file exists, 0 otherwise
  ------------------------------------------------------------------------*/
int
FileExists(char *fname)
{
  FILE *fp ;
  int old_errno;

  old_errno = errno;

  fp = fopen(fname, "r") ;
  if (fp)
    fclose(fp) ;
  else
    errno = old_errno;

  return(fp != NULL) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract just the file name (no path) from a string.

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
FileName(char *full_name)
{
  char *fname, *number, *at ;

  fname = strrchr(full_name, '/') ;
  if (!fname)
    fname = full_name ;
  else
    fname++ ;   /* skip '/' */

  if (*fname == '@')
    fname++ ;

  number = strrchr(fname, '#') ;
  if (number)
    *number = 0 ;
  at = strrchr(fname, '@') ;
  if (at)
    *at = 0 ;

  return(fname) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine the file type from the name

  Return Values:

  ------------------------------------------------------------------------*/
int
FileType(char *fname)
{
  char *dot, buf[STR_LEN], *number ;

  if (*fname == '@')
    return(LIST_FILE) ;

  strcpy(buf, fname) ;
  dot = strrchr(buf, '@') ;
  number = strchr(buf, '#') ;
  if (number)
    *number = 0 ;  /* don't consider : part of extension */
  if (!dot)
    dot = strrchr(buf, '.') ;

  if (dot)
  {
    dot++ ;
    StrUpper(buf) ;
    if (!strcmp(dot, "MAT"))
      return(MATLAB_FILE) ;
    else if (!strcmp(dot, "HIPL") ||
             !strcmp(dot, "HIPS") ||
             !strcmp(dot,"HIP"))
      return(HIPS_FILE) ;
    else if (!strcmp(dot, "LST"))
      return(LIST_FILE) ;
  }
  return(UNKNOWN_FILE) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  if a number is specified return it, otherwise return -1

  Return Values:

  ------------------------------------------------------------------------*/
int
FileNumber(char *fname)
{
  char buf[STR_LEN], *number ;
  int  num ;

  strcpy(buf, fname) ;
  number = strchr(buf, '#') ;
  if (number)
    sscanf(number+1, "%d", &num) ;
  else
    num = -1 ;
  return(num) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine the number of separate entries in a file

  Return Values:

  ------------------------------------------------------------------------*/
int
FileNumberOfEntries(char *fname)
{
  int  type, num, nentries ;
  FILE *fp ;
  char buf[STR_LEN], line[2*STR_LEN], *cp ;

  strcpy(buf, fname) ;   /* we will modify fname, don't ruin callers copy */
  fname = buf ;

  num = FileNumber(fname) ;
  if (num == -1)
  {
    type = FileType(fname) ;
    switch (type)
    {
    case LIST_FILE:
      fp = fopen(FileName(fname), "rb") ;
      if (!fp)
        ErrorReturn(-1, (ERROR_NO_FILE,
                         "FileNumberOfEntries: could not open %s",
                         FileName(fname))) ;
      cp = fgetl(line, 199, fp) ;
      nentries = 0 ;
      while (cp)
      {
        sscanf(cp, "%s", buf) ;
        num = FileNumberOfEntries(buf) ;
        nentries += num ;
        cp = fgetl(line, 199, fp) ;
      }
      fclose(fp) ;

      break ;
    case HIPS_FILE:
      nentries = ImageNumFrames(fname) ;
      break ;
    case MATLAB_FILE:
    default:
      nentries = 1 ;
      break ;
    }
  }
  else
    nentries = 1 ;

  return(nentries) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract the file name including path, removing additional
  modifiers such as '@' and '#'

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *
FileFullName(char *full_name)
{
  char *fname, *number, *at ;

  fname = full_name ;
  if (*fname == '@')
    fname++ ;

  number = strrchr(fname, '#') ;
  if (number)
    *number = 0 ;
  at = strrchr(fname, '@') ;
  if (at)
    *at = 0 ;

  return(fname) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  create a temporary filename which is not currently in use

  Return Values:
  pointer to the filename
  ------------------------------------------------------------------------*/
char *
FileTmpName(char *basename)
{
  static char fname[STR_LEN] ;
  int         i ;
  FILE        *fp ;

  if (!basename)
    basename = "tmp" ;

  i = 0 ;
  do
  {
    sprintf(fname, "%s%d", basename, i++) ;
    fp = fopen(fname, "r") ;
    if (fp)
      fclose(fp) ;
  }
  while (fp) ;

  return(fname) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  create a temporary filename which is not currently in use

  Return Values:
  pointer to the filename
  ------------------------------------------------------------------------*/
void
FileRename(char *inName, char *outName)
{
#ifndef _MSDOS
  char  cmd_string[200] ;
  sprintf(cmd_string, "mv %s %s", inName, outName) ;
  system(cmd_string) ;
#else
  rename(inName,outName);
#endif
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  calculate the distance between 2 angles.

  Return Values:
  the distance.
  ------------------------------------------------------------------------*/
float
angleDistance(float theta1, float theta2)
{
  float adist ;

  adist = (float)fabs(theta1 - theta2) ;
  if (adist >= PI)
    adist = (float)fabs(adist - 2*PI) ;

  return(adist) ;
}
/*#ifndef SunOS*/
int
stricmp(char *str1, char *str2)
{
  char buf1[STR_LEN], buf2[STR_LEN] ;

  strcpy(buf1, str1) ;
  strcpy(buf2, str2) ;
  StrUpper(buf1) ;
  StrUpper(buf2) ;
  return(strcmp(buf1, buf2)) ;
}
/*#endif*/

/*------------------------------------------------------------------------
  Parameters:

  Description:
  calculate the distance between 2 angles.

  Return Values:
  remove leading spaces from a string
  ------------------------------------------------------------------------*/
char *
StrRemoveSpaces(char *str)
{
  while (isspace((int)*str))
    str++ ;

  return(str) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  ------------------------------------------------------------------------*/
#ifdef SunOS
extern char *getcwd(char *pathname, size_t size) ;
#endif

char *
FileNameAbsolute(char *fname, char *absFname)
{
  char pathname[MAXPATHLEN] ;
  int  len ;

  if (*fname == '/')
  {
    if (absFname != fname)
      strcpy(absFname, fname) ;
  }
  else   /* not already absolute */
  {
    len = strlen(fname) ;
    if (fname[len-1] == '/')
      fname[len-1] = 0 ;
#ifndef Linux
    getcwd(pathname,MAXPATHLEN-1) ;
#else
#if 0
    getcwd(pathname, MAXPATHLEN-1) ;
#else
    sprintf(pathname, ".") ;
#endif
#endif
    sprintf(absFname, "%s/%s", pathname, fname) ;
  }
  return(absFname) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract the path name from a file name and return a pointer
  to it.

  Return Values:
  ------------------------------------------------------------------------*/
char *
FileNamePath(char *fname, char *pathName)
{
  char *slash ;

  strcpy(pathName, fname) ;
  slash = strrchr(pathName, '/') ;
  if (slash)
    *slash = 0 ;          /* remove file name */
  else
#ifndef Linux
    getcwd(pathName, MAXPATHLEN-1)  ;    /* no path at all, must be cwd */
#else
#if 0
    getcwd(pathName, MAXPATHLEN-1) ;
#else
  sprintf(pathName, ".") ;
#endif
#endif

  return(pathName) ;
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  advance a string pointer past a numeric value

  Return Values:
  ------------------------------------------------------------------------*/
char *
StrSkipNumber(char *str)
{
  while (*str && isdigit((int)*str))
    str++ ;
  if (*str == '.')  /* check for floating point # */
  {
    str++ ;
    while (*str && isdigit((int)*str))
      str++ ;
  }
  while (*str && isspace((int)*str))
    str++ ;

  return(str) ;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
float
deltaAngle(float angle1, float angle2)
{
  float   delta ;

  delta = angle1 - angle2 ;
  if (delta > M_PI)
    delta = 2.0 * M_PI - delta ;
  else if (delta < -M_PI)
    delta = -2.0f * M_PI - delta ;

  return(delta) ;
}

/* ---------------------------------------------------------
   Name: AppendString()
   Author: Douglas Greve, 9/20/99
   Purpose: Appends a string to a given source string,
   reallocating the source string in the process.
   Return:  pointer to the (reallocated) source string.
   Notes:
   1. src can be null.
   2. if app is null, then the function just returns src.
   3. if src cannot be reallocated, then the function
   forces an exit.
   ---------------------------------------------------------*/
char *AppendString(char *src, char *app)
{
  int sz1 = 0, sz2 = 0;
  char *tmp;

  if (!app) return src;
  sz2 = strlen(app);

  if (src) sz1 = strlen(src);
  else
  {
    src = (char *) calloc(sizeof(char),(sz2+1));
    memcpy(src,app,sz2);
    return(src);
  }

  tmp = (char *) realloc(src, sizeof(char)*(sz1+sz2+1));
  if (!tmp)
  {
    fprintf(stderr,"ERROR: AppendString: \n");
    fprintf(stderr,"Could not realloc\n");
    fprintf(stderr,"%s:%d\n",__FILE__,__LINE__);
    exit(1);
  }
  sprintf(tmp,"%s%s",tmp,app);

  return(tmp);
}

/////////////////////////////////////////////////////////
// IEEE754 defines
// see
// http://cch.loria.fr/documentation/IEEE754/numerical_comp_guide/index.html
//
// +INF   sign 0, exponent 255, fraction 0
// -INF   sign 1, exponent 255, fraction 0
// NaN    sign ?, exponent 255, fraction non-zero
//
// sign(1), e(23-30), f(0-22)
//
/* -1 if value is -inf, 1 if val is inf, 0 otherwise */
//
#if __BYTE_ORDER == __BIG_ENDIAN
union ieee754_float
{
  float f;

  /* This is the IEEE 754 single-precision format.  */
  struct
  {
unsigned int negative:
    1;
unsigned int exponent:
    8;
unsigned int mantissa:
    23;
  }
  ieee;
};
#else // little endian
union ieee754_float
{
  float f;

  /* This is the IEEE 754 single-precision format.  */
  struct
  {
unsigned int mantissa:
    23;
unsigned int exponent:
    8;
unsigned int negative:
    1;
  }
  ieee;
};
#endif

int devIsinf(float value)
{
  unsigned int s,e,f;

  union ieee754_float v;
  v.f = value;
  s = v.ieee.negative;
  e = v.ieee.exponent;
  f = v.ieee.mantissa;

  if (e == 255 && s == 0 && f == 0)
    return(1);

  if (e == 255 && s == 1 && f == 0)
    return(-1);

  return(0);
} /* end devIsinf() */

/* isnan non-zero if NaN, 0 otherwise */
int devIsnan(float value)
{
  unsigned int s, e, f;

  union ieee754_float v;
  v.f = value;
  s = v.ieee.negative;
  e = v.ieee.exponent;
  f = v.ieee.mantissa;

  if (e == 255 && f != 0)
    return(1);

  return(0);
} /* end devIsnan() */

/* non-zero if neither infinite nor NaN, 0 otherwise */
int devFinite(float value)
{
  if (!devIsinf(value) && !devIsnan(value))
    return(1);

  return(0);
} /* end devFinite() */

char *
FileNameRemoveExtension(char *in_fname, char *out_fname)
{
  char *dot ;

  if (out_fname != in_fname)
    strcpy(out_fname, in_fname) ;
  dot = strrchr(out_fname, '.') ;
  if (dot)
    *dot = 0 ;
  return(out_fname) ;
}
char *
FileNameExtension(char *fname, char *ext)
{
  char *dot, buf[STR_LEN] ;

  ext[0] = 0 ;
  strcpy(buf, fname) ;
  dot = strrchr(buf, '.') ;
  if (dot)
    strcpy(ext, dot+1) ;

  return(ext) ;
}

#include <glob.h>

char *
FileNameFromWildcard(char *inStr, char *outStr)
{
  char *cp ;
  glob_t  gbuf ;

  if (inStr != outStr)
    strcpy(outStr, inStr) ;
  cp = strchr(inStr, '*') ;
  if (NULL != cp)
  {
    if (glob(inStr, 0, NULL, &gbuf) == 0 && gbuf.gl_pathc > 0)
      strcpy(outStr, gbuf.gl_pathv[0]) ;
  }

  return(outStr) ;
}

// return Kbytes memory used
int getMemoryUsed()
{
#ifdef Linux
  FILE *fp = 0;
  char buf[256];
  int memused = 0;
  int numassigned = 0;
  /////////////////////////////////////////////////////////////////////////
  // Linux /proc/$pid/status file memory usage information in Kbytes
  // VmSize : virtual memory usage of entire process
  // VmRSS  : resident set currently in physical memory including
  //          code, data, stack
  // VmData : virtual memory usage of heap
  // VmStk  : virtual memory usage of stack
  // VmExe  : virtual memory usage by executable and statically linked libs
  // VmLib  : virtual memory usage by dlls loaded
  /////////////////////////////////////////////////////////////////////////
  sprintf(buf, "grep -i vmdata /proc/%d/status | cut -f 2", getpid());
  errno = 0;
  fp = popen(buf, "r");
  if (fp)
  {
    numassigned = fscanf(fp, "%d", &memused);
    if (numassigned == 1)
    {
      pclose(fp);
      return memused;
    }
    else
    {
      pclose(fp);
      errno = 0;
      fprintf(stderr, "getting memoryused failed");
      return -1;
    }
  }
  if (errno)
  {
    errno = 0;
    fprintf(stderr, "getting memoryused failed");
    return -1;
  }
  return -1;  // this should never happen
#else
  static int used = 0;
  if (!used)
  {
    fprintf(stderr, "getMemoryUsed works only under Linux\n");
    used = 1;
  }
  return -1;
#endif
}

void printMemoryUsed()
{
  printf("heap used: %d Kbytes.\n", getMemoryUsed());
}

// String copy will allocation.
char *strcpyalloc(char *str)
{
  char *cpstr;
  cpstr = (char *) calloc(strlen(str)+1,sizeof(char));
  strcpy(cpstr,str);
  return(cpstr);
}


/*---------------------------------------------------------------------
  ItemsInString() - counts the number of items in a string, which an
  "item" is one or more contiguous non-white space characters. Items
  are separated by white space as determined by isspace(). These
  include \f, \n, \r, \t and \v as well as the simple space.
  *-----------------------------------------------------------*/
int ItemsInString(char *str)
{
  int items, nthchar, len;

  len = strlen(str);
  if (len == 0) return(-1);

  items = 0;
  nthchar = 0;

  // Scroll through any white space at the beginning
  while (isspace(str[nthchar]))
  {
    nthchar ++;
    if (nthchar == len) return(0); // only whitespace
  }

  // Scroll through the rest of the string
  while (1)
  {
    items++;
    while (!isspace(str[nthchar]))
    {
      // scroll thru chars in the item = nonwhitespace
      nthchar ++;
      if (nthchar == len) return(items);
    }
    while (isspace(str[nthchar]))
    {
      // scroll thru whitespace after the item
      nthchar ++;
      if (nthchar == len)  return(items);
    }
  }

  return(items); // should never get here
}
