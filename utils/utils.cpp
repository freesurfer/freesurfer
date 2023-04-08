/**
 * @brief miscellaneous utility functions
 *
 * Among other junk, the central routine for random number generation is here.
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

#ifdef _HOME_
#define FZERO(d) (fabs(d) < 0.0000001)
#endif

/*------------------------------------------------------------------------
  HEADERS
  ------------------------------------------------------------------------*/

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h> /* msvc (dng) */
#include <unistd.h>
#include <vector>
#include <random>
#include <iostream>
#include <string>
#include <sstream>
#include <iterator>

/* This should be in ctype.h, but the compiler complains */
#ifndef Darwin
#ifndef isblank
int isblank(int c);
#endif
#endif

#include "const.h"
#include "error.h"
#include "image.h"
#include "macros.h"
#include "mghendian.h"
#include "numerics.h"
#include "proto.h"
#include "romp_support.h"
#include "utils.h"

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
static bool seedHasBeenSet = false;

static long idum = 0L, nrgcalls = 0L;

int setRandomSeed(long seed)
{
  // also seed the 'standard' random number generators: rand() and random()
  srand(seed);
  srandom(seed);

  // seed vnl_random thingy
  if (idum == seed) {
    /* If you want to reseed it to the same seed you had before, you have to run OpenRan1()
       with a different seed, then run it again with the same seed as before, otherwise
       it does not restart the generator.  This change was made by DNG on March 21, 2014.
       It might break automatic tests that try to test with the same seed. This code
       is called by MRIScomputeNormals() when there is a degenerate vertex where the
       average normal is 0.*/
    long idummy = seed - 1;
    OpenRan1(&idummy);
  }
  idum = seed;
  OpenRan1(&idum);
  nrgcalls = 1;

  seedHasBeenSet = true;
  
  return (NO_ERROR);
}

long getRandomSeed(void) { return (idum); }
long getRandomCalls(void) { return (nrgcalls); }

double randomNumber(double low, double hi)
{
  double val;

  if (low > hi) {
    val = low;
    low = hi;
    hi = val;
  }

  // check if seed not set
  if (idum == 0L)
  {
    std::string seed = getEnvironVar("FREESURFER_SEED");
    if (!seed.empty())
    {
      // set seed from env variable
      std::cout << "Setting seed for random number genererator to " << seed << std::endl;
      setRandomSeed(std::stol(seed));
    }
    else
    {
      // set seed by random time
      static int laterTime = 0;
      if (!laterTime) {
        laterTime = 1; 
        char commBuffer[1024];
        FILE* commFile = fopen("/proc/self/comm", "r");
        int commSize = 0;
        if (commFile)
        {
          commSize = fread(commBuffer, 1, 1023, commFile);
          if (commSize > 0) commSize-=1; // drop the \n
          for (int i = 0; i < commSize; i++)
          {
            if (commBuffer[i] == '/') commBuffer[i] = '@';
          }
          fclose(commFile);
        }
        commBuffer[commSize] = 0;
        fprintf(stderr, "%s supposed to be reproducible but seed not set\n", commBuffer);
      }
      idum = -1L * (long)(abs((int)time(NULL)));
    }
  }
  
  double range = hi - low;
  val = OpenRan1(&idum) * range + low;
  nrgcalls++;
  if ((val < low) || (val > hi))
  {
    ErrorPrintf(ERROR_BADPARM, "randomNumber(%2.1f, %2.1f) - %2.1f\n", (float)low, (float)hi, (float)val);
  }

  return (val);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double normAngle(double angle)
{
  while (angle > PI) angle -= 2.0 * PI;

  while (angle < -PI) angle += 2.0 * PI;

  return (angle);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double calcDeltaPhi(double phi1, double phi2)
{
  double delta_phi;

  if (phi1 < 0.0)
    phi1 += 2.0 * PI;
  else if (phi1 > 2.0 * PI)
    phi1 -= 2.0 * PI;

  if (phi2 < 0.0)
    phi2 += 2.0 * PI;
  else if (phi2 > 2.0 * PI)
    phi2 -= 2.0 * PI;

  delta_phi = (phi1 - phi2);

  while (delta_phi > PI) delta_phi -= 2.0 * PI;

  while (delta_phi < -PI) delta_phi += 2.0 * PI;

  return (delta_phi);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
double latan2(double y, double x)
{
  int oerr;
  double val;

  oerr = errno;
  if (FZERO(x) && FZERO(y))
    val = 0.0;
  else
    val = atan2(y, x);
  if (val < -PI) val += 2.0 * PI;
  if (val > PI) val -= 2.0 * PI;

  if (oerr != errno) ErrorPrintf(ERROR_BADPARM, "error %d, y %f, x %f\n", errno, (float)y, (float)x);
  return (val);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
int QuadEqual(double a1, double a2)
{
  a1 = normAngle(a1);
  a2 = normAngle(a2);
  if (fabs(a1 - a2) < RADIANS(90.0))
    return (1);
  else
    return (0);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
void fComplementCode(double *pdIn, double *pdOut, int iLen)
{
  int i;
  double d;

  //  int iFullLen = iLen * 2 ;
  for (i = 0; i < iLen; i++, pdOut++) {
    d = *pdIn++;
    /*    *pdOut = d ;*/
    pdOut[iLen] = 1.0 - d;
  }
}

#ifndef _HOME_
/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *fgetl(char *s, int n, FILE *fp)
{
  char *cp, *cp2;
  int len;

  do {
    cp = fgets(s, n, fp);
    if (!cp) return (NULL);

    while (isspace((int)*cp)) cp++;

  } while (((*cp) == '#') || ((*cp) == '\n') || ((*cp) == 0));

  for (cp2 = cp; *cp2; cp2++)
    if (*cp2 == '#') *cp2 = 0;

  len = strlen(cp);
  if (cp[len - 1] == '\n') /* strip newline */
    cp[len - 1] = 0;
  return (cp);
}
#endif

/*
 * IntSqrt(n) -- integer square root
 *
 * Algorithm due to Isaac Newton, takes log(n)/2 iterations.
 */
int IntSqrt(int n)
{
  int approx, prev;

  if (n == 0) return 0;

  if (n < 4) return 1;

  approx = n / 2; /* 1st approx */

  do {
    /*
     * f(x) = x**2 - n
     *
     * nextx = x - f(x)/f'(x)
     *   = x - (x**2-n)/(2*x) x**2 may overflow
     *   = x - (x - n/x)/2  but this will *not*
     */

    prev = approx;
    approx = prev - (prev - n / prev) / 2;

  } while (approx != prev);

  return (approx);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *StrUpper(char *str)
{
  char *cp;

  for (cp = str; *cp; cp++) *cp = (char)toupper(*cp);

  return (str);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *StrLower(char *str)
{
  char *cp;

  for (cp = str; *cp; cp++) *cp = (char)tolower(*cp);

  return (str);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  replace all occurences of the character 'csrc' in the string
  'src' with the character 'cdst' in the string dst.

  Return Values:

  ------------------------------------------------------------------------*/
char *StrReplace(const char *src, char *dst, char csrc, int cdst)
{
  const char *cp_src;
  char *cp_dst;

  for (cp_src = src, cp_dst = dst; *cp_src; cp_src++, cp_dst++) {
    if (*cp_src == csrc)
      *cp_dst = (char)cdst;
    else
      *cp_dst = *cp_src;
  }

  return (dst);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract just the file name (no path) from a string.

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *FileNameOnly(const char *full_name, char *fname)
{
  const char *slash;
  char *number, *at;

  slash = strrchr(full_name, '/');

  if (fname) {
    if (!slash)
      strcpy(fname, full_name);
    else
      strcpy(fname, slash + 1);
  }
  else {
    // cannot process in place due to con
    //
    // best solution: copy full_name to fname
    //
    fname = strcpyalloc(full_name);
  }

  number = strrchr(fname, '#');
  if (number) *number = 0;
  at = strrchr(fname, '@');
  if (at) *at = 0;

  return (fname);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine whether a given file exists or not

  Return Values:
  1 if the file exists, 0 otherwise
  ------------------------------------------------------------------------*/
int FileExists(const char *fname)
{
  FILE *fp;
  int old_errno;

  old_errno = errno;

  fp = fopen(fname, "r");
  if (fp)
    fclose(fp);
  else
    errno = old_errno;

  return (fp != NULL);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract just the file name (no path) from a string.

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *FileName(char *full_name)
{
  char *fname, *number, *at;

  fname = strrchr(full_name, '/');
  if (!fname)
    fname = full_name;
  else
    fname++; /* skip '/' */

  if (*fname == '@') fname++;

  number = strrchr(fname, '#');
  if (number) *number = 0;
  at = strrchr(fname, '@');
  if (at) *at = 0;

  return (fname);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine the file type from the name

  Return Values:

  ------------------------------------------------------------------------*/
int FileType(const char *fname)
{
  char *dot, buf[STR_LEN], *number;

  if (*fname == '@') return (LIST_FILE);

  strcpy(buf, fname);
  dot = strrchr(buf, '@');
  number = strchr(buf, '#');
  if (number) *number = 0; /* don't consider : part of extension */
  if (!dot) dot = strrchr(buf, '.');

  if (dot) {
    dot++;
    StrUpper(buf);
    if (!strcmp(dot, "MAT"))
      return (MATLAB_FILE);
    else if (!strcmp(dot, "HIPL") || !strcmp(dot, "HIPS") || !strcmp(dot, "HIP"))
      return (HIPS_FILE);
    else if (!strcmp(dot, "LST"))
      return (LIST_FILE);
    else if (!strcmp(dot, "TXT"))
      return (TEXT_FILE);
  }
  return (UNKNOWN_FILE);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  if a number is specified return it, otherwise return -1

  Return Values:

  ------------------------------------------------------------------------*/
int FileNumber(const char *fname)
{
  char buf[STR_LEN], *number;
  int num;

  strcpy(buf, fname);
  number = strchr(buf, '#');
  if (number)
    sscanf(number + 1, "%d", &num);
  else
    num = -1;
  return (num);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  determine the number of separate entries in a file

  Return Values:

  ------------------------------------------------------------------------*/
int FileNumberOfEntries(const char *fname)
{
  int type, num, nentries;
  FILE *fp;
  char buf[STR_LEN], line[2 * STR_LEN], *cp;

  strcpy(buf, fname); /* we will modify fname, don't ruin callers copy */

  num = FileNumber(buf);
  if (num == -1) {
    type = FileType(buf);
    switch (type) {
      case LIST_FILE:
        fp = fopen(FileName(buf), "rb");
        if (!fp) ErrorReturn(-1, (ERROR_NO_FILE, "FileNumberOfEntries: could not open %s", FileName(buf)));
        cp = fgetl(line, 199, fp);
        nentries = 0;
        while (cp) {
          sscanf(cp, "%s", buf);
          num = FileNumberOfEntries(buf);
          nentries += num;
          cp = fgetl(line, 199, fp);
        }
        fclose(fp);

        break;
      case UNKNOWN_FILE:  // assume it is text if it has no extension or modifieer
      case TEXT_FILE:
        fp = fopen(buf, "rb");
        if (!fp) ErrorReturn(-1, (ERROR_NO_FILE, "FileNumberOfEntries: could not open %s", buf));
        cp = fgetl(line, 199, fp);
        nentries = 0;
        while (cp) {
          nentries++;
          cp = fgetl(line, 199, fp);
        }
        fclose(fp);
        break;
      case HIPS_FILE:
        nentries = ImageNumFrames(buf);
        break;
      case MATLAB_FILE:
      default:
        nentries = -1;
        break;
    }
  }
  else
    nentries = 1;

  return (nentries);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract the file name including path, removing additional
  modifiers such as '@' and '#'

  Return Values:
  nothing.
  ------------------------------------------------------------------------*/
char *FileFullName(char *full_name)
{
  char *fname, *number, *at;

  fname = full_name;
  if (*fname == '@') fname++;

  number = strrchr(fname, '#');
  if (number) *number = 0;
  at = strrchr(fname, '@');
  if (at) *at = 0;

  return (fname);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  create a temporary filename which is not currently in use

  Return Values:
  pointer to the filename
  ------------------------------------------------------------------------*/
char *FileTmpName(const char *basename)
{
  static char fname[STR_LEN];
  int i;
  FILE *fp;

  const char *realname;
  if (!basename) {
    realname = "tmp";
  } else {
    realname = basename;
  }

  i = 0;
  do {
    sprintf(fname, "%s%d", realname, i++);
    fp = fopen(fname, "r");
    if (fp) fclose(fp);
  } while (fp);

  return (fname);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  create a temporary filename which is not currently in use

  Return Values:
  pointer to the filename
  ------------------------------------------------------------------------*/
void FileRename(const char *inName, const char *outName)
{
#ifndef _MSDOS
  char cmd_string[200];
  sprintf(cmd_string, "mv %s %s", inName, outName);
  int retval = system(cmd_string);
  if (retval) {
    fprintf(stdout, "Renaming file returned: %d", retval);
  }
#else
  rename(inName, outName);
#endif
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  calculate the distance between 2 angles.

  Return Values:
  the distance.
  ------------------------------------------------------------------------*/
float angleDistance(float theta1, float theta2)
{
  float adist;

  adist = (float)fabs(theta1 - theta2);
  if (adist >= PI) adist = (float)fabs(adist - 2 * PI);

  return (adist);
}
/*#ifndef SunOS*/
#ifndef Windows_NT
int stricmp(const char *str1, const char *str2)
{
  char buf1[STR_LEN], buf2[STR_LEN];

  strcpy(buf1, str1);
  strcpy(buf2, str2);
  StrUpper(buf1);
  StrUpper(buf2);
  return (strcmp(buf1, buf2));
}
#endif
/*#endif*/

/*------------------------------------------------------------------------
  Parameters:

  Description:
  calculate the distance between 2 angles.

  Return Values:
  remove leading spaces from a string
  ------------------------------------------------------------------------*/
char *StrRemoveSpaces(char *str)
{
  while (isspace((int)*str)) str++;

  return (str);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:

  Return Values:
  ------------------------------------------------------------------------*/
#ifdef SunOS
extern char *getcwd(char *pathname, size_t size);
#endif

char *FileNameAbsolute(const char *fname, char *absFname)
{
  char pathname[MAXPATHLEN];
  int len;

  if (*fname == '/') {
    if (absFname != fname) strcpy(absFname, fname);
  }
  else /* not already absolute */
  {
    len = strlen(fname);
    char *fn = strcpyalloc(fname);
    if (fn[len - 1] == '/') fn[len - 1] = 0;
#ifndef Linux
    getcwd(pathname, MAXPATHLEN - 1);
#else
#if 0
    getcwd(pathname, MAXPATHLEN-1) ;
#else
    sprintf(pathname, ".");
#endif
#endif
    sprintf(absFname, "%s/%s", pathname, fn);
    free(fn);
  }
  return (absFname);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  extract the path name from a file name and return a pointer
  to it.

  Return Values:
  ------------------------------------------------------------------------*/
char *FileNamePath(const char *fname, char *pathName)
{
  char *slash;

  strcpy(pathName, fname);
  slash = strrchr(pathName, '/');
  if (slash)
    *slash = 0; /* remove file name */
  else
#ifndef Linux
    getcwd(pathName, MAXPATHLEN - 1); /* no path at all, must be cwd */
#else
#if 0
    getcwd(pathName, MAXPATHLEN-1) ;
#else
    sprintf(pathName, ".");
#endif
#endif

  return (pathName);
}

/*------------------------------------------------------------------------
  Parameters:

  Description:
  advance a string pointer past a numeric value

  Return Values:
  ------------------------------------------------------------------------*/
char *StrSkipNumber(char *str)
{
  while (*str && isdigit((int)*str)) str++;
  if (*str == '.') /* check for floating point # */
  {
    str++;
    while (*str && isdigit((int)*str)) str++;
  }
  while (*str && isspace((int)*str)) str++;

  return (str);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
float deltaAngle(float angle1, float angle2)
{
  float delta;

  delta = angle1 - angle2;
  if (delta > M_PI)
    delta = 2.0f * (float)M_PI - delta;
  else if (delta < -M_PI)
    delta = -2.0f * (float)M_PI - delta;

  return (delta);
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

  if (src)
    sz1 = strlen(src);
  else {
    src = (char *)calloc(sizeof(char), (sz2 + 1));
    memmove(src, app, sz2);
    return (src);
  }

  tmp = (char *)realloc(src, sizeof(char) * (sz1 + sz2 + 1));
  if (!tmp) {
    fprintf(stderr, "ERROR: AppendString: \n");
    fprintf(stderr, "Could not realloc\n");
    fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
    exit(1);
  }
  strcat(tmp, app);

  return (tmp);
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
union ieee754_float {
  float f;

  /* This is the IEEE 754 single-precision format.  */
  struct
  {
    unsigned int negative : 1;
    unsigned int exponent : 8;
    unsigned int mantissa : 23;
  } ieee;
};
#else  // little endian
union ieee754_float {
  float f;

  /* This is the IEEE 754 single-precision format.  */
  struct
  {
    unsigned int mantissa : 23;
    unsigned int exponent : 8;
    unsigned int negative : 1;
  } ieee;
};
#endif

int devIsinf(float value)
{
  unsigned int s, e, f;

  union ieee754_float v;
  v.f = value;
  s = v.ieee.negative;
  e = v.ieee.exponent;
  f = v.ieee.mantissa;

  if (e == 255 && s == 0 && f == 0) return (1);

  if (e == 255 && s == 1 && f == 0) return (-1);

  return (0);
} /* end devIsinf() */

/* isnan non-zero if NaN, 0 otherwise */
int devIsnan(float value)
{
  unsigned int e, f;
  //  unsigned int s;

  union ieee754_float v;
  v.f = value;
  //  s = v.ieee.negative;
  e = v.ieee.exponent;
  f = v.ieee.mantissa;

  if (e == 255 && f != 0) return (1);

  return (0);
} /* end devIsnan() */

/* non-zero if neither infinite nor NaN, 0 otherwise */
int devFinite(float value)
{
  if (!devIsinf(value) && !devIsnan(value)) 
    return (1);

  return (0);
} /* end devFinite() */

char *FileNameRemoveExtension(const char *in_fname, char *out_fname)
{
  char *dot;

  if (out_fname != in_fname) strcpy(out_fname, in_fname);
  dot = strrchr(out_fname, '.');
  if (dot) *dot = 0;
  return (out_fname);
}
char *FileNameExtension(const char *fname, char *ext)
{
  char *dot, buf[STR_LEN];

  ext[0] = 0;
  strcpy(buf, fname);
  dot = strrchr(buf, '.');
  if (dot) strcpy(ext, dot + 1);

  return (ext);
}

#include <glob.h>

char *FileNameFromWildcard(const char *inStr, char *outStr)
{
  const char *cp;
  glob_t gbuf;

  if (inStr != outStr) strcpy(outStr, inStr);
  cp = strchr(inStr, '*');
  if (NULL != cp) {
    if (glob(inStr, 0, NULL, &gbuf) == 0 && gbuf.gl_pathc > 0) strcpy(outStr, gbuf.gl_pathv[0]);
  }

  return (outStr);
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
  if (fp) {
    numassigned = fscanf(fp, "%d", &memused);
    if (numassigned == 1) {
      pclose(fp);
      return memused;
    }
    else {
      pclose(fp);
      errno = 0;
      fprintf(stderr, "getting memoryused failed");
      return -1;
    }
  }
  if (errno) {
    errno = 0;
    fprintf(stderr, "getting memoryused failed");
    return -1;
  }
  return -1;  // this should never happen
#else
  static int used = 0;
  if (!used) {
    fprintf(stderr, "getMemoryUsed works only under Linux\n");
    used = 1;
  }
  return -1;
#endif
}

void printMemoryUsed() { printf("heap used: %d Kbytes.\n", getMemoryUsed()); }

// String copy will allocation.
char *strcpyalloc(const char *str)
{
  char *cpstr;
  cpstr = (char *)calloc(strlen(str) + 1, sizeof(char));
  strcpy(cpstr, str);
  return (cpstr);
}

/*---------------------------------------------------------------------
  ItemsInString() - counts the number of items in a string, which an
  "item" is one or more contiguous non-white space characters. Items
  are separated by white space as determined by isspace(). These
  include \f, \n, \r, \t and \v as well as the simple space.
  *-----------------------------------------------------------*/
int ItemsInString(const char *str)
{
  int items, nthchar, len;

  len = strlen(str);
  if (len == 0) return (-1);

  items = 0;
  nthchar = 0;

  // Scroll through any white space at the beginning
  while (isspace(str[nthchar])) {
    nthchar++;
    if (nthchar == len) return (0);  // only whitespace
  }

  // Scroll through the rest of the string
  while (1) {
    items++;
    while (!isspace(str[nthchar])) {
      // scroll thru chars in the item = nonwhitespace
      nthchar++;
      if (nthchar == len) return (items);
    }
    while (isspace(str[nthchar])) {
      // scroll thru whitespace after the item
      nthchar++;
      if (nthchar == len) return (items);
    }
  }

  return (items);  // should never get here
}

/*!
  \fn char *deblank(char *str)
  \brief removes blanks from a string.
*/
char *deblank(const char *str)
{
  char *dbstr;
  int m;

  dbstr = (char *)calloc(strlen(str) + 1, sizeof(char));

  m = 0;
  for (unsigned int n = 0; n < strlen(str); n++) {
    // printf("%d %c %d\n",n,str[n],isspace(str[n]));
    if (isspace(str[n])) continue;
    dbstr[m] = str[n];
    // printf("   %d %c\n",m,dbstr[m]);
    m++;
  }
  // printf("*%s*  *%s*\n",str,dbstr);

  return (dbstr);
}
/*!
  \fn int str_toupper(char *str)
  \brief Converts each char in str to upper case. Done in place.
*/
char *str_toupper(char *str)
{
  for (unsigned int n = 0; n < strlen(str); n++) str[n] = toupper(str[n]);
  return (0);
}

/*
 * The Intel C/C++ compiler references the routines 'ltoq' and 'qtol', which
 * don't seem to exist anywhere (not even in Google-land).  So dead-end
 * routines are declared here (instead of trying to guess how long to quad
 * should be implemented)...
 */
void __ltoq(void)
{
  printf("ERROR: Attempting usage of '__ltoq' routine!\n");
  exit(1);
}
void __qtol(void)
{
  printf("ERROR: Attempting usage of '__qtol' routine!\n");
  exit(1);
}

/*------------------------------------------------------------------*/
/*!
  \fn double sum2stddev(double xsum, double xsum2, int nx)
  \brief Computes the stddev of a list of numbers given: the sum of
  the numbers (xsum), the sum of the square of the numbers (xsum2),
  and the number in the list (nx). This allows the computation of the
  stddev with only one trip through the list.
*/
double sum2stddev(double xsum, double xsum2, int nx)
{
  double xmean, xstd;
  xmean = xsum / nx;
  xstd = sqrt((xsum2 - 2 * xmean * xsum + nx * xmean * xmean) / (nx - 1));
  // printf("%g %g %d %g %g\n",xsum,xsum2,nx,xmean,xstd);
  return (xstd);
}

/* --------------------------------------------- */
/*!
  \fn int compare_ints(const void *v1,const void *v2)
  \brief Int comparison function suitable for qsort.
*/
int compare_ints(const void *v1, const void *v2)
{
  int i1, i2;
  i1 = *((int *)v1);
  i2 = *((int *)v2);
  if (i1 < i2) return (-1);
  if (i1 > i2) return (+1);
  return (0);  // equal
}
/* --------------------------------------------- */
/*!
  \fn int compare_ints(const void *v1,const void *v2)
  \brief Float comparison function suitable for qsort.
*/
int compare_floats(const void *v1, const void *v2)
{
  float i1, i2;
  i1 = *((float *)v1);
  i2 = *((float *)v2);
  if (i1 < i2) return (-1);
  if (i1 > i2) return (+1);
  return (0);  // equal
}
/* ---------------------------------------------------*/
/*!
  \fn int nunqiue_int_list(int *idlist, int nlist)
  \brief Returns/counts the number of unique items in a
  list of integers. The list will be sorted.
*/
int nunqiue_int_list(int *idlist, int nlist)
{
  int idprev, nunique, n;

  qsort(idlist, nlist, sizeof(int), compare_ints);
  nunique = 1;
  idprev = idlist[0];
  for (n = 1; n < nlist; n++) {
    if (idprev != idlist[n]) {
      nunique++;
      idprev = idlist[n];
    }
  }
  return (nunique);
}
/* ---------------------------------------------------*/
/*!
  \fn int *unqiue_int_list(int *idlist, int nlist, int *nunique)
  \brief Returns the unique items in a list of integers.
  The list will be sorted.
*/
int *unqiue_int_list(int *idlist, int nlist, int *nunique)
{
  int n, *ulist, nthu;

  /* count number of unique elements in the list,
     this also sorts the list */
  *nunique = nunqiue_int_list(idlist, nlist);

  /* alloc the unqiue list */
  ulist = (int *)calloc(sizeof(int), *nunique);

  nthu = 0;
  ulist[nthu] = idlist[0];
  for (n = 1; n < nlist; n++) {
    if (ulist[nthu] != idlist[n]) {
      nthu++;
      ulist[nthu] = idlist[n];
    }
  }
  return (ulist);
}
/* ---------------------------------------------------*/
/*!
  \fn int most_frequent_int_list(int *idlist, int nlist)
  \brief Returns the item in the list that appears most frequently.
  If there is a tie, it returns the first sorted. The list will be
  sorted.
*/
int most_frequent_int_list(int *idlist, int nlist, int *nmax)
{
  int n, *ulist, nthu, nthumax, nunique, *nper, mostfreq;

  ulist = unqiue_int_list(idlist, nlist, &nunique);
  nper = (int *)calloc(sizeof(int), nunique);
  for (nthu = 0; nthu < nunique; nthu++) {
    for (n = 0; n < nlist; n++)
      if (idlist[n] == ulist[nthu]) nper[nthu]++;
  }

  nthumax = 0;
  *nmax = nper[nthumax];
  for (nthu = 0; nthu < nunique; nthu++) {
    if (*nmax < nper[nthu]) {
      *nmax = nper[nthu];
      nthumax = nthu;
    }
  }
  mostfreq = ulist[nthumax];
  free(ulist);
  free(nper);
  return (mostfreq);
}

/*--------------------------------------------------
  CountItemsInString() returns the number of items
  in the given string, where an item is defined as
  one or more contiguous non-blank characters. Same
  as gdfCountItemsInString().
  --------------------------------------------------*/
int CountItemsInString(const char *str)
{
  int len, n, nhits;
  len = strlen(str);
  nhits = 0;
  n = 0;
  while (n < len) {
    while (n < len && isblank(str[n])) n++;
    if (n >= len) break;
    if (str[n] == '\0' || str[n] == '\n' || str[n] == '\r') break;
    while (n < len && !isblank(str[n])) n++;
    nhits++;
  }
  return (nhits);
}

/*-------------------------------------------------------------------
  GetNthItemFromString() - extracts the nth item from a string.
  An item is defined as one or more non-white space chars. If nth
  is -1, then it returns the last item. item is a string that
  must be freed by the caller. Same as gdfGetNthItemFromString().
  ------------------------------------------------------------------*/
char *GetNthItemFromString(const char *str, int nth)
{
  std::istringstream iss(str);
  std::vector<std::string> tokens{
    std::istream_iterator<std::string>{iss},
    std::istream_iterator<std::string>{}};

  if (nth < 0) nth = tokens.size() - 1;

  if (nth >= tokens.size()) {
    fs::error() << "asking for item " << nth << ", but " << tokens.size() << " items in string";
    return nullptr;
  }

  char *item = strcpyalloc(tokens[nth].c_str());
  return item;
}

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array

                Reference:

                  Author: Wirth, Niklaus
                   Title: Algorithms + data structures = programs
               Publisher: Englewood Cliffs: Prentice-Hall, 1976
    Physical description: 366 p.
                  Series: Prentice-Hall Series in Automatic Computation


                                                                        reorders a[] (of length n)
                  returns k_th smallest in array
                  can be used to compute the median

 ---------------------------------------------------------------------------*/
float kth_smallest(float a[], int n, int k)
{
  int i, j, l, m;
  float x, t;
  int kk = k - 1;

  l = 0;
  m = n - 1;
  while (l < m) {
    x = a[kk];
    i = l;
    j = m;
    do {
      while (a[i] < x) i++;
      while (x < a[j]) j--;
      if (i <= j) {
        t = a[i];
        a[i] = a[j];
        a[j] = t;
        i++;
        j--;
      }
    } while (i <= j);
    if (j < kk) l = i;
    if (kk < i) m = j;
  }
  return a[kk];
}

/*---------------------------------------------------------------------------
  * This Quickselect routine is based on the algorithm described in
  * "Numerical recipes in C", Second Edition,
  * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 ---------------------------------------------------------------------------*/
float quick_select(float arr[], int n, int k)
{
  int low, high;
  int median;
  int middle, ll, hh;
  float t;

  low = 0;
  high = n - 1;
  median = k - 1;
  for (;;) {
    if (high <= low) /* One element only */
    {
      return arr[median];
    }
    if (high == low + 1) { /* Two elements only */
      if (arr[low] > arr[high]) {
        t = arr[low];
        arr[low] = arr[high];
        arr[high] = t;
      }
      return arr[median];
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high]) {
      t = arr[middle];
      arr[middle] = arr[high];
      arr[high] = t;
    }
    if (arr[low] > arr[high]) {
      t = arr[low];
      arr[low] = arr[high];
      arr[high] = t;
    }
    if (arr[middle] > arr[low]) {
      t = arr[middle];
      arr[middle] = arr[low];
      arr[low] = t;
    }

    /* Swap low item (now in position middle) into position (low+1) */
    t = arr[middle];
    arr[middle] = arr[low + 1];
    arr[low + 1] = t;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
      do
        ll++;
      while (arr[low] > arr[ll]);
      do
        hh--;
      while (arr[hh] > arr[low]);
      if (hh < ll) break;
      t = arr[ll];
      arr[ll] = arr[hh];
      arr[hh] = t;
    }
    /* Swap middle item (in position low) back into correct position */
    t = arr[low];
    arr[low] = arr[hh];
    arr[hh] = t;
    /* Re-set active partition */
    if (hh <= median) low = ll;
    if (hh >= median) high = hh - 1;
  }
}

/*---------------------------------------------------------------------------
// Function median:
//       compute median
//       Input float array t[] of length n
//       computation in situ, t will be reordered
 ---------------------------------------------------------------------------*/
float median(float t[], int n)
{
  float q, q2;
  if (n % 2 == 1)  // odd
  {
    q = kth_smallest(t, n, (n + 1) / 2);
    return q;
  }

  /*  else even: */

  /*  q = kth_smallest(t,n,n/2);
     q2 = kth_smallest(t,n,n/2 + 1); */

  q = quick_select(t, n, n / 2);
  q2 = quick_select(t, n, n / 2 + 1);

  return 0.5 * (q + q2);
}

/*---------------------------------------------------------------------------
// Function mad:
//       compute median absolute deviation
//       (robust estimate for sigma):
//          1.4826 med_i(|r_i - med_j(r_j)|)
//
//       Input float array a[] of length n
//       computation in situ, a will be reordered
//
 ---------------------------------------------------------------------------*/
float mad(float a[], int n)
{
  float d, mm, medi;
  int i;
  float *t = (float *)calloc(n, sizeof(float));

  d = 1.4826;

  medi = median(a, n);

  if (t == NULL) ErrorExit(ERROR_NO_MEMORY, "utils.c mad(...) could not allocate memory for t");

  for (i = 0; i < n; i++) {
    t[i] = fabs(a[i] - medi);
  }

  mm = median(t, n);
  free(t);
  return d * mm;
}

/* These two functions will do rounding. The difference between them
   is that {nint(-0.5)=-1,nint(+0.5)=+1} whereas
   {nint2(-0.5)=-1,nint(+0.5)=0} In a very specific case (mapping a
   single slice into a volume where the source is aligned with the
   target but only 0.5vox offset) this caused the output slice to have
   all 0s because nint() never allowed a legal index (was always -1 or
   +1 and never 0).
 */
int nint(const double f)  { return (f < 0 ? ((int)(f - 0.5)) : ((int)(f + 0.5))); }
int nint2(const double f) { return (f < 0 ? ((int)(f - 0.5)) : ((int)(f + 0.49999999))); }

void (*progress_callback)(int) = 0;
int global_progress_range[2] = {0, 100};

/*---------------------------------------------------------------------------
// Function SetProgressCallback:
//       set call back function to respond to progress change
//
//       input start and end as the range of progress
//       default is 0 and 100
//
 ---------------------------------------------------------------------------*/
void SetProgressCallback(void (*callback)(int), int start, int end)
{
  progress_callback = callback;
  global_progress_range[0] = start;
  global_progress_range[1] = end;
}

/*---------------------------------------------------------------------------
// Function exec_progress_callback:
//       convenient function to call progress callback function
//       and set current progress
//
//       In case of single frame volume, set frame to 0 and
//       total_frames to 1
//
 ---------------------------------------------------------------------------*/
void exec_progress_callback(int slice, int total_slices, int frame, int total_frames)
{
  if (progress_callback)
    progress_callback(global_progress_range[0] + (global_progress_range[1] - global_progress_range[0]) *
                                                     (slice + total_slices * frame) / (total_slices * total_frames));
}
int *compute_permutation(int num, int *vec)
{
  int n, index, tmp;

  if (vec == NULL) vec = (int *)calloc(num, sizeof(vec[0]));
  if (vec == NULL) ErrorExit(ERROR_NOMEMORY, "compute_permutation(%d): calloc failed", num);

  for (n = 0; n < num; n++) vec[n] = n;

  for (n = 0; n < num; n++) {
    index = (int)randomNumber(0.0, (double)(num - 0.0001));
    tmp = vec[index];
    vec[index] = vec[n];
    vec[n] = tmp;
  }
  return (vec);
}

/*!
  \fn int GetVmPeak(void)
  \brief returns Vm Peak by reading
  /proc/self/status. Note that it is only as accurate as
  /proc/self/status which is a file that needs to be updated by the
  OS; not sure how often that happens. 
*/
int GetVmPeak(void)
{
  static int u[5];
  GetMemUsage(u);
  return(u[1]);
}

/*!
  \fn int GetVmSize(void)
  \brief returns Vm Size by reading /proc/self/status. Note that it is
  only as accurate as /proc/self/status which is a file that needs to
  be updated by the OS; not sure how often that happens.
*/
int GetVmSize(void)
{
  static int u[5];
  GetMemUsage(u);
  return(u[0]);
}

/*!
  \fn int *GetMemUsage(int *u)
  \brief returns Vm Size, Peak, RSS, Data, and Stk by reading
  /proc/self/status. Note that it is only as accurate as
  /proc/self/status which is a file that needs to be updated by the
  OS; not sure how often that happens. Better to use this than
  getrusage() below.
*/
int *GetMemUsage(int *u)
{
  FILE *fp;
  static char tag[500];
  int r;

  if (u == NULL) u = (int *)calloc(sizeof(int), 5);
  fp = fopen("/proc/self/status", "r");
  if (fp == NULL) return (u);  // not there on MAC

  while (1) {
    r = fscanf(fp, "%s", tag);
    if (r == EOF) break;
    if (!strcasecmp(tag, "VmSize:")) fscanf(fp, "%d", &u[0]);
    if (!strcasecmp(tag, "VmPeak:")) fscanf(fp, "%d", &u[1]);
    if (!strcasecmp(tag, "VmRSS:")) fscanf(fp, "%d", &u[2]);
    if (!strcasecmp(tag, "VmData:")) fscanf(fp, "%d", &u[3]);
    if (!strcasecmp(tag, "VmStk:")) fscanf(fp, "%d", &u[4]);
    fgets(tag, 499, fp);
  }
  fclose(fp);
  // printf("SPRDS: ");
  // int k;
  // for(k=0; k < 5; k++) printf("%d ",u[k]);
  // printf("\n");

  return (u);
}

int PrintMemUsage(FILE *fp)
{
  static int u[5];
  GetMemUsage(u);
  fprintf(fp, "MEM: Size %d  Peak: %d   RSS: %d  Data: %d  Stk: %d\n", u[0], u[1], u[2], u[3], u[4]);
  fflush(fp);
  return (0);
}

/* ---------------------------------------------------*/
/*!
  \fn int PrintRUsage(int who, char *pre, FILE *fp)
  \brief Returns usage statistics. See man page for getrusage.
  probably want: who = RUSAGE_SELF. pre is a string that will start
  each line. It appears that this will underestimate the actual
  memory usage. GetMemUsage() is probably better.
*/
int PrintRUsage(int who, const char *pre, FILE *fp)
{
  int err;
  struct rusage u;
  err = getrusage(who, &u);
  if (err) return (err);
  // user and system usage time in seconds
  fprintf(fp, "%sutimesec    %lf\n", pre, u.ru_utime.tv_sec + u.ru_utime.tv_usec / 1000000.0);
  fprintf(fp, "%sstimesec    %lf\n", pre, u.ru_stime.tv_sec + u.ru_stime.tv_usec / 1000000.0);
  fprintf(fp, "%sru_maxrss   %ld\n", pre, u.ru_maxrss);   /* maximum resident set size */
  fprintf(fp, "%sru_ixrss    %ld\n", pre, u.ru_ixrss);    /* integral shared memory size */
  fprintf(fp, "%sru_idrss    %ld\n", pre, u.ru_idrss);    /* integral unshared data size */
  fprintf(fp, "%sru_isrss    %ld\n", pre, u.ru_isrss);    /* integral unshared stack size */
  fprintf(fp, "%sru_minflt   %ld\n", pre, u.ru_minflt);   /* page reclaims */
  fprintf(fp, "%sru_majflt   %ld\n", pre, u.ru_majflt);   /* page faults */
  fprintf(fp, "%sru_nswap    %ld\n", pre, u.ru_nswap);    /* swaps */
  fprintf(fp, "%sru_inblock  %ld\n", pre, u.ru_inblock);  /* block input operations */
  fprintf(fp, "%sru_oublock  %ld\n", pre, u.ru_oublock);  /* block output operations */
  fprintf(fp, "%sru_msgsnd   %ld\n", pre, u.ru_msgsnd);   /* messages sent */
  fprintf(fp, "%sru_msgrcv   %ld\n", pre, u.ru_msgrcv);   /* messages received */
  fprintf(fp, "%sru_nsignals %ld\n", pre, u.ru_nsignals); /* signals received */
  fprintf(fp, "%sru_nvcsw    %ld\n", pre, u.ru_nvcsw);    /* voluntary context switches */
  fprintf(fp, "%sru_nivcsw   %ld\n", pre, u.ru_nivcsw);   /* involuntary context switches */
  return (0);
}
/* ---------------------------------------------------*/
/*!
  \fn int WriteRUsage(int who, char *pre, char *fname)
  \brief Writes usage statistics to file fname. See man getrusage.
  probably want: who = RUSAGE_SELF;
  pre is a string that will start each line.
*/
int WriteRUsage(int who, const char *pre, char *fname)
{
  int err;
  FILE *fp;
  fp = fopen(fname, "w");
  if (fp == NULL) return (1);
  err = PrintRUsage(who, pre, fp);
  return (err);
}


static const float     Pi = PI;
static const float halfPi = PI/2;

static float fastApproxAtan2fWkr(float ay, float ax) {     // ax >= ay

    #define tableLast 100000
    static volatile float table[tableLast + 1];
    static volatile int   tableInited;
    if (!tableInited) {
        int index;
        for (index = 0; index <= tableLast; index++) {
            table[index] = atan2(index,tableLast);
        }
        tableInited = 1;
    }
    
    int index = (ax==0) ? 0 : (tableLast*ay)/ax;
        
    return table[index];
}

float fastApproxAtan2f(float y, float x) {
    static bool volatile once;
    static bool volatile use_atanf2;
    if (!once) {
        use_atanf2 = !!getenv("FREESURFER_fastApproxAtan2f_use_atanf2");
        once = true;
        if (use_atanf2)
            fprintf(stdout, "%s:%d %s\n", __FILE__, __LINE__,
                use_atanf2?"using slow atan2f":"using table lookup for fast but slightly inaccurate atan");
    }
    if (use_atanf2) return atan2f(y,x);
    
    float  ax = fabsf(x), ay = fabsf(y);
    double r;
    if (ax >= ay) { 
        r =          fastApproxAtan2fWkr(ay,ax);
    } else {
        r = halfPi - fastApproxAtan2fWkr(ax,ay);
    }
    if (ax == x && ay == y)          ; else
    if (ax == x && ay != y) r =    -r; else
    if (ax != x && ay == y) r =  Pi-r; else
    if (ax != x && ay != y) r = -Pi+r;
    return r; 
}


// Support for reporting logic problems and other phenomona
//
typedef struct ReportEntry { 
    struct ReportEntry* next;           // list for this line
    const char*         file;           // key
    int                 check_count; 
    int                 check_elideUntil; 
    int                 report_count; 
    int                 report_elideUntil; 
    } ReportEntry;

static ReportEntry* reportEntry(const char* file, int line) {
  static ReportEntry* entries[1000000];
  cheapAssert(line <          1000000);
  ReportEntry** ep = &entries[line];
  while (*ep && strcmp((*ep)->file,file)) ep = &(*ep)->next;
  ReportEntry* e = *ep;
  if (!e) {
    e = (ReportEntry*)calloc(1,sizeof(ReportEntry)); 
    e->file              = file;
    e->check_elideUntil  = 1;
    e->report_elideUntil = 1;
    *ep = e;
  }
  return e;
}

bool spendTimeCheckingForLogicProblem(const char* file, int line) {
  ReportEntry* e = reportEntry(file, line);
  if (!e) return false;
  e->check_count++;
  if (e->check_count < e->check_elideUntil) return false;
  e->check_elideUntil = (e->check_elideUntil < 32) ? e->check_elideUntil + 1 : e->check_elideUntil*2;
  return true;
}

LogicProblemResponse copeWithLogicProblem2(
    bool* wasReported, 
    const char* envvarFixer,
    const char* msg, 
    const char* file, int line, const char* function)
{
  static bool emitReport, getenvDone;
  if (!getenvDone) { getenvDone = true;
    emitReport = !!getenv("FREESURFER_reportLogicErrors");
  }

  LogicProblemResponse response = LogicProblemResponse_old;
  if (envvarFixer && !!getenv(envvarFixer)) response = LogicProblemResponse_fix;
  
  if (wasReported) *wasReported = false;
  
  if (!emitReport) return response;
  
  ReportEntry* e = reportEntry(file, line);
  if (!e) return response;

  e->report_count++;
  if (e->report_count < e->report_elideUntil) return response;

  e->report_elideUntil = (e->report_elideUntil < 32) ? e->report_elideUntil + 1 : e->report_elideUntil*2;

  fprintf(stdout,
    "%s:%d %s - probable logic error: %s%s\n%s%s%s",
    file,line,function, msg, 
    (response == LogicProblemResponse_fix) ? ", fixing it" : ", not fixing it",
    envvarFixer?".  Set ":"",
    envvarFixer?envvarFixer:"",
    envvarFixer?" to fix it\n":"");

  if (wasReported) *wasReported = true;

  return response;
}

/*!
  \fn double *DListStats(double *dlist, int nlist, double *stats)
  \brief Computes statistics of a list of double numbers. If stats
  is NULL, it allocs.
  stats[0] = nlist;
  stats[1] = mean;
  stats[2] = stddev;
  stats[3] = min;
  stats[4] = max;
 */
double *DListStats(double *dlist, int nlist, double *stats)
{
  int n;
  long double sum, sum2;
  double mean,min,max,stddev;
  if(nlist <= 0){
    printf("ERROR: DListStats() nlist = %d\n",nlist);
    return(NULL);
  }
  if(dlist == NULL){
    printf("ERROR: DListStats() dlist is NULL\n");
    return(NULL);
  }

  if(stats==NULL) stats = (double*)calloc(sizeof(double),5);
  sum = 0;
  sum2 = 0;
  min = dlist[0];
  max = dlist[0];
  for(n=0; n < nlist; n++){
    sum  += dlist[n];
    sum2 += (dlist[n]*dlist[n]);
    if(min > dlist[n]) min = dlist[n];
    if(max < dlist[n]) max = dlist[n];
  }
  mean = sum/nlist;
  stddev = sqrt(sum2/nlist - mean*mean);
  stats[0] = nlist;
  stats[1] = mean;
  stats[2] = stddev;
  stats[3] = min;
  stats[4] = max;
  return(stats);
}


/*!
  /brief Checks whether a string ends with a substring.
*/
bool stringEndsWith(const std::string& value, const std::string& ending) {
  if (value.length() < ending.length()) return false;
  return (0 == value.compare (value.length() - ending.length(), ending.length(), ending));
}


/*!
  /brief Returns environment variable value as string. If the variable does not
  exist, an empty string is returned. 
*/
std::string getEnvironVar(std::string const &key)
{
  char * val = getenv(key.c_str());
  return val == NULL ? std::string("") : std::string(val);
}


/*!
  /brief Checks whether the directory exists.
*/
bool directoryExists(std::string const &directory) {
  struct stat info;
  if (stat(directory.c_str(), &info) != 0) {
    return false;
  } else if (S_ISDIR(info.st_mode)) {
    return true;
  }
  return false;
}


/*!
  /brief Checks whether the directory is writable.
*/
bool directoryIsWritable(std::string const &directory) {
  if (!directoryExists(directory)) return false;
  std::string filename = directory + "/tmp.XXXXXX";
  int fid = mkstemp(&filename[0]);
  if (fid < 0) return false;
  close(fid);
  unlink(filename.c_str());
  return true;
}


/*!
  /brief Returns the base directory path for temporary files. Checks for environment
  variables `FS_TMPDIR`, `TMPDIR`, `TMP`, `TEMP`, and `TEMPDIR`. If none of those are set,
  the default is `/tmp`.
*/
std::string getBaseTempDir()
{
  // check if base temporary directory has been set in the env
  std::string basedir = "";
  std::vector<std::string> varnames = { "FS_TMPDIR", "TMPDIR", "TMP", "TEMP", "TEMPDIR" };
  for (const std::string& varname : varnames) {
    basedir = getEnvironVar(varname);
    if ((!basedir.empty()) && (directoryIsWritable(basedir))) break;
  }

  // default to /tmp
  if (basedir.empty()) basedir = "/tmp";
  if (!directoryIsWritable(basedir)) fs::fatal() << "Can not write to temporary directory " << basedir;

  return basedir;
}


/*!
  /brief Generates a unique, temporary file. Suffix is optional.
*/
std::string makeTempFile(std::string const &suffix)
{
  std::string basedir = getBaseTempDir();
  std::string filename = basedir + "/tmp.XXXXXX" + suffix;
  int fid = mkstemps(&filename[0], suffix.length());
  if (fid < 0) fs::fatal() << "Could not create temporary file in " << basedir;
  close(fid);
  return filename;
}


/*!
  /brief Generates a unique, temporary directory.
*/
std::string makeTempDir()
{
  std::string basedir = getBaseTempDir();
  std::string dirname = basedir + "/tmp.XXXXXX";
  char * res = mkdtemp(&dirname[0]);
  if (res == nullptr) fs::fatal() << "Could not create temporary directory in " << basedir;
  return dirname;
}


/*!
  \fn std::vector<unsigned int> FindMatches(std::vector<int> iv, int i)
  \brief return all indices in iv that match i. Returns a possibly
  empty unsigned int vector of indices. Note that this function can be
  overloaded to use other data types
 */
std::vector<unsigned int> FindMatches(std::vector<int> iv, int i)
{
  std::vector<unsigned int> indices;
  std::vector<int>::iterator itr = iv.begin();
  while(itr != iv.end()){
    if(*itr == i) indices.push_back(itr-iv.begin());
    itr++;
  }
  return(indices);
}
/*!
  \fn int FindFirstMatch(std::vector<int> iv, int i)
  \brief return the first index in iv that matchs i. If nothing
  matches, then returns -1. FindMatches() could be used but this is
  potentially faster when only one match is expected or only the first
  or any match is wanted.  Note that this function can be overloaded
  to use other data types
 */
int FindFirstMatch(std::vector<int> iv, int i)
{
  std::vector<int>::iterator itr = iv.begin();
  while(itr != iv.end()){
    if(*itr == i) return(itr-iv.begin());
    itr++;
  }
  return(-1);
}
