/*
   @(#)utils.c  1.12
   10/16/95
*/
/*------------------------------------------------------------------------
      File Name: utils.c

         Author: Bruce Fischl

        Created: Jan. 1994

    Description: miscellaneous utility functions

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

#include "nr.h"
#include "const.h"
#include "utils.h"
#include "proto.h"
#include "error.h"
#include "image.h"
#include "macros.h"

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


double drand48(void) ;

/*------------------------------------------------------------------------
                              FUNCTIONS
------------------------------------------------------------------------*/
/*------------------------------------------------------------------------
       Parameters:

      Description:
  
    Return Values:
        nothing.
------------------------------------------------------------------------*/
double
randomNumber(double low, double hi)
{
  double val, range ;
  static long idum = -1L ;

  if (low > hi)
  {
    val = low ;
    low = hi ;
    hi = val ;
  }

  if (idum <= 0L)     /* change seed from run to run */
    idum = -1L * (long)(abs((int)time(NULL))) ; 

  range = hi - low ;
  val = ran1(&idum) * range + low ;
  
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

#ifdef _HOME_
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

    while (isspace(*cp))
      cp++ ;

  } while (((*cp) == '#') || ((*cp) == '\n') || ((*cp) == 0)) ;

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

  } while (approx != prev);

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
FileName(char *full_name, char *fname)
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

  fp = fopen(fname, "r") ;
  if (fp)
    fclose(fp) ;

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
FileNameNoExtensions(char *full_name)
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
  char *dot, buf[100], *number ;

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
    else if (!strcmp(dot, "HIPL")||!strcmp(dot, "HIPS") || !strcmp(dot,"HIP"))
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
  char buf[100], *number ;
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
  char buf[100], line[200], *cp ;

  strcpy(buf, fname) ;   /* we will modify fname, don't ruin callers copy */
  fname = buf ;
  
  num = FileNumber(fname) ;
  if (num == -1)
  {
    type = FileType(fname) ;
    switch (type)
    {
    case LIST_FILE:
      fp = fopen(FileNameNoExtensions(fname), "rb") ;
      if (!fp)
        ErrorReturn(-1, (ERROR_NO_FILE, 
                         "FileNumberOfEntries: could not open %s",
                         FileNameNoExtensions(fname))) ;
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
  static char fname[100] ;
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
  } while (fp) ;

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
#ifndef SunOS
int
stricmp(char *str1, char *str2)
{
  char buf1[100], buf2[100] ;

  strcpy(buf1, str1) ;
  strcpy(buf2, str2) ;
  StrUpper(buf1) ;
  StrUpper(buf2) ;
  return(strcmp(buf1, buf2)) ;
}
#endif
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
  while (isspace(*str))
    str++ ;

  return(str) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
------------------------------------------------------------------------*/
#ifdef SunOS
extern char *getwd(char *pathname) ;
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
    getwd(pathname) ;
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
    getwd(pathName)  ;    /* no path at all, must be cwd */

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
  while (*str && isdigit(*str))
    str++ ;
  if (*str == '.')  /* check for floating point # */
  {
    str++ ;
    while (*str && isdigit(*str))
      str++ ;
  }
  while (*str && isspace(*str))
    str++ ;

  return(str) ;
}

