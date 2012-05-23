/**
 * @file  error.c
 * @brief error handling routines
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/05/23 17:36:46 $
 *    $Revision: 1.23 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include <hips_error.h>

#include "error.h"
#include "hips.h"
#include "proto.h"
#include "rgb_image.h"
#include "fsinit.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define ERROR_FNAME  "error.log"

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static char error_fname[100] = ERROR_FNAME ;
static int (*error_vprintf)(const char *fmt, va_list args) = vprintf ;
static int (*error_vfprintf)(FILE *fp,const char *fmt,va_list args) = vfprintf;
static void (*error_exit)(int ecode) = NULL ;
static void rgb_error(char *error_str) ;

/*-----------------------------------------------------
                      GLOBAL DATA
-------------------------------------------------------*/

int Gerror = NO_ERROR ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void
rgb_error(char *error_str)
{
  ErrorPrintf(ERROR_BADPARM, error_str) ;
  return ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ErrorInit(char *fname,
          int (*vfprint)(FILE *fp, const char *fmt, va_list args),
          int (*vprint)(const char *fmt, va_list args))
{
  FSinit() ;
  error_exit = (void *)exit;
  i_seterror(rgb_error) ;
  if (fname)
    strcpy(error_fname, fname) ;
  if (vfprint)
    error_vfprintf = vfprint ;
  if (vprint)
    error_vprintf = vprint ;

  unlink(error_fname) ; /* start with a fresh log file */
  errno = 0 ;

  /* probably should be some info into log file like date/user etc... */
  return(NO_ERROR) ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
void
ErrorExit(int ecode, const char *fmt, ...)
{
  va_list  args ;

  Gerror = ecode ;
  va_start(args, fmt) ;
  vfprintf(stderr, fmt, args) ;
  fprintf(stderr, "\n") ;
  fflush(stderr);
  fflush(stdout);
  if (errno)
    perror(NULL) ;
  if (hipserrno)
    perr(ecode, "Hips error:") ;

  if (error_exit)
    (*error_exit)(ecode) ;
  else
    exit(ecode) ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ErrorPrintf(int ecode, const char *fmt, ...)
{
  va_list  args ;
  FILE     *fp ;

  Gerror = ecode ;
  va_start(args, fmt) ;
  (*error_vfprintf)(stderr, fmt, args) ;
  fprintf(stderr, "\n") ;
  fflush(stderr);
  fflush(stdout);
  va_end(args);
  if (errno)
    perror(NULL) ;
#if 0
  if (hipserrno)
    perr(ecode, "Hips error:") ;
#endif

  va_start(args, fmt) ;
  fp = fopen(ERROR_FNAME, "a") ;
  if (fp)
  {
    (*error_vfprintf)(fp, fmt, args) ;
    fprintf(fp, "\n") ;
    fclose(fp) ;     /* close file to flush changes */
  }
  va_end(args);
  return(ecode) ;
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ErrorSetExitFunc(void (*exit_func)(int ecode))
{
  error_exit = exit_func ;
  return(1) ;
}

