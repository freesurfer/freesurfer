/*
 *       FILE NAME:   error.c
 *
 *       DESCRIPTION: error handling routines
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
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

#include <hips_error.h>

#include "error.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
void
ErrorExit(int ecode, char *fmt, ...)
{
  va_list  args ;

  va_start(args, fmt) ;
  vfprintf(stderr, fmt, args) ;
  if (errno)
    perror(NULL) ;
  if (hipserrno)
    perr(ecode, "Hips error:") ;

  exit(ecode) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
ErrorPrintf(int ecode, char *fmt, ...)
{
  va_list  args ;

  va_start(args, fmt) ;
  vfprintf(stderr, fmt, args) ;
  if (errno)
    perror(NULL) ;
  if (hipserrno)
    perr(ecode, "Hips error:") ;

  return(ecode) ;
}

