/**
 * @file  error.c
 * @brief error handling routines
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2015/07/28 21:24:19 $
 *    $Revision: 1.25 $
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

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "error.h"
#include "fsinit.h"
#include "hips.h"
#include "proto.h"
#include "rgb.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

#define ERROR_FNAME "error.log"

/*-----------------------------------------------------
                    STATIC DATA
-------------------------------------------------------*/

static char error_fname[100] = ERROR_FNAME;
static int (*error_vprintf)(const char *fmt, va_list args) = vprintf;
static int (*error_vfprintf)(FILE *fp, const char *fmt, va_list args) = vfprintf;
static void (*error_exit)(int ecode) = NULL;
static void rgb_error(char *error_str);
char *ErrorExitDoneFile = NULL;

/*-----------------------------------------------------
                      GLOBAL DATA
-------------------------------------------------------*/

int Gerror = NO_ERROR;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/

void SetErrorExitDoneFile(char *DoneFile)
{
  extern char *ErrorExitDoneFile;
  ErrorExitDoneFile = DoneFile;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static void rgb_error(char *error_str)
{
  ErrorPrintf(ERROR_BADPARM, error_str);
  return;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ErrorInit(char *fname,
              int (*vfprint)(FILE *fp, const char *fmt, va_list args),
              int (*vprint)(const char *fmt, va_list args))
{
  FSinit();
  error_exit = exit;
  i_seterror(rgb_error);
  if (fname) strcpy(error_fname, fname);
  if (vfprint) error_vfprintf = vfprint;
  if (vprint) error_vprintf = vprint;

  unlink(error_fname); /* start with a fresh log file */
  errno = 0;

  /* probably should be some info into log file like date/user etc... */
  return (NO_ERROR);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
void ErrorExit(int ecode, const char *fmt, ...)
{
  va_list args;
  extern char *ErrorExitDoneFile;

  Gerror = ecode;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  fflush(stderr);
  fflush(stdout);
  if (errno) perror(NULL);

  if (ErrorExitDoneFile != NULL) {
    // This creates a text file with the value of 1. This can
    // be used to let another process know that this process
    // is finshed and exited with an error. Good when submitting
    // to cluster
    ErrorWriteDoneFile(ErrorExitDoneFile, 1);
  }

  if (error_exit)
    (*error_exit)(ecode);
  else
    exit(ecode);
}

/*!
  \fn int ErrorWriteDoneFile(char *DoneFile, int errorcode)
  \brief This creates a text file with the contents being
  errorcode. This can be used to let another process know that this
  process is finshed and exited with or without an error. Good when
  submitting to cluster. The basic idea is to use
  SetErrorExitDoneFile(char *DoneFile) to set the file for
  ErrorExit(). If no error occurs during processing, then
  the process should run ErrorWriteDoneFile(DoneFile, 0);
 */
int ErrorWriteDoneFile(char *DoneFile, int errorcode)
{
  FILE *fp;
  if (DoneFile == NULL) return (0);
  fp = fopen(DoneFile, "w");
  if (fp == NULL) return (1);
  fprintf(fp, "%d\n", errorcode);
  fclose(fp);
  return (0);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ErrorPrintf(int ecode, const char *fmt, ...)
{
  va_list args;
  FILE *fp;

  Gerror = ecode;
  va_start(args, fmt);
  (*error_vfprintf)(stderr, fmt, args);
  fprintf(stderr, "\n");
  fflush(stderr);
  fflush(stdout);
  va_end(args);
  if (errno) perror(NULL);

  va_start(args, fmt);
  fp = fopen(ERROR_FNAME, "a");
  if (fp) {
    (*error_vfprintf)(fp, fmt, args);
    fprintf(fp, "\n");
    fclose(fp); /* close file to flush changes */
  }
  va_end(args);
  return (ecode);
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int ErrorSetExitFunc(void (*exit_func)(int ecode))
{
  error_exit = exit_func;
  return (1);
}
