/**
 * @file  error.h
 * @brief error handling prototypes
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2015/07/27 20:49:35 $
 *    $Revision: 1.21 $
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


#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <stdarg.h>

int     ErrorInit(char *fname,
                  int (*vfprint)(FILE *fp, const char *fmt, va_list args),
                  int (*vprint)(const char *fmt, va_list args)) ;
int     ErrorSetExitFunc(void (*exit_func)(int ecode)) ;
void    ErrorExit(int ecode, const char *fmt, ...) ;
int     ErrorPrintf(int ecode, const char *fmt, ...) ;
void    SetErrorExitDoneFile(char *DoneFile);
int ErrorWriteDoneFile(char *DoneFile, int errorcode);
#define ErrorReturn(ret, args)  { ErrorPrintf args ; return(ret) ; }

#define ESCAPE   ErrorExit
#define ErrorSet ErrorPrintf
/* error codes */

#define NO_ERROR              0
#define ERROR_NONE            NO_ERROR
#define ERROR_NO_FILE         -1
#define ERROR_NOFILE          ERROR_NO_FILE
#define ERROR_NO_MEMORY       -2
#define ERROR_NOMEMORY        ERROR_NO_MEMORY
#define ERROR_UNSUPPORTED     -3
#define ERROR_BADPARM         -4
#define ERROR_BAD_PARM        ERROR_BADPARM
#define ERROR_BADFILE         -5
#define ERROR_BAD_FILE        ERROR_BADFILE
#define ERROR_SIZE            -6
#define ERROR_BADLOOP         -7
#define ERROR_OUT_OF_BOUNDS   -8

extern int Gerror ;    /* global error value */

#endif
