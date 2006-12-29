/**
 * @file  error.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.15 $
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
 *       FILE NAME:   error.h
 *
 *       DESCRIPTION: error handling prototypes
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
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
void    ErrorExit(int ecode, char *fmt, ...) ;
int     ErrorPrintf(int ecode, char *fmt, ...) ;
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


extern int Gerror ;    /* global error value */

#endif
