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
