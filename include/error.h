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
#define ERROR_NO_FILE         -1
#define ERROR_NO_MEMORY       -2
#define ERROR_UNSUPPORTED     -3
#define ERROR_BADPARM         -4
#define ERROR_BADFILE         -5

#endif
