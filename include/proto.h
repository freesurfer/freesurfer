/**
 * @file  proto.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.34 $
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


#ifndef _PROTO_H_
#define _PROTO_H_

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sys/types.h>
#include <time.h>
#include <sys/timeb.h>

#if defined(SunOS) | defined(IRIX)
#include <ieeefp.h>
#endif

/*----------------- ALL PLATFORMS *--------------------*/
int stricmp(char *str1, char *str2) ;


/*----------- SunOS -----------------------*/
#ifdef SunOSX

/* should be in stdlib.h */
#ifndef EXIT_FAILURE
#define EXIT_FAILURE  1
#endif

/* needed for select */
#include <sys/types.h>
#include <sys/time.h>

int _flsbuf(unsigned char c, FILE *p) ;
int _filbuf(FILE *p) ;
int select (int width, fd_set *readfds, fd_set *writefds, fd_set *exceptfds,
            struct timeval *timeout) ;

int puts(char *s) ;
int fputs(char *s, FILE *stream) ;
int getw(FILE *stream);
int putw(int w, FILE *stream);
int fputc(int c, FILE *stream);
int fgetc(FILE *stream);
int pclose(FILE *stream);
int ftime(struct timeb *tp) ;
/* void ftime(struct timeb *tm) ;*/
int fflush(FILE *fp) ;
int scanf(const char *fmt, ...) ;
int fscanf(FILE *fp, const char *fmt, ...) ;
int fprintf(FILE *fp, const char *fmt, ...) ;
int printf(const char *fmt, ...) ;
int vsprintf(char *str, const char *fmt, va_list args) ;
int vfprintf(FILE *fp, const char *fmt, va_list args) ;
int vprintf(const char *fmt, va_list args) ;
int sscanf(char *str, const char *fmt, ...) ;
void system(char *command_string) ;
void perror(char *s) ;
int fgetc(FILE *fp) ;
time_t time(time_t *tloc) ;
void fclose(FILE *fp) ;
void rewind(FILE *fp) ;
char toupper(char c) ;
char tolower(char c) ;
int fseek(FILE *fp, long offset, int whence) ;
int fread(void *ptr, int size, int nitems, FILE *fp) ;
int fwrite(void *ptr, int size, int nitems, FILE *fp) ;
#endif   /* SunOS */

/*----------- Mac OS/X -----------------------*/
#ifdef Darwin
#include "utils.h"
void *malloc(size_t byteSize) ;
#define nint(f)  ((int) (floor(f + .5)))
/* ((int)(rint((double)f))) */
#if 0
double drand48(void);
#else
#define drand48()   randomNumber(0.0, 1.0)
#endif
void srand48(long seed);
/* #define log2(d)  log10(d)/log10(2) */
/* #define exp2(d)  pow(2.0, d) */
int ftime(struct timeb *tp) ;

#endif

/*----------- Linux -----------------------*/
#ifdef Linux
#include "macros.h"
#define nint(f)  ((int) (floor(f + .5)))
/* ((int)(rint((double)f))) */

#if 0
int getw(FILE *stream);
int putw(int w, FILE *stream);
void swab(const void *from, void *to, size_t n);
#endif

#endif   /* Linux */

#ifdef SunOS
int ftime(struct timeb *tp) ;
#define nint(f)   ((int) (floor(f + .5)))
/* ((in ((int)(rint((double)f))) */
#include <ieeefp.h>
/*void bzero(void *s, int n);*/

#endif

/*----------- IRIX (SGI) -------------------*/
#ifdef IRIX
#define nint(f)   ((int) (floor(f + .5)))
/*  ((int)(rint((double)f))) */
/*#define isnan(f)  0*/
double rint(double x) ;
#endif

/*----------- MSDOS -----------------------*/
#ifdef _MSDOS
#define nint(f)   ((int) (floor(f + .5)))
/* ((int)((double)f+(f < 0 ? -0.5 : 0.5))) */
#define isnan(f) 0
#define unlink _unlink
#define hypot  _hypot

#endif   /* MSDOS */

#endif   /* #ifndef PROTO_H */


