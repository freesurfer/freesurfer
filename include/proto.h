#ifndef PROTO_H
#define PROTO_H

#include <stdio.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/timeb.h>

#ifdef SPARC

int fflush(FILE *fp) ;
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
int stricmp(char *str1, char *str2) ;
void ftime(struct timeb *tm) ;

#else
#ifdef LINUX

#include "macros.h"

int stricmp(char *str1, char *str2) ;
#define nint(f)   ((int)(rint((double)f)))

#endif   /* LINUX */
#endif   /* SPARC */

#endif

