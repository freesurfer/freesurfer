#ifndef FIO_H
#define FIO_H

FILE  *MGHopen_file(char *fname, char *rwmode) ;
int   putf(float f, FILE *fp) ;
float getf(FILE *fp) ;

int   fread1(int *v, FILE *fp) ;
int   fread2(int *v, FILE *fp) ;
int   fread3(int *v, FILE *fp) ;
int   fread4(float *v, FILE *fp) ;
float freadFloat(FILE *fp) ;
int   freadInt(FILE *fp) ;
short freadShort(FILE *fp) ;

int   fwriteFloat(float f, FILE *fp) ;
int   fwriteShort(short s, FILE *fp) ;
int   fwriteInt(int v, FILE *fp) ;
int   fwrite1(int v,FILE *fp) ;
int   fwrite2(int v, FILE *fp) ;
int   fwrite3(int v, FILE *fp) ;
int   fwrite4(int v, FILE *fp) ;


#endif
