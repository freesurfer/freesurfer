#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fio.h"
#include "machine.h"
#include "proto.h"

FILE *MGHopen_file(char *fname, char *rwmode)
{
  FILE *f1;

  if ((f1 = fopen(fname,rwmode)) == NULL)
  {
    printf("Can't open %s\n",fname);
    exit(1);
  }
  return f1;
}

int
putf(float f, FILE *fp)
{
#ifdef Linux
  f = swapFloat(f) ;
#endif
  return(fwrite(&f,4,1,fp));
}

float
getf(FILE *fp)
{
  float f;

  fread(&f,4,1,fp);
#ifdef Linux
  f = swapFloat(f) ;
#endif
  return f;
}

int
fread1(int *v, FILE *fp)
{
  unsigned char c;
  int  ret ;

  ret = fread(&c,1,1,fp);
  *v = c;
  return(ret) ;
}


int
fread2(int *v, FILE *fp)
{
  short s;
  int   ret ;

  ret = fread(&s,2,1,fp);
#ifdef Linux
  s = swapShort(s) ;
#endif
  *v = s;
  return(ret) ;
}

int
fread3(int *v, FILE *fp)
{
  unsigned int i = 0;
  int  ret ;

  ret = fread(&i,3,1,fp);
#ifdef Linux
  i = (unsigned int)swapInt(i) ;
#endif
  *v = ((i>>8) & 0xffffff);
  return(ret) ;
}

int
fread4(float *v, FILE *fp)
{
  float f;
  int   ret ;

  ret = fread(&f,4,1,fp);
#ifdef Linux
  f = swapFloat(f) ;
#endif
  *v = f;
  return(ret) ;
}

int
fwrite1(int v,FILE *fp)
{
  unsigned char c = (unsigned char)v;

  return(fwrite(&c,1,1,fp));
}

int
fwrite2(int v, FILE *fp)
{
  short s = (short)v;

#ifdef Linux
  s = swapShort(s) ;
#endif
  return(fwrite(&s,2,1,fp));
}

int
fwrite3(int v, FILE *fp)
{
  unsigned int i = (unsigned int)(v<<8);

#ifdef Linux
  i = (unsigned int)swapInt(i) ;
#endif
  return(fwrite(&i,3,1,fp));
}

int
fwrite4(int v,FILE *fp)
{
#ifdef Linux
  v = swapInt(v) ;
#endif
  return(fwrite(&v,4,1,fp));
}

float
freadFloat(FILE *fp)
{
  float f;
  int   ret ;

  ret = fread(&f,4,1,fp);
#ifdef Linux
  f = swapFloat(f) ;
#endif
  return(f) ;
}

int
freadInt(FILE *fp)
{
  int  i, nread ;

  nread = fread(&i,sizeof(int),1,fp);
#ifdef Linux
  i = swapInt(i) ;
#endif
  return(i) ;
}

short
freadShort(FILE *fp)
{
  int   nread ;
	short s ;

  nread = fread(&s,sizeof(short),1,fp);
#ifdef Linux
  s = swapShort(s) ;
#endif
  return(s) ;
}

int
fwriteInt(int v, FILE *fp)
{
#ifdef Linux
  v = swapInt(v) ;
#endif
  return(fwrite(&v,sizeof(int),1,fp));
}

int
fwriteFloat(float f, FILE *fp)
{
#ifdef Linux
  f = swapFloat(f) ;
#endif
  return(fwrite(&f,sizeof(float),1,fp));
}

