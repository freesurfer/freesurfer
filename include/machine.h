#ifndef MACHINE_H
#define MACHINE_H

/* this is the unix version */
#define huge
#define far
#define large
#define near
#define hmemmove(dst,src,n)  memcpy(dst,src,n)

short  swapShort(short s) ;
long   swapLong(long l) ;
float   swapFloat(float l) ;
double swapDouble(double dval) ;
int    swapInt(int i) ;

#endif
