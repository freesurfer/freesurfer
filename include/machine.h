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

#ifdef Linux

#define orderIntBytes(i)     swapInt(i)
#define orderShortBytes(i)   swapShort(i)
#define orderFloatBytes(i)   swapFloat(i)
#define orderDoubleBytes(i)  swapDouble(i)
#define orderLongBytes(i)    swapLong(i)

#else

#define orderIntBytes(i)     (i)
#define orderShortBytes(i)   (i)
#define orderFloatBytes(i)   (i)
#define orderDoubleBytes(i)  (i)
#define orderLongBytes(i)    (i)


#endif

#endif
