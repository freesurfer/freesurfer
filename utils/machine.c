#include <stdio.h>

#include "machine.h"

typedef union
{
  long  l ;
  float f ;
  int   i ;
  char  buf[4] ;
  short s[2] ;
} SWAP_LONG ;

long
swapLong(long l)
{
  SWAP_LONG  sl ;
  short      s ;

  /* first swap bytes in each word */
  sl.l = l ;
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;

  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.l) ;
}

float
swapFloat(float f)
{
  SWAP_LONG  sl ;
  short      s ;

  /* first swap bytes in each word */
  sl.f = f ;
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;

  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.f) ;
}

typedef union
{
  short  s ;
  char   buf[sizeof(short)] ;
} SWAP_SHORT ;

short
swapShort(short s)
{
  SWAP_SHORT ss ;
  char       c ;

  /* first swap bytes in word */
  ss.s = s ;
  c = ss.buf[0] ;
  ss.buf[0] = ss.buf[1] ;
  ss.buf[1] = c ;

  return(ss.s) ;
}

typedef union
{
  double  d ;
  long    l[sizeof(double) / sizeof(long)] ;
} SWAP_DOUBLE ;

double
swapDouble(double d)
{
  SWAP_DOUBLE  sd ;
  long         l ;

  sd.d = d ;

  sd.l[0] = swapLong(sd.l[0]) ;
  sd.l[1] = swapLong(sd.l[1]) ;
  l = sd.l[0] ;
  sd.l[0] = sd.l[1] ;
  sd.l[1] = l ;

  return(sd.d) ;
}
int
swapInt(int i)
{
  SWAP_LONG  sl ;
  short      s ;

  /* first swap bytes in each word */
  sl.i = i ;
  sl.s[0] = swapShort(sl.s[0]) ;
  sl.s[1] = swapShort(sl.s[1]) ;

  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.i) ;
}

