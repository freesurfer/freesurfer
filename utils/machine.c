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
/*---------------------------------------------------------
  Name: int Arch486(void)
  Determines whether the current architecture is 486-based
  or not, and returns 1 if it is and 0 otherwise.
  Notes: 
    1. endianness = 1 = 486/Linux/PC
    2. endianness = 0 = Sun/IRIX
    3. This is consistent with the endianness flag in the 
       bshort/bfloat format.
  --------------------------------------------------------*/
int Arch486(void)
{
  int endian;
  short tmp = 1;
  char *ctmp;
  
  ctmp = (char *)(&tmp);
  if(*(ctmp+1) == 1) endian = 0;
  else               endian = 1;
  return(endian);
}
/*---------------------------------------------------------
  Name: ByteSwapBuf()
  Reverses the byte order of each nBytesPerItem buffer
  in Buf. nitems is the number of nBytesPerItem-byte items 
  in Buf.
  ---------------------------------------------------------*/
int ByteSwapBuf(void *Buf, long int nItems, int nBytesPerItem)
{
  switch(nBytesPerItem){
  case 1: break;
  case 2: ByteSwap2(Buf, nItems); break;
  case 4: ByteSwap4(Buf, nItems); break;
  case 8: ByteSwap8(Buf, nItems); break;
  default:
    printf("ERROR: ByteSwapBuf: nBytesPerItem = %d not supported\n",
     nBytesPerItem);
    return(1);
  }
  return(0);
}
/*---------------------------------------------------------
  Name: ByteSwap2()
  Reverses the byte order of each 2-byte buffer of buf2.
  nitems is the number of 2-byte items in buf2;
  ---------------------------------------------------------*/
int ByteSwap2(void *buf2, long int nitems)
{
  register char *cbuf, ctmp;
  register long int n;

  cbuf = (char *) buf2;
  for(n=0; n < nitems ; n+=2){
    ctmp = *cbuf;
    *cbuf = *(cbuf+1);
    *(cbuf+1) = ctmp;
    cbuf += 2;
  }
  return(0);
}
/*---------------------------------------------------------
  Name: ByteSwap4()
  Reverses the byte order of each 4-byte buffer of buf4.
  nitems is the number of 4-byte items in buf4;
  ---------------------------------------------------------*/
int ByteSwap4(void *buf4, long int nitems)
{
  register char *cbuf, ctmp;
  register long int n;

  cbuf = (char *) buf4;
  for(n=0; n < nitems ; n+=4){

    /* swap the first and fourth */
    ctmp = *cbuf;
    *cbuf = *(cbuf+3);
    *(cbuf+3) = ctmp;

    /* swap the second and third */
    ctmp = *(cbuf+1);
    *(cbuf+1) = *(cbuf+2);
    *(cbuf+2) = ctmp;

    cbuf += 4;
  }

  return(0);
}
/*---------------------------------------------------------
  Name: ByteSwap8()
  Reverses the byte order of each 8-byte buffer of buf8.
  nitems is the number of 8-byte items in buf8;
  ---------------------------------------------------------*/
int ByteSwap8(void *buf8, long int nitems)
{
  register char *cbuf, ctmp;
  register long int n;

  cbuf = (char *) buf8;
  for(n=0; n < nitems ; n+=8){

    /* swap the first and eigth */
    ctmp = *cbuf;
    *cbuf = *(cbuf+7);
    *(cbuf+7) = ctmp;

    /* swap the second and seventh */
    ctmp = *(cbuf+1);
    *(cbuf+1) = *(cbuf+6);
    *(cbuf+6) = ctmp;

    /* swap the third and sixth */
    ctmp = *(cbuf+2);
    *(cbuf+2) = *(cbuf+5);
    *(cbuf+5) = ctmp;

    /* swap the fourth and fifth */
    ctmp = *(cbuf+3);
    *(cbuf+3) = *(cbuf+4);
    *(cbuf+4) = ctmp;

    cbuf += 8;
  }

  return(0);
}

