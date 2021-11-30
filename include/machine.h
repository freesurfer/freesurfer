/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


#ifndef MACHINE_H
#define MACHINE_H

#include "mghendian.h"

/* this is the unix version */
#define huge
#define far
#define large
#define near
#define hmemmove(dst,src,n)  memcpy(dst,src,n)

// 32 bit OS long is 32bit.  64 bit OS long is 64 bit
// but int is 32 bit on both
typedef int long32;
typedef long long long64;

unsigned short swapUShort(unsigned short us);
short  swapShort(short s) ;
long32 swapLong32(long32 l);
long64 swapLong64(long64 l);
float  swapFloat(float l) ;
double swapDouble(double dval) ;
int    swapInt(int i) ;

int Arch486(void);
int ByteSwapBuf(void *Buf, long int nItems, int nBytesPerItem);
int ByteSwap2(void *buf2, long int nitems);
int ByteSwap4(void *buf4, long int nitems);
int ByteSwap8(void *buf8, long int nitems);

#if (BYTE_ORDER == LITTLE_ENDIAN)

#define orderIntBytes(i)     swapInt(i)
#define orderShortBytes(i)   swapShort(i)
#define orderUShortBytes(i)  swapUShort(i)
#define orderFloatBytes(i)   swapFloat(i)
#define orderDoubleBytes(i)  swapDouble(i)
#define orderLong32Bytes(i)  swapLong32(i)
#define orderLong64Bytes(i)  swapLong64(i)

#else

#define orderIntBytes(i)     (i)
#define orderShortBytes(i)   (i)
#define orderFloatBytes(i)   (i)
#define orderDoubleBytes(i)  (i)
#define orderLong32Bytes(i)  (i)
#define orderLong64Bytes(i)  (i)

#endif

#endif
