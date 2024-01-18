#include <string.h>
#include <stdlib.h>

#include "fsbufferio.h"
#include "machine.h"
#include "bfileio.h"
#include "tags.h"
#include "mri.h"


FSbufferIO::FSbufferIO(int size)
{
  fsbufsize = size;
  fsbuflen  = 0;

  fsbuf = (unsigned char*)malloc(fsbufsize);
  memset(fsbuf, 0, fsbufsize);
}


FSbufferIO::~FSbufferIO()
{
  free(fsbuf);
  
  fsbuf = NULL;
  fsbufsize = 0; fsbuflen = 0;
}


int FSbufferIO::write_buf(void *buf, int buflen)
{
  _checkbufsize(buflen);  
  memcpy(fsbuf+fsbuflen, buf, buflen);
  fsbuflen += buflen;

  return buflen;
}


// based on znzWriteInt()
int FSbufferIO::write_int(int v)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapInt(v);
#endif
    
  int dlen = sizeof(int);
  return write_buf(&v, dlen);
#if 0  
  _checkbufsize(dlen);  
  memcpy(fsbuf+fsbuflen, &v, dlen);
  fsbuflen += dlen;

  return dlen;
#endif  
}

// based on znzWriteLong()
int FSbufferIO::write_long(long long v)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  v = swapLong64(v);
#endif

  int dlen = sizeof(long long);
  return write_buf(&v, dlen);
#if 0
  _checkbufsize(dlen);    
  memcpy(fsbuf+fsbuflen, &v, dlen);
  fsbuflen += dlen;

  return dlen;
#endif
}

// based on znzWriteFloat()
int FSbufferIO::write_float(float f)
{
  char buf[4];
  memmove(buf, &f, 4);
  
#if (BYTE_ORDER == LITTLE_ENDIAN)
  byteswapbuffloat(buf, 1);
// f = swapFloat(f);  // old way
#endif

  int dlen = sizeof(float);
  return write_buf(buf, dlen);
#if 0  
  _checkbufsize(dlen);    
  memcpy(fsbuf+fsbuflen, buf, dlen);
  fsbuflen += dlen;

  return dlen;
#endif  
}

// based on znzWriteDouble()
int FSbufferIO::write_double(double d)
{
#if (BYTE_ORDER == LITTLE_ENDIAN)
  d = swapDouble(d);
#endif

  int dlen = sizeof(double);
  return write_buf(&d, dlen);
#if 0
  _checkbufsize(dlen);
  memcpy(fsbuf+fsbuflen, &d, dlen);
  fsbuflen += dlen;

  return dlen;
#endif  
}


int FSbufferIO::write_matrix(MATRIX *M)
{
  long long dlen = MATRIX_STRLEN;
  char matbuf[dlen];

  bzero(matbuf, dlen);
  sprintf(matbuf,
          "Matrix %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf %10lf",
          M->rptr[1][1],
          M->rptr[1][2],
          M->rptr[1][3],
          M->rptr[1][4],
          M->rptr[2][1],
          M->rptr[2][2],
          M->rptr[2][3],
          M->rptr[2][4],
          M->rptr[3][1],
          M->rptr[3][2],
          M->rptr[3][3],
          M->rptr[3][4],
          M->rptr[4][1],
          M->rptr[4][2],
          M->rptr[4][3],
          M->rptr[4][4]);

  return write_buf(matbuf, dlen);
}


int FSbufferIO::write_geom(VOL_GEOM *volgeom)
{
  int bytewritten = 0;
  bytewritten += write_int(volgeom->valid);
    
  bytewritten += write_int(volgeom->width);
  bytewritten += write_int(volgeom->height);    
  bytewritten += write_int(volgeom->depth);
    
  bytewritten += write_float(volgeom->xsize);
  bytewritten += write_float(volgeom->ysize);
  bytewritten += write_float(volgeom->zsize);
    
  bytewritten += write_float(volgeom->x_r);
  bytewritten += write_float(volgeom->x_a);
  bytewritten += write_float(volgeom->x_s);
    
  bytewritten += write_float(volgeom->y_r);
  bytewritten += write_float(volgeom->y_a);    
  bytewritten += write_float(volgeom->y_s);
    
  bytewritten += write_float(volgeom->z_r);
  bytewritten += write_float(volgeom->z_a);    
  bytewritten += write_float(volgeom->z_s);
    
  bytewritten += write_float(volgeom->c_r);
  bytewritten += write_float(volgeom->c_a);
  bytewritten += write_float(volgeom->c_s);

  bytewritten += write_buf(volgeom->fname, strlen(volgeom->fname));

  return bytewritten;
}


void FSbufferIO::_checkbufsize(int buflen)
{
  if (fsbuflen + buflen <= fsbufsize)
    return;
  
  while (fsbuflen + buflen > fsbufsize)
    fsbufsize *= 2;

  unsigned char *tmpbuf = (unsigned char*)malloc(fsbufsize);
    
  memcpy(tmpbuf, fsbuf, fsbuflen);
  free(fsbuf);
  fsbuf = tmpbuf;
}
