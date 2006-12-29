/*
 *  iflsbuf -
 *
 *        Paul Haeberli - 1984
 *
 */
#include  "rgb_image.h"

unsigned int iflsbuf(RGB_IMAGE *image, unsigned int c) {
  register unsigned short *base;
  register int n, rn;
  int size;

  if ((image->flags&_IOWRT)==0)
    return(EOF);
  if ((base=image->base)==NULL) {
    size = IBUFSIZE(image->xsize);
    if ((image->base=base=ibufalloc(image)) == NULL) {
      i_errhdlr("flsbuf: error on buf alloc\n",0,0,0,0);
      return EOF;
    }
    rn = n = 0;
  } else if ((rn = n = image->ptr - base) > 0)  {
    n = putrow(image,base,image->y,image->z);
    if (++image->y >= image->ysize) {
      image->y = 0;
      if (++image->z >= image->zsize) {
        image->z = image->zsize-1;
        image->flags |= _IOEOF;
        return -1;
      }
    }
  }
  image->cnt = image->xsize-1;
  *base++ = c;
  image->ptr = base;
  if (rn != n) {
    image->flags |= _IOERR;
    return(EOF);
  }
  return(c);
}
