/*
 *  ifilbuf -
 *
 *        Paul Haeberli - 1984
 *
 */
#include  "rgb_image.h"

int ifilbuf(RGB_IMAGE *image) {
  int size;

  if ((image->flags&_IOREAD) == 0)
    return(EOF);
  if (image->base==NULL) {
    size = IBUFSIZE(image->xsize);
    if ((image->base = ibufalloc(image)) == NULL) {
      i_errhdlr("can't alloc image buffer\n",0,0,0,0);
      return EOF;
    }
  }
  image->cnt = getrow(image,image->base,image->y,image->z);
  image->ptr = image->base;
  if (--image->cnt < 0) {
    if (image->cnt == -1) {
      image->flags |= _IOEOF;
      if (image->flags & _IORW)
        image->flags &= ~_IOREAD;
    } else
      image->flags |= _IOERR;
    image->cnt = 0;
    return -1;
  }
  if (++image->y >= image->ysize) {
    image->y = 0;
    if (++image->z >= image->zsize) {
      image->z = image->zsize-1;
      image->flags |= _IOEOF;
      return -1;
    }
  }
  return *image->ptr++ & 0xffff;
}
