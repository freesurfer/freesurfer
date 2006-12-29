/*
 *  iclose and iflush -
 *
 *        Paul Haeberli - 1984
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include  "rgb_image.h"

int iclose(RGB_IMAGE *image) {
  long tablesize, ret;

  iflush(image);
  img_optseek(image, 0);
  if (image->flags&_IOWRT) {
    if (image->dorev)
      cvtimage((long *)image);
    swapImage(image) ;
    if (img_write(image,(char *)image,sizeof(RGB_IMAGE)) != sizeof(RGB_IMAGE)) {
      i_errhdlr("iclose: error on write of image header\n",0,0,0,0);
      return EOF;
    }
    swapImage(image) ;
    if (image->dorev)
      cvtimage((long *)image);
    if (ISRLE(image->type)) {
      img_optseek(image, 512L);
      tablesize = image->ysize*image->zsize*sizeof(long);
      if (image->dorev)
        cvtlongs((long *)image->rowstart,(long)tablesize);
      if (img_write(image,(char *)(image->rowstart),tablesize) != tablesize) {
        i_errhdlr("iclose: error on write of rowstart\n",0,0,0,0);
        return EOF;
      }
      if (image->dorev)
        cvtlongs((long *)image->rowsize,tablesize);
      if (img_write(image,(char *)(image->rowsize),tablesize) != tablesize) {
        i_errhdlr("iclose: error on write of rowsize\n",0,0,0,0);
        return EOF;
      }
    }
  }
  if (image->base) {
    free(image->base);
    image->base = 0;
  }
  if (image->tmpbuf) {
    free(image->tmpbuf);
    image->tmpbuf = 0;
  }
  if (ISRLE(image->type)) {
    free(image->rowstart);
    image->rowstart = 0;
    free(image->rowsize);
    image->rowsize = 0;
  }
  ret = close(image->file);
  if (ret != 0)
    i_errhdlr("iclose: error on close of file\n",0,0,0,0);
  free(image);
  return ret;
}

int iflush(RGB_IMAGE *image) {
  unsigned short *base;

  if ( (image->flags&_IOWRT)
       && (base=image->base)!=NULL && (image->ptr-base)>0) {
    if (putrow(image, base, image->y,image->z)!=image->xsize) {
      image->flags |= _IOERR;
      return(EOF);
    }
  }
  return(0);
}
