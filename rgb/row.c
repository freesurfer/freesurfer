/*
 *  putrow, getrow -
 *
 *        Paul Haeberli - 1984
 *
 */
#include  <stdio.h>
#include  "rgb_image.h"

/* twitzel added this function */

int putrow_uc(RGB_IMAGE *image, unsigned char *buffer,
              unsigned int y, unsigned int z) {
  register unsigned char   *sptr;
  register unsigned char      *cptr;
  register unsigned int x;
  register unsigned long min, max;
  register long cnt;

  if ( !(image->flags & (_IORW|_IOWRT)) )
    return -1;
  if (image->dim<3)
    z = 0;
  if (image->dim<2)
    y = 0;
  if (ISVERBATIM(image->type)) {
    switch (BPP(image->type)) {
    case 1:
      min = image->min;
      max = image->max;
      cptr = (unsigned char *)image->tmpbuf;
      sptr = buffer;
      /* this might be uneccessary for this uc function */
      for (x=image->xsize; x--;) {
        *cptr = *sptr++;
        if (*cptr > max) max = *cptr;
        if (*cptr < min) min = *cptr;
        cptr++;
      }
      image->min = min;
      image->max = max;
      img_seek(image,y,z);
      cnt = image->xsize;
      if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
        i_errhdlr("putrow: error on write of row\n",0,0,0,0);
        return -1;
      } else
        return cnt;
      /* NOTREACHED */

    case 2:
      printf("ERROR: this rgb save function BPP=2 is not implemented\n");
      /* min = image->min;
      max = image->max;
      sptr = buffer;
      for(x=image->xsize; x--;) {
      if (*sptr > max) max = *sptr;
      if (*sptr < min) min = *sptr;
      sptr++;
      }
      image->min = min;
      image->max = max;
      img_seek(image,y,z);
      cnt = image->xsize<<1;
      if (img_write(image,(char *)(buffer),cnt) != cnt) {
      i_errhdlr("putrow: error on write of row\n",0,0,0,0);
      return -1;
      } else {
      return image->xsize;
      }
      */
      /* NOTREACHED */

    default:
      i_errhdlr("putrow: weird bpp\n",0,0,0,0);
    }
  } else if (ISRLE(image->type)) {
    printf("ERROR: compressed RGB is not implemented !\n");
    /*
    switch(BPP(image->type)) {
    case 1:
      min = image->min;
      max = image->max;
      sptr = buffer;
      for(x=image->xsize; x--;) {
    if (*sptr > max) max = *sptr;
    if (*sptr < min) min = *sptr;
    sptr++;
      }
      image->min = min;
      image->max = max;
      cnt = img_rle_compact(buffer,2,image->tmpbuf,1,image->xsize);
      img_setrowsize(image,cnt,y,z);
      img_seek(image,y,z);
      if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
    i_errhdlr("putrow: error on write of row\n",0,0,0,0);
    return -1;
      } else
    return image->xsize;
      */
    /* NOTREACHED */
    /*
    case 2:
    min = image->min;
    max = image->max;
    sptr = buffer;
    for(x=image->xsize; x--;) {
    if (*sptr > max) max = *sptr;
    if (*sptr < min) min = *sptr;
    sptr++;
    }
    image->min = min;
    image->max = max;
    cnt = img_rle_compact(buffer,2,image->tmpbuf,2,image->xsize);
    cnt <<= 1;
    img_setrowsize(image,cnt,y,z);
    img_seek(image,y,z);
    if(image->dorev)
    cvtshorts(image->tmpbuf,cnt);
    if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
    if(image->dorev)
    cvtshorts(image->tmpbuf,cnt);
    i_errhdlr("putrow: error on write of row\n",0,0,0,0);
    return -1;
    } else {
    if(image->dorev)
    cvtshorts(image->tmpbuf,cnt);
    return image->xsize;
    }
    */
    /* NOTREACHED */
    /*
    default:
    i_errhdlr("putrow: weird bpp\n",0,0,0,0);
    } */
  } else
    i_errhdlr("putrow: weird image type\n",0,0,0,0);
  return(-1);
}

int putrow(RGB_IMAGE *image, unsigned short *buffer,
           unsigned int y, unsigned int z) {
  register unsigned short   *sptr;
  register unsigned char      *cptr;
  register unsigned int x;
  register unsigned long min, max;
  register long cnt;

  if ( !(image->flags & (_IORW|_IOWRT)) )
    return -1;
  if (image->dim<3)
    z = 0;
  if (image->dim<2)
    y = 0;
  if (ISVERBATIM(image->type)) {
    switch (BPP(image->type)) {
    case 1:
      min = image->min;
      max = image->max;
      cptr = (unsigned char *)image->tmpbuf;
      sptr = buffer;
      for (x=image->xsize; x--;) {
        *cptr = *sptr++;
        if (*cptr > max) max = *cptr;
        if (*cptr < min) min = *cptr;
        cptr++;
      }
      image->min = min;
      image->max = max;
      img_seek(image,y,z);
      cnt = image->xsize;
      if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
        i_errhdlr("putrow: error on write of row\n",0,0,0,0);
        return -1;
      } else
        return cnt;
      /* NOTREACHED */

    case 2:
      min = image->min;
      max = image->max;
      sptr = buffer;
      for (x=image->xsize; x--;) {
        if (*sptr > max) max = *sptr;
        if (*sptr < min) min = *sptr;
        sptr++;
      }
      image->min = min;
      image->max = max;
      img_seek(image,y,z);
      cnt = image->xsize<<1;
      if (image->dorev)
        cvtshorts(buffer,cnt);
      if (img_write(image,(char *)(buffer),cnt) != cnt) {
        if (image->dorev)
          cvtshorts(buffer,cnt);
        i_errhdlr("putrow: error on write of row\n",0,0,0,0);
        return -1;
      } else {
        if (image->dorev)
          cvtshorts(buffer,cnt);
        return image->xsize;
      }
      /* NOTREACHED */

    default:
      i_errhdlr("putrow: weird bpp\n",0,0,0,0);
    }
  } else if (ISRLE(image->type)) {
    switch (BPP(image->type)) {
    case 1:
      min = image->min;
      max = image->max;
      sptr = buffer;
      for (x=image->xsize; x--;) {
        if (*sptr > max) max = *sptr;
        if (*sptr < min) min = *sptr;
        sptr++;
      }
      image->min = min;
      image->max = max;
      cnt = img_rle_compact(buffer,2,image->tmpbuf,1,image->xsize);
      img_setrowsize(image,cnt,y,z);
      img_seek(image,y,z);
      if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
        i_errhdlr("putrow: error on write of row\n",0,0,0,0);
        return -1;
      } else
        return image->xsize;
      /* NOTREACHED */

    case 2:
      min = image->min;
      max = image->max;
      sptr = buffer;
      for (x=image->xsize; x--;) {
        if (*sptr > max) max = *sptr;
        if (*sptr < min) min = *sptr;
        sptr++;
      }
      image->min = min;
      image->max = max;
      cnt = img_rle_compact(buffer,2,image->tmpbuf,2,image->xsize);
      cnt <<= 1;
      img_setrowsize(image,cnt,y,z);
      img_seek(image,y,z);
      if (image->dorev)
        cvtshorts(image->tmpbuf,cnt);
      if (img_write(image,(char *)(image->tmpbuf),cnt) != cnt) {
        if (image->dorev)
          cvtshorts(image->tmpbuf,cnt);
        i_errhdlr("putrow: error on write of row\n",0,0,0,0);
        return -1;
      } else {
        if (image->dorev)
          cvtshorts(image->tmpbuf,cnt);
        return image->xsize;
      }
      /* NOTREACHED */

    default:
      i_errhdlr("putrow: weird bpp\n",0,0,0,0);
    }
  } else
    i_errhdlr("putrow: weird image type\n",0,0,0,0);
  return(-1);
}

int getrow(RGB_IMAGE *image, unsigned short *buffer,
           unsigned int y, unsigned int z) {
  register short i;
  register unsigned char *cptr;
  register unsigned short *sptr;
  register short cnt;

  if ( !(image->flags & (_IORW|_IOREAD)) )
    return -1;
  if (image->dim<3)
    z = 0;
  if (image->dim<2)
    y = 0;
  img_seek(image, y, z);
  if (ISVERBATIM(image->type)) {
    switch (BPP(image->type)) {
    case 1:
      if (img_read(image,(char *)image->tmpbuf,image->xsize)
          != image->xsize) {
        i_errhdlr("getrow: error on read of row\n",0,0,0,0);
        return -1;
      } else {
        cptr = (unsigned char *)image->tmpbuf;
        sptr = buffer;
        for (i=image->xsize; i--;)
          *sptr++ = *cptr++;
      }
      return image->xsize;
      /* NOTREACHED */

    case 2:
      cnt = image->xsize<<1;
      if (img_read(image,(char *)(buffer),cnt) != cnt) {
        i_errhdlr("getrow: error on read of row\n",0,0,0,0);
        return -1;
      } else {
        if (image->dorev)
          cvtshorts(buffer,cnt);
        return image->xsize;
      }
      /* NOTREACHED */

    default:
      i_errhdlr("getrow: weird bpp\n",0,0,0,0);
      break;
    }
  } else if (ISRLE(image->type)) {
    switch (BPP(image->type)) {
    case 1:
      if ( (cnt = img_getrowsize(image)) == -1 )
        return -1;
      if ( img_read(image,(char *)(image->tmpbuf),cnt) != cnt ) {
        i_errhdlr("getrow: error on read of row\n",0,0,0,0);
        return -1;
      } else {
        img_rle_expand(image->tmpbuf,1,buffer,2);
        return image->xsize;
      }
      /* NOTREACHED */

    case 2:
      if ( (cnt = img_getrowsize(image)) == -1 )
        return -1;
      if ( cnt != img_read(image,(char *)(image->tmpbuf),cnt) ) {
        i_errhdlr("getrow: error on read of row\n",0,0,0,0);
        return -1;
      } else {
        if (image->dorev)
          cvtshorts(image->tmpbuf,cnt);
        img_rle_expand(image->tmpbuf,2,buffer,2);
        return image->xsize;
      }
      /* NOTREACHED */

    default:
      i_errhdlr("getrow: weird bpp\n",0,0,0,0);
      break;
    }
  } else
    i_errhdlr("getrow: weird image type\n",0,0,0,0);
  return(-1) ;
}
