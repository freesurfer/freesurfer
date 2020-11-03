#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

#include "rgb.h"

#define OPEN_GL_CODE  1

void isetname(RGB_IMAGE *image, const char *name) {
  strncpy(image->name,name,80-1);
}

void isetcolormap(RGB_IMAGE *image, int colormap) {
  image->colormap = colormap;
}

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

unsigned int iflsbuf(RGB_IMAGE *image, unsigned int c) {
  unsigned short *base;
  int n, rn;
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

int getpix(RGB_IMAGE *image) {
  if (--(image)->cnt>=0)
    return (int)(*(image)->ptr++);
  else
    return ifilbuf(image);
}

unsigned int putpix(RGB_IMAGE *image, unsigned int pix) {
  if (--(image)->cnt>=0)
    return (unsigned int)(*(image)->ptr++ = pix);
  else
    return iflsbuf(image,pix);
}

unsigned long img_seek(RGB_IMAGE *image, unsigned int y, unsigned int z) {
  if (img_badrow(image,y,z)) {
    i_errhdlr("img_seek: row number out of range\n", 0,0,0,0);
    return EOF;
  }
  image->x = 0;
  image->y = y;
  image->z = z;
  if (ISVERBATIM(image->type)) {
    switch (image->dim) {
    case 1:
      return img_optseek(image, 512L);
    case 2:
      return img_optseek(image,512L+(y*image->xsize)*BPP(image->type));
    case 3:
      return img_optseek(image,
                         512L+(y*image->xsize+z*image->xsize*image->ysize)*
                         BPP(image->type));
    default:
      i_errhdlr("img_seek: weird dim\n",0,0,0,0);
      break;
    }
  } else if (ISRLE(image->type)) {
    switch (image->dim) {
    case 1:
      return img_optseek(image, image->rowstart[0]);
    case 2:
      return img_optseek(image, image->rowstart[y]);
    case 3:
      return img_optseek(image, image->rowstart[y+z*image->ysize]);
    default:
      i_errhdlr("img_seek: weird dim\n",0,0,0,0);
      break;
    }
  } else
    i_errhdlr("img_seek: weird image type\n",0,0,0,0);
  return((unsigned long)-1);
}

int img_badrow(RGB_IMAGE *image, unsigned int y, unsigned int z) {
  if (y>=image->ysize || z>=image->zsize)
    return 1;
  else
    return 0;
}

int img_write(RGB_IMAGE *image, char *buffer,int count) {
  int retval;

  retval =  write(image->file,buffer,count);
  if (retval == count)
    image->offset += count;
  else
    image->offset = -1;
  return retval;
}

int img_read(RGB_IMAGE *image, char *buffer, int count) {
  int retval;

  retval =  read(image->file,buffer,count);
  if (retval == count)
    image->offset += count;
  else
    image->offset = -1;
  return retval;
}

unsigned long img_optseek(RGB_IMAGE *image, unsigned long offset) {
  if (image->offset != offset) {
    image->offset = offset;
    return ((unsigned long) lseek(image->file,offset,0));
  }
  return offset;
}

long img_getrowsize(RGB_IMAGE *image) {
  switch (image->dim) {
  case 1:
    return image->rowsize[0];
  case 2:
    return image->rowsize[image->y];
  case 3:
    return image->rowsize[image->y+image->z*image->ysize];
  }
  return(0L) ;
}

void img_setrowsize(RGB_IMAGE *image, long cnt, long y, long z) {
  int *sizeptr;

  if (img_badrow(image,y,z))
    return;
  switch (image->dim) {
  default:
  case 1:
    sizeptr = &image->rowsize[0];
    image->rowstart[0] = image->rleend;
    break;
  case 2:
    sizeptr = &image->rowsize[y];
    image->rowstart[y] = image->rleend;
    break;
  case 3:
    sizeptr = &image->rowsize[y+z*image->ysize];
    image->rowstart[y+z*image->ysize] = image->rleend;
  }
  if (*sizeptr != -1)
    image->wastebytes += *sizeptr;
  *sizeptr = cnt;
  image->rleend += cnt;
}

#define docompact               \
  while(iptr<ibufend) {           \
      sptr = iptr;            \
      iptr += 2;              \
      while((iptr<ibufend)&&((iptr[-2]!=iptr[-1])||(iptr[-1]!=iptr[0])))\
    iptr++;             \
      iptr -= 2;              \
      count = iptr-sptr;            \
      while(count) {            \
    todo = count>126 ? 126:count;         \
    count -= todo;            \
    *optr++ = 0x80|todo;          \
    while(todo--)           \
        *optr++ = *sptr++;          \
      }               \
      sptr = iptr;            \
      cc = *iptr++;           \
      while( (iptr<ibufend) && (*iptr == cc) )      \
    iptr++;             \
      count = iptr-sptr;            \
      while(count) {            \
    todo = count>126 ? 126:count;         \
    count -= todo;            \
    *optr++ = todo;           \
    *optr++ = cc;           \
      }               \
  }               \
  *optr++ = 0;

int img_rle_compact(unsigned short *expbuf, int ibpp,
                    unsigned short *rlebuf, int obpp, int cnt) {
  if (ibpp == 1 && obpp == 1) {
    unsigned char *iptr = (unsigned char *)expbuf;
    unsigned char *ibufend = iptr+cnt;
    unsigned char *sptr;
    unsigned char *optr = (unsigned char *)rlebuf;
    short todo, cc;
    long count;

    docompact;
    return optr - (unsigned char *)rlebuf;
  } else if (ibpp == 1 && obpp == 2) {
    unsigned char *iptr = (unsigned char *)expbuf;
    unsigned char *ibufend = iptr+cnt;
    unsigned char *sptr;
    unsigned short *optr = rlebuf;
    short todo, cc;
    long count;

    docompact;
    return optr - rlebuf;
  } else if (ibpp == 2 && obpp == 1) {
    unsigned short *iptr = expbuf;
    unsigned short *ibufend = iptr+cnt;
    unsigned short *sptr;
    unsigned char *optr = (unsigned char *)rlebuf;
    short todo, cc;
    long count;

    docompact;
    return optr - (unsigned char *)rlebuf;
  } else if (ibpp == 2 && obpp == 2) {
    unsigned short *iptr = expbuf;
    unsigned short *ibufend = iptr+cnt;
    unsigned short *sptr;
    unsigned short *optr = rlebuf;
    short todo, cc;
    long count;

    docompact;
    return optr - rlebuf;
  } else  {
    i_errhdlr("rle_compact: bad bpp: %d %d\n",ibpp,obpp,0,0);
    return 0;
  }
}

#define doexpand        \
  while(1) {        \
      pixel = *iptr++;      \
      if ( !(count = (pixel & 0x7f)) )  \
    return;       \
      if(pixel & 0x80) {      \
         while(count--)     \
        *optr++ = *iptr++;    \
      } else {        \
         pixel = *iptr++;     \
         while(count--)     \
        *optr++ = pixel;    \
      }         \
  }

void img_rle_expand(unsigned short *rlebuf, int ibpp,
                    unsigned short *expbuf, int obpp) {
  if (ibpp == 1 && obpp == 1) {
    unsigned char *iptr = (unsigned char *)rlebuf;
    unsigned char *optr = (unsigned char *)expbuf;
    unsigned short pixel,count;

    doexpand;
  } else if (ibpp == 1 && obpp == 2) {
    unsigned char *iptr = (unsigned char *)rlebuf;
    unsigned short *optr = expbuf;
    unsigned short pixel,count;

    doexpand;
  } else if (ibpp == 2 && obpp == 1) {
    unsigned short *iptr = rlebuf;
    unsigned char  *optr = (unsigned char *)expbuf;
    unsigned short pixel,count;

    doexpand;
  } else if (ibpp == 2 && obpp == 2) {
    unsigned short *iptr = rlebuf;
    unsigned short *optr = expbuf;
    unsigned short pixel,count;

    doexpand;
  } else
    i_errhdlr("rle_expand: bad bpp: %d %d\n",ibpp,obpp,0,0);
}

int putrow_uc(RGB_IMAGE *image, unsigned char *buffer, unsigned int y, unsigned int z) {
  unsigned char   *sptr;
  unsigned char      *cptr;
  unsigned int x;
  unsigned long min, max;
  long cnt;

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
      } else {
        return cnt;
      }
    case 2:
      printf("ERROR: this rgb save function BPP=2 is not implemented\n");
#if __GNUC__  >= 8
      [[gnu::fallthrough]];
#endif
    default:
      i_errhdlr("putrow: weird bpp\n",0,0,0,0);
    }
  } else if (ISRLE(image->type)) {
    printf("ERROR: compressed RGB is not implemented !\n");
  } else
    i_errhdlr("putrow: weird image type\n",0,0,0,0);
  return(-1);
}

int putrow(RGB_IMAGE *image, unsigned short *buffer, unsigned int y, unsigned int z) {
  unsigned short   *sptr;
  unsigned char      *cptr;
  unsigned int x;
  unsigned long min, max;
  long cnt;

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
  short i;
  unsigned char *cptr;
  unsigned short *sptr;
  short cnt;

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

RGB_IMAGE *imgopen(int, const char *, const char *,unsigned int, unsigned int, unsigned int, unsigned int, unsigned int);

RGB_IMAGE *iopen(const char *file, const char *mode, unsigned int type, unsigned int dim, unsigned int xsize, unsigned int ysize, unsigned int zsize) {
  return(imgopen(0, file, mode, type, dim, xsize, ysize, zsize));
}

RGB_IMAGE *fiopen(int f, const char *mode, unsigned int type, unsigned int dim, unsigned int xsize, unsigned int ysize, unsigned int zsize) {
  return(imgopen(f, 0, mode, type, dim, xsize, ysize, zsize));
}

RGB_IMAGE *imgopen(int f, const char *file, const char *mode, unsigned int type, unsigned int dim, unsigned int xsize, unsigned int ysize, unsigned int zsize) {
  RGB_IMAGE  *image;
  int rw;
  int tablesize;
  int i, max;

  image = (RGB_IMAGE*)calloc(1,sizeof(RGB_IMAGE));
  if (!image ) {
    i_errhdlr("iopen: error on image struct alloc\n",0,0,0,0);
    return NULL;
  }
  rw = mode[1] == '+';
  if (rw) {
    i_errhdlr("iopen: read/write mode not supported\n",0,0,0,0);
    return NULL;
  }
  if (*mode=='w') {
    if (file) {
      f = creat(file, 0666);
      if (rw && f>=0) {
        close(f);
        f = open(file, 2);
      }
    }
    if (f < 0) {
      i_errhdlr("iopen: can't open output file %s\n",0,0,0,0);
      return NULL;
    }
    image->imagic = IMAGIC;
    image->type = type;
    image->xsize = xsize;
    image->ysize = 1;
    image->zsize = 1;
    if (dim>1)
      image->ysize = ysize;
    if (dim>2)
      image->zsize = zsize;
    if (image->zsize == 1) {
      image->dim = 2;
      if (image->ysize == 1)
        image->dim = 1;
    } else {
      image->dim = 3;
    }
    image->min = 10000000;
    image->max = 0;
    isetname(image,"no name");
    image->wastebytes = 0;
    image->dorev = 0;
    swapImage(image) ;
    if (write(f,image,sizeof(RGB_IMAGE)) != sizeof(RGB_IMAGE)) {
      i_errhdlr("iopen: error on write of image header\n",0,0,0,0);
      return NULL;
    }
    swapImage(image) ;
  } else {
    if (file)
      f = open(file, rw? 2: 0);
    if (f < 0)
      return(NULL);
    if (read(f,image,sizeof(RGB_IMAGE)) != sizeof(RGB_IMAGE)) {
      i_errhdlr("iopen: error on read of image header\n",0,0,0,0);
      return NULL;
    }
    if ( ((image->imagic>>8) | ((image->imagic&0xff)<<8))
         == IMAGIC ) {
      image->dorev = 1;
      cvtimage((long *)image);
    } else
      image->dorev = 0;
    if (image->imagic != IMAGIC) {
      i_errhdlr("iopen: bad magic in image file %x\n",image->imagic,0,0,0);
      return NULL;
    }
  }
  if (rw)
    image->flags = _IORW;
  else if (*mode != 'r')
    image->flags = _IOWRT;
  else
    image->flags = _IOREAD;
  if (ISRLE(image->type)) {
    tablesize = image->ysize*image->zsize*sizeof(long);
    image->rowstart = (unsigned int *)malloc(tablesize);
    image->rowsize = (int *)malloc(tablesize);
    if ( image->rowstart == 0 || image->rowsize == 0 ) {
      i_errhdlr("iopen: error on table alloc\n",0,0,0,0);
      return NULL;
    }
    image->rleend = 512L+2*tablesize;
    if (*mode=='w') {
      max = image->ysize*image->zsize;
      for (i=0; i<max; i++) {
        image->rowstart[i] = 0;
        image->rowsize[i] = -1;
      }
    } else {
      tablesize = image->ysize*image->zsize*sizeof(long);
      lseek(f, 512L, 0);
      if (read(f,image->rowstart,tablesize) != tablesize) {
        i_errhdlr("iopen: error on read of rowstart\n",0,0,0,0);
        return NULL;
      }
      if (image->dorev)
        cvtlongs((long *)image->rowstart,tablesize);
      if (read(f,image->rowsize,tablesize) != tablesize) {
        i_errhdlr("iopen: error on read of rowsize\n",0,0,0,0);
        return NULL;
      }
      if (image->dorev)
        cvtlongs((long *)image->rowsize,tablesize);
    }
  }
  image->cnt = 0;
  image->ptr = 0;
  image->base = 0;
  if ( (image->tmpbuf = ibufalloc(image)) == 0 ) {
    i_errhdlr("iopen: error on tmpbuf alloc %d\n",image->xsize,0,0,0);
    return NULL;
  }
  image->x = image->y = image->z = 0;
  image->file = f;
  image->offset = 512L;     /* set up for img_optseek */
  lseek(image->file, 512L, 0);
  return(image);
}

unsigned short *ibufalloc(RGB_IMAGE *image) {
  return (unsigned short *)malloc(IBUFSIZE(image->xsize));
}

long reverse(unsigned long lwrd) {
  return ((lwrd>>24)    |
          (lwrd>>8 & 0xff00)   |
          (lwrd<<8 & 0xff0000) |
          (lwrd<<24)     );
}

void cvtshorts( unsigned short *buffer, long n) {
  short i;
  long nshorts = n>>1;
  unsigned short swrd;

  for (i=0; i<nshorts; i++) {
    swrd = *buffer;
    *buffer++ = (swrd>>8) | (swrd<<8);
  }
}


void cvtlongs( long *buffer, long n) {
  short i;
  long nlongs = n>>2;
  unsigned long lwrd;

  for (i=0; i<nlongs; i++) {
    lwrd = buffer[i];
    buffer[i] =     ((lwrd>>24)     |
                     (lwrd>>8 & 0xff00)  |
                     (lwrd<<8 & 0xff0000)  |
                     (lwrd<<24)    );
  }
}

void cvtimage( long *buffer) {
  cvtshorts((unsigned short *)buffer,12);
  cvtlongs(buffer+3,12);
  cvtlongs(buffer+26,4);
}

static void (*i_errfunc)(char *ebuf);

/*  error handler for the image library.  If the iseterror() routine
    has been called, sprintf's the args into a string and calls the
    error function.  Otherwise calls fprintf with the args and then
    exit.  This allows 'old' programs to assume that no errors
    ever need be worried about, while programs that know how and
    want to can handle the errors themselves.  Olson, 11/88
*/
/* most args currently used is 2 */
void i_errhdlr(const char *fmt, int a1, int a2, int a3, int a4) {
  if (i_errfunc) {
    char ebuf[2048];  /* be generous; if an error includes a
                                         pathname, the maxlen is 1024, so we shouldn't ever
                                         overflow this! */
    sprintf(ebuf, fmt, a1, a2, a3, a4);
    (*i_errfunc)(ebuf);
    return;
  }
  fprintf(stderr, fmt, a1, a2, a3, a4);
  exit(1);
}

/* this function sets the error handler for i_errhdlr */
void i_seterror(void (*func)(char *)) {
  i_errfunc = func;
}

typedef union {
  short  s ;
  char   buf[sizeof(short)] ;
} SSWAP_SHORT ;

typedef union {
  float f ;
  int   i ;
  char  buf[4] ;
  short s[2] ;
} SSWAP_LONG ;

short tmpswapShort(short s) {
  SSWAP_SHORT ss ;
  char       c ;

  /* first swap bytes in word */
  ss.s = s ;
  c = ss.buf[0] ;
  ss.buf[0] = ss.buf[1] ;
  ss.buf[1] = c ;

  return(ss.s) ;
}

int tmpswapInt(int i) {
  SSWAP_LONG  sl ;
  short      s ;

  /* first swap bytes in each word */
  sl.i = i ;
  sl.s[0] = tmpswapShort(sl.s[0]) ;
  sl.s[1] = tmpswapShort(sl.s[1]) ;

  /* now swap words */
  s = sl.s[0] ;
  sl.s[0] = sl.s[1] ;
  sl.s[1] = s ;

  return(sl.i) ;
}

void swapImage(RGB_IMAGE *image) {
#if (BYTE_ORDER == LITTLE_ENDIAN)
  image->imagic = tmpswapShort(image->imagic) ;
  image->type = tmpswapShort(image->type) ;
  image->dim = tmpswapShort(image->dim) ;
  image->xsize = tmpswapShort(image->xsize) ;
  image->ysize = tmpswapShort(image->ysize) ;
  image->zsize = tmpswapShort(image->zsize) ;
  image->min = tmpswapInt(image->min) ;
  image->max = tmpswapInt(image->max) ;
  image->wastebytes = tmpswapShort(image->wastebytes) ;
  image->colormap = tmpswapInt(image->colormap) ;
#endif
}
