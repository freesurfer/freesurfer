#ifndef RGB_H
#define RGB_H

#include <stdio.h>

/*
 *  Defines for image files . . . .
 *
 *        Paul Haeberli - 1984
 *      Look in /usr/people/4Dgifts/iristools/imgtools for example code!
 *
 */

#define IMAGIC  0732
#define CM_NORMAL 0
#define CM_DITHERED 1
#define CM_SCREEN 2
#define CM_COLORMAP 3

#define TYPEMASK 0xff00
#define BPPMASK 0x00ff
#define ITYPE_VERBATIM 0x0000
#define ITYPE_RLE 0x0100
#define ISRLE(type) (((type) & 0xff00) == ITYPE_RLE)
#define ISVERBATIM(type) (((type) & 0xff00) == ITYPE_VERBATIM)
#define BPP(type) ((type) & BPPMASK)
#define RLE(bpp) (ITYPE_RLE | (bpp))
#define UNCOMPRESSED(bpp) (ITYPE_VERBATIM | (bpp))
#define IBUFSIZE(pixels) ((pixels+(pixels>>6))<<2)
#define RLE_NOP 0x00

#define ierror(p) (((p)->flags&_IOERR)!=0)
#define ifileno(p) ((p)->file)

typedef struct {
  unsigned short  imagic;
  unsigned short  type;
  unsigned short  dim;
  unsigned short  xsize;
  unsigned short  ysize;
  unsigned short  zsize;
  unsigned int  min;
  unsigned int  max;
  unsigned int  wastebytes;
  char    name[80];
  unsigned int  colormap;
  int     file;
  unsigned short  flags;
  short   dorev;
  short   x;
  short   y;
  short   z;
  short   cnt;
  unsigned short  *ptr;
  unsigned short  *base;
  unsigned short  *tmpbuf;
  unsigned int  offset;
  unsigned int  rleend;
  unsigned int  *rowstart;
  int     *rowsize;
}
RGB_IMAGE;

RGB_IMAGE *icreate();
RGB_IMAGE *iopen(const char *file, const char *mode, unsigned int type, unsigned int dim, unsigned int xsize, unsigned int ysize, unsigned int zsize);
RGB_IMAGE *fiopen(int f, const char *mode, unsigned int type, unsigned int dim, unsigned int xsize, unsigned int ysize, unsigned int zsize);
long reverse(unsigned long lwrd);
void cvtshorts( unsigned short *buffer, long n);
void i_seterror(void (*func)(char *));
int getpix(RGB_IMAGE *image);
unsigned int putpix(RGB_IMAGE *image, unsigned int pix);
int img_read(RGB_IMAGE *image, char *buffer, int count);
int img_write(RGB_IMAGE *image, char *buffer,int count);
int img_badrow(RGB_IMAGE *image, unsigned int y, unsigned int z);
unsigned long img_seek(RGB_IMAGE *image, unsigned int y, unsigned int z);
long img_getrowsize(RGB_IMAGE *image);
void img_setrowsize(RGB_IMAGE *image, long cnt, long y, long z);
int img_rle_compact(unsigned short *expbuf, int ibpp, unsigned short *rlebuf, int obpp, int cnt);
void img_rle_expand(unsigned short *rlebuf, int ibpp, unsigned short *expbuf, int obpp);
int putrow(RGB_IMAGE *image, unsigned short *buffer, unsigned int y, unsigned int z);
int getrow(RGB_IMAGE *image, unsigned short *buffer, unsigned int y, unsigned int z);

unsigned short *ibufalloc(RGB_IMAGE *image);
int ifilbuf(RGB_IMAGE *image);
int iflush(RGB_IMAGE *image);
unsigned int iflsbuf(RGB_IMAGE *image, unsigned int c);
void isetname(RGB_IMAGE *image, const char *name);
void isetcolormap(RGB_IMAGE *image, int colormap);

int iclose(RGB_IMAGE *image);
int putrow(RGB_IMAGE *image, unsigned short *buffer,unsigned int y,unsigned int z);
int putrow_uc(RGB_IMAGE *image, unsigned char *buffer, unsigned int y, unsigned int z);
int getrow(RGB_IMAGE *image, unsigned short *buffer,unsigned int y,unsigned int z);

#define IMAGEDEF

unsigned long img_optseek(RGB_IMAGE *image, unsigned long offset);
int img_write(RGB_IMAGE *image, char *buffer,int count);
void cvtimage( long *buffer);
void cvtlongs( long *buffer, long n);
void i_errhdlr(const char *fmt, int a1, int a2, int a3, int a4);
void swapImage(RGB_IMAGE *image);

#ifndef _IOREAD
#define _IOREAD         0001
#endif
#ifndef _IOWRT
#define _IOWRT          0002
#endif
#ifndef _IOEOF
#define _IOEOF          0020
#endif
#ifndef _IOERR
#define _IOERR          0040
#endif
#ifndef _IORW
#define _IORW           0200
#endif

#endif
