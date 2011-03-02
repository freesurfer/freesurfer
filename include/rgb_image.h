/**
 * @file  rgb_image.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.12 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef __GL_RGB_IMAGE_H__
#define __GL_RGB_IMAGE_H__
#ifdef __cplusplus
extern "C"
{
#endif


  /*
   *  Defines for image files . . . .
   *
   *        Paul Haeberli - 1984
   *      Look in /usr/people/4Dgifts/iristools/imgtools for example code!
   *
   */

#include <stdio.h>

#define IMAGIC  0732

  /* colormap of images */
#define CM_NORMAL   0 /* file contains rows of values which
  * are either RGB values (zsize == 3)
  * or greyramp values (zsize == 1) */
#define CM_DITHERED   1
#define CM_SCREEN   2 /* file contains data which is a screen
  * image;
  getrow returns buffer which
               * can be displayed directly with
               * writepixels */
#define CM_COLORMAP   3 /* a colormap file */

#define TYPEMASK    0xff00
#define BPPMASK     0x00ff
#define ITYPE_VERBATIM    0x0000
#define ITYPE_RLE   0x0100
#define ISRLE(type)   (((type) & 0xff00) == ITYPE_RLE)
#define ISVERBATIM(type)  (((type) & 0xff00) == ITYPE_VERBATIM)
#define BPP(type)   ((type) & BPPMASK)
#define RLE(bpp)    (ITYPE_RLE | (bpp))
#define UNCOMPRESSED(bpp)   (ITYPE_VERBATIM | (bpp))
#define IBUFSIZE(pixels)  ((pixels+(pixels>>6))<<2)
#define RLE_NOP     0x00

#define ierror(p)   (((p)->flags&_IOERR)!=0)
#define ifileno(p)    ((p)->file)
#define getpix(p)   (--(p)->cnt>=0 ? *(p)->ptr++ : ifilbuf(p))
#define putpix(p,x)   (--(p)->cnt>=0 \
            ? ((int)(*(p)->ptr++=(unsigned)(x))) \
            : iflsbuf(p,(unsigned)(x)))

               typedef struct
               {
                 unsigned short  imagic;   /* stuff saved on disk . . */
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

                 int     file;   /* stuff used in core only */
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
                 unsigned int  rleend;   /* for rle images */
                 unsigned int  *rowstart;  /* for rle images */
                 int     *rowsize; /* for rle images */
               }
               RGB_IMAGE;

  RGB_IMAGE *icreate();
  RGB_IMAGE *iopen(char *file, char *mode, unsigned int type, unsigned int dim,
                   unsigned int xsize, unsigned int ysize, unsigned int zsize);
  RGB_IMAGE *fiopen(int f, char *mode, unsigned int type, unsigned int dim,
                    unsigned int xsize, unsigned int ysize, unsigned int zsize);
  long reverse(unsigned long lwrd)  ;
  void cvtshorts( unsigned short *buffer, long n) ;
  void i_seterror(void (*func)(char *)) ;
#undef getpix
#undef putpix
  int getpix(RGB_IMAGE *image) ;
  unsigned int putpix(RGB_IMAGE *image, unsigned int pix) ;
  int img_read(RGB_IMAGE *image, char *buffer, int count) ;
  int img_write(RGB_IMAGE *image, char *buffer,int count) ;
  int img_badrow(RGB_IMAGE *image, unsigned int y, unsigned int z) ;
  unsigned long img_seek(RGB_IMAGE *image, unsigned int y, unsigned int z) ;
  long img_getrowsize(RGB_IMAGE *image) ;
  void img_setrowsize(RGB_IMAGE *image, long cnt, long y, long z) ;
  int img_rle_compact(unsigned short *expbuf, int ibpp,
                      unsigned short *rlebuf, int obpp, int cnt) ;
  void img_rle_expand(unsigned short *rlebuf, int ibpp,
                      unsigned short *expbuf, int obpp) ;
  int putrow(RGB_IMAGE *image, unsigned short *buffer,
             unsigned int y, unsigned int z)  ;
  int getrow(RGB_IMAGE *image, unsigned short *buffer,
             unsigned int y, unsigned int z)  ;

  /*
   *
   * ...while iopen and fiopen can take an extended set of parameters, the
   * last five are optional, so a more correct prototype would be:
   *
   *
   * RGB_IMAGE *iopen(char *file, char *mode, ...);
   * RGB_IMAGE *fiopen(int f, char *mode, ...);
   */

  unsigned short *ibufalloc(RGB_IMAGE *image);
  int ifilbuf(RGB_IMAGE *image);
  int iflush(RGB_IMAGE *image);
  unsigned int iflsbuf(RGB_IMAGE *image, unsigned int c);
  void isetname(RGB_IMAGE *image, char *name);
  void isetcolormap(RGB_IMAGE *image, int colormap);

  int iclose(RGB_IMAGE *image);
  int putrow(RGB_IMAGE *image, unsigned short *buffer,unsigned int y,unsigned int z);
  int putrow_uc(RGB_IMAGE *image, unsigned char *buffer,
                unsigned int y, unsigned int z)  ;
  int getrow(RGB_IMAGE *image, unsigned short *buffer,unsigned int y,unsigned int z);

  /* removed at tosa's request
  RGB_IMAGE *iopen();
  RGB_IMAGE *icreate();

  */

  /* RKT: why is there a second definition here? why hasn't anything
     complained before? */
  /* unsigned short *ibufalloc(); */

#define IMAGEDEF    /* for backwards compatibility */
#ifdef __cplusplus
}
#endif

unsigned long img_optseek(RGB_IMAGE *image, unsigned long offset);
int img_write(RGB_IMAGE *image, char *buffer,int count) ;
void cvtimage( long *buffer) ;
void cvtlongs( long *buffer, register long n) ;
void i_errhdlr(char *fmt, int a1, int a2, int a3, int a4)  ;
void swapImage(RGB_IMAGE *image) ;

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

#endif  /* !__GL_IMAGE_H__ */
