//
// rgbutils.c
//
/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <fcntl.h>
#include <unistd.h>  /* for SEEK_ constants */

#include "hmem.h"
#include <hipl_format.h>

#include "hips.h"

#include "image.h"
#include "error.h"
#include "matrix.h"
#include "matfile.h"
#include "utils.h"
#include "macros.h"
#include "machine.h"
#include "proto.h"
#include "diag.h"
#include "canny.h"
#include "rgb_image.h"
#ifndef IRIX
#include "pgm.h"
#include "ppm.h"
#include "pbm.h"
#endif

int    RGBwrite(IMAGE *I, char *fname, int frame) ;
IMAGE *RGBReadImage(char *fname);
IMAGE *RGBReadHeader(char *fname, IMAGE *);


IMAGE *
RGBReadHeader(char *fname, IMAGE *I) {
  RGB_IMAGE *rgb;

  rgb = iopen(fname, "r", 0, 0, 0, 0, 0);

  if (!I)
    I = ImageAlloc(rgb->ysize, rgb->xsize, PFRGB, 1) ;
  else
    init_header(I, "orig", "seq", 1, "today", rgb->ysize,
                rgb->xsize, PFRGB, 1, "temp");

  iclose(rgb);

  return(I) ;
}

IMAGE *RGBReadImage(char *fname) {
  IMAGE *I;
  RGB_IMAGE *rgb;
  unsigned short rows,cols,*r,*g,*b,i,j,*tr,*tg,*tb;
  byte *iptr;

  rgb = iopen(fname, "r", 0, 0, 0, 0, 0);
  rows = rgb->ysize;
  cols = rgb->xsize;

  if (rgb->zsize>3)
    ErrorReturn(NULL, (ERROR_BAD_PARM,
                       "Too many color planes in RGBReadImage (%s)\n",fname));

  I = ImageAlloc(rows, cols, PFRGB, 1);

  if ((r = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  if ((g = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  if ((b = (unsigned short *)malloc(sizeof(unsigned short)*cols)) == NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate color buffer\n");

  iptr = I->image;

  for (i=0;i<rows;i++) {
    getrow(rgb,r,i,0); /* Red */
    getrow(rgb,g,i,1); /* Green */
    getrow(rgb,b,i,2); /* Blue */

    /* Translate color planes to RGB format */
    tr = r;
    tg = g;
    tb = b;
    for (j=0;j<cols;j++) {
      *iptr++ = *tr++;
      *iptr++ = *tg++;
      *iptr++ = *tb++;
    }
  }

  free(r);
  free(g);
  free(b);
  iclose(rgb);

  return I;
}

int
RGBwrite(IMAGE *I, char *fname, int frame) {
  RGB_IMAGE  *image ;
  int    x, y ;
  unsigned short *r ;

#ifndef Linux
  image = iopen(fname,"w",RLE(1), 2, I->cols, I->rows, 1);
#else
  image = iopen(fname,"w",UNCOMPRESSED(1), 2, I->cols, I->rows, 1);
#endif
  r = (unsigned short *)calloc(I->cols, sizeof(unsigned short)) ;
  for (y = 0 ; y < I->rows; y++) {
    for (x = 0 ; x < I->cols ; x++)
      r[x] = (unsigned short)(*IMAGEpix(I, x, y)) ;

    /* fill rbuf, gbuf, and bbuf with pixel values */
    putrow(image, r, y, 0);    /* red row */
    putrow(image, r, y, 1);    /* green row */
    putrow(image, r, y, 2);    /* blue row */
  }
  iclose(image);
  free(r) ;
  return(NO_ERROR) ;
}

