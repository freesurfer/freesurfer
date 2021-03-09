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

/* MRI I/O - routines for reading and writing large files fast */
/* 2/1/91 - AD */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "proto.h"

#include "MRIio_old.h"

static void MGHprint_error(const char *str);

static void MGHprint_error(const char *str)
{
  printf("%s", str);
  exit(0);
}

char *lmalloc(unsigned long size)
{
  char *p;

  p = (char *)malloc(size);
  if (p == NULL) MGHprint_error("Cannot malloc()\n");
  return p;
}

char *lcalloc(size_t nmemb, size_t size)
{
  char *p;

  p = (char *)calloc(nmemb, size);
  if (p == NULL) MGHprint_error("Cannot calloc()\n");
  return p;
}

void file_name(const char *fpref, char *fname, int num, const char *form)
{
  char ext[10];

  sprintf(ext, form, num);
  strcpy(fname, fpref);
  strcat(fname, ext);
}

void buffer_to_image(unsigned char *buf, unsigned char **im, int ysize, int xsize)
{
  int i, j;
  unsigned long k;
  float sum;

  k = 0;
  sum = 0;
  for (i = 0; i < ysize; i++)
    for (j = 0; j < xsize; j++) {
      im[i][j] = buf[k++];
      sum += im[i][j];
    }
  /*
    printf("avg=%f\n",sum/(ysize*xsize));
  */
}

void image_to_buffer(unsigned char **im, unsigned char *buf, int ysize, int xsize)
{
  int i, j;
  unsigned long k;

  k = 0;
  for (i = 0; i < ysize; i++)
    for (j = 0; j < xsize; j++) buf[k++] = im[i][j];
}
