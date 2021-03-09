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

/*
  Name:   corio.c
  Author: Douglas N. Greve
  Date:   12/15/00

  Purpose: handles basic input/output for the MGH COR format. In this
  format, volumes have a dimension of 256x256x256. They are saved to
  disk as 256 files, each file is a corronal slice.  Each slice is
  saved with the sagital dimension as the faster index and the axial
  dimension as the slower index. Each coronal slice forms an image
  with rows corresponding to the axial index, and columns
  corresponding to the sagital index. As the column index goes from 0
  to 255, the anatomical location changes from right to left. As the
  row index goes from 0 to 255, the anatomical location changes from
  top to bottom. As the slice goes from 0 to 255, the anatomical
  location changes from back to front. Each value is saved as an
  unsigned char. Only one value per voxel is allowed (in, nframes=1).

  column: sagital:  R:  x:  0 = right, 255 = left
  slice:  corronal: A:  y:  0 = back,  255 = front
  row:    axial:    S:  z:  0 = top,   255 = bottom

  The volumes are stored in unsigned char**, eg
     unsigned char **COR;

  Then, COR[slice] is a pointer to a corronal slice, and
  COR[slice][col+256*row] is the value of the volume at the given
  row, col, and slice.

*/
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
extern int errno;

#include "corio.h"

/*-------------------------------------*/
unsigned char **alloc_cor(void)
{
  unsigned char **COR;
  int n;

  COR = (unsigned char **)calloc(256, sizeof(unsigned char *));
  if (COR == NULL) {
    fprintf(stderr, "alloc_cor(): could not alloc unsigned char **\n");
    return (NULL);
  }

  for (n = 0; n < 256; n++) {
    COR[n] = (unsigned char *)calloc(256 * 256, sizeof(unsigned char));
    if (COR[n] == NULL) {
      fprintf(stderr, "alloc_cor(): could not alloc slice %d\n", n);
      return (NULL);
    }
  }
  return (COR);
}
/*-------------------------------------*/
int free_cor(unsigned char ***pppCOR)
{
  int n;
  unsigned char **COR;

  COR = *pppCOR;

  for (n = 0; n < 256; n++) free(COR[n]);
  free(COR);

  return (0);
}

/*-------------------------------------*/
unsigned char **ld_cor(char *cordir)
{
  unsigned char **COR;
  int n;
  char fname[1000];
  FILE *fp;
  int nread, n_to_be_read = 256 * 256;

  COR = alloc_cor();
  if (COR == NULL) return (NULL);

  for (n = 0; n < 256; n++) {
    /*fprintf(stderr,"%3d ",n);
      if((n+1)%16==0) fprintf(stderr,"\n");*/

    sprintf(fname, "%s/COR-%03d", cordir, n + 1);
    fp = fopen(fname, "r");
    if (fp == NULL) {
      perror("ld_cor()");
      fprintf(stderr, "Could not open %s\n", fname);
      free_cor(&COR);
      return (NULL);
    }

    nread = fread(COR[n], sizeof(unsigned char), n_to_be_read, fp);
    fclose(fp);
    if (nread != n_to_be_read) {
      perror("ld_cor()");
      fprintf(stderr, "Error reading %s\n", fname);
      fprintf(stderr, "nread = %d, n-to-be-read = %d\n", nread, n_to_be_read);
      free_cor(&COR);
      return (NULL);
    }
  }
  return (COR);
}
/*-------------------------------------------------*/
int cordir_iswritable(char *cordir)
{
  char tmpstr[2000];
  FILE *fp;

  sprintf(tmpstr, "%s/junk-tmp.huh", cordir);
  fp = fopen(tmpstr, "w");
  if (fp == NULL) {
    perror("");
    fprintf(stderr, "Cannot write to %s\n", cordir);
    return (0);
  }
  fclose(fp);
  unlink(tmpstr);

  return (1);
}

/*-------------------------------------------------*/
int sv_cor(unsigned char **COR, char *cordir)
{
  int n;
  char fname[1000];
  FILE *fp;
  int nwritten, n_to_be_written = 256 * 256;

  if (!cordir_iswritable(cordir)) return (1);

  for (n = 0; n < 256; n++) {
    /*fprintf(stderr,"%3d ",n);
      if((n+1)%16==0) fprintf(stderr,"\n");*/

    sprintf(fname, "%s/COR-%03d", cordir, n + 1);
    fp = fopen(fname, "w");
    if (fp == NULL) {
      perror("sv_cor()");
      fprintf(stderr, "Could not open %s\n", fname);
      return (1);
    }

    nwritten = fwrite(COR[n], sizeof(unsigned char), n_to_be_written, fp);
    fclose(fp);
    if (nwritten != n_to_be_written) {
      perror("sv_cor()");
      fprintf(stderr, "Error writing %s\n", fname);
      fprintf(stderr, "nwriten = %d, n-to-be-written = %d\n", nwritten, n_to_be_written);
      return (1);
    }
  }

  /* write the COR-.info file */
  sprintf(fname, "%s/COR-.info", cordir);
  fp = fopen(fname, "w");
  if (fp == NULL) {
    perror("sv_cor()");
    fprintf(stderr, "Could not open %s\n", fname);
    return (1);
  }
  fprintf(fp, "imnr0 1\n");
  fprintf(fp, "imnr1 256\n");
  fprintf(fp, "ptype 2\n");
  fprintf(fp, "x 256\n");
  fprintf(fp, "y 256\n");
  fprintf(fp, "fov 0.256000\n");
  fprintf(fp, "thick 0.001000\n");
  fprintf(fp, "psiz 0.001000\n");
  fprintf(fp, "locatn 0.000000\n");
  fprintf(fp, "strtx -0.128000\n");
  fprintf(fp, "endx 0.128000\n");
  fprintf(fp, "strty -0.128000\n");
  fprintf(fp, "endy 0.128000\n");
  fprintf(fp, "strtz -0.128000\n");
  fprintf(fp, "endz 0.128000\n");
  fprintf(fp, "tr 0.000000\n");
  fprintf(fp, "te 0.000000\n");
  fprintf(fp, "ti 0.000000\n");

  return (0);
}
/*-------------------------------------------------*/
unsigned char getcorval(unsigned char **COR, int row, int col, int slc)
{
  if (col < 0 || col > 255) {
    fprintf(stderr, "getcorval: col = %d, out of bounds\n", col);
    return (0);
  }
  if (row < 0 || row > 255) {
    fprintf(stderr, "getcorval: row = %d, out of bounds\n", row);
    return (0);
  }
  if (slc < 0 || slc > 255) {
    fprintf(stderr, "getcorval: slc = %d, out of bounds\n", slc);
    return (0);
  }

  return (*(COR[slc] + col + row * 256));
}
/*-------------------------------------------------*/
int setcorval(unsigned char val, unsigned char **COR, int row, int col, int slc)
{
  if (col < 0 || col > 255) {
    fprintf(stderr, "setcorval: col = %d, out of bounds\n", col);
    return (1);
  }
  if (row < 0 || row > 255) {
    fprintf(stderr, "setcorval: row = %d, out of bounds\n", row);
    return (1);
  }
  if (slc < 0 || slc > 255) {
    fprintf(stderr, "setcorval: slc = %d, out of bounds\n", slc);
    return (1);
  }

  *(COR[slc] + col + row * 256) = val;

  return (0);
}
