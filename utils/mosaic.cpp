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

/**************************************************************
  Name: mosaic.c
  Author: Douglas Greve
  Date: 09/12/01
  Description: routines to convert between mosaics and volumes

  ncvol - number of columns in the volume
  nrvol - number of rows in the volume
  nsvol - number of slices in the volume

  cvol  - a particular column in the volume
  rvol  - a particular row in the volume
  svol  - a particular slice in the volume

  ncmos - number of pixel columns in the mosaic
  nrmos - number of pixel rows in the mosaic

  nctmos - number of tile columns in the mosaic
  nrtmos - number of tile rows in the mosaic

  ctmos - a particular column tile in the mosaic
  rtmos - a particular row tile in the mosaic

  cmos - a particular column pixel in the mosaic
  rmos - a particular row pixelin the mosaic

  Notes:
    1. Indices always start at 0.

**************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef Darwin
#include <malloc.h>
#else
#include "proto.h"
#endif

#include "math.h"

#include "mosaic.h"

//#define _DEBUG

/*--------------------------------------------------------
  VolSS2MosSS() - converts from a volume subscript (CRS)
  to a mosaic subscript (CR). If the point is outside of
  the range of the mosaic, OutOfBounds will be set to 1.
  Returns:
    0 - no error
    1 - if the number of columns/rows in the mosaic is not
        an integer multiple of the number of columns/rows
        in the image
--------------------------------------------------------*/
int VolSS2MosSS(
    int cvol, int rvol, int svol, int ncvol, int nrvol, int ncmos, int nrmos, int *cmos, int *rmos, int *OutOfBounds)
{
  int nctmos;
  // int nrtmos;
  int ctmos, rtmos;

  if ((ncmos % ncvol) != 0) {
    fprintf(stderr,
            "ERROR: VolSS2MosSS: number of columns in the mosaic (%d) \n"
            "is not an integer multiple of the number of columns in \n"
            "the image (%d)\n",
            ncmos,
            ncvol);
    return (1);
  }

  if ((nrmos % nrvol) != 0) {
    fprintf(stderr,
            "ERROR: VolSS2MosSS: number of row in the mosaic (%d) \n"
            "is not an integer multiple of the number of rows in \n"
            "the image (%d)\n",
            nrmos,
            nrvol);
    return (1);
  }

  /* number of col and row tiles in the mosaic */
  nctmos = ncmos / ncvol;
  // nrtmos = nrmos / nrvol;

  /* the row and col tile of the given slice */
  rtmos = (int)(floor(svol / nctmos));
  ctmos = svol - rtmos * nctmos;

  *cmos = ctmos * ncvol + cvol;
  *rmos = rtmos * nrvol + rvol;

  if (OutOfBounds != NULL) {
    if (*cmos < 0 || *cmos > (ncmos - 1) || *rmos < 0 || *rmos > (nrmos - 1))
      *OutOfBounds = 1;
    else
      *OutOfBounds = 0;
  }

  return (0);
}
/*--------------------------------------------------------
  MosSS2VolSS() - converts from a mosaic subscript (CR)
  to a volume subscript (CRS). If the point is outside of
  the range of the volume, OutOfBounds will be set to 1.
  Returns:
    0 - no error
    1 - if the number of columns/rows in the mosaic is not
        an integer multiple of the number of columns/rows
        in the image
--------------------------------------------------------*/
int MosSS2VolSS(int cmos,
                int rmos,
                int ncmos,
                int nrmos,
                int ncvol,
                int nrvol,
                int nsvol,
                int *cvol,
                int *rvol,
                int *svol,
                int *OutOfBounds)
{
  int nctmos;
  // int nrtmos;
  int ctmos, rtmos;

  if ((ncmos % ncvol) != 0) {
    fprintf(stderr,
            "ERROR: MosSS2Vol: number of columns in the mosaic (%d) \n"
            "is not an integer multiple of the number of columns in \n"
            "the image (%d)\n",
            ncmos,
            ncvol);
    return (1);
  }

  if ((nrmos % nrvol) != 0) {
    fprintf(stderr,
            "ERROR: MosSS2Vol: number of row in the mosaic (%d) \n"
            "is not an integer multiple of the number of rows in \n"
            "the image (%d)\n",
            nrmos,
            nrvol);
    return (1);
  }

  /* number of col and row tiles in the mosaic */
  nctmos = ncmos / ncvol;
  // nrtmos = nrmos / nrvol;

  /* the row and col tile of the given point in the mosaic */
  ctmos = (int)(floor(cmos / ncvol));
  rtmos = (int)(floor(rmos / nrvol));

  *cvol = cmos - ctmos * ncvol;
  *rvol = rmos - rtmos * nrvol;
  *svol = ctmos + rtmos * nctmos;

  if (OutOfBounds != NULL) {
    if (*cvol < 0 || *cvol > (ncvol - 1) || *rvol < 0 || *rvol > (nrvol - 1) || *svol < 0 || *svol > (nsvol - 1))
      *OutOfBounds = 1;
    else
      *OutOfBounds = 0;
  }

#if 0
  printf("----------------------------------\n");
  printf("cmos = %d\n",cmos);
  printf("rmos = %d\n",rmos);
  printf("ncmos = %d\n",ncmos);
  printf("nrmos = %d\n",nrmos);
  printf("nctmos = %d\n",nctmos);
  printf("nrtmos = %d\n",nrtmos);
  printf("ctmos = %d\n",ctmos);
  printf("rtmos = %d\n",rtmos);
  printf("cvol = %d\n",*cvol);
  printf("rvol = %d\n",*rvol);
  printf("svol = %d\n",*svol);
#endif

  return (0);
}
/*----------------------------------------------------------
  CheckMosaic(void) - this is just a diagnostic that prints
  info to stdout.
  ----------------------------------------------------------*/
int CheckMosaic(void)
{
  int cvol, rvol, svol, ncvol, nrvol, nsvol;
  int ncmos, nrmos, cmos, rmos;
  int cvol2, rvol2, svol2;
  int oob;
  int ok = 0;
  int n;

  ncvol = 8;
  nrvol = 3;
  nsvol = 11;
  ncmos = 4 * ncvol;
  nrmos = 3 * nrvol;

  n = 0;
  for (svol = 0; svol < nsvol; svol++) {
    for (rvol = 0; rvol < nrvol; rvol++) {
      for (cvol = 0; cvol < ncvol; cvol++) {
        VolSS2MosSS(cvol, rvol, svol, ncvol, nrvol, ncmos, nrmos, &cmos, &rmos, &oob);
        MosSS2VolSS(cmos, rmos, ncmos, nrmos, ncvol, nrvol, nsvol, &cvol2, &rvol2, &svol2, &oob);
        ok = 0;
        if (cvol == cvol2 && rvol == rvol2 && svol == svol2) ok = 1;
        printf("%5d  %2d %2d %2d   %3d %3d    %2d %2d %2d  %d  %d\n",
               n,
               cvol,
               rvol,
               svol,
               cmos,
               rmos,
               cvol2,
               rvol2,
               svol2,
               oob,
               ok);
        n++;
      }
    }
  }
  printf("-----------------------------------\n");
  VolSS2MosSS(0, 0, 11, ncvol, nrvol, ncmos, nrmos, &cmos, &rmos, &oob);
  printf("%5d  %2d %2d %2d   %3d %3d    %2d %2d %2d  %d  %d\n",
         n,
         cvol,
         rvol,
         svol,
         cmos,
         rmos,
         cvol2,
         rvol2,
         svol2,
         oob,
         ok);

  MosSS2VolSS(cmos, rmos, ncmos, nrmos, ncvol, nrvol, nsvol, &cvol2, &rvol2, &svol2, &oob);
  printf("%5d  %2d %2d %2d   %3d %3d    %2d %2d %2d  %d  %d\n",
         n,
         cvol,
         rvol,
         svol,
         cmos,
         rmos,
         cvol2,
         rvol2,
         svol2,
         oob,
         ok);

  return (0);
}
