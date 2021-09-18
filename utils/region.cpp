/**
 * @brief utilities for the REGION data structure
 *
 */
/*
 * Original Author: Bruce Fischl
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

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "region.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int regionCornerCoords(MRI_REGION *r, int which_corner, int *px, int *py, int *pz);

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
MRI_REGION *REGIONalloc(void)
{
  MRI_REGION *r;

  r = (MRI_REGION *)calloc(1, sizeof(MRI_REGION));
  if (!r) ErrorExit(ERROR_NO_MEMORY, "REGIONalloc: could not allocate region");

  return (r);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          return the region in reg1 but not in reg2

          note that a cubic region is not a rich enough descriptor
          for this procedure, but mostly we are concerned with
          finding the portion of an advancing window which is
          new.
------------------------------------------------------*/
MRI_REGION *REGIONsubtract(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  int x1_start, x1_end, x2_end, y1_start, y1_end, y2_end, z1_start, z1_end, z2_end;

  x1_start = reg1->x;
  x1_end = reg1->x + reg1->dx - 1;
  x2_end = reg2->x + reg2->dx - 1;

  rdst->x = MAX(x1_start, x2_end);
  rdst->dx = x1_end - rdst->x + 1;
  if (rdst->dx < 0) rdst->dx = 0;

  y1_start = reg1->y;
  y1_end = reg1->y + reg1->dy - 1;
  y2_end = reg2->y + reg2->dy - 1;

  rdst->y = MAX(y1_start, y2_end);
  rdst->dy = y1_end - rdst->y + 1;
  if (rdst->dy < 0) rdst->dy = 0;

  z1_start = reg1->z;
  z1_end = reg1->z + reg1->dz - 1;
  z2_end = reg2->z + reg2->dz - 1;

  rdst->z = MAX(z1_start, z2_end);
  rdst->dz = z1_end - rdst->z + 1;
  if (rdst->dz < 0) rdst->dz = 0;

  return (rdst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONadd(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst) { return (rdst); }
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONclear(MRI_REGION *r)
{
  memset(r, 0, sizeof(*r));
  return (r);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONcopy(MRI_REGION *rsrc, MRI_REGION *rdst)
{
  if (!rdst) rdst = REGIONalloc();

  memmove(rdst, rsrc, sizeof(*rsrc));
  return (rdst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONintersect(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  int x2, y2, z2;

  rdst->x = MAX(reg1->x, reg2->x);
  rdst->y = MAX(reg1->y, reg2->y);
  rdst->z = MAX(reg1->z, reg2->z);

  x2 = MIN(reg1->x + reg1->dx, reg2->x + reg2->dx) - 1;
  y2 = MIN(reg1->y + reg1->dy, reg2->y + reg2->dy) - 1;
  z2 = MIN(reg1->z + reg1->dz, reg2->z + reg2->dz) - 1;
  rdst->dx = x2 - rdst->x + 1;
  rdst->dy = y2 - rdst->y + 1;
  rdst->dz = z2 - rdst->z + 1;
  return (rdst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONunion(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst) { return (rdst); }
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int REGIONinside(MRI_REGION *reg, int x, int y, int z)
{
  int x1, y1, z1;

  x1 = reg->x + reg->dx - 1;
  y1 = reg->y + reg->dy - 1;
  z1 = reg->z + reg->dz - 1;

  if (x < reg->x || x > x1 || y < reg->y || y > y1 || z < reg->z || z > z1) return (REGION_OUTSIDE);

  if ((x == reg->x || x == x1) && (y == reg->y || y == y1) && (z == reg->z || z == z1)) return (REGION_ON_BORDER);

  return (REGION_INSIDE);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *REGIONexpand(MRI_REGION *rsrc, MRI_REGION *rdst, int n)
{
  rdst->x = rsrc->x - n;
  rdst->y = rsrc->y - n;
  rdst->z = rsrc->z - n;
  rdst->dx = rsrc->dx + 2 * n;
  rdst->dy = rsrc->dy + 2 * n;
  rdst->dz = rsrc->dz + 2 * n;
  return (rdst);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float REGIONminCornerDistance(MRI_REGION *r1, MRI_REGION *r2)
{
  float min_dist = 10000.0f, dist, dx, dy, dz;
  int i, j, x0 = 0, y0 = 0, z0 = 0, x1 = 0, y1 = 0, z1 = 0;
  MRI_REGION r3;

  REGIONintersect(r1, r2, &r3);
  if (r3.dx > 0 && r3.dy > 0) return (0.0);

  for (i = 0; i < 8; i++) /* each corner of r1 */
  {
    regionCornerCoords(r1, i, &x0, &y0, &z0);
    for (j = 0; j < 8; j++) /* each corner of r2 */
    {
      regionCornerCoords(r2, j, &x1, &y1, &z1);
      dx = (float)(x1 - x0);
      dy = (float)(y1 - y0);
      dz = (float)(z1 - z0);
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (dist < min_dist) min_dist = dist;
    }
  }
  return (min_dist);
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int regionCornerCoords(MRI_REGION *r, int which_corner, int *px, int *py, int *pz)
{
  switch (which_corner) {
    case 0:
      *px = r->x;
      *py = r->y;
      *pz = r->z;
      break;
    case 1:
      *px = r->x + r->dx - 1;
      *py = r->y;
      *pz = r->z;
      break;
    case 2:
      *px = r->x;
      *py = r->y + r->dy - 1;
      *pz = r->z;
      break;
    case 3:
      *px = r->x;
      *py = r->y;
      *pz = r->z + r->dz - 1;
      break;
    case 4:
      *px = r->x + r->dx - 1;
      *py = r->y + r->dy - 1;
      *pz = r->z;
      break;
    case 5:
      *px = r->x + r->dx - 1;
      *py = r->y;
      *pz = r->z + r->dz - 1;
      break;
    case 6:
      *px = r->x;
      *py = r->y + r->dy - 1;
      *pz = r->z + r->dz - 1;
      break;
    case 7:
      *px = r->x + r->dx - 1;
      *py = r->y + r->dy - 1;
      *pz = r->z + r->dz - 1;
      break;
  }
  return (NO_ERROR);
}
/*!
  \fn MRI_REGION *REGIONgetBoundingBoxM(MRI *mask, int npad[6])
  \brief Determines bounding box as corners of the smallest box needed
  to fit all the non-zero voxels. npad[X] allows the BB to be expanded
  in each direction. It will expand the first dimension by npad[0]
  toward 0.  It will expand the first dimension by npad[1] toward inf.
  Etc, for the remaining dims. region->{x,y,z} is the CRS 0-based
  starting point of the box.  region->{dx,dy,dz} is the size of the
  box such that a loop would run for(c=region->x; c <
  region->x+region->dx; c++)
*/
MRI_REGION *REGIONgetBoundingBoxM(const MRI *mask, const int npad[6])
{
  int c, r, s;
  int cmin, cmax, rmin, rmax, smin, smax;
  MRI_REGION *region;
  float v;

  cmin = rmin = smin = 1000000;
  cmax = rmax = smax = 0;

  for (c = 0; c < mask->width; c++) {
    for (r = 0; r < mask->height; r++) {
      for (s = 0; s < mask->depth; s++) {
        v = MRIgetVoxVal(mask, c, r, s, 0);
        if (iszero(v)) continue;
        if (cmin > c) cmin = c;
        if (rmin > r) rmin = r;
        if (smin > s) smin = s;
        if (cmax < c) cmax = c;
        if (rmax < r) rmax = r;
        if (smax < s) smax = s;
      }
    }
  }
  region = REGIONalloc();
  region->x = MAX(cmin - npad[0], 0);
  region->y = MAX(rmin - npad[2], 0);
  region->z = MAX(smin - npad[4], 0);
  region->dx = MIN(cmax - cmin + npad[0] + npad[1], mask->width  - region->x);
  region->dy = MIN(rmax - rmin + npad[2] + npad[3], mask->height - region->y);
  region->dz = MIN(smax - smin + npad[4] + npad[5], mask->depth  - region->z);

  return (region);
}

/*!
  \fn MRI_REGION *REGIONgetBoundingBox(MRI *mask, int npad)
  \brief Determines bounding box as corners of the smallest box needed
  to fit all the non-zero voxels. If npad is non-zero then then the box
  is expanded by npad in each direction (making it 2*npad bigger in each
  dimension). region->{x,y,z} is the CRS 0-based starting point of the box.
  region->{dx,dy,dz} is the size of the box such that a loop would run
  for(c=region->x; c < region->x+region->dx; c++). Note: should change this
  to use REGIONgetBoundingBoxM()
*/
MRI_REGION *REGIONgetBoundingBox(MRI *mask, int npad)
{
  int c, r, s;
  int cmin, cmax, rmin, rmax, smin, smax;
  MRI_REGION *region;
  float v;

  //Note: should change this to use REGIONgetBoundingBoxM()

  cmin = rmin = smin = 1000000;
  cmax = rmax = smax = 0;

  for (c = 0; c < mask->width; c++) {
    for (r = 0; r < mask->height; r++) {
      for (s = 0; s < mask->depth; s++) {
        v = MRIgetVoxVal(mask, c, r, s, 0);
        if (iszero(v)) continue;
        if (cmin > c) cmin = c;
        if (rmin > r) rmin = r;
        if (smin > s) smin = s;
        if (cmax < c) cmax = c;
        if (rmax < r) rmax = r;
        if (smax < s) smax = s;
      }
    }
  }
  region = REGIONalloc();
  region->x = MAX(cmin - npad, 0);
  region->y = MAX(rmin - npad, 0);
  region->z = MAX(smin - npad, 0);
  region->dx = MIN(cmax - cmin + 2 * npad, mask->width - region->x);
  region->dy = MIN(rmax - rmin + 2 * npad, mask->height - region->y);
  region->dz = MIN(smax - smin + 2 * npad, mask->depth - region->z);

  return (region);
}

/*!
  \fn MRI_REGION *REGIONgetBoundingBoxEqOdd(MRI *mask, int npad)
  \brief Same as REGIONgetBoundingBox() but forces the output to
  have nrows=ncols and nslices to be odd. This is a customization
  to accomodate the STIR PET simulator
*/
MRI_REGION *REGIONgetBoundingBoxEqOdd(MRI *mask, int npad)
{
  MRI_REGION *region;
  int dxy,delta, delta1, delta2, isodd;

  region = REGIONgetBoundingBox(mask, npad);
  isodd = (region->dz % 2);

  printf("Volume size: %d %d %d, OrigRegion: ",mask->width,mask->height,mask->depth);
  REGIONprint(stdout,region);
 
  if(region->dx > mask->height){
    printf("ERROR: REGIONgetBoundingBoxEqOdd(): dx is bigger than MRI height\n");
    return(NULL);
  }
  if(region->dy > mask->width){
    printf("ERROR: REGIONgetBoundingBoxEqOdd(): dy is bigger than MRI width\n");
    return(NULL);
  }

  if(region->dx == region->dy && isodd) 
    return(region);

  if(!isodd){
    if(region->dz < mask->depth){
      // Add 1 to dz to make it odd
      region->dz += 1; 
    }
    else if(region->z > 0){
      // Remove 1 from z and add 1 to dz to make it odd
      region->z -= 1;
      region->dz += 1; 
    }
    else {
      printf("ERROR: REGIONgetBoundingBoxEqOdd(): cannot make z odd: z=%d, dz=%d, nz=%d\n",
	     region->z,region->dz,mask->depth);
      return(NULL);
    }
  }
  if(region->dx == region->dy)
    return(region);

  // If it gets here, then dx != dy
  dxy = MAX(region->dx,region->dy);
  if(region->x + dxy <= mask->width && region->y + dxy <= mask->height){
    region->dx = dxy;
    region->dy = dxy;
    return(region);
  }
  // Divide the difference into two (maybe equal) parts
  delta = fabs(region->dx-region->dy);
  delta1 = delta/2;
  delta2 = delta - delta1;
  if(region->dx < region->dy){
    // reduce the x start by half the difference
    if(region->x - delta1 > 0 && region->x + region->dx + delta2 <= mask->width){
      region->x  -= delta1;
      region->dx += delta;
      return(region);
    }
    // Well, that did not work, use asymmetric padding
    int pad2 = mask->width - (region->x + region->dx);
    if(region->x < pad2){
      // Add as much as possible to the left, then pad to the right
      region->x = 0;
      region->dx = region->dy;
    }
    else {
      // Add as much as possible to the right, then pad to the left
      region->dx = region->dy;
      region->x = mask->width - region->dx;
    }
    printf("Volume size: %d %d %d, xAsymRegion: ",mask->width,mask->height,mask->depth);
    REGIONprint(stdout,region);
    return(region);
  }
  if(region->dy < region->dx){
    // reduce the y start by half the difference
    if(region->y - delta1 > 0 && region->y + region->dy + delta2 <= mask->height){
      region->y  -= delta1;
      region->dy += delta;
      return(region);
    }
    int pad2 = mask->height - (region->y + region->dy);
    if(region->y < pad2){
      region->y = 0;
      region->dy = region->dx;
    }
    else {
      region->dy = region->dx;
      region->y = mask->height - region->dy;
    }
    printf("Volume size: %d %d %d, yAsymRegion: ",mask->width,mask->height,mask->depth);
    REGIONprint(stdout,region);
    return(region);
  }

  // Should probably never get here
  printf("ERROR: REGIONgetBoundingBoxEqOdd(): cannot make x=y\n");
  printf("Volume size: %d %d %d, Region: ",mask->width,mask->height,mask->depth);
  REGIONprint(stdout,region);
  return(NULL);
}

/*!
  \fn int REGIONprint(FILE *fp, MRI_REGION *r)
  \brief Prints REGION struct.
*/
int REGIONprint(FILE *fp, MRI_REGION *r)
{
  fprintf(fp, "%d %d %d  %d %d %d\n", r->x, r->y, r->z, r->dx, r->dy, r->dz);
  return (0);
}
