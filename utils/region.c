/**
 * @file  region.c
 * @brief utilities for the REGION data structure
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/12/07 21:53:31 $
 *    $Revision: 1.9 $
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

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>

#include "mri.h"
#include "macros.h"
#include "const.h"
#include "region.h"
#include "diag.h"
#include "error.h"

/*-----------------------------------------------------
                    MACROS AND CONSTANTS
-------------------------------------------------------*/

/*-----------------------------------------------------
                    STATIC PROTOTYPES
-------------------------------------------------------*/

static int regionCornerCoords(MRI_REGION *r, int which_corner, int *px,
                              int *py, int *pz) ;

/*-----------------------------------------------------
                    GLOBAL FUNCTIONS
-------------------------------------------------------*/
MRI_REGION *
REGIONalloc(void)
{
  MRI_REGION *r ;

  r = (MRI_REGION *)calloc(1, sizeof(MRI_REGION)) ;
  if (!r)
    ErrorExit(ERROR_NO_MEMORY, "REGIONalloc: could not allocate region") ;

  return(r) ;
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
MRI_REGION *
REGIONsubtract(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  int   x1_start, x1_end, x2_end, y1_start, y1_end, y2_end,
  z1_start, z1_end, z2_end ;

  x1_start = reg1->x ;
  x1_end = reg1->x + reg1->dx - 1 ;
  x2_end = reg2->x + reg2->dx - 1 ;

  rdst->x = MAX(x1_start, x2_end) ;
  rdst->dx = x1_end - rdst->x + 1 ;
  if (rdst->dx < 0)
    rdst->dx = 0 ;

  y1_start = reg1->y ;
  y1_end = reg1->y + reg1->dy - 1 ;
  y2_end = reg2->y + reg2->dy - 1 ;

  rdst->y = MAX(y1_start, y2_end) ;
  rdst->dy = y1_end - rdst->y + 1 ;
  if (rdst->dy < 0)
    rdst->dy = 0 ;

  z1_start = reg1->z ;
  z1_end = reg1->z + reg1->dz - 1 ;
  z2_end = reg2->z + reg2->dz - 1 ;

  rdst->z = MAX(z1_start, z2_end) ;
  rdst->dz = z1_end - rdst->z + 1 ;
  if (rdst->dz < 0)
    rdst->dz = 0 ;

  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONadd(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONclear(MRI_REGION *r)
{
  memset(r, 0, sizeof(*r)) ;
  return(r) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONcopy(MRI_REGION *rsrc, MRI_REGION *rdst)
{
  if (!rdst)
    rdst = REGIONalloc() ;

  memmove(rdst, rsrc, sizeof(*rsrc)) ;
  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONintersect(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  int x2, y2, z2 ;

  rdst->x = MAX(reg1->x, reg2->x) ;
  rdst->y = MAX(reg1->y, reg2->y) ;
  rdst->z = MAX(reg1->z, reg2->z) ;

  x2 = MIN(reg1->x+reg1->dx, reg2->x+reg2->dx) - 1 ;
  y2 = MIN(reg1->y+reg1->dy, reg2->y+reg2->dy) - 1 ;
  z2 = MIN(reg1->z+reg1->dz, reg2->z+reg2->dz) - 1 ;
  rdst->dx = x2 - rdst->x + 1 ;
  rdst->dy = y2 - rdst->y + 1 ;
  rdst->dz = z2 - rdst->z + 1 ;
  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONunion(MRI_REGION *reg1, MRI_REGION *reg2, MRI_REGION *rdst)
{
  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
int
REGIONinside(MRI_REGION *reg, int x, int y, int z)
{
  int x1, y1, z1 ;

  x1 = reg->x + reg->dx - 1 ;
  y1 = reg->y + reg->dy - 1 ;
  z1 = reg->z + reg->dz - 1 ;

  if (x < reg->x || x > x1 || y < reg->y || y > y1 || z < reg->z || z > z1)
    return(REGION_OUTSIDE) ;

  if ((x == reg->x || x == x1) &&
      (y == reg->y || y == y1) &&
      (z == reg->z || z == z1))
    return(REGION_ON_BORDER) ;

  return(REGION_INSIDE) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
MRI_REGION *
REGIONexpand(MRI_REGION *rsrc, MRI_REGION *rdst, int n)
{
  rdst->x = rsrc->x - n ;
  rdst->y = rsrc->y - n ;
  rdst->z = rsrc->z - n ;
  rdst->dx = rsrc->dx + 2*n ;
  rdst->dy = rsrc->dy + 2*n ;
  rdst->dz = rsrc->dz + 2*n ;
  return(rdst) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
float
REGIONminCornerDistance(MRI_REGION *r1, MRI_REGION *r2)
{
  float       min_dist = 10000.0f, dist, dx, dy, dz ;
  int         i, j, x0=0, y0=0, z0=0, x1=0, y1=0, z1=0 ;
  MRI_REGION  r3 ;

  REGIONintersect(r1, r2, &r3) ;
  if (r3.dx > 0 && r3.dy > 0)
    return(0.0) ;

  for (i = 0 ; i < 8 ; i++)   /* each corner of r1 */
  {
    regionCornerCoords(r1, i, &x0, &y0, &z0) ;
    for (j = 0 ; j < 8 ; j++)   /* each corner of r2 */
    {
      regionCornerCoords(r2, j, &x1, &y1, &z1) ;
      dx = (float)(x1 - x0) ;
      dy = (float)(y1 - y0) ;
      dz = (float)(z1 - z0) ;
      dist = sqrt(dx*dx + dy*dy + dz*dz) ;
      if (dist < min_dist)
        min_dist = dist ;
    }
  }
  return(min_dist) ;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
------------------------------------------------------*/
static int
regionCornerCoords(MRI_REGION *r, int which_corner, int *px, int *py, int *pz)
{
  switch (which_corner)
  {
  case 0:
    *px = r->x ;
    *py = r->y ;
    *pz = r->z ;
    break ;
  case 1:
    *px = r->x + r->dx - 1 ;
    *py = r->y ;
    *pz = r->z ;
    break ;
  case 2:
    *px = r->x ;
    *py = r->y + r->dy - 1 ;
    *pz = r->z ;
    break ;
  case 3:
    *px = r->x ;
    *py = r->y ;
    *pz = r->z + r->dz - 1 ;
    break ;
  case 4:
    *px = r->x + r->dx - 1 ;
    *py = r->y + r->dy - 1 ;
    *pz = r->z ;
    break ;
  case 5:
    *px = r->x + r->dx - 1 ;
    *py = r->y ;
    *pz = r->z + r->dz - 1 ;
    break ;
  case 6:
    *px = r->x ;
    *py = r->y + r->dy - 1 ;
    *pz = r->z + r->dz - 1 ;
    break ;
  case 7:
    *px = r->x + r->dx - 1 ;
    *py = r->y + r->dy - 1 ;
    *pz = r->z + r->dz - 1 ;
    break ;
  }
  return(NO_ERROR) ;
}
/*!
  \fn MRI_REGION *REGIONgetBoundingBox(MRI *mask, int npad)
  \brief Determines bounding box as corners of the smallest box needed
  to fit all the non-zero voxels. If npad is non-zero then then the box
  is expanded by npad in each direction (making it 2*npad bigger in each
  dimension).
*/
MRI_REGION *REGIONgetBoundingBox(MRI *mask, int npad)
{
  int c, r, s;
  int cmin, cmax, rmin, rmax, smin, smax;
  MRI_REGION *region;
  float v;
  
  cmin = rmin = smin = 1000000;
  cmax = rmax = smax = 0;

  for(c=0; c < mask->width; c++){
    for(r=0; r < mask->height; r++){
      for(s=0; s < mask->depth; s++){
	v = MRIgetVoxVal(mask,c,r,s,0);
	if(iszero(v)) continue;
	if(cmin > c) cmin = c;
	if(rmin > r) rmin = r;
	if(smin > s) smin = s;
	if(cmax < c) cmax = c;
	if(rmax < r) rmax = r;
	if(smax < s) smax = s;
      }
    }
  }
  region = REGIONalloc();
  region->x  = cmin - npad;
  region->y  = rmin - npad;
  region->z  = smin - npad;
  region->dx = cmax-cmin + 2*npad;
  region->dy = rmax-rmin + 2*npad;
  region->dz = smax-smin + 2*npad;

  return(region);
}

/*!
  \fn int REGIONprint(FILE *fp, MRI_REGION *r)
  \brief Prints REGION struct.
*/
int REGIONprint(FILE *fp, MRI_REGION *r)
{
  fprintf(fp,"%d %d %d  %d %d %d\n",r->x,r->y,r->z,r->dx,r->dy,r->dz);
  return(0);
}
