/*
 *       FILE NAME:   region.c
 *
 *       DESCRIPTION: utilities for the REGION  data structure
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        1/8/97
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

  memcpy(rdst, rsrc, sizeof(*rsrc)) ;
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

