/*----------------------------------------------------------------------
           File Name:  xvmri.h

           Author:

           Description:

           Conventions:

----------------------------------------------------------------------*/
/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>

#include "xvutil.h"
#include "mri.h"
#include "proto.h"
#include "xvmri.h"
#include "error.h"
#include "diag.h"
#include "histo.h"

/*----------------------------------------------------------------------
                           GLOBAL DATA
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                              FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void 
mri_event_handler(XV_FRAME *xvf,int depth,MRI *mri,Event *event,DIMAGE *dimage,
                  int view, int *px, int *py, int *pz)
{
  int       x, y, z ;
  Real      xr, yr, zr ;
  HISTOGRAM *histo ;
  float     fmin, fmax ;
  XV_FRAME  *xvnew ;

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (view)
  {
  default:
  case MRI_CORONAL:
    /*
      Z=0 is the back of the head,
      X=0 is the right side of the head
      Y=0 is the neck/brain stem
      */
    x = event_x(event) ;
    if (xvf->ydir < 0)
      y = mri->height - (event_y(event)+1) ;
    else
      y = event_y(event) ;
    z = depth - mri->imnr0 ;
    break ;
  case MRI_HORIZONTAL:
    x = event_x(event) ;
    y = depth - mri->imnr0 ;
    if (xvf->ydir < 0)
      z = mri->height - (event_y(event)+1) ;
    else
      z = event_y(event) ;
    break ;
  case MRI_SAGITAL:
    x = depth - mri->imnr0 ;
    if (xvf->ydir < 0)
      y = mri->height - (event_y(event)+1) ;
    else
      y = event_y(event) ;
    z = event_x(event) ;
    break ;
  }

  if (x <= 0)
    x = 0 ;
  else if (x >= mri->width)
    x = mri->width-1 ;
  if (y <= 0)
    y = 0 ;
  else if (y >= mri->height)
    y = mri->height-1 ;
  if (z <= 0)
    z = 0 ;
  else if (z >= mri->depth)
    z = mri->depth-1 ;

  MRIvoxelToWorld(mri, (Real)x, (Real)y, (Real)z, &xr, &yr, &zr) ;
  switch (event_id(event))
  {
  case LOC_DRAG:
    if (!event_left_is_down(event))
      return ;
  case MS_LEFT:
    XVprintf(xvf, 0, "(%d,%d,%d) --> %d",x,y,z,MRIvox(mri, x, y, z));
    break ;
  default:
    switch ((char)event->ie_code)
    {
    case 'h':
      if (event_is_up(event))
      {
        MRIvalRange(mri, &fmin, &fmax) ;
        histo = MRIhistogram(mri, (int)(fmax - fmin + 1)) ;
        xvnew = XValloc(1, 1, 2, 200, 200, "histogram tool", NULL) ;
        XVshowHistogram(xvnew, 0, histo) ;
        xv_main_loop(xvnew->frame);
        HISTOfree(&histo) ;
        break ;
      }
    }
    break ;
  }
  if (px)
    *px = x ;
  if (py)
    *py = y ;
  if (pz)
    *pz = z ;
}
void
XVMRIdrawPoint(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
               int x,int y,int z,int color)
{
  int xi, yi ;

  switch (view)
  {
  default:
  case MRI_CORONAL:
    /*
      Z=0 is the back of the head,
      X=0 is the right side of the head
      Y=0 is the neck/brain stem
      */
    xi = x ;
    yi = y ;
    break ;
  case MRI_HORIZONTAL:
    xi = x ;
    yi = z ;
    break ;
  case MRI_SAGITAL:
    xi = z ;
    yi = y ;
    break ;
  }
  XVdrawPoint(xvf, which, xi, yi, color) ;
}

void
XVMRIdrawRegion(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
                MRI_REGION *reg, int color)
{
  int xi, yi, dx, dy ;

  switch (view)
  {
  default:
  case MRI_CORONAL:
    /*
      Z=0 is the back of the head,
      X=0 is the right side of the head
      Y=0 is the neck/brain stem
      */
    xi = reg->x ;
    yi = reg->y ;
    dx = reg->dx ;
    dy = reg->dy ;
    break ;
  case MRI_HORIZONTAL:
    xi = reg->x ;
    yi = reg->z ;
    dx = reg->dx ;
    dy = reg->dz ;
    break ;
  case MRI_SAGITAL:
    xi = reg->z ;
    dx = reg->dz ;
    yi = reg->y ;
    dy = reg->dy ;
    break ;
  }
  XVdrawBox(xvf, which, xi, yi, dx, dy, color) ;
}

