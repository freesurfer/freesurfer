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
                  int view)
{
  int       x, y, z ;
  Real      xr, yr, zr ;
  MRI_HISTO *mrih ;
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
    y = mri->height - (event_y(event)+1) ;
    z = depth - mri->imnr0 ;
    break ;
  case MRI_HORIZONTAL:
    x = event_x(event) ;
    y = depth - mri->imnr0 ;
    z = mri->height - (event_y(event)+1) ;
    break ;
  case MRI_SAGITAL:
    x = depth - mri->imnr0 ;
    y = mri->height - (event_y(event)+1) ;
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
        mrih = MRIhistogram(mri, (int)(fmax - fmin + 1)) ;
        xvnew = XValloc(1, 1, 2, 200, 200, "histogram tool", NULL) ;
        XVMRIshowHistogram(xvnew, 0, mrih) ;
        xv_main_loop(xvnew->frame);
        /*        MRIdumpHistogram(mrih, stderr) ;*/
        MRIfreeHistogram(&mrih) ;
        break ;
      }
    }
    break ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVMRIshowHistogram(XV_FRAME *xvf, int which, MRI_HISTO *mrih)
{
  int  binno, max_count, nbins ;

  nbins = mrih->nbins ;
  max_count = 0 ;
  for (binno = 1 ; binno < nbins ; binno++)  /* ignore 0th bin */
    if (mrih->counts[binno] > max_count)
      max_count = mrih->counts[binno] ;

  fprintf(stderr, "max_count = %d\n", max_count) ;

  return(NO_ERROR) ;
}
