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
                           CONSTANTS
----------------------------------------------------------------------*/

#define MAX_IMAGES  30

/*----------------------------------------------------------------------
                           GLOBAL DATA
----------------------------------------------------------------------*/

IMAGE          *Idisplay[MAX_IMAGES] = { NULL } ;
MRI            *mris[MAX_IMAGES] ;
int            mri_views[MAX_IMAGES] ;
int            mri_depths[MAX_IMAGES] ;
int            mri_frames[MAX_IMAGES] ;
int            which_click = -1 ;

/*----------------------------------------------------------------------
                           STATIC DATA
----------------------------------------------------------------------*/

static XV_FRAME       *xvf ;
static Menu           view_menu ;
static char           view_str[100] ;
static Panel_item     view_panel ;

static int            x_click ;
static int            y_click ;
static int            z_click ;

static int            talairach = 0 ; /* show image or Talairach coords */

/*----------------------------------------------------------------------
                        STATIC PROTOTYPES
----------------------------------------------------------------------*/

static void viewMenuItem(Menu menu, Menu_item menu_item) ;
static IMAGE *get_next_slice(IMAGE *Iold, int which, int dir) ;
static void repaint_handler(XV_FRAME *xvf, DIMAGE *dimage) ;

/*----------------------------------------------------------------------
                              FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void 
mri_event_handler(XV_FRAME *xvf, Event *event,DIMAGE *dimage, 
                  int *px, int *py, int *pz)
{
  int       x, y, z, which, depth, view ;
  Real      xr, yr, zr, xt, yt, zt ;
  HISTOGRAM *histo ;
  float     fmin, fmax ;
  XV_FRAME  *xvnew ;
  MRI       *mri ;

  which = dimage->which ;
  mri = mris[which] ;
  depth = mri_depths[which] ;
  view = mri_views[which] ;

  /* click can occur in the middle of other stuff (sort of asynchonous) */
  if (!mri || !mri->slices)  
    return ;

/*
  The convention  is  that  positive xspace coordinates run
  from the patient's  left  side  to  right  side,  positive
  yspace  coordinates run from patient posterior to anterior
  and positive zspace coordinates run from inferior to superior.
 */   
  switch (mri_views[which])
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
    y = depth ;
    if (xvf->ydir < 0)
      z = mri->height - (event_y(event)+1) ;
    else
      z = event_y(event) ;
    break ;
  case MRI_SAGITAL:
    x = depth ;
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

  if (talairach)
    MRIvoxelToTalairachVoxel(mri, (Real)x, (Real)y, (Real)z, &xt, &yt, &zt) ;
  MRIvoxelToWorld(mri, (Real)x, (Real)y, (Real)z, &xr, &yr, &zr) ;
  switch (event_id(event))
  {
  case LOC_DRAG:
    if (!event_left_is_down(event))
      return ;
  case MS_LEFT:
    switch (mri->type)
    {
    case MRI_UCHAR:
      if (talairach)
        XVprintf(xvf, 0, "T: (%d,%d,%d) --> %d",
                 nint(xt),nint(yt),nint(zt),MRIvox(mri, x, y, z));
      else
        XVprintf(xvf, 0, "(%d,%d,%d) --> %d",x,y,z,MRIvox(mri, x, y, z));
      break ;
    case MRI_FLOAT:
      if (talairach)
        XVprintf(xvf, 0, "T: (%d,%d,%d) --> %2.3f",
                 nint(xt),nint(yt),nint(zt),MRIFvox(mri, x, y, z));
      else
        XVprintf(xvf, 0, "(%d,%d,%d) --> %2.3f",x,y,z,MRIFvox(mri, x, y, z));
      break ;
    }
    break ;
  default:
    if (event_is_up(event)) switch ((char)event->ie_code)
    {
    case 'T':
      talairach = 1 ;
      break ;
    case 't':
      talairach = 0 ;
      break ;
    case 'h':
      MRIvalRange(mri, &fmin, &fmax) ;
      histo = MRIhistogram(mri, (int)(fmax - fmin + 1)) ;
      xvnew = XValloc(1, 1, 2, 200, 200, "histogram tool", NULL) ;
      XVshowHistogram(xvnew, 0, histo) ;
      xv_main_loop(xvnew->frame);
      HISTOfree(&histo) ;
      break ;
    }
    break ;
  }

  if (event_is_up(event) && (event_id(event) == MS_LEFT))
  {
    int view, old_which ;

    if (which_click != which)   /* erase old point */
    {
      old_which = which_click ;
      which_click = dimage->which ;
      XVrepaintImage(xvf, old_which) ;
      view = mri_views[which] ;
      sprintf(view_str, "view: %s", 
              view == MRI_CORONAL ? "CORONAL" :
              view == MRI_SAGITAL ? "SAGITAL" : "HORIZONTAL") ;
      xv_set(view_panel, PANEL_LABEL_STRING, view_str, NULL) ;
    }

    x_click = x ;
    y_click = y ;
    z_click = z ;
    XVrepaintImage(xvf, which_click) ;
  }
  if (px)
    *px = x ;
  if (py)
    *py = y ;
  if (pz)
    *pz = z ;

}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVMRIdrawPoint(XV_FRAME *xvf, int which, int view, int depth, MRI *mri,
               int x,int y,int z,int color)
{
  int xi, yi ;

  switch (mri_views[which])
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
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
XVMRIshowFrame(XV_FRAME *xvf, MRI *mri, int which, int slice,int frame)
{
  IMAGE  *I ;
  float  mag ;

  if (frame < 0)
    frame = 0 ;

  if (which_click == which && mri != mris[which])
    which_click = -1 ;  /* a new MR image shown */

  if (!mri)
    return(NULL) ;

  mri_frames[which] = frame ;

  if (slice < 0)  /* set slice to middle of slice direction */
  {
    switch (mri_views[which])
    {
    case MRI_CORONAL:
      slice = (mri->imnr0 + mri->imnr1) / 2 ;
      break ;
    case MRI_SAGITAL:
      slice = mri->width / 2 ;
      break ;
    case MRI_HORIZONTAL:
      slice = mri->height / 2 ;
      break ;
    }
  }

  mri_depths[which] = slice ;
  switch (mri_views[which])
  {
  case MRI_CORONAL:
    slice -= mri->imnr0 ;
    break ;
  case MRI_SAGITAL:
  case MRI_HORIZONTAL:
    break ;
  }

  
  I = MRItoImageView(mri, Idisplay[which], slice, mri_views[which], frame) ;
  if (!I)
    return(NULL) ;


  mag = MIN((float)xvf->orig_disp_rows / (float)I->rows,
            (float)xvf->orig_disp_cols / (float)I->cols) ;

  XVsetImageSize(xvf, which, nint((float)I->rows*mag), 
                 nint((float)I->cols*mag));
  XVresize(xvf) ;

  /* must be done before XVshowImage to draw point properly */
  if (which_click < 0)  /* reset current click point */
  {
    which_click = which ;
    z_click = slice ;
    y_click = mri->height / 2 ;
    x_click = mri->width / 2 ;
#if 0
    XVMRIdrawPoint(xvf, which, mri_views[which], 0, mri, x_click,
                   y_click, z_click, XXOR);
#endif
  }
  XVshowImage(xvf, which, I, 0) ;


  if (Idisplay[which] && (I != Idisplay[which]))
    ImageFree(&Idisplay[which]) ;

  Idisplay[which] = I ;

  mris[which] = mri ;
  return(I) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
XVMRIshow(XV_FRAME *xvf, MRI *mri, int which, int slice)
{
  return(XVMRIshowFrame(xvf, mri, which, slice, 0)) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int  
XVMRIinit(XV_FRAME *xvf_init, int view_row, int view_col)
{
  int i ;

  for (i = 0 ; i < MAX_IMAGES ; i++)
    mri_views[i] = MRI_CORONAL ;

  xvf = xvf_init ;
  XVsetDepthFunc(xvf, get_next_slice) ;
  XVsetPrintStatus(xvf, 0) ;
  XVsetYDir(xvf, 1) ;
  sprintf(view_str, "view: CORONAL") ;
  view_menu = (Menu)
    xv_create((Xv_opaque)NULL, MENU,
              MENU_NOTIFY_PROC,    viewMenuItem,
              MENU_STRINGS,        "CORONAL", "SAGITAL", "HORIZONTAL", NULL,
              NULL) ;

  view_panel = (Panel_item)
    xv_create(xvf->panel, PANEL_BUTTON,
              PANEL_LABEL_STRING,    view_str,
              XV_X,                  view_col,
              XV_Y,                  view_row,
              PANEL_ITEM_MENU,       view_menu,
              NULL) ;

  XVsetRepaintHandler(repaint_handler) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void 
viewMenuItem(Menu menu, Menu_item menu_item)
{
  char   *menu_str ;
  MRI    *mri, *mri2 ;
  int    slice, which, view, which2, sync ;
  DIMAGE *dimage, *dimage2 ;

  which = which_click ;
  if (which_click < 0)   /* no current window */
    which = 0 ;

  mri = mris[which] ;
  if (!mri)              /* no mri in current window */
    return ;

  menu_str = (char *)xv_get(menu_item, MENU_STRING) ;

  if (!stricmp(menu_str, "CORONAL"))
  {
    view = MRI_CORONAL ;
    slice = z_click + mri->imnr0 ;
  }
  else if (!stricmp(menu_str, "SAGITAL"))
  {
    view = MRI_SAGITAL ;
    slice = x_click ;
  }
  else
  {
    view = MRI_HORIZONTAL ;
    slice = y_click ;
  }

  sprintf(view_str, "view: %s", menu_str) ;
  xv_set(view_panel, PANEL_LABEL_STRING, view_str, NULL) ;
  XVMRIsetView(xvf, which, view) ;

  dimage = XVgetDimage(xvf, which, DIMAGE_IMAGE) ;
  sync = dimage->sync ;
  if (sync)
  {
    for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
    {
      if (which2 == which)
        continue ;
      dimage2 = XVgetDimage(xvf, which2, DIMAGE_IMAGE) ;
      mri2 = mris[which2] ;
      if (dimage2 && mri2 && (dimage2->sync == sync))
      {
        XVMRIsetView(xvf, which2, view) ;
        XVMRIshowFrame(xvf, mri2, which2, slice, mri_frames[which2]) ;
      }
    }
  }

  /* will reset syncs */
  XVMRIshowFrame(xvf, mri, which, slice, mri_frames[which]) ;  
  if (sync)  /* if they were synced, reinstate it */
    XVsyncAll(xvf, which) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static IMAGE *
get_next_slice(IMAGE *Iold, int which, int dir)
{
  IMAGE *I ;
  MRI   *mri ;
  int   depth, offset ;

  if (!Iold)
    return(NULL) ;

  mri = mris[which] ;
  if (mri)
  {
    depth = mri_depths[which] + dir ;
    if (mri_views[which] == MRI_CORONAL)
      offset = mri->imnr0 ;
    else
      offset = 0 ;

    I = MRItoImageView(mri, Idisplay[which], depth-offset, mri_views[which],
                       mri_frames[which]);
    
    if (!I)
      I = Idisplay[which] ;         /* failed to get next slice */
    else
    {
      mri_depths[which] = depth ;   /* new depth is valid */
      if (Idisplay[which] && (Idisplay[which] != I))
        ImageFree(&Idisplay[which]) ;
      Idisplay[which] = I ;
    }
  }
  else
    I = Iold ;

  return(I) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVMRIshowAll(XV_FRAME *xvf)
{
  int i ;

  for (i = 0 ; i < MAX_IMAGES ; i++)
    XVMRIshowFrame(xvf, mris[i], i, mri_depths[i], mri_frames[i]) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
repaint_handler(XV_FRAME *xvf, DIMAGE *dimage)
{
  if (dimage->which == which_click && (mris[which_click] != NULL))
    XVMRIdrawPoint(xvf, which_click, mri_views[which_click], 0, 
                   mris[which_click], x_click, y_click, z_click, XRED) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVMRIfree(MRI **pmri, int which)
{
  MRI *mri ;

  mri = *pmri ;
  mris[which] = NULL ;
  MRIfree(pmri) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVMRIsetView(XV_FRAME *xvf, int which, int view)
{
  mri_views[which] = view ;
  return(NO_ERROR) ;
}
