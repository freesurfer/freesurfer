/**
 * @file  xvmri.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.42 $
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
#include "utils.h"
#include "mri_circulars.h"

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
int            mri_slices[MAX_IMAGES] ;
int            x_click ;
int            y_click ;
int            z_click ;

int            brush = 0 ;


/*----------------------------------------------------------------------
                           STATIC DATA
----------------------------------------------------------------------*/

static XV_FRAME       *xvf ;
static Menu           view_menu ;
static char           view_str[100] ;
static Panel_item     view_panel ;

static char           image_path[100] = "." ;

static int            talairach = 0 ; /* show image or Talairach coords */
static int            fixed [MAX_IMAGES] = { 0 } ;
static MRI_SURFACE    *mri_surface = NULL ;
static int            dont_redraw = 0 ;
static int            show_three_views = 0 ;

/*----------------------------------------------------------------------
                        STATIC PROTOTYPES
----------------------------------------------------------------------*/

static void viewMenuItem(Menu menu, Menu_item menu_item) ;
static IMAGE *get_next_slice(IMAGE *Iold, int which, int dir) ;
static void repaint_handler(XV_FRAME *xvf, DIMAGE *dimage) ;
static int mri_write_func(Event *event, DIMAGE *dimage, char *fname) ;
#if 0
static int xvmriRepaintValue(XV_FRAME *xvf, int which, int x, int y, int z) ;
#endif

/*----------------------------------------------------------------------
                              FUNCTIONS
----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
mri_event_handler(XV_FRAME *xvf, Event *event,DIMAGE *dimage,
                  int *px, int *py, int *pz)
{
  int       x, y, z, which, depth, view, frame, xi, yi, zi, xk, yk, zk, slice ;
  Real      xr, yr, zr, xt, yt, zt, xv, yv, zv, xtv, ytv, ztv ;
  float     xf, yf, zf, xft, yft, zft ;
  MRI       *mri ;
  char      fname[100] ;
  FILE      *fp ;
  BUFTYPE   val, old_val ;
  static int repaint_needed = 0 ;

  which = dimage->which ;
  mri = mris[which] ;
  depth = mri_depths[which] ;
  frame = mri_frames[which] ;
  view = mri_views[which] ;
  slice = mri_slices[which] ;

  /* click can occur in the middle of other stuff (sort of asynchonous) */
  if (!mri || !mri->slices)
    return(ERROR_BAD_PARM) ;

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
  case MRI_SAGITTAL:
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

  if (event_is_up(event) && repaint_needed)
  {
    repaint_needed = 0 ;
    XVMRIredisplayFrame(xvf, mri, which, mri_depths[which], mri_frames[which]);
  }
  if (talairach)
  {
    MRIvoxelToTalairach(mri, (Real)x, (Real)y, (Real)z, &xt, &yt, &zt) ;
    MRIvoxelToTalairachVoxel(mri, (Real)x, (Real)y, (Real)z, &xtv,&ytv,&ztv);
  }
  MRIvoxelToWorld(mri, (Real)x, (Real)y, (Real)z, &xr, &yr, &zr) ;

  if (px)
    *px = x ;
  if (py)
    *py = y ;
  if (pz)
    *pz = z ;

  if (event_shift_is_down(event) && !event_is_ascii(event))
  {
    if (event_right_is_down(event))
      val = 255 ;
    else if (event_middle_is_down(event))
      val = 1 ;
    else
      val = 0 ;

    for (zk = z - brush ; zk <= z+brush ; zk++)
    {
      zi = mri->zi[zk] ;
      for (yk = y - brush ; yk <= y+brush ; yk++)
      {
        yi = mri->yi[yk] ;
        for (xk = x - brush ; xk <= x+brush ; xk++)
        {
          xi = mri->xi[xk] ;
          if (val)
          {
            repaint_needed = 1 ;
            switch (mri->type)
            {
            default:
            case MRI_UCHAR:
              old_val = MRIseq_vox(mri, xi, yi, zi, mri_frames[which])  ;
              MRIseq_vox(mri, xi, yi, zi, mri_frames[which]) = val ;
              break ;
            case MRI_FLOAT:
              old_val = (BUFTYPE)MRIFseq_vox(mri, xi,yi,zi,mri_frames[which]);
              MRIFseq_vox(mri, xi, yi, zi, mri_frames[which]) = (float)val ;
              break ;
            }
          }
        }
      }
    }
  }
  else switch (event_id(event))
    {
    case LOC_DRAG:
      if (!event_left_is_down(event))
        break ;
    case MS_LEFT:
      switch (mri->type)
      {
      case MRI_UCHAR:
        if (talairach)
          XVprintf(xvf, 0, "T: (%d,%d,%d) --> %d",
                   nint(xt),nint(yt),nint(zt),
                   MRIseq_vox(mri, nint(xtv), nint(ytv), nint(ztv),frame));
        else
          XVprintf(xvf, 0, "(%d,%d,%d) --> %d",x,y,z,
                   MRIseq_vox(mri,x,y,z,frame));
        break ;
      case MRI_FLOAT:
        if (talairach)
          XVprintf(xvf, 0, "T: (%d,%d,%d) --> %2.3f",
                   nint(xt),nint(yt),nint(zt),
                   MRIFseq_vox(mri, nint(xtv), nint(ytv), nint(ztv),frame));
        else
          XVprintf(xvf, 0, "(%d,%d,%d) --> %2.3f",x,y,z,
                   MRIFseq_vox(mri, x, y, z,frame));
        break ;
      }
      break ;
    default:
      if (event_is_up(event)) switch ((char)event->ie_code)
        {
        case 'B':
          brush++ ;
          break ;
        case 'b':
          if (brush > 0)
            brush-- ;
          break ;
        case 'S':   /* change all views and slices to be the same */
          XVMRIsetView(xvf, which, mri_views[which]) ;
          break ;
#if 0
        case '0':
          fprintf(stderr, "turning off voxel (%d, %d, %d)\n", x, y, z) ;
          switch (mri->type)
          {
          case MRI_UCHAR:
            MRIseq_vox(mri, x, y, z, mri_frames[which]) = 1 ;
            break ;
          case MRI_FLOAT:
            MRIseq_vox(mri, x, y, z, mri_frames[which]) = 0.0f ;
            break ;
          }
          XVMRIshowFrame(xvf, mri, which, mri_slices[which], mri_frames[which]) ;
          break ;
        case '1':
          fprintf(stderr, "turning on voxel (%d, %d, %d)\n", x, y, z) ;
          switch (mri->type)
          {
          case MRI_UCHAR:
            MRIseq_vox(mri, x, y, z, mri_frames[which]) = 255 ;
            break ;
          case MRI_FLOAT:
            MRIseq_vox(mri, x, y, z, mri_frames[which]) = 1.0f ;
            break ;
          }
          XVMRIshowFrame(xvf, mri, which, mri_slices[which], mri_frames[which]) ;
          break ;
#endif
        case 'G':
        case 'g':
          /* look in 4 places for edit.dat - same dir as image, tmp/edit.dat
             ../tmp and ../../tmp
             */
          sprintf(fname, "%s/edit.dat", image_path) ;
          if (!FileExists(fname))
          {
            sprintf(fname, "%s/../tmp/edit.dat", image_path) ;
            if (!FileExists(fname))
            {
              sprintf(fname, "%s/../../tmp/edit.dat", image_path) ;
              if (!FileExists(fname))
              {
                sprintf(fname, "%s/tmp/edit.dat", image_path) ;
                if (!FileExists(fname))
                {
                  XVprintf(xvf,0,"could not find edit.dat from %s", image_path);
                  return(ERROR_NO_FILE) ;
                }
              }
            }
          }
          fp = fopen(fname, "r") ;
          if (fscanf(fp, "%f %f %f", &xf, &yf, &zf) != 3)
          {
            XVprintf(xvf, 0, "could not scan coordinates out of %s", fname) ;
            fclose(fp) ;
            return(ERROR_BAD_FILE) ;
          }
          if (fscanf(fp, "%f %f %f", &xft, &yft, &zft) != 3)
          {
            XVprintf(xvf,0,"could not scan Talairach coordinates out of %s",
                     fname);
            fclose(fp) ;
            return(ERROR_BAD_FILE) ;
          }
          fclose(fp) ;
          if (talairach)
            MRItalairachToVoxel(mri, (Real)xft,(Real)yft,(Real)zft,&xv,&yv,&zv);
          else
            MRIworldToVoxel(mri, (Real)xf, (Real)yf, (Real)zf, &xv, &yv, &zv) ;
          XVMRIsetPoint(xvf, which, nint(xv), nint(yv), nint(zv)) ;
          XVprintf(xvf, 0, "current point: (%d, %d, %d) --> (%d, %d, %d)",
                   nint(xf), nint(yf), nint(zf), nint(xv), nint(yv), nint(zv)) ;
#if 0
          fprintf(stderr, "read (%2.3f, %2.3f, %2.3f) and (%2.3f, %2.3f, %2.3f)\n",
                  xf, yf, zf, xft, yft, zft) ;
          fprintf(stderr, "voxel (%d, %d, %d)\n", nint(xv), nint(yv), nint(zv)) ;
#endif
          break ;
        case 'W':
        case 'w':
          /* look in 4 places for edit.dat - same dir as image, tmp/edit.dat
             ../tmp and ../../tmp
             */
          sprintf(fname, "%s/edit.dat", image_path) ;
          if (!FileExists(fname))
          {
            sprintf(fname, "%s/../tmp/edit.dat", image_path) ;
            if (!FileExists(fname))
            {
              sprintf(fname, "%s/../../tmp/edit.dat", image_path) ;
              if (!FileExists(fname))
              {
                sprintf(fname, "%s/tmp/edit.dat", image_path) ;
                if (!FileExists(fname))
                {
                  XVprintf(xvf,0,"could not find edit.dat from %s", image_path) ;
                  return(ERROR_BAD_FILE) ;
                }
              }
            }
          }
          fp = fopen(fname, "w") ;
          MRIvoxelToWorld(mri, (Real)x_click, (Real)y_click, (Real)z_click,
                          &xr, &yr, &zr) ;
          MRIvoxelToTalairach(mri, (Real)x_click, (Real)y_click, (Real)z_click,
                              &xt, &yt, &zt) ;
          fprintf(fp, "%f %f %f\n", (float)xr, (float)yr, (float)zr) ;
          fprintf(fp, "%f %f %f\n", (float)xt, (float)yt, (float)zt) ;
          fclose(fp) ;
#if 0
          fprintf(stderr, "wrote (%2.3f, %2.3f, %2.3f) and (%2.3f, %2.3f, %2.3f)\n",
                  xr, yr, zr, xt, yt, zt) ;
          fprintf(stderr, "voxel (%d, %d, %d)\n", x_click, y_click, z_click) ;
#endif
          break ;
        case 'x':
        case 'X':
          XVMRIsetView(xvf, which, MRI_SAGITTAL) ;
          break ;
        case 'y':
        case 'Y':
          XVMRIsetView(xvf, which, MRI_HORIZONTAL) ;
          break ;
        case 'z':
        case 'Z':
          XVMRIsetView(xvf, which, MRI_CORONAL) ;
          break ;
        case 'T':
          talairach = 1 ;
          break ;
        case 't':
          talairach = 0 ;
          break ;
#if 0
        case 'h':
        {
          HISTOGRAM *histo ;
          float     fmin, fmax ;
          XV_FRAME  *xvnew ;

          MRIvalRange(mri, &fmin, &fmax) ;
          histo = MRIhistogram(mri, (int)(fmax - fmin + 1)) ;
          xvnew = XValloc(1, 1, 2, 200, 200, "histogram tool", NULL) ;
          XVshowHistogram(xvnew, 0, histo) ;
          xv_main_loop(xvnew->frame);
          HISTOfree(&histo) ;
        }
#else
        case 'h':
#endif
          break ;
        }
      break ;
    }

#if 0
  if (event_is_up(event) && (event_id(event) == MS_LEFT))
#else
  if (event_id(event) == MS_LEFT)
#endif
  {
    int view, old_which ;

    if (which_click != which)   /* erase old point */
    {
      if (!event_is_up(event))
        dont_redraw = 1 ;  /* inhibit 1 redraw of surface */
      old_which = which_click ;
      which_click = dimage->which ;
      XVrepaintImage(xvf, old_which) ;
      view = mri_views[which] ;
      sprintf(view_str, "view: %s",
              view == MRI_CORONAL ? "CORONAL" :
              view == MRI_SAGITTAL ? "SAGITTAL" : "HORIZONTAL") ;
      xv_set(view_panel, PANEL_LABEL_STRING, view_str, NULL) ;
    }

    if (!event_is_up(event))
      dont_redraw = 1 ;  /* inhibit 1 redraw of surface */
    XVMRIsetPoint(xvf, which_click, x, y, z) ;
  }
  return(NO_ERROR) ;
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
  case MRI_SAGITTAL:
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
int
XVMRIdraw3DPoint(XV_FRAME *xvf, int which, int x,int y,int z,int color)
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
    if (z != mri_slices[which])
      return(0) ;
    break ;
  case MRI_HORIZONTAL:
    xi = x ;
    yi = z ;
    if (y != mri_slices[which])
      return(0) ;
    break ;
  case MRI_SAGITTAL:
    xi = z ;
    yi = y ;
    if (x != mri_slices[which])
      return(0) ;
    break ;
  }
  XVdrawPoint(xvf, which, xi, yi, color) ;
  return(1) ;
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
  case MRI_SAGITTAL:
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
  float  mag, fmin, fmax ;

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
    case MRI_SAGITTAL:
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
    if (slice >= mri->depth)
      slice = mri->depth-1 ;
    if (slice < 0)
      slice = 0 ;
    break ;
  case MRI_SAGITTAL:
    if (slice >= mri->width)
      slice = mri->width-1 ;
    break ;
  case MRI_HORIZONTAL:
    if (slice >= mri->height)
      slice = mri->height-1 ;
    break ;
  }


  mri_slices[which] = slice ;
  I = MRItoImageView(mri, Idisplay[which], slice, mri_views[which], frame) ;
  if (!I)
    return(NULL) ;


  mag = MIN((float)xvf->orig_disp_rows / (float)I->rows,
            (float)xvf->orig_disp_cols / (float)I->cols) ;

  XVsetImageSize(xvf, which, nint((float)I->rows*mag),
                 nint((float)I->cols*mag));
  /*  XVresize(xvf) ;*/

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
  MRIvalRange(mri, &fmin, &fmax) ;
  XVshowImageRange(xvf, which, I, 0, fmin, fmax) ;


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
XVMRIshowRange(XV_FRAME *xvf, MRI *mri, int which, int slice,
               float fmin, float fmax)
{
  IMAGE  *I ;
  float  mag ;

  if (which_click == which && mri != mris[which])
    which_click = -1 ;  /* a new MR image shown */

  if (!mri)
    return(NULL) ;

  mri_frames[which] = 0 ;

  if (slice < 0)  /* set slice to middle of slice direction */
  {
    switch (mri_views[which])
    {
    case MRI_CORONAL:
      slice = (mri->imnr0 + mri->imnr1) / 2 ;
      break ;
    case MRI_SAGITTAL:
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
    if (slice >= mri->depth)
      slice = mri->depth-1 ;
    if (slice < 0)
      slice = 0 ;
    break ;
  case MRI_SAGITTAL:
    if (slice >= mri->width)
      slice = mri->width-1 ;
    break ;
  case MRI_HORIZONTAL:
    if (slice >= mri->height)
      slice = mri->height-1 ;
    break ;
  }


  mri_slices[which] = slice ;
  I = MRItoImageView(mri, Idisplay[which], slice, mri_views[which], 0) ;
  if (!I)
    return(NULL) ;


  mag = MIN((float)xvf->orig_disp_rows / (float)I->rows,
            (float)xvf->orig_disp_cols / (float)I->cols) ;

  XVsetImageSize(xvf, which, nint((float)I->rows*mag),
                 nint((float)I->cols*mag));
  /*  XVresize(xvf) ;*/

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
  XVshowImageRange(xvf, which, I, 0, fmin, fmax) ;


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
  XVsetWriteFunc(xvf, "WRITE MRI", "MR image file name", mri_write_func) ;
  sprintf(view_str, "view: CORONAL") ;
  view_menu = (Menu)
              xv_create((Xv_opaque)NULL, MENU,
                        MENU_NOTIFY_PROC,    viewMenuItem,
                        MENU_STRINGS,        "CORONAL", "SAGITTAL", "HORIZONTAL", NULL,
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
  MRI    *mri ;
  int    which, view ;

  which = which_click ;
  if (which_click < 0)   /* no current window */
    which = 0 ;

  mri = mris[which] ;
  if (!mri)              /* no mri in current window */
    return ;

  menu_str = (char *)xv_get(menu_item, MENU_STRING) ;

  if (!stricmp(menu_str, "CORONAL"))
    view = MRI_CORONAL ;
  else if (!stricmp(menu_str, "SAGITTAL"))
    view = MRI_SAGITTAL ;
  else
    view = MRI_HORIZONTAL ;

  XVMRIsetView(xvf, which, view) ;
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

  if (fixed [which])
    return(Iold) ;

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
      mri_slices[which] += dir ;   /* new depth is valid */
      if (Idisplay[which] && (Idisplay[which] != I))
        ImageFree(&Idisplay[which]) ;
      Idisplay[which] = I ;
    }

    /* if current image is changing depth, then change the location of
       the current click to agree with the current depth.
       */
    if (which_click == which)
    {
      switch (mri_views[which])
      {
      case MRI_CORONAL:
        z_click = depth - mri->imnr0 ;
        break ;
      case MRI_SAGITTAL:
        x_click = depth ;
        break ;
      case MRI_HORIZONTAL:
        y_click = depth ;
        break ;
      }
      dont_redraw = 2 ;  /* inhibit 1 redraw of surface */
      XVMRIsetPoint(xvf, which_click, x_click, y_click, z_click) ;
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
  int which = dimage->which ;
  MRI *mri ;

  mri = mris[which] ;
  if (!mri)
    return ;

  if (mri_surface && !dont_redraw)
  {
    VERTEX  *v, *vn ;
    int     vno, n ;
    Real    xv, yv, zv, slice, xv2, yv2, zv2, dx, dy ;

    for (vno = 0 ; vno < mri_surface->nvertices ; vno++)
    {
      v = &mri_surface->vertices[vno] ;
      if (v->ripflag)
        continue ;
      MRIworldToVoxel(mri, v->x, v->y, v->z, &xv, &yv, &zv) ;

      switch (mri_views[which])
      {
      default:
      case MRI_CORONAL:
        /*
          Z=0 is the back of the head,
          X=0 is the right side of the head
          Y=0 is the neck/brain stem
          */
        slice = mri_depths[which] - mri->imnr0 ;
        if ((zv > slice-1) && (zv < slice+1))
        {
          for (n = 0 ; n < v->vnum ; n++)
          {
            vn = &mri_surface->vertices[v->v[n]] ;
            MRIworldToVoxel(mri,vn->x, vn->y, vn->z, &xv2, &yv2, &zv2) ;
            if ((zv2 > slice-.5) && (zv2 < slice+.5))
            {
              dx = xv2-xv ;
              dy = yv2-yv ;
              XVdrawLinef(xvf, which, xv, yv, dx, dy, XCYAN) ;
            }
          }
        }
        break ;
      case MRI_SAGITTAL:
        /*
          Z=0 is the back of the head,
          X=0 is the right side of the head
          Y=0 is the neck/brain stem
          */
        slice = mri_slices[which] ;
        if ((xv > slice-.5) && (xv < slice+.5))
        {
          for (n = 0 ; n < v->vnum ; n++)
          {
            vn = &mri_surface->vertices[v->v[n]] ;
            MRIworldToVoxel(mri, vn->x, vn->y, vn->z, &xv2, &yv2, &zv2);
            if ((xv2 > slice-.5) && (xv2 < slice+.5))
            {
              dx = zv2-zv ;
              dy = yv2-yv ;
              XVdrawLinef(xvf, which, zv, yv, dx, dy, XCYAN) ;
            }
          }
        }
        break ;
      case MRI_HORIZONTAL:
        /*
          Z=0 is the back of the head,
          X=0 is the right side of the head
          Y=0 is the neck/brain stem
          */
        slice = mri_depths[which] ;
        if ((yv > slice-.5) && (yv < slice+.5))
        {
          for (n = 0 ; n < v->vnum ; n++)
          {
            vn = &mri_surface->vertices[v->v[n]] ;
            MRIworldToVoxel(mri, vn->x, vn->y, vn->z, &xv2, &yv2, &zv2) ;
            if ((yv2 > slice-.5) && (yv2 < slice+.5))
            {
              dx = xv2-xv ;
              dy = zv2-zv ;
              XVdrawLinef(xvf, which, xv, zv, dx, dy, XCYAN) ;
            }
          }
        }
        break ;
      }
    }
  }
  dont_redraw = 0 ;
  if (dimage->sync || (which == which_click))
    XVMRIdrawPoint(xvf, which, mri_views[which], 0,
                   mri, x_click, y_click, z_click, XRED) ;
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
  int     slice, which2, offset, slice2, sync ;
  DIMAGE  *dimage, *dimage2 ;
  MRI     *mri, *mri2 ;
  char    *menu_str ;
  float   xsize, ysize, zsize ;

  if (!mris[which])
    return(NO_ERROR) ;

  mri = mris[which] ;

  switch (view)
  {
  default:
  case MRI_CORONAL:
    slice = z_click + mri->imnr0 ;
    menu_str = "CORONAL" ;
    break ;
  case MRI_SAGITTAL:
    slice = x_click ;
    menu_str = "SAGITTAL" ;
    break ;
  case MRI_HORIZONTAL:
    slice = y_click ;
    menu_str = "HORIZONTAL" ;
    break ;
  }

  dimage = XVgetDimage(xvf, which, DIMAGE_IMAGE) ;

  /* check to see if showing a zoomed image, if so, change region */
  if (/*(view != mri_views[which]) && */(dimage->dx > 0))
  {
    switch (view)
    {
    default:
    case MRI_CORONAL:
      slice = z_click + mri->imnr0 ;
      dimage->x0 = x_click - dimage->dx/2 ;
      dimage->y0 = mri->height - (y_click + dimage->dy/2) ;
      menu_str = "CORONAL" ;
      break ;
    case MRI_SAGITTAL:
      slice = x_click ;
      menu_str = "SAGITTAL" ;
      dimage->x0 = z_click - dimage->dx/2 ;
      dimage->y0 = mri->height - (y_click + dimage->dy/2) ;
      break ;
    case MRI_HORIZONTAL:
      dimage->x0 = x_click - dimage->dx/2 ;
      dimage->y0 = mri->depth - (z_click + dimage->dy/2) ;
      slice = y_click ;
      menu_str = "HORIZONTAL" ;
      break ;
    }
  }

  if (which == which_click)
  {
    sprintf(view_str, "view: %s", menu_str) ;
    xv_set(view_panel, PANEL_LABEL_STRING, view_str, NULL) ;
  }

  sync = dimage->sync ;
  if (sync && !show_three_views)
  {
    for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
    {
      mri2 = mris[which2] ;
      if (which2 == which || fixed [which2] || !mri2)
          continue ;
      dimage2 = XVgetDimage(xvf, which2, DIMAGE_IMAGE) ;
      xsize = mri2->xsize / mri->xsize ;
      ysize = mri2->ysize / mri->ysize ;
      zsize = mri2->zsize / mri->zsize ;
      if (dimage2 && mri2)
      {
        switch (view)
        {
        case MRI_CORONAL:
          offset = mri->zstart - mri2->zstart ;
          slice2 = slice - mri->imnr0 ;  /* turn it into an index */
          slice2 = (slice2+offset) * zsize + mri2->imnr0 ;
          break ;
        case MRI_SAGITTAL:
          offset = mri->xstart - mri2->xstart ;
          slice2 = (slice+offset) * xsize ;
          break ;
        case MRI_HORIZONTAL:
          offset = mri->ystart - mri2->ystart ;
          slice2 = (slice+offset) * ysize ;
          break ;
        default:
          slice2 = slice ;
        }
        dimage2->dx = dimage->dx ;
        dimage2->dy = dimage->dy ;
        dimage2->x0 = dimage->x0 ;
        dimage2->y0 = dimage->y0 ;
        if (view == mri_views[which2])
          XVMRIredisplayFrame(xvf, mri2, which2, slice2, mri_frames[which2]) ;
        else
        {
          mri_views[which2] = view ;
          XVMRIshowFrame(xvf, mri2, which2, slice2, mri_frames[which2]) ;
        }
      }
    }
  }

  /* will reset syncs */
  if (view != mri_views[which])
  {
    mri_views[which] = view ;
    XVMRIshowFrame(xvf, mri, which, slice, mri_frames[which]) ;
    if (show_three_views)
      XVsyncAll(xvf, 0) ;
#if 0
    if (sync)  /* if they were synced, reinstate it */
      XVsyncAll(xvf, which) ;
#endif
  }
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
mri_write_func(Event *event, DIMAGE *dimage, char *fname)
{
  int  which ;
  MRI  *mri ;

  which = dimage->which ;
  mri = mris[which] ;
  if (!mri)
    return(ERROR_BADPARM) ;
  MRIwrite(mri, fname) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVMRIsetPoint(XV_FRAME *xvf, int which, int x, int y, int z)
{
  DIMAGE   *dimage, *dimage2 ;
  int      which2, x2, y2, z2, slice, old_redraw = dont_redraw ;
  MRI      *mri, *mri2 ;
  char     fmt[150], title[50], buf[100] ;
  Real     xr, yr, zr ;

  if (which != which_click)   /* erase old point */
  {
    int old_which = which_click ;

    which_click = which ;  /* must do before repaint so point isn't drawn */
#if 0
    xvmriRepaintValue(xvf, old_which, x, y, z) ;
#else
    XVrepaintImage(xvf, old_which) ;
#endif
  }
  else which_click = which ;
  dimage = XVgetDimage(xvf, which, DIMAGE_IMAGE) ;
  mri = mris[which] ;
  x_click = x ;
  y_click = y ;
  z_click = z ;

  switch (mri_views[which])   /* make sure correct slice is shown */
  {
  default:
  case MRI_CORONAL:
    slice = z + mri->imnr0 ;
    break ;
  case MRI_SAGITTAL:
    slice = x ;
    break ;
  case MRI_HORIZONTAL:
    slice = y ;
    break ;
  }

  if (slice != mri_depths[which])
    XVMRIshowFrame(xvf, mri, which, slice, mri_frames[which]) ;
  else
    XVrepaintImage(xvf, which) ;

  if (dimage->sync)
  {
    for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
    {
      dimage2 = XVgetDimage(xvf, which2, DIMAGE_IMAGE) ;
      mri2 = mris[which2] ;
      if (dimage2 /* && (which2 != which) */ && mri2)
      {
        MRIvoxelToVoxel(mri, mri2, (Real)x, (Real)y, (Real)z, &xr, &yr, &zr);
        x2 = nint(xr) ;
        y2 = nint(yr) ;
        z2 = nint(zr) ;
        if (x2 < 0 || x2 >= mri2->width || y2 < 0 || y2 >= mri2->height ||
            z2 < 0 || z2 >= mri2->depth)
          continue ;
#if 0
        x2 = MAX(0, x2) ;
        y2 = MAX(0, y2) ;
        z2 = MAX(0, z2) ;
        x2 = MIN(mri2->width-1, x2) ;
        y2 = MIN(mri2->height-1, y2) ;
        z2 = MIN(mri2->depth-1, z2) ;
#endif
        XVgetTitle(xvf, which2, title, 0) ;
        sprintf(fmt, "%%10.10s: (%%3d, %%3d) --> %%2.%dlf\n",xvf->precision);
        switch (mri2->type)
        {
        case MRI_UCHAR:
          sprintf(buf, "%d", MRIseq_vox(mri2, x2, y2,z2,mri_frames[which2]));
          break ;
        case MRI_FLOAT:
          sprintf(buf,"%2.3f",MRIFseq_vox(mri2,x2,y2,z2,mri_frames[which2]));
          break ;
        }
        XVshowImageTitle(xvf, which2, "%s (%s)", title, buf) ;
        if (which2 == which)
          continue ;
        /*
          this dont redraw stuff is such a hack, but I can't think
          of an easy way to prevent tons of irritating redraws
          otherwise.
          */
        if (old_redraw < 2 || which2 > which)  /* will get redrawn itself */
          dont_redraw = old_redraw ;
        XVrepaintImage(xvf, which2) ;
      }
    }
  }
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              set the current path. This will be used to try and
              pass data back and forth between applications like
              surfer.
----------------------------------------------------------------------*/
int
XVMRIsetImageName(XV_FRAME *xvf, char *image_name)
{
  FileNamePath(image_name, image_path) ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Redisplay just the point (x,y,z)
----------------------------------------------------------------------*/
#if 0
static int
xvmriRepaintValue(XV_FRAME *xvf, int which, int x, int y, int z)
{
  MRI    *mri ;
  IMAGE  *I ;
  int    xp, yp ;
  DIMAGE *dimage ;
  float  fmin, fmax ;

  dimage = XVgetDimage(xvf, which, DIMAGE_IMAGE) ;
  I = Idisplay[which] ;
  mri = mris[which] ;

  switch (mri_views[which])
  {
  default:
  case MRI_CORONAL:
    xp = x ;
    yp = y ;
    break ;
  case MRI_SAGITTAL:
    xp = z - mri->imnr0 ;
    yp = y ;
    break ;
  case MRI_HORIZONTAL:
    xp = x ;
    yp = z-mri->imnr0 ;
    break ;
  }
  yp = I->rows - (yp+1) ;
  switch (I->pixel_format)
  {
  case PFBYTE:
    *IMAGEpix(I, xp, yp) = MRIseq_vox(mri,x, y,z, mri_frames[which]);
    break ;
  case PFFLOAT:
    *IMAGEFpix(I, xp, yp) = MRIFseq_vox(mri,x, y,z, mri_frames[which]);
    break ;
  }
  MRIvalRange(mri, &fmin, &fmax) ;
  XVshowImageRange(xvf, which, I, 0, fmin, fmax) ;
  return(NO_ERROR) ;
}
#endif
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
IMAGE *
XVMRIredisplayFrame(XV_FRAME *xvf, MRI *mri, int which, int slice,int frame)
{
  IMAGE  *I ;
  float  fmin, fmax ;

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
    case MRI_SAGITTAL:
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
    if (slice >= mri->depth)
      slice = mri->depth-1 ;
    if (slice < 0)
      slice = 0 ;
    break ;
  case MRI_SAGITTAL:
    if (slice >= mri->width)
      slice = mri->width-1 ;
    break ;
  case MRI_HORIZONTAL:
    if (slice >= mri->height)
      slice = mri->height-1 ;
    break ;
  }

  mri_slices[which] = slice ;
  I = MRItoImageView(mri, Idisplay[which], slice, mri_views[which], frame) ;
  if (!I)
    return(NULL) ;


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

  MRIvalRange(mri, &fmin, &fmax) ;
  XVshowImageRange(xvf, which, I, 0, fmin, fmax) ;


  if (Idisplay[which] && (I != Idisplay[which]))
    ImageFree(&Idisplay[which]) ;

  Idisplay[which] = I ;

  mris[which] = mri ;
  return(I) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Allow user to specify some MRs that shouldn't move in depth
             with the rest.
----------------------------------------------------------------------*/
int
XVMRIfixDepth(XV_FRAME *xvf, int which, int fix)
{
  fixed [which] = fix ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
             Allow user to specify some MRs that shouldn't move in depth
             with the rest.
----------------------------------------------------------------------*/
int
XVMRIspecifySurface(XV_FRAME *xvf, MRI_SURFACE *mrisurf)
{
  mri_surface = mrisurf ;
  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              Show three views of the same MRI at the same time,
              similar to syncing, but don't sync the views.
----------------------------------------------------------------------*/
int
XVMRIshowThreeViews(XV_FRAME *xvf)
{
  XVsyncAll(xvf, 0) ;
  show_three_views = 1 ;

  return(NO_ERROR) ;
}

