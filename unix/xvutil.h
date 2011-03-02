/**
 * @file  xvutil.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.38 $
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


#ifndef XVUTIL_H
#define XVUTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xproto.h>
#include <X11/Xutil.h>

#include <xview/xv_c_types.h>
#include <xview/xview.h>
#include <xview/frame.h>
#include <xview/panel.h>
#include <xview/canvas.h>
#include <xview/scrollbar.h>
#include <xview/cms.h>
#include <xview/xv_xrect.h>
#include <xview/notify.h>

#include "image.h"
#include "const.h"
#include "histo.h"

#if 1
#ifdef STDC_SWITCHED
#undef STDC_SWITCHED
#define __STDC__
#endif
#endif

#define MSGS      3

/* colors for XVdraw... calls */
#define XXOR      -1
#define XRED      0
#define XWHITE    1
#define XGREEN    2
#define XBLUE     3
#define XCYAN     4
#define XYELLOW   5
#define XPURPLE   6

#define MAX_COLORS       64

typedef struct
{
  Canvas      canvas ;
  Panel_item  title_item ;
  char        title_string[STR_LEN] ;
  IMAGE       *dispImage ;    /* image displayed on screen */
  IMAGE       *sourceImage ;  /* copy of source image passed to XVshowImage */
  IMAGE       *oSourceImage ; /* source image passed to XVshowImage */
  XImage      *ximage ;
  Window      window ;
  GC          clearGC, greenGC, blueGC, redGC,
  xorGC, whiteGC, cyanGC, yellowGC, purpleGC ;
  float       xscale ;        /* scale of dispImage relative to sourceImage */
  float       yscale ;        /* scale of dispImage relative to sourceImage */
  int         rescale_range ;
  int         frame ;       /* frame of sourceImage which is shown*/
  int         used ;        /* whether this dimage is in use */
  int         entered ;     /* to prevent re-entrancy in repaint proc. */
  int         which ;       /* index into table in XV_FRAME structure */
  int         x, y ;        /* position on canvas */
  int         bshift ;      /* for modifying brightness of displayed image */
  float       gamma[MAX_COLORS] ;  /* for gamma correction */
  float       zoom ;        /* current zoom scale */
  int         x0 ;          /* x zoom origin */
  int         y0 ;          /* y zoom origin */
  int         dx ;          /* width of zoomed image */
  int         dy ;          /* height of zoomed image */

  /* the next 4 are temporary copies of the prior 4 while drawing is underway */
  int         x1 ;
  int         y1 ;          /* y box origin */
  int         dx1 ;         /* width of zoomed image */
  int         dy1 ;         /* height of zoomed image */

  IMAGE       *zoomImage ;  /* zoomed image */
  int         sync ;        /* sync zooming and stuff with other images */
  HISTOGRAM   *histo ;
  float       fmin ;        /* fixed minimum for display scaling */
  float       fmax ;        /* fixed maximum for display scaling */
}
DISPLAY_IMAGE, DIMAGE ;

typedef struct
{
  Frame           frame ;
  Panel           panel ;
  Display         *display ;
  int             screen ;
  unsigned long   white_pixel ;
  unsigned long   black_pixel ;
  unsigned long   red_pixel ;
  unsigned long   blue_pixel ;
  unsigned long   green_pixel ;
  unsigned long   cyan_pixel ;
  unsigned long   yellow_pixel ;
  unsigned long   purple_pixel ;
  Cms             cms ;
  GC              gc ;
  int             rows ;
  int             cols ;
  int             button_rows ;
  int             display_rows ;
  int             display_cols ;
  char            msg_str[MSGS][STR_LEN] ;
  Panel_item      msg_item[MSGS] ;
  DIMAGE          **dimages ;
  int             min_panel_width ;
  int             precision ;
  int             ydir ;        /* direction of y coord relative to xview */
  int             noprint ;     /* 1 --> don't do XVprintf in event_handler */
  int         rescale ;     /* rescale range with all images of this format */
  char            fname_prompt[30] ;
  char            file_name[100] ;
  Frame           fname_frame ;
  Panel           fname_panel ;
  Panel_item      fname_panel_item ;
  int             (*fname_func)(char *fname) ;
  IMAGE           *(*get_next_image)(IMAGE *Iold, int which, int dir) ;
  int             orig_disp_rows ;  /* original value of display_rows */
  int             orig_disp_cols ;  /* original values of display_cols */
  int             (*write_func)(Event *event, DIMAGE *dimage, char *str) ;
}
XV_FRAME ;

int  XVsetWriteFunc(XV_FRAME *xvf, char *frame_name, char *prompt_str,
                    int (*write_func)(Event *event,DIMAGE *dimage,char *str));
int  XVshowHistogram(XV_FRAME *xvf, int which, HISTOGRAM *mrih) ;
XV_FRAME *XValloc(int rows, int cols, int button_rows, int display_rows,
                  int display_cols, char *name, Notify_value (*poll)(void)) ;
int XVprintf(XV_FRAME *xvf, int which, ...) ;
void XVclearImage(XV_FRAME *xvf, int which, int dotitle) ;
void XVclearImageTitle(XV_FRAME *xvf, int which) ;
int  XVshowImageTitle(XV_FRAME *xvf, int which, ...) ;
void XVshowImage(XV_FRAME *xvf, int which, IMAGE *image, int frame) ;
void XVshowImageRange(XV_FRAME *xvf, int which, IMAGE *image, int frame,
                      float fmin, float fmax) ;
void XVdrawBox(XV_FRAME *xvf, int which, int x, int y, int dx, int dy,
               int color) ;
void XVdrawLine(XV_FRAME *xvf, int which, int x, int y, int dx, int dy,
                int color) ;
void XVdrawLinef(XV_FRAME *xvf, int which, float x, float y, float dx,
                 float dy,  int color) ;
void XVdrawArrow(XV_FRAME *xvf, int which, int x, int y, float dx, float dy,
                 int color) ;
void XVdrawPoint(XV_FRAME *xvf, int which, int x, int y, int color) ;
void XVsetParms(void (*event_handler)(Event *event, DIMAGE *dimage)) ;
void XVsetKBhandler(void (*kb_handler)(Event *event, DIMAGE *dimage)) ;
void XVsetRepaintHandler(void(*repaint_handler)(XV_FRAME *xvf,DIMAGE *dimage));
void XVshowVectorImage(XV_FRAME *xvf, int which, int x0, int y0,
                       int width, int height, int color, IMAGE *image) ;
void XVsetQuitFunc(void (*quit_func)(void)) ;
void XVrepaintImage(XV_FRAME *xvf, int which) ;
int XVsetMessagePosition(XV_FRAME *xvf, int which, int col, int row) ;
int XVsetImageSize(XV_FRAME *xvf, int which, int rows, int cols) ;
void XVsetMinPanelWidth(XV_FRAME *xvf, int min_panel_width) ;
int XVresize(XV_FRAME *xvf) ;
int XVchangeDisplaySize(XV_FRAME *xvf) ;
void buttonQuit(Panel_item item, Event *event) ;
int XVaddImageCol(XV_FRAME *xvf) ;
int XVdeleteImageCol(XV_FRAME *xvf) ;
int XVsetPrecision(XV_FRAME *xvf, int precision) ;
int XVbrighten(XV_FRAME *xvf, int which, int offset) ;
int XVgamma(XV_FRAME *xvf, int which, float beta) ;
int XVzoom(XV_FRAME *xvf, int which, float zoom) ;
int XVsync(XV_FRAME *xvf, int which, int sync) ;
int XVdoSync(XV_FRAME *xvf, int which) ;
int XVsyncAll(XV_FRAME *xvf, int which) ;
int XVunsyncAll(XV_FRAME *xvf, int which) ;
int XVshowAllSyncedImages(XV_FRAME *xvf, int which) ;
int XVsetPrintStatus(XV_FRAME *xvf, int status) ;
int XVsetYDir(XV_FRAME *xvf, int ydir) ;
int XVshowAll(XV_FRAME *xvf) ;
int XVgetFileName(XV_FRAME *xvf, char *default_fname,
                  int (*fname_func)(char *fname), ...) ;
int XVsetDepthFunc(XV_FRAME *xvf,
                   IMAGE *(*get_image)(IMAGE *Iold, int which, int dir)) ;
DIMAGE *XVgetDimage(XV_FRAME *xvf, int which, int type) ;
char   *XVgetTitle(XV_FRAME *xvf,int which, char *title, int with_value) ;

#define WINDOW_PAD          3


#ifdef IRIX
#define ROW_HEIGHT           40
#define CHAR_WIDTH           13
#else
#define ROW_HEIGHT           30
#define CHAR_WIDTH           8
#endif

#define FIRST_BUTTON_ROW     ROW_HEIGHT
#define SECOND_BUTTON_ROW    FIRST_BUTTON_ROW+ROW_HEIGHT
#define THIRD_BUTTON_ROW     SECOND_BUTTON_ROW+ROW_HEIGHT
#define FOURTH_BUTTON_ROW    THIRD_BUTTON_ROW+ROW_HEIGHT
#define FIFTH_BUTTON_ROW     FOURTH_BUTTON_ROW+ROW_HEIGHT
#define SIXTH_BUTTON_ROW     FIFTH_BUTTON_ROW+ROW_HEIGHT
#define SEVENTH_BUTTON_ROW   SIXTH_BUTTON_ROW+ROW_HEIGHT
#define EIGHTH_BUTTON_ROW    SEVENTH_BUTTON_ROW+ROW_HEIGHT
#define LAST_BUTTON_ROW      EIGHTH_BUTTON_ROW


#define DIMAGE_UNUSED        0
#define DIMAGE_IMAGE         1
#define DIMAGE_HISTOGRAM     2
#define DIMAGE_ALLOC         3
#define DIMAGE_ALLOCATED     4

#endif
