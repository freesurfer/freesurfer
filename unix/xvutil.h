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


typedef struct
{
  Canvas      canvas ;
  Panel_item  title_item ;
  char        title_string[STR_LEN] ;
  IMAGE       *dispImage ;
  IMAGE       *sourceImage ;
  XImage      *ximage ;
  Window      window ;
  GC          clearGC, greenGC, blueGC, redGC, 
              xorGC, whiteGC, cyanGC, yellowGC, purpleGC ;
  float       xscale ;
  float       yscale ;
  int         rescale_range ;
  int         frame ;
  int         used ;
  int         entered ;     /* to prevent re-entrancy in repaint proc. */
  int         which ;
  int         x, y ;        /* position on canvas */
  int         bshift ;      /* for modifying brightness of displayed image */
} DISPLAY_IMAGE, DIMAGE ;

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
} XV_FRAME ;


XV_FRAME *XValloc(int rows, int cols, int button_rows, int display_rows, 
                  int display_cols, char *name, Notify_value (*poll)(void)) ;
int XVprintf(XV_FRAME *xvf, int which, ...) ;
void XVclearImage(XV_FRAME *xvf, int which, int dotitle) ;
void XVclearImageTitle(XV_FRAME *xvf, int which) ;
int  XVshowImageTitle(XV_FRAME *xvf, int which, ...) ;
void XVshowImage(XV_FRAME *xvf, int which, IMAGE *image, int frame) ;
void XVdrawBox(XV_FRAME *xvf, int which, int x, int y, int dx, int dy, 
               int color) ;
void XVdrawLine(XV_FRAME *xvf, int which, int x, int y, int dx, int dy,  
               int color) ;
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
void buttonQuit(Panel_item item, Event *event) ;
int XVaddImageCol(XV_FRAME *xvf) ;
int XVdeleteImageCol(XV_FRAME *xvf) ;
int XVsetPrecision(XV_FRAME *xvf, int precision) ;
int XVbrighten(XV_FRAME *xvf, int which, int offset) ;

#define WINDOW_PAD          3


#ifdef IRIX
#define ROW_HEIGHT           40
#define CHAR_WIDTH           16
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


#endif
