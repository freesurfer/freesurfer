/**
 * @file  xwin.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.3 $
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


/*
   @(#)xwin.h  1.2
   8/10/95
*/
/*------------------------------------------------------------------------
      File Name:  xwin.h

         Author:  Bruce Fischl

        Created:  Jan. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef XWIN_H
#define XWIN_H

#include <stdarg.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/keysymdef.h>

#ifndef sparc
#define sparc 1
#endif

typedef struct
{
  Display             *display;
  Window              window;
  int                 xsize,ysize;
  XSizeHints          hint;
  int                 depth,screen;
  GC                  black,white,flip ;
  XFontStruct         *font_struct;
  unsigned long       charwidth, charheight;
  float               scale ;
}
xwindow_type;

#define xBLACK  0
#define xWHITE  1
#define xFLIP   2
#define xRED    3
#define xGREEN  4
#define FULLCIRCLE 23040

#if ANSI
void xInit(xwindow_type *xwindow, int xsize, int ysize) ;
void xEnd(xwindow_type *xwindow) ;
void xNextEvent(xwindow_type *xwindow, XEvent *event) ;
void xDrawLine(xwindow_type *xwindow,int x0, int y0, int x1, int y1,
               int color, int style) ;
void xDrawLines(xwindow_type *xwindow, XSegment *segs, int nsegs, int color) ;

void xDrawCircle(xwindow_type *xwindow, int x0, int y0, int radius,int color);
void xDrawPoint(xwindow_type *xwindow, int x, int y, int color) ;
xwindow_type
*xNewWindow(xwindow_type *parent, int xpos, int ypos, int xsize,
            int ysize, char *name, int border, int event_mask) ;

void xFreeWindow(xwindow_type *xwindow) ;
void xShowWindowAttributes(xwindow_type *xwindow, FILE *fp) ;
void xShowEventMask(unsigned long eventmask, FILE *fp) ;
void xShowEventType(unsigned long event, FILE *fp) ;
int  xUseGrayscale(Display *display, int screen) ;
void xShowEvent(unsigned long event, FILE *fp) ;
void xShowDisplayInfo(xwindow_type *xwin, FILE *fp) ;
int xprintf(xwindow_type *xwindow, int x, int y, char *fmt, ...) ;

#endif

#endif
