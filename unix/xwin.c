/**
 * @file  xwin.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.5 $
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
   @(#)xwin.c 1.2
   8/10/95
*/
/*------------------------------------------------------------------------
      File Name: xwin.c

         Author: Bruce Fischl

        Created: Jan. 1994

    Description:

------------------------------------------------------------------------*/

#include <stdio.h>
#define sparc 1
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "xwin.h"
#include "diag.h"
#include "proto.h"

#define BORDERWIDTH     10
#define FONT            "-adobe-courier-bold-r-normal--*-100-*-*-*-*-iso8859-1"
#define NCOLORS         255

void open_xwindow(xwindow_type *xwindow, char displayname[], int xpos,int ypos,
                  int xsize, int ysize,char title[],long eventmask) ;
void close_xwindow(xwindow_type *xwindow) ;

/* colors for gray scale */
static XColor  xcolors[NCOLORS] ;
static int ncolors = 0 ;   /* actual # of colors allocated */

void
xEnd(xwindow_type *xwindow)
{
  if (!xwindow) return ;

  XFlush(xwindow->display) ;
  close_xwindow(xwindow) ;
}

void
xInit(xwindow_type *xwin, int xsize, int ysize)
{
  long                  eventmask ;

  if (!xwin) return ;

  eventmask = FocusChangeMask | EnterWindowMask | LeaveWindowMask |
              ButtonReleaseMask | ButtonPressMask | KeyReleaseMask | KeyPressMask |
              ExposureMask | ButtonMotionMask;
  if (!xwin->scale) xwin->scale = 1 ;
  open_xwindow(xwin,"",-1,-1, xsize, ysize,
               "DEBUGGING IMAGE", eventmask);
  XSetWindowBorderWidth(xwin->display, xwin->window, BORDERWIDTH) ;
  XSetLineAttributes(xwin->display,xwin->black,
                     3, LineSolid, CapRound, JoinRound );
  XSetLineAttributes(xwin->display,xwin->white,
                     3, LineSolid, CapRound, JoinRound );
  XSetLineAttributes(xwin->display,xwin->flip,
                     3, LineSolid, CapRound, JoinRound );

  if (!ncolors) xUseGrayscale(xwin->display, xwin->screen) ;
  XSync(xwin->display, 0) ;
}
/*************************************************************************/
/**     open_xwindow                                                    **/
/** arguments:                                                          **/
/**     struct xwindow_type *xwindow    xwindow structure               **/
/**     int xpos,ypos                   window position- negative value **/
/**                                       means position with mouse     **/
/**     int xsize,ysize                 window dimensions               **/
/**     char title[]                    window and icon title           **/
/**     long eventmask                  event mask e.g. ButtonPressMask **/
/**                                                     KeyPressMask    **/
/**                                                     ExposureMask    **/
/**     Steve Lehar (slehar@park.bu.edu)                                **/
/**     Greg Lesher (lesher@park.bu.edu)                                **/
/*************************************************************************/
void
open_xwindow(xwin,displayname,xpos,ypos,xsize,ysize,title,eventmask)
xwindow_type *xwin;
char displayname[];
int xpos,ypos,xsize,ysize;
char title[];
long eventmask;
{
  GC create_gc();
  XWMHints xwmh;
  unsigned long valuemask ;
  XSetWindowAttributes xswa ;
  XVisualInfo xinfo ;
  int ret ;

  /**** open display and screen ****/
  xwin->display = XOpenDisplay(displayname);
  if (xwin->display == NULL)
  {
    printf("ERROR open_xwindow: can't open display %s\n",displayname);
    exit(1);
  }
  xwin->screen  = DefaultScreen(xwin->display);

  /**** set window position and size hints for window manager ****/
  xwin->hint.x      = xpos;
  xwin->hint.y      = ypos;
  xwin->hint.width  = xsize;
  xwin->hint.height = ysize;
  xwin->xsize  = xsize;
  xwin->ysize = ysize;
  if (xpos < 0 || ypos < 0)
    xwin->hint.flags = ( PPosition | PSize);
  else
    xwin->hint.flags = (USPosition | PSize);

#if 0
  /**** create the x window ****/
  xwin->window = XCreateSimpleWindow(xwin->display,
                                     DefaultRootWindow(xwin->display),
                                     xwin->hint.x,     xwin->hint.y,
                                     xwin->hint.width, xwin->hint.height,
                                     BORDERWIDTH,
                                     BlackPixel(xwin->display,xwin->screen),
                                     WhitePixel(xwin->display,xwin->screen));
#else

  ret = XMatchVisualInfo(xwin->display, DefaultScreen(xwin->display),
                         8, GrayScale, &xinfo) ;

  /*  if (ret == False)  always use default visual */
  xinfo.visual = DefaultVisual(xwin->display, xwin->screen) ;
  valuemask = CWEventMask ;
  xswa.event_mask = eventmask ;
  xwin->window = XCreateWindow(xwin->display,
                               DefaultRootWindow(xwin->display),
                               xwin->hint.x,     xwin->hint.y,
                               xwin->hint.width, xwin->hint.height,
                               BORDERWIDTH, CopyFromParent, InputOutput,
                               xinfo.visual,
                               valuemask, &xswa) ;
  xwmh.flags = InputHint ;
  xwmh.input = True ;
  XSetWMHints(xwin->display, xwin->window, &xwmh) ;

#endif

  XSetStandardProperties(xwin->display, xwin->window, title, title,
                         None, NULL, None, &xwin->hint);

  /**** make window sensitive selected events ****/
  XSelectInput(xwin->display, xwin->window, eventmask);

  /**** create graphics contexts (drawing colors & properties) ****/
  xwin->black = XCreateGC(xwin->display, xwin->window, 0, 0);
  xwin->white = XCreateGC(xwin->display, xwin->window, 0, 0);
  xwin->flip  = XCreateGC(xwin->display, xwin->window, 0, 0);

  /**** set gc drawing functions ****/
  XSetFunction(xwin->display, xwin->black, GXcopy);
  XSetFunction(xwin->display, xwin->white, GXcopy);
  XSetFunction(xwin->display, xwin->flip,  GXinvert);

  /**** set gc drawing colors ****/
  XSetForeground(xwin->display, xwin->black,
                 BlackPixel(xwin->display,xwin->screen));
  XSetBackground(xwin->display, xwin->black,
                 WhitePixel(xwin->display,xwin->screen));
  XSetForeground(xwin->display, xwin->white,
                 WhitePixel(xwin->display,xwin->screen));
  XSetBackground(xwin->display, xwin->white,
                 BlackPixel(xwin->display,xwin->screen));

  /**** load font ****/
  xwin->font_struct = XLoadQueryFont(xwin->display, FONT);
  if (xwin->font_struct == 0)
  {
    printf("font %s not found, using default\n",FONT);
    xwin->font_struct =
      XQueryFont(xwin->display,
                 XGContextFromGC(DefaultGC(xwin->display,xwin->screen)));
  }
  else
  {
    XSetFont(xwin->display, xwin->black, xwin->font_struct->fid);
    XSetFont(xwin->display, xwin->white, xwin->font_struct->fid);
    XSetFont(xwin->display, xwin->flip,  xwin->font_struct->fid);
  }

  /**** get font character dimensions ****/
  XGetFontProperty(xwin->font_struct, XA_QUAD_WIDTH, &xwin->charwidth);
  XGetFontProperty(xwin->font_struct, XA_CAP_HEIGHT, &xwin->charheight);

  /**** map window to screen (make it visible) ****/
  XMapRaised(xwin->display, xwin->window);
}

/*************************************************************************/
/**     close_xwindow                                                   **/
/**  Close the x window                                                 **/
/**     Steve Lehar (slehar@park.bu.edu)                                **/
/*************************************************************************/
void
close_xwindow(xwin)
xwindow_type *xwin;
{
  XDestroyWindow(xwin->display,xwin->window);
  XCloseDisplay(xwin->display);
}

void
xNextEvent(xwindow_type *xwin, XEvent *event)
{
  if (!xwin)
  {
    static int callno = 0 ;

    switch (callno)
    {
    case 0:
      break ;
    case 1:
      event->type = ButtonPress ;
      event->xbutton.button = Button1 ;
      event->xbutton.x = xwin->hint.width / 3 ;
      event->xbutton.y = 0 ;
      break ;
    default:
      break ;
    }
    return ;
  }

#if 0
  /* KLUGE a GraphicsExpose event */
  /* to keep XNextEvent from blocking */
  XCopyArea(xwin->display,xwin->window,xwin->window,
            xwin->white,0,0,1,1,0,0);
#else
  if (XEventsQueued(xwin->display, QueuedAfterReading))
#endif
  XNextEvent(xwin->display, event) ;
  event->xbutton.x -= xwin->xsize / 2 ;
  event->xbutton.y = xwin->ysize / 2 - event->xbutton.y ;
}

void
xDrawLine(xwindow_type *xwin, int x0, int y0, int x1, int y1, int color,
          int style)
{
  int ox, oy ;
  GC  gc ;

  if (!xwin) return ;

  gc = xwin->black ;  /* default */
  x0 *= xwin->scale ;
  y0 *= xwin->scale ;
  x1 *= xwin->scale ;
  y1 *= xwin->scale ;

  ox = xwin->xsize / 2 ;
  oy = xwin->ysize / 2 ;
#if 0
  x0 += ox ;
  x1 += ox ;
  y1 = oy - y1 ;
  y0 = oy - y0 ;
#endif

  switch (style)
  {
  case PLOT_LINES_DASHED:
    style = LineOnOffDash ;
    break ;
  case PLOT_LINES_DOUBLE_DASHED:
    style = LineDoubleDash ;
  default:
    style = LineSolid ;
    break ;
  }

  switch (color)
  {
  case xBLACK:
    gc = xwin->black ;
    break ;
  case xFLIP:
    gc = xwin->flip ;
    break ;
  case xWHITE:
    gc = xwin->white ;
    break ;
  default:
    break ;
  }

  XSetLineAttributes(xwin->display, gc, 0, style, CapRound, JoinBevel) ;
  XDrawLine(xwin->display,xwin->window, gc, x0, y0, x1, y1) ;
}
void
xDrawLines(xwindow_type *xwin, XSegment *segs, int nsegs, int color)
{
  int i ;
  GC  gc ;

  if (!xwin) return ;

  /* scale according to window size */
  for (i = 0 ; i < nsegs ; i++)
  {
    segs[i].x1 *= xwin->scale ;
    segs[i].y1 *= xwin->scale ;
    segs[i].x2 *= xwin->scale ;
    segs[i].y2 *= xwin->scale ;
  }

  switch (color)
  {
  case xBLACK:
    gc = xwin->black ;
    break ;
  case xFLIP:
    gc = xwin->flip ;
    break ;
  case xWHITE:
    gc = xwin->white ;
    break ;
  default:
    gc = xwin->black ;
    break ;
  }
  XDrawSegments(xwin->display, xwin->window, gc, segs, nsegs) ;
}

void
xDrawCircle(xwindow_type *xwin, int x0, int y0, int radius, int color)
{
  int ox, oy ;

  if (!xwin) return ;

  x0 *= xwin->scale ;
  y0 *= xwin->scale ;
  radius *= xwin->scale ;

  ox = xwin->xsize / 2 ;
  oy = xwin->ysize / 2 ;
#if 0
  x0 += ox ;
  y0 = oy - y0 ;
#endif

  switch (color)
  {
  case xBLACK:
    XDrawArc(xwin->display, xwin->window, xwin->black,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  case xWHITE:
    XDrawArc(xwin->display, xwin->window, xwin->white,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  case xFLIP:
    XDrawArc(xwin->display, xwin->window, xwin->flip,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  default:
    break ;
  }
#if 0
  XSync(xwin->display, 0) ;
#endif
}

#if 0
void
xDrawCircles(xwindow_type *xwin, XArc *xcircles, int ncircles, int color)
{
  int i ;

  if (!xwin) return ;

  for (i = 0 ; i < ncircles ; i++)
  {
    x0 *= xwin->scale ;
    y0 *= xwin->scale ;
    radius *= xwin->scale ;
  }

  switch (color)
  {
  case xBLACK:
    XDrawArc(xwin->display, xwin->window, xwin->black,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  case xWHITE:
    XDrawArc(xwin->display, xwin->window, xwin->white,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  case xFLIP:
    XDrawArc(xwin->display, xwin->window, xwin->flip,
             x0 - radius, y0 -  radius, 2*radius, 2*radius, 0, FULLCIRCLE) ;
    break ;
  default:
    break ;
  }
}
#endif

void
xFreeWindow(xwindow_type *xwin)
{
  if (xwin->black) XFreeGC(xwin->display, xwin->black) ;
  if (xwin->white) XFreeGC(xwin->display, xwin->white) ;
  if (xwin->flip) XFreeGC(xwin->display, xwin->flip) ;
  XDestroyWindow(xwin->display, xwin->window) ;
  free(xwin) ;
}

xwindow_type *
xNewWindow(xwindow_type *parent, int xpos, int ypos, int xsize, int ysize,
           char *name, int border, int event_mask)
{
  xwindow_type *xwin ;
  Display      *display ;
  int          screen ;
  Window       pwindow ;

  xwin = (xwindow_type *)calloc(1, sizeof(xwindow_type)) ;
  xwin->scale = 1 ;

  if (!parent)
  {
    display = XOpenDisplay("") ;
    screen = DefaultScreen(display) ;
    pwindow = DefaultRootWindow(display) ;
  }
  else
  {
    display = parent->display ;
    screen = parent->screen ;
    pwindow = parent->window ;
  }
  xwin->display = display ;
  xwin->screen = screen ;
  xwin->hint.x = xpos ;
  xwin->hint.y = ypos ;
  xwin->xsize = xwin->hint.width = xsize ;
  xwin->ysize = xwin->hint.height = ysize ;
  xwin->hint.flags = (PPosition | PSize) ;
  xwin->screen = screen ;
  xwin->display = display ;

  xwin->window = XCreateSimpleWindow(display,
                                     pwindow,
                                     xwin->hint.x, xwin->hint.y,
                                     xwin->hint.width, xwin->hint.height,
                                     border,
                                     BlackPixel(display, screen),
                                     WhitePixel(display, screen));
  XSetStandardProperties(display, xwin->window, name, name,
                         None, NULL, None, &xwin->hint) ;


  /* create graphics contexts */
  xwin->white = XCreateGC(display, xwin->window, 0, 0) ;
  xwin->black = XCreateGC(xwin->display, xwin->window, 0, 0);
  XSetFunction(xwin->display, xwin->white, GXcopy);
  XSetFunction(xwin->display, xwin->black, GXcopy);

  XSetBackground(display, xwin->white, BlackPixel(display, screen)) ;
  XSetForeground(display, xwin->white, WhitePixel(display, screen)) ;
  XSetForeground(xwin->display, xwin->black,
                 BlackPixel(xwin->display,xwin->screen));
  XSetBackground(xwin->display, xwin->black,
                 WhitePixel(xwin->display,xwin->screen));

  XSelectInput(display, xwin->window,
               ButtonReleaseMask | ButtonPressMask | KeyPressMask |
               ExposureMask | KeyReleaseMask) ;

  /**** load font ****/
  xwin->font_struct = XLoadQueryFont(xwin->display, FONT);
  if (xwin->font_struct == 0)
  {
    printf("font %s not found, using default\n",FONT);
    xwin->font_struct =
      XQueryFont(xwin->display,
                 XGContextFromGC(DefaultGC(xwin->display,xwin->screen)));
  }
  else
  {
    XSetFont(xwin->display, xwin->white, xwin->font_struct->fid);
#if 0
    XSetFont(xwin->display, xwin->black, xwin->font_struct->fid);
    XSetFont(xwin->display, xwin->flip,  xwin->font_struct->fid);
#endif
  }

  /**** get font character dimensions ****/
  XGetFontProperty(xwin->font_struct, XA_QUAD_WIDTH, &xwin->charwidth);
  XGetFontProperty(xwin->font_struct, XA_CAP_HEIGHT, &xwin->charheight);

  return(xwin) ;
}
void
xDrawPoint(xwindow_type *xwin, int x, int y, int color)
{
  GC gc ;
  int  xstart, ystart, xend, yend ;

  xstart = x * xwin->scale ;
  xend = xstart + xwin->scale - 1 ;
  ystart = y * xwin->scale ;
  yend = ystart + xwin->scale - 1 ;

  switch (color)
  {
  default:
  case xBLACK:
    gc = xwin->black ;
    break ;
  case xWHITE:
    gc = xwin->white ;
    break ;
  case xFLIP:
    gc = xwin->flip ;
  }

  for (x = xstart ; x <= xend ; x++)
    for (y = ystart ; y <= yend ; y++)
      XDrawPoint(xwin->display, xwin->window, gc, x, y) ;
}
void
xShowWindowAttributes(xwindow_type *xwin, FILE *fp)
{
  XWindowAttributes xattribs ;

  XGetWindowAttributes(xwin->display, xwin->window, &xattribs) ;

  fprintf(fp, "xwin @%lx, x %d, y %d, width %d, height %d, border %d\n",
          (long)xwin, xattribs.x, xattribs.y, xattribs.width, xattribs.height,
          xattribs.border_width) ;
  fprintf(fp, "class %s, event_mask %lx:\n",
          xattribs.class == InputOnly ? "InputOnly" : "InputOutput",
          xattribs.your_event_mask) ;
  xShowEventMask(xattribs.your_event_mask, fp) ;

#if 0
  fprintf(fp, "all events %lx:\n", xattribs.all_event_masks) ;
  xShowEventMask(xattribs.all_event_masks, fp) ;
#endif
}

void
xShowEventMask(unsigned long eventmask, FILE *fp)
{
  if (eventmask & KeyPressMask)
    fprintf(fp, "\tKeyPress\n") ;
  if (eventmask & KeyReleaseMask)
    fprintf(fp, "\tKeyRelease\n") ;
  if (eventmask & ButtonPressMask)
    fprintf(fp, "\tButtonPress\n") ;
  if (eventmask & ButtonReleaseMask)
    fprintf(fp, "\tButtonRelease\n") ;
  if (eventmask & EnterWindowMask)
    fprintf(fp, "\tEnterWindow\n") ;
  if (eventmask & LeaveWindowMask)
    fprintf(fp, "\tLeaveWindow\n") ;
  if (eventmask & PointerMotionMask)
    fprintf(fp, "\tPointerMotion\n") ;
  if (eventmask & PointerMotionHintMask)
    fprintf(fp, "\tPointerMotionHint\n") ;
  if (eventmask & Button1MotionMask)
    fprintf(fp, "\tButton1Motion\n") ;
  if (eventmask & Button2MotionMask)
    fprintf(fp, "\tButton2Motion\n") ;
  if (eventmask & Button3MotionMask)
    fprintf(fp, "\tButton3Motion\n") ;
  if (eventmask & Button4MotionMask)
    fprintf(fp, "\tButton4Motion\n") ;
  if (eventmask & Button5MotionMask)
    fprintf(fp, "\tButton5Motion\n") ;
  if (eventmask & ButtonMotionMask)
    fprintf(fp, "\tButtonMotion\n") ;
  if (eventmask & KeymapStateMask)
    fprintf(fp, "\tKeymapState\n") ;
  if (eventmask & ExposureMask)
    fprintf(fp, "\tExposure\n") ;
  if (eventmask & VisibilityChangeMask)
    fprintf(fp, "\tVisibilityChange\n") ;
  if (eventmask & StructureNotifyMask)
    fprintf(fp, "\tStructureNotify\n") ;
  if (eventmask & ResizeRedirectMask)
    fprintf(fp, "\tResizeRedirect\n") ;
  if (eventmask & SubstructureNotifyMask)
    fprintf(fp, "\tSubstructureNotify\n") ;
  if (eventmask & SubstructureRedirectMask)
    fprintf(fp, "\tSubstructureRedirect\n") ;
  if (eventmask & FocusChangeMask)
    fprintf(fp, "\tFocusChange\n") ;
  if (eventmask & PropertyChangeMask)
    fprintf(fp, "\tPropertyChange\n") ;
  if (eventmask & ColormapChangeMask)
    fprintf(fp, "\tColormapChange\n") ;
  if (eventmask & OwnerGrabButtonMask)
    fprintf(fp, "\tOwnerGrabButton\n") ;
}
void
xShowEvent(unsigned long event, FILE *fp)
{
  switch (event)
  {
  case KeyPress:
    fprintf(fp, "KeyPress\n") ;
    break ;
  case KeyRelease:
    fprintf(fp, "KeyRelease\n") ;
    break ;
  case ButtonPress:
    fprintf(fp, "ButtonPress\n") ;
    break ;
  case ButtonRelease:
    fprintf(fp, "ButtonRelease\n") ;
    break ;
  case MotionNotify:
    fprintf(fp, "MotionNotify\n") ;
    break ;
  case EnterNotify:
    fprintf(fp, "EnterNotify\n") ;
    break ;
  case LeaveNotify:
    fprintf(fp, "LeaveNotify\n") ;
    break ;
  case FocusIn:
    fprintf(fp, "FocusIn\n") ;
    break ;
  case FocusOut:
    fprintf(fp, "FocusOut\n") ;
    break ;
  case KeymapNotify:
    fprintf(fp, "KeymapNotify\n") ;
    break ;
  case Expose:
    fprintf(fp, "Expose\n") ;
    break ;
  case GraphicsExpose:
    fprintf(fp, "GraphicsExpose\n") ;
    break ;
  case NoExpose:
    fprintf(fp, "NoExpose\n") ;
    break ;
  case VisibilityNotify:
    fprintf(fp, "VisibilityNotify\n") ;
    break ;
  case CreateNotify:
    fprintf(fp, "CreateNotify\n") ;
    break ;
  case DestroyNotify:
    fprintf(fp, "DestroyNotify\n") ;
    break ;
  case UnmapNotify:
    fprintf(fp, "UnmapNotify\n") ;
    break ;
  case MapNotify:
    fprintf(fp, "MapNotify\n") ;
    break ;
  case MapRequest:
    fprintf(fp, "MapRequest\n") ;
    break ;
  case ReparentNotify:
    fprintf(fp, "ReparentNotify\n") ;
    break ;
  case ConfigureNotify:
    fprintf(fp, "ConfigureNotify\n") ;
    break ;
  case ConfigureRequest:
    fprintf(fp, "ConfigureRequest\n") ;
    break ;
  case GravityNotify:
    fprintf(fp, "GravityNotify\n") ;
    break ;
  case ResizeRequest:
    fprintf(fp, "ResizeRequest\n") ;
    break ;
  case CirculateNotify:
    fprintf(fp, "CirculateNotify\n") ;
    break ;
  case CirculateRequest:
    fprintf(fp, "CirculateRequest\n") ;
    break ;
  case PropertyNotify:
    fprintf(fp, "PropertyNotify\n") ;
    break ;
  case SelectionClear:
    fprintf(fp, "SelectionClear\n") ;
    break ;
  case SelectionRequest:
    fprintf(fp, "SelectionRequest\n") ;
    break ;
  case SelectionNotify:
    fprintf(fp, "SelectionNotify\n") ;
    break ;
  case ColormapNotify:
    fprintf(fp, "ColormapNotify\n") ;
    break ;
  case ClientMessage:
    fprintf(fp, "ClientMessage\n") ;
    break ;
  case MappingNotify:
    fprintf(fp, "MappingNotify\n") ;
    break ;
  case LASTEvent:
    fprintf(fp, "LASTEvent\n") ;
    break ;
  }
}
void
xShowDisplayInfo(xwindow_type *xwin, FILE *fp)
{
  int   /* ret, */ major_ver, minor_rev, release, cells, planes ;
  char *display_name, *vendor_name ;
  Visual *visual ;

  display_name = DisplayString(xwin->display) ;

  cells = DisplayCells(xwin->display, xwin->screen) ;
  planes = DisplayPlanes(xwin->display, xwin->screen) ;
  major_ver = ProtocolVersion(xwin->display) ;
  minor_rev = ProtocolRevision(xwin->display) ;
  release = VendorRelease(xwin->display) ;
  vendor_name = ServerVendor(xwin->display) ;

  fprintf(fp, "x version %d.%d.%d from %s\n", major_ver, minor_rev,
          release, vendor_name) ;
  fprintf(fp, "display %s, %d planes, %d cells\n", display_name,  planes,
          cells) ;
  visual = DefaultVisual(xwin->display, xwin->screen) ;
  fprintf(fp, "visual class ") ;
  switch (visual->class)
  {
  case PseudoColor:
    fprintf(fp, "%s\n", "PseudoColor") ;
    break ;
  case DirectColor:
    fprintf(fp, "%s\n", "DirectColor") ;
    break ;
  case GrayScale:
    fprintf(fp, "%s\n", "GrayScale") ;
    break ;
  case StaticColor:
    fprintf(fp, "%s\n", "StaticColor") ;
    break ;
  case TrueColor:
    fprintf(fp, "%s\n", "TrueColor") ;
    break ;
  case StaticGray:
    fprintf(fp, "%s\n", "StaticGray") ;
    break ;
  }

#if 0
  if (ret)
    fprintf(fp, "GrayScale available\n") ;
  else
    fprintf(fp, "no GrayScale available\n") ;
#endif
}

#define NPLANES 0
#define SHIFT 8


int
xUseGrayscale(Display *display, int screen)
{
  int            i, result, cindex ;
  double         cno, slotsPerColor ;
  unsigned long  plane_masks[NPLANES+10], colors[NCOLORS+10] ;
  Colormap       cmap ;

  cmap = DefaultColormap(display, screen) ;

  ncolors = NCOLORS ;
  do
  {
    result = XAllocColorCells(display, cmap,
                              False, plane_masks, NPLANES,
                              colors, ncolors) ;
    if (result == False)
      ncolors-- ;
  }
  while (result == False) ;

  if (result == False)
  {
    fprintf(stderr, "could not allocate color cells\n") ;
    exit(1) ;
  }

  slotsPerColor = (double)ncolors / (double)NCOLORS ;
  for (cno = 0.0, i = 0 ; i < NCOLORS ; i++, cno += slotsPerColor)
  {
    cindex = nint(cno) ;
    xcolors[i].pixel = colors[cindex] ;
    xcolors[i].flags = DoRed | DoGreen | DoBlue ;
    xcolors[i].red = cindex << SHIFT ;
    xcolors[i].green = cindex << SHIFT ;
    xcolors[i].blue = cindex << SHIFT ;
  }
  XStoreColors(display, cmap, xcolors, NCOLORS) ;
  printf("%d gray scale values allocated\n", ncolors) ;

  return(ncolors) ;
}
#define MAX_STRING 250

int
xprintf(xwindow_type *xwin, int x, int y, char *fmt, ...)
{
  va_list args ;
  char    str[MAX_STRING] ;
  int     len ;


  va_start(args, fmt) ;
#if 0
  xwin = va_arg(args, xwindow_type *) ;
  x = va_arg(args, int) ;
  y = va_arg(args, int) ;
  fmt = va_arg(args, char *) ;
#endif

  /* this will clear the whole window from x, y down and to the right */
  XClearArea(xwin->display, xwin->window, 0, 0, 0, 0, False) ;

  vsprintf(str, fmt, args) ;
  len = strlen(str) ;
  XDrawString(xwin->display, xwin->window, xwin->black, x, y, str, len);
  /*  XFlush(xwin->display) ;*/
  va_end(args) ;
  return(len) ;
}
