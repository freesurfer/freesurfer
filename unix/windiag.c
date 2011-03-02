/**
 * @file  windiag.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:41 $
 *    $Revision: 1.10 $
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
   @(#)window.c 1.1
   2/28/94
*/
/*------------------------------------------------------------------------
      File Name: window.c

         Author: Bruce Fischl

        Created: Jan. 1994

    Description:

------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdarg.h>
#include <math.h>

#include <xview/notify.h>

#include "windiag.h"
#include "error.h"
#include "proto.h"
#include "xwin.h"
#include "image.h"
#include "xvutil.h"
#include "thread.h"


/*------------------------------------------------------------------------
                            CONSTANTS
------------------------------------------------------------------------*/

#define MAX_WINDOWS   10
#define MIN_VAL       0
#define MAX_VAL       255

/*------------------------------------------------------------------------
                            STRUCTURES
------------------------------------------------------------------------*/

typedef struct
{
  XV_FRAME     *xvf ;
  xwindow_type *xwin ;
  double       dXscale ;
  double       dYscale ;
  double       dXmin ;
  double       dYmin ;
  double       dXmax ;
  double       dYmax ;
}
DIAG_WINDOW ;


/*------------------------------------------------------------------------
                            GLOBAL DATA
------------------------------------------------------------------------*/

static int inited = 0 ;

/*------------------------------------------------------------------------
                            STATIC DATA
------------------------------------------------------------------------*/

static DIAG_WINDOW window_table[MAX_WINDOWS] ;
static DIAG_WINDOW *HandleToPtr(int iWin) ;
static void  event_handler(Event *event, DIMAGE *dimage) ;

/*------------------------------------------------------------------------
                            STATIC PROTOTYPES
------------------------------------------------------------------------*/

static int            winNewHandle(void) ;
static Notify_value   winPoll(void) ;
static void           winThread(int iTid, void *parm) ;

/*------------------------------------------------------------------------
                              FUNCTIONS
------------------------------------------------------------------------*/

/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinShow(int iWin)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  XMapRaised(pwin->xvf->display, pwin->xwin->window) ;
  XFlush(pwin->xvf->display) ;
  return(1) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
#define MAX_RANGE   200
int
WinShowImage(int iWin, IMAGE *image, int which)
{
  DIAG_WINDOW       *pwin ;

  pwin = HandleToPtr(iWin) ;
  XVshowImage(pwin->xvf, which, image, 0) ;
  ThreadYield() ;
  return(NO_ERROR) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinSetName(int win, int which, char *fmt, ...)
{
  char   name[100] ;
  DIAG_WINDOW *pwin ;
  va_list args ;

  va_start(args, fmt) ;
  /*  fmt = va_arg(args, char *) ;*/
  vsprintf(name, fmt, args) ;

  pwin = HandleToPtr(win) ;
  XVshowImageTitle(pwin->xvf, which, name) ;

  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinClear(int iWin)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;

  /* this will clear the whole window from x, y down and to the right */
  XClearArea(pwin->xvf->display, pwin->xwin->window, 0, 0, 0, 0, False) ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinFree(int iWin)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  xFreeWindow(pwin->xwin) ;
  pwin->xwin = NULL ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinFlush(int iWin)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  XSync(pwin->xvf->display, 0) ;
  XFlush(pwin->xvf->display) ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinAlloc(char *pcName, int iXpos, int iYpos, int iWidth, int iHeight)
{
  int    iWin ;
  DIAG_WINDOW *pwin ;

  iWin = winNewHandle() ;
  if (iWin < 0)
    return(ErrorPrintf(ERROR_NO_MEMORY,
                       "WinAlloc: could not allocate new window"));

  pwin = HandleToPtr(iWin) ;
#if 1
  pwin->xwin = xNewWindow(NULL, iXpos, iYpos, iWidth, iHeight, pcName, 0, 0) ;
#else
  pwin->xwin = xNewWindow(NULL, iXpos, iYpos, iWidth, iHeight, pcName, 0, 0) ;
#endif

  return(iWin) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
static Notify_value
winPoll(void)
{
  static int ncalls = 0 ;

  if (ncalls++ == 4)
    ThreadSignal(0, SIG_ALL) ;

  ThreadYield() ;
  /*  ThreadSleep(TID_SELF, 100) ;*/
  return(NOTIFY_DONE) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinCreate(char *pcName, int button_rows, int image_rows, int image_cols,
          int rows, int cols)
{
  int    iWin ;
  DIAG_WINDOW *pwin ;

  iWin = winNewHandle() ;
  if (iWin < 0)
    return(ErrorPrintf(ERROR_NO_MEMORY,
                       "WinAlloc: could not allocate new window"));

  pwin = HandleToPtr(iWin) ;
  pwin->xvf = XValloc(rows, cols, button_rows, image_rows, image_cols,
                      pcName, winPoll) ;
  XVsetParms(event_handler) ;

  if (!inited)
  {
    ThreadInit(0, MAX_WINDOWS, 10*1024, 1) ;
    inited = 1 ;
  }

  ThreadStart("window", winThread, pwin, MIN_PRIORITY) ;

  ThreadSuspend(TID_SELF, SIG_ANY) ;
  WinFlush(iWin) ;
  return(iWin) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
------------------------------------------------------------------------*/
static void
winThread(int iTid, void *parm)
{
  DIAG_WINDOW  *pwin ;

  pwin = (DIAG_WINDOW *)parm ;
  fprintf(stderr, "entering main loop\n") ;
  xv_main_loop(pwin->xvf->frame) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
static int
winNewHandle(void)
{
  int iWin ;

  for (iWin = 0 ; iWin < MAX_WINDOWS ; iWin++)
    if (!window_table[iWin].xwin)
      return(iWin) ;

  return(-1) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
static DIAG_WINDOW *
HandleToPtr(int iWin)
{
  if ((iWin < 0) || (iWin >= MAX_WINDOWS))
  {
    ErrorPrintf(ERROR_NO_MEMORY, "HandleToPtr(%d): bad handle", iWin) ;
    return(NULL) ;
  }
  return(&window_table[iWin]) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinDrawLine(int iWin, int x0, int y0, int x1, int y1, int color, int style)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;

  /* invert coordinate system so increasing y goes up */
  y0 = pwin->xwin->ysize - y0 ;
  y1 = pwin->xwin->ysize - y1 ;
  xDrawLine(pwin->xwin, x0, y0, x1, y1, color, style) ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
#define MAX_STRING 250
int
WinPrintf(int iWin, int x, int y, char *fmt, ...)
{
  va_list args ;
  char    str[MAX_STRING] ;
  int     len ;
  DIAG_WINDOW  *pwin ;

  va_start(args, fmt) ;
  pwin = HandleToPtr(iWin) ;

  /* AHH! what coordinate system to use? */
  /*  y = pwin->xwin->ysize - y ;*/
  fmt = va_arg(args, char *) ;

  vsprintf(str, fmt, args) ;
  len = strlen(str) ;
  XDrawString(pwin->xvf->display, pwin->xwin->window, pwin->xwin->black,
              x, y, str, len);
  /*  XFlush(pwin->xvf->display) ;*/
  va_end(args) ;
  return(len) ;
}

int
WinClearArea(int iWin, int x0, int y0, int width, int height)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  XClearArea(pwin->xvf->display, pwin->xwin->window,x0,y0,width,height,False);
  return(1) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinDrawCircle(int iWin, int x0, int y0, int radius, int color)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;

  /* invert coordinate system so increasing y goes up */
  /* AHH! what coordinate system to use? */
  /*  y0 = pwin->xwin->ysize - y0 ;*/
  xDrawCircle(pwin->xwin, x0, y0, radius, color) ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinSetScale(int iWin, double dXscale, double dYscale)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;

  pwin->dXscale = dXscale ;
  pwin->dYscale = dYscale ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinGetScale(int iWin, double *pdXscale, double *pdYscale)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  *pdXscale = pwin->dXscale ;
  *pdYscale = pwin->dYscale ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinSetRange(int iWin, double dXmin, double dXmax, double dYmin, double dYmax)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  pwin->dXmin = dXmin ;
  pwin->dXmax = dXmax ;
  pwin->dYmin = dYmin ;
  pwin->dYmax = dYmax ;
  return(0) ;
}
/*------------------------------------------------------------------------
       Parameters:

      Description:

    Return Values:
            0 on success, < 0 otherwise.

------------------------------------------------------------------------*/
int
WinGetRange(int iWin, double *pdXmin, double *pdXmax, double *pdYmin,
            double *pdYmax)
{
  DIAG_WINDOW *pwin ;

  pwin = HandleToPtr(iWin) ;
  *pdXmin = pwin->dXmin ;
  *pdXmax = pwin->dXmax ;
  *pdYmin = pwin->dYmin ;
  *pdYmax = pwin->dYmax ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
event_handler(Event *event, DIMAGE *dimage)
{
#if 0
  int    x, y ;
  float  scale, dx, dy ;

  x = event_x(event) ;
  y = event_y(event) ;

  scale = dimage->scale ;
  if (event_id(event) == MS_MIDDLE)
  {
    if (event_is_down(event))
    {
      if (processed)
      {
        dx = *HIPSFpix(offsetImage, x, y) ;
        dy = *HIPSFseq_pix(offsetImage, x, y, 1) ;
        XVprintf(xvf, 0, "offset at (%d, %d) = (%2.3f, %2.3f)", x, y, dx, dy) ;
      }
      else
        dx = dy = 0.0f ;

      showOffsetArea(x,y, dx, dy) ;
      switch (filter_type)
      {
      case FILTER_EXPONENTIAL:
        calculateExponentialKernel(gradImage, x, y, filter_size,
                                   filter_parm, kernelImage, dx, dy) ;
        break ;
      case FILTER_EXP_LAPLACIAN:
        calculateExponentialKernel(laplacianImage, x, y, filter_size,
                                   filter_parm, kernelImage, dx, dy) ;
        break ;
      case FILTER_EXP_SUM:
        tmpImage = HipsAdd(laplacianImage, gradImage, tmpImage) ;
        calculateExponentialKernel(tmpImage, x, y, filter_size,
                                   filter_parm, kernelImage, dx, dy) ;
        break ;
      default:
        return ;
      }

      XVshowImage(xvf, FILTER_IMAGE, kernelImage, 0) ;
      XVshowImageTitle(xvf, FILTER_IMAGE, "filter at (%d, %d)    ", x, y) ;
    }
  }
#endif
}
