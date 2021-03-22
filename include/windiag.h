/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
   @(#)window.h 1.4
   8/10/95
*/
/*------------------------------------------------------------------------
      File Name:  window.h

         Author:  Bruce Fischl

        Created:  Feb. 1994

    Description:

------------------------------------------------------------------------*/
#ifndef WINDOW_H
#define WINDOW_H


#include <stdarg.h>

#include "image.h"

int  WinAlloc(char *pcName, int iXpos, int iYpos, int iWidth, int iHeight) ;
int  WinShow(int iWin) ;
int  WinFree(int iWin) ;
int  WinDrawLine(int iWin, int x0, int y0, int x1, int y1,int color,int style);
int  WinClear(int iWin) ;
int  WinClearArea(int iWin, int x0, int y0, int width, int height) ;
int  WinFlush(int iWin) ;
int  WinDrawCircle(int iWin, int x0, int y0, int radius, int color) ;
int  WinSetScale(int iWin, double dXscale, double dYscale) ;
int  WinGetScale(int iWin, double *pdXscale, double *pdYscale) ;
int  WinSetRange(int iWin,double dXmin,double dXmax,double dYmin,double dYmax);
int  WinGetRange(int iWin, double *pdXmin, double *pdXmax, double *pdYmin,
                 double *pdYmax);

int  WinCreate(char *pcName, int button_rows, int image_rows,
               int image_cols, int rows, int cols) ;
int  WinPrintf(int iWin, int x, int y, char *fmt, ...) ;
int  WinShowImage(int iWin, IMAGE *image, int which) ;
int  WinSetName(int win, int which, char *fmt, ...) ;

#endif
