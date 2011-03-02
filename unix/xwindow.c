/**
 * @file  xwindow.c
 * @brief XWindow crap
 *
 */
/*
 * Original Author: Bruce Fischl
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

#ifdef HAVE_OPENGL

/*############################################################################*/
/* orig: window.c: marty: tk=>tko,TK_=>TKO_ in all to avoid tcl/tk */
/*                        add definitions of predeclared (event.c,cursor.c) */
/*                        add 3 functions from getset.c */
/*                        InitWindow: add KeyRelease Event */
/*                        wrote tkoInitPixmap() for offscreen render */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "proto.h"
#include "xwindow.h"  /* marty: = tk.h + private.h */

/******************************************************************************/

int print_all_glx_visuals = 0;  /* marty */

Display *xDisplay = 0;
int xScreen = 0;
Window wRoot = 0;
Atom deleteWindowAtom;
WINDOW_REC w = {
                 0, 0, 300, 300, TKO_RGB|TKO_SINGLE|TKO_DIRECT, TKO_MINIMUM_CRITERIA
               };
float colorMaps[] = {
                      0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000,
                      0.000000, 1.000000, 0.333333, 0.776471, 0.443137, 0.556863,
                      0.443137, 0.556863, 0.219608, 0.666667, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.039216, 0.078431, 0.117647, 0.156863,
                      0.200000, 0.239216, 0.278431, 0.317647, 0.356863, 0.400000,
                      0.439216, 0.478431, 0.517647, 0.556863, 0.600000, 0.639216,
                      0.678431, 0.717647, 0.756863, 0.800000, 0.839216, 0.878431,
                      0.917647, 0.956863, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 0.000000, 0.000000,
                      1.000000, 1.000000, 0.000000, 0.000000, 1.000000, 1.000000,
                      0.333333, 0.443137, 0.776471, 0.556863, 0.443137, 0.219608,
                      0.556863, 0.666667, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.039216, 0.078431, 0.117647, 0.156863, 0.200000, 0.239216,
                      0.278431, 0.317647, 0.356863, 0.400000, 0.439216, 0.478431,
                      0.517647, 0.556863, 0.600000, 0.639216, 0.678431, 0.717647,
                      0.756863, 0.800000, 0.839216, 0.878431, 0.917647, 0.956863,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.141176, 0.282353, 0.427451,
                      0.568627, 0.713726, 0.854902, 1.000000, 0.000000, 0.141176,
                      0.282353, 0.427451, 0.568627, 0.713726, 0.854902, 1.000000,
                      0.000000, 0.141176, 0.282353, 0.427451, 0.568627, 0.713726,
                      0.854902, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 0.333333, 0.443137,
                      0.443137, 0.219608, 0.776471, 0.556863, 0.556863, 0.666667,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.666667, 0.333333,
                      0.666667, 0.333333, 0.666667, 0.333333, 0.039216, 0.078431,
                      0.117647, 0.156863, 0.200000, 0.239216, 0.278431, 0.317647,
                      0.356863, 0.400000, 0.439216, 0.478431, 0.517647, 0.556863,
                      0.600000, 0.639216, 0.678431, 0.717647, 0.756863, 0.800000,
                      0.839216, 0.878431, 0.917647, 0.956863, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
                      0.000000, 0.000000, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.247059, 0.247059, 0.247059, 0.247059, 0.247059, 0.247059,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.498039, 0.498039,
                      0.498039, 0.498039, 0.498039, 0.498039, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 0.749020, 0.749020, 0.749020, 0.749020,
                      0.749020, 0.749020, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                      1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
                    };

float tkoRGBMap[8][3] = {
                          {
                            0, 0, 0
                          },
                          {
                            1, 0, 0
                          },
                          {
                            0, 1, 0
                          },
                          {
                            1, 1, 0
                          },
                          {
                            0, 0, 1
                          },
                          {
                            1, 0, 1
                          },
                          {
                            0, 1, 1
                          },
                          {
                            1, 1, 1
                          }
                        };

/******************************************************************************/

void tkoCloseWindow(void)
{

  if (xDisplay)
  {
    cursorNum = 0;

    ExposeFunc = 0;
    ReshapeFunc = 0;
    IdleFunc = 0;
    DisplayFunc = 0;
    KeyDownFunc = 0;
    MouseDownFunc = 0;
    MouseUpFunc = 0;
    MouseMoveFunc = 0;

    glFlush();
    glFinish();
    if (TKO_HAS_OVERLAY(w.type))
    {
      XDestroyWindow(xDisplay, w.wOverlay);
      glXDestroyContext(xDisplay, w.cOverlay);
      XFreeColormap(xDisplay, w.cMapOverlay);
      XFree((char *)w.vInfoOverlay);
    }
    XDestroyWindow(xDisplay, w.wMain);
    glXDestroyContext(xDisplay, w.cMain);
    XFreeColormap(xDisplay, w.cMapMain);
    XFree((char *)w.vInfoMain);
    XCloseDisplay(xDisplay);
    xDisplay = 0;
  }
}

/******************************************************************************/

static int ErrorHandler(Display *xDisplay, XErrorEvent *event)
{
  char buf[80];

  printf("\nReceived X error!\n");
  printf("\tError code   : %d\n", event->error_code);
  printf("\tRequest code : %d\n", event->request_code);
  printf("\tMinor code   : %d\n\n", event->minor_code);
  XGetErrorText(xDisplay, event->error_code, buf, 80);
  printf("\tError text : '%s'\n\n", buf);
  return 0;
}

/******************************************************************************/

void tkoInitDisplayMode(GLenum type)
{

  w.type = type;
}

/******************************************************************************/

void tkoInitDisplayModePolicy(GLenum type)
{

  w.dmPolicy = type;
}

/******************************************************************************/

GLenum tkoInitDisplayModeID(GLint id)
{
  XVisualInfo sampleVis;
  int nvis;

  if (!xDisplay && !tkoInitDisplay()) return(GL_FALSE);

  sampleVis.visualid = id;
  w.vInfoMain = XGetVisualInfo(xDisplay, VisualIDMask, &sampleVis, &nvis);

  if (w.vInfoMain)
    return GL_TRUE;
  else
    return GL_FALSE;
}

/******************************************************************************/

void tkoInitPosition(int x, int y, int width, int height)
{

  w.x = x;
  w.y = y;
  w.w = width;
  w.h = height;
}

/******************************************************************************/

static XVisualInfo *FindBestMainVisual(GLenum type)
{
  int list[32], i;

  i = 0;

  list[i++] = GLX_LEVEL;
  list[i++] = 0;

  if (TKO_IS_DOUBLE(type))
  {
    list[i++] = GLX_DOUBLEBUFFER;
  }

  if (TKO_IS_RGB(type))
  {
    list[i++] = GLX_RGBA;
    list[i++] = GLX_RED_SIZE;
    list[i++] = 1;
    list[i++] = GLX_GREEN_SIZE;
    list[i++] = 1;
    list[i++] = GLX_BLUE_SIZE;
    list[i++] = 1;
    if (TKO_HAS_ALPHA(type))
    {
      list[i++] = GLX_ALPHA_SIZE;
      list[i++] = 1;
    }
    if (TKO_HAS_ACCUM(type))
    {
      list[i++] = GLX_ACCUM_RED_SIZE;
      list[i++] = 1;
      list[i++] = GLX_ACCUM_GREEN_SIZE;
      list[i++] = 1;
      list[i++] = GLX_ACCUM_BLUE_SIZE;
      list[i++] = 1;
      if (TKO_HAS_ALPHA(type))
      {
        list[i++] = GLX_ACCUM_ALPHA_SIZE;
        list[i++] = 1;
      }
    }
  }
  else if (TKO_IS_INDEX(type))
  {
    list[i++] = GLX_BUFFER_SIZE;
    list[i++] = 1;
  }

  if (TKO_HAS_DEPTH(type))
  {
    list[i++] = GLX_DEPTH_SIZE;
    list[i++] = 1;
  }

  if (TKO_HAS_STENCIL(type))
  {
    list[i++] = GLX_STENCIL_SIZE;
    list[i++] = 1;
  }

  list[i] = (int)None;

  return glXChooseVisual(xDisplay, xScreen, list);
}

/******************************************************************************/

static XVisualInfo *FindExactMainVisual(GLenum type)
{
  int i, nvis, val, rval, gval, bval, aval;
  XVisualInfo *vis_list, *this_vis, *best_vis, sampleVis;
  int this_score, best_score;

  /* Get list of visuals for this screen */
  sampleVis.screen = xScreen;
  vis_list = XGetVisualInfo( xDisplay, VisualScreenMask, &sampleVis, &nvis);

  /*
   * Loop through the visuals; find first one that matches the attr
   * specified in type
   */
  best_score = -1;
  best_vis = NULL;
  for ( i = 0; i < nvis; i++ )
  {
    this_vis = &vis_list[i];

    /* marty */
    glXGetConfig(xDisplay, this_vis, GLX_USE_GL, &val);
    if (print_all_glx_visuals && val)
    {
      printf("vis=%d:\n",i);
      glXGetConfig(xDisplay, this_vis, GLX_LEVEL, &val);
      printf("  GLX_LEVEL (> 0 is overlay)=%d\n",val);
      glXGetConfig( xDisplay, this_vis, GLX_DOUBLEBUFFER, &val);
      printf("  GLX_DOUBLEBUFFER flag=%d\n",val);
      glXGetConfig(xDisplay, this_vis, GLX_RGBA, &val);
      printf("  GLX_RGBA flag=%d\n",val);
      glXGetConfig(xDisplay, this_vis, GLX_BUFFER_SIZE, &val);
      printf("  GLX_BUFFER_SIZE bits=%d\n",val);
      glXGetConfig(xDisplay, this_vis, GLX_DEPTH_SIZE, &val);
      printf("  GLX_DEPTH_SIZE bits=%d\n",val);
      printf("  colormapsize=%d\n",this_vis->colormap_size);
    }

    /* Visual must be supported by GLX */
    if ( glXGetConfig(xDisplay, this_vis, GLX_USE_GL, &val) ) continue;
    if ( !val ) continue;

    /* Visual must be in main planes which is level 0 */
    glXGetConfig(xDisplay, this_vis, GLX_LEVEL, &val);
    if ( val != 0 ) continue;

    /* Color Index or RGBA? It must match the requested value */
    glXGetConfig(xDisplay, this_vis, GLX_RGBA, &val);
    if ( TKO_IS_RGB(type) && !val ) continue;
    if ( TKO_IS_INDEX(type) && val ) continue;

    /* Double buffered or Single buffered? */
    glXGetConfig( xDisplay, this_vis, GLX_DOUBLEBUFFER, &val);
    if ( TKO_IS_DOUBLE(type) && !val ) continue;
    if ( TKO_IS_SINGLE(type) && val ) continue;

    /* If accum requested then accum rgb size must be > 0 */
    /* If alpha requested then alpha size must be > 0 */
    /* if accum & alpha requested then accum alpha size must be > 0 */
    if ( TKO_IS_RGB(type) )
    {
      glXGetConfig(xDisplay, this_vis, GLX_ACCUM_RED_SIZE, &rval);
      glXGetConfig(xDisplay, this_vis, GLX_ACCUM_GREEN_SIZE, &gval);
      glXGetConfig(xDisplay, this_vis, GLX_ACCUM_BLUE_SIZE, &bval);
      glXGetConfig(xDisplay, this_vis, GLX_ACCUM_ALPHA_SIZE, &aval);
      if ( TKO_HAS_ACCUM(type) )
      {
        if ( rval <= 0 || gval <= 0 || bval <= 0 ) continue;
      }
      else
      {
        if ( rval > 0 || gval > 0 || bval > 0 || aval > 0 ) continue;
      }

      glXGetConfig(xDisplay, this_vis, GLX_ALPHA_SIZE, &val);
      if ( TKO_HAS_ALPHA(type) )
      {
        if ( val <= 0 ) continue;
        if ( TKO_HAS_ACCUM(type) && aval <= 0 ) continue;
      }
      else
      {
        if ( val > 0 ) continue;
      }

    }

    /* Check depth buffer */
    glXGetConfig(xDisplay, this_vis, GLX_DEPTH_SIZE, &val);
    if ( TKO_HAS_DEPTH(type) )
    {
      if ( val <= 0 ) continue;
    }
    else
    {
      if ( val > 0 ) continue;
    }

    /* Check stencil buffer */
    glXGetConfig( xDisplay, this_vis, GLX_STENCIL_SIZE, &val);
    if ( TKO_HAS_STENCIL(type) )
    {
      if ( val <= 0 ) continue;
    }
    else
    {
      if ( val > 0 ) continue;
    }

    glXGetConfig(xDisplay, this_vis, GLX_BUFFER_SIZE, &this_score);

    if (this_score > best_score )
    {
      best_score = this_score;
      best_vis = this_vis;
    }

  }

  if ( best_vis )
  {
    sampleVis.visualid = best_vis->visualid;
    sampleVis.screen = xScreen;
    if ( nvis > 0 ) XFree((char *)vis_list);
    return XGetVisualInfo(xDisplay, VisualIDMask|VisualScreenMask,
                          &sampleVis, &nvis);
  }
  else
  {
    if ( nvis > 0 ) XFree((char *)vis_list);
    return None;
  }

}

static XVisualInfo *FindOverlayVisual(void)
{
  int list[3];

  list[0] = GLX_LEVEL;
  list[1] = 1;
  list[2] = (int)None;

  return glXChooseVisual(xDisplay, xScreen, list);
}

static GLenum GetMainWindowType(XVisualInfo *vi)
{
  GLenum mask;
  int x, y, z;

  mask = 0;

  glXGetConfig(xDisplay, vi, GLX_DOUBLEBUFFER, &x);
  if (x)
  {
    mask |= TKO_DOUBLE;
  }
  else
  {
    mask |= TKO_SINGLE;
  }

  glXGetConfig(xDisplay, vi, GLX_RGBA, &x);
  if (x)
  {
    mask |= TKO_RGB;
    glXGetConfig(xDisplay, vi, GLX_ALPHA_SIZE, &x);
    if (x > 0)
    {
      mask |= TKO_ALPHA;
    }
    glXGetConfig(xDisplay, vi, GLX_ACCUM_RED_SIZE, &x);
    glXGetConfig(xDisplay, vi, GLX_ACCUM_GREEN_SIZE, &y);
    glXGetConfig(xDisplay, vi, GLX_ACCUM_BLUE_SIZE, &z);
    if (x > 0 && y > 0 && z > 0)
    {
      mask |= TKO_ACCUM;
    }
  }
  else
  {
    mask |= TKO_INDEX;
  }

  glXGetConfig(xDisplay, vi, GLX_DEPTH_SIZE, &x);
  if (x > 0)
  {
    mask |= TKO_DEPTH;
  }

  glXGetConfig(xDisplay, vi, GLX_STENCIL_SIZE, &x);
  if (x > 0)
  {
    mask |= TKO_STENCIL;
  }

  if (glXIsDirect(xDisplay, w.cMain))
  {
    mask |= TKO_DIRECT;
  }
  else
  {
    mask |= TKO_INDIRECT;
  }

  return mask;
}

static int WaitForMainWindow(Display *d, XEvent *e, char *arg)
{

  if (e->type == MapNotify && e->xmap.window == w.wMain)
  {
    return GL_TRUE;
  }
  else
  {
    return GL_FALSE;
  }
}

#if 0
static int WaitForOverlayWindow(Display *d, XEvent *e, char *arg)
{

  if (e->type == MapNotify && e->xmap.window == w.wOverlay)
  {
    return GL_TRUE;
  }
  else
  {
    return GL_FALSE;
  }
}
#endif

GLenum tkoInitDisplay(void)
{
  int erb, evb;

  xDisplay = XOpenDisplay(0);
  if (!xDisplay)
  {
    fprintf(stderr, "Can't connect to xDisplay!\n");
    return GL_FALSE;
  }
  if (!glXQueryExtension(xDisplay, &erb, &evb))
  {
    fprintf(stderr, "No glx extension!\n");
    return GL_FALSE;
  }
  xScreen = DefaultScreen(xDisplay);
  wRoot = RootWindow(xDisplay, xScreen);
  XSetErrorHandler(ErrorHandler);
  return(GL_TRUE);
}

GLenum tkoInitPixmap(long x, long y)   /* marty */
{
  static Pixmap pixmap;
  static GLXPixmap glxpixmap;
  static GLXContext cx;
  GLenum overlayFlag;

  if (!xDisplay && !tkoInitDisplay()) return(GL_FALSE);

  if (TKO_HAS_OVERLAY(w.type))
  {
    overlayFlag = GL_TRUE;
  }
  else
  {
    overlayFlag = GL_FALSE;
  }
  w.type &= ~TKO_OVERLAY;

  if (w.dmPolicy == TKO_MINIMUM_CRITERIA)
    w.vInfoMain = FindBestMainVisual(w.type);
  else if (w.dmPolicy == TKO_EXACT_MATCH)
    w.vInfoMain = FindExactMainVisual(w.type);
  if (!w.vInfoMain)
  {
    fprintf(stderr, "Window type not found!\n");
    return GL_FALSE;
  }

  cx = glXCreateContext(xDisplay, w.vInfoMain, None, GL_FALSE);
  if (!cx)
  {
    fprintf(stderr, "could not create a drawing context!\n");
    return GL_FALSE;
  }
  pixmap = XCreatePixmap(xDisplay, RootWindow(xDisplay, xScreen),
                         x, y, w.vInfoMain->depth);
  glxpixmap = glXCreateGLXPixmap(xDisplay, w.vInfoMain, pixmap);
  glXMakeCurrent(xDisplay, glxpixmap, cx) ;
  return GL_TRUE;
}

GLenum tkoInitWindow(char *title)
{
  XSetWindowAttributes wa;
  XTextProperty tp;
  XSizeHints sh;
  XEvent e;
  GLenum overlayFlag;

  if (!xDisplay && !tkoInitDisplay()) return(GL_FALSE);

  if (TKO_HAS_OVERLAY(w.type))
  {
    overlayFlag = GL_TRUE;
  }
  else
  {
    overlayFlag = GL_FALSE;
  }
  w.type &= ~TKO_OVERLAY;

  if (w.dmPolicy == TKO_MINIMUM_CRITERIA)
    w.vInfoMain = FindBestMainVisual(w.type);
  else if (w.dmPolicy == TKO_EXACT_MATCH)
    w.vInfoMain = FindExactMainVisual(w.type);
  if (!w.vInfoMain)
  {
    fprintf(stderr, "Window type not found!\n");
    return GL_FALSE;
  }

  w.cMain = glXCreateContext(xDisplay, w.vInfoMain, None,
                             (TKO_IS_DIRECT(w.type))?GL_TRUE:GL_FALSE);
  if (!w.cMain)
  {
    fprintf(stderr, "Can't create a context!\n");
    return GL_FALSE;
  }

  w.type = GetMainWindowType(w.vInfoMain);

  if (TKO_IS_INDEX(w.type))
  {
    if (w.vInfoMain->class != StaticColor &&
        w.vInfoMain->class != StaticGray)
    {
      w.cMapMain = XCreateColormap(xDisplay, wRoot, w.vInfoMain->visual,
                                   AllocAll);
    }
    else
    {
      w.cMapMain = XCreateColormap(xDisplay, wRoot, w.vInfoMain->visual,
                                   AllocNone);
    }
  }
  else
  {
    w.cMapMain = XCreateColormap(xDisplay, wRoot, w.vInfoMain->visual,
                                 AllocNone);
  }
  tkoSetRGBMap(256, colorMaps);
  wa.colormap = w.cMapMain;
  wa.background_pixmap = None;
  wa.border_pixel = 0;
  /* marty: added KeyReleaseMask */
  wa.event_mask = StructureNotifyMask | ExposureMask | KeyPressMask |
                  KeyReleaseMask | ButtonPressMask | ButtonReleaseMask | PointerMotionMask;
  w.wMain = XCreateWindow(xDisplay, wRoot, w.x, w.y, w.w, w.h, 0,
                          w.vInfoMain->depth, InputOutput,
                          w.vInfoMain->visual,
                          CWBackPixmap|CWBorderPixel|CWEventMask|CWColormap,
                          &wa);

  XStringListToTextProperty(&title, 1, &tp);
  sh.flags = USPosition | USSize;
  XSetWMProperties(xDisplay, w.wMain, &tp, &tp, 0, 0, &sh, 0, 0);
  deleteWindowAtom = XInternAtom(xDisplay, "WM_DELETE_WINDOW", False);
  XSetWMProtocols(xDisplay, w.wMain, &deleteWindowAtom, 1);
  XMapWindow(xDisplay, w.wMain);
  drawAllowFlag = GL_FALSE;
  XIfEvent(xDisplay, &e, WaitForMainWindow, 0);

  if (overlayFlag == GL_TRUE)
  {
    w.vInfoOverlay = FindOverlayVisual();
    if (w.vInfoOverlay)
    {
      w.cOverlay = glXCreateContext(xDisplay, w.vInfoOverlay, None,
                                    GL_TRUE);
      w.cMapOverlay = XCreateColormap(xDisplay, wRoot,
                                      w.vInfoOverlay->visual, AllocNone);
      tkoSetOverlayMap(256, colorMaps);
      wa.colormap = w.cMapOverlay;
      wa.background_pixmap = None;
      wa.border_pixel = 0;
      w.wOverlay = XCreateWindow(xDisplay, w.wMain, 0, 0, w.w, w.h, 0,
                                 w.vInfoOverlay->depth, InputOutput,
                                 w.vInfoOverlay->visual,
                                 CWBackPixmap|CWBorderPixel|CWColormap,
                                 &wa);
      XMapWindow(xDisplay, w.wOverlay);
      XSetWMColormapWindows(xDisplay, w.wMain, &w.wOverlay, 1);
      w.type |= TKO_OVERLAY;
    }
    else
    {
      fprintf(stderr, "Can't create a overlay plane!\n");
    }
  }

  if (!glXMakeCurrent(xDisplay, w.wMain, w.cMain))
  {
    fprintf(stderr, "Can't make window current drawable!\n");
    return GL_FALSE;
  }
  XFlush(xDisplay);

  return GL_TRUE;
}

/******************************************************************************/

void tkoQuit(void)
{

  tkoCloseWindow();
  exit(0);
}

/******************************************************************************/

void tkoSwapBuffers(void)
{

  if (xDisplay)
  {
    glXSwapBuffers(xDisplay, w.wMain);
  }
}

/******************************************************************************/


/*############################################################################*/
/* marty: add definitions from event.c (pre-declared in header) */

#include <X11/keysym.h>

void (*ExposeFunc)(int, int) = 0;
void (*ReshapeFunc)(int, int) = 0;
void (*DisplayFunc)(void) = 0;
GLenum (*KeyDownFunc)(int, GLenum) = 0;
GLenum (*MouseDownFunc)(int, int, GLenum) = 0;
GLenum (*MouseUpFunc)(int, int, GLenum) = 0;
GLenum (*MouseMoveFunc)(int, int, GLenum) = 0;
void (*IdleFunc)(void) = 0;
int lastEventType = -1;
GLenum drawAllowFlag;


/*############################################################################*/
/* marty: add definition from cursor.c (pre-declared in header) */

int cursorNum = 0;


/*############################################################################*/
/* marty: add three functions from getset.c */

/******************************************************************************/
void tkoSetOneColor(int index, float r, float g, float b)
{
  XColor c;
  int rShift, gShift, bShift;

  switch (w.vInfoMain->class)
  {
  case DirectColor:
    rShift = ffs((unsigned int)w.vInfoMain->red_mask) - 1;
    gShift = ffs((unsigned int)w.vInfoMain->green_mask) - 1;
    bShift = ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
    c.pixel = ((index << rShift) & w.vInfoMain->red_mask) |
              ((index << gShift) & w.vInfoMain->green_mask) |
              ((index << bShift) & w.vInfoMain->blue_mask);
    c.red = (unsigned short)(r * 65535.0 + 0.5);
    c.green = (unsigned short)(g * 65535.0 + 0.5);
    c.blue = (unsigned short)(b * 65535.0 + 0.5);
    c.flags = DoRed | DoGreen | DoBlue;
    XStoreColor(xDisplay, w.cMapMain, &c);
    break;
  case GrayScale:
  case PseudoColor:
    if (index < w.vInfoMain->colormap_size)
    {
      c.pixel = index;
      c.red = (unsigned short)(r * 65535.0 + 0.5);
      c.green = (unsigned short)(g * 65535.0 + 0.5);
      c.blue = (unsigned short)(b * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, w.cMapMain, &c);
    }
    break;
  }

  XSync(xDisplay, 0);
}

/******************************************************************************/
void tkoSetOverlayMap(int size, float *rgb)
{
  XColor c;
  unsigned long *buf;
  int max, i;

  if (w.vInfoOverlay->class == PseudoColor)
  {
    max = (size > w.vInfoOverlay->colormap_size) ?
          w.vInfoOverlay->colormap_size : size;
    buf = (unsigned long *)calloc(max, sizeof(unsigned long));
    XAllocColorCells(xDisplay, w.cMapOverlay, True, NULL, 0, buf, max-1);
    for (i = 1; i < max; i++)
    {
      c.pixel = i;
      c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
      c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
      c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, w.cMapOverlay, &c);
    }
    free(buf);
  }

  XSync(xDisplay, 0);
}

/******************************************************************************/
void tkoSetRGBMap(int size, float *rgb)
{
  XColor c;
  int rShift, gShift, bShift, max, i;

  switch (w.vInfoMain->class)
  {
  case DirectColor:
    max = (size > w.vInfoMain->colormap_size) ? w.vInfoMain->colormap_size
          : size;
    for (i = 0; i < max; i++)
    {
      rShift = ffs((unsigned int)w.vInfoMain->red_mask) - 1;
      gShift = ffs((unsigned int)w.vInfoMain->green_mask) - 1;
      bShift = ffs((unsigned int)w.vInfoMain->blue_mask) - 1;
      c.pixel = ((i << rShift) & w.vInfoMain->red_mask) |
                ((i << gShift) & w.vInfoMain->green_mask) |
                ((i << bShift) & w.vInfoMain->blue_mask);
      c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
      c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
      c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, w.cMapMain, &c);
    }
    break;
  case GrayScale:
  case PseudoColor:
    max = (size > w.vInfoMain->colormap_size) ? w.vInfoMain->colormap_size
          : size;
    for (i = 0; i < max; i++)
    {
      c.pixel = i;
      c.red = (unsigned short)(rgb[i] * 65535.0 + 0.5);
      c.green = (unsigned short)(rgb[size+i] * 65535.0 + 0.5);
      c.blue = (unsigned short)(rgb[size*2+i] * 65535.0 + 0.5);
      c.flags = DoRed | DoGreen | DoBlue;
      XStoreColor(xDisplay, w.cMapMain, &c);
    }
    break;
  }

  XSync(xDisplay, 0);
}
/******************************************************************************/


#endif // ifdef HAVE_OPENGL
