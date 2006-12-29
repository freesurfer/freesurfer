/*      @(#)scrollbar.h 1.35 91/03/19 */

/*
 * (c) Copyright 1989 Sun Microsystems, Inc. Sun design patents
 * pending in the U.S. and foreign countries. See LEGAL NOTICE
 * file for terms of the license.
 */

#ifndef xview_scrollbar_DEFINED
#define xview_scrollbar_DEFINED

/*
 * Module: scrollbar.h
 * Library: libxview.a
 *
 * Level: public
 *
 * Description: Describes attributes for scrollbar
 *
 */

/*
 ***********************************************************************
 *   Include Files
 ***********************************************************************
 */

#include <xview/pkg.h>
#include <xview/window.h>
#include <xview/attrol.h>

/*
 ***********************************************************************
 *   Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines
 */

#define SCROLLBAR      &xv_scrollbar_pkg

#define SCROLLBAR_FIRST      vuid_first(SCROLL_DEVID) /* 32256 */
#define SCROLLBAR_REQUEST    (SCROLLBAR_FIRST+0)         /* 32256 */

/* SunView1 Compatiblity */
#define SCROLL_LINE_HEIGHT SCROLLBAR_PIXELS_PER_UNIT

/*
 * PRIVATE #defines
 */

#define SCROLLBAR_ATTR(type, ordinal) ATTR(ATTR_PKG_SCROLLBAR, type, ordinal)

/*
 ***********************************************************************
 *  Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

typedef Xv_opaque  Scrollbar;

typedef struct
{
  Xv_window_struct parent_data;
  Xv_opaque  private_data;
}
Xv_scrollbar;

typedef enum {
  /*
   * Public Attributes
   */
  SCROLLBAR_COMPUTE_SCROLL_PROC = SCROLLBAR_ATTR(ATTR_FUNCTION_PTR,  6),
  SCROLLBAR_DIRECTION  = SCROLLBAR_ATTR(ATTR_ENUM,   8),
  SCROLLBAR_INACTIVE  = SCROLLBAR_ATTR(ATTR_BOOLEAN,   13),
  SCROLLBAR_LAST_VIEW_START = SCROLLBAR_ATTR(ATTR_INT,  10),
  SCROLLBAR_MENU   = SCROLLBAR_ATTR(ATTR_OPAQUE,  11),
  SCROLLBAR_NORMALIZE_PROC = SCROLLBAR_ATTR(ATTR_FUNCTION_PTR,  5),
  SCROLLBAR_NOTIFY_CLIENT  = SCROLLBAR_ATTR(ATTR_OPAQUE,   9),
  SCROLLBAR_OBJECT_LENGTH  = SCROLLBAR_ATTR(ATTR_INT,   1),
  SCROLLBAR_OVERSCROLL  = SCROLLBAR_ATTR(ATTR_INT,  12),
  SCROLLBAR_PAGE_LENGTH  = SCROLLBAR_ATTR(ATTR_INT,   4),
  SCROLLBAR_PERCENT_OF_DRAG_REPAINTS = SCROLLBAR_ATTR(ATTR_INT,       14),
  SCROLLBAR_PIXELS_PER_UNIT = SCROLLBAR_ATTR(ATTR_INT,   0),
  SCROLLBAR_SPLITTABLE  = SCROLLBAR_ATTR(ATTR_BOOLEAN,   7),
  SCROLLBAR_VIEW_START  = SCROLLBAR_ATTR(ATTR_INT,   2),
  SCROLLBAR_VIEW_LENGTH  = SCROLLBAR_ATTR(ATTR_INT,   3)
} Scrollbar_attribute;


typedef enum {
  /*
   * absolute motion
   */
  SCROLLBAR_ABSOLUTE,
  /*
   * forward motions
   */
  SCROLLBAR_POINT_TO_MIN,
  SCROLLBAR_PAGE_FORWARD,
  SCROLLBAR_LINE_FORWARD,
  /*
   * backward motions
   */
  SCROLLBAR_MIN_TO_POINT,
  SCROLLBAR_PAGE_BACKWARD,
  SCROLLBAR_LINE_BACKWARD,
  /*
   * first last motions
   */
  SCROLLBAR_TO_END,
  SCROLLBAR_TO_START,
  SCROLLBAR_PAGE_ALIGNED,
  /*
   * no scrolling
   */
  SCROLLBAR_NONE
} Scroll_motion;

typedef enum {
  SCROLLBAR_VERTICAL,
  SCROLLBAR_HORIZONTAL
} Scrollbar_setting;

/*
 ***********************************************************************
 *    Globals
 ***********************************************************************
 */

extern Xv_pkg   xv_scrollbar_pkg;

/*
 * Public functions
 */

EXTERN_FUNCTION (void scrollbar_default_compute_scroll_proc, (Scrollbar sb, int pos, int length, Scroll_motion motion, unsigned long *offset, unsigned long *object_length));

EXTERN_FUNCTION (void scrollbar_paint, (Scrollbar sb));

/*
 * XView-private functions
 */
EXTERN_FUNCTION (int scrollbar_width_for_scale, (Window_rescale_state scale));

/*
 * For SunView 1 Compatibility
 */
EXTERN_FUNCTION (Scrollbar scrollbar_create, (Attr_attribute attr1, DOTDOTDOT));
EXTERN_FUNCTION (int scrollbar_set, (Scrollbar sb, DOTDOTDOT));
EXTERN_FUNCTION (Xv_opaque scrollbar_get, (Scrollbar sb, Scrollbar_attribute attr));
EXTERN_FUNCTION (int scrollbar_destroy, (Scrollbar sb));
EXTERN_FUNCTION (void scrollbar_scroll_to, (Scrollbar sb, unsigned long new_start));

#endif  /* ~xview_scrollbar_DEFINED */
