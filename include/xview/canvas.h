/*      @(#)canvas.h 20.35 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_canvas_DEFINED
#define xview_canvas_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */
#include <xview/openwin.h>
#include <xview/pixwin.h>
#include <xview/win_input.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define CANVAS			&xv_canvas_pkg
#define CANVAS_VIEW		&xv_canvas_view_pkg
#define CANVAS_PAINT_WINDOW	&xv_canvas_pw_pkg

#define	CANVAS_PIXWIN		CANVAS_NTH_PAINT_WINDOW, 0
#define	CANVAS_AUTO_CLEAR	OPENWIN_AUTO_CLEAR

/*
 * Some useful macros
 */
#define	canvas_pixwin(canvas)	\
	((Pixwin *)xv_get(canvas, CANVAS_NTH_PAINT_WINDOW, 0))
#define	canvas_paint_window(canvas) \
	((Xv_Window)xv_get(canvas, CANVAS_NTH_PAINT_WINDOW, 0))

#define	CANVAS_EACH_PAINT_WINDOW(canvas, pw)	\
   {int _pw_cnt = 0; \
   while (((pw) = (Xv_Window) xv_get((canvas), CANVAS_NTH_PAINT_WINDOW, _pw_cnt++)) != NULL) { \

#define	CANVAS_END_EACH	}}

/*
 * PRIVATE #defines 
 */

#define	CANVAS_ATTR(type, ordinal)	ATTR(ATTR_PKG_CANVAS, type, ordinal)
#define	CANVAS_VIEW_ATTR(type, ordinal)	ATTR(ATTR_PKG_CANVAS_VIEW,type, ordinal)
#define	CANVAS_PAINT_ATTR(type, ordinal) ATTR(ATTR_PKG_CANVAS_PAINT_WINDOW,type, ordinal)
#define CANVAS_ATTR_LIST(ltype, type, ordinal) \
    CANVAS_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))


/*
 * For SunView 1 compatibility 
 */
#define	CANVAS_TYPE	ATTR_PKG_CANVAS
#define CANVAS_MARGIN	CANVAS_VIEW_MARGIN

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structs
 ***********************************************************************
 */

typedef Xv_opaque	Canvas;
typedef Xv_opaque	Canvas_view;
typedef Xv_opaque	Canvas_paint_window;

/*
 * Enumerations 
 */

typedef enum {
	/*
 	 * Public attributes 
	 */
    CANVAS_AUTO_EXPAND		= CANVAS_ATTR(ATTR_BOOLEAN, 	  1),
    CANVAS_AUTO_SHRINK		= CANVAS_ATTR(ATTR_BOOLEAN, 	  5),
    CANVAS_FIXED_IMAGE		= CANVAS_ATTR(ATTR_BOOLEAN, 	 10),
    CANVAS_HEIGHT		= CANVAS_ATTR(ATTR_Y, 		 15),
    CANVAS_MIN_PAINT_HEIGHT 	= CANVAS_ATTR(ATTR_Y, 		 20),
    CANVAS_MIN_PAINT_WIDTH  	= CANVAS_ATTR(ATTR_X, 		 25),
    CANVAS_NTH_PAINT_WINDOW	= CANVAS_ATTR(ATTR_OPAQUE,  	 30),
    CANVAS_REPAINT_PROC		= CANVAS_ATTR(ATTR_FUNCTION_PTR, 35),
    CANVAS_RESIZE_PROC		= CANVAS_ATTR(ATTR_FUNCTION_PTR, 40),
    CANVAS_RETAINED		= CANVAS_ATTR(ATTR_BOOLEAN, 	 45),
    CANVAS_VIEW_MARGIN		= CANVAS_ATTR(ATTR_INT, 	 50),
    CANVAS_VIEWABLE_RECT 	= CANVAS_ATTR(ATTR_RECT_PTR, 	 55),
    CANVAS_WIDTH		= CANVAS_ATTR(ATTR_X, 		 60),
    CANVAS_X_PAINT_WINDOW	= CANVAS_ATTR(ATTR_BOOLEAN,      65),
    CANVAS_PAINTWINDOW_ATTRS	= CANVAS_ATTR_LIST(ATTR_RECURSIVE, ATTR_AV, 70),
    CANVAS_NO_CLIPPING		= CANVAS_ATTR(ATTR_BOOLEAN,      75),
    CANVAS_CMS_REPAINT          = CANVAS_ATTR(ATTR_BOOLEAN,      80)
} Canvas_attribute;


typedef enum {
    CANVAS_VIEW_PAINT_WINDOW = CANVAS_VIEW_ATTR(ATTR_OPAQUE, 1),
    CANVAS_VIEW_CANVAS_WINDOW = CANVAS_VIEW_ATTR(ATTR_OPAQUE, 2)
} Canvas_view_attribute;

typedef enum {
    CANVAS_PAINT_CANVAS_WINDOW = CANVAS_PAINT_ATTR(ATTR_OPAQUE, 1),
    CANVAS_PAINT_VIEW_WINDOW   = CANVAS_PAINT_ATTR(ATTR_OPAQUE, 2)
} Canvas_paint_attribute;

/*
 * Structures 
 */

typedef struct {
    Xv_openwin	parent_data;
    Xv_opaque	private_data;
} Xv_canvas;

typedef struct {
    Xv_window_struct	parent_data;
    Xv_opaque		private_data;
} Xv_canvas_view;

typedef struct {
    Xv_window_struct	parent_data;
    Xv_opaque		private_data;
} Xv_canvas_pw;


/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

/*
 * Variables 
 */
extern Xv_pkg	xv_canvas_pkg;
extern Xv_pkg	xv_canvas_view_pkg;
extern Xv_pkg	xv_canvas_pw_pkg;

/*
 * Functions
 */
EXTERN_FUNCTION (Event * canvas_event, (Canvas canvas_obj, Event *event));
EXTERN_FUNCTION (Event * canvas_window_event, (Canvas canvas_obj, Event *event));

#endif /* ~xview_canvas_DEFINED */
