/* @(#)fullscreen.h 20.28 91/09/14 	 */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_fullscreen_DEFINED
#define xview_fullscreen_DEFINED

#include <xview/win_input.h>	/* needed for SV1 compat routines */
#include <xview/rect.h>		/* needed for SV1 compat routines */
#include <xview/cursor.h>	/* needed for SV1 compat routines */
#include <xview/window.h>	/* needed for SV1 compat routines */
#include <xview/generic.h>

/*
 **********************************************************************
 * PUBLIC #defines
 **********************************************************************
 */

#ifndef XV_ATTRIBUTES_ONLY
#define	FULLSCREEN		&xv_fullscreen_pkg
#endif				/* ~XV_ATTRIBUTES_ONLY */

/*
 * PRIVATE #defines
 */

#define	FULLSCREEN_ATTR(type, ordinal)	ATTR(ATTR_PKG_FULLSCREEN, type, ordinal)

/*
 * **********************************************************************
 * Typedefs, Enumerations, and Structures 
 * **********************************************************************
 */

#ifndef XV_ATTRIBUTES_ONLY

/*
 * Public typedefs
 */

typedef Xv_opaque Fullscreen;

/*
 * Public structures
 */

/*
 * struct fullscreen is For SunView 1 fullscreen compatibility only
 */

struct fullscreen {
    int             fs_windowfd;
    Rect            fs_screenrect;
    Xv_Window       fs_rootwindow;
    Inputmask       inputmask;	/* current mask */
    Xv_Cursor       cursor;	/* current cursor */
};

typedef struct {
    Xv_generic_struct parent;
    Xv_opaque       private_data;
    Xv_embedding    embedding_data;
    struct fullscreen fullscreen_struct;
} Xv_fullscreen;

#endif				/* ~XV_ATTRIBUTES_ONLY */

typedef enum {
    FULLSCREEN_SYNCHRONOUS  = 0,
    FULLSCREEN_ASYNCHRONOUS = 1
} Fullscreen_grab_mode;

typedef enum {
    /*
     * Public Attributes
     */
    FULLSCREEN_CURSOR_WINDOW =    FULLSCREEN_ATTR(ATTR_OPAQUE, 5),    /* CG- */
    FULLSCREEN_INPUT_WINDOW =     FULLSCREEN_ATTR(ATTR_OPAQUE, 10),   /* CG- */
    FULLSCREEN_PAINT_WINDOW =     FULLSCREEN_ATTR(ATTR_OPAQUE, 15),   /* CG- */
    FULLSCREEN_RECT =             FULLSCREEN_ATTR(ATTR_NO_VALUE, 20), /* -G- */
    FULLSCREEN_SYNC =             FULLSCREEN_ATTR(ATTR_INT, 25),      /* CGS */
    FULLSCREEN_ALLOW_SYNC_EVENT = FULLSCREEN_ATTR(ATTR_NO_VALUE, 30), /* --S */
    FULLSCREEN_ALLOW_EVENTS =     FULLSCREEN_ATTR(ATTR_OPAQUE, 32),   /* --S */
    FULLSCREEN_GRAB_KEYBOARD =    FULLSCREEN_ATTR(ATTR_BOOLEAN, 35),  /* CGS */
    FULLSCREEN_GRAB_POINTER =     FULLSCREEN_ATTR(ATTR_BOOLEAN, 40),  /* CGS */
    FULLSCREEN_GRAB_SERVER =      FULLSCREEN_ATTR(ATTR_BOOLEAN, 45),  /* CGS */
    FULLSCREEN_KEYBOARD_GRAB_PTR_MODE = 
				  FULLSCREEN_ATTR(ATTR_OPAQUE, 50),   /* CGS */
    FULLSCREEN_KEYBOARD_GRAB_KBD_MODE = 
				  FULLSCREEN_ATTR(ATTR_OPAQUE, 55),   /* CGS */
    FULLSCREEN_POINTER_GRAB_PTR_MODE = 
				  FULLSCREEN_ATTR(ATTR_OPAQUE, 60),   /* CGS */
    FULLSCREEN_POINTER_GRAB_KBD_MODE = 
				  FULLSCREEN_ATTR(ATTR_OPAQUE, 65),   /* CGS */
    FULLSCREEN_OWNER_EVENTS =     FULLSCREEN_ATTR(ATTR_BOOLEAN, 60)  /* CGS */
} Fullscreen_attr;

#define FULLSCREEN_ALLOW_ASYNC_EVENTS	FULLSCREEN_SYNC, FALSE	      /* --S */

/*
 * **********************************************************************
 * 
 * 
 * Globals **********************************************************************
 * 
 */

#ifndef XV_ATTRIBUTES_ONLY

/*
 * PUBLIC variables
 */

extern Xv_pkg   xv_fullscreen_pkg;

/*
 * PUBLIC functions
 */

EXTERN_FUNCTION (struct fullscreen *fullscreen_init, (Xv_Window window));
EXTERN_FUNCTION (int	fullscreen_set_cursor, (struct fullscreen *fs, Xv_Cursor cursor));
EXTERN_FUNCTION (int	fullscreen_set_inputmask, (struct fullscreen *fs, Inputmask *im));
EXTERN_FUNCTION (int	fullscreen_destroy, (struct fullscreen *fs));

/*
 * Global variables for debugging purposes
 */
extern int	fullscreendebug;
extern int	fullscreendebugserver;
extern int	fullscreendebugptr;
extern int	fullscreendebugkbd;


#endif				/* ~XV_ATTRIBUTES_ONLY */

#endif				/* ~xview_fullscreen_DEFINED */
