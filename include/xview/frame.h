/*	@(#)frame.h 20.67 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_frame_DEFINED
#define xview_frame_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/window.h>
#include <xview/attrol.h>
#include <X11/X.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define FRAME			FRAME_BASE
#define FRAME_BASE		&xv_frame_base_pkg
#define FRAME_CMD		&xv_frame_cmd_pkg
#define FRAME_PROPS		FRAME_CMD
#define FRAME_HELP		&xv_frame_help_pkg
#define FRAME_CLASS		&xv_frame_class_pkg

#define	ROOT_FRAME		((Frame)0)

#define FRAME_SHOW_HEADER	FRAME_SHOW_LABEL
#define FRAME_FOCUS_UP_WIDTH	13
#define FRAME_FOCUS_UP_HEIGHT	13
#define FRAME_FOCUS_RIGHT_WIDTH		13
#define FRAME_FOCUS_RIGHT_HEIGHT	13

/*
 * Utility Macros 
 */

#define frame_fit_all(frame) 					\
{ 								\
    Xv_Window win; 						\
    int n = 0; 							\
    while (win = window_get(frame, FRAME_NTH_SUBWINDOW, n++, 0))\
	window_fit(win); 					\
    window_fit(frame); 						\
}

#define frame_done_proc(frame) 					\
   (((void (*)())window_get(frame, FRAME_DONE_PROC))(frame))

#define frame_default_done_proc(frame) 				\
   (((void (*)())window_get(frame, FRAME_DEFAULT_DONE_PROC))(frame))

/*
 * PRIVATE #defines 
 */

#define	FRAME_ATTR(type, ordinal)	ATTR(ATTR_PKG_FRAME, type, ordinal)
#define FRAME_ATTR_LIST(ltype, type, ordinal) \
	FRAME_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))

/*
 * BUG: maybe these should be attributes 
 */


/*
 * width of border around a frame 
 */
#define FRAME_BORDER_WIDTH      (0)
/*
 * width of border around a subwindow 
 */
#define FRAME_SUBWINDOW_SPACING (FRAME_BORDER_WIDTH)


/*
 * PUBLIC #defines 
 * For SunView 1 Compatibility
 */

#define FRAME_TYPE		ATTR_PKG_FRAME

#define	FRAME_ARGS		XV_INIT_ARGS
#define	FRAME_ARGC_PTR_ARGV	XV_INIT_ARGC_PTR_ARGV
#define	FRAME_CMDLINE_HELP_PROC	XV_USAGE_PROC
#define	FRAME_LABEL		XV_LABEL
#ifdef OW_I18N
#define	FRAME_LABEL_WCS		XV_LABEL_WCS
#endif
#define FRAME_OPEN_RECT		WIN_RECT

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures	
 ***********************************************************************
 */

typedef	Xv_opaque 	Frame;
typedef Xv_opaque	Frame_cmd;
typedef Xv_opaque	Frame_props;
typedef Xv_opaque	Frame_help;

typedef struct {
    Xv_window_struct	parent_data;
    Xv_opaque		private_data;
} Xv_frame_class;

typedef struct {
    Xv_frame_class	parent_data;
    Xv_opaque		private_data;
} Xv_frame_base;

typedef Xv_frame_base 	Xv_frame_cmd;
typedef Xv_frame_base 	Xv_frame_props;
typedef Xv_frame_base 	Xv_frame_help;

typedef enum {
    FRAME_FOCUS_UP,
    FRAME_FOCUS_RIGHT
} Frame_focus_direction;

typedef enum {    /* do not change this order */
    FRAME_CMD_PIN_OUT,
    FRAME_CMD_PIN_IN
} Frame_cmd_pin_state;

typedef struct frame_accelerator {
    short	    code;   /* = event->ie_code */
    KeySym	    keysym; /* X keysym */
    void	  (*notify_proc)(); /* accelerator notify proc */
    Xv_opaque	    data;   /* opaque data handle */
    struct frame_accelerator *next;
} Frame_accelerator;
 
typedef enum {
	/*
	 * PUBLIC attributes 
	 */
	FRAME_BACKGROUND_COLOR	= FRAME_ATTR(ATTR_SINGLE_COLOR_PTR,	   5),
	FRAME_BUSY		= FRAME_ATTR(ATTR_BOOLEAN,		  10),
	FRAME_CLOSED	= FRAME_ATTR(ATTR_BOOLEAN,		  15),
	FRAME_CLOSED_RECT	= FRAME_ATTR(ATTR_RECT_PTR,		  20),
	FRAME_WM_COMMAND_STRINGS= FRAME_ATTR_LIST(ATTR_NULL, ATTR_STRING, 21),
	FRAME_WM_COMMAND_ARGC_ARGV
				= FRAME_ATTR(ATTR_INT_PAIR, 		  22),
	FRAME_WM_COMMAND_ARGV	= FRAME_ATTR(ATTR_OPAQUE,		  23),
	FRAME_WM_COMMAND_ARGC	= FRAME_ATTR(ATTR_INT,			  24),
	FRAME_CMD_PANEL		= FRAME_ATTR(ATTR_OPAQUE,		  25),
	FRAME_CURRENT_RECT	= FRAME_ATTR(ATTR_RECT_PTR,		  35),
	FRAME_OLD_RECT          = FRAME_ATTR(ATTR_RECT_PTR,               36),
	FRAME_DEFAULT_DONE_PROC	= FRAME_ATTR(ATTR_FUNCTION_PTR,		  40),
	FRAME_DONE_PROC		= FRAME_ATTR(ATTR_FUNCTION_PTR,		  45),
	FRAME_FOCUS_WIN		= FRAME_ATTR(ATTR_INT_PAIR,		 165),
	FRAME_FOCUS_DIRECTION	= FRAME_ATTR(ATTR_ENUM,			 170),
	FRAME_FOREGROUND_COLOR	= FRAME_ATTR(ATTR_SINGLE_COLOR_PTR,	  50),
	FRAME_ICON		= FRAME_ATTR(ATTR_OPAQUE,		  55),
	FRAME_INHERIT_COLORS	= FRAME_ATTR(ATTR_BOOLEAN,		  60),
	FRAME_LEFT_FOOTER	= FRAME_ATTR(ATTR_STRING,		  65),
#ifdef OW_I18N
	FRAME_LEFT_FOOTER_WCS	= FRAME_ATTR(ATTR_WSTRING,		  66),
#endif
	FRAME_NEXT_PANE		= FRAME_ATTR(ATTR_NO_VALUE,		  67),
	FRAME_NO_CONFIRM	= FRAME_ATTR(ATTR_BOOLEAN,		  70),
	FRAME_NTH_SUBFRAME	= FRAME_ATTR(ATTR_INT,			  75),
	FRAME_NTH_SUBWINDOW	= FRAME_ATTR(ATTR_INT,			  80),
	FRAME_PREVIOUS_ELEMENT	= FRAME_ATTR(ATTR_NO_VALUE,		  81),
	FRAME_PREVIOUS_PANE	= FRAME_ATTR(ATTR_NO_VALUE,		  82),
	FRAME_PROPERTIES_PROC	= FRAME_ATTR(ATTR_FUNCTION_PTR,		  85),
	FRAME_CMD_PUSHPIN_IN	= FRAME_ATTR(ATTR_BOOLEAN,		 105),
	FRAME_CMD_DEFAULT_PIN_STATE = FRAME_ATTR(ATTR_ENUM,		 106),
	FRAME_CMD_PIN_STATE	= FRAME_ATTR(ATTR_ENUM,		 	 107),
	FRAME_RIGHT_FOOTER	= FRAME_ATTR(ATTR_STRING,		 115),
#ifdef OW_I18N
	FRAME_RIGHT_FOOTER_WCS	= FRAME_ATTR(ATTR_WSTRING,		 116),
#endif
	FRAME_SHOW_FOOTER	= FRAME_ATTR(ATTR_BOOLEAN,		 125),
	FRAME_SHOW_RESIZE_CORNER = FRAME_ATTR(ATTR_BOOLEAN,		 128),
	FRAME_SHOW_LABEL	= FRAME_ATTR(ATTR_BOOLEAN,		 130),
	FRAME_GROUP_LEADER	= FRAME_ATTR(ATTR_BOOLEAN,		 135),
	FRAME_MIN_SIZE		= FRAME_ATTR(ATTR_INT_PAIR,	 	 136),
	FRAME_MAX_SIZE		= FRAME_ATTR(ATTR_INT_PAIR,	 	 137),
	/*
	 * PRIVATE attributes 
	 */
#ifdef	OW_I18N
	FRAME_INPUT_WINDOW	= FRAME_ATTR(ATTR_OPAQUE,		 139),
#endif	/* OW_I18N */
	FRAME_NEXT_CHILD	= FRAME_ATTR(ATTR_OPAQUE,		 140),
	FRAME_PREVIOUS_CHILD	= FRAME_ATTR(ATTR_OPAQUE,		 142),
	FRAME_SCALE_STATE	= FRAME_ATTR(ATTR_INT,			 145),
	FRAME_SUBWINDOWS_ADJUSTABLE	
				= FRAME_ATTR(ATTR_BOOLEAN,		 150),
        FRAME_COUNT             = FRAME_ATTR(ATTR_INT,                   160),
	FRAME_FOCUS_UP_IMAGE	= FRAME_ATTR(ATTR_OPAQUE,		 175),
	FRAME_FOCUS_RIGHT_IMAGE	= FRAME_ATTR(ATTR_OPAQUE,		 180),
	FRAME_FOCUS_GC		= FRAME_ATTR(ATTR_OPAQUE,		 185),
	FRAME_ORPHAN_WINDOW	= FRAME_ATTR(ATTR_INT,			 190),
	FRAME_FOOTER_WINDOW	= FRAME_ATTR(ATTR_BOOLEAN,               195),
#ifdef OW_I18N
	FRAME_IMSTATUS_WINDOW	= FRAME_ATTR(ATTR_BOOLEAN,               196),
#endif	
	FRAME_ACCELERATOR	= FRAME_ATTR(ATTR_INT_TRIPLE,		 200),
	FRAME_X_ACCELERATOR	= FRAME_ATTR(ATTR_INT_TRIPLE,		 205),
#ifdef OW_I18N
	FRAME_LEFT_IMSTATUS_WCS	= FRAME_ATTR(ATTR_WSTRING,		 210),
	FRAME_LEFT_IMSTATUS     = FRAME_ATTR(ATTR_STRING, 		 215),
	FRAME_RIGHT_IMSTATUS_WCS= FRAME_ATTR(ATTR_WSTRING,		 220),
        FRAME_RIGHT_IMSTATUS    = FRAME_ATTR(ATTR_STRING, 		 225),
	FRAME_SHOW_IMSTATUS	= FRAME_ATTR(ATTR_BOOLEAN,               230),
	FRAME_CMD_POINTER_WARP	= FRAME_ATTR(ATTR_BOOLEAN,		 240),
#endif
	FRAME_COMPOSE_STATE	= FRAME_ATTR(ATTR_BOOLEAN,               235)
	
} Frame_attribute;

#define	FRAME_PROPS_PUSHPIN_IN	FRAME_CMD_PUSHPIN_IN
#define	FRAME_PROPS_PANEL	FRAME_CMD_PANEL

/*
 *  values for scale state
 */
#define Frame_rescale_state	Window_rescale_state
#define FRAME_SCALE_SMALL	WIN_SCALE_SMALL
#define FRAME_SCALE_MEDIUM	WIN_SCALE_MEDIUM
#define FRAME_SCALE_LARGE	WIN_SCALE_LARGE
#define FRAME_SCALE_EXTRALARGE	WIN_SCALE_EXTRALARGE

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

extern Xv_pkg	xv_frame_class_pkg;
extern Xv_pkg	xv_frame_base_pkg;
extern Xv_pkg	xv_frame_cmd_pkg;
extern Xv_pkg	xv_frame_props_pkg;
extern Xv_pkg	xv_frame_help_pkg;

/*
 * XView Private functions
 */
EXTERN_FUNCTION (void frame_cmdline_help, (char *name));
EXTERN_FUNCTION (void frame_grant_extend_to_edge, (Frame frame, int to_right));
EXTERN_FUNCTION (void frame_kbd_use, (Frame frame, Xv_Window sw, Xv_Window pw));
EXTERN_FUNCTION (void frame_kbd_done, (Frame frame, Xv_Window sw));

#endif /* xview_frame_DEFINED */

