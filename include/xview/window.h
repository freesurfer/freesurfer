/*      @(#)window.h 20.88 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_window_DEFINED
#define xview_window_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/generic.h>
#include <xview/server.h>
#include <xview/screen.h>
#include <xview/drawable.h>
#include <xview/win_input.h>
#include <xview/rect.h>
#include <X11/Xlib.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define	WINDOW				&xv_window_pkg

#define WIN_EXTEND_TO_EDGE 		-1
#define WIN_DEFAULT_BORDER_WIDTH	1
#define WIN_MESSAGE_DATA_SIZE	        20

#define window_main_loop	xv_main_loop

typedef enum {                  /* values for WIN_SCALE_STATE */
	WIN_SCALE_SMALL,
	WIN_SCALE_MEDIUM,
	WIN_SCALE_LARGE,
	WIN_SCALE_EXTRALARGE
} Window_rescale_state;

#ifndef ForgetGravity
/* Bit Gravity */
 
#define ForgetGravity           0
#define NorthWestGravity        1
#define NorthGravity            2
#define NorthEastGravity        3
#define WestGravity             4
#define CenterGravity           5
#define EastGravity             6
#define SouthWestGravity        7
#define SouthGravity            8
#define SouthEastGravity        9
#define StaticGravity           10

 /* Window gravity + bit gravity above */

#define UnmapGravity            0
#endif


/*
 * Useful conversion macros
 */
#define XV_SCREEN_FROM_WINDOW(window) \
	(Xv_Screen) xv_get(window, XV_SCREEN)

#define XV_SERVER_FROM_WINDOW(window) \
	(Xv_Server) xv_get(xv_get(window, XV_SCREEN), SCREEN_SERVER)

#define XV_DISPLAY_FROM_WINDOW(window) \
	(Display *) xv_get(xv_get(xv_get(window, XV_SCREEN), \
		SCREEN_SERVER), XV_DISPLAY)

/*
 * Window-fitting macros 
 */
#define window_fit(win) \
   (void)window_set(win, WIN_FIT_HEIGHT, 0, WIN_FIT_WIDTH, 0, 0)

#define window_fit_height(win) \
   (void)window_set(win, WIN_FIT_HEIGHT, 0, 0)

#define window_fit_width(win) \
   (void)window_set(win, WIN_FIT_WIDTH, 0, 0)

/*
 * PRIVATE #defines
 */

/*
 * Attribute macros 
 */
#define	WIN_ATTR(type, ordinal)		ATTR(ATTR_PKG_WIN, type, ordinal)
#define WIN_ATTR_LIST(ltype, type, ordinal) \
		WIN_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))
/*
 * Fake types -- This should be resolved someday 
 */
#define ATTR_IMASK			ATTR_OPAQUE

/*
 * WIN_RECT_INFO flags for package implementors 
 */
#define	WIN_X_SET			0x1
#define	WIN_Y_SET			0x2
#define	WIN_WIDTH_SET			0x4
#define	WIN_HEIGHT_SET			0x8

/*
 * PUBLIC #defines
 * For SunView 1 Compatibility Only 
 */

#define WIN_X				XV_X
#define WIN_Y				XV_Y
#define WIN_WIDTH			XV_WIDTH
#define WIN_HEIGHT			XV_HEIGHT
#define WIN_FONT			XV_FONT
#define WIN_DEVICE_NAME			XV_XNAME
#define WIN_DEVICE_NUMBER		XV_DEVICE_NUMBER
#define WIN_TOP_MARGIN			XV_TOP_MARGIN
#define WIN_BOTTOM_MARGIN		XV_BOTTOM_MARGIN
#define WIN_LEFT_MARGIN			XV_LEFT_MARGIN
#define WIN_RIGHT_MARGIN		XV_RIGHT_MARGIN
#define WIN_NAME			XV_NAME
#define WIN_OWNER			XV_OWNER
#define WIN_FD				XV_SELF
#define WIN_PIXWIN			XV_SELF
#define WIN_RECT			XV_RECT

#define WINDOW_TYPE			ATTR_PKG_WIN

#define	WIN_CONSUME_KBD_EVENT		WIN_CONSUME_EVENT
#define	WIN_IGNORE_KBD_EVENT		WIN_IGNORE_EVENT
#define	WIN_CONSUME_KBD_EVENTS		WIN_CONSUME_EVENTS	
#define	WIN_IGNORE_KBD_EVENTS		WIN_IGNORE_EVENTS
#define	WIN_PICK_INPUT_MASK		WIN_INPUT_MASK
#define	WIN_CONSUME_PICK_EVENT		WIN_CONSUME_EVENT
#define	WIN_IGNORE_PICK_EVENT		WIN_IGNORE_EVENT
#define	WIN_CONSUME_PICK_EVENTS		WIN_CONSUME_EVENTS
#define	WIN_IGNORE_PICK_EVENTS		WIN_IGNORE_EVENTS

#define WIN_NOTIFY_EVENT_PROC		WIN_NOTIFY_SAFE_EVENT_PROC

#ifdef OW_I18N
#define WIN_IM_PREEDIT_START 		WIN_IC_PREEDIT_START
#define	WIN_IM_PREEDIT_DRAW  		WIN_IC_PREEDIT_DRAW
#define	WIN_IM_PREEDIT_DONE  		WIN_IC_PREEDIT_DONE
#define	WIN_IM_STATUS_START  		WIN_IC_STATUS_START
#define	WIN_IM_STATUS_DRAW   		WIN_IC_STATUS_DRAW
#define	WIN_IM_STATUS_DONE   		WIN_IC_STATUS_DONE
#define	WIN_IM_LUC_START     		WIN_IC_LUC_START
#define	WIN_IM_LUC_DRAW      		WIN_IC_LUC_DRAW
#define	WIN_IM_LUC_DONE      		WIN_IC_LUC_DONE
#define	WIN_IM_LUC_PROCESS   		WIN_IC_LUC_PROCESS
#endif /* OW_I18N */


/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

typedef	Xv_opaque 		Xv_Window;
typedef	Xv_opaque 		Xv_window;
typedef	Attr_pkg	 	Window_type;

typedef enum {
	/*
	 * Public Attributes 
	 */
	WIN_ALARM		= WIN_ATTR(ATTR_NO_VALUE, 	  1),	/* --S */
	WIN_BACK		= WIN_ATTR(ATTR_NO_VALUE,	  2),   /* -S- */
	WIN_BACKGROUND_PIXMAP	= WIN_ATTR(ATTR_OPAQUE,		  3),
	WIN_BELOW		= WIN_ATTR(ATTR_OPAQUE,		  4),	/* --S */
	WIN_CLIENT_DATA		= WIN_ATTR(ATTR_OPAQUE,		  8),
	WIN_COLUMNS		= WIN_ATTR(ATTR_INT,		 12),
	WIN_COLUMN_GAP		= WIN_ATTR(ATTR_INT,		 16),
	WIN_COLUMN_WIDTH	= WIN_ATTR(ATTR_INT,		 20),
#ifdef OW_I18N
        WIN_IC_PREEDIT_START    = WIN_ATTR(ATTR_INT_PAIR,        21),   /* C-S */
        WIN_IC_PREEDIT_DRAW     = WIN_ATTR(ATTR_INT_PAIR,        22),   /* C-S */
        WIN_IC_PREEDIT_DONE     = WIN_ATTR(ATTR_INT_PAIR,        23),   /* C-S */
#endif /* OW_I18N */
	WIN_CONSUME_EVENT	= WIN_ATTR(ATTR_ENUM,		 24),	/* --S */
#ifdef OW_I18N
        WIN_IC_STATUS_START     = WIN_ATTR(ATTR_INT_PAIR,        25),   /* C-S */
        WIN_IC_STATUS_DRAW      = WIN_ATTR(ATTR_INT_PAIR,        26),   /* C-S */
        WIN_IC_STATUS_DONE      = WIN_ATTR(ATTR_INT_PAIR,        27),   /* C-S */
#endif /* OW_I18N */
	WIN_CONSUME_EVENTS	
			= WIN_ATTR_LIST(ATTR_NULL, ATTR_ENUM,	 28),	/* --S */
#ifdef OW_I18N
        WIN_IC_LUC_START        = WIN_ATTR(ATTR_INT_PAIR,        29),   /* C-S */
        WIN_IC_LUC_DRAW         = WIN_ATTR(ATTR_INT_PAIR,        30),   /* C-S */
        WIN_IC_LUC_DONE         = WIN_ATTR(ATTR_INT_PAIR,        31),   /* C-S */
#endif /* OW_I18N */
	WIN_CURSOR		= WIN_ATTR(ATTR_CURSOR_PTR,	 32),
#ifdef OW_I18N
        WIN_IC_LUC_PROCESS      = WIN_ATTR(ATTR_INT_PAIR,        33),
#endif /* OW_I18N */

	WIN_DEPTH		= WIN_ATTR(ATTR_INT,		102),	/* CG- */
	/* WIN_DEPTH_V2 is needed to keep backwards compatibility with
	 * V2.  WIN_DEPTH was mistakenly defined as and ATTR_NO_VALUE,
	 * which meant that despite the documentation you couldn't
	 * use it in the creation of a window.  Thus we created
	 * a new WIN_DEPTH that is of the correct type, and moved
	 * the old WIN_DEPTH to WIN_DEPTH_V2.  This can be removed
	 * when it has sufficiently aged.
	 */
	WIN_DEPTH_V2		= WIN_ATTR(ATTR_NO_VALUE,	 36),	/* -G- */

	WIN_DESIRED_HEIGHT	= WIN_ATTR(ATTR_INT, 		 40),
	WIN_DESIRED_WIDTH	= WIN_ATTR(ATTR_INT, 		 44),
	WIN_ERROR_MSG		= WIN_ATTR(ATTR_STRING,		 48),
#ifdef  OW_I18N
        WIN_ERROR_MSG_WCS       = WIN_ATTR(ATTR_STRING,          49),
#endif  /* OW_I18N */
	WIN_EVENT_PROC		= WIN_ATTR(ATTR_FUNCTION_PTR,	 52),
	WIN_FIT_HEIGHT		= WIN_ATTR(ATTR_Y,		 60),
	WIN_FIT_WIDTH		= WIN_ATTR(ATTR_X,		 64),
	WIN_FRONT		= WIN_ATTR(ATTR_NO_VALUE,	 66),   /* -S- */
	WIN_GLYPH_FONT		= WIN_ATTR(ATTR_OPAQUE,		 67),
	WIN_GRAB_ALL_INPUT	= WIN_ATTR(ATTR_BOOLEAN,	 68),	/* --S */
	WIN_HORIZONTAL_SCROLLBAR= WIN_ATTR(ATTR_OPAQUE,		 72),
#ifdef OW_I18N
        WIN_IC                  = WIN_ATTR(ATTR_OPAQUE,          74),   /* -GS */
        WIN_USE_IM              = WIN_ATTR(ATTR_BOOLEAN,         75),   /* CG- */
#endif /* OW_I18N */
	WIN_IGNORE_EVENT	= WIN_ATTR(ATTR_ENUM,		 76),	/* --S */
	WIN_IGNORE_EVENTS
			= WIN_ATTR_LIST(ATTR_NULL, ATTR_ENUM,	 80),	/* --S */
#ifdef OW_I18N
	WIN_IC_COMMIT		= WIN_ATTR(ATTR_STRING,		 81),  /* --G */
	WIN_IC_COMMIT_WCS	= WIN_ATTR(ATTR_STRING,		 82),  /* --G */
	WIN_IC_CONVERSION	= WIN_ATTR(ATTR_BOOLEAN,	 83),  /* -SG */
#endif /* OW_I18N */
	WIN_INPUT_MASK		= WIN_ATTR(ATTR_IMASK,		 84),
#ifdef OW_I18N
	WIN_IC_RESET		= WIN_ATTR(ATTR_NO_VALUE,	 85),  /* -S- */
#endif /* OW_I18N */
	WIN_IS_CLIENT_PANE	= WIN_ATTR(ATTR_NO_VALUE,	 88),	
	WIN_MENU		= WIN_ATTR(ATTR_OPAQUE,		 92),
	WIN_MOUSE_XY		= WIN_ATTR(ATTR_XY,		 96),	/* -GS */
	WIN_PARENT		= WIN_ATTR(ATTR_OPAQUE,		100),	/* -GS */
	WIN_PERCENT_HEIGHT	= WIN_ATTR(ATTR_INT,		104),
	WIN_PERCENT_WIDTH	= WIN_ATTR(ATTR_INT, 		108),
	WIN_RIGHT_OF		= WIN_ATTR(ATTR_OPAQUE,		116),	/* --S */
	WIN_ROWS		= WIN_ATTR(ATTR_INT,		120),
	WIN_ROW_GAP		= WIN_ATTR(ATTR_INT,		124),
	WIN_ROW_HEIGHT		= WIN_ATTR(ATTR_INT,		128),
	WIN_SCREEN_RECT		= WIN_ATTR(ATTR_NO_VALUE,	132),	/* -G- */
	WIN_SET_FOCUS		= WIN_ATTR(ATTR_NO_VALUE,	228),	/* --S */
	WIN_VERTICAL_SCROLLBAR	= WIN_ATTR(ATTR_OPAQUE,		140),
	WIN_MESSAGE_TYPE	= WIN_ATTR(ATTR_OPAQUE, 	141),
	WIN_MESSAGE_FORMAT	= WIN_ATTR(ATTR_INT, 		142),
	WIN_MESSAGE_DATA	= WIN_ATTR(ATTR_OPAQUE,		143),
	/*
	 * Private Attributes 
	 */
	WIN_ALARM_DATA		= WIN_ATTR(ATTR_OPAQUE, 	144),	/* -G- */
	WIN_BORDER		= WIN_ATTR(ATTR_BOOLEAN, 	148),
	WIN_FINDINTERSECT	= WIN_ATTR(ATTR_XY, 		152),	/* -G- */
	WIN_FRAME		= WIN_ATTR(ATTR_OPAQUE, 	156),
	WIN_INPUT_ONLY		= WIN_ATTR(ATTR_NO_VALUE,	160),
	WIN_IS_IN_FULLSCREEN_MODE = WIN_ATTR(ATTR_INT, 		164),
	WIN_IS_ROOT		= WIN_ATTR(ATTR_NO_VALUE,	168),
	WIN_KBD_FOCUS		= WIN_ATTR(ATTR_BOOLEAN, 	172),
	WIN_LAYOUT_PROC		= WIN_ATTR(ATTR_FUNCTION_PTR,	176),
	WIN_LOCKDATA		= WIN_ATTR(ATTR_NO_VALUE, 	180), 	/* --S */
	WIN_MAP			= WIN_ATTR(ATTR_BOOLEAN,	184),
	WIN_NOTIFY_SAFE_EVENT_PROC 	= WIN_ATTR(ATTR_FUNCTION_PTR,	192),
	WIN_NOTIFY_IMMEDIATE_EVENT_PROC	= WIN_ATTR(ATTR_FUNCTION_PTR,	193),
	WIN_OUTER_RECT		= WIN_ATTR(ATTR_RECT_PTR, 	200),	
	WIN_RECT_INFO		= WIN_ATTR(ATTR_INT, 		204),
	WIN_RETAINED		= WIN_ATTR(ATTR_BOOLEAN, 	208),
	WIN_SELECTION_WINDOW	= WIN_ATTR(ATTR_NO_VALUE, 	212),	/* --S */
	WIN_TOP_LEVEL		= WIN_ATTR(ATTR_BOOLEAN, 	216),
	WIN_TOP_LEVEL_NO_DECOR	= WIN_ATTR(ATTR_BOOLEAN, 	220),
	WIN_TRANSPARENT		= WIN_ATTR(ATTR_NO_VALUE,	223),	/* C-- */
	WIN_SAVE_UNDER		= WIN_ATTR(ATTR_BOOLEAN, 	226),
	WIN_REMOVE_CARET	= WIN_ATTR(ATTR_NO_VALUE,	227),   /* --S */
	WIN_X_PAINT_WINDOW	= WIN_ATTR(ATTR_BOOLEAN,        229),
	WIN_INHERIT_COLORS	= WIN_ATTR(ATTR_BOOLEAN,        230),
	WIN_CMS			= WIN_ATTR(ATTR_OPAQUE,		231),
        WIN_DYNAMIC_VISUAL      = WIN_ATTR(ATTR_BOOLEAN,        232),
	WIN_CMS_CHANGE		= WIN_ATTR(ATTR_NO_VALUE,	241),
	WIN_COLOR_INFO		= WIN_ATTR(ATTR_OPAQUE, 	242),
	WIN_CMD_LINE		= WIN_ATTR(ATTR_STRING, 	244),
	WIN_NO_CLIPPING		= WIN_ATTR(ATTR_BOOLEAN,	245),
	WIN_ADD_DROP_INTEREST   = WIN_ATTR(ATTR_OPAQUE,         246),
	WIN_DELETE_DROP_INTEREST  = WIN_ATTR(ATTR_OPAQUE,       247),
	WIN_ADD_DROP_ITEM       = WIN_ATTR(ATTR_OPAQUE,         252),
	WIN_DELETE_DROP_ITEM    = WIN_ATTR(ATTR_LONG,           253),

	/*
         * Public Attributes
         */
	WIN_X_CLIP_RECTS	= WIN_ATTR(ATTR_OPAQUE,		233),  /* -G- */
	WIN_CMS_DATA		= WIN_ATTR(ATTR_OPAQUE,		235),
	WIN_CMS_NAME		= WIN_ATTR(ATTR_STRING,		236),
	WIN_BIT_GRAVITY         = WIN_ATTR(ATTR_INT,            237),
	WIN_WINDOW_GRAVITY      = WIN_ATTR(ATTR_INT,            238),
	WIN_FOREGROUND_COLOR	= WIN_ATTR(ATTR_INT,            239),
	WIN_BACKGROUND_COLOR	= WIN_ATTR(ATTR_INT,            240),
	WIN_X_COLOR_INDICES	= WIN_ATTR(ATTR_OPAQUE,         243),
	WIN_CONSUME_X_EVENT_MASK = WIN_ATTR(ATTR_LONG,		248),
	WIN_IGNORE_X_EVENT_MASK = WIN_ATTR(ATTR_LONG,		249),
	WIN_X_EVENT_MASK 	= WIN_ATTR(ATTR_LONG,		250),
	WIN_COLLAPSE_EXPOSURES  = WIN_ATTR(ATTR_BOOLEAN,	251),
	WIN_SOFT_FNKEY_LABELS   = WIN_ATTR(ATTR_STRING, 	203),

	/*
	 * Public Attributes 
	 * For SunView 1 Compatibility
	 */
	WIN_TYPE		= WIN_ATTR(ATTR_ENUM,		224)
} Window_attribute;
#define WIN_SHOW	XV_SHOW

typedef enum {
	WIN_NULL_VALUE = 0,
	WIN_NO_EVENTS,
	WIN_UP_EVENTS,
	WIN_ASCII_EVENTS,
	WIN_UP_ASCII_EVENTS,
	WIN_MOUSE_BUTTONS,
	WIN_IN_TRANSIT_EVENTS,
	WIN_LEFT_KEYS,
	WIN_TOP_KEYS,
	WIN_RIGHT_KEYS,
	WIN_META_EVENTS,
	WIN_UP_META_EVENTS,
	/*
 	 * semantic event classes 
 	 */
	WIN_SUNVIEW_FUNCTION_KEYS,
	WIN_EDIT_KEYS,
	WIN_MOTION_KEYS,
	WIN_TEXT_KEYS
} Window_input_event;

typedef enum {
	WIN_CREATE, 
	WIN_INSERT,
	WIN_REMOVE,
	WIN_DESTROY,
	WIN_GET_RIGHT_OF, 
	WIN_GET_BELOW, 
	WIN_ADJUST_RECT, 
	WIN_GET_X, 
	WIN_GET_Y, 
	WIN_GET_WIDTH, 
	WIN_GET_HEIGHT,
	WIN_GET_RECT, 
	WIN_LAYOUT,
	WIN_INSTALL
} Window_layout_op;

typedef struct {
	Xv_drawable_struct	parent_data;
	Xv_opaque		private_data;
} Xv_window_struct;

typedef struct window_rescale_rect_obj {
    Rect        old_rect;
    Rect        new_rect;
    int         width_change, height_change,x_change,y_change;
    int         adjusted;
    Xv_Window   sw;
/* relationships */
} Window_rescale_rect_obj;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

extern Xv_pkg		xv_window_pkg;

/*
 * PUBLIC functions 
 */
EXTERN_FUNCTION (int window_done, (Xv_Window window));
EXTERN_FUNCTION (void xv_main_loop, (Xv_Window window));
EXTERN_FUNCTION (int window_read_event, (Xv_Window window, Event *event));
EXTERN_FUNCTION (void window_bell, (Xv_Window window));
EXTERN_FUNCTION (Xv_public int xv_rows, (Xv_Window window, int rows));
EXTERN_FUNCTION (Xv_public int xv_row, (Xv_Window window, int row));
EXTERN_FUNCTION (Xv_public int xv_cols, (Xv_Window window, int cols));
EXTERN_FUNCTION (Xv_public int xv_col, (Xv_Window window, int col));
EXTERN_FUNCTION (Xv_public int xv_send_message, (Xv_window window, Xv_opaque addresse, char *msg_type, int format, Xv_opaque *data, int len));

/*
 * PRIVATE functions
 */

EXTERN_FUNCTION (void window_set_cache_rect, (Xv_Window window, Rect *rect));

/*
 * The Xv_opaque type for the rect object list is defined in a private 
 * include file windowimpl.h
 */
EXTERN_FUNCTION (Window_rescale_rect_obj *window_create_rect_obj_list, (int num_elems));
EXTERN_FUNCTION (void window_destroy_rect_obj_list, (Window_rescale_rect_obj *rect_obj_list));
EXTERN_FUNCTION (void window_add_to_rect_list, (Window_rescale_rect_obj *rect_obj_list, Xv_Window window, Rect *rect, int i));
EXTERN_FUNCTION (int window_rect_equal_ith_obj, (Window_rescale_rect_obj *rect_obj_list, Rect *rect, int i));
EXTERN_FUNCTION (void window_set_client_message, ( Xv_Window window, XClientMessageEvent *msg));
EXTERN_FUNCTION (Xv_opaque * xv_get_selected_windows, (Xv_object window));

/*
 * PUBLIC functions 
 * For SunView 1 Compatibility
 */

EXTERN_FUNCTION (Xv_Window window_create, (Xv_Window window, Xv_pkg *pkg, DOTDOTDOT));
EXTERN_FUNCTION (Xv_opaque window_get, (Xv_Window window, Window_attribute attr, DOTDOTDOT));
EXTERN_FUNCTION (int window_set, (Xv_Window window, DOTDOTDOT));
EXTERN_FUNCTION (int window_destroy, (Xv_Window window));


/*
 * PRIVATE functions 
 * For SunView 1 Compatibility 
 */

EXTERN_FUNCTION (void window_rc_units_to_pixels, (Xv_Window win, Attr_avlist attr));

/*
 * Some wmgr stuff that needs to be here for the split libs.
 * This should be moved out as soon as all the pushpin stuff in moved
 * out of the intrinsic layer.  [csk 3/23/89]
 */

/* value for pushpin_initial_state */
#ifndef WMPushpinIsOut
#define WMPushpinIsOut  0
#endif /* WMPushpinIsOut */
#ifndef WMPushpinIsIn
#define WMPushpinIsIn   1
#endif /* WMPushpinIsIn */

#endif /* ~xview_window_DEFINED */
