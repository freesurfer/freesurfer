/*	@(#)cursor.h 20.34 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef  xview_cursor_DEFINED
#define  xview_cursor_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <stdio.h>
#include <xview/generic.h>
#include <xview/pkg.h>
#include <sys/types.h>
#include <pixrect/pixrect.h>

#ifdef __STDC__ 
#ifndef CAT
#define CAT(a,b)        a ## b 
#endif 
#endif
#include <pixrect/memvar.h>

#include <xview/rect.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

#define CURSOR				&xv_cursor_pkg

#define	CURSOR_ATTR(type, ordinal)	ATTR(ATTR_PKG_CURSOR, type, ordinal)

/*
 * For Sunview 1 compatibility
 */
#define CURSOR_MAX_IMAGE_BYTES		32	/* max. # of image bytes */
#define CURSOR_MAX_IMAGE_WORDS		16

#define CUR_MAXIMAGEBYTES		CURSOR_MAX_IMAGE_BYTES
#define CUR_MAXIMAGEWORDS		CURSOR_MAX_IMAGE_WORDS

#define NOFONTCURSOR			-1

/*
 * OPEN LOOK Cursor character definitions
 */
/* Source characters */
#define OLC_BASIC_PTR		0
#define OLC_BASIC_MASK_PTR	1
#define OLC_MOVE_PTR		2
#define OLC_MOVE_MASK_PTR	3
#define OLC_COPY_PTR		4
#define OLC_COPY_MASK_PTR	5
#define OLC_BUSY_PTR		6
#define OLC_BUSY_MASK_PTR	7
#define OLC_STOP_PTR		8 /* Note: Not an official OPEN LOOK cursor */
#define OLC_STOP_MASK_PTR	9
#define OLC_PANNING_PTR		10
#define OLC_PANNING_MASK_PTR	11
#define OLC_NAVIGATION_LEVEL_PTR	12
#define OLC_NAVIGATION_LEVEL_MASK_PTR	13

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef Xv_opaque	Xv_Cursor;
typedef Xv_opaque	Xv_cursor;

typedef enum {
    CURSOR_NEUTRAL = 0,
    CURSOR_ACCEPT = 1,
    CURSOR_REJECT = 2
} Cursor_drag_state;

typedef enum {
    CURSOR_MOVE = 0,
    CURSOR_DUPLICATE = 1
} Cursor_drag_type;

typedef enum {
	/*
 	 * Public Attributes 
 	 */
    CURSOR_BACKGROUND_COLOR	= CURSOR_ATTR(ATTR_SINGLE_COLOR_PTR,    27),
    CURSOR_DRAG_STATE		= CURSOR_ATTR(ATTR_ENUM,		2),
    CURSOR_DRAG_TYPE		= CURSOR_ATTR(ATTR_ENUM,		3),
    CURSOR_FOREGROUND_COLOR	= CURSOR_ATTR(ATTR_SINGLE_COLOR_PTR,    26),
    CURSOR_IMAGE		= CURSOR_ATTR(ATTR_PIXRECT_PTR,	  	1),
    CURSOR_MASK_CHAR		= CURSOR_ATTR(ATTR_INT,		  	5),
    CURSOR_OP			= CURSOR_ATTR(ATTR_INT,		 	10),
    CURSOR_SRC_CHAR		= CURSOR_ATTR(ATTR_INT,		 	15),
    CURSOR_STRING		= CURSOR_ATTR(ATTR_STRING,		16),
#ifdef OW_I18N
    CURSOR_STRING_WCS		= CURSOR_ATTR(ATTR_WSTRING,	 	17),
#endif
    CURSOR_XHOT			= CURSOR_ATTR(ATTR_INT,		 	20),
    CURSOR_YHOT			= CURSOR_ATTR(ATTR_INT,		 	25),
	/*
	 * OPEN LOOK cursors.  Used as index to XV_KEY_DATA on server.
	 */
    CURSOR_BASIC_PTR		= CURSOR_ATTR(ATTR_OPAQUE,		30),
    CURSOR_MOVE_PTR		= CURSOR_ATTR(ATTR_OPAQUE,		35),
    CURSOR_COPY_PTR		= CURSOR_ATTR(ATTR_OPAQUE,	 	40),
    CURSOR_BUSY_PTR		= CURSOR_ATTR(ATTR_OPAQUE,		45),
    CURSOR_PANNING_PTR		= CURSOR_ATTR(ATTR_OPAQUE,		50),
    CURSOR_NAVIGATION_LEVEL_PTR	= CURSOR_ATTR(ATTR_OPAQUE,		55),
    CURSOR_STOP_PTR		= CURSOR_ATTR(ATTR_OPAQUE,		60),
    CURSOR_QUESTION_MARK_PTR	= CURSOR_ATTR(ATTR_OPAQUE,		65)
} Cursor_attribute;
#define CURSOR_SHOW_CURSOR		XV_SHOW

typedef struct {
	Xv_generic_struct	parent_data;
	Xv_opaque		private_data;
}   Xv_cursor_struct;


/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

extern Xv_pkg		xv_cursor_pkg;

/*
 * Public Functions 
 */

EXTERN_FUNCTION (void		cursor_set_cursor, (Xv_object window, Xv_Cursor cursor));

/*
 * For Sunview 1 compatibility 
 */
EXTERN_FUNCTION (Xv_Cursor cursor_create, (Attr_attribute attr1, DOTDOTDOT));
EXTERN_FUNCTION (void		cursor_destroy, (Xv_Cursor cursor));
EXTERN_FUNCTION (Xv_Cursor	cursor_copy, (Xv_Cursor cursor));
EXTERN_FUNCTION (Xv_opaque	cursor_get, (Xv_Cursor cursor, Cursor_attribute  attr));
EXTERN_FUNCTION (int		cursor_set, (Xv_Cursor cursor, DOTDOTDOT));

#endif	/* xview_cursor_DEFINED */
