/*
 * @(#)tty.h 20.17 91/09/14 SMI
 *
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_tty_DEFINED
#define xview_tty_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/window.h>
#include <xview/openwin.h>
#include <xview/pkg.h>
#include <xview/attrol.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define TTY			&xv_tty_pkg

/*
 * Data type declaration for ttysw view 
 */
#define TTY_VIEW_TYPE		ATTR_PKG_TTY_VIEW
#define TTY_VIEW		&xv_tty_view_pkg

/*
 * PRIVATE #defines 
 */

#define TTY_ATTR(type, ordinal)	ATTR(ATTR_PKG_TTY, type, ordinal)
#define	ATTR_BUF_LEN_USED	ATTR_TYPE(ATTR_BASE_OPAQUE, 3)
#define	TTY_ARGV_DO_NOT_FORK	-1
#define	TTY_INFINITY		((long)0x77777777)

/*
 * PUBLIC #defines
 * For SunView 1 Compatibility Only
 */

/*
 * Data type declaration for ttysw folio 
 */
#define TTY_TYPE		ATTR_PKG_TTY

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures 
 ***********************************************************************
 */

typedef Xv_opaque	Tty;
typedef Xv_opaque	Tty_view;

typedef struct {
	Xv_openwin	parent_data;
	Xv_opaque	private_data;
} Xv_tty;

typedef struct {
	Xv_window_struct	parent_data;
	Xv_opaque		private_data;
} Xv_tty_view;
 
typedef enum {
	/*
	 * Public attributes 
	 */
	TTY_ARGV		= TTY_ATTR(ATTR_OPAQUE, 	 1),
	TTY_CONSOLE		= TTY_ATTR(ATTR_BOOLEAN,	 5),
	TTY_INPUT		= TTY_ATTR(ATTR_BUF_LEN_USED,	10),
	TTY_OUTPUT		= TTY_ATTR(ATTR_BUF_LEN_USED,	15),
	TTY_PAGE_MODE		= TTY_ATTR(ATTR_BOOLEAN,	20),
	TTY_QUIT_ON_CHILD_DEATH
				= TTY_ATTR(ATTR_BOOLEAN,	25),
	/*
	 * Private attributes 
	 */
	TTY_BOLDSTYLE		= TTY_ATTR(ATTR_INT,		30),
	TTY_BOLDSTYLE_NAME	= TTY_ATTR(ATTR_STRING,		35),
	TTY_INVERSE_MODE	= TTY_ATTR(ATTR_INT,		40),
	TTY_PID			= TTY_ATTR(ATTR_INT,		45),
	TTY_PTY_FD		= TTY_ATTR(ATTR_INT,		50),	/* --G */
	TTY_TTY_FD		= TTY_ATTR(ATTR_INT,		60),	/* --G */
	TTY_UNDERLINE_MODE	= TTY_ATTR(ATTR_INT,		65)
} Tty_attribute;

#undef ATTR_BUF_LEN_USED

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

extern  Xv_pkg		xv_tty_pkg;
extern  Xv_pkg		xv_tty_view_pkg;

/*
 * 		Escape sequences recognized by TTY subwindows
 *
 *      \E[1t           - open
 *      \E[2t           - close (become iconic)
 *      \E[3t           - move, with interactive feedback
 *      \E[3;TOP;LEFTt  - move, TOP LEFT in pixels
 *      \E[4t           - stretch, with interactive feedback
 *      \E[4;ROWS;COLSt - stretch, ROWS COLS in pixels
 *      \E[5t           - top (expose)
 *      \E[6t           - bottom (hide)
 *      \E[7t           - refresh
 *      \E[8;ROWS;COLSt - stretch, ROWS COLS in characters
 *      \E[11t          - report open or iconic, sends \E[1t or \E[2t
 *      \E[13t          - report position, sends \E[3;TOP;LEFTt
 *      \E[14t          - report size in pixels, sends \E[8;ROWS;COLSt
 *      \E[18t          - report size in chars, sends \E[4;ROWS;COLSt
 *      \E[20t          - report icon label, sends \E]Llabel\E\
 *      \E[21t          - report tool label, sends \E]llabel\E\
 *      \E]l<text>\E\   - set tool label to <text>
 *      \E]I<file>\E\   - set icon file to <file>
 *      \E]L<label>\E\  - set icon label to <label>
 */

#endif /* ~xview_tty_DEFINED  */
