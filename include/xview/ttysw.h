/*	@(#)ttysw.h 20.15 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 * A tty subwindow is a subwindow type that is used to provide a
 * terminal emulation for teletype based programs.
 *
 * The caller of ttysw_start typically waits for the child process to die
 * before exiting.
 *
 */

#ifndef xview_ttysw_DEFINED
#define xview_ttysw_DEFINED

#include <xview/xv_c_types.h>
#include <xview/tty.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PRIVATE #defines 
 */

/* 
 * options - controlled by ttysw_getopt(), ttysw_setopt
 * The values of the #define's are wired into aliases, shell-scripts,
 * etc. and should not be changed!
 */
#define	TTYOPT_PAGEMODE			1
#define TTYOPT_SELSVC			3
#define TTYOPT_TEXT			4	/* TERMSW */

/*
 * styles for rendering boldface characters 
 */
#define TTYSW_BOLD_NONE			0x0
#define TTYSW_BOLD_OFFSET_X		0x1
#define TTYSW_BOLD_OFFSET_Y		0x2
#define TTYSW_BOLD_OFFSET_XY		0x4
#define TTYSW_BOLD_INVERT		0x8
#define TTYSW_BOLD_MAX			0x8

/*
 * Modes for invert and underline 
 */
#define TTYSW_ENABLE			0x0
#define TTYSW_DISABLE			0x1
#define TTYSW_SAME_AS_BOLD		0x2

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

typedef caddr_t	Ttysubwindow;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

EXTERN_FUNCTION (int ttysw_input, (Tty ttysw, char *addr, int len));
EXTERN_FUNCTION (int ttysw_output, (Tty ttysw, char *addr, int len));


#ifdef _OTHER_TTYSW_FUNCTIONS

/*
 * C Library routines specifically related to ttysw subwindow functions.
 */

/*
 * PRIVATE functions 
 */

EXTERN_FUNCTION (void ttysw_done, (Tty ttysw));
EXTERN_FUNCTION (void ttysw_setopt, (Tty ttysw, int opt, int on));
EXTERN_FUNCTION (int ttysw_getopt, (Tty ttysw, int opt));

/*
 * PUBLIC functions
 * for compatibility with pre-SunView 1 code
 */
EXTERN_FUNCTION (void ttysw_becomeconsole, (Tty ttysw));

#endif /* _OTHER_TTYSW_FUNCTIONS */

#endif /* ~xview_ttysw_DEFINED */
