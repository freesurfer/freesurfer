/*	@(#)win_env.h 20.14 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_win_environ_DEFINED
#define xview_win_environ_DEFINED

#include <xview/xv_c_types.h>

/*
 ***********************************************************************
 *                      Definitions and Macros
 ***********************************************************************
 */

/*
 * PRIVATE #defines
 */

/*
 * The window library uses the ENVIRONMENT to pass a limited amount of
 * window information to newly created child processes.
 * A collection of conventions dictate the use of these values.
 */

#define WE_PARENT               "WINDOW_PARENT"
#define WE_INITIALDATA          "WINDOW_INITIALDATA"
#define WE_GFX                  "WINDOW_GFX"
#define WE_ME                   "WINDOW_ME"

/*
 ***********************************************************************
 *                              Globals
 ***********************************************************************
 */
 
/*
 * C Library routines specifically related to ENVIRONMENT conventions.
 */
 
/*
 * PUBLIC functions
 */
 
/*
 * Get/set window that can be taken over by graphics programs.
 */
#ifdef _OTHER_WIN_ENV_FUNCTIONS

EXTERN_FUNCTION (void we_setgfxwindow, (char *windevname));
EXTERN_FUNCTION (int we_getgfxwindow, (char *windevname)); 
 
/*
 * PUBLIC functions provided only for
 * compatibility with pre-SunView code
 */
 
/*
 * Get/set parent of window that is being created.
 */
EXTERN_FUNCTION (void we_setparentwindow, (char *windevname));
EXTERN_FUNCTION (int we_getparentwindow, (char *windevname));

#endif /* _OTHER_WIN_ENV_FUNCTIONS */
 
#endif /* ~xview_win_environ_DEFINED */

