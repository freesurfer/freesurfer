/*	@(#)sun.h 20.20 91/09/14 SMI */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef sunwindow_sun_DEFINED
#define sunwindow_sun_DEFINED

#include <xview/base.h>
#include <sys/types.h>
#include <stdio.h> 

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif

#ifndef Bool_DEFINED
#define Bool int
#define True 1
#define False 0
#endif

#define strequal(s1, s2) (strcmp(s1, s2) == 0)
#define STRDUP(str)	 (strcpy(xv_malloc((unsigned) strlen(str) + 1), str))
/*
 * Get some storage and copy a string. Note that str is evaluated twice,
 * so no side effects.
 */
#define ord(e)  ((int) (e))
/*
 * Make an enum usable as an integer
 * Successor or predecessor macros might also be convenient, but
 * much less so.
 */

#define FOREVER for (;;)

#ifndef MIN
#define MIN(x, y) ( ((x) < (y)) ? (x) : (y) )
#endif

#ifndef MAX
#define MAX(x, y) ( ((x) > (y)) ? (x) : (y) )
#endif

#endif
