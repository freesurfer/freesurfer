/*	@(#)attrol.h 20.11 88/09/05	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_attrol_DEFINED
#define	xview_attrol_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/base.h>
#include <xview/attr.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * These are the packages that are part of the XView OL library.
 *
 * IMPORTANT NOTE:  The attr id numbers start where the "intrinsics"
 *		    attr ids leave off.  The range of valid values
 *		    for these objects is [64..128) where 
 *		    ATTR_PKG_LAST_VALUE must be less than 128 and all
 *		    ATTR_PKG_ ids defined here must fit in the range 
 *		    (ATTR_PKG_LAST_VALUE..128).
 *
 *		    Be sure to check the value of ATTR_PKG_LAST_VALUE
 *		    when adding any new packages.
 */

#define ATTR_PKG_CANVAS		     (ATTR_PKG_LAST_VALUE +  1)
#define ATTR_PKG_ENTITY   	     (ATTR_PKG_LAST_VALUE +  2)
#define ATTR_PKG_TERMSW		     (ATTR_PKG_LAST_VALUE +  3)
#define ATTR_PKG_FRAME		     (ATTR_PKG_LAST_VALUE +  4)	
#define ATTR_PKG_ICON		     (ATTR_PKG_LAST_VALUE +  5)
#define ATTR_PKG_MENU		     (ATTR_PKG_LAST_VALUE +  6)
#define ATTR_PKG_PANEL		     (ATTR_PKG_LAST_VALUE +  7)
#define ATTR_PKG_OPENWIN             (ATTR_PKG_LAST_VALUE +  8)
#define ATTR_PKG_TEXTSW		     (ATTR_PKG_LAST_VALUE +  9)
#define ATTR_PKG_TTY		     (ATTR_PKG_LAST_VALUE + 10)
#define ATTR_PKG_NOTICE		     (ATTR_PKG_LAST_VALUE + 11)
#define ATTR_PKG_HELP		     (ATTR_PKG_LAST_VALUE + 12)
#define ATTR_PKG_TEXTSW_VIEW	     (ATTR_PKG_LAST_VALUE + 13)
#define ATTR_PKG_PANEL_VIEW	     (ATTR_PKG_LAST_VALUE + 14)
#define ATTR_PKG_CANVAS_VIEW         (ATTR_PKG_LAST_VALUE + 15)
#define ATTR_PKG_CANVAS_PAINT_WINDOW (ATTR_PKG_LAST_VALUE + 16)
#define ATTR_PKG_TTY_VIEW 	     (ATTR_PKG_LAST_VALUE + 17)
#define ATTR_PKG_TERMSW_VIEW 	     (ATTR_PKG_LAST_VALUE + 18)
#define ATTR_PKG_SCROLLBAR 	     (ATTR_PKG_LAST_VALUE + 19)

/* See REMIND in attr.h before adding any new pkgs. */

#endif /* ~xview_attrol_DEFINED */

