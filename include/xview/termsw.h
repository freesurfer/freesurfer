/*
 * @(#)termsw.h 20.17 91/09/14 SMI      
 *
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_termsw_DEFINED
#define xview_termsw_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#ifdef XV_ATTRIBUTES_ONLY
#include <xview/attrol.h>

#else
#include <xview/textsw.h>

#endif /* XV_ATTRIBUTES_ONLY */

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#ifndef XV_ATTRIBUTES_ONLY

#define TERMSW     		     	&xv_termsw_pkg
#define TERMSW_VIEW			&xv_termsw_view_pkg

#endif	/* ~XV_ATTRIBUTES_ONLY */

/*
 * PRIVATE #defines 
 */

#define TERMSW_TYPE      		ATTR_PKG_TERMSW
#define TERMSW_VIEW_TYPE		ATTR_PKG_TERMSW_VIEW
#define	TERMSW_ATTR(type, ordinal)	ATTR(ATTR_PKG_TERMSW, type, ordinal)

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

#ifndef XV_ATTRIBUTES_ONLY

typedef Xv_opaque		Termsw;
typedef Xv_opaque  		Termsw_view;

typedef struct {
	/*
	 * This isn't really a textsw, only shares few attrs 
	 */
        Xv_textsw		parent_data; 
        Xv_opaque    		private_data;
	Xv_opaque		private_text;
	Xv_opaque		private_tty;
} Xv_termsw;

typedef struct {
	/*
	 * This isn't really a textsw view, only shares few attrs 
	 */
        Xv_textsw_view    	parent_data; 
        Xv_opaque  		private_data;
        Xv_opaque		private_text;
	Xv_opaque		private_tty;
} Xv_termsw_view;

#endif /* ~XV_ATTRIBUTES_ONLY */

typedef enum {
	TERMSW_MODE_TYPE,
	TTYSW_MODE_TYPE
} Termsw_mode;

typedef enum {
	/*
 	 * Public attributes. 
	 */
	TERMSW_MODE 	= TERMSW_ATTR(ATTR_INT,  	1)	/* --G */
} Termsw_attribute;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

extern  Xv_pkg			xv_termsw_pkg;
extern  Xv_pkg       		xv_termsw_view_pkg;

#endif /* ~xview_termsw_DEFINED */
