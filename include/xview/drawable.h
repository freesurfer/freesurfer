/*      @(#)drawable.h 20.14 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_drawable_DEFINED
#define xview_drawable_DEFINED

/*
 * Interface to generic attributes of Drawable objects, where "drawable"
 * is defined by X server.  This is currently an implementation concept.
 */

/*
 ***********************************************************************
 *			Include files
 ***********************************************************************
 */

#include <xview/generic.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PRIVATE #defines 
 */

#define XV_DRAWABLE_OBJECT		&xv_drawable_pkg

#define DRAWABLE_ATTR(type, ordinal)	ATTR(ATTR_PKG_DRAWABLE, type, ordinal)

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef Xv_opaque    Xv_Drawable;
typedef Xv_opaque    Xv_drawable;

typedef enum {
	/*
	 * Private Attributes
	 */
    DRAWABLE_INFO	= DRAWABLE_ATTR(ATTR_OPAQUE,	100)
} Drawable_attr;

typedef struct {			/* For sub-pkg implementors only */
    Xv_generic_struct	parent_data;
    Xv_opaque		private_data;
} Xv_drawable_struct;

/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

extern Xv_pkg			xv_drawable_pkg;

#endif	/* ~xview_drawable_DEFINED */
