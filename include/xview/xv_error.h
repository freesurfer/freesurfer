/*
 * @(#)xv_error.h 1.25 91/09/14 SMI
 *
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_xv_error_DEFINED
#define xview_xv_error_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/attr.h>
#ifndef sunwindow_sun_DEFINED
#include <xview/sun.h>
#endif

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PRIVATE #defines 
 */

#define ERROR_ATTR(type, ordinal)      ATTR(ATTR_PKG_ERROR, type, ordinal)
#define ERROR_ATTR_LIST(ltype, type, ordinal) \
		ERROR_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))

#define ERROR_MAX_STRING_SIZE	512


/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

/*
 * Data Structures
 */

typedef enum {
	ERROR_RECOVERABLE,
	ERROR_NON_RECOVERABLE
} Error_severity;

typedef enum {
	ERROR_SYSTEM,
	ERROR_SERVER,
	ERROR_TOOLKIT,
	ERROR_PROGRAM
} Error_layer;

typedef enum {
	ERROR_BAD_ATTR	 	= ERROR_ATTR(ATTR_OPAQUE,		 1),
	ERROR_BAD_VALUE 	= ERROR_ATTR(ATTR_OPAQUE_PAIR,		 3),
	ERROR_CANNOT_GET	= ERROR_ATTR(ATTR_OPAQUE,		 6),
	ERROR_CANNOT_SET	= ERROR_ATTR(ATTR_OPAQUE,		 9),
	ERROR_CREATE_ONLY	= ERROR_ATTR(ATTR_OPAQUE,		12),
	ERROR_INVALID_OBJECT	= ERROR_ATTR(ATTR_STRING,		15),
	ERROR_LAYER		= ERROR_ATTR(ATTR_ENUM,			18),
	ERROR_PKG		= ERROR_ATTR(ATTR_OPAQUE,		21),
	ERROR_SERVER_ERROR	= ERROR_ATTR(ATTR_OPAQUE,		23),
	ERROR_SEVERITY		= ERROR_ATTR(ATTR_ENUM,			24),
	ERROR_STRING		= ERROR_ATTR(ATTR_STRING,		27)
} Error_attr;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * Functions 
 */

EXTERN_FUNCTION (int	xv_error, (Xv_object object, DOTDOTDOT));
EXTERN_FUNCTION (int 	xv_error_default, (Xv_object object, Attr_avlist avlist));
EXTERN_FUNCTION (char *	xv_error_format, (Xv_object object, Attr_avlist avlist));

#endif /* ~xview_xv_error_DEFINED */
