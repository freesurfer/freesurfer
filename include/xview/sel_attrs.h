/*	@(#)sel_attrs.h 20.20 91/09/14		*/

#ifndef	xview_selection_attributes_DEFINED
#define	xview_selection_attributes_DEFINED

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/attr.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PRIVATE #defines 
 */

/*
 *	Common requests a client may send to a selection-holder
 */
#define ATTR_PKG_SELN	ATTR_PKG_SELN_BASE

#define SELN_ATTR(type, n)	ATTR(ATTR_PKG_SELN, type, n)

#define SELN_ATTR_LIST(list_type, type, n)	\
	ATTR(ATTR_PKG_SELN, ATTR_LIST_INLINE(list_type, type), n)

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

/*
 * Public Enumerations 
 */

/*
 *	Attributes of selections
 * 	The numbering scheme has to match the scheme of Sunview selection_attributes.h
 */
typedef enum	{
	/*
	 * Public Attributes 
	 */
	SELN_REQ_BYTESIZE	= SELN_ATTR(ATTR_INT,		         1),
	SELN_REQ_COMMIT_PENDING_DELETE	
				= SELN_ATTR(ATTR_NO_VALUE,	         65),
	SELN_REQ_CONTENTS_ASCII	= SELN_ATTR_LIST(ATTR_NULL, ATTR_CHAR,  2),

#ifdef OW_I18N
	SELN_REQ_CONTENTS_WCS   = SELN_ATTR_LIST(ATTR_NULL, ATTR_WCHAR, 202),
        SELN_REQ_CHARSIZE       = SELN_ATTR(ATTR_INT,                   204),
        SELN_REQ_FILE_NAME_WCS  = SELN_ATTR_LIST(ATTR_NULL, ATTR_WCHAR, 205),
#endif /*OW_I18N*/

	SELN_REQ_CONTENTS_PIECES= SELN_ATTR_LIST(ATTR_NULL, ATTR_CHAR,  3),
	SELN_REQ_DELETE		= SELN_ATTR(ATTR_NO_VALUE,	        66),
	SELN_REQ_END_REQUEST	= SELN_ATTR(ATTR_NO_VALUE,	        253),
	SELN_REQ_FAILED		= SELN_ATTR(ATTR_INT,		        255),
	SELN_REQ_FAKE_LEVEL	= SELN_ATTR(ATTR_INT,		        98),
	SELN_REQ_FILE_NAME	= SELN_ATTR_LIST(ATTR_NULL, ATTR_CHAR,  9),
	SELN_REQ_FIRST		= SELN_ATTR(ATTR_INT,		        4),
	SELN_REQ_FIRST_UNIT	= SELN_ATTR(ATTR_INT,		        5),
	SELN_REQ_LAST		= SELN_ATTR(ATTR_INT,		        6),
	SELN_REQ_LAST_UNIT	= SELN_ATTR(ATTR_INT,		        7),
	SELN_REQ_LEVEL		= SELN_ATTR(ATTR_INT,		        8),
	SELN_REQ_RESTORE	= SELN_ATTR(ATTR_NO_VALUE,	        67),
	SELN_REQ_SET_LEVEL	= SELN_ATTR(ATTR_INT,		        99),
	SELN_REQ_UNKNOWN	= SELN_ATTR(ATTR_INT,		        254),
	SELN_REQ_YIELD		= SELN_ATTR(ATTR_ENUM,		        97),
	/*
	 * Private Attributes 
	 */
#ifdef OW_I18N
        SELN_REQ_CONTENTS_CT    = SELN_ATTR_LIST(ATTR_NULL, ATTR_CHAR,  203),
#endif /*OW_I18N*/

	SELN_AGENT_INFO		= SELN_ATTR(ATTR_OPAQUE,                100),
	SELN_REQ_FUNC_KEY_STATE	= SELN_ATTR(ATTR_INT,		 	101),
	SELN_REQ_SELECTED_WINDOWS= SELN_ATTR_LIST(ATTR_NULL, ATTR_INT, 	102),
	SELN_REQ_CONTENTS_OBJECT= SELN_ATTR_LIST(ATTR_NULL, ATTR_CHAR, 	103),
	SELN_REQ_OBJECT_SIZE	= SELN_ATTR(ATTR_INT, 			104),
	SELN_REQ_IS_READONLY	= SELN_ATTR(ATTR_BOOLEAN,	       105),
	SELN_TRACE_ACQUIRE	= SELN_ATTR(ATTR_BOOLEAN,	       193),
	SELN_TRACE_DONE		= SELN_ATTR(ATTR_BOOLEAN,	       194),
	SELN_TRACE_DUMP		= SELN_ATTR(ATTR_ENUM,		       200),
	SELN_TRACE_HOLD_FILE	= SELN_ATTR(ATTR_BOOLEAN,	       195),
	SELN_TRACE_INFORM	= SELN_ATTR(ATTR_BOOLEAN,	       196),
	SELN_TRACE_INQUIRE	= SELN_ATTR(ATTR_BOOLEAN,	       197),
	SELN_TRACE_STOP		= SELN_ATTR(ATTR_BOOLEAN,	       199),
	SELN_TRACE_YIELD	= SELN_ATTR(ATTR_BOOLEAN,	       198)
}	Seln_attribute;

/* Meta-levels available for use with SELN_REQ_FAKE/SET_LEVEL.
 *	SELN_LEVEL_LINE is "text line bounded by newline characters,
 *			    including only the terminating newline"
 */
typedef enum {
	SELN_LEVEL_FIRST	= 0x40000001,
	SELN_LEVEL_LINE		= 0x40000101,
	SELN_LEVEL_ALL		= 0x40008001,
	SELN_LEVEL_NEXT		= 0x4000F001,
	SELN_LEVEL_PREVIOUS	= 0x4000F002
}	Seln_level;

#endif /* ~xview_selection_attributes_DEFINED */
