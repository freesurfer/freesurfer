#ifndef lint
#ifdef sccs
static char     sccsid[] = "@(#)sel_pkg.h 1.5 90/11/13";
#endif
#endif

/*
 *	(c) Copyright 1990 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_selection_DEFINED
#define xview_selection_DEFINED

/*
 ***********************************************************************
 *		Include Files
 ***********************************************************************
 */

#include <sys/time.h>
#include <X11/X.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>
#include <X11/Xutil.h>
#include <xview/generic.h>
#include <errno.h>

/*
 ***********************************************************************
 *		Definitions and Macros
 ***********************************************************************
 */

/*
 * Package definitions
 */
#define SELECTION	&xv_sel_pkg
#define SELECTION_OWNER &xv_sel_owner_pkg
#define SELECTION_REQUESTOR &xv_sel_requestor_pkg
#define	SELECTION_ITEM	&xv_sel_item_pkg


/*
 * Various utility macros 
 */


/*
 * Public constants 
 */


/*
 * PRIVATE #defines 
 */

#define	SEL_ATTR(type, ordinal)	ATTR(ATTR_PKG_SELECTION, type, ordinal)
#define SEL_ATTR_LIST(ltype, type, ordinal) \
	SEL_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))


/* Errors */
#define SEL_ERROR              -1

#define SEL_BEGIN_MULTIPLE      2
#define SEL_END_MULTIPLE        4
#define SEL_MULTIPLE            8

#define SEL_BAD_TIME            0
#define SEL_BAD_WIN_ID          1
#define SEL_BAD_PROPERTY        2
#define SEL_BAD_CONVERSION      3
#define SEL_TIMEDOUT            4
#define SEL_PROPERTY_DELETED    5
#define SEL_BAD_PROPERTY_EVENT  6 

#define SEL_INCREMENT        2
#define SEL_EMPTY            1
#define SEL_PROPERTY_DATA    1

#define SEL_BUSY             1
#define SEL_LOSE             2
#define SEL_LOCAL_PROCESS    4
#define SEL_ADD_PROP_NOTIFY  8
#define SEL_INTERNAL_ERROR   16

#define OLD_SEL_CLIENT       2
#define NEW_SEL_CLIENT       4
    

#define MAX_NUM_INCR         20

#define MAX_SEL_BUFF_SIZE(dpy) ( XMaxRequestSize(dpy) )
#define BYTE_SIZE( len, format ) ( ( len * format ) >> 3 )


/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

/*
 * Typedefs 
 */
typedef	Xv_opaque 		Selection;
typedef	Selection 		Selection_owner;
typedef	Selection	 	Selection_requestor;
typedef	Xv_opaque 		Selection_item;

/*
 * Enumerations 
 */
typedef enum {

	/*
	 * Public Attributes 
	 */
	/* Common Selection package attributes */
	SEL_DATA		= SEL_ATTR(ATTR_OPAQUE,			   5),
	SEL_TYPE		= SEL_ATTR(ATTR_LONG,			  10),
	SEL_TYPE_NAME		= SEL_ATTR(ATTR_STRING,			  15),

	/* Selection object attributes */
	SEL_RANK		= SEL_ATTR(ATTR_LONG,			  20),
	SEL_RANK_NAME		= SEL_ATTR(ATTR_STRING,			  25),
	SEL_TIME		= SEL_ATTR(ATTR_OPAQUE,			  30),
	SEL_TIMEOUT_VALUE	= SEL_ATTR(ATTR_INT,			  35),

	/* Selection_owner object attributes */
	SEL_CONVERT_PROC	= SEL_ATTR(ATTR_FUNCTION_PTR,		  40),
	SEL_DONE_PROC		= SEL_ATTR(ATTR_FUNCTION_PTR,		  45),
	SEL_FIRST_ITEM		= SEL_ATTR(ATTR_OPAQUE,			  50),
	SEL_LOSE_PROC		= SEL_ATTR(ATTR_FUNCTION_PTR,		  55),
	SEL_NEXT_ITEM		= SEL_ATTR(ATTR_OPAQUE,			  60),
	SEL_OWN			= SEL_ATTR(ATTR_BOOLEAN,		  65),
	SEL_PROP_INFO		= SEL_ATTR(ATTR_OPAQUE,			  70),
	SEL_PROPERTY		= SEL_ATTR(ATTR_LONG,			  75),

	/* Selection_requestor object attributes */
	SEL_REPLY_PROC		= SEL_ATTR(ATTR_FUNCTION_PTR,		  80),
	SEL_TYPES		= SEL_ATTR_LIST(ATTR_NULL, ATTR_LONG,	  85),
	SEL_TYPE_NAMES		= SEL_ATTR_LIST(ATTR_NULL, ATTR_STRING,	  90),
	SEL_PROP_DATA		= SEL_ATTR(ATTR_OPAQUE,			  95),
	SEL_PROP_FORMAT		= SEL_ATTR(ATTR_INT,			  100),
	SEL_PROP_LENGTH		= SEL_ATTR(ATTR_LONG,			  105),
	SEL_PROP_TYPE		= SEL_ATTR(ATTR_LONG,			  110),
	SEL_PROP_TYPE_NAME	= SEL_ATTR(ATTR_STRING,			  115),
	SEL_TYPE_INDEX		= SEL_ATTR(ATTR_INT,			  120),
	SEL_APPEND_TYPES	= SEL_ATTR_LIST(ATTR_NULL, ATTR_LONG,	  125),
	SEL_APPEND_TYPE_NAMES	= SEL_ATTR_LIST(ATTR_NULL, ATTR_STRING,	  130),
	/* Selection_item object attributes */
	SEL_COPY		= SEL_ATTR(ATTR_BOOLEAN,		  135),
	SEL_FORMAT		= SEL_ATTR(ATTR_INT,			  140),
	SEL_LENGTH		= SEL_ATTR(ATTR_LONG,			  145)

	/*
	 * Private Attributes 
	 */
} Selection_attr;

/*
 * Structures 
 */

typedef struct {
    Xv_generic_struct	parent_data;
    Xv_opaque		private_data;
} Xv_sel;

typedef struct {
    Xv_sel		parent_data;
    Xv_opaque		private_data;
} Xv_sel_owner;

typedef struct {
    Xv_sel		parent_data;
    Xv_opaque		private_data;
} Xv_sel_requestor;

typedef struct {
    Xv_sel_owner	parent_data;
    Xv_opaque		private_data;
} Xv_sel_item;


/*
 ***********************************************************************
 *		Globals
 ***********************************************************************
 */

/*
 * 	Package Structures 
 */
extern Xv_pkg 		xv_sel_pkg;
extern Xv_pkg 		xv_sel_owner_pkg;
extern Xv_pkg 		xv_sel_requestor_pkg;
extern Xv_pkg 		xv_sel_item_pkg;

/*
 * 	Public Functions 
 */

EXTERN_FUNCTION (Bool sel_convert_proc, (Selection_owner sel_owner, Atom * type, Xv_opaque *data, unsigned long *length, int *format));
EXTERN_FUNCTION (void sel_post_request, (Selection_requestor sel_req));


typedef struct sel_prop_info {
    Xv_opaque	    data;
    int		    format;	/* data element size: 8, 16 or 32 bits */
    unsigned long   length;	/* nbr of elements in data */
    Atom            type;
    char	    *typeName;	
} Sel_prop_info;


typedef  struct sel_compat_info{
	Window   owner;
	Atom     selection;
	int      clientType;
	struct sel_compat_info *next;
} Sel_cmpat_info;

extern XContext  cmpatCtx;

#endif	/* ~xview_selection_DEFINED */
