/*      @(#)generic.h 20.40 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_generic_DEFINED
#define xview_generic_DEFINED

/*
 * Generic attributes fall into two classes:
 *	(1) Truly generic, implemented by attr_xxx.o or generic.o, use the
 * package ATTR_PKG_GENERIC, shared with attr.h.
 *	(2) Common but not truly generic, implemented by .o's spread
 * across many sub-systems, use the package ATTR_PKG_SV, shared with xview.h.
 * Many of these common attributes pertain to server properties and thus only
 * apply to objects with a window server component.
 *
 * Implementation dependent notes on generic X attributes:
 *	XV_XNAME has the format
 * "<host name>:<display number in decimal>:<xid in decimal>".
 *	XV_DEVICE_NUMBER is the XID of the underlying X object.  XV_XID is
 * provided when a piece of code wants to emphasize that the "X id" is what
 * is needed, rather than an abstract "tree link".
 * 	Most of these attributes are only supported on Drawable objects,
 * but some, like XV_XID, are supported by all objects that have direct
 * underlying X components, e.g. Fonts.
 */

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */


#include <xview/pkg_public.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define XV_GENERIC_OBJECT	&xv_generic_pkg

/* XV_COPY is "magic" package xv_create checks for to distinguish
 * creation of a new object from creation of a copy of an existing object.
 */
#define XV_COPY			(Xv_pkg *)1

/*
 * Accelerator for XV_HELP and HELP_STRING_FILENAME
 */
#define XV_HELP_DATA 		XV_KEY_DATA, XV_HELP
#define HELP_STRING_FILENAME	XV_KEY_DATA, XV_HELP_STRING_FILENAME

#define XV_XID			XV_DEVICE_NUMBER

#define	XV_DUMMY_WINDOW		0x77777777

/*
 * Focus Client Rank.  Value is of type Xv_focus_rank.
 * Referred to in the Mouseless Model Specification as "Focus Client Classes".
 * Transient and Limited focus classes are Window Manager objects.
 * An Ordinary Focus client has an XV_FOCUS_RANK of XV_FOCUS_SECONDARY.
 * A First Class Focus client has an XV_FOCUS_RANK of XV_FOCUS_PRIMARY.
 */
#define XV_FOCUS_RANK		XV_KEY_DATA, XV_FOCUS_RANK_KEY

#define	XV_RC_SPECIAL		0x20000
#define	XV_RESET_REF_COUNT	XV_REF_COUNT, XV_RC_SPECIAL
#define	XV_INCREMENT_REF_COUNT	XV_REF_COUNT, XV_RC_SPECIAL+1
#define	XV_DECREMENT_REF_COUNT	XV_REF_COUNT, XV_RC_SPECIAL-1


/*
 * PRIVATE #defines 
 */
#define	XV_ATTR(type, ordinal)	ATTR(ATTR_PKG_SV, type, ordinal)

#define XV_ATTR_LIST(ltype, type, ordinal) \
        XV_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))


/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef enum {
    XV_FOCUS_SECONDARY = 0,	/* default value (a.k.a., Ordinary Focus) */
    XV_FOCUS_PRIMARY = 1	/* a.k.a., First Class Focus */
} Xv_focus_rank;

/*
 * WARNING: GENERIC_ATTR shared with attr.h (which claims [0..64)) 
 */
typedef enum {
	/*
	 * PUBLIC and Generic 
	 */
	/*
	 * For "contexts": key & data (& optional destroy for data) 
	 */
	XV_KEY_DATA		= GENERIC_ATTR(ATTR_INT_PAIR,	 64),
	XV_KEY_DATA_COPY_PROC	= GENERIC_ATTR(ATTR_OPAQUE_PAIR, 65),
	XV_KEY_DATA_REMOVE	= GENERIC_ATTR(ATTR_INT,	 66), /* -S- */
	XV_KEY_DATA_REMOVE_PROC = GENERIC_ATTR(ATTR_OPAQUE_PAIR, 67),
	/*
	 * For "reference count" on shared objects, e.g. fonts & menus 
	 */
	XV_REF_COUNT		= GENERIC_ATTR(ATTR_INT,	 68),
	/*
	 * Type of object 
	 */
	XV_TYPE 		= GENERIC_ATTR(ATTR_OPAQUE,	 69), /* --G */
	XV_IS_SUBTYPE_OF 	= GENERIC_ATTR(ATTR_OPAQUE,	 70), /* --G */
	/*
	 * Miscellaneous 
	 */
	XV_LABEL		= GENERIC_ATTR(ATTR_STRING,	 71),
	XV_NAME			= GENERIC_ATTR(ATTR_STRING,	 72),
	XV_STATUS		= GENERIC_ATTR(ATTR_INT,	 73),
	XV_STATUS_PTR		= GENERIC_ATTR(ATTR_OPAQUE,	 74),
	XV_HELP			= GENERIC_ATTR(ATTR_STRING,	 80),
	XV_HELP_STRING_FILENAME	= GENERIC_ATTR(ATTR_STRING,	 82),
	XV_SHOW			= GENERIC_ATTR(ATTR_BOOLEAN,	 81),
#ifdef OW_I18N
        XV_LABEL_WCS            = GENERIC_ATTR(ATTR_WSTRING,     164),
        XV_NAME_WCS             = GENERIC_ATTR(ATTR_WSTRING,     161),
        XV_HELP_WCS             = GENERIC_ATTR(ATTR_WSTRING,     162),
        XV_HELP_STRING_FILENAME_WCS
                                = GENERIC_ATTR(ATTR_WSTRING,     163),
#endif /* OW_I18N */
	/*
	 * Required by package implementations, used only by xv_create 
	 */
	XV_COPY_OF		= GENERIC_ATTR(ATTR_OPAQUE,	 75), /* -S- */
	XV_END_CREATE		= GENERIC_ATTR(ATTR_NO_VALUE,	 76), /* -S- */
	/*
	 * To simplify SunView1.X compatibility 
	 */
	XV_SELF			= GENERIC_ATTR(ATTR_OPAQUE,	 77), /* --G */
	/*
	 * Managing (usually containing) object 
	 */
	XV_OWNER		= GENERIC_ATTR(ATTR_OPAQUE,	 78),
    	/*
	 * Required by package implementations, used only by xv_find 
	 */
	XV_AUTO_CREATE		= GENERIC_ATTR(ATTR_INT,	 79), /* C-- */
	/*
	 * PUBLIC but only Common 
	 */
	/*
	 * For layout 
	 */
	XV_FONT			= XV_ATTR(ATTR_OPAQUE,		 64),
	XV_MARGIN		= XV_ATTR(ATTR_INT,		 65),
	XV_LEFT_MARGIN		= XV_ATTR(ATTR_INT,		 66),
	XV_TOP_MARGIN		= XV_ATTR(ATTR_INT,		 67),
	XV_RIGHT_MARGIN		= XV_ATTR(ATTR_INT,		 68),
	XV_BOTTOM_MARGIN	= XV_ATTR(ATTR_INT,		 69),
	/*
	 * Origin is usually parent's most upper-left coord inside margins 
	 */
	XV_X			= XV_ATTR(ATTR_X,		 70),
	XV_Y			= XV_ATTR(ATTR_Y,		 71),
	XV_WIDTH		= XV_ATTR(ATTR_X,		 72),
	XV_HEIGHT		= XV_ATTR(ATTR_Y,		 73),
	XV_RECT			= XV_ATTR(ATTR_RECT_PTR,	 74),
	/*
	 * Server specific or dependent 
	 */
	XV_XNAME		= XV_ATTR(ATTR_STRING,		 96), /* C-G */
	XV_DEVICE_NUMBER	= XV_ATTR(ATTR_LONG,		 97), /* C-G */
	XV_ROOT			= XV_ATTR(ATTR_OPAQUE,		 98), /* --G */
	XV_VISUAL		= XV_ATTR(ATTR_OPAQUE,          125), /* C-G */
	XV_VISUAL_CLASS		= XV_ATTR(ATTR_INT,		117), /* C-G */
	XV_DEPTH		= XV_ATTR(ATTR_INT,		126), /* C-G */
	XV_DISPLAY		= XV_ATTR(ATTR_OPAQUE,		110), /* --G */
	XV_SCREEN		= XV_ATTR(ATTR_OPAQUE,		116), /* --G */
	/*
	 * Mouseless Model support
	 */
	XV_FOCUS_ELEMENT	= XV_ATTR(ATTR_INT,		118), /* -S- */
	XV_FOCUS_RANK_KEY	= XV_ATTR(ATTR_ENUM,		119), /* CSG */

	/*
	 * Added to support the Xrm resource database
	 */
	XV_USE_DB		= XV_ATTR_LIST(ATTR_RECURSIVE, ATTR_AV,	120),
	XV_INSTANCE_NAME	= XV_ATTR(ATTR_STRING,		125),
	XV_INSTANCE_QLIST	= XV_ATTR(ATTR_OPAQUE,		130),
	XV_QUARK		= XV_ATTR(ATTR_OPAQUE,		135),
	XV_USE_INSTANCE_RESOURCES=XV_ATTR(ATTR_OPAQUE,		140),

#ifdef OW_I18N
	/*
	 * The I18N Level 4 attribute XV_IM goes here:
	 */
	 XV_IM                   = XV_ATTR(ATTR_OPAQUE,          150),
#endif /* OW_I18N */

	/*
	 * Added to support locale announcement
	 */
	XV_LC_BASIC_LOCALE	= XV_ATTR(ATTR_STRING, 		155),
	XV_LC_DISPLAY_LANG	= XV_ATTR(ATTR_STRING, 		160),
	XV_LC_INPUT_LANG	= XV_ATTR(ATTR_STRING, 		165),
	XV_LC_NUMERIC		= XV_ATTR(ATTR_STRING, 		170),
	XV_LC_TIME_FORMAT	= XV_ATTR(ATTR_STRING,		175),
 	XV_LOCALE_DIR		= XV_ATTR(ATTR_STRING, 		180),
	XV_USE_LOCALE		= XV_ATTR(ATTR_BOOLEAN,		185),

	/*
	 * PRIVATE now, but ... 
	 */
	XV_GC			= XV_ATTR(ATTR_OPAQUE,		113)  /* --G */
} Xv_generic_attr;

/*
 * Generic package definition	
 */
typedef struct {
    Xv_base	parent_data;
    Xv_opaque	private_data;
} Xv_generic_struct;

typedef enum {
    
    /* XXX: This is a hack to support V3 applications dropping on V2 apps. 
     * You ask what does XV_INIT_* have to do with dnd, well...We needed
     * some way to tell the difference between a V2 client using V3
     * shared libs and a V3 client using V3 shared libs.  By changing
     * the value of the attr XV_INIT_ARGS{_PTR_ARGV} we could detect this.
     * A V2 client using V3 libs will still have the old attr values
     * compiled into them, while V3 client compiled under V3 header files
     * will use the new value for these attrs. It's a hack, but it works.
     * We can remove this for V4 as we only need to support one release
     * back.  Everything marked with DND_HACK can be removed generic.h
     * xv_init.c, frame_init.c
     */

    /* XXX: This is a hack to support the navigator which shipped on
     * pre-V3 FCS.  The navigator has a bug in which it expects XView to
     * send a WIN_RESIZE event when certain window objects are created.
     * XView no longer sends bogus WIN_RESIZE events.  This breaks the
     * navigator.  This hack can be removed in a future release of XView.
     */

    /* Private Attributes */

    /* DND_HACK begin */
    /* The following two attributes are private to XView.  Do not use them as
     * they will be removed in a future release of the toolkit.
     */
    XV_INIT_ARGS_FOR_DND     = XV_ATTR(ATTR_INT_PAIR,         	3),
    XV_INIT_ARGC_PTR_ARGV_FOR_DND = XV_ATTR(ATTR_INT_PAIR,	6),  /* -S- */
    /* DND_HACK end */

    /* NAV_HACK begin */
    /* The following two attributes are private to XView.  Do not use them as
     * they will be removed in a future release of the toolkit.
     */
    XV_INIT_ARGS_NAV          = XV_ATTR(ATTR_INT_PAIR,         	4),
    XV_INIT_ARGC_PTR_ARGV_NAV = XV_ATTR(ATTR_INT_PAIR,         	7),  /* -S- */
    /* NAV_HACK end */

    /* Public Attributes */

    XV_INIT_ARGS             = XV_ATTR(ATTR_INT_PAIR,         	5),
    XV_INIT_ARGC_PTR_ARGV    = XV_ATTR(ATTR_INT_PAIR,         	8),  /* -S- */
    XV_USAGE_PROC       = XV_ATTR(ATTR_FUNCTION_PTR,     	9),  /* -S- */
    XV_ERROR_PROC       = XV_ATTR(ATTR_FUNCTION_PTR,    	12),
    XV_X_ERROR_PROC	= XV_ATTR(ATTR_FUNCTION_PTR,    	15)
} Xv_attr;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

extern Xv_pkg		xv_generic_pkg;

/*
 * PUBLIC functions 
 */

EXTERN_FUNCTION (Xv_object xv_init, (Attr_attribute attr1, DOTDOTDOT));
EXTERN_FUNCTION (Attr_attribute xv_unique_key, (void));

#endif /* ~xview_generic_DEFINED */
