/*    @(#)attr.h 20.42 91/06/17  */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_attr_DEFINED
#define	xview_attr_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <sys/types.h>
#include <xview/base.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * This macro is machine dependent in that it assumes
 * the cardinality will be in the lower 5 bits of the type-cardinality
 * pair.
 */
#define	ATTR_TYPE(base_type, cardinality) \
    (((((unsigned)(base_type)) & 0x07f) << 5) |	\
     (((unsigned)(cardinality)) & 0x0f))

/*
 * Note that this macro is machine dependent in that it assumes the
 * base_type-cardinality pair will be in the lower 13 bits of the 
 * list_type-base_cardinality pair.
 */
#define	ATTR_LIST_OF(list_type, list_ptr_type, base_cardinality) \
    (((((unsigned)(list_type)) & 0x3) << 14) | \
     (((unsigned)(list_ptr_type) & 0x1) << 13) | \
     (((unsigned)(base_cardinality)) & 0x3fff))

#define	ATTR_LIST_INLINE(list_type, base_cardinality)	\
    ATTR_LIST_OF(list_type, ATTR_LIST_IS_INLINE, base_cardinality)

#define	ATTR_LIST_PTR(list_type, base_cardinality)	\
    ATTR_LIST_OF(list_type, ATTR_LIST_IS_PTR, base_cardinality)

/*
 * An attribute is composed of 
 * pkg     ( 8 bits): id of package that uses the attribute
 * ordinal ( 8 bits): ordinal number of the attribute
 * type    (16 bits): list type, list pointer type, base type, and cardinality
 */
#define	ATTR(pkg, type, ordinal)	\
    ( ((((unsigned)(pkg))	& 0x7f) << 24) | \
      ((((unsigned)(ordinal))	& 0xff) << 16) | \
       (((unsigned)(type))	& 0xffef)	)


/*
 *  The base type space potentially runs from 0 to 127 inclusive.  This is
 * subdivided as follows:
 *	[0..32)	[ATTR_BASE_UNUSED_FIRST..ATTR_BASE_UNUSED_LAST]
 *		Reserved for non-Sun packages.
 *     [32..64)	Reserved for future use.
 *    [64..128)	[ATTR_BASE_FIRST..ATTR_BASE_LAST] is currently used by
 *		Sun packages
 *		(ATTR_BASE_LAST..128) is for future use by Sun packages.
 */
#define	ATTR_BASE_UNUSED_FIRST		0
#define	ATTR_BASE_UNUSED_LAST		31
#define	ATTR_BASE_FIRST			64

#define	ATTR_PKG_UNUSED_FIRST		1
#define ATTR_PKG_UNUSED_LAST		31
#define	ATTR_PKG_FIRST			64

/*
 * Values from the Attr_pkg enumeration which are required for some
 * of the following #defines. This allows those defines to remain in this
 * section rather than below the enumeration.
 */
#define ATTR_PKG_ZERO_VALUE		0
#define ATTR_PKG_GENERIC_VALUE		ATTR_PKG_FIRST
#define ATTR_PKG_LAST_VALUE		ATTR_PKG_FIRST + 14

/*
 * ATTR_STANDARD_SIZE is large enough to allow for 
 * most attribute-value lists.
 */
#define	ATTR_STANDARD_SIZE		250

/*
 * Generic attribute definition	
 */
#define GENERIC_ATTR(type, ordinal)	\
				ATTR(ATTR_PKG_GENERIC_VALUE, type, ordinal)

/*
 * Note that in these macros, attr must not be 
 * in a register or be a constant.
 * Since this is deemed too restrictive, we use
 * shifting & masking instead.  
 * #define ATTR_UNION(attr)	((Attr_union *) &((Attr_attribute) (attr)))
 * #define ATTR_INFO(attr)		(ATTR_UNION(attr)->info)
 * #define ATTR_CODE(attr)		(ATTR_UNION(attr)->code)
 * #define ATTR_LIST_TYPE(attr)		(ATTR_INFO(attr).list_type)
 * #define ATTR_LIST_PTR_TYPE(attr)	(ATTR_INFO(attr).list_ptr_type)
 * #define ATTR_BASE_TYPE(attr)		(ATTR_INFO(attr).base_type)
 * #define ATTR_CARDINALITY(attr)	(ATTR_INFO(attr).cardinality)
 * #define ATTR_CONSUMED(attr)		(ATTR_INFO(attr).consumed)
 * #define ATTR_PKG(attr)		(ATTR_INFO(attr).pkg)
 * #define ATTR_ORDINAL(attr)		(ATTR_INFO(attr).ordinal)
 */

#define	ATTR_CODE(attr)		((unsigned) (attr))

#define	ATTR_PKG(attr) \
	((Attr_pkg) ((ATTR_CODE(attr) >> 24) & 0xFF))

#define	ATTR_ORDINAL(attr)	\
	((unsigned) ((ATTR_CODE(attr) >> 16) & 0xFF))

#define	ATTR_LIST_TYPE(attr)	\
	((Attr_list_type) ((ATTR_CODE(attr) >> 14) & 0x3))

#define	ATTR_LIST_PTR_TYPE(attr)	\
	((Attr_list_ptr_type) ((ATTR_CODE(attr) >> 13) & 0x1))

#define	ATTR_BASE_TYPE(attr)	\
	((Attr_base_type) ((ATTR_CODE(attr) >> 5) & 0x7F))

/*
 * This macro is different from ATTR_BASE_TYPE in that it returns
 * the 'type' of the attribute, and not just the base type.
 * The base type is one of the things that is used to specify 'type'.
 *
 * It is also different from ATTR_TYPE, which returns the 
 * <base type><spare1><cardinality> bits.
 *
 * Note: Bit 4 i.e. the spare1 bit is zeroed.
 */
#define ATTR_WHICH_TYPE(attr)	\
	((Attr_base_cardinality) ((ATTR_CODE(attr)) & 0xffef))

#define	ATTR_CARDINALITY(attr)	\
	((unsigned) (ATTR_CODE(attr) & 0xF))

#define ATTR_CONSUMED(attr)	((unsigned) ((ATTR_CODE(attr) >> 12) & 0x1))

#define	attr_skip(attr, argv)	\
    ((ATTR_LIST_TYPE((attr)) == ATTR_NONE) \
	? (Attr_avlist) (argv) + ATTR_CARDINALITY((attr)) \
	: attr_skip_value((Attr_attribute)(attr), (argv)))

#define	attr_next(attrs)	attr_skip((*(attrs)), ((attrs)+1))

#define attr_make( listhead, listlen, valist )   \
           attr_make_count( listhead, listlen, valist, NULL )

/*
 * Character unit support 
 * Provided for SunView 1 compatibility.
 */
#ifndef	lint
#define attr_replace_cu(avlist, font, lmargin, tmargin, rgap) 	\
    attr_rc_units_to_pixels(avlist, xv_get(font, FONT_DEFAULT_CHAR_WIDTH), 	\
	xv_get(font, FONT_DEFAULT_CHAR_HEIGHT), lmargin, tmargin, 0, rgap)

#define attr_cu_to_x(encoded_value, font, left_margin) 		\
    attr_rc_unit_to_x(encoded_value, xv_get(font, FONT_DEFAULT_CHAR_WIDTH), left_margin, 0)

#define attr_cu_to_y(encoded_value, font, top_margin, row_gap) 	\
    attr_rc_unit_to_y(encoded_value, xv_get(font, FONT_DEFAULT_CHAR_HEIGHT), top_margin, \
	row_gap)
#endif	/* lint */

#define	ATTR_CU_TAG		0x80000000
#define ATTR_PIXEL_OFFSET	0x00008000
#define	ATTR_CU_MASK		0xC0000000

#define	ATTR_CU_TYPE(n)					\
    ((Attr_cu_type) ((n) & (unsigned) (ATTR_CU_LENGTH)))

#define	ATTR_CU(unit, n)				\
   (((unsigned)(unit)) | (((n) & 0x1FFF) << 16) |	\
    ATTR_CU_TAG | ATTR_PIXEL_OFFSET)

/*
 * attr_is_cu(n) returns non-zero if n has 
 * been encoded using ATTR_CU()
 */
#define	attr_is_cu(n)		(((n) & ATTR_CU_MASK) == ATTR_CU_TAG)

/*
 * Following are useful for multi-pass avlist processing. 
 */
#define	ATTR_NOP(attr)		(ATTR_CODE(attr) | (0x1 << 12))
#define	ATTR_CONSUME(attr)	(attr) = ((Xv_opaque)ATTR_NOP(attr))

/*
 * For Sunview 1 compatibility
 */

#define	ATTR_PIXWIN_PTR	ATTR_OPAQUE

/*
 * Macros for position including
 * margins.
 */
#define	ATTR_COL(n)		ATTR_CU(ATTR_CU_POSITION, n)
#define	ATTR_ROW(n)		ATTR_CU(ATTR_CU_POSITION, n)

/*
 * Macros for length excluding
 * margins.
 */
#define	ATTR_COLS(n)		ATTR_CU(ATTR_CU_LENGTH, n)
#define	ATTR_ROWS(n)		ATTR_CU(ATTR_CU_LENGTH, n)
#define	ATTR_CHARS(n)		ATTR_CU(ATTR_CU_LENGTH, n)
#define	ATTR_LINES(n)		ATTR_CU(ATTR_CU_LENGTH, n)

/*
 ***********************************************************************
 *			Typedefs, enumerations, and structs
 ***********************************************************************
 */

/*
 * Attr_avlist is not an array of Attr_attributes, because it is an array
 * of intermixed attributes and values.
 */
typedef unsigned long	 Attr_attribute;	/* 32 bit quantity */
typedef Attr_attribute	*Attr_avlist;

/*
 * Enumerations 
 */

typedef enum {
    ATTR_CU_POSITION	= 0x0,			/* bit 29 is off */
    ATTR_CU_LENGTH	= 0x20000000		/* bit 29 is on */
} Attr_cu_type;

typedef enum {
    ATTR_LIST_IS_INLINE	= 0,
    ATTR_LIST_IS_PTR	= 1
} Attr_list_ptr_type;

typedef enum {
    /* Note that ATTR_NONE must have a value of zero,
     * since a no-list type is assumed for each of the
     * types in Attr_base_cardinality.
     */
    ATTR_NONE		= 0,
    ATTR_RECURSIVE	= 1,
    ATTR_NULL		= 2,
    ATTR_COUNTED	= 3
} Attr_list_type;

/*
 * NOTE: The base type numbers have to be EXACTLY the same as SunView1 in order
 *	  to support cut and paste between SunView1 and XView windows.
 *	  Nothing changes!
 */
typedef enum {
    ATTR_BASE_NO_VALUE		= ATTR_BASE_FIRST + 17,
	/*
	 * Fundamental C types. 
	 */
    ATTR_BASE_INT		= ATTR_BASE_FIRST,
    ATTR_BASE_LONG		= ATTR_BASE_FIRST + 24,
    ATTR_BASE_SHORT		= ATTR_BASE_FIRST + 25,
    ATTR_BASE_ENUM		= ATTR_BASE_FIRST + 9,
    ATTR_BASE_CHAR		= ATTR_BASE_FIRST + 10,
#ifdef OW_I18N
    ATTR_BASE_WCHAR		= ATTR_BASE_FIRST + 26,
#define ATTR_BASE_CHAR_WC	ATTR_BASE_WCHAR
#endif    
    ATTR_BASE_STRING		= ATTR_BASE_FIRST + 11,
#ifdef OW_I18N
    ATTR_BASE_WSTRING		= ATTR_BASE_FIRST + 27,
#define	ATTR_BASE_STRING_WCS	ATTR_BASE_WSTRING
#endif    
    ATTR_BASE_FUNCTION_PTR	= ATTR_BASE_FIRST + 19,
	/*
	 * Derivative C types. 
	 */
    ATTR_BASE_BOOLEAN		= ATTR_BASE_FIRST + 8,
    ATTR_BASE_OPAQUE		= ATTR_BASE_FIRST + 16,
	/*
	 * Special coordinate types; look in attr_cu.c for the details. 
	 */
    ATTR_BASE_X			= ATTR_BASE_FIRST + 2,
    ATTR_BASE_INDEX_X		= ATTR_BASE_FIRST + 3,
    ATTR_BASE_Y			= ATTR_BASE_FIRST + 4,
    ATTR_BASE_INDEX_Y		= ATTR_BASE_FIRST + 5,
    ATTR_BASE_XY		= ATTR_BASE_FIRST + 6,
    ATTR_BASE_INDEX_XY		= ATTR_BASE_FIRST + 7,	
	/*
	 * Pointer types. 
	 */
    ATTR_BASE_PIXRECT_PTR	= ATTR_BASE_FIRST + 12,
    ATTR_BASE_PIXFONT_PTR	= ATTR_BASE_FIRST + 13,
    ATTR_BASE_RECT_PTR		= ATTR_BASE_FIRST + 15,
    ATTR_BASE_AV		= ATTR_BASE_FIRST + 18,
    ATTR_BASE_ICON_PTR		= ATTR_BASE_FIRST + 20,
    ATTR_BASE_SINGLE_COLOR_PTR	= ATTR_BASE_FIRST + 21,
    ATTR_BASE_CURSOR_PTR	= ATTR_BASE_FIRST + 22,
#ifdef OW_I18N
    ATTR_BASE_LAST		= ATTR_BASE_FIRST + 27
#else
    ATTR_BASE_LAST		= ATTR_BASE_FIRST + 25
#endif
} Attr_base_type;

/* Clients of the attribute value package should use
 * Attr_base_cardinality elements to define the base type
 * and cardinality of their attributes.
 */  
typedef enum {
    ATTR_NO_VALUE		= ATTR_TYPE(ATTR_BASE_NO_VALUE, 	0),
    ATTR_INT			= ATTR_TYPE(ATTR_BASE_INT, 		1),
    ATTR_INT_PAIR		= ATTR_TYPE(ATTR_BASE_INT, 		2),
    ATTR_INT_TRIPLE		= ATTR_TYPE(ATTR_BASE_INT, 		3),
    ATTR_LONG			= ATTR_TYPE(ATTR_BASE_LONG, 		1),
    ATTR_SHORT			= ATTR_TYPE(ATTR_BASE_SHORT, 		1),
    ATTR_ENUM			= ATTR_TYPE(ATTR_BASE_ENUM, 		1),
    ATTR_CHAR			= ATTR_TYPE(ATTR_BASE_CHAR, 		1),
#ifdef OW_I18N
    ATTR_WCHAR			= ATTR_TYPE(ATTR_BASE_WCHAR, 		1),
#define ATTR_CHAR_WC		ATTR_WCHAR
#endif
    ATTR_STRING			= ATTR_TYPE(ATTR_BASE_STRING, 		1),
#ifdef OW_I18N
    ATTR_WSTRING		= ATTR_TYPE(ATTR_BASE_WSTRING, 		1),
#define	ATTR_STRING_WCS		ATTR_WSTRING
#endif
    ATTR_FUNCTION_PTR		= ATTR_TYPE(ATTR_BASE_FUNCTION_PTR, 	1),
    ATTR_BOOLEAN		= ATTR_TYPE(ATTR_BASE_BOOLEAN, 		1),
    ATTR_OPAQUE			= ATTR_TYPE(ATTR_BASE_OPAQUE, 		1),
    ATTR_OPAQUE_PAIR		= ATTR_TYPE(ATTR_BASE_OPAQUE, 		2),
    ATTR_OPAQUE_TRIPLE		= ATTR_TYPE(ATTR_BASE_OPAQUE, 		3),
    ATTR_X			= ATTR_TYPE(ATTR_BASE_X, 		1),
    ATTR_INDEX_X		= ATTR_TYPE(ATTR_BASE_INDEX_X, 		2),
    ATTR_Y			= ATTR_TYPE(ATTR_BASE_Y, 		1),
    ATTR_INDEX_Y		= ATTR_TYPE(ATTR_BASE_INDEX_Y, 		2),
    ATTR_XY			= ATTR_TYPE(ATTR_BASE_XY, 		2),
    ATTR_INDEX_XY		= ATTR_TYPE(ATTR_BASE_INDEX_XY, 	3),
    ATTR_PIXRECT_PTR		= ATTR_TYPE(ATTR_BASE_PIXRECT_PTR, 	1),
    ATTR_PIXFONT_PTR		= ATTR_TYPE(ATTR_BASE_PIXFONT_PTR, 	1),
    ATTR_RECT_PTR		= ATTR_TYPE(ATTR_BASE_RECT_PTR, 	1),
    ATTR_AV			= ATTR_TYPE(ATTR_BASE_AV, 		1),
    ATTR_ICON_PTR		= ATTR_TYPE(ATTR_BASE_ICON_PTR, 	1),
    ATTR_SINGLE_COLOR_PTR	= ATTR_TYPE(ATTR_BASE_SINGLE_COLOR_PTR, 1),
    ATTR_CURSOR_PTR		= ATTR_TYPE(ATTR_BASE_CURSOR_PTR, 	1)
} Attr_base_cardinality;

/*
 *  The package id space potentially runs from 0 to 255 inclusive.  This is
 * subdivided as follows:
 *	     0	NEVER a valid package id.
 *	[1..32)	[ATTR_PKG_UNUSED_FIRST..ATTR_PKG_UNUSED_LAST]
 *		Reserved for non-Sun packages.
 *     [32..64)	Reserved for future use.
 *    [64..128)	[ATTR_PKG_FIRST..ATTR_PKG_LAST] is currently used by
 *		Sun packages
 *		(ATTR_PKG_LAST..128) is for future use by Sun packages.
 *   [128..256)	Reserved for future use.
 */
typedef enum {
    ATTR_PKG_ZERO		= ATTR_PKG_ZERO_VALUE,
    ATTR_PKG_GENERIC		= ATTR_PKG_GENERIC_VALUE,

    ATTR_PKG_CURSOR		= ATTR_PKG_FIRST +  1,
    ATTR_PKG_DRAWABLE		= ATTR_PKG_FIRST +  2,
    ATTR_PKG_FONT		= ATTR_PKG_FIRST +  3,
    ATTR_PKG_IMAGE		= ATTR_PKG_FIRST + 4,
    ATTR_PKG_SERVER_IMAGE	= ATTR_PKG_FIRST + 5,
    ATTR_PKG_SCREEN		= ATTR_PKG_FIRST + 6,
    ATTR_PKG_SELN_BASE		= 71,	/* ATTR_PKG_FIRST +  7 */
	/* ATTR_PKG_SELN_BASE must be 71, as it is known to 3.X and 4.X code.
	 * In fact, the layout of the bits in an attribute is known, and also
	 * cannot change without breaking communication between SunView 1 and
	 * XView selections.
	 */
    ATTR_PKG_SERVER		= ATTR_PKG_FIRST + 8,
    ATTR_PKG_WIN		= ATTR_PKG_FIRST + 9,
    ATTR_PKG_SV			= ATTR_PKG_FIRST + 10,
    ATTR_PKG_FULLSCREEN		= ATTR_PKG_FIRST + 11,
    ATTR_PKG_ERROR		= ATTR_PKG_FIRST + 12,
    ATTR_PKG_CMS		= ATTR_PKG_FIRST + 13,
    ATTR_PKG_DND		= ATTR_PKG_FIRST + 14,

    /* REMIND: ATTR_PKG_SELECTION should be ATTR_PKG_FIRST+15 put this
	       will cause the pkg values for the OL pkgs to change.  This
	       would break binary compatibility.  So we put the new intrinsic
	       pkg after the OL pkgs.  When we can again break binary 
	       compatibility, we should change this and add some space
	       between the intrinsic pkgs and the OL pkgs.
	       Remove comment in attrol.h when this is fixed.
     */
    ATTR_PKG_SELECTION		= ATTR_PKG_LAST_VALUE + 20,

    ATTR_PKG_LAST		= ATTR_PKG_LAST_VALUE
	/* 
	 * Change ATTR_PKG_LAST_VALUE to be EQUAL to the last legal pkg id.
	 * The procedure counter(), called by attr_make, aborts if 
	 * PKG_ID > ATTR_PKG_LAST 
	 * PKG name should also be added to attr_names[] in attr.c 
	 */
} Attr_pkg;

/*
 * Generic attributes: ATTR_PKG_GENERIC is shared with
 * generic.h [64..128).
 */
typedef enum {
    ATTR_LIST = GENERIC_ATTR(ATTR_LIST_PTR(ATTR_RECURSIVE, ATTR_NO_VALUE), 0),
    ATTR_NOP0 = GENERIC_ATTR(ATTR_NO_VALUE,			16),
    ATTR_NOP1 = GENERIC_ATTR(ATTR_OPAQUE,			17),
    ATTR_NOP2 = GENERIC_ATTR(ATTR_TYPE(ATTR_BASE_OPAQUE, 2),	18),
    ATTR_NOP3 = GENERIC_ATTR(ATTR_TYPE(ATTR_BASE_OPAQUE, 3),	19),
    ATTR_NOP4 = GENERIC_ATTR(ATTR_TYPE(ATTR_BASE_OPAQUE, 4),	20)
} Attr_generic;

/*
 * Structs 
 */

typedef union {
    struct {
	ENUM_BITFIELD(Attr_pkg) 		pkg             : 8;
	unsigned				ordinal		: 8;
        ENUM_BITFIELD(Attr_list_type)           list_type       : 2;
        ENUM_BITFIELD(Attr_list_ptr_type)       list_ptr_type   : 1;
	unsigned				consumed	: 1;
        ENUM_BITFIELD(Attr_base_type)           base_type       : 7;
	unsigned				spare1		: 1; /*unused*/
	unsigned				cardinality	: 4;
    } 			info;
    Attr_attribute	code;
} Attr_union;

/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

/*
 * Public Functions
 */
EXTERN_FUNCTION (Attr_avlist attr_create_list, (Attr_attribute attr1, DOTDOTDOT));

/*
 * Private Functions
 */
EXTERN_FUNCTION (int attr_copy, (Attr_avlist *source, Attr_avlist *dest));
EXTERN_FUNCTION (int attr_count, (Attr_avlist count));
EXTERN_FUNCTION (char * attr_name, (Attr_attribute attr));
EXTERN_FUNCTION (Attr_avlist attr_skip_value, (Attr_attribute attr, Attr_avlist avlist));
EXTERN_FUNCTION (int attr_rc_unit_to_x, (unsigned int encoded_value, int col_width,  int left_margin, int col_gap));
EXTERN_FUNCTION (int attr_rc_unit_to_y, (unsigned int encoded_value, int row_height, int top_margin,  int row_gap));
EXTERN_FUNCTION (void attr_rc_units_to_pixels, (Attr_avlist avlist, int col_width, int row_height, int left_margin, int top_margin, int col_gap, int row_gap));

#endif /* xview_attr_DEFINED */
