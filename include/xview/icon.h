/*	@(#)icon.h 20.13 90/03/13 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_icon_DEFINED
#define xview_icon_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/window.h>
#include <xview/attrol.h>
#include <xview/svrimage.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define ICON 			&xv_icon_pkg
#define ICON_TYPE  		ATTR_PKG_ICON

#define ICON_WIDTH		XV_WIDTH
#define ICON_HEIGHT		XV_HEIGHT
#define ICON_LABEL		XV_LABEL
#ifdef OW_I18N
#define ICON_LABEL_WCS          XV_LABEL_WCS
#endif
#define ICON_FONT		XV_FONT

#define icon_attr_next(attr) 	(Icon_attribute *)attr_next((Xv_opaque *)attr)

#define	ICON_DEFAULT_WIDTH	64
#define	ICON_DEFAULT_HEIGHT	64


/*
 * PRIVATE #defines 
 */

#define ICON_ATTR(type, ordinal)       ATTR(ATTR_PKG_ICON, type, ordinal)
#define ICON_ATTR_LIST(ltype, type, ordinal) \
        ICON_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structures
 ***********************************************************************
 */

typedef Xv_opaque      Icon;

typedef enum {
   ICON_IMAGE		  = ICON_ATTR(ATTR_PIXRECT_PTR,    5),
   ICON_IMAGE_RECT	  = ICON_ATTR(ATTR_RECT_PTR, 	  10),
   ICON_LABEL_RECT  	  = ICON_ATTR(ATTR_RECT_PTR, 	  15),
   ICON_TRANSPARENT       = ICON_ATTR(ATTR_BOOLEAN,       20),
   ICON_MASK_IMAGE	  = ICON_ATTR(ATTR_PIXRECT_PTR,   25),
   ICON_TRANSPARENT_LABEL = ICON_ATTR(ATTR_STRING,        30)
#ifdef OW_I18N
                                                             ,
   ICON_TRANSPARENT_LABEL_WCS
                          = ICON_ATTR(ATTR_WSTRING,        35)
#endif /* OW_I18N */
} Icon_attribute;

typedef struct {
	Xv_window_struct	 parent_data;
	Xv_opaque		 private_data;
}   Xv_icon;

/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

/*
 * Variables 
 */
extern Xv_pkg		xv_icon_pkg;

/*
 * PUBLIC Functions 
 * For Sunview 1 Compatibility
 */

EXTERN_FUNCTION (Icon icon_create, (Attr_attribute attr1, DOTDOTDOT));
EXTERN_FUNCTION (Xv_opaque icon_get, (Icon icon_public, Icon_attribute attr));
EXTERN_FUNCTION (int icon_set, (Icon icon, DOTDOTDOT));
EXTERN_FUNCTION (int icon_destroy, (Icon icon));
EXTERN_FUNCTION (void icon_display, (Icon icon, int x, int y));

#endif /* ~xview_icon_DEFINED */
