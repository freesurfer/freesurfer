/*	@(#)svrimage.h 20.23 91/09/14 SMI	*/

#ifndef xview_server_image_DEFINED
#define xview_server_image_DEFINED

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 ***********************************************************************
 *				Include Files
 ***********************************************************************
 */

#include <xview/generic.h>
#include <xview/pkg.h>
#include <xview/drawable.h>
#include <sys/types.h>				/* needed for pixrect.h */
#include <pixrect/pixrect.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define	SERVER_IMAGE 			&xv_server_image_pkg

#define	server_image_begin_bits(name)	unsigned short name[] = {
#define	server_image_end_bits(name)	};

/*
 * PRIVATE #defines 
 */

#define SERVER_IMAGE_ATTR(type, ordinal)	\
				ATTR(ATTR_PKG_SERVER_IMAGE, type, ordinal)

#define SERVER_IMAGE_TYPE  	ATTR_PKG_SERVER_IMAGE

#define SERVER_IMAGE_ATTR_LIST(ltype, type, ordinal) \
        SERVER_IMAGE_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

typedef Xv_opaque	 	Server_image;

typedef enum {
	/*
	 * Public attributes 
	 */
	SERVER_IMAGE_DEPTH	= 	SERVER_IMAGE_ATTR(ATTR_INT, 	1),
	SERVER_IMAGE_BITS	= 	SERVER_IMAGE_ATTR(ATTR_OPAQUE, 	2),
	SERVER_IMAGE_X_BITS	=	SERVER_IMAGE_ATTR(ATTR_OPAQUE, 	3),
	SERVER_IMAGE_COLORMAP   = 	SERVER_IMAGE_ATTR(ATTR_STRING,  4),
	SERVER_IMAGE_BITMAP_FILE=       SERVER_IMAGE_ATTR(ATTR_STRING,  5),
	SERVER_IMAGE_PIXMAP     =       SERVER_IMAGE_ATTR(ATTR_OPAQUE,  6),
	SERVER_IMAGE_SAVE_PIXMAP=       SERVER_IMAGE_ATTR(ATTR_BOOLEAN, 7),
#ifdef  OW_I18N
        SERVER_IMAGE_BITMAP_FILE_WCS=   SERVER_IMAGE_ATTR(ATTR_STRING,  8),
#endif
	SERVER_IMAGE_CMS	=	SERVER_IMAGE_ATTR(ATTR_OPAQUE,  9)
} Server_image_attribute;

typedef struct {
	Xv_drawable_struct	parent_data;
	Xv_opaque		private_data;
	Xv_embedding		embedding_data;
	Pixrect			pixrect;
}   Xv_server_image;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

Xv_public Xv_pkg  xv_server_image_pkg;

/*
 * PUBLIC functions 
 */
EXTERN_FUNCTION (int server_image_rop, (Xv_opaque dest, int dx,	int dy, int dw, int dh, unsigned long op, Xv_opaque src, int sx, int sy));
EXTERN_FUNCTION (int 	server_image_stencil, (Xv_opaque dest, int dx, int dy, int dw, int dh, unsigned long op, Xv_opaque st, int stx, int dty, Xv_opaque src, int sx, int sy));
EXTERN_FUNCTION (int 	server_image_destroy, (Pixrect *pr));
EXTERN_FUNCTION (int 	server_image_get, (Xv_opaque dest, int x, int y));
EXTERN_FUNCTION (int 	server_image_put, (Xv_opaque dest, int x, int y, int value));
EXTERN_FUNCTION (int 	server_image_vector, (Xv_opaque dest, int x0, int y0, int x1, int y1, int op, int value));
EXTERN_FUNCTION (Pixrect *server_image_region, (Xv_opaque dest, int x, int y, int w, int h));
EXTERN_FUNCTION (int 	server_image_colormap, (Xv_opaque dest, int index, int count, unsigned char *red, unsigned char *green, unsigned char *blue));

#endif /* xview_server_image_DEFINED */
