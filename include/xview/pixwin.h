/*	@(#)pixwin.h 20.25 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_pixwin_DEFINED
#define	xview_pixwin_DEFINED

/*
 ***********************************************************************
 *                      Include Files
 ***********************************************************************
 */
#include <pixrect/pixrect.h>
#include <xview/rect.h>
#include <xview/rectlist.h>
#include <xview/base.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */
#define PIX_MAX_PLANE_GROUPS  	1
#define	PW_PIXEL_CACHE_NULL	((Pw_pixel_cache *)0)

#define XV_DEFAULT_FG_BG        0
#define XV_INVERTED_FG_BG       1
#define XV_USE_OP_FG            0 
#define XV_USE_CMS_FG           1

#define PW_ROP_NULL_SRC 0
#define PW_VECTOR       1
#define PW_LINE         2
#define PW_POLYLINE     3
#define PW_TEXT         4
#define PW_STENCIL      5
#define PW_POLYGON2     6
#define PW_POLYPOINT    7
#define PW_ROP          8
#define PW_REPLROP      9
#define PW_NUM_OPS      10

/*
 *    SV1 compatibility definitions.
 */
#define pw_copy(dpw, dx, dy, w, h, op, spw, sx, sy) \
		xv_rop(dpw, dx, dy, w, h, op, (Pixrect *)spw, sx, sy)

#define pw_read(dpr, dx, dy, w, h, op, spw, sx, sy) \
		xv_read(dpr, dx, dy, w, h, op, spw, sx, sy)

#define pw_replrop(dpw, xw, yw, width, height, op, pr, xr, yr) \
		xv_replrop(dpw, xw, yw, width, height, op, pr, xr, yr)

#define pw_rop(dpw, dx, dy, w, h, op, sp, sx, sy) \
	    	xv_rop(dpw, dx, dy, w, h, op, (Pixrect *)sp, sx, sy)

#define pw_setcmsname(pw, name) xv_set(pw, WIN_CMS_NAME, name, 0)

#define pw_stencil(dpw, x, y, w, h, op, stpr, stx, sty, spr, sy, sx) \
		xv_stencil(dpw, x, y, w, h, op, stpr, stx, sty, spr, sy, sx)

#define pw_text(window, xbasew, ybasew, op, pixfont, str) \
		xv_text(window, xbasew, ybasew, op, pixfont, str)

#define pw_ttext(pw, xbasew, ybasew, op, pixfont, str) \
		xv_ttext(pw, xbasew, ybasew, op, pixfont, str)

#define pw_vector(window, x0, y0, x1, y1, op, cms_index) \
		xv_vector(window, x0, y0, x1, y1, op, cms_index)

#define pw_write(dpw, dx, dy, w, h, op, spr, sx, sy) \
	   	xv_rop(dpw, dx, dy, w, h, op, (Pixrect *)spr, sx, sy)

#define pw_writebackground(dpw, dx, dy, w, h, op) \
		xv_rop(dpw, dx, dy, w, h, op, (Pixrect *)NULL, 0, 0)

/* 
 *	Obsolete SV1 pixwin functions.
 */
#define pw_batch(pw, type)
#define pw_batch_off(pw)
#define pw_batch_on(pw)
#define pw_close(pw)
#define pw_dbl_access(pw)
#define pw_dbl_flip(pw)
#define pw_dbl_release(pw)
#define pw_dbl_get(pw, attrs)
#define pw_dbl_set(pw, attrs)
#define pw_get_region_rect(pw, r)
#define pw_get_x_offset(pw) 
#define pw_get_y_offset(pw) 
#define pw_lock(pixwin, rect)
#define pw_region(pw, x, y, w, h)
#define pw_set_region_rect(pw, r, use_same_pr)
#define pw_set_xy_offset(pw, x_offset, y_offset)
#define pw_set_x_offset(pw, x_offset)
#define pw_set_y_offset(pw, y_offset)
#define pw_show(pw)
#define pw_unlock(pixwin)
#define pw_use_fast_monochrome(pw)


/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */
typedef struct pw_pixel_cache {
	Rect r;
	struct pixrect * plane_group[PIX_MAX_PLANE_GROUPS];
} Pw_pixel_cache;

typedef	struct	pixwin {
	char	dummy;		/* dummy field for compatibility */
} Pixwin;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC functions 
 */

EXTERN_FUNCTION (Pw_pixel_cache * pw_save_pixels, (Xv_opaque pw, Rect *rect));
EXTERN_FUNCTION (void pw_restore_pixels, (Xv_opaque ps, Pw_pixel_cache *pc));
 
/*
 * For SunView Pixrect/Pixwin Graphics Compatibility
 */

EXTERN_FUNCTION (int xv_read, (Pixrect *pr, int op, int x, int y, int width, int height, Xv_opaque window, int sx, int sy));
EXTERN_FUNCTION (int xv_replrop, (Xv_opaque window, int op, int xw, int yw, int width, int height, Pixrect *pr, int xr, int yr)); 
EXTERN_FUNCTION (int xv_rop, (Xv_opaque window, int op, int x, int y, int width, int height, Pixrect *pr, int xr, int yr));
EXTERN_FUNCTION (int xv_stencil, (Xv_opaque window, int op, int dx, int dy, int width, int height, Pixrect *stpr, int stx, int sty, Pixrect *spr, int sx, int sy));
EXTERN_FUNCTION (int xv_text, (Xv_opaque window, int op, int xbasew, int ybasew, Xv_opaque font, char *str));
EXTERN_FUNCTION (int xv_ttext, (Xv_opaque window, int xbasew, int ybasew, int op, Xv_opaque font, char *str));
EXTERN_FUNCTION (int xv_vector, (Xv_opaque window, int x0, int y0, int x1, int y1, int op, int cms_index));
EXTERN_FUNCTION (int pw_batchrop, (Pixwin *pw, int x, int y, int op, struct pr_prpos *sbp, int m));
EXTERN_FUNCTION (int pw_get, (Xv_opaque xv_drawable, int x, int y));
EXTERN_FUNCTION (int pw_put, (Xv_opaque pw, int x, int y, int val));
EXTERN_FUNCTION (int pw_putattributes, (Xv_opaque pw, int *planes));
EXTERN_FUNCTION (int pw_getattributes, (Xv_opaque pw, int *planes));
EXTERN_FUNCTION (int pw_putcolormap, (Xv_opaque pw, int index, int count, unsigned char *red, unsigned char *green, unsigned char *blue));
EXTERN_FUNCTION (int pw_getcolormap, (Xv_opaque pw, int index, int count, unsigned char *red, unsigned char *green, unsigned char *blue));
EXTERN_FUNCTION (struct pixfont * pw_pfsysopen, (void));
EXTERN_FUNCTION (int pw_pfsysclose, (void));
 
#endif	/* ~xview_pixwin_DEFINED */
