/*	@(#)rect.h 20.19 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 * Defines the interface to the geometric object
 * called a Rect which is a rectangle.
 */

#ifndef xview_rect_DEFINED
#define xview_rect_DEFINED

#include <xview/xv_c_types.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines
 */

#define	RECT_NULL	((Rect *)0)

#ifndef coord
#define coord	short
#endif

/*
 * Rectangle sort ordering.
 */
#define	RECTS_TOPTOBOTTOM	0
#define	RECTS_BOTTOMTOTOP	1
#define	RECTS_LEFTTORIGHT	2
#define	RECTS_RIGHTTOLEFT	3

#define	RECTS_UNSORTED		-1
#define	RECTS_SORTS		4

/*
 * Rect Geometry macros
 */
#define	rect_right(rect)  	((rect)->r_left+(rect)->r_width-1)
#define	rect_bottom(rect) 	((rect)->r_top+(rect)->r_height-1)

#define rect_print(rect) 							\
        (void)fprintf(stderr,"[left: %d, top: %d, width: %d, height: %d]\n",	\
            (rect)->r_left, (rect)->r_top, (rect)->r_width, (rect)->r_height)

#define rect_marginadjust(r,m) 				\
	{ (r)->r_left-=(m);(r)->r_top-=(m); 		\
	 (r)->r_width+=(m)+(m);(r)->r_height+=(m)+(m);}

#define rect_borderadjust(r,m) 				\
	{ (r)->r_width+=(m)+(m);(r)->r_height+=(m)+(m);}

#define rect_construct(r,x,y,w,h) \
	{(r)->r_left=(x);(r)->r_top=(y);(r)->r_width=(w);(r)->r_height=(h);}

/*
 * Rect Predicate macros
 */
#define rect_equal(r1,r2) \
	((r1)->r_left==(r2)->r_left && (r1)->r_width==(r2)->r_width && \
	 (r1)->r_top==(r2)->r_top && (r1)->r_height==(r2)->r_height)

#define rect_sizes_differ(r1, r2) \
        ((r1)->r_width != (r2)->r_width || (r1)->r_height != (r2)->r_height)

#define rect_isnull(r)		((r)->r_width == 0 || (r)->r_height == 0)

#define rect_includespoint(r,x,y) \
	((x) >= (r)->r_left && (y) >= (r)->r_top && \
	 (x)<(r)->r_left+(r)->r_width && (y)<(r)->r_top+(r)->r_height)

#define rect_includesrect(r1, r2) \
	((r1)->r_left <= (r2)->r_left && (r1)->r_top <= (r2)->r_top && \
	 (r1)->r_left+(r1)->r_width >= (r2)->r_left+(r2)->r_width && \
	 (r1)->r_top+(r1)->r_height >= (r2)->r_top+(r2)->r_height)

#define rect_intersectsrect(r1,r2) \
	((r1)->r_left<(r2)->r_left+(r2)->r_width && \
	 (r1)->r_top<(r2)->r_top+(r2)->r_height &&  \
	 (r2)->r_left<(r1)->r_left+(r1)->r_width && \
	 (r2)->r_top<(r1)->r_top+(r1)->r_height)

/*
 * Rect Transformation macros used for passing rects up/down embedded
 * coordinate systems.
 */
#define	rect_passtoparent(x,y,rect) \
	{(rect)->r_left=(rect)->r_left+(x); (rect)->r_top=(rect)->r_top+(y);}

#define	rect_passtochild(x,y,rect) \
	{(rect)->r_left=(rect)->r_left-(x); (rect)->r_top=(rect)->r_top-(y);}

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

/*
 * PUBLIC structures 
 */

typedef struct rect {
	coord	r_left, r_top;
	short	r_width, r_height;
} Rect;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC variables 
 */
extern	struct rect 	rect_null;

/*
 * PUBLIC Functions
 */

EXTERN_FUNCTION (struct rect rect_bounding, (Rect *r1, Rect *r2));
EXTERN_FUNCTION (unsigned rect_clipvector, (Rect *rect, int *x1arg, int *y1arg, int *x2arg, int *y2arg));
EXTERN_FUNCTION (unsigned rect_order, (Rect *rl, Rect *r2, int sortorder));
EXTERN_FUNCTION (int rect_right_of, (Rect *rect1, Rect *rect2));
EXTERN_FUNCTION (int rect_below, (Rect *rect1, Rect *rect2));
#ifdef _OTHER_RECT_FUNCTIONS

EXTERN_FUNCTION (void rect_intersection, (Rect *rl, Rect *r2, Rect *r));
EXTERN_FUNCTION (int rect_distance, (Rect *rect, int x, int y, int * x_used, int * y_used));
#endif /* _OTHER_RECT_FUNCTIONS */

#endif /* ~xview_rect_DEFINED */
