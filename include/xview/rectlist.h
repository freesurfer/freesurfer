/*	@(#)rectlist.h 20.15 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_rectlist_DEFINED
#define xview_rectlist_DEFINED

#include <xview/xv_c_types.h>

/*
 * Defines the interface to the data structure called 
 * a rectlist which is a list of rectangles.
 *
 * A rectlist has an offset (rl_x, rl_y) assoicated with it that is used
 * to efficiently pass this data structure from one embedded coordinate
 * system to another without changing the offsets of all the rectangles in
 * the list.
 *
 * Rl_bound is the rectangle that bounds all the rectangles in the list
 * and is maintained to facilitate efficient rectangle algebra.
 */


/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define	RECTNODE_NULL	((Rectnode *)0)
#define	RECTLIST_NULL	((Rectlist *)0)

/*
 * Rectlist Transformation macros used for passing rectlists up/down embedded
 * coordinate systems.
 */
#define	rl_passtoparent(x,y,rl) \
	{(rl)->rl_x=(rl)->rl_x+(x); (rl)->rl_y=(rl)->rl_y+(y);}

#define	rl_passtochild(x,y,rl) \
	{(rl)->rl_x=(rl)->rl_x-(x); (rl)->rl_y=(rl)->rl_y-(y);}

/*
 * Rectlist Offset Adjustment macros
 */
#define	rl_rectoffset(rl,r1,r) \
	{*(r)= *(r1); (r)->r_left+=(rl)->rl_x; (r)->r_top+=(rl)->rl_y;}

#define	rl_coordoffset(rl,x,y) {*(x)-=(rl)->rl_x;*(y)-=(rl)->rl_y;}


/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

/*
 * PUBLIC structures 
 */

typedef	struct	rectnode {
	struct	rectnode *rn_next;	/* Pointer to next rectnode */
	struct	rect rn_rect;
} Rectnode;


typedef	struct	rectlist {
	coord	rl_x, rl_y;		/* Offset to apply to each rect
					   in list including bound */
	struct	rectnode *rl_head;	/* Pointer to first rectnode */
	struct	rectnode *rl_tail;	/* Pointer to last rectnode */
	struct	rect rl_bound;		/* Describes bounding rect of 
					   all rects in list */
} Rectlist;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC functions 
 */
extern  struct rectlist rl_null;

EXTERN_FUNCTION (unsigned rl_empty, (Rectlist *rl));
EXTERN_FUNCTION (unsigned rl_equal, (Rectlist *rl1, Rectlist *rl2));
EXTERN_FUNCTION (unsigned rl_boundintersectsrect, (Rect *r, Rectlist *rl));
EXTERN_FUNCTION (unsigned rl_rectintersects, (Rect *r, Rectlist *rl));
EXTERN_FUNCTION (unsigned rl_equalrect, (Rect *r, Rectlist *rl));
EXTERN_FUNCTION (unsigned rl_includespoint, (Rectlist *r, int x, int y));

#ifdef xview_other_rl_funcs

EXTERN_FUNCTION (void rl_rectintersection, (Rect *r, Rectlist *rl1, Rectlist *rl));
EXTERN_FUNCTION (void rl_rectunion, (Rect *r, Rectlist *rl1, Rectlist *rl));
EXTERN_FUNCTION (void rl_rectdifference, (Rect *r, Rectlist *rl1, Rectlist *rl));
EXTERN_FUNCTION (void rl_intersection, (Rectlist *rl1, Rectlist *rl2, Rectlist *rl));
EXTERN_FUNCTION (void rl_sort, (Rectlist *rl1, Rectlist *rl, int sortorder));
EXTERN_FUNCTION (void rl_union, (Rectlist *rl1, Rectlist *rl2, Rectlist *rl));
EXTERN_FUNCTION (void rl_difference, (Rectlist *rl1, Rectlist *rl2, Rectlist *rl));
EXTERN_FUNCTION (void rl_initwithrect, (Rect *r, Rectlist *rl));
EXTERN_FUNCTION (void rl_copy, (Rectlist *rl1, Rectlist *rl));
EXTERN_FUNCTION (void rl_free, (Rectlist *rl));
EXTERN_FUNCTION (void rl_coalesce, (Rectlist *rl));
EXTERN_FUNCTION (void rl_normalize, (Rectlist *rl));
EXTERN_FUNCTION (void rl_print, (Rectlist *rl, char *tag));

#endif /* xview_other_rl_funcs */

#endif	/* ~xview_rectlist_DEFINED */
