/*	@(#)cms_grays.h 20.11 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_cms_grays_DEFINED
#define xview_cms_grays_DEFINED

/*
 * Definition of the colormap segment CMS_GRAYS, a collection of grays.
 */

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */
#define	CMS_GRAYS	"grays"
#define	CMS_GRAYSSIZE	8

#define	WHITE	0
#define	GRAY(i)	((BLACK)*(i)/100)
#define	BLACK	7

#define	cms_grayssetup(r,g,b) 						\
	{ unsigned char v = 0, vi; 					\
	  int	i, gi; 							\
	  vi = 255/BLACK; 						\
	  for (i=BLACK-1;i>WHITE;i--) { /* Dark (small v)->light (large v) */ \
		v += vi; 						\
		gi = GRAY(100/8*i); 					\
		(r)[gi] = v; (g)[gi] = v; (b)[gi] = v;  		\
	  } 								\
	  (r)[WHITE] = 255;	(g)[WHITE] = 255;	(b)[WHITE] = 255; \
	  (r)[BLACK] = 0;	(g)[BLACK] = 0;		(b)[BLACK] = 0; \
	}

#endif /* ~xview_cms_grays_DEFINED */
