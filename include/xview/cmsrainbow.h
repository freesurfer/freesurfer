/*	@(#)cmsrainbow.h 20.11 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 * Definition of the colormap segment CMS_RAINBOW,
 * a small collection of colors of the rainbow.
 */

#ifndef xview_cms_rainbow_DEFINED
#define xview_cms_rainbow_DEFINED

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */
#define	CMS_RAINBOW		"rainbow"
#define	CMS_RAINBOWSIZE	8

#define	WHITE			0
#define	RED			1
#define	ORANGE			2
#define	YELLOW			3
#define	GREEN			4
#define	BLUE			5
#define	INDIGO			6
#define	VIOLET			7

#define	cms_rainbowsetup(r,g,b) 					\
	(r)[WHITE] = 255;	(g)[WHITE] = 255;	(b)[WHITE] = 255; \
	(r)[RED] = 255;		(g)[RED] = 0;		(b)[RED] = 0; 	\
	(r)[ORANGE] = 192;	(g)[ORANGE] = 64;	(b)[ORANGE] = 0; \
	(r)[YELLOW] = 128;	(g)[YELLOW] = 128;	(b)[YELLOW] = 0; \
	(r)[GREEN] = 0;		(g)[GREEN] = 255;	(b)[GREEN] = 0; \
	(r)[BLUE] = 0;		(g)[BLUE] = 0;		(b)[BLUE] = 255; \
	/* \
	 * The rule for indigo is B > R & SQRT(B**2+R**2) < .5 \
	 * where 0.0<=B|R<=1.0).  Trying R=.25 and B=.3. \
	 */ \
	(r)[INDIGO] = 64;	(g)[INDIGO] = 0;	(b)[INDIGO] = 76; \
	/* \
	 * The rule for violet is R > B & SQRT(B**2+R**2) > .5 \
	 * where 0.0<=B|R<=1.0).  Trying R=.5 and B=.7. \
	 */ \
	(r)[VIOLET] = 128;	(g)[VIOLET] = 0;	(b)[VIOLET] = 178;


#endif	/* ~xview_cms_rainbow_DEFINED */
