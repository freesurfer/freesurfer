/*	@(#)win_struct.h 20.11 91/09/14 SMI	*/

#ifndef xview_win_struct_DEFINED
#define xview_win_struct_DEFINED

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

/* Host name maximum size is 31 + 1 for colon + 2 for display number
 * + 10 for xid + 1 for null.
 */
#define	WIN_NAMESIZE		45

/*
 * Link names.
 */
#define	WL_PARENT		0
#define	WL_OLDERSIB		1
#define	WL_YOUNGERSIB		2
#define	WL_OLDESTCHILD		3
#define	WL_YOUNGESTCHILD	4

#define	WL_ENCLOSING		WL_PARENT
#define	WL_COVERED		WL_OLDERSIB
#define	WL_COVERING		WL_YOUNGERSIB
#define	WL_BOTTOMCHILD		WL_OLDESTCHILD
#define	WL_TOPCHILD		WL_YOUNGESTCHILD

#define	WIN_NULLLINK		-1

#endif /* ~xview_win_struct_DEFINED */
