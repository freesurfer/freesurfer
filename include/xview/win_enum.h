/*      @(#)win_enum.h 20.12 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifdef TEST_WIN_ENUM_DEPENDENCIES
#ifndef xview_win_enum_DEFINED
#define xview_win_enum_DEFINED	1

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/rect.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

/*
 * Values for 'flags' field in struct Win_enum_node 
 */
#define WIN_NODE_INSERTED     0x1
#define WIN_NODE_OPEN         0x2
#define WIN_NODE_IS_ROOT      0x4

/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 */

typedef int	  Window_handle;

typedef enum win_enumerator_result { 
	Enum_Normal, 
	Enum_Succeed, 
	Enum_Fail 
} Win_enum_result, (*Enumerator)();

/*
 * For a fast window enumerator in user-space 
 */

typedef struct win_enum_node	{
	unsigned char	me;
	unsigned char	parent;
	unsigned char	upper_sib;
	unsigned char	lowest_kid;
	unsigned int	flags;
	Rect		open_rect;
	Rect		icon_rect;
}   Win_enum_node;

typedef struct win_tree_layer	{
	unsigned int	bytecount;
	Win_enum_node  *buffer;
}   Win_tree_layer;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC functions 
 */


#endif /* ~xview_win_enum_DEFINED */
#endif /* TEST_WIN_ENUM_DEPENDENCIES */
