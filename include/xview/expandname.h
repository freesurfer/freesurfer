/*	@(#)expandname.h 20.13 91/09/14 SMI	*/
 
/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_expand_name_DEFINED
#define xview_expand_name_DEFINED

#include <xview/xv_c_types.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */

#define NONAMES	((struct namelist*) 0)

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

struct namelist	{
	unsigned	 count;
	char		*names[1];	/* LIE --it's actually names[count]
					 * followed by the strings themselves
					 */
};

/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

/*
 * Public Functions
 *
 *	struct namelist *
 *	xv_expand_name(str)
 *	    char	*str;
 *
 *	void
 *	free_namelist(nl)
 *	    struct namelist *nl;
 *
 *  xv_expand_name returns a pointer to a struct namelist containing
 *	the client's shell's idea of what its argument expands to.
 *	If str contains no shell metacharacters, the shell is never
 *	consulted, and the argument string is simply returned
 *	in a singleton namelist.
 * 
 *	In case of errors, xv_expand_name() writes an error message to
 *	stderr, and returns NONAMES.
 *
 *  free_namelist
 *	The struct namelist is dynamically allocated, and should be
 *	freed after use by a call to free_namelist().
 */

EXTERN_FUNCTION (struct 	namelist *xv_expand_name, (char *str));
EXTERN_FUNCTION (void 		free_namelist, (struct namelist *nl));

#endif	/* ~xview_expand_name_DEFINED */
