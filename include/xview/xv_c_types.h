/*	@(#)xv_c_types.h 1.8 91/09/14 SMI      */

/*
 *      (c) Copyright 1989 Sun Microsystems, Inc. Sun design patents
 *      pending in the U.S. and foreign countries. See LEGAL NOTICE
 *      file for terms of the license.
 */
#ifndef XV_C_TYPES_H
#define XV_C_TYPES_H

#ifndef C_VARIETIES_H
#define C_VARIETIES_H

/*
 *	This file defines some macros that are used to make code
 *	portable between the major C dialects currently in use at
 *	Sun.  As of 12/88, these include Sun C (a lot like K&R C),
 *	ANSI C, and C++.
 *
 * external functions:
 *	To declare an external function, invoke the EXTERN_FUNCTION
 *	macro; the macro's first parameter should be the function's
 *	return type and function name, and the second macro parameter
 *	should be the parenthesized list of function arguments (or an
 *	ellipsis - DOTDOTDOT macro should be used to indicate the 
 *      ellipsis as explained later in this file - if the arguments are 
 *      unspecified or of varying number or type, or the uppercase word 
 *      "_VOID_" if the function takes no arguments).  Some examples:
 *
 *	    EXTERN_FUNCTION( void printf, (char *, DOTDOTDOT) );
 *	    EXTERN_FUNCTION( int fread, (char*, int, int, FILE*) );
 *	    EXTERN_FUNCTION( int getpid, (_VOID_) );
 *	    EXTERN_FUNCTION( int atoi, (CONST char *str) );
 *
 *	Note that to be ANSI-C conformant, one should put "," at the end
 *	first argument of printf() declaration.
 *
 * structure tags:
 *	In order to handle cases where a structure tag has the same name
 *	as a type, the STRUCT_TAG macro makes the tag disappear in C++.
 *	An example (from <sys/types.h>):
 *
 *	    typedef struct STRUCT_TAG(fd_set) { ... } fd_set;
 *
 * enum bitfields:
 *	In K&R C as interpreted at UCB, bitfields may be declared to
 *	be of an enumerated type.  Neither ANSI C nor C++ permit this,
 *	so the ENUM_BITFIELD macro replaces the enum declaration with
 *	"unsigned".   An example (from <sunwindow/attr.h>):
 *
 *	    struct {
 *		ENUM_BITFIELD( Attr_pkg )	pkg		: 8;
 *		unsigned			ordinal		: 8;
 *		ENUM_BITFIELD( Attr_list_type )	list_type	: 8;
 *		...
 *	    };
 *
 * enum type specifier:
 *	In K&R C, it is OK to use "enum xyz" as return type but in C++,
 * 	one should use "xyz".  ENUM_TYPE macro is used to handle this.
 *
 *		ENUM_TYPE(enum, xyz) (*func) (...);
 *
 * "struct s s;":
 *	C++ does not allow this sort of name conflict, since struct tags are
 *	in the same namespace as variables.  In this case, we use the
 *	NAME_CONFLICT macro to prepend (for C++) an underscore to the
 *	variable (or struct member) name.  E.g. from <pixrect/pixrect.h>:
 *
 *	    typedef struct pixrect {
 *		struct	pixrectops *pr_ops;
 *		struct	pr_size NAME_CONFLICT(pr_size);
 *	    } Pixrect;
 *	    #define pr_height	NAME_CONFLICT(pr_size).y
 *	    #define pr_width	NAME_CONFLICT(pr_size).x
 *
 *	Note that no spaces are allowed within the parentheses in the
 *	invocation of NAME_CONFLICT.
 *
 * Pointers to functions declared as struct members:
 *	Instead of getting picky about the types expected by struct
 *	members which are pointers to functions, we use DOTDOTDOT to
 *	tell C++ not to be so uptight:
 *
 *	    struct pixrectops {
 *		    int	(*pro_rop)( DOTDOTDOT );
 *		    int	(*pro_stencil)( DOTDOTDOT );
 *		    int	(*pro_batchrop)( DOTDOTDOT );
 *		    . . .
 *	    };
 *
 */


/* Which type of C/C++ compiler are we using? */

#if defined(__cplusplus)
    /*
     * Definitions for C++ 2.0 and later require extern "C" { decl; }
     */
#   define EXTERN_FUNCTION( rtn, args ) extern "C" { rtn args; }
#   define STRUCT_TAG( tag_name ) /* the tag disappears */
#   define ENUM_BITFIELD( enum_type ) unsigned
#   define ENUM_TYPE( enum_sp, enum_ty ) enum_ty

#if defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus)
#   define NAME_CONFLICT( name ) _##name
#else
#   define NAME_CONFLICT( name ) _/**/name
#endif

#   define DOTDOTDOT ...
#   define _VOID_ /* anachronism */
#   define CONST const

/*
 * This is not necessary for 2.0 since 2.0 has corrected the void (*) () problem
 */
typedef void (*_PFV_)();
typedef int (*_PFI_)();

#elif defined(c_plusplus)
    /*
     * Definitions for C++ 1.2
     */
#   define EXTERN_FUNCTION( rtn, args ) rtn args
#   define STRUCT_TAG( tag_name )  /* the tag disappears */
#   define ENUM_BITFIELD( enum_type ) unsigned
#   define ENUM_TYPE( enum_sp, enum_ty ) enum_ty
#   define NAME_CONFLICT( name ) _/**/name
#   define DOTDOTDOT ...
#   define _VOID_ /* anachronism */
#   define CONST const 

typedef void (*_PFV_)();
typedef int (*_PFI_)();

#elif defined(__STDC__)
    /*
     * Definitions for ANSI C
     */
#   define EXTERN_FUNCTION( rtn, args ) rtn args
#   define STRUCT_TAG( tag_name ) tag_name
#   define ENUM_BITFIELD( enum_type ) unsigned
#   define ENUM_TYPE( enum_sp, enum_ty ) enum_sp enum_ty
#   define NAME_CONFLICT( name ) name
#   define DOTDOTDOT ...
#   define _VOID_ void
#   define CONST   

#else
    /*
     * Definitions for Sun/K&R C -- ignore function prototypes,
     * but preserve tag names and enum bitfield declarations.
     */
#   define EXTERN_FUNCTION( rtn, args ) rtn()
#   define STRUCT_TAG( tag_name ) tag_name
#   define ENUM_BITFIELD( enum_type ) enum_type
#   define ENUM_TYPE( enum_sp, enum_ty ) enum_sp enum_ty
#   define NAME_CONFLICT( name ) name
#   define DOTDOTDOT
#   define _VOID_
    /* VOID is only used where it disappears anyway */
#   define CONST   

#endif /* Which type of C/C++ compiler are we using? */

#endif /* ~defined(C_VARIETIES_H) */

#endif /* defined(XV_C_TYPES_H) */

