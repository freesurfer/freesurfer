/*      @(#)pkg_public.h 20.24 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_pkg_public_DEFINED
#define xview_pkg_public_DEFINED


/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/pkg.h>
#include <xview/xv_error.h>


/*
 ***********************************************************************
 *		Typedefs, Enumerations, and Structures
 ***********************************************************************
 *
 * SunView pkg. definition	
 */

/*
 * PRIVATE structures for pkg implementors only  
 */

/*
 * Last field before "embedded" struct in an "embedding object". 
 */
typedef long unsigned	 Xv_embedding;


/*
 * Base instance for all objects	
 */
typedef struct {
    long unsigned	 seal;	/* Has "special" value meaning "am object" */
    Xv_pkg		*pkg;   /* Always points to pkg chain for an object */
} Xv_base;

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC General interface functions	
 */
EXTERN_FUNCTION (Xv_object xv_create, (Xv_opaque owner, Xv_pkg *pkg, DOTDOTDOT));
EXTERN_FUNCTION (Xv_object xv_find, (Xv_opaque owner, Xv_pkg *pkg, DOTDOTDOT));
EXTERN_FUNCTION (Xv_opaque xv_set, (Xv_opaque object, DOTDOTDOT));
EXTERN_FUNCTION (Xv_opaque xv_get, (Xv_opaque object, Attr_attribute attr, DOTDOTDOT));
EXTERN_FUNCTION (int xv_destroy_safe, (Xv_object object));
EXTERN_FUNCTION (int xv_destroy_check, (Xv_object object));
EXTERN_FUNCTION (int xv_destroy, (Xv_object object));
EXTERN_FUNCTION (int xv_destroy_immediate,(Xv_object object));

/*
 * PRIVATE functions for pkg implementors only  
 */

EXTERN_FUNCTION (Xv_opaque xv_object_to_standard, (Xv_object object, const char *caller));

#if !(defined(__STDC__) || defined(__cplusplus) || defined(c_plusplus))
#define const 
#endif
extern const char *xv_notptr_str;
#define XV_OBJECT_SEAL          0xF0A58142
#define XV_OBJECT_TO_STANDARD(_passed_object, _caller, _object)\
{\
      if (!_passed_object) {\
        xv_error(NULL, ERROR_INVALID_OBJECT,xv_notptr_str,\
                 ERROR_STRING, _caller,\
                 0);\
        _object = ((Xv_opaque)0);\
    }\
    else\
       _object = (((Xv_base *)_passed_object)->seal == XV_OBJECT_SEAL) ? _passed_object : xv_object_to_standard(_passed_object, _caller);\
}

#endif /* xview_pkg_public_DEFINED */
