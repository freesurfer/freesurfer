/*      @(#)db.h 50.12 91/09/14 SMI */
/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_db_h_DEFINED
#define	xview_db_h_DEFINED

#include <xview/generic.h>
#include <xview/attr.h>
#include <X11/Xlib.h>
#include <X11/Xresource.h>

/* type converters */
#define db_cvt_int_to_string(from_value, to_value) \
            (void) sprintf(to_value, "%d", (int)from_value)

#define db_cvt_bool_to_string(from_value, to_value) \
            (void) strcpy(to_value, ((int)from_value) ? "True" : "False")

#define db_cvt_char_to_string(from_value, to_value) \
            to_value[0] = (char)from_value, to_value[1] = '\0'

#endif	/* ~xview_db_h_DEFINED */
