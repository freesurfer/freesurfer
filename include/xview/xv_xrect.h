/*      @(#)xv_xrect.h 1.13 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef _xv_xrect_h_already_included
#define _xv_xrect_h_already_included
#include <X11/Xlib.h>

/*
 * To enable applications to get the current clipping list
 * for direct X graphics.
 */
#define XV_MAX_XRECTS 32
typedef struct {
        XRectangle      rect_array[XV_MAX_XRECTS];
        int             count;
} Xv_xrectlist;

#endif /* _xv_xrect_h_already_included */
