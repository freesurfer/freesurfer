/*      @(#)window_hs.h 20.19 91/09/14	 */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*  WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 *  Added 6/91:
 *  This header file is left here for compatibility only, it will be removed
 *  in some future release of the toolkit.  
 *  WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
 */

#ifndef	xview_window_hs_DEFINED
#define	xview_window_hs_DEFINED	

/*
 * Include this header file to get all sunwindow related header files.
 */

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <sys/time.h>
#include <xview/notify.h>

#ifndef	pixrect_hs_DEFINED
#define	pixrect_hs_DEFINED
/* <pixrect/pixrect_hs.h> without frame buffer variable include files */
#include <sys/types.h>
#include <pixrect/pixrect.h>
#include <pixrect/pr_dblbuf.h>
#include <pixrect/pr_line.h>
#include <pixrect/pr_planegroups.h>
#include <pixrect/pr_util.h>
#include <pixrect/traprop.h>

#ifdef __STDC__ 
#ifndef CAT
#define CAT(a,b)        a ## b 
#endif 
#endif
#include <pixrect/memvar.h>

#include <pixrect/pixfont.h>

#ifndef XV_OS_SVR4
#include <rasterfile.h>
#endif

#include <pixrect/pr_io.h>
#endif	/* pixrect_hs_DEFINED */

#include <xview/rect.h>
#include <xview/rectlist.h>
#include <xview/pixwin.h>
#include <xview/win_struct.h>
#include <xview/win_screen.h>
#include <xview/win_input.h>
#include <xview/win_notify.h>

#endif	/* ~xview_window_hs_DEFINED */
