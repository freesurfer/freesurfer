/*      @(#)xview_xvin.h 1.13 91/09/14 SMI      */

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_xview_xvin_DEFINED
#define xview_xview_xvin_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <signal.h>

#include <sys/types.h>
#include <pixrect/pixrect.h>
#include <pixrect/pr_planegroups.h>
#include <pixrect/pr_util.h>

#ifdef __STDC__
#ifndef CAT
#define CAT(a,b)        a ## b
#endif
#endif
#include <pixrect/memvar.h>

#include <pixrect/pixfont.h>
#include <pixrect/traprop.h>
#include <pixrect/pr_line.h>

#if defined(__cplusplus) || defined(__STDC__)
#include <stdlib.h>
#endif /* __cplusplus || __STDC__ */

#include <xview/xv_c_types.h>   
#include <xview/generic.h>
#include <xview/server.h>
#include <xview/screen.h>

#include <xview/notify.h>
#include <xview/pixwin.h>
#include <xview/win_input.h>

#endif /* ~xview_xview_xvin_DEFINED */
