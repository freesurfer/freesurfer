/*      @(#)xv_i18n.h 1.5 91/09/14; SMI */
/*
 *      (c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *      pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *      file for terms of the license.
 */

#ifndef xv_i18n_h_DEFINED
#define xv_i18n_h_DEFINED

#ifdef OW_I18N

/*
 * This is public header file for the all XView I18N programs as well
 * as XView itself.  This file should provide the following
 * information in the all platforms.
 */
/*
 *	1. definition for the "wchar_t".
 *	2. definition for the MB_CUR_MAX, MB_LEN_MAX.
 *	3. definition for the setlocale() function and related
 *	    #define.
 *	4. definition for the wide character and multibyte
 *	   functions (ie, wscpy, mblen....).
 *	5. definieion for the wide character classification
 *	   functions.
 *	6. include the all i18n specific Xlib include files.
 */
/*
 * Also, because of different platform may have a different naming
 * schema for the wide character functions, you may want to provide
 * the macros to adpat to the specific platform.  The current code is
 * using AT&T MNLS and/or J/ALE naming schema (wsXXX, ie. wscpy).
 */

#ifdef SVR4


#	include <stdlib.h>		/* #2 (MB_CUR_MAX) */
#	include <limits.h>		/* #2 (MB_LEN_MAX) */
#	include <widec.h>		/* #1, #4 */
#	include <locale.h>		/* #3 */
#	include <wctype.h>		/* #5 */

#	include <X11/Xlib.h>		/* We should remove this - Ako */
#	include <xim/immgr.h>		/* #6 */
#	include <xim/immgr_ic.h>	/* #6 */
#	include <xim/immgr_cb.h>	/* #6 */
#	include <mltext/XFontSet.h>	/* #6 */

#else /* SVR4 */

	/*
	 * In case of the SunOS 4.X/JLE 1.X.
	 */
#	include <stdlib.h>		/* #2 (MB_CUR_MAX) */
#	include <limits.h>		/* #2 (MB_LEN_MAX) */
#	include <widec.h>		/* #1, #4 */
#	include <locale.h>		/* #3 */
#	include <wctype.h>		/* #5 */

	/*
	 * The following will only valid for the Sun's implementation
	 * of the i18n Xlib level functions.
	 */
#	include <X11/Xlib.h>		/* We should remove this - Ako */
#	include <xim/immgr.h>		/* #6 */
#	include <xim/immgr_ic.h>	/* #6 */
#	include <xim/immgr_cb.h>	/* #6 */
#	include <mltext/XFontSet.h>	/* #6 */

#endif /* SVR4 */

#endif /* OW_I18N */

#endif /* xv_i18n_h_DEFINED */
