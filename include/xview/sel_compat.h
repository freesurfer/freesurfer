/*	@(#)sel_compat.h 1.16 91/09/14		*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef	xview_selection_compat_DEFINED
#define	xview_selection_compat_DEFINED

#include <xview/server.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines
 * For Compatibility with SunView 1 
 */

#define seln_acquire(seln_client, asked)  	\
		selection_acquire(xv_default_server, seln_client, asked)
#define seln_done(seln_client, rank)  		\
		selection_done(xv_default_server, seln_client, rank)
#define seln_inform(seln_client, which, down)  	\
		selection_inform(xv_default_server, seln_client, which, down)
#define seln_get_function_state(func)  		\
		selection_get_function_state(server, func)
#define seln_inquire(which)			\
		selection_inquire(xv_default_server, which)
#define seln_clear_functions()   		\
		selection_clear_functions(xv_default_server)
#define	seln_request(holder, buffer) 		\
		selection_request(xv_default_server, holder, buffer)
#define seln_send_yield(rank, holder) 		\
		selection_send_yield(xv_default_server, rank, holder)
#define seln_function_pending()			\
		selection_function_pending(xv_default_sever)
#define seln_report_event(seln_client, event)	\
		selection_report_event(xv_default_server, seln_client, event)

#endif ~xview_selection_compat_DEFINED
