/*
 * @(#)win_notify.h 20.15 91/09/14 SMI
 *
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/*
 * SunView related notification definitions (also see notify.h).
 */
#ifndef win_notify_DEFINED
#define win_notify_DEFINED

#include <xview/xv_c_types.h>

/*
 ***********************************************************************
 *				Globals
 ***********************************************************************
 */

/*
 * PUBLIC Functions 
 */

/*
 * Posting of client events to window notifier clients 
 */

EXTERN_FUNCTION (Notify_error 	win_post_id, (Notify_client client, int id, Notify_event_type when));
		
EXTERN_FUNCTION (Notify_error 	win_post_id_and_arg, (Notify_client client, int id, Notify_event_type when, Notify_arg arg, Notify_copy copy_func, Notify_release release_func));

EXTERN_FUNCTION (Notify_error 	win_post_event, ( Notify_client client, Event *event, Notify_event_type when));

EXTERN_FUNCTION (Notify_error 	win_post_event_arg, ( Notify_client client, Event * event, Notify_event_type when, Notify_arg arg, Notify_copy copy_func, Notify_release release_func));

/*
 * Utilities to call if posting with win_post_id_and_arg or win_post_event_arg
 */

EXTERN_FUNCTION (Notify_arg win_copy_event, ( Notify_client client, Notify_arg arg, Event **event_ptr));

EXTERN_FUNCTION (void win_free_event, ( Notify_client client, Notify_arg arg, Event *event));
					
#endif /* win_notify_DEFINED */
