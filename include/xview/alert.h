/*	@(#)alert.h 20.15 91/09/14	*/

/* 
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

/* SunView 1.x "alerts" are now called "notices" in SunXView.x.  Below 
 * is for backwards compatibility only.  All usages of alert.h should 
 * be discontinued as further releases may not support this interface.
 * Not all "alert" attributes apply in SunXView.x, so full compatiblity
 * does not exist.
 * Include notice.h definitions and attributes instead.
 */

#ifndef xview_alert_DEFINED
#define xview_alert_DEFINED

/*
 ***********************************************************************
 *			Include Files
 ***********************************************************************
 */

#include <xview/notice.h>

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines
 */ 
#define Alert_attribute			Notice_attribute

#define ALERT_NO_BEEPING		NOTICE_NO_BEEPING
#define ALERT_MESSAGE_STRINGS		NOTICE_MESSAGE_STRINGS
#define ALERT_MESSAGE_STRINGS_ARRAY_PTR NOTICE_MESSAGE_STRINGS_ARRAY_PTR
#define ALERT_MESSAGE_FONT		NOTICE_FONT
#define ALERT_BUTTON_YES		NOTICE_BUTTON_YES
#define ALERT_BUTTON_NO			NOTICE_BUTTON_NO
#define ALERT_BUTTON			NOTICE_BUTTON
#define ALERT_TRIGGER			NOTICE_TRIGGER

/*
 * Useful constants 
 */
#define ALERT_YES			 1
#define ALERT_NO			 0
#define ALERT_FAILED			-1
#define ALERT_TRIGGERED			-2

#define alert_prompt			notice_prompt

/*
 * PRIVATE #defines 
 */
#define alert_attr_next(attr) (Notice_attribute *)attr_next((caddr_t *)attr)

#endif /* ~xview_alert_DEFINED */
