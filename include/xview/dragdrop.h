#ifndef lint
#ifdef sccs
static char     sccsid[] = "@(#)dragdrop.h 1.18 91/09/14";
#endif
#endif

/*
 *      (c) Copyright 1990 Sun Microsystems, Inc. Sun design patents 
 *      pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *      file for terms of the license.
 */

#ifndef xview_dnd_DEFINED
#define xview_dnd_DEFINED

#include <xview/attr.h>
#include <xview/sel_pkg.h>
#include <xview/xv_c_types.h>

/*
 * Private 
 */

#define DND_ATTR(type, ordinal)      ATTR(ATTR_PKG_DND, type, ordinal)
#define DND_ATTR_LIST(ltype, type, ordinal) \
        DND_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))

/*
 * Public 
 */

#define DRAGDROP		&xv_dnd_pkg
#define DROP_SITE_ITEM		&xv_drop_site_item
typedef Xv_opaque		Xv_drag_drop, Drag_drop, Dnd;
typedef Xv_opaque		Xv_drop_site, Drop_site_item;
	
				/* Error codes */
#define	DND_SUCCEEDED		XV_OK
#define	DND_ERROR		XV_ERROR
#define	DND_ILLEGAL_TARGET	2
#define	DND_TIMEOUT		3
#define	DND_SELECTION		4
#define	DND_ROOT		5
#define	DND_ABORTED		6 /* STOP key pressed */

		 	/* Drag type */
typedef enum {
	DND_MOVE=0,
	DND_COPY
} DndDragType;

#define DND_ENTERLEAVE		(1<<0)
#define DND_MOTION		(1<<1)
#define DND_DEFAULT_SITE        (1<<2)

#define DND_RECT_SITE		0
	/* XXX: Currently, only rect sites are supported. */
#define DND_WINDOW_SITE		1

#define DND_VERSION		0

#define dnd_is_local(event)	(event_flags(event) & DND_LOCAL)
#define dnd_is_forwarded(event)	(event_flags(event) & DND_FORWARDED)

/* Due to the unclear future for a standard drag and drop protocol, do not
 * rely on being able to access the drop site id in the X event.   It may
 * not always be there.  You've been warned! 
 */
#define dnd_site_id(event)      (event->ie_xevent->xclient.data.l[3])

			/* Private defines */
#define DND_MOVE_FLAG           (1<<0)
#define DND_ACK_FLAG            (1<<1)
#define DND_TRANSIENT_FLAG      (1<<2)
#define DND_FORWARDED_FLAG      (1<<3)
                        /* ie_flags */
#define DND_LOCAL               (1<<0)
#define DND_FORWARDED           (1<<1)


/* 
 * Public attributes
 */
typedef enum {
	/* Public Drag & Drop pkg attrs */
	DND_TYPE			= DND_ATTR(ATTR_SHORT,		1),
	DND_CURSOR			= DND_ATTR(ATTR_OPAQUE,		5),
	DND_X_CURSOR			= DND_ATTR(ATTR_LONG,		10),
	DND_ACCEPT_CURSOR		= DND_ATTR(ATTR_OPAQUE,		15),
	DND_ACCEPT_X_CURSOR		= DND_ATTR(ATTR_LONG,		20),
	DND_TIMEOUT_VALUE		= DND_ATTR(ATTR_LONG,		25),
	/* Private Drop Site Item attrs */
	DROP_SITE_SIZE			= DND_ATTR(ATTR_INT,		95),
	/* Public Drop Site Item attrs */
#ifdef WINDOW_SITES
	DROP_SITE_TYPE			= DND_ATTR(ATTR_INT,		100),
#endif /* WINDOW_SITES */
	DROP_SITE_ID			= DND_ATTR(ATTR_LONG,		105),
	DROP_SITE_EVENT_MASK		= DND_ATTR(ATTR_INT,		110),
	DROP_SITE_REGION		= DND_ATTR(ATTR_OPAQUE,		115),
	DROP_SITE_REGION_PTR		= DND_ATTR(ATTR_OPAQUE, 	120),
	DROP_SITE_DELETE_REGION		= DND_ATTR(ATTR_OPAQUE,		125),
	DROP_SITE_DELETE_REGION_PTR	= DND_ATTR(ATTR_OPAQUE, 	130),
	DROP_SITE_DEFAULT		= DND_ATTR(ATTR_BOOLEAN,	135)
} Dnd_attribute;

typedef struct {
	Xv_sel_owner		parent_data;
	Xv_opaque		private_data;
} Xv_dnd_struct, Xv_drop_site_struct;

/*
 * Public Functions 
 */
extern Xv_pkg	xv_dnd_pkg;
extern Xv_pkg	xv_drop_site_item;
EXTERN_FUNCTION (int            dnd_send_drop, (Xv_opaque object));
EXTERN_FUNCTION (Xv_opaque      dnd_decode_drop, (Xv_opaque object, Event *event));
EXTERN_FUNCTION (void      	dnd_done, (Xv_opaque object));

/* The function xv_decode_drop is obsolete, please use the new V3 
   drap and drop interface. The function is still provided for backwards
   compatibility for V2 programs.
*/
EXTERN_FUNCTION( int xv_decode_drop, ( Event *ev, char *buffer, unsigned int bsize ));

#endif /* ~xview_dnd_DEFINED */
