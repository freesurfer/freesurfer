/*	@(#)cms.h 20.12 89/08/17 SMI	*/
/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_cms_h_DEFINED
#define	xview_cms_h_DEFINED

#include <sys/types.h>

/*
 ***********************************************************************
 *                      Include files
 ***********************************************************************
 */

#include <xview/generic.h>


/*
 ***********************************************************************
 *                      Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines
 */
#define CMS			&xv_cms_pkg
#define Xv_Singlecolor		Xv_singlecolor
#define CMS_NAME		XV_NAME
#define XV_DEFAULT_CMS_SIZE	2

#define CMS_CONTROL_BG1		0
#define CMS_CONTROL_BG2		1
#define CMS_CONTROL_BG3		2
#define CMS_CONTROL_HIGHLIGHT	3
#define CMS_CONTROL_COLORS	4

/*
 * Private #defines
 */
#define CMS_ATTR(type, ordinal)	    ATTR(ATTR_PKG_CMS, type, ordinal)
#define CMS_ATTR_LIST(ltype, type, ordinal) \
	CMS_ATTR(ATTR_LIST_INLINE((ltype), (type)), (ordinal))
#define CMS_CONTROL_CMS_NAME(name, visual, color) \
	    sprintf(name, "xv_ctl_%x_%d%d%d", visual, \
		    color.red, color.green, color.blue);
#define CMS_STATUS(bits, status)    ((bits) & (1 << (status)))

#define CMS_STATUS_DEFAULT	0
#define CMS_STATUS_CONTROL	1
#define CMS_STATUS_FRAME	2

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef Xv_opaque	Cms;

typedef enum Xv_cmstype	{
    XV_STATIC_CMS    = 1, 
    XV_DYNAMIC_CMS   = 2
} Cms_type;

typedef enum {
    /*
     * Public attributes.
     */
    CMS_TYPE		= CMS_ATTR(ATTR_ENUM,	    	    	5),
    CMS_SIZE    	= CMS_ATTR(ATTR_INT, 	    	    	10),
    CMS_COLORS  	= CMS_ATTR(ATTR_SINGLE_COLOR_PTR,   	15),
    CMS_NAMED_COLORS  	= CMS_ATTR_LIST(ATTR_NULL, ATTR_STRING,	16),
    CMS_X_COLORS  	= CMS_ATTR(ATTR_OPAQUE, 		17),
    CMS_INDEX   	= CMS_ATTR(ATTR_INT,        	    	20),
    CMS_COLOR_COUNT   	= CMS_ATTR(ATTR_INT,        	    	25),
    CMS_SCREEN		= CMS_ATTR(ATTR_OPAQUE,		    	30),
    CMS_INDEX_TABLE	= CMS_ATTR(ATTR_OPAQUE,             	35),
    CMS_FOREGROUND_PIXEL= CMS_ATTR(ATTR_INT,                	40),
    CMS_BACKGROUND_PIXEL= CMS_ATTR(ATTR_INT,                	45),
    CMS_PIXEL		= CMS_ATTR(ATTR_INT,                	50),
    CMS_CMS_DATA	= CMS_ATTR(ATTR_OPAQUE,		    	55),
    CMS_CONTROL_CMS	= CMS_ATTR(ATTR_BOOLEAN,            	60),
    /*
     * Private attributes.
     */
    CMS_DEFAULT_CMS	= CMS_ATTR(ATTR_BOOLEAN,	    	70),
    CMS_FRAME_CMS	= CMS_ATTR(ATTR_BOOLEAN,	    	75),
    CMS_CMAP_ID		= CMS_ATTR(ATTR_OPAQUE,		    	80),
    CMS_CMAP_TYPE	= CMS_ATTR(ATTR_ENUM,		    	85),
    CMS_STATUS_BITS	= CMS_ATTR(ATTR_INT,		    	90)
} Cms_attribute;

typedef struct {
    Xv_generic_struct   parent_data;
    Xv_opaque           private_data;
} Xv_cms_struct;


typedef struct xv_cmsdata {
    Cms_type		 type;
    short                size;
    short                index;
    short                rgb_count;
    unsigned char       *red;
    unsigned char       *green;
    unsigned char       *blue;
} Xv_cmsdata;

typedef struct xv_singlecolor {
    u_char 	red, green, blue;
} Xv_singlecolor;


/* 
 * do we really need this ? check and remove.
 */
typedef struct xv_color_info {
    Cms     		cms;
    int			cms_fg;	/* CMS index of the foreground color */
    int			cms_bg; /* CMS index of the background color */
} Xv_Color_info;

/*
 *    Private typedefs
 */
typedef enum {
    XV_STATIC_CMAP    = 1, 
    XV_DYNAMIC_CMAP   = 2
} Cmap_type;

/*
 ***********************************************************************
 *                              Globals
 ***********************************************************************
 */
extern Xv_pkg	xv_cms_pkg;

#endif	/* ~xview_cms_h_DEFINED */
