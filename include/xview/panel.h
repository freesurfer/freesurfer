/*	@(#)panel.h 20.89 91/09/14 SMI	*/

/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_panel_DEFINED
#define xview_panel_DEFINED

/*
 ***********************************************************************
 *		Include Files
 ***********************************************************************
 */

#include <xview/canvas.h>
#include <xview/frame.h>

/*
 ***********************************************************************
 *		Definitions and Macros
 ***********************************************************************
 */

/*
 * PUBLIC #defines 
 */
#define PANEL 		&xv_panel_pkg
#define PANEL_ABBREV_MENU_BUTTON &xv_panel_ambtn_pkg 
#define	PANEL_BUTTON	&xv_panel_button_pkg
#define	PANEL_CHOICE	&xv_panel_choice_pkg
#define PANEL_DROP_TARGET &xv_panel_drop_pkg
#define	PANEL_GAUGE	&xv_panel_gauge_pkg
#define	PANEL_ITEM	&xv_panel_item_pkg
#define	PANEL_LIST	&xv_panel_list_pkg
#define	PANEL_MESSAGE	&xv_panel_message_pkg
#define PANEL_MULTILINE_TEXT &xv_panel_multiline_text_pkg
#define PANEL_NUMERIC_TEXT  &xv_panel_num_text_pkg
#define	PANEL_SLIDER	&xv_panel_slider_pkg
#define	PANEL_TEXT	&xv_panel_text_pkg
#define PANEL_VIEW	&xv_panel_view_pkg

#define SCROLLABLE_PANEL &xv_scrollable_panel_pkg

#define PANEL_CHOICE_STACK	\
			PANEL_CHOICE, 	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT

#define PANEL_CYCLE	PANEL_CHOICE,	PANEL_DISPLAY_LEVEL,	PANEL_CURRENT

#define	PANEL_TOGGLE	PANEL_CHOICE, 	PANEL_CHOOSE_ONE, 	FALSE
#define PANEL_CHECK_BOX	PANEL_TOGGLE, 	PANEL_FEEDBACK, 	PANEL_MARKED

#define PANEL_CHECK_BOX_VALUE		PANEL_TOGGLE_VALUE


/*
 * Various utility macros 
 */
#define	panel_get_value(ip) 		xv_get((ip), PANEL_VALUE)
#define	panel_set_value(ip, val)	xv_set((ip), PANEL_VALUE, (val), 0)
#ifdef OW_I18N
#define panel_get_value_wcs(ip)         xv_get((ip), PANEL_VALUE_WCS)
#define panel_set_value_wcs(ip, val)    xv_set((ip), PANEL_VALUE_WCS, (val), 0)
#endif  /* OW_I18N */

/* Note: In PANEL_EACH_ITEM, we need "_next" since the current item
 * (ip) may be destroyed from within the for loop.
 */
#define	PANEL_EACH_ITEM(panel, ip)			\
   { Panel_item _next; 					\
   for ((ip) = xv_get((panel), PANEL_FIRST_ITEM); 	\
	(ip);				 		\
	(ip) = _next) { 				\
	_next = xv_get((ip), PANEL_NEXT_ITEM);		\
	if (xv_get(ip, PANEL_ITEM_OWNER))		\
	  continue;					\
	{

#define	PANEL_END_EACH	}}}

/*
 * Miscellaneous constants 
 */
#define PANEL_ITEM_X_START	4	/* offset from left edge */
#define PANEL_ITEM_Y_START	4	/* offset from top edge */

/* Panel defined events.
 * These are given to the Panel's or Panel_item's event proc
 */
#define	PANEL_EVENT_FIRST	vuid_first(PANEL_DEVID)		/* 32000 */
#define	PANEL_EVENT_CANCEL	(PANEL_EVENT_FIRST + 0)		/* 32000 */

/*
 * PRIVATE #defines 
 */

#define	PANEL_ATTR(type, ordinal)	ATTR(ATTR_PKG_PANEL, type, ordinal)

/*
 * panel specific attribute types 
 */
#define	PANEL_INDEX_STRING		ATTR_TYPE(ATTR_BASE_UNUSED_FIRST,     2)
#define	PANEL_INDEX_PIXRECT_PTR		ATTR_TYPE(ATTR_BASE_UNUSED_FIRST + 1, 2)
#define	PANEL_INDEX_BOOLEAN		ATTR_TYPE(ATTR_BASE_UNUSED_FIRST + 2, 2)
#define	PANEL_INDEX_FONT		ATTR_TYPE(ATTR_BASE_UNUSED_FIRST + 3, 2)
#define	PANEL_INDEX_CLIENT_DATA		ATTR_TYPE(ATTR_BASE_UNUSED_FIRST + 4, 2)


#define PANEL_FONT		WIN_FONT
#define PANEL_TYPE		ATTR_PKG_PANEL

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef struct {
    Xv_canvas		parent_data;
    Xv_opaque		private_data;
} Xv_panel_or_item;	/* scrollable_panel, panel or item */


/*
 * Typedefs 
 */
typedef	Xv_opaque 		Panel;
typedef	Xv_opaque 		Scrollable_panel;
typedef	Xv_opaque	 	Panel_item;
typedef	Xv_opaque 		Panel_attribute_value;
typedef Xv_panel_or_item	Xv_panel;
typedef Xv_panel_or_item	Xv_item;
typedef Xv_opaque		Panel_view;
typedef Xv_item			Xv_panel_message;

typedef Panel_item		Panel_abbrev_menu_button_item;
typedef Panel_item		Panel_button_item;
typedef Panel_item              Panel_choice_item;
typedef Panel_item		Panel_drop_target_item;
typedef Panel_item		Panel_gauge_item;
typedef Panel_item              Panel_list_item;
typedef Panel_item              Panel_message_item; 
typedef Panel_item              Panel_slider_item;
typedef Panel_item              Panel_text_item;
typedef Panel_item              Panel_numeric_text_item;
typedef Panel_item		Panel_multiline_text_item;

/*
 * Enumerations 
 */
typedef enum {
	/*********************
	 * Public Attributes 
	 *********************/

	/* Panel and Panel_item attributes */
	PANEL_ACCEPT_KEYSTROKE	= PANEL_ATTR(ATTR_BOOLEAN,		   1),
	PANEL_CLIENT_DATA	= PANEL_ATTR(ATTR_OPAQUE,		  36),
	PANEL_EVENT_PROC	= PANEL_ATTR(ATTR_FUNCTION_PTR,		  44),
	PANEL_GINFO		= PANEL_ATTR(ATTR_OPAQUE,		  45),
	    /* PANEL_ITEM_X and PANEL_ITEM_Y must still be supported in addition
	     * to XV_X and XV_Y because they can be used on the panel. When used
	     * on the panel, PANEL_ITEM_X != XV_X.
	     */
	PANEL_ITEM_X		= PANEL_ATTR(ATTR_X,			  63),
	PANEL_ITEM_Y		= PANEL_ATTR(ATTR_Y,			  64),
	PANEL_ITEM_X_GAP	= PANEL_ATTR(ATTR_X,			  65),
	PANEL_ITEM_Y_GAP	= PANEL_ATTR(ATTR_Y,			  66),
	PANEL_LAYOUT		= PANEL_ATTR(ATTR_ENUM,			  82),
	PANEL_PAINT		= PANEL_ATTR(ATTR_ENUM,			 160),

	/* Panel attributes */
	PANEL_BACKGROUND_PROC	= PANEL_ATTR(ATTR_FUNCTION_PTR,		   2),
	PANEL_BLINK_CARET	= PANEL_ATTR(ATTR_BOOLEAN,		   4),
	PANEL_BOLD_FONT		= PANEL_ATTR(ATTR_PIXFONT_PTR,		   6),
	PANEL_CARET_ITEM	= PANEL_ATTR(ATTR_OPAQUE,		   8),
	PANEL_CURRENT_ITEM	= PANEL_ATTR(ATTR_OPAQUE,		   9),
	PANEL_DEFAULT_ITEM	= PANEL_ATTR(ATTR_OPAQUE,		  38),
	PANEL_FOCUS_PW		= PANEL_ATTR(ATTR_OPAQUE,		  39),
	PANEL_EXTRA_PAINT_HEIGHT = PANEL_ATTR(ATTR_Y, 			  46),
	PANEL_EXTRA_PAINT_WIDTH = PANEL_ATTR(ATTR_X,			  48),
	PANEL_FIRST_ITEM	= PANEL_ATTR(ATTR_OPAQUE,		  52),
	PANEL_FIRST_PAINT_WINDOW = PANEL_ATTR(ATTR_OPAQUE,		  35),
	PANEL_ITEM_X_POSITION	= PANEL_ATTR(ATTR_INT,			  47),
	PANEL_ITEM_Y_POSITION	= PANEL_ATTR(ATTR_INT,			  49),
	PANEL_NO_REDISPLAY_ITEM	= PANEL_ATTR(ATTR_BOOLEAN,		  43),
	PANEL_PRIMARY_FOCUS_ITEM = PANEL_ATTR(ATTR_OPAQUE,		  51),
	PANEL_REPAINT_PROC	= PANEL_ATTR(ATTR_FUNCTION_PTR,		 164),
	PANEL_STATUS		= PANEL_ATTR(ATTR_OPAQUE,		 166),

	/* Panel_item attributes */
	PANEL_BUSY		= PANEL_ATTR(ATTR_BOOLEAN,		   7),
	PANEL_CHILD_CARET_ITEM	= PANEL_ATTR(ATTR_OPAQUE,		  81),
	PANEL_INACTIVE		= PANEL_ATTR(ATTR_BOOLEAN,		  54),
	PANEL_ITEM_CLASS	= PANEL_ATTR(ATTR_OPAQUE,		  56),
	PANEL_ITEM_COLOR	= PANEL_ATTR(ATTR_INT,			  58),
	PANEL_ITEM_CREATED	= PANEL_ATTR(ATTR_BOOLEAN,		  59),
	PANEL_ITEM_DEAF		= PANEL_ATTR(ATTR_BOOLEAN,		  57),
	PANEL_ITEM_LABEL_RECT	= PANEL_ATTR(ATTR_RECT_PTR,		  55),
	PANEL_ITEM_MENU		= PANEL_ATTR(ATTR_OPAQUE,		  60),
	PANEL_ITEM_NTH_WINDOW	= PANEL_ATTR(ATTR_OPAQUE,		  77),
	PANEL_ITEM_NWINDOWS	= PANEL_ATTR(ATTR_INT,			  79),
	PANEL_ITEM_RECT		= PANEL_ATTR(ATTR_RECT_PTR,		  62),
	PANEL_ITEM_VALUE_RECT	= PANEL_ATTR(ATTR_RECT_PTR,		  71),
	PANEL_ITEM_WANTS_ADJUST	= PANEL_ATTR(ATTR_BOOLEAN,		  61),
	PANEL_ITEM_WANTS_ISO	= PANEL_ATTR(ATTR_BOOLEAN,		  69),
	PANEL_LABEL_BOLD	= PANEL_ATTR(ATTR_BOOLEAN, 		  68),
	PANEL_LABEL_FONT	= PANEL_ATTR(ATTR_PIXFONT_PTR,		  70),
	PANEL_LABEL_IMAGE	= PANEL_ATTR(ATTR_PIXRECT_PTR,		  72),
	PANEL_LABEL_STRING	= PANEL_ATTR(ATTR_STRING,		  74),
	PANEL_LABEL_WIDTH	= PANEL_ATTR(ATTR_INT,			  76),
	PANEL_LABEL_X		= PANEL_ATTR(ATTR_X,			  78),
	PANEL_LABEL_Y		= PANEL_ATTR(ATTR_Y,			  80),
	PANEL_NEXT_COL		= PANEL_ATTR(ATTR_INT,			 145),
	PANEL_NEXT_ITEM		= PANEL_ATTR(ATTR_OPAQUE,		 146),
	PANEL_NEXT_ROW		= PANEL_ATTR(ATTR_INT,			 147),
	PANEL_NOTIFY_PROC	= PANEL_ATTR(ATTR_FUNCTION_PTR,		 154),
	PANEL_NOTIFY_STATUS	= PANEL_ATTR(ATTR_INT,			 156),
	PANEL_OPS_VECTOR	= PANEL_ATTR(ATTR_OPAQUE,		 157),
	PANEL_VALUE_X		= PANEL_ATTR(ATTR_X,			 190),
	PANEL_VALUE_Y		= PANEL_ATTR(ATTR_Y,			 192),

	/* Panel_choice_item attributes */
	PANEL_CHOICES_BOLD	= PANEL_ATTR(ATTR_BOOLEAN,		  10),
	PANEL_CHOICE_COLOR	= PANEL_ATTR(ATTR_INT_PAIR,		  11),
	PANEL_CHOICE_FONT	= PANEL_ATTR(PANEL_INDEX_FONT,		  12),
	PANEL_CHOICE_FONTS	=
        	PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXFONT_PTR), 14),
	PANEL_CHOICE_IMAGE	= PANEL_ATTR(PANEL_INDEX_PIXRECT_PTR, 	  16),
	PANEL_CHOICE_IMAGES	= 
        	PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXRECT_PTR), 18),
	PANEL_CHOICE_NCOLS      = PANEL_ATTR(ATTR_INT,			  19),
	PANEL_CHOICE_NROWS      = PANEL_ATTR(ATTR_INT,			  21),
	PANEL_CHOICE_RECT	= PANEL_ATTR(ATTR_INT,			  73),
	PANEL_CHOICE_STRING	= PANEL_ATTR(PANEL_INDEX_STRING,	  20), 
	PANEL_CHOICE_STRINGS	=
       		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_STRING),	  22),
	PANEL_CHOICE_X		= PANEL_ATTR(ATTR_INDEX_X,		  26),
	PANEL_CHOICE_XS	= PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_X), 28),
	PANEL_CHOICE_Y		= PANEL_ATTR(ATTR_INDEX_Y,		  30),
	PANEL_CHOICE_YS	= PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_Y), 32),
	PANEL_DEFAULT_VALUE	= PANEL_ATTR(ATTR_OPAQUE,		  40),
	PANEL_DISPLAY_LEVEL	= PANEL_ATTR(ATTR_ENUM,			  42),
	PANEL_FEEDBACK		= PANEL_ATTR(ATTR_ENUM,			  50),
	PANEL_MARK_IMAGE	= PANEL_ATTR(PANEL_INDEX_PIXRECT_PTR,	 124),
	PANEL_MARK_IMAGES		= 
        	PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXRECT_PTR),126),
	PANEL_MARK_X		= PANEL_ATTR(ATTR_INDEX_X,		 128),
	PANEL_MARK_XS	= PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_X),130),
	PANEL_MARK_Y		= PANEL_ATTR(ATTR_INDEX_Y,		 132),
	PANEL_MARK_YS	= PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_Y),134),
	PANEL_NCHOICES		= PANEL_ATTR(ATTR_INT,			 143),
	PANEL_NOMARK_IMAGE	= PANEL_ATTR(PANEL_INDEX_PIXRECT_PTR,	 148),
	PANEL_NOMARK_IMAGES		=
        	PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXRECT_PTR),150),
	PANEL_TOGGLE_VALUE	= PANEL_ATTR(PANEL_INDEX_BOOLEAN,	 178),

	/* Panel_choice_item, Panel_gauge_item, Panel_multiline_text_item,
	 * Panel_numeric_text_item, Panel_slider_item and Panel_text_item
	 * attributes.
	 */
	PANEL_VALUE		= PANEL_ATTR(ATTR_INT, 			 180),

	/* Panel_choice_item and Panel_list_item attributes */
	PANEL_CHOOSE_NONE	= PANEL_ATTR(ATTR_BOOLEAN,		  33),
	PANEL_CHOOSE_ONE	= PANEL_ATTR(ATTR_BOOLEAN,		  34),

	/* Panel_drop_target_item attributes */
	PANEL_DROP_BUSY_GLYPH	= PANEL_ATTR(ATTR_OPAQUE,		   3),
	PANEL_DROP_DND		= PANEL_ATTR(ATTR_OPAQUE,		   5),
	PANEL_DROP_FULL		= PANEL_ATTR(ATTR_BOOLEAN,		  13),
	PANEL_DROP_GLYPH	= PANEL_ATTR(ATTR_OPAQUE,		  15),
	PANEL_DROP_HEIGHT	= PANEL_ATTR(ATTR_INT,			  17),
	PANEL_DROP_SEL_REQ	= PANEL_ATTR(ATTR_OPAQUE,		  25),
	PANEL_DROP_SITE_DEFAULT	= PANEL_ATTR(ATTR_BOOLEAN,		  83),
	PANEL_DROP_WIDTH	= PANEL_ATTR(ATTR_INT,			  27),

	/* Panel_gauge_item attributes */
	PANEL_GAUGE_WIDTH	= PANEL_ATTR(ATTR_INT,			  53),

	/* Panel_gauge_item, Panel_numeric_text_item and
	 * Panel_slider_item attributes
	 */
	PANEL_MAX_VALUE		= PANEL_ATTR(ATTR_INT,			 138),
	PANEL_MIN_VALUE		= PANEL_ATTR(ATTR_INT,			 144),

	/* Panel_gauge_item and Panel_slider_item attributes */
	PANEL_DIRECTION		= PANEL_ATTR(ATTR_ENUM,			  41),
	PANEL_MAX_TICK_STRING	= PANEL_ATTR(ATTR_STRING,		 137),
	PANEL_MIN_TICK_STRING	= PANEL_ATTR(ATTR_STRING,		 140),
	PANEL_SHOW_RANGE	= PANEL_ATTR(ATTR_BOOLEAN,		 172),
	PANEL_TICKS		= PANEL_ATTR(ATTR_INT,			 177),

	/* Panel_list_item attributes */
	PANEL_LIST_CLIENT_DATA	= PANEL_ATTR(PANEL_INDEX_CLIENT_DATA,	  88),
	PANEL_LIST_CLIENT_DATAS	=
       		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_OPAQUE),	  90),
	PANEL_LIST_DELETE	= PANEL_ATTR(ATTR_INT,			  92),
	PANEL_LIST_DELETE_ROWS	= PANEL_ATTR(ATTR_INT_PAIR,		  93),
	PANEL_LIST_DELETE_SELECTED_ROWS = PANEL_ATTR(ATTR_NO_VALUE,	  97),
	PANEL_LIST_DISPLAY_ROWS	= PANEL_ATTR(ATTR_INT,			  94),
	PANEL_LIST_FIRST_SELECTED = PANEL_ATTR(ATTR_NO_VALUE,		  95),
	PANEL_LIST_FONT		= PANEL_ATTR(PANEL_INDEX_FONT,		  96),
	PANEL_LIST_FONTS		=
       		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXFONT_PTR), 98),
	PANEL_LIST_GLYPH	= PANEL_ATTR(PANEL_INDEX_PIXRECT_PTR,	 100),
	PANEL_LIST_GLYPHS	=
       		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_PIXRECT_PTR),102),
	PANEL_LIST_INSERT	= PANEL_ATTR(ATTR_INT,			 106),
	PANEL_LIST_INSERT_DUPLICATE = PANEL_ATTR(ATTR_BOOLEAN,		 105),
	PANEL_LIST_INSERT_GLYPHS = PANEL_ATTR(ATTR_OPAQUE_PAIR,		 107),
	PANEL_LIST_INSERT_STRINGS = PANEL_ATTR(ATTR_OPAQUE_PAIR,	 109),
	PANEL_LIST_MODE		= PANEL_ATTR(ATTR_ENUM,			 119),
	PANEL_LIST_NEXT_SELECTED = PANEL_ATTR(ATTR_INT,			  29),
	PANEL_LIST_NROWS	= PANEL_ATTR(ATTR_INT,			 108),
	PANEL_LIST_ROW_HEIGHT	= PANEL_ATTR(ATTR_INT,			 110),
	PANEL_LIST_SCROLLBAR	= PANEL_ATTR(ATTR_OPAQUE,		 111),
	PANEL_LIST_SELECT	= PANEL_ATTR(ATTR_INT_PAIR,	 	 112),
	PANEL_LIST_SELECTED	= PANEL_ATTR(ATTR_INT,			 113),
	PANEL_LIST_SORT		= PANEL_ATTR(ATTR_ENUM,			 117),
	PANEL_LIST_STRING	= PANEL_ATTR(PANEL_INDEX_STRING,	 114),
	PANEL_LIST_STRINGS	=
       		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_STRING),	 116),
	PANEL_LIST_TITLE	= PANEL_ATTR(ATTR_STRING,		  31),
	PANEL_LIST_WIDTH	= PANEL_ATTR(ATTR_INT,			 122),

	/* Panel_list_item, Panel_multiline_text_item,
	 * Panel_numeric_text_item, Panel_slider_item and
	 * Panel_text_item attributes
	 */
	PANEL_READ_ONLY		= PANEL_ATTR(ATTR_BOOLEAN,		 162),

	/* Panel_multiline_text_item attributes */
	PANEL_LINE_BREAK_ACTION	= PANEL_ATTR(ATTR_ENUM,			  84),

	/* Panel_multiline_text_item, Panel_numeric_text_item,
	 * Panel_slider_item and Panel_text_item attributes
	 */
	PANEL_NOTIFY_LEVEL	= PANEL_ATTR(ATTR_ENUM,			 152),
	PANEL_VALUE_DISPLAY_LENGTH	= PANEL_ATTR(ATTR_INT,		 182),

	/* Panel_multiline_text_item, Panel_numeric_text_item and
	 * Panel_text_item attributes
	 */
	PANEL_NOTIFY_STRING	= PANEL_ATTR(ATTR_STRING,		 158),
	PANEL_VALUE_DISPLAY_WIDTH = PANEL_ATTR(ATTR_INT,		 183),
	PANEL_VALUE_STORED_LENGTH	= PANEL_ATTR(ATTR_INT,		 186),

	/* Panel_numeric_text_item and Panel_slider_item attributes */
	PANEL_JUMP_DELTA	= PANEL_ATTR(ATTR_INT,			  67),
	PANEL_VALUE_FONT	= PANEL_ATTR(ATTR_PIXFONT_PTR,		 184),

	/* Panel_numeric_text_item and Panel_text_item attributes */
	PANEL_MASK_CHAR		= PANEL_ATTR(ATTR_CHAR,			 136),
	PANEL_VALUE_UNDERLINED	= PANEL_ATTR(ATTR_BOOLEAN,		 188),

	/* Panel_slider_item attributes */
	PANEL_MAX_VALUE_STRING	= PANEL_ATTR(ATTR_STRING,		 139),
	PANEL_MIN_VALUE_STRING	= PANEL_ATTR(ATTR_STRING,		 142),
	PANEL_SHOW_VALUE	= PANEL_ATTR(ATTR_BOOLEAN,		 174),
	PANEL_SLIDER_END_BOXES	= PANEL_ATTR(ATTR_BOOLEAN,		 175),
	PANEL_SLIDER_WIDTH	= PANEL_ATTR(ATTR_X,			 176),

	/* Panel_text_item attributes */
	PANEL_TEXT_SELECT_LINE	= PANEL_ATTR(ATTR_NO_VALUE,		 210),

	/*********************
	 * Private Attributes 
	 *********************/
	PANEL_CHOICE_OFFSET	= PANEL_ATTR(ATTR_INT,			 196),
	PANEL_ITEM_OWNER	= PANEL_ATTR(ATTR_OPAQUE,		 197),
	PANEL_LABEL_BOXED	= PANEL_ATTR(ATTR_BOOLEAN,		 198),
	PANEL_LABEL_INVERTED	= PANEL_ATTR(ATTR_BOOLEAN,		 200),
	PANEL_LIST_ROW_HIGHLIGHTED = PANEL_ATTR(ATTR_BOOLEAN,		 202),
	PANEL_MENU_ITEM		= PANEL_ATTR(ATTR_BOOLEAN,		 204),
	PANEL_NUM_TXT		= PANEL_ATTR(ATTR_OPAQUE,		 205),
	PANEL_TEXT_CURSOR	= PANEL_ATTR(ATTR_OPAQUE,		 206)

#ifdef OW_I18N
	/****************************************
	 * Internationalization Public Attributes
	 ****************************************/
                                                                             ,
	PANEL_LABEL_STRING_WCS	= PANEL_ATTR(ATTR_STRING,                 75),
	PANEL_NOTIFY_PROC_WCS	= PANEL_ATTR(ATTR_FUNCTION_PTR,          155),
	PANEL_CHOICE_STRING_WCS = PANEL_ATTR(PANEL_INDEX_STRING,          23),
	PANEL_CHOICE_STRINGS_WCS = 
		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_STRING),      24),
	PANEL_VALUE_WCS         = PANEL_ATTR(ATTR_INT,                   181),
	PANEL_LIST_INSERT_STRINGS_WCS = 
		PANEL_ATTR(ATTR_OPAQUE_PAIR,	 165),
	PANEL_LIST_STRING_WCS   = PANEL_ATTR(PANEL_INDEX_STRING,         115),
	PANEL_LIST_STRINGS_WCS  =
		PANEL_ATTR(ATTR_LIST_INLINE(ATTR_NULL, ATTR_STRING),     118),
	PANEL_LIST_TITLE_WCS	= PANEL_ATTR(ATTR_STRING,		 163),
	PANEL_MASK_CHAR_WC	= PANEL_ATTR(ATTR_CHAR,			 161),
	PANEL_NOTIFY_STRING_WCS = PANEL_ATTR(ATTR_STRING,                159),
	PANEL_MAX_TICK_STRING_WCS = PANEL_ATTR(ATTR_STRING,		 141),
	PANEL_MIN_TICK_STRING_WCS = PANEL_ATTR(ATTR_STRING,		 149),
	PANEL_MAX_VALUE_STRING_WCS = PANEL_ATTR(ATTR_STRING,		 151),
	PANEL_MIN_VALUE_STRING_WCS = PANEL_ATTR(ATTR_STRING,		 153),

	/******************************************
	 * Internationalization Private Attributes
	 ******************************************/

	PANEL_TEXT_IGNORE_IM	= PANEL_ATTR(ATTR_INT,			 167)
#endif /* OW_I18N */
} Panel_attr;

#define PANEL_PARENT_PANEL	XV_OWNER
#define PANEL_SHOW_ITEM		XV_SHOW
#define PANEL_DISPLAY_ROWS	PANEL_LIST_DISPLAY_ROWS

typedef enum {
    PANEL_LIST_OP_DESELECT,
    PANEL_LIST_OP_SELECT,
    PANEL_LIST_OP_VALIDATE,
    PANEL_LIST_OP_DELETE
} Panel_list_op;

typedef enum {
    PANEL_LIST_READ,
    PANEL_LIST_EDIT
} Panel_list_mode;

typedef struct panel_ops {
   void   (*panel_op_handle_event)();
   void   (*panel_op_begin_preview)();
   void   (*panel_op_update_preview)();
   void   (*panel_op_cancel_preview)();
   void   (*panel_op_accept_preview)();
   void   (*panel_op_accept_menu)();
   void   (*panel_op_accept_key)();
   void	  (*panel_op_clear)();
   void   (*panel_op_paint)();
   void	  (*panel_op_resize)();
   void   (*panel_op_remove)();
   void   (*panel_op_restore)();
   void   (*panel_op_layout)();
   void   (*panel_op_accept_kbd_focus)();
   void   (*panel_op_yield_kbd_focus)();
   void	   *panel_op_extension;		/* for future extensions */
} Panel_ops;

typedef struct panel_pw_struct {
	Xv_Window	pw;
	Xv_Window	view;
	struct panel_pw_struct	*next;
} Panel_paint_window;

/*
 * 		Panel status flags.  (Used in panel->status)
 *
 * Note: These values are readable but not settable at the Panel_item level.
 */
typedef struct panel_status {
    unsigned blinking : 1;
    unsigned current_item_active : 1;
    unsigned destroying : 1;
    unsigned focus_item_set : 1;
	/* The keyboard focus item has been set when the panel paint window
	 * did not have the keyboard focus.  When the panel paint window
	 * receives the KBD_USE event, don't change the kbd_focus to the
	 * primary focus item.  Instead, just clear the flag.
	 */
    unsigned has_input_focus : 1;
    unsigned mouseless : 1;
    unsigned nonstd_cursor : 1;	/* cursor is not basic pointer */
    unsigned painted : 1;
    unsigned pointer_grabbed : 1;
    unsigned quick_move : 1;	/* quick move pending */
    unsigned select_displays_menu : 1;
    unsigned three_d : 1;
} Panel_status;


/*
 * Values for LEVEL attributes 
 */
typedef enum {
    PANEL_CLEAR,        	/* painting */
    PANEL_NO_CLEAR,     	/* painting */
    PANEL_NONE,			/* text notify, menu, display, painting */
    PANEL_ALL,			/* text notify, slider notify, display */
    PANEL_NON_PRINTABLE,	/* text notify */
    PANEL_SPECIFIED,		/* text notify */
    PANEL_CURRENT,		/* display */
    PANEL_DONE,			/* slider notify */
    PANEL_MARKED,		/* feedback */
    PANEL_VERTICAL,		/* layout, slider direction */
    PANEL_HORIZONTAL,		/* layout, slider direction */
    PANEL_INVERTED,		/* feedback */
    /*
     * values returnable by notify routines 
     */
    PANEL_INSERT,
    PANEL_NEXT,
    PANEL_PREVIOUS,
    /*
     * mouse state 
     */
    PANEL_NONE_DOWN,		/* no buttons are down */
    PANEL_LEFT_DOWN,		/* left button is down */
    PANEL_MIDDLE_DOWN,		/* middle button is down */
    PANEL_RIGHT_DOWN,		/* right button is down */
    PANEL_CHORD_DOWN,		/* chord of buttons are down */
    /*
     * PANEL_LIST sort keys
     */
    PANEL_FORWARD,		/* forward alphabetic sort */
    PANEL_REVERSE,		/* reverse alphabetic sort */
    /*
     * PANEL_LINE_BREAK_ACTION values
     */
    PANEL_WRAP_AT_CHAR,
    PANEL_WRAP_AT_WORD
} Panel_setting;

/*
 * Types of items
 */
typedef enum { 
   PANEL_MESSAGE_ITEM,
   PANEL_BUTTON_ITEM, 
   PANEL_CHOICE_ITEM, 
   PANEL_TOGGLE_ITEM,
   PANEL_TEXT_ITEM,
   PANEL_NUMERIC_TEXT_ITEM,
   PANEL_SLIDER_ITEM,
   PANEL_LIST_ITEM,
   PANEL_GAUGE_ITEM,
   PANEL_ABBREV_MENU_BUTTON_ITEM,
   PANEL_EXTENSION_ITEM,
   PANEL_MULTILINE_TEXT_ITEM,
   PANEL_DROP_TARGET_ITEM
} Panel_item_type;
    
typedef struct {
    Xv_canvas_view	parent_data;
    Xv_opaque		private_data;
} Xv_panel_view;

/*
 * Structures 
 */

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_choice;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_gauge;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_slider;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_text;

typedef struct {
    Xv_item             parent_data;
    Xv_opaque           private_data;
} Xv_panel_num_text;

typedef struct {
    Xv_item             parent_data;
    Xv_opaque           private_data;
} Xv_panel_multiline_text;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_ambtn;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_button;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_drop;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_list;

typedef struct {
    Xv_item		parent_data;
    Xv_opaque		private_data;
} Xv_panel_extension_item;

/*
 * For SunView 1 compatibility 
 */
typedef Panel_attr	Panel_attribute;

/*
 ***********************************************************************
 *		Globals
 ***********************************************************************
 */

/*
 * 	Variables 
 */
extern Xv_pkg 		xv_panel_pkg;
extern Xv_pkg		xv_panel_ambtn_pkg;
extern Xv_pkg		xv_panel_button_pkg;
extern Xv_pkg		xv_panel_choice_pkg;
extern Xv_pkg		xv_panel_drop_pkg;
extern Xv_pkg		xv_panel_gauge_pkg;
extern Xv_pkg		xv_panel_item_pkg;
extern Xv_pkg		xv_panel_list_pkg;
extern Xv_pkg		xv_panel_message_pkg;
extern Xv_pkg           xv_panel_multiline_text_pkg;
extern Xv_pkg           xv_panel_num_text_pkg;
extern Xv_pkg		xv_panel_slider_pkg;
extern Xv_pkg		xv_panel_text_pkg;
extern Xv_pkg		xv_panel_view_pkg;
extern Xv_pkg		xv_scrollable_panel_pkg;

/*
 * 	Public Functions 
 */

/*
 * Panel routines 
 */
EXTERN_FUNCTION (Panel_item 	panel_advance_caret, (Panel panel));
EXTERN_FUNCTION (Panel_item 	panel_backup_caret, (Panel panel));
EXTERN_FUNCTION (void		panel_show_focus_win, (Panel_item item, Frame frame, int x, int y));

/*
 * event mapping routines 
 */
EXTERN_FUNCTION (int panel_handle_event, (Panel_item item, Event *event));
EXTERN_FUNCTION (void panel_default_handle_event, (Panel_item item, Event *event));
EXTERN_FUNCTION (int panel_cancel, (Panel_item item, Event *event));

/*
 * Panel_item action routines 
 */
EXTERN_FUNCTION (int panel_begin_preview, (Panel_item item, Event * event));
EXTERN_FUNCTION (int panel_update_preview, (Panel_item item, Event *event));
EXTERN_FUNCTION (int panel_accept_preview, (Panel_item item, Event *event));
EXTERN_FUNCTION (int panel_cancel_preview, (Panel_item item, Event *event));
EXTERN_FUNCTION (int panel_accept_menu, (Panel_item item, Event *event));
EXTERN_FUNCTION (int panel_accept_key, (Panel_item item, Event *event));

/*
 * utilities 
 */
EXTERN_FUNCTION (Panel_setting panel_text_notify, (Panel_item item, Event *event));
EXTERN_FUNCTION (void panel_paint_label, (Panel_item item));
EXTERN_FUNCTION (struct pixrect *panel_button_image, (Panel panel, char *string, int width, Xv_opaque font));
EXTERN_FUNCTION (void panel_default_clear_item, (Panel_item item));

/* routines to translate event coordinates
 * Note that struct inputevent * is the same as
 * Event *, this is used for compatibility with previous
 * releases.
 */
EXTERN_FUNCTION (Event * panel_window_event, (Panel panel, Event *event));
EXTERN_FUNCTION (Event * panel_event, (Panel panel, Event *event));

/*
 * For SunView 1 Compatibility Only
 */

/*
 * Panel & Panel_item routines in XView should use xv_ routines instead.
 */
EXTERN_FUNCTION (Panel_attribute_value panel_get, (Panel panel, Panel_attr attr, DOTDOTDOT));
EXTERN_FUNCTION (int panel_set, (Panel panel, DOTDOTDOT));
EXTERN_FUNCTION (int panel_paint, (Panel panel, Panel_setting flag));
EXTERN_FUNCTION (int panel_free, (Panel panel));
EXTERN_FUNCTION (int panel_destroy_item, (Panel_item item));
EXTERN_FUNCTION (Panel_item panel_create_item, (Panel panel, Xv_pkg *item_type, DOTDOTDOT));

#endif	/* ~xview_panel_DEFINED */
