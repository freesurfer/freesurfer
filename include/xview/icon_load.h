/*	@(#)icon_load.h 20.16 91/09/14 SMI	*/
 
/*
 *	(c) Copyright 1989 Sun Microsystems, Inc. Sun design patents 
 *	pending in the U.S. and foreign countries. See LEGAL NOTICE 
 *	file for terms of the license.
 */

#ifndef xview_icon_load_DEFINED
#define xview_icon_load_DEFINED

#include <xview/svrimage.h>

#ifdef notdef
/*
 * The format of an ASCII representation of an icon file is defined by:
 *	<icon_file>	::= <header> <icon_image_bits>
 *	<header>	::= <header_start> <header_body> <header_end>
 *	<header_start>	::= /* Format_version=1
 *	<header_end>	::= */
 /* 	<header_body>	::= <header_options> <header_misc> 
 *	<header_options>::= | <header_option> <header_options> 
 * 	<header_option> ::= Depth=1
 *			  | Height=<Number>
 *			  | Width=<Multiple_of_16>
 *			  | Valid_bits_per_item=<Sixteen_or_ThirtyTwo>
 *	<header_misc>	::= <Anything_except_header_option_or_header_end>
 *	<icon_image_bits>
 *			::= <icon_image_value>
 * 			  | <icon_image_bits> , <icon_image_value>
 *
 * Default values for header parameters are:
 *	Depth			 1
 *	Height			64
 * 	Width			64
 *	Valid_bits_per_item	16
 * 
 * A sample icon file follows:
 *
 * Format_version=1, Width=16, Height=16, Depth=1, Valid_bits_per_item=16
 * This file is the template for all images in the cursor/icon library.
 * The first line contains the information needed to properly interpret the
 *   actual bits, which are expected to be used directly by software that
 *   wants to do compile-time binding to an image via a #include.
 * The actual bits must be specified in hex.
 * The default interpretation of the bits below is specified by the
 *   behavior of mpr_static.
 * Note that Valid_bits_per_item uses the least-significant bits.
 * Description: An empty (clear) cursor
 * Background: White
 *
 *	0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000,
 *	0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000
 *
 *
 * The utilities available to read icon files are:
 *
 *	int
 *	icon_load(icon, from_file, error_msg);
 *		register Icon	icon;
 *		char			*from_file, *error_msg;
 * Loads the pixrect for the icon, and then copies information from the pixrect
 * to initialize the icon.  Information which cannot be specified in the
 * current format is set to default values.
 * See the following description of icon_load_mpr for details about the
 * arguments other than icon.
 * Returns 0 iff it managed to successfully initialize the icon.
 *
 * 	int
 *	  icon_init_from_pr(icon, pr)
 *		register Icon	icon;
 *		register struct pixrect	*pr;
 * Initializes the icon struct from the specified pr.
 * Information, such as the font for any icon text, which cannot be ascertained
 * from a pixrect, is set to default values.
 *
 * The return value is meaningless.
 *
 *	struct pixrect *
 *	icon_load_mpr(from_file, error_msg)
 *		char		*from_file, *error_msg;
 * Loads the icon named by from_file into a dynamically allocated memory
 *   pixrect.
 * If NULL is returned, which happens if the file cannot be accessed or if the
 *   file is not in a valid format, an appropriate error message is placed
 *   in error_msg, which should be at least IL_ERRORMSG_SIZE long.
 *
 *	FILE *
 *	icon_open_header(from_file, error_msg, info)
 *		char				*from_file, *error_msg;
 *		register icon_header_handle	 info;
 * icon_open_header() allows a client to preserve extra descriptive material
 *   when rewriting an icon file.  It is also the front-end routine used by
 *   icon_load_mpr (and thus icon_load).
 * If NULL is returned, which happens if the file cannot be accessed or if the
 *   file is not in a valid format, an appropriate error message is placed
 *   in error_msg, which should be at least IL_ERRORMSG_SIZE long.
 * info is mostly filled in with the information from the header options.  The
 * exception is that info->last_param_pos is filled in with the position
 * immediately after the last header option that was read.
 * The returned FILE * is left positioned at the end of the header, thus
 * ftell(icon_open_header()) indicates where the actual bits of the image
 * are expected to begin.  Hence, the characters in the range
 * [info->last_param_pos..ftell(icon_open_header())) should contain all extra
 * descriptive material contained in the icon file's header.
 *
 * icon_read_pr is not a public routine
 *	int
 *	icon_read_pr(fd, header, pr)
 *		register FILE			*fd;
 * 		register icon_header_handle	 header;
 * 		register struct pixrect		*pr;
 * Reads the actual bit pattern for the pixrect's image data from the specified
 *  fd using the parameters in the specified header.
 * The fd is expected to be positioned to the beginning of the image data
 *  (usually by a prior call to icon_open_header).
 * The fd is left positioned at the end of the image data, so that by doing an
 *  ftell(fd) before and after calling icon_read_pr, a caller can write out a
 *  modified version of the image data while keeping intact the surrounding
 *  commentary, declarations, etc.
 * The return value is meaningless.
 * CAVEAT: currently pr must be a memory pixrect.
 *
 */
#endif /* notdef */

/*
 ***********************************************************************
 *			Definitions and Macros
 ***********************************************************************
 */

#define IL_ERRORMSG_SIZE   256

/*
 ***********************************************************************
 *		Typedefs, enumerations, and structs
 ***********************************************************************
 */

typedef struct icon_header_object {
	int	depth,
		height,
		format_version,
		valid_bits_per_item,
		width;
	long	last_param_pos;
} Xv_icon_header_info;

typedef Xv_icon_header_info	*icon_header_handle;
typedef Xv_icon_header_info	icon_header_object;

/*
 ***********************************************************************
 *			Globals
 ***********************************************************************
 */

/*
 * Public Functions
 */
EXTERN_FUNCTION (int icon_load, (Icon icon, char *from_file, char *error_msg));
EXTERN_FUNCTION (int icon_init_from_pr, (Icon icon, struct pixrect *pr));
EXTERN_FUNCTION (struct pixrect *icon_load_mpr, (char *from_file, char *error_msg));
EXTERN_FUNCTION (Server_image	icon_load_svrim, (char *from_file, char *error_msg));
EXTERN_FUNCTION (FILE *icon_open_header, (char *from_file, char *error_msg, Xv_icon_header_info *info));

#endif /* xview_icon_load_DEFINED */
