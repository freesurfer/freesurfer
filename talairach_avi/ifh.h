/************************************************************************************************/
/* Copyright 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007	*/
/* Washington University Mallinckrodt Institute of Radiology. All Rights Reserved.		*/
/* This software may not be reproduced, copied, or distributed without written			*/
/* permission of Washington University. For further information contact A. Z. Snyder.		*/
/************************************************************************************************/
/*$Header: /space/repo/1/dev/dev/talairach_avi/ifh.h,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: ifh.h,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.3  2007/05/03  22:17:09  avi
 * linux gcc v3 compliant
 *
 * Revision 1.2  2006/03/16  06:33:26  avi
 * radical pruning of unneccesary fields
 * add endian field
 *
 * Revision 1.1  2004/03/03  02:26:41  avi
 * Initial revision
 **/

#ifndef IFH_INCLUDED
#define IFH_INCLUDED
/*________________________________________________________________________________ 
File:		ifh.h

Description:	General definition and includes for the image list checking software

Author:		Tom Yang

Date:		05/03/95

History:	Created by Tom Yang on May 3, 1995.
		Modified (center[3], mmppix[3]) Dec 10, 1998, AZS
________________________________________________________________________________*/ 
typedef struct {
	char	interfile[32];			/* <ASCII> interfile label */
	char	version_of_keys[32];		/* <ASCII> "3.3" */
	char	conversion_program[256];	/* <ASCII> program that generated the data */
	char	name_of_data_file[256];		/* <ASCII> image data filename */
	char	number_format[32];		/* <ASCII> "float" */
	char    imagedata_byte_order[32];       /* <ASCII> byte order ["bigendian" | "littleendian"] of the image data */
	int	number_of_bytes_per_pixel;	/* <Numeric> [1 | 2 | 4] */
	int	number_of_dimensions;		/* <Numeric> 4 */
	int	matrix_size[4];			/* <Numeric> image array dimensions */
	int	orientation;			/* <Numeric> 2 = Transverse; 3 = Coronal; 4 = sagittal */
	float	scaling_factor[4];		/* <Numeric> (x, y, z) (mm/pixel) */
	float	mmppix[3];			/* <Numeric> (x, y, z) */
	float	center[3];			/* <Numeric> (x, y, z) */
} IFH;

#endif /* IFH_INCLUDED */
