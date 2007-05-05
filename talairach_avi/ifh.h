/**
 * @file  ifh.h
 * @brief general definition and includes for the image list checking software
 */
/*
 * Original Author: Tom Yang on May 3, 1995
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/05/05 00:00:06 $
 *    $Revision: 1.2 $
 *
 * Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
 * Washington University, Mallinckrodt Institute of Radiology.
 * All Rights Reserved.
 *
 * This software may not be reproduced, copied, or distributed without 
 * written permission of Washington University. For further information 
 * contact A. Z. Snyder.
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#ifndef IFH_INCLUDED
#define IFH_INCLUDED

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
