/**
 * @file  hipspub.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/*
 * Copyright (c) 1992 Michael Landy
 *
 * Disclaimer:  No guarantees of performance accompany this software,
 * nor is any responsibility assumed on the part of the authors.  All the
 * software has been tested extensively and every effort has been made to
 * insure its reliability.
 */

/*
 * hipspub.h - definitions of raster HIPS formats
 *
 * Michael Landy - 9/30/92
 */

/*
 * Pixel Format Codes (raster formats only)
 */

#define PFBYTE  0 /* Bytes interpreted as unsigned integers */
#define PFSHORT  1 /* Short integers (2 bytes) */
#define PFINT  2 /* Integers (4 bytes) */
#define PFFLOAT  3 /* Float's (4 bytes)*/
#define PFCOMPLEX  4 /* 2 Float's interpreted as (real,imaginary) */
#define PFASCII  5 /* ASCII rep, with linefeeds after each row */
#define PFDOUBLE  6 /* Double's (8 byte floats) */
#define PFDBLCOM  7 /* Double complex's (2 Double's) */
#define PFMSBF  30 /* packed, most-significant-bit first */
#define PFLSBF  31 /* packed, least-significant-bit first */
#define PFSBYTE  32 /* signed bytes */
#define PFUSHORT 33 /* unsigned shorts */
#define PFUINT  34 /* unsigned ints */
#define PFRGB  35 /* RGBRGBRGB bytes */
#define PFRGBZ  36 /* RGB0RGB0RGB0 bytes */

#define LINELENGTH 200  /* max characters per line in header vars */

#define MSBFIRST 1 /* HIPS-1 bit_packing value for PFMSBF */
#define LSBFIRST 2 /* HIPS-1 bit_packing value for PFLSBF */
