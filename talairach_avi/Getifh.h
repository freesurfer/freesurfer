/**************************************************************************************/
/* Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006                           */
/* Washington University, Mallinckrodt Institute of Radiology.                        */
/* All Rights Reserved.                                                               */
/* This software may not be reproduced, copied, or distributed without written        */
/* permission of Washington University. For further information contact A. Z. Snyder. */
/**************************************************************************************/
/*$Header: /space/repo/1/dev/dev/talairach_avi/Getifh.h,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: Getifh.h,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.1  2006/09/23  05:31:16  avi
 * Initial revision
 **/

#include <ifh.h>

/************************************************************************/
/* read IFH file associated with imgfile and populate the IFH structure */
/************************************************************************/
extern int Getifh 	(char *imgfile, IFH *ifhdr);

/*****************************************************************************************************/
/* write IFH file from IFH strucrure. argument control sets the value of the image_byte_order field. */
/*****************************************************************************************************/
extern int Writeifh 	(char *program, char *outfile, IFH *ifhdr, char control);

/*****************************************************************************/
/* write IFH file from the arguments passed. omits mmppix and center values. */
/*****************************************************************************/
extern int writeifhe 	(char *program, char *outfile, int *imgdim, float *voxdim, int orient, char control); 

/****************************************************************/
/* writeifhe() functionality including mmppix and center values */
/****************************************************************/
extern int writeifhmce 	(char *program, char *outfile, int *imgdim, float *voxdim, int orient,float *mmppix, float *center, char control);

/************************************************************************************************/
/* writeifhmce() functionality excluding control. Output ifh imagedata byte order is CPU-endian */
/************************************************************************************************/
extern int writeifhmc	(char *program, char *outfile, int *imgdim, float *voxdim, int orient, float *mmppix, float *center);

