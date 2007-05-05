/**
 * @file  conc.h
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
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

#define MAXL	256

typedef struct {
	char	program[MAXL];
	char	lstroot[MAXL], lstfile[MAXL];
	char	outroot[MAXL], outfile[MAXL];
	char	control;
	int	imgdim[4], vdim, orient, isbig, osbig;
	float	voxdim[3], mmppix[3], center[3];
	FILE	*lstfp, *imgfp, *outfp;
	char	**imgfile0;		/* input  4dfp filenames */
	char	**imgfile1;		/* output 4dfp filenames */
	int	rivol, wivol, *nvol;
	int	rifile, wifile, rnfile, wnfile, imgfp_open;
} CONC_BLOCK;

int split		(char *string, char *srgv[], int maxp);
void conc_init_quiet	(CONC_BLOCK *conc_block, char *program);
void conc_init		(CONC_BLOCK *conc_block, char *program);
void conc_open_quiet	(CONC_BLOCK *conc_block, char *lstfile);
void conc_open		(CONC_BLOCK *conc_block, char *lstfile);
void conc_read_vol	(CONC_BLOCK *conc_block, float *imgt);
void conc_write_vol	(CONC_BLOCK *conc_block, float *imgt);
void conc_new		(CONC_BLOCK *conc_block, char *trailer);
void conc_newe		(CONC_BLOCK *conc_block, char *trailer, char control);
int conc_ifh_hdr_rec	(CONC_BLOCK *conc_block, int argc, char *argv[], char *recstr);
void conc_free		(CONC_BLOCK *conc_block);
void conc_rewind	(CONC_BLOCK *conc_block);
