/**
 * @file  ifh2hdr.c
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL	256

void getrange (char *string, float *minval, float *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atof (string);
		*maxval = atof (str + 2);
	} else {
		*minval = 0.0;
		*maxval = atof (string);
	}
}

extern int	Inithdr (struct dsr *phdr, int *imgdim, float *voxdim, char *proto_imgfile);

static char rcsid[] = "$Id: ifh2hdr.c,v 1.2 2007/05/05 00:00:06 nicks Exp $";
int main (int argc, char *argv[]) {
	FILE		*fp;
	struct dsr	hdr;
	IFH		ifh;

	char		filespc[MAXL], imgroot[MAXL];
	float		voxsiz[3];
	float           fmin, fmax;
	int		imgdim[4];

/***********/
/* utility */
/***********/
	char		*str, command[MAXL], program[MAXL];
	int 		c, i, k;

/*********/
/* flags */
/*********/
	int		verbose = 0;

	int		isbig, swab_flag, range_flag = 0;

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while ((c = *str++)) switch (c) {
				case 'r':	getrange (str, &fmin, &fmax);
						*str = '\0'; range_flag++;		break;
				case 'v':	verbose++;				break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) file>\n", program);
		printf ("e.g.,\t%s vc654_mpr_atl -r-500to1500\n", program);
		printf ("\toption\n");
		printf ("\t-r<flt>[to<flt>]\tset range\n");
		printf ("N.B.:\t%s preserves the endian state of the input ifh\n", program);
		exit (1);
	}

	sprintf (filespc, "%s.4dfp.ifh", imgroot);
	if (Getifh (filespc, &ifh)) errr (program, filespc);
	for (k = 0; k < 3; k++) voxsiz[k] = ifh.scaling_factor[k];
	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];

	Inithdr (&hdr, imgdim, voxsiz, "");
	hdr.dime.datatype = 16;			/* float */
	hdr.dime.bitpix = 32;
	hdr.hist.orient = ifh.orientation - 2;
	if (range_flag) {
		hdr.dime.glmin = fmin;
		hdr.dime.glmax = fmax;
	}

	isbig = !strcmp (ifh.imagedata_byte_order, "bigendian");
	swab_flag = (CPU_is_bigendian() != 0) != (isbig != 0);
	if (swab_flag) swab_hdr (&hdr);

	sprintf (filespc, "%s.4dfp.hdr", imgroot);
	printf ("Writing: %s\n", filespc);
	if (!(fp = fopen (filespc, "wb"))
	|| fwrite (&hdr, sizeof (struct dsr), 1, fp) != 1) errw (program, filespc);
	fclose (fp);

	exit (0);
}
