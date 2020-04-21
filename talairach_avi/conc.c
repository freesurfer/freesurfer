/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
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
#include <math.h>
#include <float.h>
#include <string.h>
#include <ctype.h>
#include <Getifh.h>
#include <endianio.h>
#include <rec.h>
#include <conc.h>

int split (char *string, char *srgv[], int maxp) {
	int	i, m;
	char	*ptr;

	if ((ptr = strchr (string, '#'))) *ptr = '\0';
	i = m = 0;
	while (m < maxp) {
		while (!isgraph ((int) string[i]) && string[i]) i++;
		if (!string[i]) break;
		srgv[m++] = string + i;
		while (isgraph ((int) string[i])) i++;
		if (!string[i]) break;
		string[i++] = '\0';
	}
	return m;
}

/********************/
/* global variables */
/********************/

void conc_init_quiet (CONC_BLOCK *conc_block, char *program) {
	strcpy (conc_block->program, program);
	conc_block->rnfile = conc_block->wnfile = conc_block->imgfp_open = 0;
	*conc_block->lstroot = *conc_block->lstfile = *conc_block->outroot = *conc_block->outfile = '\0';
}

void conc_init (CONC_BLOCK *conc_block, char *program) {
	printf ("%s\n", "freesurfer conc.c");
	strcpy (conc_block->program, program);
	conc_block->rnfile = conc_block->wnfile = conc_block->imgfp_open = 0;
	*conc_block->lstroot = *conc_block->lstfile = *conc_block->outroot = *conc_block->outfile = '\0';
}

void conc_open_quiet (CONC_BLOCK *conc_block, char *lstfile) {
	char		*ptr, string[MAXL], *srgv[MAXL];
	int		i, k, status = 0;
	IFH		ifh;

	getroot (lstfile, conc_block->lstroot);
	sprintf (conc_block->lstfile, "%s.conc", conc_block->lstroot);
	if (!(conc_block->lstfp = fopen (conc_block->lstfile, "r"))) errr (conc_block->program, conc_block->lstfile);
	while (fgets (string, MAXL, conc_block->lstfp)) {
		if (!(ptr = strstr (string, "number_of_files:"))) continue;
		conc_block->rnfile = atoi (ptr + 16);
	}
	if (!(conc_block->imgfile0 = (char **) malloc (conc_block->rnfile * sizeof (char *)))) errm (conc_block->program);
	if (!(conc_block->nvol     = (int *)   malloc (conc_block->rnfile * sizeof (int))))    errm (conc_block->program);
	for (i = 0; i < conc_block->rnfile; i++) {
		if (!(conc_block->imgfile0[i] = (char *) malloc (MAXL * sizeof (char)))) errm (conc_block->program);
	}
	rewind (conc_block->lstfp);
	i = 0; while (fgets (string, MAXL, conc_block->lstfp)) {
		if ((ptr = strstr (string, "file:"))) {
			if (i == conc_block->rnfile) errr (conc_block->program, conc_block->lstfile);
			split (string, srgv, MAXL);
			strcpy (conc_block->imgfile0[i++], srgv[0] + 5);
		}
	}
	fclose (conc_block->lstfp);
	if (i != conc_block->rnfile) errr (conc_block->program, conc_block->lstfile);

	for (conc_block->imgdim[3] = i = 0; i < conc_block->rnfile; i++) {
		if (Getifh (conc_block->imgfile0[i], &ifh)) errr (conc_block->program, conc_block->imgfile0[i]);
		if (!i) {
			conc_block->orient = ifh.orientation;
			conc_block->isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
			for (conc_block->vdim = 1, k = 0; k < 3; k++) {
				conc_block->vdim *= conc_block->imgdim[k] = ifh.matrix_size[k];
				conc_block->voxdim[k] = ifh.scaling_factor[k];
				conc_block->mmppix[k] = ifh.mmppix[k];
				conc_block->center[k] = ifh.center[k];
			}
		} else {
			status |= (conc_block->orient != ifh.orientation);
			status |= (conc_block->isbig != strcmp (ifh.imagedata_byte_order, "littleendian"));
			for (k = 0; k < 3; k++) {
				status |= (conc_block->imgdim[k] != ifh.matrix_size[k]);
				status |= (fabs (conc_block->voxdim[k] - ifh.scaling_factor[k]) > 1.e-5);
				status |= (fabs (conc_block->mmppix[k] - ifh.mmppix[k]) > 1.e-5);
				status |= (fabs (conc_block->center[k] - ifh.center[k]) > 1.e-5);
			}
		}
		if (status) break;
		conc_block->imgdim[3] += conc_block->nvol[i] = ifh.matrix_size[3];
	}
	if (status) {
		fprintf (stderr, "%s: %s %s dimensional or byte order inconsistency\n",
			conc_block->program, conc_block->imgfile0[0], conc_block->imgfile0[i]);
		exit (-1);
	}
	conc_block->rifile = conc_block->rivol = 0;
}

void conc_open (CONC_BLOCK *conc_block, char *lstfile) {
	char		*ptr, string[MAXL], *srgv[MAXL];
	int		i, k, status = 0;
	IFH		ifh;

	getroot (lstfile, conc_block->lstroot);
	sprintf (conc_block->lstfile, "%s.conc", conc_block->lstroot);
	if (!(conc_block->lstfp = fopen (conc_block->lstfile, "r"))) errr (conc_block->program, conc_block->lstfile);
	while (fgets (string, MAXL, conc_block->lstfp)) {
		if (!(ptr = strstr (string, "number_of_files:"))) continue;
		conc_block->rnfile = atoi (ptr + 16);
	}
	printf ("number of files = %d\n", conc_block->rnfile);
	if (!(conc_block->imgfile0 = (char **) malloc (conc_block->rnfile * sizeof (char *)))) errm (conc_block->program);
	if (!(conc_block->nvol     = (int *)   malloc (conc_block->rnfile * sizeof (int))))    errm (conc_block->program);
	for (i = 0; i < conc_block->rnfile; i++) {
		if (!(conc_block->imgfile0[i] = (char *) malloc (MAXL * sizeof (char)))) errm (conc_block->program);
	}
	rewind (conc_block->lstfp);
	i = 0; while (fgets (string, MAXL, conc_block->lstfp)) {
		if ((ptr = strstr (string, "file:"))) {
			if (i == conc_block->rnfile) errr (conc_block->program, conc_block->lstfile);
			split (string, srgv, MAXL);
			strcpy (conc_block->imgfile0[i++], srgv[0] + 5);
		}
	}
	fclose (conc_block->lstfp);
	if (i != conc_block->rnfile) errr (conc_block->program, conc_block->lstfile);

	for (conc_block->imgdim[3] = i = 0; i < conc_block->rnfile; i++) {
		if (Getifh (conc_block->imgfile0[i], &ifh)) errr (conc_block->program, conc_block->imgfile0[i]);
		if (!i) {
			conc_block->orient = ifh.orientation;
			conc_block->isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
			for (conc_block->vdim = 1, k = 0; k < 3; k++) {
				conc_block->vdim *= conc_block->imgdim[k] = ifh.matrix_size[k];
				conc_block->voxdim[k] = ifh.scaling_factor[k];
				conc_block->mmppix[k] = ifh.mmppix[k];
				conc_block->center[k] = ifh.center[k];
			}
		} else {
			status |= (conc_block->orient != ifh.orientation);
			status |= (conc_block->isbig != strcmp (ifh.imagedata_byte_order, "littleendian"));
			for (k = 0; k < 3; k++) {
				status |= (conc_block->imgdim[k] != ifh.matrix_size[k]);
				status |= (fabs (conc_block->voxdim[k] - ifh.scaling_factor[k]) > 1.e-5);
				status |= (fabs (conc_block->mmppix[k] - ifh.mmppix[k]) > 1.e-5);
				status |= (fabs (conc_block->center[k] - ifh.center[k]) > 1.e-5);
			}
		}
		if (status) break;
		conc_block->imgdim[3] += conc_block->nvol[i] = ifh.matrix_size[3];
		printf ("%-5d%s%5d\n", i + 1, conc_block->imgfile0[i], conc_block->nvol[i]);
	}
	if (status) {
		fprintf (stderr, "%s: %s %s dimensional or byte order inconsistency\n",
			conc_block->program, conc_block->imgfile0[0], conc_block->imgfile0[i]);
		exit (-1);
	}
	printf ("orient   %d\n", conc_block->orient);
	printf ("imgdim %10d%10d%10d%10d\n",
			conc_block->imgdim[0], conc_block->imgdim[1], conc_block->imgdim[2], conc_block->imgdim[3]);
	printf ("voxdim %10.6f%10.6f%10.6f\n", conc_block->voxdim[0], conc_block->voxdim[1], conc_block->voxdim[2]);
	printf ("mmppix %10.6f%10.6f%10.6f\n", conc_block->mmppix[0], conc_block->mmppix[1], conc_block->mmppix[2]);
	printf ("center %10.4f%10.4f%10.4f\n", conc_block->center[0], conc_block->center[1], conc_block->center[2]);
	conc_block->rifile = conc_block->rivol = 0;
}

void conc_rewind (CONC_BLOCK *conc_block) {
	if (!strlen (conc_block->lstfile) || !conc_block->rnfile) {
		fprintf (stderr, "%s: conc_rewind ignored as read not initialized\n", conc_block->program);
		return;
	}
	if (conc_block->imgfp_open) fclose (conc_block->imgfp);
	conc_block->rivol = 0;
	conc_block->rifile = 0;
}

void conc_read_vol (CONC_BLOCK *conc_block, float *imgt) {
	if (!strlen (conc_block->lstfile) || !conc_block->rnfile) {
		fprintf (stderr, "%s: conc read not initialized\n", conc_block->program);
		exit (-1);
	}
	if (!conc_block->rivol) {
		if (!(conc_block->imgfp = fopen (conc_block->imgfile0[conc_block->rifile], "rb")))
			errr (conc_block->program, conc_block->imgfile0[conc_block->rifile]);
		conc_block->imgfp_open++;
	}
	if (eread (imgt, conc_block->vdim, conc_block->isbig, conc_block->imgfp))
		errr (conc_block->program, conc_block->imgfile0[conc_block->rifile]);
	if (++conc_block->rivol == conc_block->nvol[conc_block->rifile]) {
		fclose (conc_block->imgfp);
		conc_block->rivol = conc_block->imgfp_open = 0;
		if (++conc_block->rifile == conc_block->rnfile) conc_block->rifile = 0;
	}
}

void conc_write_vol (CONC_BLOCK *conc_block, float *imgt) {
	if (!strlen (conc_block->outfile)) {
		fprintf (stderr, "%s: conc output not initialized\n", conc_block->program);
		exit (-1);
	}
	if (!conc_block->wivol) {
		if (!(conc_block->outfp = fopen (conc_block->imgfile1[conc_block->wifile], "wb")))
		      errw (conc_block->program, conc_block->imgfile1[conc_block->wifile]);
	}
	if (ewrite (imgt, conc_block->vdim, conc_block->control, conc_block->outfp))
		errw (conc_block->program, conc_block->imgfile1[conc_block->wifile]);
	if (++conc_block->wivol == conc_block->nvol[conc_block->wifile]) {
		fclose (conc_block->outfp);
		conc_block->wivol = 0;
		if (++conc_block->wifile == conc_block->wnfile) conc_block->wifile = 0;
	}
}

void conc_newe (CONC_BLOCK *conc_block, char *trailer, char control) {
	conc_new (conc_block, trailer);
	conc_block->control = control;
}

void conc_new (CONC_BLOCK *conc_block, char *trailer) {
	char	imgroot[MAXL];
	int		i;

	if (!strlen (conc_block->lstfile) || !conc_block->rnfile) {
		fprintf (stderr, "%s: conc_new not possible because conc read not initialized\n", conc_block->program);
		exit (-1);
	}
	conc_block->wnfile = conc_block->rnfile;
	if (!(conc_block->imgfile1 = (char **) malloc (conc_block->wnfile * sizeof (char *)))) errm (conc_block->program);
	sprintf (conc_block->outfile, "%s_%s.conc", conc_block->lstroot, trailer);
	if (!(conc_block->outfp = fopen (conc_block->outfile, "w"))) errw (conc_block->program, conc_block->outfile);
	fprintf (conc_block->outfp, "number_of_files: %d\n", conc_block->wnfile);
	for (i = 0; i < conc_block->wnfile; i++) {
		if (!(conc_block->imgfile1[i] = (char *) malloc (MAXL * sizeof (char)))) errm (conc_block->program);
		getroot (conc_block->imgfile0[i], imgroot);
		sprintf (conc_block->imgfile1[i], "%s_%s.4dfp.img", imgroot, trailer);
		fprintf (conc_block->outfp, "\tfile:%s\n", conc_block->imgfile1[i]);
	}
	if (fclose (conc_block->outfp)) errw (conc_block->program, conc_block->outfile);
	conc_block->wifile = conc_block->wivol = 0;
/***************************************************************/
/* default output input-endian overridden by using conc_newe() */
/***************************************************************/
	conc_block->control = (conc_block->isbig) ? 'b' : 'l';
}

int conc_ifh_hdr_rec (CONC_BLOCK *conc_block, int argc, char *argv[], char *recstr) {
	char	command[MAXL];
	int		i, status = 0;
	int		imgdim[4];

	if (!conc_block->wnfile) {
		fprintf (stderr, "%s: conc output not initialized\n", conc_block->program);
		exit (-1);
	}

	switch (conc_block->control) {
		case 'b': case 'B': conc_block->osbig = 1; break;
		case 'l': case 'L': conc_block->osbig = 0; break;
		default: conc_block->osbig = CPU_is_bigendian(); break;
	}

	for (i = 0; i < 3; i++) imgdim[i] = conc_block->imgdim[i];
	for (i = 0; i < conc_block->wnfile; i++) {
		imgdim[3] = conc_block->nvol[i];
		writeifhmce (conc_block->program, conc_block->imgfile1[i], imgdim, conc_block->voxdim,
			conc_block->orient, conc_block->mmppix, conc_block->center, conc_block->control);
		sprintf (command, "ifh2hdr %s", conc_block->imgfile1[i]); status |= system (command);
		startrece (conc_block->imgfile1[i], argc, argv, recstr, conc_block->control);
		catrec (conc_block->imgfile0[i]);
		endrec ();
	}
	return status;
}

void conc_free (CONC_BLOCK *conc_block) {
	int	i;

	if (conc_block->rnfile) {
		for (i = 0; i < conc_block->rnfile; i++) free (conc_block->imgfile0[i]);
		free (conc_block->imgfile0);
		free (conc_block->nvol);
	}
	if (conc_block->wnfile) {
		for (i = 0; i < conc_block->wnfile; i++) free (conc_block->imgfile1[i]);
		free (conc_block->imgfile1);
	}
}
