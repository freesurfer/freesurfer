/**
 * @file  endianio.c
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/01/04 03:17:25 $
 *    $Revision: 1.3.6.1 $
 *
 * Copyright 1999 - 2011
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ANALYZE.h>
#include <Getifh.h>

#define MAXL		256

static char rcsid[] = "$Id: endianio.c,v 1.3.6.1 2012/01/04 03:17:25 nicks Exp $";
void endianio_rcs (void) {printf ("%s\n", rcsid);}

void swab2 (char *a) {
	char	t;

	t = a[0]; a[0] = a[1]; a[1] = t;
}

void swab4 (char *a) {
	char	t;

	t = a[0]; a[0] = a[3]; a[3] = t;
	t = a[1]; a[1] = a[2]; a[2] = t;
}

void swab_hdr (struct dsr *phdr) {
	int	i;
	float	*fptr;

	swab4 ((char *) &phdr->hk.sizeof_hdr);
	swab4 ((char *) &phdr->hk.extents);
	swab2 ((char *) &phdr->hk.session_error);
	for (i = 0; i < 8; i++) swab2 ((char *) &phdr->dime.dim[i]);
	swab2 ((char *) &phdr->dime.datatype);
	swab2 ((char *) &phdr->dime.bitpix);
	swab2 ((char *) &phdr->dime.dim_un0);
	for (i = 0; i < 8; i++) swab4 ((char *) &phdr->dime.pixdim[i]);
	fptr = &phdr->dime.funused8; for (i = 0; i < 8; i++) swab4 ((char *) (fptr + i));
	swab4 ((char *) &phdr->dime.compressed);
	swab4 ((char *) &phdr->dime.verified);
	swab4 ((char *) &phdr->dime.glmax);
	swab4 ((char *) &phdr->dime.glmin);
	swab4 ((char *) &phdr->hist.views);
	swab4 ((char *) &phdr->hist.vols_added);
	swab4 ((char *) &phdr->hist.start_field);
	swab4 ((char *) &phdr->hist.field_skip);
	swab4 ((char *) &phdr->hist.omax);
	swab4 ((char *) &phdr->hist.omin);
	swab4 ((char *) &phdr->hist.smax);
	swab4 ((char *) &phdr->hist.smin);
}

int CPU_is_bigendian () {
	union {
		float		q;
		unsigned char	c[4];
	} fone;
	fone.q = 1.0;
	return fone.c[0];
}

void errf (char* program) {
	fprintf (stderr, "%s: illegal bytes/word\n", program);
	exit (-1);
}

void errm (char* program) {
	fprintf (stderr, "%s: memory allocation error\n", program);
	exit (-1);
}

void errr (char* program, char* filespc) {
	fprintf (stderr, "%s: %s read error\n", program, filespc);
	exit (-1);
}

void errw (char* program, char* filespc) {
	fprintf (stderr, "%s: %s write error\n", program, filespc);
	exit (-1);
}

void getroot (char *filespc, char *imgroot) {
	char	*str;
	strcpy (imgroot, filespc);
	while ((str = strrchr (imgroot, '.'))) {
			if (!strcmp (str, ".rec"))	*str = '\0';
		else	if (!strcmp (str, ".img"))	*str = '\0';
		else	if (!strcmp (str, ".ifh"))	*str = '\0';
		else	if (!strcmp (str, ".4dfp"))	*str = '\0';
		else	if (!strcmp (str, ".hdr"))	*str = '\0';
		else	if (!strcmp (str, ".conc"))	*str = '\0';
		else	break;
	}
}

int eread (float *imgt, int n, int isbig, FILE *fp) {
	int	i, swab_flag;

	swab_flag = (CPU_is_bigendian() != 0) != (isbig != 0);
	if (0) printf ("eread swab_flag=%d\n", swab_flag);
	if (fread (imgt, sizeof (float), n, fp) != n) return -1;
	if (swab_flag) for (i = 0; i < n; i++) swab4 ((char *) (imgt + i));

	return 0;
}

int ewrite (float *imgt, int n, char control, FILE *fp) {
	int	i, swab_flag;
	float	f;

	swab_flag = ((CPU_is_bigendian() != 0) && (control == 'l' || control == 'L'))
		 || ((CPU_is_bigendian() == 0) && (control == 'b' || control == 'B'));

	if (0) printf ("ewrite swab_flag=%d\n", swab_flag);
	if (swab_flag) {
		for (i = 0; i < n; i++) {
			f = imgt[i];
			swab4 ((char *) &f);
			if (fwrite (&f, sizeof (float), 1, fp) != 1) return -1;
		}
		return 0;
	} else {
		return (fwrite (imgt, sizeof (float), n, fp) != n);
	}
}

int gread (char *imgt, size_t bytes, int n, FILE *fp, int isbig) {
	int	i, swab_flag;

	swab_flag = (CPU_is_bigendian() != 0) != (isbig != 0);
	if (0) printf ("gread swab_flag=%d\n", swab_flag);

	if (fread (imgt, bytes, n, fp) != n) return -1;
	if (swab_flag) for (i = 0; i < n; i++) switch (bytes) {
		case 2: swab2 (imgt + 2*i); break;
		case 4: swab4 (imgt + 4*i); break;
		default: errf ("gread");
	}
	return 0;
}

int gwrite (char *imgt, size_t bytes, int n, FILE *fp, char control) {
	int	i, swab_flag;
	char	*imgb;		/* i/o buffer */
	int	status;

	swab_flag = ((CPU_is_bigendian() != 0) && (control == 'l' || control == 'L'))
		 || ((CPU_is_bigendian() == 0) && (control == 'b' || control == 'B'));
	if (0) printf ("gwrite swab_flag=%d\n", swab_flag);

	status = 0;
	if (swab_flag) {
		if (!(imgb = malloc (bytes*n))) errm ("gwrite");
		for (i = 0; i < bytes*n; i++) imgb[i] = imgt[i];
		for (i = 0; i < n; i++) switch (bytes) {
			case 2:	swab2 (imgb + 2*i); break;
			case 4:	swab4 (imgb + 4*i); break;
			default: errf ("gwrite");
		}
	} else {
		imgb = imgt;
	}
	status = fwrite (imgb, bytes, n, fp) != n;
	if (swab_flag) free (imgb);
	return status;
}

void load_4dfp_frame (char *fileroot, int *imgdim, int frame, int isbig, float *fimg) {
	FILE		*fp;
	static char	subr[] = "load_4dfp_frame";
	char		filespc[4*MAXL];
	int		vdim, i, swab_flag;

	vdim = imgdim[0]*imgdim[1]*imgdim[2];

	sprintf (filespc, "%s.4dfp.img", fileroot);
	fprintf (stdout, "Reading: %s frame %d\n", filespc, frame + 1);
	if (!(fp = fopen (filespc, "rb")) || fseek (fp, (long) frame * vdim * sizeof (float), SEEK_SET)
	|| fread (fimg, sizeof (float), vdim, fp) != vdim || fclose (fp)) errr (subr, filespc);

	swab_flag = (CPU_is_bigendian() != 0) != (isbig != 0);
	if (swab_flag) for (i = 0; i < vdim; i++) swab4 ((char *) (fimg + i));
}

int get_4dfp_dimoe (char *fileroot, int *imgdim, float *voxsiz, int *orient, int *isbig) {
	IFH		ifh;
	static char	subr[] = "get_4dfp_dimoe";
	char		filespc[4*MAXL];
	int		i, k, status;
	char		*TCS[3] = {"T", "C", "S"};

	getroot (fileroot, filespc);
	if ((i = strlen (filespc)) + 10 > 4*MAXL) {
		fprintf (stdout, "%s: %s filename too long\n", subr, fileroot);
		return -1;
	}

	status = 0;
	strcat (filespc, ".4dfp.img");
	if (Getifh (filespc, &ifh)) errr (subr, filespc);
	for (k = 0; k < 3; k++) voxsiz[k] = ifh.scaling_factor[k];
	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];

	fprintf (stdout, "%s\n", fileroot);
	fprintf (stdout, "%10d%10d%10d%10d\n", imgdim[0], imgdim[1], imgdim[2], imgdim[3]);
	fprintf (stdout, "%10f%10f%10f\n", voxsiz[0], voxsiz[1], voxsiz[2]);
	*orient = ifh.orientation;
	if (*orient < 2 || *orient > 4) {
		fprintf (stderr, "%s warning: %s illegal orientation (%d)\n", subr, fileroot, ifh.orientation);
		status = -2;
	} else {
		fprintf (stdout, "orientation %s byte_order %s\n", TCS[*orient - 2], ifh.imagedata_byte_order);
	}
	*isbig = !strcmp (ifh.imagedata_byte_order, "bigendian");
	return status;
}

int get_4dfp_dimoe_quiet (char *fileroot, int *imgdim, float *voxsiz, int *orient, int *isbig) {
	IFH		ifh;
	static char	subr[] = "get_4dfp_dimoe_quiet";
	char		filespc[4*MAXL];
	int		i, k, status;

	getroot (fileroot, filespc);
	if ((i = strlen (filespc)) + 10 > 4*MAXL) {
		fprintf (stderr, "%s: %s filename too long\n", subr, fileroot);
		return -1;
	}

	status = 0;
	strcat (filespc, ".4dfp.img");
	if (Getifh (filespc, &ifh)) errr (subr, filespc);
	for (k = 0; k < 3; k++) voxsiz[k] = ifh.scaling_factor[k];
	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];

	*orient = ifh.orientation;
	if (*orient < 2 || *orient > 4) {
		fprintf (stderr, "%s warning: %s illegal orientation (%d)\n", subr, fileroot, ifh.orientation);
		status = -2;
	}
	*isbig = !strcmp (ifh.imagedata_byte_order, "bigendian");
	return status;
}

#define N	1000000
void eread_ewrite_test () {
	FILE	*fp;
	char	program[] = "eread_ewrite_test", filespc[] = "/data/petsun23/testx";
	int	i, n = N;
	float	*x, *y;

	printf ("CPU_is_bigendian=%d\n", CPU_is_bigendian ());

	if (!(x = (float *) malloc (n * sizeof (float)))) errm (program);
	if (!(y = (float *) malloc (n * sizeof (float)))) errm (program);
	for (i = 0; i < n; i++) {
		x[i] = -1. + 2.0*rand ()/(float) RAND_MAX;
		if (i < 20) printf ("%10d%10.4f\n", i, x[i]);
	}

	printf ("Writing: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "wb")))		errw (program, filespc);
	if (ewrite (x, n, 'b', fp) || fclose (fp))	errw (program, filespc);
	printf ("Reading: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "rb")))		errr (program, filespc);
	if (eread (y, n, 1, fp) || fclose (fp))		errr (program, filespc);
	for (i = 0; i < n; i++) assert (y[i] == x[i]);

	printf ("Writing: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "wb")))		errw (program, filespc);
	if (ewrite (x, n, 'l', fp) || fclose (fp))	errw (program, filespc);
	printf ("Reading: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "rb")))		errr (program, filespc);
	if (eread (y, n, 0, fp) || fclose (fp))		errr (program, filespc);
	for (i = 0; i < n; i++) assert (y[i] == x[i]);

	printf ("Done\n"); fflush (stdout);
	free (x); free (y);
}

/* main (int argc, char *argv[]) {	*/
void gread_gwrite_test () {
	FILE	*fp;
	char	program[] = "gread_ewrite_test", filespc[] = "/data/petsun23/testx";
	int	i, n = N;
	float	*x, *y;

	printf ("CPU_is_bigendian=%d\n", CPU_is_bigendian ());

	if (!(x = (float *) malloc (n * sizeof (float)))) errm (program);
	if (!(y = (float *) malloc (n * sizeof (float)))) errm (program);
	for (i = 0; i < n; i++) {
		x[i] = -1. + 2.0*rand ()/(float) RAND_MAX;
		if (i < 20) printf ("%10d%10.4f\n", i, x[i]);
	}

	printf ("Writing: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "wb")))					errw (program, filespc);
	if (gwrite ((char *) x, sizeof (float), n, fp, 'b') || fclose (fp))	errw (program, filespc);
	printf ("Reading: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "rb")))					errr (program, filespc);
	if (gread  ((char *) y, sizeof (float), n, fp, 1) || fclose (fp))	errr (program, filespc);
	if (0) for (i = 0; i < n; i++) printf ("%10f%10f\n", x[i], y[i]);
	for (i = 0; i < n; i++) assert (y[i] == x[i]);

	printf ("Writing: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "wb")))					errw (program, filespc);
	if (gwrite ((char *) x, sizeof (float), n, fp, 'l') || fclose (fp))	errw (program, filespc);
	printf ("Reading: %s\n", filespc); fflush (stdout);
	if (!(fp = fopen (filespc, "rb")))					errr (program, filespc);
	if (gread  ((char *) y, sizeof (float), n, fp, 0) || fclose (fp))	errr (program, filespc);
	if (0) for (i = 0; i < n; i++) printf ("%10f%10f\n", x[i], y[i]);
	for (i = 0; i < n; i++) assert (y[i] == x[i]);

	printf ("Done\n"); fflush (stdout);
	free (x); free (y);
}
