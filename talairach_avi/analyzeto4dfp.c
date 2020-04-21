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
#include <string.h>
#include <math.h>
#include <float.h>
#include <endianio.h>	/* includes <ANALYZE.h>	*/
#include <Getifh.h>
#include <rec.h>

#define MAXL	256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

extern void flipx (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void flipy (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */
extern void flipz (float *imgf, int *pnx, int* pny, int *pnz);	/* cflip.c */

int main (int argc, char *argv[]) {
	FILE		*fpimg, *fpout;
	struct dsr	hdr;					/* ANALYZE hdr */
	char		*str, command[MAXL], program[MAXL];
	char		imgroot[MAXL], outroot[MAXL];
	char		imgfile[MAXL], outfile[MAXL];

/****************/
/* image arrays */
/****************/
	float		fmin = +FLT_MAX;
	float		fmax = -FLT_MAX;
	float		*imgf;
	unsigned char	*imgu;
	short		*imgi;
	float		ROIScaleFactor=0.0;
	float		voxsiz[3];
	int		imgdim[4], vdim, bytepix, orient=0;
	char		control = '\0';

/***********/
/* utility */
/***********/
	int 		c, i, k;

/*********/
/* flags */
/*********/
	int		status = 0;
	int		Weber_flag = 0;
	int		swab_flag = 0;
	int		scale_flag = 0;
	int		O_flag = 0;
	int		xflag = 0, yflag = 0, zflag = 0;

	fprintf (stdout, "%s\n", "freesurfer analyzeto4dfp.c");
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while ((c = *str++)) switch (c) {
				case 's': scale_flag++;				break;
				case 'x': xflag++;				break;
				case 'y': yflag++;				break;
				case 'z': zflag++;				break;
				case '@': control = *str++;			*str = '\0'; break;
				case 'O': orient = atoi (str++); O_flag++;	*str = '\0'; break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
		}
	}	
	if (k < 1) {
		fprintf (stderr, "Usage: %s <analyze_image>\n", program);
		fprintf (stderr, "\toption\n");
		fprintf (stderr, "\t-s	apply ROIScaleFactor\n");
		fprintf (stderr, "\t-x	flip first  axis\n");
		fprintf (stderr, "\t-y	flip second axis\n");
		fprintf (stderr, "\t-z	flip third  axis\n");
		fprintf (stderr, "\t-@<b|l>	toutput big or little endian (default CPU endian)\n");
		fprintf (stderr, "\t-O<int>	supply orientation code (in range [0-5])\n");
		exit (1);
	}

	strcpy (outroot, imgroot);
	while ((str = strrchr (outroot, '.'))) {
		if (!strcmp (str, ".4dint")) *str = '\0';
		else break;
	}

 	sprintf (imgfile, "%s.hdr", imgroot);
	if (!(fpimg = fopen (imgfile, "rb"))
	|| fread (&hdr, sizeof (struct dsr), 1, fpimg) != 1
	|| fclose (fpimg)) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);
	if (hdr.hk.sizeof_hdr != sizeof (struct dsr)) {
		printf ("converting byte order\n");
		swab_hdr (&hdr);
		swab_flag++;
	}
	printf ("header size %d bytes\n", hdr.hk.sizeof_hdr);
	i = (unsigned long) &hdr.dime.datatype	- (unsigned long) &hdr;
	printf ("hdr.dime.datatype\toffset=%d	value=%d\n", i, hdr.dime.datatype);
	i = (unsigned long) &hdr.dime.bitpix	- (unsigned long) &hdr;
	printf ("hdr.dime.bitpix\t\toffset=%d	value=%d\n", i,	hdr.dime.bitpix);

	k = (hdr.dime.datatype ==  2 && hdr.dime.bitpix == 8)	/* unsigned char */
	||  (hdr.dime.datatype ==  4 && hdr.dime.bitpix == 16)	/* signed short int */
	||  (hdr.dime.datatype == 16 && hdr.dime.bitpix == 32);	/* float */
	if (!k) {
		fprintf (stderr, "%s: unrecognized data type %s\n", program, imgfile);
		exit (-1);
	}
	bytepix = hdr.dime.bitpix / 8;

/*******************************/
/* determine image orientation */
/*******************************/
	i = (unsigned long) &hdr.hist.orient	- (unsigned long) &hdr;
	printf ("hdr.hist.orient\t\toffset=%d	value=%d\n", i,	hdr.hist.orient);
	if (O_flag) {
		hdr.hist.orient = orient;
	} else {
		orient = hdr.hist.orient;
	}
	if (orient < 0	|| orient > 5) {
		fprintf (stderr, "%s: unrecognized orientation (%d)\n", program, orient);
		exit (-1);
	}
	if (orient > 2 && orient < 6) {		/* unflip according to Darren Weber */
		orient -= 3; hdr.hist.orient = orient;
		Weber_flag = 1;
	}

	printf ("dimensionality%6d\n", hdr.dime.dim[0]);
	if (hdr.dime.dim[0] < 3) {
		fprintf (stderr, "input image dimensionality (%d) must be at least 3\n", hdr.dime.dim[0]);
		exit (-1);
	}
	if (hdr.dime.dim[0] == 3) hdr.dime.dim[4] = 1;
	printf ("dimensions");
	for (i = 0; i < 4; i++) printf ("%10d", hdr.dime.dim[i + 1]); printf ("\n");
	for (i = 0; i < 4; i++) imgdim[i] =  hdr.dime.dim[i + 1];
	for (i = 0; i < 3; i++) voxsiz[i] =  hdr.dime.pixdim[i + 1];
	vdim = imgdim[0] * imgdim[1] * imgdim[2];

	if (!(imgu = (unsigned char *) malloc (vdim * bytepix))) errm (program);
	if (!(imgi = 	     (short *) malloc (vdim * bytepix))) errm (program);
	if (!(imgf = (float *) malloc (vdim * sizeof (float))))  errm (program);

/***************************/
/* determine output endian */
/***************************/
	if (!control) control = (CPU_is_bigendian()) ? 'b' : 'l';

/*************************/
/* open input and output */
/*************************/
 	sprintf (imgfile, "%s.img", imgroot);
	if (!(fpimg = fopen (imgfile, "rb"))) errr (program, imgfile);
	printf ("Reading: %s\n", imgfile);
 	sprintf (outfile, "%s.4dfp.img", outroot);
	if (!(fpout = fopen (outfile, "wb"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);

/***********************/
/* process all volumes */
/***********************/
	for (k = 0; k < imgdim[3]; k++) {
		switch (hdr.dime.datatype) {
		case 16:
			if (fread (imgf, bytepix, vdim, fpimg) != vdim) errr (program, imgfile);
			if (swab_flag) for (i = 0; i < vdim; i++) swab4 ((char *) &imgf[i]);
			break;
		case 4:
			if (fread (imgi, bytepix, vdim, fpimg) != vdim) errr (program, imgfile);
			if (swab_flag) for (i = 0; i < vdim; i++) swab2 ((char *) &imgi[i]);
			for (i = 0; i < vdim; i++) imgf[i] = (float) imgi[i];
			break;
		case 2:
			if (fread (imgu, bytepix, vdim, fpimg) != vdim) errr (program, imgfile);
			for (i = 0; i < vdim; i++) imgf[i] = (float) imgu[i];
			break;
		}
		if (scale_flag) {
			ROIScaleFactor = hdr.dime.funused9;
			hdr.dime.funused9 = 0.;
		}
		if (Weber_flag) flipy (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		switch (orient) {	/* 0 = transverse; 1 = coronal; 2 = sagittal */
			case 2:	flipx (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
			case 1:	flipz (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
			case 0:	flipy (imgf, imgdim + 0, imgdim + 1, imgdim + 2);	break;
			default:							break;
		}
		if (xflag)	flipx (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		if (yflag)	flipy (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		if (zflag)	flipz (imgf, imgdim + 0, imgdim + 1, imgdim + 2);
		for (i = 0; i < vdim; i++) {
			if (scale_flag) imgf[i] *= ROIScaleFactor;
			if (imgf[i] > fmax) fmax = imgf[i];
			if (imgf[i] < fmin) fmin = imgf[i];
		}
		if (ewrite (imgf, vdim, control, fpout)) errw (program, outfile);
	}
	fclose (fpimg);
	fclose (fpout);

/**************/
/* create ifh */
/**************/
	writeifhe (program, outfile, imgdim, voxsiz, orient + 2, control);

/*******/
/* hdr */
/*******/
	sprintf (command, "ifh2hdr %s -r%dto%d", outroot, (int) (fmin - 0.5), (int) (fmax + 0.5));
	printf ("%s\n", command); status = system (command);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, "freesurfer analyzeto4dfp.c", control);
	if (swab_flag) {
		sprintf (command, "Byte order swapped\n");
		printrec (command);
	}
	if (Weber_flag) {
		sprintf (command, "Orientation of second axis flipped according to Fred Weber\n");
		printrec (command);
	}
	if (xflag) {
		sprintf (command, "First  axis flipped\n");
		printrec (command);
	}
	if (yflag) {
		sprintf (command, "Second axis flipped\n");
		printrec (command);
	}
	if (zflag) {
		sprintf (command, "Third  axis flipped\n");
		printrec (command);
	}
	if (scale_flag) {
		sprintf (command, "Voxel values multiplied by ROIScaleFactor=%.6e\n", ROIScaleFactor);
		printrec (command);
	}
	catrec   (imgfile);
	endrec ();

	free (imgu);
	free (imgi);
	free (imgf);
	exit (status);
}
