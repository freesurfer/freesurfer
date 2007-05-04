/*$Header: /space/repo/1/dev/dev/talairach_avi/zero_lt_4dfp.c,v 1.1 2007/05/04 22:34:03 nicks Exp $*/
/*$Log: zero_lt_4dfp.c,v $
/*Revision 1.1  2007/05/04 22:34:03  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.11  2007/05/01  04:39:58  avi
 * Solaris10 and gccc v3 compliant
 *
 * Revision 1.10  2005/10/25  03:32:15  avi
 * allow first argumnet to be negative (not pseudo-option)
 *
 * Revision 1.9  2004/11/05  22:25:09  rsachs
 * Removed 'Get4dfpDimN'. Installed 'errm,'errr','errw','getroot','get_4dfp_dimo'.
 *
 * Revision 1.8  2001/08/02  01:11:51  avi
 * correct usage
 *
 * Revision 1.7  2001/08/02  00:39:42  avi
 * code standardization
 *
 * Revision 1.6  1998/12/13  00:51:22  avi
 * correct error computing dimension
 *
 * Revision 1.5  1998/12/13  00:18:57  avi
 * correct usage
 *
 * Revision 1.4  1998/10/12  22:15:07  mcavoy
 * Revision 1.2  1998/10/12  21:57:58  mcavoy
 * Revision 1.1  1998/05/08  22:32:04  tscull
 * Revision 1.2  1997/11/26  20:04:53  tscull
 * Initial revision
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <unistd.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>

#define MAXL 256

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

static char rcsid[]= "$Id: zero_lt_4dfp.c,v 1.1 2007/05/04 22:34:03 nicks Exp $";
int main (int argc, char **argv) {
/*************/
/* image I/O */
/*************/
	FILE            *fp_img, *fp_out;
	IFH		ifh;
	char            imgfile[MAXL], imgroot[MAXL];
	char            outfile[MAXL], outroot[MAXL] = "";

/**************/
/* processing */
/**************/
	int             imgdim[4], dimension, orient, isbig;
	float           voxdim[3];
	float           *imgr;
	float		thresh=0.0;
	char		control ='\0';

/***********/
/* utility */
/***********/
	int             c, i, k;
	char            *ptr, command[MAXL], program[MAXL];

/*********/
/* flags */
/*********/
        int             status = 0;

	printf ("%s\n", rcsid);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (i > 1 && *argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while ((c = *ptr++)) switch (c) {
				case '@': control = *ptr++;		*ptr = '\0'; break;
			}
		} else switch (k) {
			case 0: thresh = atof (argv[i]);    k++;	break;
			case 1: getroot (argv[i], imgroot); k++;	break;
			case 2: getroot (argv[i], outroot); k++;	break;
		}
	}
	if (k < 2) {
		printf ("Usage:	%s <flt> <file_4dfp> [outroot]\n", program);
		printf ("	option\n");
		printf ("	-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("e.g.,	%s 90 pt349_study9to9\n", program);
		printf ("e.g.,	%s 90 pt349_study9to9 pt349_study9to9z\n", program);
		printf ("N.B.:	default output 4dfp root is <file_4dfp>\"z\"\n");
		exit (1);
	}

	sprintf (imgfile, "%s.4dfp.img", imgroot);
/***************************************/
/* create output filename if not given */
/***************************************/
	if (!strlen (outroot)) sprintf (outroot, "%sz", imgroot);
	sprintf (outfile, "%s.4dfp.img", outroot);

/*****************************/
/* get 4dfp input dimensions */
/*****************************/
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgroot);
	if (Getifh (imgfile, &ifh)) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];

/*****************/
/* alloc buffers */
/*****************/
	if (!(imgr = (float *) malloc (dimension * sizeof (float)))) errm (program);

/***********/
/* process */
/***********/
	if (!(fp_img = fopen (imgfile, "rb"))) errr (program, imgfile);
	if (!(fp_out = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);
	for (k = 0; k < imgdim[3]; k++) {
		if (gread  ((char *) imgr, sizeof (float), dimension, fp_img, isbig)) errr (program, imgfile);
		for (i = 0; i < dimension; i++) if (imgr[i] < thresh) imgr[i] = 0.0;
		if (gwrite ((char *) imgr, sizeof (float), dimension, fp_out, control)) errw (program, outfile);
	}
	fclose (fp_img);
	fclose (fp_out);

/***************/
/* ifh and hdr */
/***************/
	if (Writeifh (program, outfile, &ifh, control)) errw (program, outfile);
	sprintf (command, "ifh2hdr %s", outroot); printf ("%s\n", command);
	status = system (command);

/*******************/
/* create rec file */
/*******************/
	startrece (outfile, argc, argv, rcsid, control);
	catrec    (imgfile);
        endrec    ();

        free (imgr);
	exit (status);
}
