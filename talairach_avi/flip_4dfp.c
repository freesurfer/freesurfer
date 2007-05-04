/*$Header: /space/repo/1/dev/dev/talairach_avi/flip_4dfp.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: flip_4dfp.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.5  2007/05/02  01:28:30  avi
 * endian gcc v3 compliant
 *
 * Revision 1.4  2004/11/16  22:26:42  rsachs
 * Installed 'setprog'. Replaced 'Get4dfpDimN' with 'get_4dfp_dimo'.
 *
 * Revision 1.3  2004/02/19  01:04:23  avi
 * eliminate calls to FORTRAN dependent flip routines
 *
 * Revision 1.2  1999/02/03  06:40:14  avi
 * correct output name generation
 *
 * Revision 1.1  1999/01/24  07:07:06  avi
 * Initial revision
 **/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <endianio.h>
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

static char rcsid[] = "$Id: flip_4dfp.c,v 1.1 2007/05/04 22:33:59 nicks Exp $";
int main (int argc, char *argv[]) {
	FILE		*imgfp, *outfp;
	IFH		ifh;
	char		imgfile[MAXL], outfile[MAXL], imgroot[MAXL], outroot[MAXL] = "";
	char		*ptr, command[MAXL], program[MAXL];
	char		control ='\0';

        int  		imgdim[4];
        float		voxdim[3];	
	float		*imgt;

	int		c, i, k;
       	int		dimension, orient, isbig;
	int		status = 0;
	int		xflag = 0, yflag = 0, zflag = 0;

	setprog (program, argv);	
	printf ("%s\n", rcsid);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) if (*argv[i] == '-') {
		strcpy (command, argv[i]); ptr = command;
		while ((c = *ptr++)) switch (c) {
			case 'x':	xflag++;	break;
			case 'y':	yflag++;	break;
			case 'z':	zflag++;	break;
			case '@': control = *ptr++;	*ptr = '\0'; break;
		}
	} else switch (k) {
		case 0:	getroot (argv[i], imgroot); k++; break;
		case 1:	getroot (argv[i], outroot); k++; break;
	}
	if (k < 1 || !(xflag || yflag || zflag)) {
		printf ("Usage:\t%s <(4dfp) image> [(4dfp) output]\n", program);
		printf ("e.g.,\t%s -yz vc345 vc345_flipyz\n", program);
		printf ("\toption\n");
		printf ("\t-x	flip x\n");
		printf ("\t-y	flip y\n");
		printf ("\t-z	flip z\n");
		printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
		printf ("N.B.:	default output fileroot = <image>_flip[xyz]\n");
		exit (1);
	}

/***************************************/
/* create output root if not specified */
/***************************************/
	if (!strlen (outroot)) {  
		sprintf (outroot, "%s_flip", imgroot);
		if (xflag) strcat (outroot, "x");
		if (yflag) strcat (outroot, "y");
		if (zflag) strcat (outroot, "z");
	}

/*****************************/
/* get input 4dfp dimensions */
/*****************************/
	sprintf (imgfile, "%s.4dfp.img", imgroot);
	if (get_4dfp_dimoe (imgfile, imgdim, voxdim, &orient, &isbig) < 0) errr (program, imgroot);
	if (Getifh (imgfile, &ifh)) errr (program, imgroot);
	if (!control) control = (isbig) ? 'b' : 'l';
	dimension = imgdim[0] * imgdim[1] * imgdim[2];
	if (!(imgt = (float *) malloc (dimension * sizeof (float)))) errm (program);

	if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
	sprintf (outfile, "%s.4dfp.img", outroot);
	if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
	fprintf (stdout, "Reading: %s\n", imgfile);
	fprintf (stdout, "Writing: %s\n", outfile);

/************/
/* process  */
/************/
	for (k = 0; k < imgdim[3]; k++) {
		if (gread  ((char *) imgt, sizeof (float), dimension, imgfp, isbig))   errr (program, imgfile);
		if (xflag) flipx (imgt, imgdim + 0, imgdim + 1, imgdim + 2);
		if (yflag) flipy (imgt, imgdim + 0, imgdim + 1, imgdim + 2);
		if (zflag) flipz (imgt, imgdim + 0, imgdim + 1, imgdim + 2);
		if (gwrite ((char *) imgt, sizeof (float), dimension, outfp, control)) errw (program, outfile);
	}
	fclose (imgfp);
	fclose (outfp);

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

	free (imgt);
	exit (status);
}
