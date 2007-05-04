/*$Header: /space/repo/1/dev/dev/talairach_avi/gauss_4dfp.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: gauss_4dfp.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.15  2007/04/23  02:37:19  avi
 * gcc v3 compliant (filter subroutines converted to C)
 * remove f_init() and f_exit()
 *
 * Revision 1.14  2006/09/25  18:59:24  avi
 * correct bug computing conc outfile
 *
 * Revision 1.13  2006/09/25  18:34:31  avi
 * Solaris 10
 *
 * Revision 1.12  2005/12/06  06:37:37  avi
 * conc file capability
 *
 * Revision 1.11  2005/12/02  06:56:48  avi
 * better usage
 *
 * Revision 1.10  2005/07/02  03:11:57  avi
 * report current volume to stdout
 *
 * Revision 1.9  2004/10/08  18:08:45  rsachs
 * Installed 'errm','errr','errw','getroot','setprog'. Replaced 'Get4dfpDi'
 *
 * Revision 1.8  1999/01/18  03:23:54  avi
 * eliminate #include <mri/mri.h>
 * initialize outroot to ""
 *
 * Revision 1.6  1998/12/03  00:28:02  avi
 * -d option
 *
 * Revision 1.4  1998/05/20  07:29:03  avi
 * new rec subroutines
 *
 * Revision 1.3  1998/05/20  07:18:43  avi
 * clean code
 * -w option
 *
 * Revision 1.2  1998/04/17  17:12:41  tscull
 * set debug to false
 *
 * Revision 1.1  1998/03/18  18:25:56  tscull
 * Initial revision
 *
 * Revision 1.3  1997/10/02  18:14:48  tscull
 * modified usage text to be more complete and easy to read
 *
 * Revision 1.2  1997/10/02  17:56:38  tscull
 * frequency to text more compact
 *
 * Revision 1.1  1997/09/30  21:00:19  tscull
 * Initial revision
 **/
/****************************************************************
  Description:	This program filters a 4dfp image volume
		by using the Gaussian filter.

  History:	Created by Tom Yang and Avi Snyder on 12/17/92. 
                Originally for ECAT images
		Modified by AZS on 10/25/95.
		Converted to 4dfp input by Tom Cull 9/30/97.
*****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <math.h>
#include <endianio.h>
#include <Getifh.h>
#include <rec.h>
#include <conc.h>

#define MAXL      256	

/*************/
/* externals */
/*************/
extern int	npad_ (int *n, int *margin);									/* FORTRAN librms */
extern void	imgpad_   (float *imag, int *nx, int *ny, int *nz, float *imgp, int *nxp, int *nyp, int *nzp);	/* FORTRAN librms */
extern void	imgdap_   (float *imag, int *nx, int *ny, int *nz, float *imgp, int *nxp, int *nyp, int *nzp);	/* FORTRAN librms */
extern void	gauss3d   (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);		/* cgauss3d.c */
extern void	gauss3dd  (float *imag, int *nx, int *ny, int *nz, float *cmppix, float *fhalf);		/* cgauss3dd.c */

void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void usage (char* program) {
	printf ("Usage:\t%s <4dfp|conc input> f_half [outroot]\n", program);
	printf (" e.g.,\t%s pt349_study9to9 0.1\n", program);
	printf (" e.g.,\t%s p1234ho5 0.7 p1234ho5_g7\n", program);
	printf ("\toptions\n");
	printf ("\t-@<b|l>\toutput big or little endian (default input endian)\n");
	printf ("\t-w\t(wrap) suppress x and y padding\n");
	printf ("\t-d\tdifferentiate\n");
	printf ("N.B.:	f_half is half frequency in 1/cm\n");
	printf ("N.B.:	default output root is <inroot>_g<10*f_half>\n");
	printf ("N.B.:	FWHM*f_half = (2ln2/pi) = 0.4412712\n");
	printf ("N.B.:	conc files must have extension \"conc\"\n");
	printf ("N.B.:	user outroot specification not possible with conc files\n");
	exit (1);
}

static char rcsid[] = "$Id: gauss_4dfp.c,v 1.1 2007/05/04 22:33:59 nicks Exp $";
int main (int argc, char **argv) {
	CONC_BLOCK	conc_block;			/* conc i/o control block */
	FILE 		*imgfp=NULL, *outfp=NULL;
	IFH		ifh;
	char 		imgroot[MAXL], imgfile[MAXL];
	char 		outroot[MAXL] = "", outfile[MAXL], trailer[MAXL];

        int  		imgdim[4], isbig;
        float	 	voxdim[3];	
	float		cmppix[3], f0;
	float		*imgt, *imgp;
	int		nx, ny, nz;
	int		nxp, nyp, nzp;
	int		margin, vdim;
	char		control = '\0';

/***********/
/* utility */
/***********/
	char 		command[MAXL], program[MAXL], *ptr;
	float		val;
	int		c, i, k;

/*********/
/* flags */
/*********/
	int		conc_flag = 0;
	int		status = 0;
	int		wrap_flag = 0;
	int		diff_flag = 0;

	fprintf (stdout, "%s\n", rcsid);
	setprog (program, argv);
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]);
			ptr = command;
			while ((c = *ptr++)) switch (c) {
				case 'w': wrap_flag++;		break;
				case 'd': diff_flag++;		break;
				case '@': control = *ptr++;	*ptr = '\0'; break;
			}
		}
		else switch (k) {
			case 0:	getroot (argv[i], imgroot);
				conc_flag = (strstr (argv[i], ".conc") == argv[i] + strlen (imgroot));
								k++; break;
			case 1:	f0 = atof (argv[i]);		k++; break;
			case 2: if (conc_flag) usage (program);
				getroot (argv[i], outroot);	k++; break;
		}	
	}
	if (k < 2) usage (program);

/***************************/
/* compute outroot trailer */
/***************************/
	sprintf (trailer, "%sg%d", (diff_flag) ? "d" : "", (int) (10.0*f0 + 0.5));

/*******************************************/
/* get 4dfp dimensions and open read/write */
/*******************************************/
	if (conc_flag) {
		conc_init (&conc_block, program);
		conc_open (&conc_block, imgroot);
		strcpy (imgfile, conc_block.lstfile);
		for (k = 0; k < 4; k++) imgdim[k] = conc_block.imgdim[k];
		for (k = 0; k < 3; k++) voxdim[k] = conc_block.voxdim[k];
		isbig = conc_block.isbig;
		strcpy (outfile, conc_block.outfile);
	} else {
		sprintf (imgfile, "%s.4dfp.img", imgroot);
		if (Getifh (imgfile, &ifh)) errr (program, imgfile);
		for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];
		for (k = 0; k < 3; k++) voxdim[k] = ifh.scaling_factor[k];
		isbig = strcmp (ifh.imagedata_byte_order, "littleendian");
		if (!(imgfp = fopen (imgfile, "rb"))) errr (program, imgfile);
		if (!strlen (outroot)) sprintf (outroot, "%s_%s", imgroot, trailer);
		sprintf (outfile, "%s.4dfp.img", outroot);
		if (!(outfp = fopen (outfile, "wb"))) errw (program, outfile);
	}
	if (!control) control = (isbig) ? 'b' : 'l';
	if (conc_flag) {conc_newe (&conc_block, trailer, control); strcpy (outfile, conc_block.outfile);}
	printf ("Reading: %s\n", imgfile);
	printf ("Writing: %s\n", outfile);
	nx = imgdim[0];
	ny = imgdim[1];
	nz = imgdim[2];
	vdim = nx * ny * nz;
	for (k = 0; k < 3; k++) cmppix[k] = voxdim[k] / 10.0;

/********************/
/* allocate buffers */
/********************/
	if (wrap_flag) {
		nxp = nx;
		nyp = ny;
	} else {
		val = (0.5 + (2.0 * 0.1874 / (cmppix[0] * f0))); margin = val; nxp = npad_ (&nx, &margin);
		val = (0.5 + (2.0 * 0.1874 / (cmppix[1] * f0))); margin = val; nyp = npad_ (&ny, &margin);
	}
	val = (0.5 + (4.0 * 0.1874 / (cmppix[2] * f0))); margin = val; nzp = npad_ (&nz, &margin);
	printf ("image dimensions %d %d %d padded to %d %d %d\n", nx, ny, nz, nxp, nyp, nzp);
	imgt = (float *) malloc (vdim * sizeof (float));
	imgp = (float *) calloc (nxp * nyp * nzp, sizeof (float));
	if (!imgp || !imgt) errm (program);

	printf ("processing volume");
	for (k = 0; k < imgdim[3]; k++) {printf(" %d", k + 1); fflush (stdout);
		if (conc_flag) {
			conc_read_vol (&conc_block, imgt);
		} else {
			if (eread (imgt, vdim, isbig, imgfp)) errr (program, imgfile);
		}
		imgpad_ (imgt, &nx, &ny, &nz, imgp, &nxp, &nyp, &nzp);
		if (diff_flag) {
			gauss3dd (imgp, &nxp, &nyp, &nzp, cmppix, &f0);
		} else {
			gauss3d  (imgp, &nxp, &nyp, &nzp, cmppix, &f0);
		}
		imgdap_ (imgt, &nx, &ny, &nz, imgp, &nxp, &nyp, &nzp);
		if (conc_flag) {
			conc_write_vol (&conc_block, imgt);
		} else {
			if (ewrite (imgt, vdim, control, outfp)) errw (program, outfile);
		}
	}
	printf("\n");

/***************/
/* ifh hdr rec */
/***************/
	if (conc_flag) {
		status |= conc_ifh_hdr_rec (&conc_block, argc, argv, rcsid);
		conc_free (&conc_block);
	} else {
		if (fclose (imgfp)) errr (program, imgfile);
		if (fclose (outfp)) errw (program, outfile);
		if (Writeifh (program, outfile, &ifh, control)) errw (program, outroot);
		sprintf (command, "ifh2hdr %s", outroot);
		status |= system (command);
	}
	startrece (outfile, argc, argv, rcsid, control);
	catrec (imgfile);
	endrec ();

	free (imgp);
	free (imgt);
	exit (status);
}
