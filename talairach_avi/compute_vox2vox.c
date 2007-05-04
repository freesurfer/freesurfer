/*$Header: /space/repo/1/dev/dev/talairach_avi/compute_vox2vox.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: compute_vox2vox.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.1  2007/05/03  00:53:52  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <endianio.h>
#include <Getifh.h>
#include <librms.h>

#define MAXL		256

/*************/
/* externals */
/*************/
extern void	f_init (void), f_exit (void);	/* FORTRAN i/o */
extern void	t4_init_ (float *t4);					/* t4_sub.f */
extern void	t4_read_ (char *t4file, float *t4);			/* t4_sub.f */
extern void	vrtflip_ (int *iori, int *imgdim, float *centeri, float *mmppixi, float *centert, float *mmppixt);	/* ft4imgo.f */

/***********/
/* globals */
/***********/
static char rcsid[] = "$Id: compute_vox2vox.c,v 1.1 2007/05/04 22:33:59 nicks Exp $";
static char program[MAXL];

void	t4list (FILE *fp, float *t4) {
	int	i;

	for (i = 0; i < 4; i++) {
		fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n", (t4 + i)[0], (t4 + i)[4], (t4 + i)[8], (t4 + i)[12]);
	}
}

void write_command_line (FILE *outfp, int argc, char *argv[]) {
	int		i;

	fprintf (outfp, "# %s", program);
	for (i = 1; i < argc; i++) fprintf (outfp, " %s", argv[i]);
	fprintf (outfp, "\n# %s\n", rcsid);
}

int main (int argc, char *argv[]) {
	FILE		*fp;
	char		srcroot[MAXL], tarroot[MAXL], t4file[MAXL], outfile[MAXL];
	IFH		ifhsrc, ifhtar;

/***********/
/* utility */
/***********/
	char		*ptr, command[MAXL];
	int		c, i, k;
	int		status = 0;
	int		four = 4;

/***************/
/* computation */
/***************/
	float		f2c[16] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,-1.0,-1.0,-1.0,1.};
	float		c2f[16] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,+1.0,+1.0,+1.0,1.};
	float		vox2ras[16], ras2vox[16], vox2voxf[16], vox2voxc[16], tempt4[16], trnsinf[16];
	float		t4[16];			/* affine warp */
	float		mmppixs[3], centers[3], mmppixt[3], centert[3];

	printf ("%s\n", rcsid);
	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; else ptr++;
	strcpy (program, ptr);
	f_init ();					/* open FORTRAN I/O */
/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); ptr = command;
			while ((c = *ptr++)) switch (c) {
			}
		} else switch (k) {
			case 2:	getroot (argv[i], tarroot);	k++; break;
			case 1:	strcpy  (t4file,  argv[i]);	k++; break;
			case 0:	getroot (argv[i], srcroot);	k++; break;
		}	
	}
	if (k < 3) {
		printf ("Usage:\t%s <(4dfp) source> <t4file> <(4dfp) target>\n", program);
		printf ("\toption\n");
		exit (1);
	}

	printf ("Reading: %s.4dfp.ifh\n", srcroot);
	if (Getifh (srcroot, &ifhsrc)) errr (program, srcroot);
/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip_ (&ifhsrc.orientation, ifhsrc.matrix_size, ifhsrc.center, ifhsrc.mmppix, centers, mmppixs);
	vrt2img_ (mmppixs, centers, ras2vox);

/***************/
/* read t4file */
/***************/
	printf ("Reading: %s\n", t4file);
	t4_read_ (t4file, t4);

	printf ("Reading: %s.4dfp.ifh\n", tarroot);
	if (Getifh (tarroot, &ifhtar)) errr (program, tarroot);
/*****************************************/
/* virtual flip instead of x4dfp2ecat () */
/*****************************************/
	vrtflip_ (&ifhtar.orientation, ifhtar.matrix_size, ifhtar.center, ifhtar.mmppix, centert, mmppixt);
	img2vrt_ (mmppixt, centert, vox2ras);

/**********************/
/* compose transforms */
/**********************/
	matmul_ (t4,       vox2ras, trnsinf,  &four);
	matmul_ (ras2vox,  trnsinf, vox2voxf, &four);
	matmul_ (vox2voxf, c2f,     tempt4,   &four);
	matmul_ (f2c,      tempt4,  vox2voxc, &four);

/****************/
/* write output */
/****************/
	if (!(ptr = strrchr (t4file, '/'))) ptr = t4file; else ptr++;
	sprintf (outfile, "%s_vox2vox.txt", t4file);
	if (!(fp = fopen (outfile, "w"))) errw (program, outfile);
	printf ("Writing: %s\n", outfile);
	fprintf (fp, "# AZS\n");
	write_command_line (fp, argc, argv);
	t4list (fp, vox2voxc);
	if (fclose (fp)) errw (program, outfile);

	f_exit ();				/* close FORTRAN I/O */
	exit (status);
} 
