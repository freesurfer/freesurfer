/**
 * @file  imgreg_4dfp.c
 * @brief compute image-image registration (t4file)
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/08/04 04:19:24 $
 *    $Revision: 1.3 $
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <endianio.h>
#include <Getifh.h>

#define MAXL	256
#define FIND	4096

/*************/
/* externals */
/*************/
extern void	f_init (void), f_exit (void);	/* FORTRAN i/o */
extern void	tparam2warp_  (int *mode, float *params, float *t4);	/* t4_sub.f */
extern void	t4file2param_ (int *mode, char *t4file, float *params);	/* t4_sub.f */
extern void	ffind_   (float *img1, short *msk1, int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1,
                      float *img2, short *msk2, int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2, float *params, int *mode);
extern void	fimgreg_ (float *img1, short *msk1, int *nx1, int *ny1, int *nz1, float *mmppix1, float *center1,
                      float *img2, short *msk2, int *nx2, int *ny2, int *nz2, float *mmppix2, float *center2, float *params, int *mode);
extern void	flipx (float *imag, int *nx, int *ny, int *nz);		/* cflip.c */
extern void	flipz (float *imag, int *nx, int *ny, int *nz);		/* cflip.c */
extern int 	x4dfp2ecat (float *imag, int *dim, int orientation);	/* below */	
	
void setprog (char *program, char **argv) {
	char *ptr;

	if (!(ptr = strrchr (argv[0], '/'))) ptr = argv[0]; 
	else ptr++;
	strcpy (program, ptr);
}

void read_file_float (char *filename, float *stack, int dimension, char *program, int isbig) {
	FILE *fp;
 
	if (!(fp = fopen (filename, "rb"))
      ||  eread (stack, dimension, isbig, fp)
      ||  fclose (fp)) errr (program, filename);
}

static char rcsid[] = "$Id: imgreg_4dfp.c,v 1.3 2007/08/04 04:19:24 nicks Exp $";
int main (int argc, char **argv) {
/************/
/* imag I/O */
/************/
	FILE		*fp;
	IFH		ifh[2], ifhm;
	char		command[MAXL];
	char		t4file[MAXL], program[MAXL];
	char		filespc[MAXL], imgroot[2][MAXL], mskroot[2][MAXL];

/**************/
/* processing */
/**************/
	float		param[13], t4[16];
	float		voxdim[3];
	float		*imag[2];
	short		*mask[2];
	int		mode;
	int		imgdim[4], mskdim[4], isbig, isbigm;
	int		nx[2], ny[2], nz[2], dimension, orientation, orientm;

/***********/
/* utility */
/***********/
	int		i, j, k;
 
	f_init ();			/* initialize FORTRAN I/O */
	fprintf (stdout, "%s\n", rcsid);
  fflush (stdout);
	setprog (program, argv);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') strcpy (command, argv[i]);
    else switch (k) {
    case 0: getroot (argv[i], imgroot[0]);		k++; break;
    case 1: getroot (argv[i], mskroot[0]); 		k++; break;
    case 2: getroot (argv[i], imgroot[1]); 		k++; break;
    case 3: getroot (argv[i], mskroot[1]); 		k++; break;
    case 4: strcpy (t4file,      argv[i]); 		k++; break;
    case 5: mode = atoi (argv[i]);			k++; break;
		}
	}
	if (k < 6) {
		printf ("Usage:\t%s target_imag target_mask source_imag source_mask t4file mode\n", program);
		printf ("or:\t%s target_imag        none source_imag source_mask t4file mode\n", program);
		printf ("or:\t%s target_imag        none source_imag        none t4file mode\n", program);
		exit (1);
	}

/*********/
/* start */
/*********/
	for (j = 0; j < 2; j++) {
    if (strlen(imgroot[j])==0) {
      fprintf (stderr, "imgreg_4dfp: imgroot[%d] is empty!\n", j);
      exit (-1);
    }
		sprintf (filespc, "%s.4dfp.img", imgroot[j]);
		if (get_4dfp_dimoe (filespc, imgdim, voxdim, &orientation, &isbig) < 0) errr (program, filespc);
		nx[j] = imgdim[0];
		ny[j] = imgdim[1];
		nz[j] = imgdim[2];
    Getifh (filespc, ifh + j);
		printf ("Reading image: %s\n", filespc);
		printf ("dimensions:%9d%10d%10d\n", nx[j], ny[j], nz[j]);
		printf ("mmppix:   %10.4f%10.4f%10.4f\n", ifh[j].mmppix[0], ifh[j].mmppix[1], ifh[j].mmppix[2]);
		printf ("center:   %10.4f%10.4f%10.4f\n", ifh[j].center[0], ifh[j].center[1], ifh[j].center[2]);

		dimension = nx[j] * ny[j] * nz[j];
		imag[j] = (float *) malloc (dimension * sizeof (float));
		mask[j]	= (short *) malloc (dimension * sizeof (short));
		if (!imag[j] || !mask[j]) errm (program);

		if (strcmp (mskroot[j], "none")) {
			sprintf (filespc, "%s.4dfp.img", mskroot[j]);
			get_4dfp_dimoe (filespc, mskdim, voxdim, &orientm, &isbigm);
      Getifh (filespc, &ifhm);
			k = orientation - ifhm.orientation;
			for (i = 0; i < 3; i++) k |= (mskdim[i] != imgdim[i]);
			for (i = 0; i < 3; i++) k |= (fabs (ifhm.mmppix[i] - ifh[j].mmppix[i]) > 0.0001);
			if (k) {
				fprintf (stderr, "imgreg_4dfp: %s %s mismatch\n", imgroot[j], mskroot[j]);
				exit (-1);
			}
			printf ("Reading mask: %s", filespc);
			read_file_float (filespc, imag[j], dimension, program, isbigm);
			x4dfp2ecat (imag[j], imgdim, orientation);
			for (k = i = 0; i < dimension; i++) 
        if ((mask[j][i] = (short) imag[j][i])) k++;
			printf (" (%d pixels)\n", k);
			for (i = 0; i < dimension; i++) mask[j][i] = (short) imag[j][i];
		} else {
			for (i = 0; i < dimension; i++) mask[j][i] = 1;
		}
		sprintf (filespc, "%s.4dfp.img", imgroot[j]);
		read_file_float (filespc, imag[j], dimension, program, isbig);
		if (x4dfp2ecat (imag[j], imgdim, orientation)) {
			fprintf (stderr, "fimgreg: %s orientation=%d not valid\n", imgroot[j], orientation);
			exit (-1);
		}
	}

/****************/
/* access check */
/****************/
	if (!access (t4file, R_OK) && access (t4file, W_OK)) {
		fprintf (stderr, "%s has read but not write permission\n", t4file);
		exit (-1);
 	}

/****************************/
/* read t4file if it exists */
/****************************/
	t4file2param_ (&mode, t4file, param);

/***********/
/* compute */
/***********/
	if (mode & FIND) {
		ffind_   (imag[0], mask[0], &nx[0], &ny[0], &nz[0], ifh[0].mmppix, ifh[0].center,
              imag[1], mask[1], &nx[1], &ny[1], &nz[1], ifh[1].mmppix, ifh[1].center, param, &mode);
 	} else {
		fimgreg_ (imag[0], mask[0], nx + 0, ny + 0, nz + 0, ifh[0].mmppix, ifh[0].center,
              imag[1], mask[1], nx + 1, ny + 1, nz + 1, ifh[1].mmppix, ifh[1].center, param, &mode); 
 	}

/****************/
/* write t4file */
/****************/
 	tparam2warp_ (&mode, param, t4);
	fp = fopen (t4file, "w");
	for (k = 0; k < argc; k++) fprintf (fp, "%s ", argv[k]);
	fprintf (fp, "\n%s\nt4\n", rcsid);
	for (k = 0; k < 4; k++) fprintf (fp, "%10.6f%10.6f%10.6f%10.4f\n",t4[0+k],t4[4+k],t4[8+k],t4[12+k]);
	if (mode & 256) fprintf (fp, "scale:    %10.6f\n", param[12]);
	fclose (fp);

	for (j = 0; j < 2; j++) {free (imag[j]); free (mask[j]);}
	f_exit ();		/* close FORTRAN I/O */
	exit (0);
}

int x4dfp2ecat (float *imag, int *dim, int orientation) {
	switch (orientation) {
  case 2:	flipx (imag, dim+0, dim+1, dim+2);	/* transverse */
    flipz (imag, dim+0, dim+1, dim+2);
    break;
  case 3:	flipx (imag, dim+0, dim+1, dim+2);	/* coronal */
    break;
  case 4: break;					/* sagittal */
  default: return -1;				/* none of the above */
	}
	return 0;
}
