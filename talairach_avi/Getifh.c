/**
 * @file  Getifh.c
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/08/04 02:23:35 $
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <ifh.h>
#include <endianio.h>

#define MAXL	256

static char	rcsid[] = "$Id: Getifh.c,v 1.3 2007/08/04 02:23:35 nicks Exp $";
void Getifh_rcs (void) {printf ("%s\n", rcsid);}

int Getifh (char *imgfile, IFH *ifhdr) {
	FILE	*fp;
	char	*str, ifhfile[MAXL], line[MAXL], parameter[MAXL];
	int 	i;
	int	mmppix_flag = 0, center_flag = 0, endian_flag = 0;
	int	debug = 0;

	getroot (imgfile, ifhfile);
	strcat (ifhfile, ".4dfp.ifh");
	if (!(fp = fopen (ifhfile, "r"))) {
		fprintf (stderr, "Getifh: %s read error, imgfile='%s'\n", 
             ifhfile, imgfile);
		return -1;
	}
	if (0) printf ("Reading: %s\n", ifhfile);	/* stdout messages can be a problem */ 
	memset (ifhdr, '\0', sizeof (IFH));
	memset (line, '\0', MAXL);
	while (fgets (line, MAXL, fp)) {
		if ((str = strrchr (line, '\n'))) *str ='\0';
		if ((str = strchr  (line, '#')))  *str ='\0';
		if (!(str = strstr (line, ":="))) continue;
		*str = '\0';		/* terminate keyword part of line */
		str += 2;
		while (isspace (*str)) str++;
		strcpy (parameter, str);

		for (i = 0; i < strlen (line); i++) line[i] = tolower (line[i]);
		if (debug) printf ("%s\t%s\n", line, parameter);

		if (strstr (line, "version of keys")) {
			strncpy (ifhdr->version_of_keys,	parameter, 31);
		}
		if (strstr (line, "conversion program")) {
			strncpy (ifhdr->conversion_program,	parameter, 255);
		}
		if (strstr (line, "name of data file")) {
			strncpy (ifhdr->name_of_data_file,	parameter, 255);
		}
		if (strstr (line, "number format")) {
			strncpy (ifhdr->number_format,		parameter, 31);
		}
		if (strstr (line, "imagedata byte order")) {
			strncpy (ifhdr->imagedata_byte_order,	parameter, 31);
			endian_flag++;
		}
		if (strstr (line, "number of bytes per pixel")) {
			ifhdr->number_of_bytes_per_pixel = atoi (parameter);
		}
		if (strstr (line, "number of dimensions")) {
			ifhdr->number_of_dimensions = atoi (parameter);
		}
		if (strstr (line, "matrix size")) {
			if (!(str = strchr (line, '['))) goto ERR;
			sscanf (str, "[%d]", &i);
			if (i < 1 || i > 4) goto ERR;		
			ifhdr->matrix_size[i-1] = atoi (parameter);
		}
		if (strstr (line, "orientation")) {
			ifhdr->orientation = atoi (parameter);
		}
		if (strstr (line, "scaling factor (mm/pixel)")) {
			if (!(str = strchr (line, '['))) goto ERR;
			sscanf (str, "[%d]", &i);
			if (i < 1 || i > 4) goto ERR;		
			ifhdr->scaling_factor[i-1] = atof (parameter);
		}
		if (strstr (line, "center")) {
			sscanf (parameter, "%f %f %f", ifhdr->center+0, ifhdr->center+1, ifhdr->center+2);
			center_flag++;
		}
		if (strstr (line, "mmppix")) {
			sscanf (parameter, "%f %f %f", ifhdr->mmppix+0, ifhdr->mmppix+1, ifhdr->mmppix+2);
			mmppix_flag++;
		}
	}
	fclose (fp);

	if (!endian_flag) {
		strcpy (ifhdr->imagedata_byte_order, "bigendian");
	}
	if (!mmppix_flag) {
		ifhdr->mmppix[0] =  ifhdr->scaling_factor[0];
		ifhdr->mmppix[1] = -ifhdr->scaling_factor[1];
		ifhdr->mmppix[2] = -ifhdr->scaling_factor[2];
	}
	if (!center_flag) {
		ifhdr->center[0] = ifhdr->mmppix[0] * (float) (ifhdr->matrix_size[0] - ifhdr->matrix_size[0]/2);
		ifhdr->center[1] = ifhdr->mmppix[1] * (float) (1 + ifhdr->matrix_size[1]/2);
		ifhdr->center[2] = ifhdr->mmppix[2] * (float) (1 + ifhdr->matrix_size[2]/2);
	}

	return 0;
ERR:	if ((str = strrchr (parameter, '\n'))) *str = '\0'; 
	fprintf (stderr, ">>> %s := %s <<<\n", line, parameter);
	return -1;
}

int Writeifh (char *program, char *outfile, IFH *ifhdr, char control) {
	FILE		*ifhfp;
	char		ifhfile[MAXL];
	int		i, osbig;

	osbig = (CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B');

	getroot (outfile, ifhfile);
	strcat (ifhfile, ".4dfp.ifh");
	if (!(ifhfp = fopen (ifhfile, "w"))) return -1;
	printf ("Writing: %s\n", ifhfile);

	fprintf (ifhfp, "INTERFILE	:=\n");
	fprintf (ifhfp, "version of keys	:= %s\n", ifhdr->version_of_keys);
	fprintf (ifhfp, "number format		:= %s\n", ifhdr->number_format);
	fprintf (ifhfp, "conversion program	:= %s\n", program);
	fprintf (ifhfp, "name of data file	:= %s\n", outfile);
	fprintf (ifhfp, "number of bytes per pixel	:= %d\n", ifhdr->number_of_bytes_per_pixel);
	fprintf (ifhfp, "imagedata byte order	:= %s\n", (osbig) ? "bigendian" : "littleendian");
	fprintf (ifhfp, "orientation		:= %d\n", ifhdr->orientation);
	fprintf (ifhfp, "number of dimensions	:= %d\n", ifhdr->number_of_dimensions);
for (i = 0; i < 4; i++) {
	fprintf (ifhfp, "matrix size [%d]	:= %d\n", i + 1, ifhdr->matrix_size[i]);
}
for (i = 0; i < 3; i++) {
	fprintf (ifhfp, "scaling factor (mm/pixel) [%d]	:= %f\n", i + 1, ifhdr->scaling_factor[i]);
}
	fprintf (ifhfp, "mmppix	:= %10.6f%10.6f%10.6f\n", ifhdr->mmppix[0], ifhdr->mmppix[1], ifhdr->mmppix[2]);
	fprintf (ifhfp, "center	:= %10.4f%10.4f%10.4f\n", ifhdr->center[0], ifhdr->center[1], ifhdr->center[2]);

	if (fclose (ifhfp)) return -1;
	return 0;
}

int writeifhe (char *program, char *outfile, int *imgdim, float *voxdim, int orient, char control) {
	FILE		*ifhfp;
	char		ifhfile[MAXL];
	int		osbig;

	osbig = (CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B');
	getroot (outfile, ifhfile);
	strcat (ifhfile, ".4dfp.ifh");
	if (!(ifhfp = fopen (ifhfile, "w"))) return -1;
	printf ("Writing: %s\n", ifhfile);

	fprintf (ifhfp, "INTERFILE	:=\n");
	fprintf (ifhfp, "version of keys	:= 3.3\n");
	fprintf (ifhfp, "number format		:= float\n");
	fprintf (ifhfp, "conversion program	:= %s\n", program);
	fprintf (ifhfp, "name of data file	:= %s\n", outfile);
	fprintf (ifhfp, "number of bytes per pixel	:= %d\n", 4);
	fprintf (ifhfp, "imagedata byte order	:= %s\n", (osbig) ? "bigendian" : "littleendian");
	fprintf (ifhfp, "orientation		:= %d\n", orient);
	fprintf (ifhfp, "number of dimensions	:= %d\n", 4);
	fprintf (ifhfp, "matrix size [1]	:= %d\n", imgdim[0]);
	fprintf (ifhfp, "matrix size [2]	:= %d\n", imgdim[1]);
	fprintf (ifhfp, "matrix size [3]	:= %d\n", imgdim[2]);
	fprintf (ifhfp, "matrix size [4]	:= %d\n", imgdim[3]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [1]	:= %f\n", voxdim[0]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [2]	:= %f\n", voxdim[1]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [3]	:= %f\n", voxdim[2]);

	if (fclose (ifhfp)) return -1;
	return 0;
}

int writeifhmc (char *program, char *outfile, int *imgdim, float *voxdim, int orient, float *mmppix, float *center) {
	FILE		*ifhfp;
	char		*ptr, imgroot[MAXL];

	getroot (outfile, imgroot);
	strcat (imgroot, ".4dfp.ifh");
	if (!(ifhfp = fopen (imgroot, "w"))) errw (program, imgroot);
	fprintf (ifhfp, "INTERFILE	:=\n");
	fprintf (ifhfp, "version of keys	:= 3.3\n");
	fprintf (ifhfp, "number format		:= float\n");
	fprintf (ifhfp, "conversion program	:= %s\n", program);
	if (!(ptr = strrchr (outfile, '/'))) ptr = outfile; else ptr++;
	fprintf (ifhfp, "name of data file	:= %s\n", ptr);
	fprintf (ifhfp, "number of bytes per pixel	:= %d\n", 4);
	fprintf (ifhfp, "imagedata byte order	:= %s\n", (CPU_is_bigendian ()) ? "bigendian" : "littleendian");
	fprintf (ifhfp, "orientation		:= %d\n", orient);
	fprintf (ifhfp, "number of dimensions   := %d\n", 4);
	fprintf (ifhfp, "matrix size [1]	:= %d\n", imgdim[0]);
	fprintf (ifhfp, "matrix size [2]	:= %d\n", imgdim[1]);
	fprintf (ifhfp, "matrix size [3]	:= %d\n", imgdim[2]);
	fprintf (ifhfp, "matrix size [4]	:= %d\n", imgdim[3]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [1]	:= %f\n", voxdim[0]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [2]	:= %f\n", voxdim[1]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [3]	:= %f\n", voxdim[2]);
	fprintf (ifhfp, "mmppix	:= %10.6f%10.6f%10.6f\n", mmppix[0], mmppix[1], mmppix[2]);
	fprintf (ifhfp, "center	:= %10.4f%10.4f%10.4f\n", center[0], center[1], center[2]);
	if (fclose (ifhfp)) return -1;
	return 0;
}

int writeifhmce (char *program, char *outfile, int *imgdim, float *voxdim, int orient,
		float *mmppix, float *center, char control) {
	FILE		*ifhfp;
	char		ifhfile[MAXL];
	int		osbig;

	osbig = (CPU_is_bigendian ()) ? !(control == 'l' || control == 'L') : (control == 'b' || control == 'B');

	getroot (outfile, ifhfile);
	strcat (ifhfile, ".4dfp.ifh");
	if (!(ifhfp = fopen (ifhfile, "w"))) return -1;
	printf ("Writing: %s\n", ifhfile);

	fprintf (ifhfp, "INTERFILE	:=\n");
	fprintf (ifhfp, "version of keys	:= 3.3\n");
	fprintf (ifhfp, "number format		:= float\n");
	fprintf (ifhfp, "conversion program	:= %s\n", program);
	fprintf (ifhfp, "name of data file	:= %s\n", outfile);
	fprintf (ifhfp, "number of bytes per pixel	:= %d\n", 4);
	fprintf (ifhfp, "imagedata byte order	:= %s\n", (osbig) ? "bigendian" : "littleendian");
	fprintf (ifhfp, "orientation		:= %d\n", orient);
	fprintf (ifhfp, "number of dimensions	:= %d\n", 4);
	fprintf (ifhfp, "matrix size [1]	:= %d\n", imgdim[0]);
	fprintf (ifhfp, "matrix size [2]	:= %d\n", imgdim[1]);
	fprintf (ifhfp, "matrix size [3]	:= %d\n", imgdim[2]);
	fprintf (ifhfp, "matrix size [4]	:= %d\n", imgdim[3]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [1]	:= %f\n", voxdim[0]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [2]	:= %f\n", voxdim[1]);
	fprintf (ifhfp, "scaling factor (mm/pixel) [3]	:= %f\n", voxdim[2]);
	fprintf (ifhfp, "mmppix	:= %10.6f%10.6f%10.6f\n", mmppix[0], mmppix[1], mmppix[2]);
	fprintf (ifhfp, "center	:= %10.4f%10.4f%10.4f\n", center[0], center[1], center[2]);

	if (fclose (ifhfp)) return -1;
	return 0;
}
