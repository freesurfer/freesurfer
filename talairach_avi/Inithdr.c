/**
 * @brief Initialize ANALYZE header using prototype ANALYZE image file
 *
 */
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

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <ANALYZE.h>
#include <unistd.h>
#include <pwd.h>

#include "rec.h"

#define MAXL		256

void Inithdr_rcs (void) {printf ("%s\n", "freesurfer Inithdr.c");}
int Inithdr (struct dsr *phdr, int *imgdim, float *voxdim, char *proto_imgfile) {
	FILE		*fp;
	char		*str, string[MAXL], proto_hdr[MAXL];
	int		j, k, status;
	int		debug = 1;
	struct passwd	*pw;

	memset (phdr, '\0', sizeof (struct dsr));
	phdr->hk.sizeof_hdr = sizeof (struct dsr); 	/* required */
	phdr->hk.extents = 16384;			/* recommended */
	phdr->hk.regular = 'r';				/* required */
	phdr->dime.datatype = 4;			/* short int */
	phdr->dime.bitpix = 16;

	status = 0;
	if (strlen (proto_imgfile)) {
		strcpy (proto_hdr, proto_imgfile);
		while ((str = strrchr (proto_hdr, '.'))) {
			if (!strcmp (str, ".rec")) *str = '\0';
			if (!strcmp (str, ".img")) *str = '\0';
			if (!strcmp (str, ".hdr")) *str = '\0';
			else break;
		}
		strcat (proto_hdr, ".hdr");
		if (!(fp = fopen (proto_hdr, "rb")) || fread (phdr, sizeof (struct dsr), 1, fp) != 1) {
			if (debug) fprintf (stderr, "Inithdr: cannot read %s\n", proto_hdr);
			status = 1;
		} else {
			fclose (fp);
		}
		str = strrchr (proto_imgfile, '/');
		if (str) str++; else str = proto_imgfile;
		strncpy (phdr->hist.descrip, str, 79);
	}

	phdr->dime.dim[0] = 4;				/* 4 dimensions */
	phdr->dime.dim[1] = imgdim[0];
	phdr->dime.dim[2] = imgdim[1];
	phdr->dime.dim[3] = imgdim[2];
	phdr->dime.dim[4] = imgdim[3];
	phdr->dime.pixdim[1] = voxdim[0];
	phdr->dime.pixdim[2] = voxdim[1];
	phdr->dime.pixdim[3] = voxdim[2];

	strcpy (string, current_date_time());
	string[24] = '\0';
	if (debug) printf ("%s\n", string);
	for (j = k = 0; k < 10; k++) if (string[k] != ' ') phdr->hist.exp_date[j++] = string[k];
	strncpy (phdr->hist.exp_time,  string + 11, 9);
	strncpy (phdr->hist.generated, string + 20, 9);

	pw = getpwuid(geteuid());
	if ((pw != NULL) && (pw->pw_name != NULL)) {
		strcpy (string, pw->pw_name);
	} else {
		strcpy (string, "UNKNOWN");
	}
	strncpy (phdr->hist.originator, string, 9);
	return status;
}
