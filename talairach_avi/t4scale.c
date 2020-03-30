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
#include <math.h>
#include <string.h>
#include <stdio.h>

#define MAXL		256

float t4scale (char *t4file) {
	FILE		*fp;
	float		q;
	char		*ptr, string[MAXL];
	
	q = 1.0;
	if (!(fp = fopen (t4file, "r"))) {
		fprintf (stderr, "t4scale: %s read error\n", t4file);
		exit (-1);
	}
	while (fgets (string, MAXL, fp)) if ((ptr = strstr (string, "scale:"))) q = atof (ptr + 6);
	fclose (fp);
	return q;
}
