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
#include <unistd.h>		/* R_OK and W_OK */
#include <time.h>
#include <pwd.h>
#include <sys/utsname.h>
#include <endianio.h>

#include "rec.h"


#define  MAXL	1024		/* accommodate very long commands */

/********************/
/* global variables */
/********************/
static char recfile[MAXL] = "";

void rec_rcsid (void) {printf ("%s\n", "freesurfer rec.c");}

const char* current_date_time() {
  time_t tt = time(&tt);
  const char* time_str = ctime(&tt);
  const char* override_time_str =
    getenv("FREESURFER_REPLACEMENT_FOR_CREATION_TIME_STRING");
  if (override_time_str) time_str = override_time_str;
  return time_str;
}

int get_machine_info (char *string) {
	struct utsname u_name;
	int rtn = 1, r;
	
	r = uname (&u_name);	
	if (r >= 0) {
		sprintf (string, "%s %s %s %s", u_name.nodename, u_name.sysname, u_name.release, u_name.machine);
 	} else {
		rtn = 0;	
 	}
 	return rtn; 
}

int printrec (char *string) {
 	FILE		*recfp;

	if (!(recfp = fopen (recfile, "a"))) {
	 	fprintf (stderr, "printrec: %s write error\n", recfile);
		return -1;
	}
	fprintf (recfp, "%s", string);
	fclose (recfp);
	return 0;
}

int catrec (char *file) {
	FILE		*recfp;
	char		filerec[MAXL], command[MAXL];
	int		k, debug = 0;
	int		isimg;

	if (access (recfile, W_OK)) {
		fprintf (stderr, "catrec: recfile not initialized\n");
		return -1;
	}
	strcpy (filerec, file);
	k = strlen (filerec);
	isimg = (!strcmp (filerec + k - 4, ".img"))
	     || (!strcmp (filerec + k - 4, ".trk"))
	     || (!strcmp (filerec + k - 5, ".conc"));
	if (isimg) strcat (filerec, ".rec");
	if (access (filerec, R_OK)) {	
		if (!(recfp = fopen (recfile, "a"))) {
	 		fprintf (stderr, "printrec: %s write error\n", recfile);
			return -1;
		}
		fprintf (recfp, "%s not found\n", filerec);
		fclose (recfp);
		return 1;
	} else {
		printf ("including: %s in %s\n", filerec, recfile);
		sprintf (command, "cat -s %s >> %s", filerec, recfile);
		if (debug) fprintf (stderr, "%s\n", command);
		return system (command);
	}
}

int startrec (char *outfile, int argc, char *argv[], char *rcsid) {
	extern void	get_time_usr (char *string);
	FILE		*recfp;
	char 		*str, string[MAXL];
	int 		k;

	strcpy (string, outfile);
	while ((str = strrchr (string, '.'))) {
		if (!strcmp (str, ".rec")) *str = '\0';
		else break;
	}
	strcpy (recfile, string);
	strcat (recfile, ".rec");
	if (!(recfp = fopen (recfile, "w"))) {
		fprintf (stderr, "startrec: %s write error\n", recfile);
		return -1;
	}
	fprintf (recfp, "rec %s", string);
	get_time_usr (string);
	fprintf (recfp, "  %s\n", string);
	for (k = 0; k < argc; k++) fprintf (recfp, "%s ", argv[k]);
	fprintf (recfp, "\n%s\n", rcsid); 
	fclose (recfp);
	return 0;
}

int startrecl (char *outfile, int argc, char *argv[], char *rcsid) {
	extern void	get_time_usr (char *string);
	FILE		*recfp;
	char 		*str, string[MAXL];
	int 		k;

	strcpy (string, outfile);
	while ((str = strrchr (string, '.'))) {
		if (!strcmp (str, ".rec")) *str = '\0';
		else break;
	}
	strcpy (recfile, string);
	strcat (recfile, ".rec");
	if (!(recfp = fopen (recfile, "w"))) {
		fprintf (stderr, "startrecl: %s write error\n", recfile);
		return -1;
	}
	fprintf (recfp, "rec %s", string);
	get_time_usr (string);
	fprintf (recfp, "  %s\n", string);
	for (k = 0; k < argc - 1; k++) fprintf (recfp, "\t%s\t\\\n", argv[k]);
	fprintf (recfp, "\t%s\n", argv[k]);
	fprintf (recfp, "%s\n", rcsid); 
	fclose (recfp);
	return 0;
}

int startrece (char *outfile, int argc, char *argv[], char *rcsid, char control) {
	extern void	get_time_usr (char *string);
	FILE		*recfp;
	char 		*str, string[MAXL];
	int 		k, osbig;

	strcpy (string, outfile);
	while ((str = strrchr (string, '.'))) {
		if (!strcmp (str, ".rec")) *str = '\0';
		else break;
	}
	strcpy (recfile, string);
	strcat (recfile, ".rec");
	if (!(recfp = fopen (recfile, "w"))) {
		fprintf (stderr, "startrece: %s write error\n", recfile);
		return -1;
	}
	fprintf (recfp, "rec %s", string);
	get_time_usr (string);
	fprintf (recfp, "  %s", string);
	if (get_machine_info (string)) {
	 	fprintf (recfp, "@%s\n", string);
	} else {
		fprintf (recfp, "\n");
	}	
	for (k = 0; k < argc; k++) fprintf (recfp, "%s ", argv[k]);
	fprintf (recfp, "\n%s\n", rcsid);
	switch (control) {
		case 'b': case 'B': osbig = 1; break;
		case 'l': case 'L': osbig = 0; break;
		default: osbig = CPU_is_bigendian(); break;
	}
	fprintf (recfp, "%s\n", ((osbig) ? "bigendian" : "littleendian")); 
	fclose (recfp);
	return 0;
}

int startrecle (char *outfile, int argc, char *argv[], char *rcsid, char control) {
	extern void	get_time_usr (char *string);
	FILE		*recfp;
	char 		*str, string[MAXL];
	int 		k, osbig;

	strcpy (string, outfile);
	while ((str = strrchr (string, '.'))) {
		if (!strcmp (str, ".rec")) *str = '\0';
		else break;
	}
	strcpy (recfile, string);
	strcat (recfile, ".rec");
	if (!(recfp = fopen (recfile, "w"))) {
		fprintf (stderr, "startrecle: %s write error\n", recfile);
		return -1;
	}
	fprintf (recfp, "rec %s", string);
	get_time_usr (string);
	fprintf (recfp, "  %s", string);
	if (get_machine_info (string)) {
	 	fprintf (recfp, "@%s\n", string);
	} else {
		fprintf (recfp,"\n");
	}	
	for (k = 0; k < argc - 1; k++) fprintf (recfp, "\t%s\t\\\n", argv[k]);
	fprintf (recfp, "\t%s\n", argv[k]);
	fprintf (recfp, "%s\n", rcsid); 
	switch (control) {
		case 'b': case 'B': osbig = 1; break;
		case 'l': case 'L': osbig = 0; break;
		default: osbig = CPU_is_bigendian(); break;
	}
	fprintf (recfp, "%s\n", ((osbig) ? "bigendian" : "littleendian")); 
	fclose (recfp);
	return 0;
}

int endrec (void) {
	extern void	get_time_usr (char *string);
	FILE		*recfp;
	char 		string[MAXL];
	
	if (!(recfp = fopen (recfile, "a"))) {
		fprintf (stderr, "endrec: recfile write error\n");
		return -1;
	}
	get_time_usr (string);
	fprintf (recfp, "endrec %s\n", string);
	fclose (recfp);
	return 0;
}

void get_time_usr (char *string) {

	struct passwd	*pw;

	strcpy (string, current_date_time());
	string[24] = '\0';
	strcat (string, "  ");
	pw = getpwuid(geteuid());
	if ((pw != NULL) && (pw->pw_name != NULL)) {
		strcat (string, pw->pw_name);
	} else {
		strcat (string, "UNKNOWN");
	}
}
