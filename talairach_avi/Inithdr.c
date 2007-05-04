/*$Header: /space/repo/1/dev/dev/talairach_avi/Inithdr.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: Inithdr.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.7  2007/05/04  03:08:38  avi
 * linux gcc compliant
 * cuserid -> getpwuid
 *
 * Revision 1.6  2005/12/16  02:41:13  avi
 * #include ANALYZE.h -> generic; point to correct file in ifh2hdr.mak
 *
 * Revision 1.5  2003/09/25  04:01:44  avi
 * get_date_log () -> in-line call to time() and cuserid() => eliminate need for libmri
 *
 * Revision 1.4  1999/06/24  06:33:01  avi
 * allow proto_imgfile to be (char *NULL)
 *
 * Revision 1.3  1998/12/30  04:53:05  avi
 * eliminate message in absence of proto_header
 *
 * Revision 1.1  1998/05/21  16:26:56  tscull
 * Initial revision
 **/
/*_________________________________________________________________________
  File:		Inithdr

  Usage:	Inithdr (&hdr, imgdim, voxdim, proto_imgfile)

  Description:	Initialize ANALYZE header using prototype ANALYZE image file

  Author:	AZS

  History:	05/09/98
___________________________________________________________________________*/

#include <string.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <ANALYZE.h>
#include <unistd.h>
#include <pwd.h>

#define MAXL		256

static char rcsid[] = "$Id: Inithdr.c,v 1.1 2007/05/04 22:33:59 nicks Exp $";
void Inithdr_rcs (void) {printf ("%s\n", rcsid);}
int Inithdr (struct dsr *phdr, int *imgdim, float *voxdim, char *proto_imgfile) {
	FILE		*fp;
	char		*str, string[MAXL], proto_hdr[MAXL];
	int		j, k, status;
	int		debug = 1;
	time_t		time_sec;
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

	time (&time_sec);
	strcpy (string, ctime (&time_sec));
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
