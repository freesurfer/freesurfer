/*$Header: /space/repo/1/dev/dev/talairach_avi/ifh2hdr.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: ifh2hdr.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.9  2007/02/28  05:35:19  avi
 * Solaris 10
 *
 * Revision 1.8  2006/03/26  00:10:17  avi
 * eliminate redundant screen message
 *
 * Revision 1.7  2006/03/25  04:30:55  avi
 * use endianio.c subroutines
 *
 * Revision 1.6  2006/03/24  06:10:16  avi
 * preserve endian state of input
 *
 * Revision 1.5  2005/12/16  02:38:09  avi
 * #include ifh.h and ANALYZE.h -> generic; point to correct file in ifh2hdr.mak
 *
 * Revision 1.4  2003/09/25  03:06:25  avi
 * eliminate call to Get4dfpDimN
 *
 * Revision 1.3  2003/03/15  04:28:43  avi
 * modernize error routines
 *
 * Revision 1.2  1999/11/23  02:48:03  avi
 * include ifh.orientation in output header
 *
 * Revision 1.1  1999/08/28  23:18:15  avi
 * Initial revision
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Getifh.h>
#include <endianio.h>

#define MAXL	256

void getrange (char *string, float *minval, float *maxval) {
        char	*str;

	str = strstr (string, "to");
	if (str) {
		*str = '\0';
		*minval = atof (string);
		*maxval = atof (str + 2);
	} else {
		*minval = 0.0;
		*maxval = atof (string);
	}
}

extern int	Inithdr (struct dsr *phdr, int *imgdim, float *voxdim, char *proto_imgfile);

static char rcsid[] = "$Id: ifh2hdr.c,v 1.1 2007/05/04 22:33:59 nicks Exp $";
int main (int argc, char *argv[]) {
	FILE		*fp;
	struct dsr	hdr;
	IFH		ifh;

	char		filespc[MAXL], imgroot[MAXL];
	float		voxsiz[3];
	float           fmin, fmax;
	int		imgdim[4];

/***********/
/* utility */
/***********/
	char		*str, command[MAXL], program[MAXL];
	int 		c, i, k;

/*********/
/* flags */
/*********/
	int		verbose = 0;

	int		isbig, swab_flag, range_flag = 0;

	printf ("%s\n", rcsid);
	if (!(str = strrchr (argv[0], '/'))) str = argv[0]; else str++;
	strcpy (program, str);

/************************/
/* process command line */
/************************/
	for (k = 0, i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			strcpy (command, argv[i]); str = command;
			while ((c = *str++)) switch (c) {
				case 'r':	getrange (str, &fmin, &fmax);
						*str = '\0'; range_flag++;		break;
				case 'v':	verbose++;				break;
			}
		} else switch (k) {
			case 0:	getroot (argv[i], imgroot);	k++; break;
		}	
	}
	if (k < 1) {
		printf ("Usage:\t%s <(4dfp) file>\n", program);
		printf ("e.g.,\t%s vc654_mpr_atl -r-500to1500\n", program);
		printf ("\toption\n");
		printf ("\t-r<flt>[to<flt>]\tset range\n");
		printf ("N.B.:\t%s preserves the endian state of the input ifh\n", program);
		exit (1);
	}

	sprintf (filespc, "%s.4dfp.ifh", imgroot);
	if (Getifh (filespc, &ifh)) errr (program, filespc);
	for (k = 0; k < 3; k++) voxsiz[k] = ifh.scaling_factor[k];
	for (k = 0; k < 4; k++) imgdim[k] = ifh.matrix_size[k];

	Inithdr (&hdr, imgdim, voxsiz, "");
	hdr.dime.datatype = 16;			/* float */
	hdr.dime.bitpix = 32;
	hdr.hist.orient = ifh.orientation - 2;
	if (range_flag) {
		hdr.dime.glmin = fmin;
		hdr.dime.glmax = fmax;
	}

	isbig = !strcmp (ifh.imagedata_byte_order, "bigendian");
	swab_flag = (CPU_is_bigendian() != 0) != (isbig != 0);
	if (swab_flag) swab_hdr (&hdr);

	sprintf (filespc, "%s.4dfp.hdr", imgroot);
	printf ("Writing: %s\n", filespc);
	if (!(fp = fopen (filespc, "wb"))
	|| fwrite (&hdr, sizeof (struct dsr), 1, fp) != 1) errw (program, filespc);
	fclose (fp);

	exit (0);
}
