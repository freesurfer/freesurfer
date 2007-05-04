/*$Header: /space/repo/1/dev/dev/talairach_avi/conc.h,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: conc.h,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.6  2006/09/23  22:59:22  avi
 * add split() prototype
 *
 * Revision 1.5  2006/09/23  05:07:03  avi
 * eliminate prototypes for subroutines moved elsewhere
 *
 * Revision 1.4  2006/09/23  04:20:53  avi
 * provision for endian control
 *
 * Revision 1.3  2006/08/02  23:57:56  avi
 * newly defined conc_rewind()
 *
 * Revision 1.2  2005/09/16  03:34:13  avi
 * conc_init_quiet() conc_open_quiet()
 *
 * Revision 1.1  2004/11/27  05:44:52  avi
 * Initial revision
 **/

#define MAXL	256

typedef struct {
	char	program[MAXL];
	char	lstroot[MAXL], lstfile[MAXL];
	char	outroot[MAXL], outfile[MAXL];
	char	control;
	int	imgdim[4], vdim, orient, isbig, osbig;
	float	voxdim[3], mmppix[3], center[3];
	FILE	*lstfp, *imgfp, *outfp;
	char	**imgfile0;		/* input  4dfp filenames */
	char	**imgfile1;		/* output 4dfp filenames */
	int	rivol, wivol, *nvol;
	int	rifile, wifile, rnfile, wnfile, imgfp_open;
} CONC_BLOCK;

int split		(char *string, char *srgv[], int maxp);
void conc_init_quiet	(CONC_BLOCK *conc_block, char *program);
void conc_init		(CONC_BLOCK *conc_block, char *program);
void conc_open_quiet	(CONC_BLOCK *conc_block, char *lstfile);
void conc_open		(CONC_BLOCK *conc_block, char *lstfile);
void conc_read_vol	(CONC_BLOCK *conc_block, float *imgt);
void conc_write_vol	(CONC_BLOCK *conc_block, float *imgt);
void conc_new		(CONC_BLOCK *conc_block, char *trailer);
void conc_newe		(CONC_BLOCK *conc_block, char *trailer, char control);
int conc_ifh_hdr_rec	(CONC_BLOCK *conc_block, int argc, char *argv[], char *recstr);
void conc_free		(CONC_BLOCK *conc_block);
void conc_rewind	(CONC_BLOCK *conc_block);
