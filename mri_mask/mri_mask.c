#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"

void usage(int exit_val);
MRI *fliplr(MRI *src);

char *Progname;


int main(int argc, char *argv[])
{

  MRI *mri_src, *mri_mask, *mri_dst ;
  int nargs;

  Progname = argv[0];

  nargs = handle_version_option (argc, argv, "$Id: mri_mask.c,v 1.3 2003/03/28 18:14:47 kteich Exp $");
  argc -= nargs ;
  if (1 == argc)
    exit (0);

  if(argc != 4)
    usage(1);

	mri_src = MRIread(argv[1]) ;
	if (!mri_src)
		ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
							Progname, argv[1]) ;
	mri_mask = MRIread(argv[2]) ;
	if (!mri_mask)
		ErrorExit(ERROR_BADPARM, "%s: could not read mask volume %s",
							Progname, argv[1]) ;
	mri_dst = MRImask(mri_src, mri_mask, NULL, 0, 0) ;
	if (!mri_dst)
		ErrorExit(Gerror, "%s: stripping failed", Progname) ;

	printf("writing masked volume to %s...\n", argv[3]) ;
  MRIwrite(mri_dst, argv[3]);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <in vol> <mask vol> <out vol>\n", Progname);
  fprintf(fout, "this program applies a mask volume (typically skull stripped)\n") ;

  exit(exit_val);

}  /*  end usage()  */
/*  EOF  */
