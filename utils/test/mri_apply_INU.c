/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


/* Try to estimate gain field from a pair of volumes went through NU_correct;
 *  and then apply the gain field to the current input volume
 */

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

static char *fname_before = NULL; /* filename for template volume before N3 */
static char *fname_after = NULL; /* filename for template volume after N3 */


void usage(int exit_val);

static int debug_flag = 0;

const char *Progname;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{

  char **av;
  MRI *mri_before, *mri_after, *mri_in, *mri_out;
  int ac, nargs;
  int  width, height, depth, x, y, z, f,nframes ;
  double v_before, v_after, v_in, v_out;
  double gain;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_apply_INU");
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }


  if (argc  !=  3)
    usage(1);

  if (fname_before == NULL  || fname_after == NULL)
  {
    printf("Use options to specify template volumes for gain field computation\n");
    usage(1);
  }

  mri_in = MRIread(argv[1]) ;
  if (!mri_in)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s",
              Progname, argv[1]) ;

  mri_before = MRIread(fname_before) ;
  if (!mri_before)
    ErrorExit(ERROR_BADPARM, "%s: could not read tempate volume %s",
              Progname, fname_before) ;

  mri_after = MRIread(fname_after) ;
  if (!mri_after)
    ErrorExit(ERROR_BADPARM, "%s: could not read tempate volume %s",
              Progname, fname_after) ;

  if ((mri_in->width != mri_after->width) ||
      (mri_in->height != mri_after->height) ||
      (mri_in->depth != mri_after->depth) ||
      (mri_in->width != mri_before->width) ||
      (mri_in->height != mri_before->height) ||
      (mri_in->depth != mri_before->depth)

     )
    ErrorExit(ERROR_BADPARM, "%s: three input volumes have different sizes \n", Progname);

  mri_out = MRIclone(mri_in, NULL) ;

  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;
  nframes = mri_in->nframes ;
  if (nframes == 0) nframes = 1;

  for (f = 0 ; f < nframes ; f++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          v_in = (double) MRIgetVoxVal(mri_in,x,y,z,f);
          v_before = (double) MRIgetVoxVal(mri_before,x,y,z,f);
          v_after = (double) MRIgetVoxVal(mri_after,x,y,z,f);
          gain = v_after/(v_before + 1e-15);
          v_out = gain*v_in;

          MRIsetVoxVal(mri_out,x,y,z,f,(float)v_out);
        }
      }
    }
  }

  printf("writing INU-corrected volume to %s...\n", argv[2]) ;
  MRIwrite(mri_out, argv[2]);

  MRIfree(&mri_in);
  MRIfree(&mri_out);
  MRIfree(&mri_before);
  MRIfree(&mri_after);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <input vol> <corrected vol>\n", Progname);
  fprintf(fout, "this program estimates gain field from volumes given in options, and apply it to input \n") ;
  fprintf(fout, "Options are (actually they have to be used): \n") ;
  fprintf(fout, "\t -before %%s: template before nu_correct \n") ;
  fprintf(fout, "\t -after %%s: template before nu_correct \n") ;

  exit(exit_val);

}  /*  end usage()  */
/*  EOF  */

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    debug_flag = 1;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "before"))
  {
    fname_before  = argv[2];
    printf("Use file %s for template volume before N3 \n", fname_before) ;
    nargs = 1;
  }
  else if (!stricmp(option, "after"))
  {
    fname_after  = argv[2];
    printf("Use file %s for template volume after N3 \n", fname_after) ;
    nargs = 1;
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      usage(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
