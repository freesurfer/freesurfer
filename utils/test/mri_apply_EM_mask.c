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


/* Clear voxels labelled as dura by mri_ms_EM (hard seg)
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

static char *fname_dura = NULL; /* filename for dura membership function */

// static float threshold = 40;

void usage(int exit_val);

static int debug_flag = 0;

const char *Progname;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{

  char **av;
  MRI *mri_dura, *mri_in, *mri_out;
  int ac, nargs;
  int  width, height, depth, x, y, z, f,nframes ;
  double v_in, v_out;
  int dura_label;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_apply_EM_mask");
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

  if (fname_dura == NULL)
  {
    printf("Use options to specify dura membership function volume \n");
    usage(1);
  }

  mri_in = MRIread(argv[1]) ;
  if (!mri_in)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s",
              Progname, argv[1]) ;

  mri_dura = MRIread(fname_dura) ;
  if (!mri_dura)
    ErrorExit(ERROR_BADPARM, "%s: could not read dura volume %s",
              Progname, fname_dura) ;

  if ((mri_in->width != mri_dura->width) ||
      (mri_in->height != mri_dura->height) ||
      (mri_in->depth != mri_dura->depth)
     )
    ErrorExit(ERROR_BADPARM, "%s: the input volumes have different sizes \n", Progname);

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
          dura_label = (int) MRIgetVoxVal(mri_dura,x,y,z,f);
          if (dura_label == 2) v_out = 0;
          else v_out = v_in;

          MRIsetVoxVal(mri_out,x,y,z,f,(float)v_out);
        }
      }
    }
  }

  printf("writing dura-removed volume to %s...\n", argv[2]) ;
  MRIwrite(mri_out, argv[2]);

  MRIfree(&mri_in);
  MRIfree(&mri_out);
  MRIfree(&mri_dura);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <input vol> <corrected vol>\n", Progname);
  fprintf(fout, "this program estimates dura from given membership function, and clear it from input volume \n") ;
  fprintf(fout, "Options are (actually -dura has to be used): \n") ;
  fprintf(fout, "\t -dura %%s: membership function for dura \n") ;

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
  else if (!stricmp(option, "dura"))
  {
    fname_dura  = argv[2];
    printf("Use file %s for hard EM-segmentation volume (dura_label = 2) \n", fname_dura) ;
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
