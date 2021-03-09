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


/* Extract a subvolume from original one
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

void usage(int exit_val);

static int debug_flag = 0;

static int start_flag = 0;
static int size_flag = 0;

static int start_x, start_y, start_z, dx, dy, dz;

const char *Progname;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{

  char **av;
  MRI *mri_in, *mri_out;
  int ac, nargs;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_extract");
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

  if (start_flag == 0)
  {
    printf("Use -start # # # flasg to specify starting coordinates for cropping\n");
    exit(0);
  }
  if (size_flag == 0)
  {
    printf("use -size # # # flasg to specify desired final size\n");
    exit(0);
  }

  mri_in = MRIread(argv[1]) ;
  if (!mri_in)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s",
              Progname, argv[1]) ;

  mri_out = MRIextract(mri_in, NULL, start_x, start_y, start_z, dx, dy, dz) ;

  printf("writing cropped volume to %s...\n", argv[2]) ;
  MRIwrite(mri_out, argv[2]);

  MRIfree(&mri_in);
  MRIfree(&mri_out);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <input vol> <out vol>\n", Progname);
  fprintf(fout, "this program crops the input volume to a new size \n") ;
  fprintf(fout, "Options are (actually they have to be used): \n") ;
  fprintf(fout, "\t -start (x, y, z): starting point \n") ;
  fprintf(fout, "\t -size (width, height, depth): final size \n") ;

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
  else if (!stricmp(option, "start"))
  {
    start_x = atoi(argv[2]) ;
    start_y = atoi(argv[3]) ;
    start_z = atoi(argv[4]) ;
    start_flag = 1;
    nargs = 3 ;
    printf("starting voxel is (%d, %d, %d)...\n", start_x, start_y, start_z) ;
  }
  else if (!stricmp(option, "size"))
  {
    dx = atoi(argv[2]) ;
    dy = atoi(argv[3]) ;
    dz = atoi(argv[4]) ;
    size_flag = 1;
    nargs = 3 ;
    printf("output size is (%d, %d, %d)...\n", dx, dy, dz) ;

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
