/**
 * @file  mri_interpolate.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/10/06 01:18:12 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "tags.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"

char *Progname ;

static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;

static int navgs = 50 ;

/*
   command line consists of two inputs:

   argv[1]  - input volume
   argv[2]  - output volume
*/


int
main(int argc, char *argv[]) {
  char         *in_fname, *out_fname, **av ;
  MRI          *mri_in, *mri_ctrl, *mri_out;
  int          ac, nargs=0, x, y, z ;
  char         cmdline[CMD_LINE_LEN] ;
  Real         val;

  make_cmd_version_string
    (argc, argv,
     "$Id: mri_interpolate.c,v 1.5 2011/10/06 01:18:12 fischl Exp $",
     "$Name:  $",cmdline);
  setRandomSeed(-1L) ;
  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (argc < 2)
    ErrorExit
    (ERROR_BADPARM,
     "usage: %s [<options>] <in volume> <out volume>\n", Progname) ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  mri_in=MRIread(in_fname);
  MRIaddCommandLine(mri_in, cmdline) ;

  mri_ctrl = MRIclone(mri_in,NULL);
  for (x = 0 ; x < mri_in->width; x++)
    for (y = 0 ; y < mri_in->height; y++)
      for (z = 0 ; z < mri_in->depth; z++)
      {
        val = MRIgetVoxVal(mri_in, x, y, z, 0);
        if (!FZERO(val))
          MRIsetVoxVal(mri_ctrl,x,y,z, 0, CONTROL_MARKED);
      }
  mri_out = MRIbuildVoronoiDiagram(mri_in, mri_ctrl, NULL) ;
  MRIsoapBubble(mri_out, mri_ctrl, mri_out, navgs, 1) ;
  MRIwrite(mri_out,out_fname);
  exit(0) ;
  return(0) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!strcmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else switch (*option) {
  case 'A':
    navgs = atoi(argv[2]);
    nargs = 1;
    printf("using %d soap bubble averages\n",navgs);
    break;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(int code) {
  printf("usage: %s <in volume> <out volume>\n",Progname) ;
  exit(code) ;
}
