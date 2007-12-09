/**
 * @file  mris_density.c
 * @brief fills the interior of a surface.
 *
 * program to fill the interior of a surface at an arbitrary resolution
 * specified by the user.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/09 22:45:02 $
 *    $Revision: 1.2.2.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

static const char vcid[] = 
"$Id: mris_fill.c,v 1.2.2.1 2007/12/09 22:45:02 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static double resolution = .25 ;
static int conform = 0 ;

int
main(int argc, char *argv[]) 
{
  char          **av, *out_fname, *in_fname ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  MRI           *mri_interior ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_fill.c,v 1.2.2.1 2007/12/09 22:45:02 nicks Exp $", "$Name:  $",
   cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_fill.c,v 1.2.2.1 2007/12/09 22:45:02 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 3)
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  fprintf(stderr, "reading surface from %s...\n", in_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, 
              "%s: could not read surface file %s",
              Progname, in_fname) ;

  mri_interior = MRISfillInterior(mris, resolution, NULL) ;

  if (conform)
  {
    MRI *mri_tmp, *mri_tmp2 ;
    mri_tmp = MRIalloc(256, 256, 256, MRI_UCHAR) ;
    MRIsetResolution(mri_tmp, 1.0, 1.0, 1.0) ;
    mri_tmp->xstart = mri_tmp->ystart = mri_tmp->zstart = -mri_tmp->width/2;
    mri_tmp->xend =   mri_tmp->yend =   mri_tmp->zend =   mri_tmp->width/2;
    mri_tmp->x_r = -1.0;  mri_tmp->x_a =  0.0; mri_tmp->x_s =  0.0;
    mri_tmp->y_r =  0.0;  mri_tmp->y_a =  0.0; mri_tmp->y_s = -1.0;
    mri_tmp->z_r =  0.0;  mri_tmp->z_a =  1.0;  mri_tmp->z_s =  0.0;
    mri_tmp->c_r = mris->vg.c_r ;
    mri_tmp->c_a = mris->vg.c_a ;
    mri_tmp->c_s = mris->vg.c_s ;
    mri_tmp2 = MRIresample(mri_interior, mri_tmp, SAMPLE_NEAREST);
    MRIfree(&mri_interior) ; MRIfree(&mri_tmp) ;
    mri_interior = mri_tmp2 ;
  }
  MRIaddCommandLine(mri_interior, cmdline) ;
  fprintf(stderr, "writing filled volume to %s...\n", out_fname) ;
  MRIwrite(mri_interior, out_fname) ;
  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version")){
    print_version() ;
  } else switch (toupper(*option)) {
  case 'R':
    resolution = (double)atof(argv[2]) ;
    nargs = 1 ;
    printf("setting resolution for intermediate calculations to %2.4f\n", 
           resolution) ;
    break ;
  case 'C':
    conform = 1 ;
    printf("conforming output volume\n") ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void) 
{
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) 
{
  printf("usage: %s [options] <input surface> <output volume>\n",
         Progname) ;
}

static void
print_help(void) 
{
  print_usage() ;
  printf("\nThis program floodfills the interior of a surface and writes\n"
         "the results into a volume of the specified resolution.\n") ;
  printf("\nvalid options are:\n\n") ;
  printf("\t-r <resolution>: set the resolution of the output volume"
         " (default = %2.3f mm/voxel)\n", resolution) ;
  printf("\t-c               'conform' the volume before writing\n") ;
  exit(1) ;
}

static void
print_version(void) 
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
