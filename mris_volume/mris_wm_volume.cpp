/**
 * @brief computes the interior volume of the ?h.white surfaces 
 * (excluding subcortical labels)
 *
 * computes the interior volume of the ?h.white surfaces (excluding 
 * subcortical labels), uses the aseg.mgz and the ?h.white surfaces.
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;

static void usage_exit(int code) ;
static char sdir[STRLEN] = "" ;
static const char *white_name = "white" ;
static const char *aseg_name = "aseg.mgz" ;

static double resolution = 1.0/4.0 ;

int
main(int argc, char *argv[]) 
{
  char          **av, *hemi, fname[STRLEN], *cp, *subject;
  int           ac, nargs, eno, nv, nf, ne;
  MRI_SURFACE   *mris;
  int           msec, minutes, seconds;
  Timer start ;
  double        total_volume ;
  MRI           *mri_aseg ;

  nargs = handleVersionOption(argc, argv, "mris_wm_volume");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(sdir)) 
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  subject = argv[1] ;
  hemi = argv[2] ;
  if (argc < 3)
    usage_exit(1) ;

  //  printf("command line parsing finished\n");

  /*** Read in the input surfaces ***/
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, white_name) ;
  mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;

  eno = MRIScomputeEulerNumber(mris, &nv, &nf, &ne) ;
  if (eno != 2)
    ErrorExit(ERROR_BADPARM, "%s: surface %s has an incorrect topology (eno=%d)",
              Progname, fname, eno) ;
  /*** Read in the aseg volume ***/
  sprintf(fname, "%s/%s/mri/%s", sdir, subject, aseg_name) ;
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, 
              "%s: could not read segmentation volume %s", 
              Progname, fname) ;

  total_volume = MRIScomputeWhiteVolume(mris, mri_aseg, resolution) ;
  printf("%g\n", total_volume);
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;


  MRISfree(&mris);
  MRIfree(&mri_aseg) ;

  exit(0);
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

  if (!stricmp(option, "white"))
  {
    white_name = argv[2] ;
    fprintf(stderr, "using %s as white surface name\n", white_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "aseg"))
  {
    aseg_name = argv[2] ;
    fprintf(stderr, "using %s as aseg volume name\n", aseg_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    fprintf(stderr, "using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  }
  else switch (*option) {
  default:
    printf("unknown option %s\n", argv[1]);
    exit(1);
    break;
  }

  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s <subject> <hemi>\n", Progname) ;
  printf("\t This program computes the volume of the enclosed ?h.white "
         "surface,\n\tignoring non-wm voxels in the aseg.\n");
  printf("\t use -v option to output more messages\n");
  printf("  -SDIR SUBJECTS_DIR \n");
  printf("  -white whitesurfname \n");
  printf("  -aseg  asegname \n");
  exit(code) ;
}

