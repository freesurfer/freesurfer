/**
 * @file  label2patch.c
 * @brief convert a label into a patch suitable for flattening
 *
 * convert a label into a patch suitable for flattening
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/03/23 14:37:42 $
 *    $Revision: 1.4 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void print_usage(void) ;

char *Progname ;

static int verbose = 0 ;


static char subjects_dir[STRLEN] = "" ;

static int nerode = 0 ;
static int ndilate = 0 ;
static int nclose = 0 ;


int
main(int argc, char *argv[]) {
  int          ac, nargs ;
  char         **av, *cp, surf_name[100], *hemi, *subject_name, *label_name,
               *out_fname ;
  MRI_SURFACE  *mris ;
  LABEL        *label ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: label2patch.c,v 1.4 2011/03/23 14:37:42 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* read in command-line options */
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    print_usage() ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  label_name = argv[3] ;
  out_fname = argv[4] ;

  if (strlen(subjects_dir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "no subjects directory in environment.\n") ;
    strcpy(subjects_dir, cp) ;
  }

  sprintf(surf_name,"%s/%s/surf/%s.inflated",subjects_dir,subject_name,hemi);
  fprintf(stderr, "reading %s...\n", surf_name) ;
  mris = MRISread(surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s\n",
              surf_name) ;
  MRIScomputeMetricProperties(mris) ;
  label = LabelRead(subject_name, label_name) ;
  if (ndilate)
    LabelDilate(label, mris, ndilate) ;
  if (nerode)
    LabelErode(label, mris, nerode) ;
  if (nclose)
  {
    LabelDilate(label, mris, nclose) ;
    LabelErode(label, mris, nclose) ;
  }

  LabelRipRestOfSurface(label, mris) ;
  MRISripFaces(mris) ;
  MRISwritePatch(mris, out_fname) ;
  if (verbose)
    fprintf(stderr, "done.\n") ;
  exit(0) ;
  return(0) ;  /* ansi */
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
  if (!stricmp(option, "dilate"))
  {
    ndilate = atoi(argv[2]) ;
    nargs = 1;
    printf("dilating label %d times...\n", ndilate) ;
  }
  else if (!stricmp(option, "erode"))
  {
    nerode = atoi(argv[2]) ;
    nargs = 1;
    printf("eroding label %d times...\n", nerode) ;
  }
  else if (!stricmp(option, "close"))
  {
    nclose = atoi(argv[2]) ;
    nargs = 1;
    printf("closing label %d times...\n", nclose) ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(subjects_dir, argv[2]) ;
    printf("using SUBJECTS_DIR = %s\n", subjects_dir) ;
    nargs = 1;
  }
  else switch (toupper(*option)) {
  case 'V':
    verbose = !verbose ;
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
print_usage(void) {
  printf("usage: %s [options] <subject name> <hemi> <label file name> <output patch file>...\n",Progname);
  printf("where valid options are:\n") ;
  printf("\t-dilate <n>  : dilate the label <n> times before creating the patch\n") ;
  printf("\t-erode <n>   : erode the label <n> times before creating the patch\n") ;
  printf("\t-close <n>   : close the label <n> times before creating the patch\n") ;
  printf("\t-sdir <path> : use <path> as the SUBJECTS_DIR instead of the environment variable\n") ;
  exit(1) ;
}
