/**
 * @file  mris_euler_number.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:10 $
 *    $Revision: 1.5 $
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

static char vcid[] = "$Id: mris_euler_number.c,v 1.5 2006/12/29 02:09:10 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static float curv_thresh = 2.0f ;
static int patch_flag = 0 ;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, fname[100] ;
  int          ac, nargs, nvertices, nfaces, nedges, eno, dno ;
  MRI_SURFACE  *mris ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_euler_number.c,v 1.5 2006/12/29 02:09:10 nicks Exp $", "$Name:  $");
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

  if (argc < 2)
    usage_exit() ;

  in_fname = argv[1] ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
  fprintf(stderr, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
          nvertices, nedges, nfaces, eno, 1-eno/2) ;

  fprintf(stderr, "      F =2V-4:          %d %s= %d-4 (%d)\n",
          nfaces, nfaces == 2*nvertices-4 ? "" : "!", 2*nvertices,
          2*nvertices-4-nfaces) ;
  fprintf(stderr, "      2E=3F:            %d %s= %d (%d)\n",
          2*nedges, 2*nedges == 3*nfaces ? "" : "!", 3*nfaces,
          2*nedges-3*nfaces) ;

  dno = MRIStopologicalDefectIndex(mris) ;
  fprintf(stderr, "\ntotal defect index = %d\n", dno) ;

  if (patch_flag) {
    MRISremoveTopologicalDefects(mris, curv_thresh) ;
    fprintf(stderr, "\nafter editing:\n") ;

    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stderr, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;

    fprintf(stderr, "      F =2V-4:          %d %s= %d-4 (%d)\n",
            nfaces, nfaces == 2*nvertices-4 ? "" : "!", 2*nvertices,
            2*nvertices-4-nfaces) ;
    fprintf(stderr, "      2E=3F:            %d %s= %d (%d)\n",
            2*nedges, 2*nedges == 3*nfaces ? "" : "!", 3*nfaces,
            2*nedges-3*nfaces) ;

    dno = MRIStopologicalDefectIndex(mris) ;
    fprintf(stderr, "total defect index = %d\n", dno) ;

    sprintf(fname, "%s.edit", in_fname) ;
    fprintf(stderr, "writing out patched surface to %s\n", fname) ;
    MRISwritePatch(mris, fname) ;
  }
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
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
    case 'P':
      patch_flag = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      curv_thresh = (float)atof(argv[2]) ;
      nargs = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr, "usage: %s [options] <input surface file>\n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will compute the euler number of a cortical surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

