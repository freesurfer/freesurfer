/**
 * @file  mris_thickness.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/09/17 00:44:24 $
 *    $Revision: 1.14 $
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

static char vcid[] = "$Id: mris_thickness.c,v 1.14 2007/09/17 00:44:24 fischl Exp $";

int main(int argc, char *argv[]) ;

int  MRISmeasureDistanceBetweenSurfaces(MRI_SURFACE *mris, MRI_SURFACE *mris2, int signed_dist) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static char pial_name[100] = "pial" ;
static char white_name[100] = WHITE_MATTER_NAME ;
static int write_vertices = 0 ;

static int nbhd_size = 2 ;
static float max_thick = 5.0 ;
static char *osurf_fname = NULL ;
static int signed_dist = 0 ;
static char sdir[STRLEN] = "" ;

int
main(int argc, char *argv[]) {
  char          **av, *out_fname, *sname, *cp, fname[STRLEN], *hemi ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_thickness.c,v 1.14 2007/09/17 00:44:24 fischl Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit() ;

  sname = argv[1] ;
  hemi = argv[2] ;
  out_fname = argv[3] ;
  if (!strlen(sdir)) {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }


#if 0
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAY_MATTER_NAME) ;
#else
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, pial_name) ;
#endif
  if (!FileExists(fname))
    sprintf(fname, "%s/%s/surf/%s.gray", sdir, sname, hemi) ;

  fprintf(stderr, "reading gray matter surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  if (osurf_fname) {
    MRI_SURFACE *mris2 ;
    mris2 = MRISread(osurf_fname) ;
    if (mris2 == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read 2nd surface from %s", Progname, osurf_fname) ;
    MRISmeasureDistanceBetweenSurfaces(mris, mris2, signed_dist) ;
    fprintf(stderr, "writing surface distance to curvature file %s...\n", out_fname) ;
    MRISwriteCurvature(mris, out_fname) ;
    exit(0) ;
  }

  if (MRISreadOriginalProperties(mris, white_name) != NO_ERROR)
    ErrorExit(Gerror, "%s: could not read white matter surface", Progname) ;
  fprintf(stderr, "measuring gray matter thickness...\n") ;

  if (write_vertices) {
    MRISfindClosestOrigVertices(mris, nbhd_size) ;
  } else {
    MRISmeasureCorticalThickness(mris, nbhd_size, max_thick) ;
  }

#if 0
  sprintf(fname, "%s/%s/surf/%s", sdir, sname, out_fname) ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
#endif
  fprintf(stderr, "writing %s to curvature file %s...\n",
          write_vertices ? "vertex correspondence" :
          "thickness", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;
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
    print_usage() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "pial")) {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "reading pial surface from file named %s\n", pial_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "white")) {
    strcpy(white_name, argv[2]) ;
    fprintf(stderr,  "reading white matter surface from file named %s\n", white_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max")) {
    max_thick = atof(argv[2]) ;
    fprintf(stderr,  "limiting maximum cortical thickness to %2.2f mm.\n",
            max_thick) ;
    nargs = 1 ;
  } else if (!stricmp(option, "osurf")) {
    osurf_fname = argv[2] ;
    signed_dist = 0 ;
    fprintf(stderr,  "measuring distance between input surface and %s\n", osurf_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "nsurf")) {
    osurf_fname = argv[2] ;
    signed_dist = 1 ;
    fprintf(stderr,  "measuring signed distance between input surface and %s\n", osurf_fname) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'V':
      write_vertices = 1 ;
      printf("writing vertex correspondences instead of thickness\n") ;
      nargs =  0 ;
      break ;
    case 'N':
      nbhd_size = atoi(argv[2]) ;
      fprintf(stderr, "using neighborhood size=%d\n", nbhd_size) ;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <subject name> <hemi> <thickness file>\n",
          Progname) ;
  print_help() ;
}

static void
print_help(void) {
  fprintf(stderr,
          "\nThis program measures the thickness of the cortical surface and\n"
          "writes the resulting scalar field into a 'curvature' file "
          "<thickness file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-max <max>\t use <max> to threshold thickness (default=5mm)\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

#include "mrishash.h"
#define MAX_DIST 10
int
MRISmeasureDistanceBetweenSurfaces(MRI_SURFACE *mris, MRI_SURFACE *mris2, int signed_dist) {
  int    vno ;
  VERTEX *v1, *v2 ;
  MRIS_HASH_TABLE *mht ;
  double           dx, dy, dz ;

  mht = MHTfillVertexTableRes(mris2, NULL, CURRENT_VERTICES, MAX_DIST) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v1 = &mris->vertices[vno] ;
    if (v1->ripflag)
      continue ;
    v2 = MHTfindClosestVertex(mht, mris2, v1) ;
    if (v2 == NULL) {
      v1->curv = MAX_DIST ;
      continue ;
    }
    dx = v1->x-v2->x ;
    dy = v1->y-v2->y ;
    dz = v1->z-v2->z ;
    v1->curv = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (signed_dist) {
      double dot ;

      dot = dx*v1->nx + dy*v1->ny + dz*v1->nz ;
      if (dot < 0)
        v1->curv *= -1 ;
    }
  }

  MHTfree(&mht) ;
  return(NO_ERROR) ;
}

