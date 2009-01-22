/**
 * @file  mris_thickness.c
 * @brief program for computing thickness of the cerebral cortex from 
 *  previously generated surfaces
 *
 * See (Fischl and Dale, 2000, PNAS)
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2009/01/22 02:47:49 $
 *    $Revision: 1.15 $
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
#include "icosahedron.h"

static char vcid[] = "$Id: mris_thickness.c,v 1.15 2009/01/22 02:47:49 fischl Exp $";

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
static char *sphere_name = "sphere" ;
static int signed_dist = 0 ;
static char sdir[STRLEN] = "" ;
static int new_thick = 0 ;
static INTEGRATION_PARMS parms ;

int
main(int argc, char *argv[]) {
  char          **av, *out_fname, *sname, *cp, fname[STRLEN], *hemi ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_thickness.c,v 1.15 2009/01/22 02:47:49 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  // for variational thickness estimation
  parms.dt = 0.2 ;
  parms.momentum = .5;
  parms.l_nlarea = 1 ;
  parms.l_thick_min = 1 ;
  parms.l_thick_parallel = 1;
  parms.tol = 1e-1 ;

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

  if (new_thick)
  {
    char              *cp, surf_fname[STRLEN], fname[STRLEN] ; ;
    MRI_SURFACE       *mris_ico ;
    
    if (parms.base_name[0] == 0) {
      
      FileNameOnly(out_fname, fname) ;
      cp = strchr(fname, '.') ;
      if (cp)
        strcpy(parms.base_name, cp+1) ;
      else
        strcpy(parms.base_name, "sphere") ;
      cp = strrchr(parms.base_name, '.') ;
      if (cp)
        *cp = 0 ;
    }

    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    if (MRISreadVertexPositions(mris, sphere_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, sphere_name) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;

    // read in icosahedron data (highly tessellated one)
    cp = getenv("FREESURFER_HOME");
    if (cp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: FREESURFER_HOME not defined in environment", cp) ;
    sprintf(surf_fname,"%s/lib/bem/ic7.tri",cp);
    mris_ico = MRISread(surf_fname) ;
    if (!mris_ico)
      ErrorExit(ERROR_NOFILE, "%s: could not open surface file %s",Progname, surf_fname) ;
    MRISscaleBrain(mris_ico, mris_ico, mris->radius/mris_ico->radius) ;
    MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char tmp[STRLEN] ;
      FileNameRemoveExtension(out_fname, tmp) ;
      
      sprintf(fname, "%s.correspondence.init", tmp) ;
      printf("writing initial correspondences to %s\n", fname) ;
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRISwrite(mris, fname) ;

      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
    }
    MRISminimizeThicknessFunctional(mris, &parms, max_thick, mris_ico) ;

    if (Gdiag & DIAG_WRITE)
    {
      char tmp[STRLEN] ;
      FileNameRemoveExtension(out_fname, tmp) ;
      
      sprintf(fname, "%s.correspondence.final", tmp) ;
      printf("writing final correspondences to %s\n", fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRISwrite(mris, fname) ;
    }
      
  }
  else if (write_vertices) {
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
  } else if (!stricmp(option, "THICK_PARALLEL")) {
    parms.l_thick_parallel = atof(argv[2]) ;
    printf("using parallel thickness coefficient %2.3f\n", parms.l_thick_parallel);
    nargs = 1 ;
  } else if (!stricmp(option, "THICK_MIN")) {
    parms.l_thick_min = atof(argv[2]) ;
    printf("using min length thickness coefficient %2.3f\n", parms.l_thick_min);
    nargs = 1 ;
  } else if (!stricmp(option, "DT")) {
    parms.dt = atof(argv[2]) ;
    printf("using time step %2.3f\n", parms.dt);
    nargs = 1 ;
  } else if (!stricmp(option, "tol")) {
    parms.tol = atof(argv[2]) ;
    printf("using tol %e\n", parms.tol);
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
  } else if (!stricmp(option, "new")) {
    new_thick = 1 ;
    fprintf(stderr,  "using variational thickness measurement\n") ;
  } else if (!stricmp(option, "nsurf")) {
    osurf_fname = argv[2] ;
    signed_dist = 1 ;
    fprintf(stderr,  "measuring signed distance between input surface and %s\n", osurf_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "vno")) {
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,  "debugging vertex %d\n", Gdiag_no) ;
  } else switch (toupper(*option)) {
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    Gdiag |= DIAG_WRITE ;
    printf("setting write iterations to %d\n", parms.write_iterations) ;
    break ;
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
    case 'M':
      parms.momentum = atof(argv[2]) ;
      fprintf(stderr, "using momentum %2.3f\n", parms.momentum) ;
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

