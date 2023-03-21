/**
 * @brief expand a surface outwards by a specified amount
 *
 * Expands a surface (typically ?h.white) outwards while maintaining smoothness
 * and self-intersection constraints.
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
#include "macros.h"
#include "error.h"
#include "fsinit.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "randomfields.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static INTEGRATION_PARMS parms ;
static int use_thickness = 0 ;
static int nsurfaces = 1 ;
static const char *thickness_name = "thickness" ;
static const char *pial_name = "pial" ;
static char *tmap_fname = NULL ;
static char *tmap_write_fname = NULL ;
static float tmap_std = 0.0 ;
static float tmap_min = .25 ;
static float tmap_max = .75 ;
static int tmap_avgs = 0 ;
static int nbrs = 2 ;

static char *orig_name = NULL ;
LABEL *label=NULL;

int
main(int argc, char *argv[])
{
  int          nargs ;
  char         *in_fname, *out_fname, fname[STRLEN], *cp ;
  int          msec, minutes, seconds ;
  Timer start ;
  float        mm_out ;
  MRI_SURFACE  *mris ;

  parms.l_spring = .05;
  parms.l_location = 1 ;
  // parms.l_curv = 1.0 ;
  parms.n_averages = 16 ;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.min_averages = 0 ;
  parms.l_surf_repulse = .0 ;
  parms.l_repulse = 0.025 ;
  parms.dt = 0.25 ;

  nargs = handleVersionOption(argc, argv, "mris_expand");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  FSinit() ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 4)
    usage_exit(1) ;

  in_fname = argv[1] ;
  mm_out = atof(argv[2]) ;
  out_fname = argv[3] ;

  FileNameOnly(out_fname, fname) ;
  cp = strchr(fname, '.') ;
  if (cp)
    strcpy(parms.base_name, cp+1) ;
  else
    strcpy(parms.base_name, "expanded") ;

  if (Gdiag & DIAG_WRITE)
  {
    char log_fname[STRLEN] ;
    int req = snprintf(log_fname, STRLEN, "%s.log", fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    INTEGRATION_PARMS_openFp(&parms, log_fname, "w") ;
    if (parms.fp)
      printf("writing log results to %s\n", log_fname) ;
    else
      ErrorExit(ERROR_BADPARM, "%s: could not open log file %s", Progname, log_fname) ;
  }
  if (use_thickness)
  {
    if (use_thickness > 0)
      printf("expanding surface %s by %2.1f%% of thickness "
	     "and writing it to %s\n",
	     in_fname, 100.0*mm_out, out_fname) ;
    else
      printf("expanding surface %s by percentages in %s "
	     "and writing it to %s\n",
	     in_fname, tmap_fname, out_fname) ;
  }
  else
    printf("expanding surface %s by %2.1f mm and writing it to %s\n",
           in_fname, mm_out, out_fname) ;
  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: MRISread(%s) failed", Progname, in_fname);

  if(label)
    LabelRip(mris,label,1);
  
  if (nbrs > 1)
  {
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
  }
  if (use_thickness)
  {
    printf("reading thickness %s...\n", thickness_name) ;
    if (MRISreadCurvatureFile(mris, thickness_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not load thickness file", Progname) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    if (MRISreadVertexPositions(mris, pial_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read pial vertex positions\n",
                Progname) ;
    if (MRISreadCanonicalCoordinates(mris, "sphere") != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "") ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    MRISripZeroThicknessRegions(mris) ;
  }
  MRIScomputeMetricProperties(mris) ;
  if (orig_name)
    MRISreadOriginalProperties(mris, orig_name) ;
  else
    MRISstoreMetricProperties(mris) ;

  if (parms.mri_brain && FZERO(mris->vg.xsize))
    MRIScopyVolGeomFromMRI(mris, parms.mri_brain) ;
  if (tmap_fname != NULL)
  {
    parms.mri_dtrans = MRIread(tmap_fname);
    if (parms.mri_dtrans == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read tmap vertex percentages from %s\n", Progname, tmap_fname) ;
    if (parms.mri_dtrans->width != mris->nvertices)
      ErrorExit(ERROR_NOFILE,
                "%s: tmap width %d != mris->nvertices %d in %s\n", 
		Progname, parms.mri_dtrans->width, mris->nvertices, tmap_fname) ;
  }
  else if (!FZERO(tmap_std))    // create a random map of distances
  {
    RFS    *rfs = RFspecInit(0, NULL);
    int    vno ;

    parms.mri_dtrans = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;
    if (parms.mri_dtrans == NULL)
      ErrorExit(ERROR_NOFILE,
                "%s: could not allocate tmap vertex percentages from %s\n", Progname, tmap_fname) ;

    printf("creating random tmap distances\n") ;
    rfs->name = strcpyalloc("gaussian");
    rfs->params[0] = mm_out ;         // mean is halfway through the ribbon
    rfs->params[1] = tmap_std;     // std

    for (vno = 0 ; vno < mris->nvertices ; vno++)
      mris->vertices[vno].val = RFdrawVal(rfs);

    MRISaverageVals(mris, tmap_avgs) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      float dist = mris->vertices[vno].val ;

      dist = MIN(MAX(dist, tmap_min), tmap_max) ;
      if (mris->vertices[vno].ripflag)
	dist = 0 ;

      if (vno == Gdiag_no)
	printf("vno %d: val %f, ripflag %d\n", vno, dist, mris->vertices[vno].ripflag) ;

      mris->vertices[vno].val = dist ;
      MRIsetVoxVal(parms.mri_dtrans, vno, 0, 0, 0, dist) ;
    }
    if (tmap_write_fname != NULL)
    {
      printf("writing random tmap to %s\n", tmap_write_fname) ;
      MRISwriteValues(mris, tmap_write_fname) ;
    }
  }
  else
    parms.mri_dtrans = NULL ;
  MRISexpandSurface(mris, mm_out, &parms, use_thickness, nsurfaces);
#if 0
  if (navgs > 0)
    MRISaverageVertexPositions(mris, navgs) ;
#endif
  if (nsurfaces == 1)
  {
    printf("writing expanded surface to %s...\n", out_fname) ;
    MRISwrite(mris, out_fname) ;
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stdout, "surface expansion took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;

  exit(0) ;
  return(0) ;
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "thickness"))
  {
    use_thickness = 1 ;
    printf("using distance as a %% of thickness\n") ;
  }
  else if (!stricmp(option, "label"))
  {
    label = LabelRead(NULL,argv[2]);
    if(label==NULL) exit(1);
    printf("Ripping vertices outside of label %s\n",argv[2]);
    nargs = 1;
  }
  else if (!stricmp(option, "convex"))
  {
    sscanf(argv[2], "%f", &parms.l_convex) ;
    nargs = 1 ;
    fprintf(stderr, "using l_convex = %2.3f\n", parms.l_convex) ;
  }
  else if (!stricmp(option, "norm"))
  {
    sscanf(argv[2], "%f", &parms.l_norm) ;
    nargs = 1 ;
    fprintf(stderr, "using l_norm = %2.3f\n", parms.l_norm) ;
  }
  else if (!stricmp(option, "max_spring"))
  {
    sscanf(argv[2], "%f", &parms.l_max_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_max_spring = %2.3f\n", parms.l_max_spring) ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    fprintf(stderr,  "using neighborhood size = %d\n", nbrs) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "thickness_name"))
  {
    thickness_name = argv[2] ;
    printf("using thickness file %s\n", thickness_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "navgs"))
  {
    parms.n_averages = atof(argv[2]) ;
    parms.min_averages = atoi(argv[3]) ;
    printf("using n_averages %d --> %d\n", parms.n_averages, parms.min_averages) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "intensity"))
  {
    parms.target_intensity = atof(argv[2]) ;
    parms.mri_brain = MRIread(argv[3]) ;
    if (parms.mri_brain == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load target intensity volume %s\n", argv[3]) ;
    printf("cropping target locations to at intensity %2.0f in %s\n", parms.target_intensity, argv[3]) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "curv"))
  {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "location"))
  {
    parms.l_location = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_location = %2.3f\n", parms.l_location) ;
  }
  else if (!stricmp(option, "nspring"))
  {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  }
  else if (!stricmp(option, "angle"))
  {
    parms.l_angle = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_angle = %2.3f\n", parms.l_angle) ;
  }
  else if (!stricmp(option, "pangle"))
  {
    parms.l_pangle = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_angle = %2.3f\n", parms.l_pangle) ;
  }
  else if (!stricmp(option, "spring_norm"))
  {
    parms.l_spring_norm = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_spring_norm = %2.3f\n", parms.l_spring_norm) ;
  }
  else if (!stricmp(option, "nltspring"))
  {
    parms.l_nltspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nltspring = %2.3f\n", parms.l_nltspring) ;
  }
  else if (!stricmp(option, "tspring"))
  {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  }
  else if (!stricmp(option, "surf_repulse"))
  {
    parms.l_surf_repulse = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_surf_repulse = %2.3f\n", parms.l_surf_repulse) ;
  }
  else if (!stricmp(option, "wd"))
  {
    tmap_write_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "writing random tmap to %s\n", tmap_write_fname) ;
  }
  else if (!stricmp(option, "pial"))
  {
    pial_name = argv[2] ;
    printf("reading pial surface from %s\n", pial_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "tmap"))
  {
    use_thickness = -1 ;
    if (!stricmp(argv[2], "random"))
    {
      tmap_std = atof(argv[3]);
      tmap_min = atof(argv[4]);
      tmap_max = atof(argv[5]);
      tmap_avgs = atoi(argv[6]);
      printf("creating random tmap in [%2.2f, %2.2f] with std %2.2f and %d averages\n",
	     tmap_min, tmap_max, tmap_std, tmap_avgs) ;
      nargs = 5 ;
    }
    else
    {
      tmap_fname = argv[2] ;
      printf("reading thickness target percent map from %s\n", tmap_fname) ;
      nargs = 1 ;
    }
  }
  else switch (toupper(*option))
    {
    case 'O':
      orig_name = argv[2] ;
      nargs = 1 ;
      printf("reading original metric properties from %s\n", orig_name) ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      break ;
    case 'R':
      parms.l_repulse = atof(argv[2]) ;
      fprintf(stderr, "l_repulse = %2.3f\n", parms.l_repulse) ;
      nargs = 1 ;
      break ;
    case 'S':
      parms.l_spring = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting spring term to %2.2f\n", parms.l_spring) ;
      break ;
    case 'T':
      parms.dt = atof(argv[2]) ;
      nargs = 1 ;
      printf("setting dt = %2.2f\n", parms.dt) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("writing snapshots of expansion every %d iterations\n",
             parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'A':
      parms.smooth_averages = atoi(argv[2]) ;
      nargs = 1 ;
      printf("smoothing surface with %d iterations after expansion\n",
             parms.smooth_averages) ;
      break ;
    case 'N':   // how many surfaces to write out
      nsurfaces = atoi(argv[2]) ;
      nargs = 1 ;
      printf("writing %d surfaces during expansion\n", nsurfaces) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("Usage: %s [options] <input surface> <mm> <output surface>\n",
         Progname) ;
  printf("  Example: mris_expand -thickness lh.white 0.5 lh.graymid\n");
  printf("  Example: mris_expand -label labelfile lh.white 0.5 lh.graymid\n");
  printf("  Example: mris_expand -tmap thickness_pct_target.mgz lh.white 0.5 lh.graymid\n");
  printf("     use a prespecified map of percent thickness to compute the target locations for expansion\n");
  printf("  Example: mris_expand -tmap random 2 .25 .75 100 -wd tmap.mgz lh.white 0.5 lh.graymid\n");
  printf("     creates a random target distance map with gaussian sampling (mean=.5, std=2) and cropping to .25/.75\n");
  printf("     and spatial averaging 100 times (cropping is average averaging). The map will be written to \n");
  printf("     file tmap.mgz\n");
  exit(code) ;
}
