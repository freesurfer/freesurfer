/**
 * @File  mris_thickness.c
 * @brief program for computing thickness of the cerebral cortex from 
 *  previously generated surfaces
 *
 * See (Fischl and Dale, 2000, PNAS)
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
 * Reporting: freesurfer@nmr.mgh.harvard.include
	 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "timer.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "icosahedron.h"
#include "label.h"


int main(int argc, char *argv[]) ;

static int fill_thickness_holes(MRI_SURFACE *mris, LABEL *cortex_label, LABEL *fsaverage_label) ;
int  MRISmeasureDistanceBetweenSurfaces(MRI_SURFACE *mris, MRI_SURFACE *mris2, int signed_dist) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
static char pial_name[100] = "pial" ;
static char white_name[100] = WHITE_MATTER_NAME ;
static int write_vertices = 0 ;

static int nbhd_size = 2 ;
static float max_thick = 5.0 ;
static char *osurf_fname = NULL ;
static const char *sphere_name = "sphere" ;
static int signed_dist = 0 ;
static char sdir[STRLEN] = "" ;
static int fmin_thick = 0 ;
static float laplace_res = 0.5 ;
static int laplace_thick = 0 ;
static INTEGRATION_PARMS parms ;

static char *long_fname = NULL ;

static LABEL *cortex_label = NULL ;
static LABEL *fsaverage_label = NULL ;
MRI *GetThickFromSeg(MRIS *surf, LABEL *label, MRI *vol, double dmax, double ddelta, MRI* thick);

#include "voxlist.h"
#include "mrinorm.h"
int
main(int argc, char *argv[]) {
  char          *out_fname, *sname, *cp, fname[STRLEN], *hemi ;
  int           nargs, msec ;
  MRI_SURFACE   *mris ;
  Timer then ;

  nargs = handleVersionOption(argc, argv, "mris_thickness");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  then.reset() ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  // for variational thickness estimation
  parms.dt = 0.1 ;
  parms.remove_neg = 1 ;
  parms.momentum = .1; 
  parms.niterations = 1000 ;
  parms.l_nlarea = 0 ;
  parms.l_thick_min = 1 ;
  parms.l_thick_spring = 0 ;
  parms.l_ashburner_triangle = 1 ;
  parms.l_ashburner_lambda = .1 ;
//  parms.l_tspring = .25;
  parms.l_thick_normal = 1;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.tol = 1e-3 ;

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

  int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", sdir, sname, hemi, pial_name) ;
  if (req >= STRLEN) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading pial surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;

  MRISresetNeighborhoodSize(mris, nbhd_size) ;
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

  if (laplace_thick)
  {
    MRI *mri_laplace ;

    {
      MRI *mri_dist_white, *mri_dist_pial ;
      int nwmissing, npmissing, nmissing, vno ;
      double xv, yv, zv, white_val, pial_val ;
      VERTEX *v ;
      FILE   *fp ;

      MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
      mri_dist_pial = MRIScomputeDistanceToSurface(mris, NULL, laplace_res) ;

      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;  // actually white
      MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
      mri_dist_white = MRIScomputeDistanceToSurface(mris, NULL, laplace_res) ;
      for (nmissing = nwmissing = npmissing = vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        MRISsurfaceRASToVoxel(mris, mri_dist_pial, v->pialx, v->pialy, v->pialz, &xv, &yv, &zv);
        MRIsampleVolumeFrameType(mri_dist_pial, xv, yv, zv, 0, SAMPLE_NEAREST, &pial_val) ;
        MRISsurfaceRASToVoxel(mris, mri_dist_white, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv);
        MRIsampleVolumeFrameType(mri_dist_white, xv, yv, zv, 0, SAMPLE_NEAREST, &white_val) ;
        if (fabs(white_val) > laplace_res)
          nwmissing++ ;
        if (fabs(pial_val) > laplace_res)
          npmissing++ ;
        if ((fabs(pial_val) > laplace_res) || (fabs(white_val) > laplace_res))
          nmissing++ ;
      }
      printf("%d of %d pial surface nodes not resolved - %2.3f %%\n",
             npmissing, mris->nvertices, 100.0*npmissing/mris->nvertices) ;
      printf("%d of %d gray/white surface nodes not resolved - %2.3f %%\n",
             nwmissing, mris->nvertices, 100.0*nwmissing/mris->nvertices) ;
      MRIfree(&mri_dist_white) ; MRIfree(&mri_dist_pial) ;
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      fp = fopen("laplace_missing.txt", "a") ;
      fprintf(fp, "%f %d %d %f %d %f %d %f\n", 
              laplace_res,  mris->nvertices,
              nwmissing, 100.0*nwmissing/mris->nvertices,
              npmissing, 100.0*npmissing/mris->nvertices,
              nmissing, 100.0*nmissing/mris->nvertices) ;
      fclose(fp) ;
    }
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
    mri_laplace = MRISsolveLaplaceEquation(mris, NULL, laplace_res) ;
    MRISmeasureLaplaceStreamlines(mris, mri_laplace, NULL, NULL) ;
  }
  else if (fmin_thick)
  {
    char              *cp, surf_fname[STRLEN], fname[STRLEN] ; ;
    
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

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      char tmp[STRLEN] ;
      FileNameRemoveExtension(out_fname, tmp) ;
      
      int req = snprintf(fname, STRLEN, "%s.correspondence.init", tmp) ;
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing initial correspondences to %s\n", fname) ;
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
      MRISwrite(mris, fname) ;

      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      MRIScomputeMetricProperties(mris) ;
    }
    if (Gdiag & DIAG_WRITE)
    {
      char tmp[STRLEN] ;
      int  vno ;
      VERTEX *v ;
      FileNameRemoveExtension(out_fname, tmp) ;
      

      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
      MRISsaveVertexPositions(mris, TMP_VERTICES) ;
      MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
      MRISfindClosestPialVerticesCanonicalCoords(mris, mris->nsize) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
        {
          v->nx = v->ny = v->nz = 0 ;
          continue ;
        }
        v->nx = v->x - v->whitex ; 
        v->ny = v->y - v->whitey ; 
        v->nz = v->z - v->whitez ; 
      }
      int req = snprintf(fname, STRLEN, "%s.normals.init.mgz", tmp) ;
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing initial surface normals to %s\n", fname) ;
      MRISwriteNormals(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;  // restore spherical canonical coords
    }
    
    //    MRISripZeroThicknessRegions(mris) ;
    MRISstoreRipFlags(mris) ;
#if 0
    {
      static INTEGRATION_PARMS nparms ;

      nparms.dt = 0.9 ;
      nparms.momentum = .5; 
      nparms.niterations = 1000 ;
      nparms.l_nlarea = 0 ;
      nparms.l_thick_min = 0 ;
      nparms.l_thick_parallel = 0;
      nparms.l_thick_normal = 1;
      nparms.l_tspring = 0;
      nparms.tol = 4e-1 ;
      MRISminimizeThicknessFunctional(mris, &nparms, max_thick) ;
      printf("------- normal vector field computed - minimizing full functional -----------\n") ;
      parms.start_t = nparms.start_t ;
    }
#endif
    MRISminimizeThicknessFunctional(mris, &parms, max_thick) ;

    if (Gdiag & DIAG_WRITE)
    {
      char tmp[STRLEN] ;
      int  vno ;
      VERTEX *v ;
      FileNameRemoveExtension(out_fname, tmp) ;
      
      int req = snprintf(fname, STRLEN, "%s.normals.mgz", tmp) ;
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing final surface normals to %s\n", fname) ;
      MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
      MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
        {
          v->nx = v->ny = v->nz = 0 ;
          continue ;
        }
        v->nx = v->tx - v->whitex ; 
        v->ny = v->ty - v->whitey ; 
        v->nz = v->tz - v->whitez ; 
      }
      MRISwriteNormals(mris, fname) ;
      MRISrestoreVertexPositions(mris, TMP2_VERTICES) ;
    }
    if (long_fname)
    {
      char line[STRLEN], subject[STRLEN], fname[STRLEN], base_name[STRLEN], *cp,
        tmp[STRLEN], out_fname_only[STRLEN] ;
      int  vno ;
      FILE *fp ;
      VERTEX *v ;
      MHT   *mht ;

      MRIScopyCurvatureToImagValues(mris) ; // save base thickness 
      // to lookup closest face
      mht = MHTcreateFaceTable_Resolution(mris, CANONICAL_VERTICES, 1.0); 

      fp = fopen(long_fname, "r") ;
      if (fp == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not open timepoint file %s\n", Progname,long_fname) ;
      strcpy(tmp, long_fname) ;
      cp = strrchr(tmp, '/') ;
      if (cp == NULL)
        ErrorExit(ERROR_BADPARM, "could not read trailing / from %s", tmp) ;
      *cp = 0 ;
      cp = strrchr(tmp, '/') ;
      if (cp == NULL)
        cp = tmp-1 ;
      strcpy(base_name, cp+1) ;
      do
      {
        if (fgetl(line, STRLEN-1, fp) == NULL)
          break ;
        sscanf(line, "%s", subject) ;
        printf("processing longitudinal subject %s\n", subject) ;
        int req = snprintf(fname, STRLEN, "%s/%s.long.%s/surf/%s.%s", sdir, subject, base_name, hemi,
                white_name) ;
        if (req >= STRLEN) {
          std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
        }
        if (MRISreadWhiteCoordinates(mris, fname) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
        req = snprintf(fname, STRLEN, "%s/%s.long.%s/surf/%s.%s", sdir, subject, base_name, hemi,
                pial_name) ;
        if (req >= STRLEN) {
          std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
        }
        if (MRISreadPialCoordinates(mris, fname) != NO_ERROR)
          ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, fname) ;
        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          float   xw, yw, zw, xp, yp, zp, thick ;
          
          v = &mris->vertices[vno] ;
          if (vno == Gdiag_no)
            DiagBreak() ;
          if (v->ripflag)
          {
            v->tx = v->whitex ; v->ty = v->whitey ; v->tz = v->whitez ;
            continue ;
          }
          MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw) ;
          MRISsampleFaceCoordsCanonical(mht, mris, v->x, v->y, v->z, 
                                        PIAL_VERTICES, &xp, &yp, &zp);
          thick = sqrt(SQR(xp-xw) + SQR(yp-yw) + SQR(zp-zw)) ;
          v->curv = thick ; v->tx = xp ; v->ty = yp ; v->tz = zp ;
        }
        FileNameOnly(out_fname, out_fname_only) ;
        req = snprintf(fname, STRLEN, "%s/%s.long.%s/surf/%s", sdir, subject, base_name, out_fname_only);
        if (req >= STRLEN) {
          std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
        }
        printf("writing thickness estimate to %s\n", fname) ;
        MRISwriteCurvature(mris, fname) ;
      } while (strlen(line) > 0) ;
      MHTfree(&mht) ;
      MRIScopyImagValuesToCurvature(mris) ; // restore base  thickness 
    }
  }
  else if (write_vertices) {
    MRISfindClosestOrigVertices(mris, nbhd_size) ;
  } else {
    MRISmeasureCorticalThickness(mris, nbhd_size, max_thick) ;
  }
 
  if (cortex_label)  // fill in thickness in holes in the cortex label where it isn't to be trusted
  {
    fill_thickness_holes(mris, cortex_label, fsaverage_label) ;
  }


#if 0
  sprintf(fname, "%s/%s/surf/%s", sdir, sname, out_fname) ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
#endif
  fprintf(stderr, "writing %s to curvature file %s...\n",
          write_vertices ? "vertex correspondence" :
          "thickness", out_fname) ;
  MRISwriteCurvature(mris, out_fname) ;
  msec = then.milliseconds() ;
  fprintf(stderr,"thickness measurement took %2.1f minutes\n", (float)msec/(60*1000.0f));
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
  else if (!stricmp(option, "fill_holes"))
  {
    cortex_label = LabelRead(NULL, argv[2]) ;
    if (cortex_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", argv[2]) ;
    fsaverage_label = LabelRead(NULL, argv[3]) ;
    if (fsaverage_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", argv[3]) ;
    nargs = 2 ;
    printf("filling holes in label %s that are not present in %s\n", argv[2], argv[3]) ;

  }
  else if (!stricmp(option, "long"))
  {
    long_fname = argv[2] ;
    nargs = 1 ;
    printf("computing longitudinal thickness from time points found in %s\n",
           long_fname) ;
  }
  else if (!stricmp(option, "pial")) {
    strcpy(pial_name, argv[2]) ;
    fprintf(stderr,  "reading pial surface from file named %s\n", pial_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "optimal")) {
    parms.integration_type = INTEGRATE_LM_SEARCH ;
    fprintf(stderr,  "using line search minimization\n") ;
  } else if (!stricmp(option, "momentum")) {
    parms.integration_type = INTEGRATE_MOMENTUM ;
    fprintf(stderr,  "using gradient descent with momentum minimization\n") ;
  } else if (!stricmp(option, "ic")) {
    parms.ico_order = atoi(argv[2]) ;
    fprintf(stderr,  "using %dth order icosahedron\n", parms.ico_order) ;
    nargs = 1 ;
  } else if (!stricmp(option, "SDIR")) {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "THICK_PARALLEL")) {
    parms.l_thick_parallel = atof(argv[2]) ;
    printf("using parallel thickness coefficient %2.3f\n", parms.l_thick_parallel);
    nargs = 1 ;
  } else if (!stricmp(option, "THICK_SPRING")) {
    parms.l_thick_spring = atof(argv[2]) ;
    printf("using spring thickness coefficient %2.3f\n", parms.l_thick_spring);
    nargs = 1 ;
  } else if (!stricmp(option, "THICK_MIN")) {
    parms.l_thick_min = atof(argv[2]) ;
    printf("using min length thickness coefficient %2.3f\n", parms.l_thick_min);
    nargs = 1 ;
  } else if (!stricmp(option, "spring")) {
    parms.l_tspring = atof(argv[2]) ;
    printf("using spring coefficient %2.3f\n", parms.l_tspring);
    nargs = 1 ;
  } else if (!stricmp(option, "neg")) {
    parms.remove_neg = 0 ;
    printf("allowing negative  vertices to occur during integration\n") ;
  } else if (!stricmp(option, "noneg")) {
    parms.remove_neg = 1 ;
    printf("not allowing negative  vertices to occur during integration\n") ;
  } else if (!stricmp(option, "nlarea")) {
    parms.l_nlarea = atof(argv[2]) ;
    printf("using nonlinear area coefficient %2.3f\n", parms.l_nlarea);
    nargs = 1 ;
  } else if (!stricmp(option, "triangle")) {
    parms.l_ashburner_triangle = atof(argv[2]) ;
    parms.l_ashburner_lambda = atof(argv[3]) ;
    printf("using Ashburner, 1999, triangle regularization with l=%2.2f and lambda=%2.1f\n", 
           parms.l_ashburner_triangle, parms.l_ashburner_lambda);
    nargs = 2 ;
  } else if (!stricmp(option, "THICK_NORMAL")) {
    parms.l_thick_normal = atof(argv[2]) ;
    printf("using normal thickness coefficient %2.3f\n", parms.l_thick_normal);
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
  } else if (!stricmp(option, "new") || !stricmp(option, "fmin") || !stricmp(option, "variational")  || !stricmp(option, "vector")) {
    fmin_thick = 1 ;
    fprintf(stderr,  "using variational thickness measurement\n") ;
  } else if (!stricmp(option, "laplace") || !stricmp(option, "laplacian")) {
    laplace_thick = 1 ;
    laplace_res = atof(argv[2]) ;
    fprintf(stderr,  "using Laplacian thickness measurement with PDE resolution = %2.3fmm\n",laplace_res) ;
    nargs = 1 ;
  } else if (!stricmp(option, "nsurf")) {
    osurf_fname = argv[2] ;
    signed_dist = 1 ;
    fprintf(stderr,  "measuring signed distance between input surface and %s\n", osurf_fname) ;
    nargs = 1 ;
  } 
  else if (!stricmp(option, "vno")) {
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,  "debugging vertex %d\n", Gdiag_no) ;
  } 
  else if (!stricmp(option, "thickness-from-seg")) {
    // surf label seg.mgz dmaxmm (eg, 6) ddeltamm (eg, .01) output.mgz
    MRIS *surf = MRISread(argv[2]);
    if(!surf) exit(1);
    LABEL *mylabel = LabelRead(NULL,argv[3]);
    if(!mylabel) exit(1);
    MRI *seg = MRIread(argv[4]);
    if(!seg) exit(1);
    double dmax,ddelta;
    sscanf(argv[5],"%lf",&dmax);
    sscanf(argv[6],"%lf",&ddelta);
    MRI *mri2 = GetThickFromSeg(surf, mylabel, seg, dmax, ddelta, NULL);
    MRIwrite(mri2,argv[7]);
    exit(0);
  } 
  else switch (toupper(*option)) {
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
  fprintf(stderr, "-fill_holes <cortex label> <fsaverage cortex label> fill in thickness in holes in the cortex label\n");
  printf("\n");
  printf("-thickness-from-seg surf label seg.mgz dmaxmm (eg, 6) ddeltamm (eg, .01) output.mgz\n");
  printf("  This is a stand-alone option that allows thickness to be computed from a surface created around\n");
  printf("  a segmentation by tracing a ray from one side to the other along the normal to the surface.\n");
  printf("  The calculations are limited to the label. The output is a surface overlay with the values\n");
  printf("  of the thickness at each vertex in the label. dmax is the maximum distance searched.\n");
  printf("  ddelta is the step size.\n");
  printf("\n\n") ;
  printf("-vector  compute the thickness using a variationally derived vector field instead of shortest distance\n") ;

  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
fill_thickness_holes(MRI_SURFACE *mris, LABEL *cortex_label, LABEL *fsaverage_label)
{
  int    vno ;
  VERTEX *v ;

  LabelMarkSurface(cortex_label, mris) ;
  LabelAddToMark(fsaverage_label, mris, 2) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked <= 1)  // not in fsaverage surface
    {
      v->marked = 0 ;   // not marked in either
      continue ;
    }

    v->val = v->curv ;  // put thickness into val field for soap bubble below
    if (v->marked == 2)  // marked only in fsaverage surface, don't make it a fixed point
      v->marked = 0 ;
    else
      v->marked = 1 ;
  }

  MRISsoapBubbleVals(mris, 500); 
  MRISclearMarks(mris) ;
  LabelMarkSurface(cortex_label, mris) ;
  LabelAddToMark(fsaverage_label, mris, 2) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked < 2)
      v->val = 0 ;   // fix everything not in the fsaverage label to have 0 thickness
  }
  MRIScopyValuesToCurvature(mris) ;
  return(NO_ERROR) ;
}

MRI *GetThickFromSeg(MRIS *surf, LABEL *label, MRI *seg, double dmax, double ddelta, MRI* thick)
{
  if(thick==NULL){
    thick = MRIalloc(surf->nvertices, 1, 1, MRI_FLOAT);
  }
  printf("GetThickFromLabel(): dmax = %g, ddelta = %g\n",dmax,ddelta);

  MATRIX *CRS, *RAS;
  CRS = MatrixAlloc(4,1,MATRIX_REAL);
  CRS->rptr[4][1] = 1;
  RAS = MatrixAlloc(4,1,MATRIX_REAL);
  RAS->rptr[4][1] = 1;
  MATRIX *vox2tkras = MRIxfmCRS2XYZtkreg(seg);
  MATRIX *tkras2vox = MatrixInverse(vox2tkras,NULL);
  MatrixFree(&vox2tkras);

  int n;
  double dsum = 0;
  if(label){
    for(n=0; n < label->n_points; n++){
      int vno = label->lv[n].vno;
      VERTEX *v;
      v = &(surf->vertices[vno]);
      double d;
      for(d=1; d <= dmax; d += ddelta){
	RAS->rptr[1][1] = v->x - d*v->nx;
	RAS->rptr[2][1] = v->y - d*v->ny;
	RAS->rptr[3][1] = v->z - d*v->nz;
	CRS = MatrixMultiply(tkras2vox,RAS,CRS);
	int c,r,s;
	c = nint(CRS->rptr[1][1]);
	r = nint(CRS->rptr[2][1]);
	s = nint(CRS->rptr[3][1]);
	if(c < 0 || c >= seg->width) continue;
	if(r < 0 || r >= seg->height) continue;
	if(s < 0 || s >= seg->depth) continue;
	double val = MRIgetVoxVal(seg,c,r,s,0);
	if(val > 0.5){
	  //hit = 1;
	  continue;
	}
	// val <= 0.5 -- this voxel not a hit
	break;
      }
      MRIsetVoxVal(thick,vno,0,0,0,d);
      dsum += d;
    }
    printf("%g %d   %g\n",dsum,label->n_points,dsum/label->n_points);
  }
  return(thick);
}
