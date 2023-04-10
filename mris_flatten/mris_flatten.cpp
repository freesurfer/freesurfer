/**
 * @brief flatten a surface patch
 *
 * "Cortical Surface-Based Analysis II: Inflation, Flattening, and a
 * Surface-Based Coordinate System", Fischl, B., Sereno, M.I., Dale, A.M.
 * (1999) NeuroImage, 9(2):195-207.
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "cma.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "utils.h"
#include "version.h"
#include "fastmarching.h"
#include "mri2.h"
#include "mrishash.h"


int main(int argc, char *argv[]) ;


static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
MRIS *SurfCopyCoords = NULL;

const char *Progname ;

static char *synth_name = NULL ;
static INTEGRATION_PARMS  parms ;
#define BASE_DT_SCALE     1.0
static float base_dt_scale = BASE_DT_SCALE ;
static char *label_fname = NULL ;
static int nbrs = 2 ;
static int do_inflate = 0 ;
static double disturb = 0 ;
static int randomly_flatten = 0 ;
static int   nospring = 0 ;
static float scale = 3 ;
static int   max_passes = 1 ;

static int sphere_flag = 0 ;
static int plane_flag = 0 ;
static int dilate = 0 ;
static int dilate_label = 0 ; // how many times to dilate label after reading

static int one_surf_flag = 0 ;
static const char *original_surf_name = SMOOTH_NAME ;
static const char *original_unfold_surf_name = ORIG_NAME ;
static float rescale = 1.0f ;

static MRI *mri_overlay ;  // if "flattening" an overlay with an existing flatmap

static LABEL *label_overlay = NULL ;

static double
rectangle_error(MRI_SURFACE *mris, double xmin, double ymin, double xmax, double ymax)
{
  int    vno ;
  VERTEX *v ;
  double max_dist, x0, y0, x, y, dx, dy, xtop, ytop, xbottom, ybottom ;

  // find point that is closest to (or outside) of the boundary
  // negative distances mean inside the rectangle
  // note this is only for a rectangle aligned with the cardinal axes
  x0 = (xmin + xmax) / 2 ; y0 = (ymin + ymax) / 2 ; // center of rectangle
  xtop = xmax - x0 ; ytop = ymax - y0 ;
  xbottom = xmin - x0 ; ybottom = ymin - y0 ;
  max_dist = -MAX(xmax-xmin, ymax-ymin) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    x = v->x - x0 ; y = v->y - y0 ;
    if (x > 0)
      dx = x - xtop ;
    else
      dx = xbottom - x ;
    if (y > 0)
      dy = y - ytop ;
    else
      dy = ybottom - y ;
    if (vno == 0 || dx > max_dist)
      max_dist = dx ;
    if (dy > max_dist)
      max_dist = dy ;
  }
  return(max_dist) ;
}

static int
find_biggest_inscribed_rectangle(MRI_SURFACE *mris, double *pxmin, double *pymin, double *pxmax, double *pymax)
{
  double   x0, y0, xmin, ymin, xmax, ymax, error ;
  int      nv, vno ;
  VERTEX   *v ;

  xmin = ymin = 1000000 ;
  xmax = ymax = -xmin ;
  for (x0 = y0 = 0.0, nv = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    nv++ ;
    x0 += v->x ; y0 += v->y ;
    if (v->x > xmax)
      xmax = v->x ;
    if (v->y > ymax)
      ymax = v->y ;
    if (v->x < xmin)
      xmin = v->x ;
    if (v->y < ymin)
      ymin = v->y ;
  }
  x0 /= nv ; y0 /= nv ;
  error = rectangle_error(mris, xmin, ymin, xmax, ymax) ;
  if (DIAG_VERBOSE_ON)
    printf("centroid of flatmap found at (%2.1f, %2.1f), bounding box (%2.1f --> %2.1f, %2.1f --> %2.1f) : %2.2f\n", 
           x0, y0, xmin, xmax, ymin, ymax, error) ;
  *pxmax = xmax ; *pxmin = xmin ; *pymin = ymin ; *pymax = ymax ;
  return(NO_ERROR) ;
}
MRI *
MRIflattenOverlay(MRI_SURFACE *mris, MRI *mri_overlay, MRI *mri_flat, double res, LABEL *label_overlay,
                  MRI **pmri_vertices)
{
  double   xmin, ymin, xmax, ymax, fdist, lambda[3], xf, yf, val0, val1, val2, val ;
  int      width, height, x, y, fno, ret, z ;
  MHT      *mht ;
  FACE     *face;
  MRI      *mri_vertices ;

  find_biggest_inscribed_rectangle(mris, &xmin, &ymin, &xmax, &ymax) ;
  width = (int)ceil((xmax-xmin)/res) ; width = (int)(floor(width/2.0)*2.0+1) ;   xmax=xmin+width;
  height = (int)ceil((ymax-ymin)/res) ;height = (int)(floor(height/2.0)*2.0+1) ; ymax=ymin+height;

  // 1st frame is correlations and 2nd is vertex #
  mri_vertices = MRIalloc(width, height, 1, MRI_FLOAT) ;
  MRIsetValues(mri_vertices, -1) ;
  mri_flat = MRIalloc(width, height, mri_overlay->nframes, MRI_FLOAT) ;
  mri_vertices->xstart = mri_flat->xstart = xmin ; mri_vertices->xend = mri_flat->xend = xmax ;
  mri_vertices->ystart = mri_flat->ystart = ymin ; mri_vertices->yend = mri_flat->yend = ymax ;
  mri_vertices->zstart = mri_flat->zstart = 0 ; 
  mri_vertices->zend = mri_flat->zend = mri_overlay->nframes-1 ;
  mri_vertices->c_r = mri_flat->c_r = xmin ;  mri_vertices->c_a = mri_flat->c_a = ymin ; 
  mri_vertices->c_s = mri_flat->c_s = 0 ;
  MRIsetResolution(mri_flat, res, res, 1) ;
  MRIsetResolution(mri_vertices, res, res, 1) ;
  if (label_overlay)  // constrain processing to only this label
    LabelRipRestOfSurface(label_overlay, mris) ;
  mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, 1.0) ;
  for (x = 0 ; x < width; x++)
    for (y = 0 ; y < height ; y++)
    {
      xf = x*res + xmin ;  yf = y*res + ymin ;   // back to flattened coords
      MHTfindClosestFaceGeneric(mht, mris, xf, yf, 0.0, 10*res, 2, 1, &face, &fno, &fdist) ;
      if (fno >= 0)  // otherwise this point is not in a face
      {
        ret = face_barycentric_coords(mris, fno, CURRENT_VERTICES, xf, yf, 0, &lambda[0], &lambda[1],&lambda[2]); 
        if (ret >= 0)
        {
          if (lambda[0] > lambda[1])
          {
            if (lambda[0] > lambda[2])
            {
              if (face->v[0] == Gdiag_no)
                DiagBreak() ;
              MRIsetVoxVal(mri_vertices, x, y, 0, 0, face->v[0]) ;
            }
            else
            {
              if (face->v[2] == Gdiag_no)
                DiagBreak() ;
              MRIsetVoxVal(mri_vertices, x, y, 0, 0, face->v[2]) ;
            }
          }
          else
          {
            if (lambda[1] > lambda[2])
            {
              if (face->v[1] == Gdiag_no)
                DiagBreak() ;
              MRIsetVoxVal(mri_vertices, x, y, 0, 0, face->v[1]) ;
            }
            else
            {
              if (face->v[2] == Gdiag_no)
                DiagBreak() ;
              MRIsetVoxVal(mri_vertices, x, y, 0, 0, face->v[2]) ;
            }
          }

          for (z = 0 ;z < mri_flat->depth ; z++)
          {
            val0 = MRIgetVoxVal(mri_overlay, face->v[0], 0, 0, z) ;
            val1 = MRIgetVoxVal(mri_overlay, face->v[1], 0, 0, z) ;
            val2 = MRIgetVoxVal(mri_overlay, face->v[2], 0, 0, z) ;
            val = lambda[0]*val0 + lambda[1]*val1 + lambda[2]*val2 ;
            MRIsetVoxVal(mri_flat, x, y, z, 0, val) ;
          }
        }
        else if (fabs(xf) < 10 && fabs(yf) < 10)
        {
          MHTfindClosestFaceGeneric(mht, mris, xf, yf, 0.0, 1000, -1, 1, &face, &fno, &fdist) ;
          printf("(%d, %d) --> %f %f unmapped (goes to face %d, v (%d, %d, %d) if projected\n",
                 x, y, xf, yf, fno, face->v[0], face->v[1], face->v[2]) ;
          DiagBreak() ;
        }
      }
    }

  if (pmri_vertices)
    *pmri_vertices = mri_vertices ;
  MHTfree(&mht) ;
  return(mri_flat) ;
}


int
main(int argc, char *argv[])
{
  char         **av, in_surf_fname[STRLEN], *in_patch_fname, *out_patch_fname,
  fname[STRLEN], path[STRLEN], *cp, hemi[10] ;
  int          ac, nargs ;
  MRI_SURFACE  *mris ;
  MRI          *mri_vertices ;

  nargs = handleVersionOption(argc, argv, "mris_flatten");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= (DIAG_SHOW | DIAG_WRITE) ;
  // memset(&parms, 0, sizeof(parms)) ;
  parms.dt = .1 ;
  parms.projection = PROJECT_PLANE ;
  parms.tol = 0.2 ;
  parms.n_averages = 1024 ;
  parms.l_dist = 1.0 ;
  parms.l_nlarea = 1.0 ;
  parms.niterations = 40 ;
  parms.area_coef_scale = 1.0 ;
  parms.dt_increase = 1.01 /* DT_INCREASE */;
  parms.dt_decrease = 0.98 /* DT_DECREASE*/ ;
  parms.error_ratio = 1.03 /*ERROR_RATIO */;
  parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
  parms.momentum = 0.9 ;
  parms.desired_rms_height = -1.0 ;
  parms.base_name[0] = 0 ;
  parms.nbhd_size = 7 ;    /* out to 7-connected neighbors */
  parms.max_nbrs = 12 ;    /* 12 at each distance */
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    print_help() ;

  parms.base_dt = base_dt_scale * parms.dt ;
  in_patch_fname = argv[1] ;
  out_patch_fname = argv[2] ;
  FileNamePath(in_patch_fname, path) ;
  cp = strrchr(in_patch_fname, '/') ;
  if (!cp)
    cp = in_patch_fname ;
  cp = strchr(cp, '.') ;
  if (cp)
  {
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
  }
  else
    strcpy(hemi, "lh") ;
  if (one_surf_flag) {
    int req = snprintf(in_surf_fname, STRLEN, "%s", in_patch_fname) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  } else {
    int req = snprintf(in_surf_fname, STRLEN, "%s/%s.%s", path, hemi, original_surf_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
  }

  if (parms.base_name[0] == 0)
  {
    FileNameOnly(out_patch_fname, fname) ;
    cp = strchr(fname, '.') ;
    if (cp)
      strcpy(parms.base_name, cp+1) ;
    else
      strcpy(parms.base_name, "flattened") ;
  }

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;

  {
    int vno;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
      if (mris->vertices_topology[vno].vnum == 0)
	mris->vertices[vno].ripflag = 1 ;
  }
  if (sphere_flag)
  {
    MRIScenter(mris, mris) ;
    mris->radius = MRISaverageRadius(mris) ;
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  }

  if (Gdiag_no >= 0)
  {
    int n ;
    printf("vertex %d has %d nbrs before patch:\n",
           Gdiag_no, mris->vertices_topology[Gdiag_no].vnum) ;
    for (n = 0 ; n < mris->vertices_topology[Gdiag_no].vnum ; n++)
      printf("\t%d\n", mris->vertices_topology[Gdiag_no].v[n]) ;
  }
  if (one_surf_flag)  /* only have the 1 surface - no patch file */
  {
    mris->patch = 1 ;
    mris->status = MRIS_PATCH ;
    if (!FEQUAL(rescale,1))
    {
      MRISscaleBrain(mris, mris, rescale) ;
      MRIScomputeMetricProperties(mris) ;
    }
    MRISstoreMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

  } 
  else
  {
    MRISsetNeighborhoodSizeAndDist(mris, mris->vertices_topology[0].nsizeMax);
    MRISresetNeighborhoodSize(mris, mris->vertices_topology[0].nsizeMax) ; // set back to max
    if (label_fname) // read in a label instead of a patch
    {
      LABEL *area ;
      area = LabelRead(NULL, label_fname) ;
      if (area == NULL)
        ErrorExit(ERROR_BADPARM, "%s: could not read label file %s",
                  Progname, label_fname) ;

      LabelDilate(area, mris, dilate_label, CURRENT_VERTICES) ;
      MRISclearMarks(mris) ;
      LabelMark(area, mris) ;
      MRISripUnmarked(mris) ;
      MRISsetRipInFacesWithRippedVertices(mris);
      mris->patch = 1 ;
      mris->status = MRIS_CUT ;
      LabelFree(&area) ;
      printf("%d valid vertices (%2.1f %% of total)\n",
             MRISvalidVertices(mris), 
             100.0*MRISvalidVertices(mris)/mris->nvertices) ;
    }
    else
    {
      if (MRISreadPatch(mris, in_patch_fname) != NO_ERROR)
        ErrorExit(ERROR_BADPARM, "%s: could not read patch file %s",
                  Progname, in_patch_fname) ;
      if (dilate)
      {
        printf("dilating patch %d times\n", dilate) ;
        MRISdilateRipped(mris, dilate) ;
        printf("%d valid vertices (%2.1f %% of total)\n",
               MRISvalidVertices(mris), 100.0*MRISvalidVertices(mris)/mris->nvertices) ;
      }
    }
    MRISremoveRipped(mris) ;
    MRISupdateSurface(mris) ;
#if 0
    mris->nsize = 1 ; // before recalculation of 2 and 3-nbrs
    {
      int vno ;
      VERTEX *v ;
      for (vno= 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        v->vtotal = v->vnum ;
        v->nsize = 1 ;
      }
    }
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
#endif
  }
  if(SurfCopyCoords) MRIScopyCoords(mris,SurfCopyCoords);

  if (Gdiag_no >= 0)
    printf("vno %d is %sin patch\n", Gdiag_no,
           mris->vertices[Gdiag_no].ripflag ? "NOT " : "") ;

  if (Gdiag_no >= 0 && mris->vertices[Gdiag_no].ripflag == 0)
  {
    int n ;
    printf("vertex %d has %d nbrs after patch:\n",
           Gdiag_no, mris->vertices_topology[Gdiag_no].vnum) ;
    for (n = 0 ; n < mris->vertices_topology[Gdiag_no].vnum ; n++)
      printf("\t%d\n", mris->vertices_topology[Gdiag_no].v[n]) ;
  }
  fprintf(stderr, "reading original vertex positions...\n") ;
  if (!FZERO(disturb))
    mrisDisturbVertices(mris, disturb) ;
  if (parms.niterations > 0)
  {
    MRISresetNeighborhoodSize(mris, nbrs) ;

    if (!FZERO(parms.l_unfold) || !FZERO(parms.l_expand))
    {
      static INTEGRATION_PARMS p2 ;
      int req = snprintf(in_surf_fname, STRLEN, "%s/%s.%s", path, hemi, original_surf_name) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (stricmp(original_unfold_surf_name,"none") == 0)
      {
        printf("using current position of patch as initial position\n") ;
        MRISstoreMetricProperties(mris) ;  /* use current positions */
      }
      else if (!sphere_flag && !one_surf_flag)
        MRISreadOriginalProperties(mris, original_unfold_surf_name) ;
      INTEGRATION_PARMS_copy(&p2, &parms) ;
      p2.l_dist = 0 ;
      p2.niterations = 100 ;
      p2.nbhd_size = p2.max_nbrs = 1 ;
      p2.n_averages = 0 ;
      p2.write_iterations = parms.write_iterations > 0 ? 25 : 0 ;
      p2.tol = -1 ;
      p2.dt = 0.5 ;
      p2.l_area = 0.0 ;
      p2.l_spring = 0.9 ;
      p2.l_convex = 0.9 ;
      p2.momentum = 0 ;
      p2.integration_type = INTEGRATE_MOMENTUM ;
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
#if 0
      p2.flags |= IPFLAG_NO_SELF_INT_TEST ;
      printf("expanding surface....\n") ;
      MRISexpandSurface(mris, 4.0, &p2) ;  // push it away from fissure
#endif
      p2.niterations = 100 ;
      MRISunfold(mris, &p2, 1) ;
      p2.niterations = 300 ;
      p2.l_unfold *= 0.25 ;
      MRISunfold(mris, &p2, 1) ;
      p2.l_unfold *= 0.25 ;
      MRISunfold(mris, &p2, 1) ;
#if 0
      printf("smoothing unfolded surface..\n");
      p2.niterations = 200 ;
      p2.l_unfold = 0 ;  // just smooth it
      MRISunfold(mris, &p2, max_passes) ;
#endif
      parms.start_t = p2.start_t ;
      parms.l_unfold = parms.l_convex = parms.l_boundary = parms.l_expand=0 ;
      MRIfree(&parms.mri_dist) ;
    }

    int req = snprintf(in_surf_fname, STRLEN, "%s/%s.%s", path, hemi, original_surf_name) ;
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!sphere_flag && !one_surf_flag)
      MRISreadOriginalProperties(mris, original_surf_name) ;
    if (randomly_flatten)
      MRISflattenPatchRandomly(mris) ;
    else
      MRISflattenPatch(mris) ;

    /* optimize metric properties of flat map */
    fprintf(stderr,"minimizing metric distortion induced by projection...\n");
    MRISscaleBrain(mris, mris, scale) ;
    MRIScomputeMetricProperties(mris) ;
    MRISunfold(mris, &parms, max_passes) ;
    MRIScenter(mris, mris) ;
    fprintf(stderr, "writing flattened patch to %s\n", out_patch_fname) ;
    MRISwritePatch(mris, out_patch_fname) ;
  }

  if (plane_flag || sphere_flag)
  {
    char fname[STRLEN] ;
    FILE *fp ;

#if 0
    sprintf(fname, "%s.%s.out",
            mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
            parms.base_name);
#else
    sprintf(fname, "flatten.log") ;
#endif
    fp = fopen(fname, "a") ;

    if (plane_flag)
      MRIScomputeAnalyticDistanceError(mris, MRIS_PLANE, fp) ;
    else if (sphere_flag)
      MRIScomputeAnalyticDistanceError(mris, MRIS_SPHERE, fp) ;
    fclose(fp) ;
  }

  if (mri_overlay)
  {
    MRI  *mri_flattened ;
    char fname[STRLEN] ;

    // if it is NxNx1x1 reshape it to be Nx1x1xN
    if ( mri_overlay->width == mri_overlay->height &&
       mri_overlay->depth == 1 &&
       mri_overlay->nframes == 1)
    {
      MRI *mri_tmp ;
      printf("reshaping to move 2nd dimension to time\n") ;
      mri_tmp = mri_reshape( mri_overlay, mri_overlay->width, 1, 1, mri_overlay->height);
      MRIfree( &mri_overlay );
      mri_overlay = mri_tmp;
    }

    // put in some special code that knows about icosahedra
    if (mris->nvertices == 163842 ||  // ic7
        mris->nvertices == 40962 ||  // ic6
        mris->nvertices == 10242 ||  // ic5
        mris->nvertices == 2562)  // ic4
    {
      int nvals, start_index, end_index ;
      MRI *mri_tmp ;
      
      printf("cross-hemispheric correlation matrix detected, reshaping...\n") ;
      nvals = mri_overlay->width * mri_overlay->height * mri_overlay->depth ;
      if (nvals == 2*mris->nvertices)   // it's a corr matrix for both hemis
      {
        if (mris->hemisphere == LEFT_HEMISPHERE || mris->hemisphere == RIGHT_HEMISPHERE)
        {
          if (mris->hemisphere == LEFT_HEMISPHERE)
          {
            start_index = 0 ; 
            end_index = mris->nvertices-1 ;
          }
          else
          {
            start_index = mris->nvertices ; 
            end_index = 2*mris->nvertices-1 ;
          }
          mri_tmp = MRIextract(mri_overlay, NULL, start_index, 0, 0, mris->nvertices, 1, 1) ;
          MRIfree(&mri_overlay) ;
          mri_overlay = mri_tmp;
        }
        else // both hemis
        {
        }
      }
    }
    
    printf("resampling overlay (%d x %d x %d x %d) into flattened coordinates..\n",
           mri_overlay->width, mri_overlay->height, mri_overlay->depth, mri_overlay->nframes) ;
    if (synth_name)
    {
      LABEL *area_lh, *area_rh ;
      char  fname[STRLEN], path[STRLEN], fname_no_path[STRLEN] ;
      int   vno, n, vno2, n2 ;

      MRIsetValues(mri_overlay, 0) ;
      FileNameOnly(synth_name, fname_no_path) ;
      FileNamePath(synth_name, path) ;
      int req = snprintf(fname, STRLEN, "%s/lh.%s", path, fname_no_path) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      area_lh = LabelRead(NULL, fname) ;
      if (area_lh == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",
                  Progname,fname) ;
      req = snprintf(fname, STRLEN, "%s/rh.%s", path, fname_no_path) ;
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      area_rh = LabelRead(NULL, fname) ;
      if (area_rh == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",
                  Progname,fname) ;
#if 0
      for (n = 0 ; n < area_lh->n_points ; n++)
      {
        vno = area_lh->lv[n].vno ;
        MRIsetVoxVal(mri_overlay, vno, 0, 0, vno, 1) ;
	printf("synthesizing map with vno %d: (%2.1f, %2.1f)\n", vno, mris->vertices[vno].x, mris->vertices[vno].y) ;
        break ;
      }
#else
      for (n = 0 ; n < area_lh->n_points ; n++)
      {
        vno = area_lh->lv[n].vno ;
        if (vno >= 0)
        {
          for (n2 = 0 ; n2 < area_lh->n_points ; n2++)
          {
            vno2 = area_lh->lv[n2].vno ;
            if (vno2 >= 0)
              MRIsetVoxVal(mri_overlay, vno, 0, 0, vno2, 1) ;
          }
          for (n2 = 0 ; n2 < area_rh->n_points ; n2++)
          {
            vno2 = area_rh->lv[n2].vno ;
            if (vno2 >= 0)
              MRIsetVoxVal(mri_overlay, vno, 0, 0, mris->nvertices+vno2, 1) ;
          }
        }
      }
#endif
    }

    mri_flattened = MRIflattenOverlay(mris, mri_overlay, NULL, 1.0, label_overlay, &mri_vertices) ;
    printf("writing flattened overlay to %s\n", out_patch_fname) ;
    MRIwrite(mri_flattened, out_patch_fname) ;
    MRIfree(&mri_flattened) ;

    FileNameRemoveExtension(out_patch_fname, fname) ;
    strcat(fname, ".vnos.mgz") ;
    printf("writing flattened vertex #s to %s\n", fname) ;
    MRIwrite(mri_vertices, fname) ;
    MRIfree(&mri_vertices) ;
  }
#if 0
  sprintf(fname, "%s.area_error", out_fname) ;
  printf("writing area errors to %s\n", fname) ;
  MRISwriteAreaError(mris, fname) ;
  sprintf(fname, "%s.angle_error", out_fname) ;
  printf("writing angle errors to %s\n", fname) ;
  MRISwriteAngleError(mris, fname) ;
  MRISfree(&mris) ;
#endif

  exit(0) ;
  return(0) ;  /* for ansi */
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
  float f ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "synth"))
  {
    synth_name = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "overlay"))
  {
    mri_overlay = MRIread(argv[2]) ;
    nargs = 1 ;
    if (mri_overlay == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read overlay from %s", argv[2]) ;
    parms.niterations = 0 ;   // this will disable the actual flattening
  }
  else if (!stricmp(option, "label_overlay") || !stricmp(option, "overlay_label"))
  {
    label_overlay = LabelRead(NULL, argv[2]) ;
    nargs = 1 ;
    if (label_overlay == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label overlay from %s", Progname,argv[2]) ;
  }
  else if (!stricmp(option, "norand"))
  {
    setRandomSeed(0L) ;
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number genererator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sphere"))
  {
    sphere_flag = 1 ;
  }
  else if (!stricmp(option, "plane"))
  {
    plane_flag = 1 ;
  }
  else if (!stricmp(option, "rescale"))
  {
    rescale = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "rescaling brain by %2.3f\n", rescale) ;
  }
  else if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dilating cuts %d times\n", dilate) ;
  }
  else if (!stricmp(option, "copy-coords"))
  {
    SurfCopyCoords = MRISread(argv[2]);
    if(SurfCopyCoords == NULL) exit(1);
    nargs = 1 ;
    fprintf(stderr, "copying coords from %s\n", argv[2]);
  }
  else if (!stricmp(option, "dist"))
  {
    sscanf(argv[2], "%f", &parms.l_dist) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.3f\n", parms.l_dist) ;
  }
  else if (!stricmp(option, "expand"))
  {
    sscanf(argv[2], "%f", &parms.l_expand) ;
    nargs = 1 ;
    printf("setting l_expand = %2.3f\n", parms.l_expand) ;
  }
  else if (!stricmp(option, "unfold"))
  {
    MRI *mri_kernel, *mri_tmp ;

    sscanf(argv[2], "%f", &parms.l_unfold) ;
    mri_tmp = MRIread(argv[3]) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read distance map %s...\n",
                Progname, argv[3]) ;

    mri_kernel = MRIgaussian1d(1.0, -1) ;
    parms.mri_dist = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
    MRIfree(&mri_kernel) ;
    MRIfree(&mri_tmp) ;
    nargs = 2 ;
    fprintf(stderr, "l_unfold = %2.3f\n", parms.l_unfold) ;
  }
  else if (!stricmp(option, "boundary"))
  {
    sscanf(argv[2], "%f", &parms.l_boundary) ;
    nargs = 1 ;
    fprintf(stderr, "l_boundary = %2.3f\n", parms.l_boundary) ;
  }
  else if (!stricmp(option, "lm"))
  {
    parms.integration_type = INTEGRATE_LINE_MINIMIZE ;
    fprintf(stderr, "integrating with line minimization\n") ;
  }
  else if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    parms.base_dt = base_dt_scale*parms.dt ;
    nargs = 1 ;
    fprintf(stderr, "momentum with dt = %2.2f\n", parms.dt) ;
  }
  else if (!stricmp(option, "curv"))
  {
    sscanf(argv[2], "%f", &parms.l_curv) ;
    nargs = 1 ;
    fprintf(stderr, "using l_curv = %2.3f\n", parms.l_curv) ;
  }
  else if (!stricmp(option, "nospring"))
  {
    nospring = 1 ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "nlarea"))
  {
    sscanf(argv[2], "%f", &parms.l_nlarea) ;
    nargs = 1 ;
    fprintf(stderr, "using l_nlarea = %2.3f\n", parms.l_nlarea) ;
  }
  else if (!stricmp(option, "boundary"))
  {
    sscanf(argv[2], "%f", &parms.l_boundary) ;
    nargs = 1 ;
    fprintf(stderr, "using l_boundary = %2.3f\n", parms.l_boundary) ;
  }
  else if (!stricmp(option, "adaptive"))
  {
    parms.integration_type = INTEGRATE_ADAPTIVE ;
    fprintf(stderr, "using adaptive time step integration\n") ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "spring"))
  {
    sscanf(argv[2], "%f", &parms.l_spring) ;
    nargs = 1 ;
    fprintf(stderr, "using l_spring = %2.3f\n", parms.l_spring) ;
  }
  else if (!stricmp(option, "name"))
  {
    strcpy(parms.base_name, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using base name = %s\n", parms.base_name) ;
  }
  else if (!stricmp(option, "angle"))
  {
    sscanf(argv[2], "%f", &parms.l_angle) ;
    nargs = 1 ;
    fprintf(stderr, "using l_angle = %2.3f\n", parms.l_angle) ;
  }
  else if (!stricmp(option, "area"))
  {
    sscanf(argv[2], "%f", &parms.l_area) ;
    nargs = 1 ;
    fprintf(stderr, "using l_area = %2.3f\n", parms.l_area) ;
  }
  else if (!stricmp(option, "tol"))
  {
    if (sscanf(argv[2], "%e", &f) < 1)
      ErrorExit(ERROR_BADPARM, "%s: could not scan tol from %s",
                Progname, argv[2]) ;
    parms.tol = (double)f ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.2e\n", (float)parms.tol) ;
  }
  else if (!stricmp(option, "error_ratio"))
  {
    parms.error_ratio = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "error_ratio=%2.3f\n", parms.error_ratio) ;
  }
  else if (!stricmp(option, "dt_inc"))
  {
    parms.dt_increase = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_increase=%2.3f\n", parms.dt_increase) ;
  }
  else if (!stricmp(option, "complete"))
  {
    parms.complete_dist_mat = 1 ;
    fprintf(stderr, "using complete distance matrix\n") ;
  }
  else if (!stricmp(option, "vnum") || (!stricmp(option, "distances")))
  {
    parms.nbhd_size = atof(argv[2]) ;
    parms.max_nbrs = atof(argv[3]) ;
    nargs = 2 ;
    fprintf(stderr, "sampling %d neighbors out to a distance of %d mm\n",
            parms.max_nbrs, parms.nbhd_size) ;
  }
  else if (!stricmp(option, "dt_dec"))
  {
    parms.dt_decrease = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt_decrease=%2.3f\n", parms.dt_decrease) ;
  }
  else if (!stricmp(option, "ou"))
  {
    original_unfold_surf_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr,"reading original unfolding surface from %s...\n",
            original_unfold_surf_name);
  }
  else if (!stricmp(option, "as"))
  {
    parms.area_coef_scale = atof(argv[2]) ;
    printf("setting area coef scale to %2.3f\n",
           parms.area_coef_scale) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'P':
      max_passes = atoi(argv[2]) ;
      fprintf(stderr, "limitting unfolding to %d passes\n", max_passes) ;
      nargs = 1 ;
      break ;
    case '1':
      one_surf_flag = 1 ;  /* patch is only surface file */
      break ;
    case 'O':
      original_surf_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr,"reading original surface from %s...\n",
              original_surf_name);
      break ;
    case 'B':
      base_dt_scale = atof(argv[2]) ;
      parms.base_dt = base_dt_scale*parms.dt ;
      nargs = 1;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'L':
      label_fname = argv[2] ;
      dilate_label = atof(argv[3]) ;
      nargs = 2;
      printf("loading label %s and dilating it %d times before flattening\n",
             label_fname, dilate_label) ;
      break ;
    case 'D':
      disturb = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "disturbing vertex positions by %2.3f\n",
              (float)disturb) ;
      break ;
    case 'R':
      randomly_flatten = !randomly_flatten ;
      fprintf(stderr, "using random placement for flattening.\n") ;
      break ;
    case 'M':
      parms.integration_type = INTEGRATE_MOMENTUM ;
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "momentum = %2.2f\n", (float)parms.momentum) ;
      break ;
    case 'S':
      scale = atof(argv[2]) ;
      fprintf(stderr, "scaling brain by = %2.3f\n", (float)scale) ;
      nargs = 1 ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      sscanf(argv[2], "%d", &parms.write_iterations) ;
      nargs = 1 ;
      fprintf(stderr, "using write iterations = %d\n",
              parms.write_iterations) ;
      break ;
    case 'I':
      do_inflate = 1 ;
      fprintf(stderr, "inflating brain...\n") ;
      break ;
    case 'A':
      sscanf(argv[2], "%d", &parms.n_averages) ;
      nargs = 1 ;
      fprintf(stderr, "using n_averages = %d\n", parms.n_averages) ;
      break ;
    case 'N':
      sscanf(argv[2], "%d", &parms.niterations) ;
      nargs = 1 ;
      fprintf(stderr, "using niterations = %d\n", parms.niterations) ;
      break ;
    default:
    case 'H':
    case '?':
    case 'U':
      print_help() ;
      break ;
    }

  return(nargs) ;
}

static void
print_usage(void)
{
  fprintf(stderr,
          "Usage: %s [options] <input patch> <output patch>\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\nThis program will flatten a surface patch\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, " -w <# iterations>\n\t"
          "write out the surface every # of iterations.\n") ;
  fprintf(stderr,
          " -distances <nbhd size> <# of vertices at each distance>\n\t"
          "specify size of neighborhood and number of vertices at each\n\t"
          "distance to be used in the optimization.\n") ;
  fprintf(stderr,
          " -dilate <# of dilations>\n\t"
          "specify the number of times to dilate the ripped edges to ensure a clean cut\n") ;
  fprintf(stderr,
          " -norand\n\t"
          "set the random seed to 0 so that flattening is repeatable\n") ;
  fprintf(stderr,
          " -seed <random seed>\n\t"
          "set the random seed to a specific value so that flattening is repeatable\n") ;
  printf(" -copy-coords surf : copy xyz coords from surface before flattening\n");
  printf("    This allows creating cuts on, eg, white surface, but flattening the inflated\n");
  printf("    eg, mris_flatten -copy-coords lh.inflated -dilate 1 lh.white.cut lh.inflated.cut.flatten\n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

