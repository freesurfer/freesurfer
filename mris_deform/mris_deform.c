/**
 * @file  mris_deform.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:10 $
 *    $Revision: 1.4 $
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "mrisurf.h"
#include "histo.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) ;
static int compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms,
                                      VERTEX_INFO *vi, float sigma) ;
static int compute_target_labels(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms,
                                 VERTEX_INFO *vi) ;


char *Progname ;

static int min_averages = 0 ;
static int max_averages = 4 ;
static float sigma = 2.0 ;
static int vavgs = 0 ;
static int nbrs = 2 ;

static void usage_exit(int code) ;

static INTEGRATION_PARMS parms ;
static double externalGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) ;

int
main(int argc, char *argv[]) {
  char        *out_fname, **av ;
  int         ac, nargs, i ;
  MRI         *mri_intensities, *mri_labels, *mri_kernel, *mri_smooth=NULL ;
  MRI_SURFACE *mris ;
  int         msec, minutes, seconds, n_averages ;
  float        current_sigma ;
  struct timeb start ;
  char         cmdline[CMD_LINE_LEN], *cp ;
  HISTOGRAM   *histos[MAX_LABEL+1] ;
  VERTEX_INFO *vi ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.integration_type = INTEGRATE_MOMENTUM ;

  // parms.l_nspring = .5; parms.l_tspring = 1; parms.l_curv = 1.0 ;
  parms.l_spring = 1;
  parms.l_curv = 1.0 ;

  // parms.l_intensity = 0.1 ;
  parms.l_repulse = 0 ;
  parms.l_external = 1 ;
  parms.n_averages = 4 ;
  parms.niterations = 1000 ;
  // parms.l_surf_repulse = .1 ;
  parms.dt = parms.base_dt = 0.5 ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_deform.c,v 1.4 2006/12/29 02:09:10 nicks Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_deform.c,v 1.4 2006/12/29 02:09:10 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read input surface %s\n", Progname,argv[1]);
  if (IS_QUADRANGULAR(mris))
    MRISremoveTriangleLinks(mris) ;

  cp = strstr(out_fname, "lh.") ;
  if (cp == NULL)
    cp = strstr(out_fname, "rh.") ;
  if (cp == NULL)
    FileNameExtension(out_fname,parms.base_name) ;  // remove hemi (e.g. lh.)
  else
    strcpy(parms.base_name, cp+3) ;

  mri_labels = MRIread(argv[2]) ;
  if (mri_labels == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume from %s\n", Progname,argv[2]);

  mri_intensities = MRIread(argv[3]) ;
  if (mri_intensities == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read intensity volume from %s\n", Progname,argv[3]);

  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  if (nbrs > 1)
    MRISsetNeighborhoodSize(mris, nbrs) ;
  build_label_histograms(mri_labels, mri_intensities, histos) ;

  vi = (VERTEX_INFO *)calloc(mris->nvertices, sizeof(VERTEX_INFO)) ;
  if (vi == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d vertex info table", Progname, mris->nvertices) ;
  parms.user_parms = (void *)vi ;
  compute_target_labels(mris, mri_labels, mri_intensities, histos, vi) ;
  current_sigma = sigma ;
  gMRISexternalGradient = externalGradient ;
  gMRISexternalSSE = externalSSE ;
  gMRISexternalRMS = externalRMS ;

  for (n_averages = max_averages, i = 0 ;
       n_averages >= min_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    parms.sigma = current_sigma ;
    mri_kernel = MRIgaussian1d(current_sigma, 100) ;
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
            current_sigma) ;
    if (!mri_smooth)
      mri_smooth = MRIclone(mri_intensities, NULL) ;
    MRIconvolveGaussian(mri_intensities, mri_smooth, mri_kernel) ;
    MRIfree(&mri_kernel) ;
    compute_target_intensities(mris, mri_labels, mri_intensities, histos, vi, current_sigma) ;
    if (vavgs > 0)
      MRISaverageVals(mris, vavgs) ;
    parms.n_averages = n_averages ;
    MRISpositionSurface(mris, mri_intensities, mri_intensities, &parms) ;
    if (!n_averages)
      break ;
  }

  MRISaddCommandLine(mris, cmdline) ;
  MRISwrite(mris, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("surface deformation took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "dt")) {}
  else if (!stricmp(option, "tspring")) {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  } else if (!stricmp(option, "vavgs")) {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing values for %d iterations\n", vavgs) ;
  } else if (!stricmp(option, "intensity")) {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.3f\n", parms.l_intensity) ;
  } else if (!stricmp(option, "grad")) {
    parms.l_grad = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_grad = %2.3f\n", parms.l_grad) ;
  } else if (!stricmp(option, "nspring")) {
    parms.l_nspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nspring = %2.3f\n", parms.l_nspring) ;
  } else if (!stricmp(option, "curv")) {
    parms.l_curv = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_curv = %2.3f\n", parms.l_curv) ;
  } else switch (toupper(*option)) {
    case 'A':
      max_averages = atoi(argv[2]) ;
      printf("using %d smoothing steps as max surface smoothing\n", max_averages) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
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
      printf("writing snapshots of deformation every %d iterations\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case '?':
    case 'U':
      usage_exit(0) ;
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
usage_exit(int code) {
  printf("usage: %s [options] <input surface> <label vol> <intensity vol> <output surface>\n", Progname) ;
  exit(code) ;
}

static int
build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) {
  int labels[MAX_LABEL+1], x, y, z, l ;

  memset(labels, 0, sizeof(labels)) ;


  for (x = 0 ; x < mri_labels->width ; x++)
    for (y = 0 ; y < mri_labels->height ; y++)
      for (z = 0 ; z < mri_labels->depth ; z++) {
        l = nint(MRIgetVoxVal(mri_labels, x, y, z, 0)) ;
        if (l == 0)
          continue ;
        if (labels[l] == 0)  // first time
        {
          char fname[STRLEN] ;
          histos[l] = MRIhistogramLabel(mri_intensities, mri_labels, l, 50) ;
          HISTOmakePDF(histos[l], histos[l]) ;
          sprintf(fname, "label%d.plt", l) ;
          HISTOplot(histos[l], fname) ;
        }

        labels[l] = 1 ;
      }

  return(NO_ERROR) ;
}

#define PTHRESH  0.005

static int
compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms, VERTEX_INFO *vi,float sigma) {
  VERTEX   *v ;
  int      vno, bin_no ;
  Real     xv, yv, zv, val, mag, mag0, nx, ny, nz, xk, yk, zk, xi, yi, zi, max_mag ;
  HISTOGRAM *h_in, *h_out ;

  sigma = 0 ;
  // first sample initiali directional derivatives in normal direction and filter them to get consistent sign
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;


    /* search locally for a voxel that has large intensity gradients, pointing in the right
      general direction with intensities that could be in the right class
    */
    MRISvertexToVoxel(mris, v, mri_intensities, &xv, &yv, &zv) ;
    MRISvertexNormalInVoxelCoords(mris, mri_labels, vno, &nx, &ny, &nz) ;
    MRIsampleVolumeDerivativeScale(mri_intensities, xv, yv, zv, nx, ny, nz, &mag0,sigma) ;
    vi[vno].mag = mag0 ;

    vi[vno].tx = v->x ;
    vi[vno].ty = v->y ;
    vi[vno].tz = v->z ;
    v->val = mag0 ;
  }
#if 0
  MRISmedianFilterVals(mris, 2) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    vi[vno].mag = v->val ;
  }
#endif
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, mri_intensities, &xv, &yv, &zv) ;

    /* search locally for a voxel that has large intensity gradients, pointing in the right
      general direction with intensities that could be in the right class
    */
    v->val = val ;
    v->val2 = 2 ;
    MRISvertexNormalInVoxelCoords(mris, mri_labels, vno, &nx, &ny, &nz) ;
    mag0 = vi[vno].mag ;
    max_mag = fabs(mag0) ;

    h_out = vi[vno].h_out ;
    h_in = vi[vno].h_in ;
    if (h_out == NULL || h_in == NULL)
      continue ;

    for (xk = -1 ; xk <= 1 ; xk += 0.25)
      for (yk = -1 ; yk <= 1 ; yk += 0.25)
        for (zk = -1 ; zk <= 1 ; zk += 0.25) {
          xi = xv+xk ;
          yi = yv+yk ;
          zi = zv+zk ;

          // check 'outside' value
          MRIsampleVolume(mri_intensities, xi+0.5*nx, yi+0.5*ny, zi+0.5*nz, &val);
          bin_no = nint((float)(val - h_out->min) / (float)h_out->bin_size) ;
          if (bin_no < 0 || bin_no >= h_out->nbins || h_out->counts[bin_no] < PTHRESH)
            continue ;  // can't be here

          // check 'inside' value
          MRIsampleVolume(mri_intensities, xi-0.5*nx, yi-0.5*ny, zi-0.5*nz, &val);
          bin_no = nint((float)(val - h_in->min) / (float)h_in->bin_size) ;
          if (bin_no < 0 || bin_no >= h_in->nbins || h_in->counts[bin_no] < PTHRESH)
            continue ;  // can't be here
          MRIsampleVolumeDerivativeScale(mri_intensities, xi, yi, zi, nx, ny, nz, &mag, sigma) ;
          if (fabs(mag) > max_mag && (mag*mag0 > 0)) {
            double xw, yw, zw ;

            MRIsampleVolume(mri_intensities, xi, yi, zi, &val);
            v->val = val ;
            max_mag = fabs(mag) ;
            if (mris->useRealRAS)
              MRIvoxelToWorld(mri_intensities, xi, yi, zi, &xw, &yw, &zw) ;
            else
              MRIvoxelToSurfaceRAS(mri_intensities, xi, yi, zi, &xw, &yw, &zw) ;
            vi[vno].tx = xw ;
            vi[vno].ty = yw ;
            vi[vno].tz = zw ;
          }
        }
    if (vno == Gdiag_no)
      printf("v %d: target (%2.1f, %2.1f, %2.1f)\n", vno, vi[vno].tx, vi[vno].ty, vi[vno].tz) ;
    v->marked = 1 ;
  }

  return(NO_ERROR) ;
}

static int
compute_target_labels(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms, VERTEX_INFO *vi) {
  VERTEX   *v ;
  int      vno, l_out, l_in ;
  Real     xv, yv, zv, val, nx, ny, nz ;
  HISTOGRAM *h_in, *h_out ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, mri_intensities, &xv, &yv, &zv) ;
    MRIsampleVolume(mri_intensities, xv, yv, zv, &val);
    v->val = val ;
    v->val2 = 2 ;

    /* search locally for a voxel that has large intensity gradients, pointing in the right
      general direction with intensities that could be in the right class
    */
    MRISvertexNormalInVoxelCoords(mris, mri_labels, vno, &nx, &ny, &nz) ;
    MRIsampleVolumeType(mri_labels, xv+1*nx, yv+1*ny, zv+1*nz, &val, SAMPLE_NEAREST);
    l_out = nint(val) ;
    MRIsampleVolumeType(mri_labels, xv-1*nx, yv-1*ny, zv-1*nz, &val, SAMPLE_NEAREST);
    l_in = nint(val) ;
    h_out = histograms[l_out] ;
    h_in = histograms[l_in] ;
    vi[vno].h_in = h_in ;
    vi[vno].h_out = h_out ;
    vi[vno].l_in = l_in ;
    vi[vno].l_out = l_out ;
  }

  return(NO_ERROR) ;
}

static double
externalGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  int         vno ;
  VERTEX      *v ;
  double      dx, dy, dz, sse, dot ;
  VERTEX_INFO *vi ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  vi = (VERTEX_INFO *)parms->user_parms ;
  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    dx = vi[vno].tx - v->x ;
    dy = vi[vno].ty - v->y ;
    dz = vi[vno].tz - v->z ;

#if 1
    dot = (dx * v->nx + dy * v->ny + dz * v->nz) ;
    dx = v->nx*dot ;
    dy = v->ny*dot ;
    dz = v->nz*dot ;
#endif
    if (vno == Gdiag_no) {
      printf("l_extern: v %d, target (%2.1f, %2.1f, %2.1f), current (%2.1f, %2.1f, %2.1f), dX=(%2.1f, %2.1f, %2.1f)\n",
             vno, vi[vno].tx, vi[vno].ty, vi[vno].tz, v->x, v->y, v->z, dx, dy, dz) ;
    }
    v->dx += parms->l_external*dx ;
    v->dy += parms->l_external*dy ;
    v->dz += parms->l_external*dz ;
    sse += dx*dx + dy*dy * dz+dz ;
  }

  return(sse * parms->l_external) ;
}
static double
externalSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) {
  int         vno ;
  VERTEX      *v ;
  double      dx, dy, dz, sse ;
  VERTEX_INFO *vi ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  vi = (VERTEX_INFO *)parms->user_parms ;
  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    dx = vi[vno].tx - v->x ;
    dy = vi[vno].ty - v->y ;
    dz = vi[vno].tz - v->z ;

    if (vno == Gdiag_no) {
      printf("l_extern: v %d, target (%2.1f, %2.1f, %2.1f), current (%2.1f, %2.1f, %2.1f), dX=(%2.1f, %2.1f, %2.1f)\n",
             vno, vi[vno].tx, vi[vno].ty, vi[vno].tz, v->x, v->y, v->z, dx, dy, dz) ;
    }
    sse += dx*dx + dy*dy * dz+dz ;
  }

  return(sse * parms->l_external) ;
}

static double
externalRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) {
  return(sqrt(externalSSE(mris, parms)/mris->nvertices)) ;
}

