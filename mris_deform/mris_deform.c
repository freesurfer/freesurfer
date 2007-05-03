/**
 * @file  mris_deform.c
 * @brief generic surface deformation binary.
 *
 * binary for deforming a surface of a voxel segmentation to more smoothly and
 * accurately represent the border.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/05/03 11:52:43 $
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
#include "transform.h"
#include "gca.h"
#include "gcamorph.h"
#include "tritri.h"  // for CROSS3 definition

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) ;
static MRI *compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms,
                                      VERTEX_INFO *vi, float sigma,
                                      TRANSFORM *transform, GCA *gca, int label) ;
static int compute_target_labels(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms,
                                 VERTEX_INFO *vi) ;


char *Progname ;

static char *gca_fname = NULL ;
static int min_averages = 0 ;
static int max_averages = 4 ;
static float sigma = 1.0 ;
static int vavgs = 0 ;
static int nbrs = 2 ;
static int target_label = -1 ;
static int renormalize_gca = 1 ;
static int read_ll = 0 ;

static void usage_exit(int code) ;

static char *gca_write_fname = NULL ;
static INTEGRATION_PARMS parms ;
static double externalGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) ;


int
main(int argc, char *argv[]) {
  char        *out_fname, **av ;
  int         ac, nargs, i ;
  MRI         *mri_intensities, *mri_labels/*, *mri_kernel, *mri_smooth=NULL*/, *mri_ll ;
  MRI_SURFACE *mris ;
  int         msec, minutes, seconds, n_averages ;
  float        current_sigma ;
  struct timeb start ;
  char         cmdline[CMD_LINE_LEN], *cp ;
  HISTOGRAM   *histos[MAX_LABEL+1] ;
  VERTEX_INFO *vi ;
  TRANSFORM   *transform ;
  GCA         *gca ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.integration_type = INTEGRATE_MOMENTUM ;

  // parms.l_nspring = .5; parms.l_tspring = 1; parms.l_curv = 1.0 ;
  parms.l_spring = .1;
  parms.l_nlspring = 1 ;
  parms.rmin = 0.5 ;
  parms.rmax = 5 ;
  parms.l_curv = .0 ;

  // parms.l_intensity = 0.1 ;
  parms.l_repulse = 0 ;
  parms.check_tol = 1 ;  // don't use intensity rms in surface deformation, use sse
  parms.tol = 0.01 ;
  parms.l_external = 1 ;
  parms.n_averages = 4 ;
  parms.niterations = 1000 ;
  // parms.l_surf_repulse = .1 ;
  parms.dt = parms.base_dt = 0.5 ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_deform.c,v 1.5 2007/05/03 11:52:43 fischl Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_deform.c,v 1.5 2007/05/03 11:52:43 fischl Exp $", "$Name:  $");
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

  if (argc < 6)
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

  transform = TransformRead(argv[3]) ;
  if (transform == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read transform from %s\n", Progname,argv[3]);

  if (gca_fname)
    cp = gca_fname ;
  else if (transform->type == MORPH_3D_TYPE)
    cp = ((GCA_MORPH *)transform->xform)->atlas.fname ;
  else
    ErrorExit(ERROR_BADPARM, "%s: must specify GCA with -gca", Progname) ;

  printf("reading atlas from %s\n", cp) ;
  gca = GCAread(cp) ;
  if (gca == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read atlas from %s\n", 
              Progname,cp);
  GCAregularizeConditionalDensities(gca, 0.5) ;

  mri_intensities = MRIread(argv[4]) ;
  if (mri_intensities == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read intensity volume from %s\n", Progname,argv[4]);

  TransformInvert(transform, mri_intensities) ;
  if (renormalize_gca)
  {
    GCAmapRenormalizeWithAlignment(gca, mri_intensities, transform,
                                   NULL, "", NULL, 0) ;
    if (gca_write_fname)
      GCAwrite(gca, gca_write_fname) ;
  }
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
  MRISaverageVertexPositions(mris, 10) ;

  for (n_averages = max_averages, i = 0 ;
       n_averages >= min_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    parms.sigma = current_sigma ;
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
            current_sigma) ;
#if 0
    mri_kernel = MRIgaussian1d(current_sigma, 100) ;
    if (!mri_smooth)
      mri_smooth = MRIclone(mri_intensities, NULL) ;
    MRIconvolveGaussian(mri_intensities, mri_smooth, mri_kernel) ;
    MRIfree(&mri_kernel) ;
#endif

    if (i == 0)
    {
      if (read_ll)
        mri_ll = MRIread("ll.mgz") ;
      else
        mri_ll = compute_target_intensities(mris, mri_labels, mri_intensities, histos, vi, current_sigma, 
                                        transform, gca, target_label) ;
      parms.mri_ll = mri_ll ;
    }
    if (vavgs > 0)
      MRISaverageVals(mris, vavgs) ;
    {
      char valname[STRLEN] ;
      strcpy(valname, out_fname) ;
      strcat(valname, ".targets") ;
      printf("writing target intensities to %s\n", valname) ;
      MRIScopyValuesToCurvature(mris) ;
      MRISwriteCurvature(mris, valname) ;
    }
    parms.n_averages = n_averages ;
    MRISpositionSurface(mris, mri_intensities, mri_intensities, &parms) ;
    //    parms.l_external /= 2;
    parms.rmin = 2 ;
    parms.rmax = 10 ;
    MRISpositionSurface(mris, mri_intensities, mri_intensities, &parms) ;
    //    parms.l_external *= 2 ;
    if (!n_averages)
      break ;
  }

  MRISaddCommandLine(mris, cmdline) ;
  printf("writing output surface to %s\n", out_fname) ;
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
  } else if (!stricmp(option, "nlspring")) {
    parms.l_nlspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nlspring = %2.3f\n", parms.l_nlspring) ;
  } else if (!stricmp(option, "vavgs")) {
    vavgs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "smoothing values for %d iterations\n", vavgs) ;
  } else if (!stricmp(option, "DEBUG_VOXEL")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx,Gy,Gz) ;
  } else if (!stricmp(option, "tol")) {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %2.3e\n", parms.tol) ;
  } else if (!stricmp(option, "write_gca")) {
    gca_write_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "writing renormalized GCA to %s\n", gca_fname) ;
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
  case 'I':
    parms.l_external = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting intensity term to %2.3f\n", parms.l_external) ;
    break ;
  case 'G':
    gca_fname = argv[2] ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    break ;
  case 'R':
    read_ll = 1 ;
    printf("reading ll volume from ll.mgz\n") ;
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
  case 'L':
    target_label = atoi(argv[2]) ;
    printf("using %s (%d) as target label\n", cma_label_to_name(target_label), target_label) ;
    nargs = 1 ;
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
build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) 
{
  int labels[MAX_LABEL+1], x, y, z, l ;

  memset(labels, 0, sizeof(labels)) ;

  for (x = 0 ; x < mri_labels->width ; x++)
    for (y = 0 ; y < mri_labels->height ; y++)
      for (z = 0 ; z < mri_labels->depth ; z++) {
        l = nint(MRIgetVoxVal(mri_labels, x, y, z, 0)) ;
        if (l == 0 || l != 17)
          continue ;
        if (labels[l] == 0)  // first time
        {
          HISTO *h ;
          char fname[STRLEN] ;

          h = MRIhistogramLabel(mri_intensities, mri_labels, l, 50) ;
          histos[l] = HISTOsmooth(h, NULL, 2) ;
          HISTOmakePDF(histos[l], histos[l]) ;
          sprintf(fname, "label%d.plt", l) ;
          HISTOplot(histos[l], fname) ;
        }

        labels[l] = 1 ;
      }

  return(NO_ERROR) ;
}

#define PTHRESH  0.005

static MRI *
compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms, VERTEX_INFO *vi,float sigma, 
                           TRANSFORM *transform, GCA *gca, int target_label) {
  int      nin, nout ;
  Real     val, nx, ny, nz, xi, yi, zi, d, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2, mag ;
  MRI      *mri_ll ;
  MRI      *mri_grad, *mri_dx, *mri_dy, *mri_dz, *mri_mask, 
           *mri_tmp, *mri_dist, *mri_pin, *mri_pout ;
  int       x, y, z, pad ;
  MATRIX    *m_vox2vox, *m_vox2vox_vector ;
  VECTOR    *v1, *v2 ;

#define WSIZE 3

  mri_tmp = MRISfillInterior(mris, mri_labels->xsize/8, NULL) ;
  pad = ceil(3/mri_tmp->xsize) ;
  mri_mask = MRIextractRegionAndPad(mri_tmp, NULL, NULL, pad) ;
  mri_dist = MRIdistanceTransform(mri_mask, NULL, 1, 10, DTRANS_MODE_SIGNED) ;
  MRIwrite(mri_dist, "d.mgz") ;

  mri_grad = MRIsobel(mri_dist, NULL, NULL) ;
  mri_dx = MRIxDerivative(mri_dist, NULL) ;
  mri_dy = MRIyDerivative(mri_dist, NULL) ;
  mri_dz = MRIzDerivative(mri_dist, NULL) ;
  MRIfree(&mri_grad) ; MRIfree(&mri_tmp) ; 
  mri_ll = MRIclone(mri_dist, NULL) ;
  mri_pin = MRIclone(mri_ll, NULL) ;
  mri_pout = MRIclone(mri_ll, NULL) ;

  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_dist, "d.mgz") ; MRIwrite(mri_mask, "m.mgz") ;
    MRIwrite(mri_dx, "dx.mgz") ; MRIwrite(mri_dy, "dy.mgz") ; 
    MRIwrite(mri_dz, "dz.mgz") ; 
  }

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_mask, mri_intensities) ;
  m_vox2vox_vector = MatrixCopy(m_vox2vox, NULL) ;
  *MATRIX_RELT(m_vox2vox_vector, 1, 4) = 0.0 ;
  *MATRIX_RELT(m_vox2vox_vector, 2, 4) = 0.0 ;
  *MATRIX_RELT(m_vox2vox_vector, 3, 4) = 0.0 ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;

  for (x = 0 ; x < mri_ll->width; x++)
  {
    for (y = 0 ; y < mri_ll->height ; y++)
    {
      for (z = 0 ; z < mri_ll->depth ; z++)
      {
        double xw, yw, zw, pout, pin, p ;

        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
        if (val > pad)
          continue ;
        nx = MRIgetVoxVal(mri_dx, x, y, z, 0) ;
        ny = MRIgetVoxVal(mri_dy, x, y, z, 0) ;
        nz = MRIgetVoxVal(mri_dz, x, y, z, 0) ;
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
          continue ;
        V3_X(v1) = nx ; V3_Y(v1) = ny ; V3_Z(v1) = nz ;
        MatrixMultiply(m_vox2vox_vector, v1, v2) ; // to aseg/norm coords
        nx = V3_X(v2) ; ny = V3_Y(v2) ; nz = V3_Z(v2) ; 
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        nx /= mag ; ny /= mag ; nz /= mag ;

        // build two tangent vectors
        e1x = -ny ; e1y = -nz ; e1z = nx ;  // non-parallel vector
        CROSS3(e2x, e2y, e2z, nx, ny, nz, e1x, e1y, e1z) ;
        CROSS3(e1x, e1y, e1z, nx, ny, nz, e2x, e2y, e2z) ;
        mag = sqrt(e1x*e1x+e1y*e1y+e1z*e1z) ;
        if (FZERO(mag))
          mag = 1e-5 ;
        e1x /= mag ; e1y /= mag ; e1z /= mag ;
        mag = sqrt(e2x*e2x+e2y*e2y+e2z*e2z) ;
        if (FZERO(mag))
          mag = 1e-5 ;
        e2x /= mag ; e2y /= mag ; e2z /= mag ;

        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xi = V3_X(v2) ; yi = V3_Y(v2) ; zi = V3_Z(v2) ; 

        for (nin=0, pin = 0.0, d = -WSIZE ; d < 0 ; d += 0.5)
        {
          if (DZERO(d))
            continue ;
          xw = xi+d*nx ; yw = yi+d*ny ; zw = zi+d*nz ;
          MRIsampleVolume(mri_intensities, xw, yw, zw, &val);  // for debugging
          p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
          
          if (DZERO(p))
            p = 1e-5;
          nin++ ;
          pin += log(p) ;
        }
        
        d = 0.5 ;  // distance to sample outwards - d1 and d2 will be in tangent directions
        for (nout=0, pout = 0.0, d1 = -1 ; d1 <= 1 ; d1 += 0.5)
          for (d2 = -1 ; d2 <= 1 ; d2 += 0.5)
          {
            xw = xi+ d*nx + d1*e1x + d2*e2x ; yw = yi + d*ny + d1*e1y + d2*e2y ; zw = zi + d*nz + d1*e1z + d2*e1z ;
            MRIsampleVolume(mri_intensities, xw, yw, zw, &val);
            p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
            
            if (p < .25)  // don't care how unlikely it is if its past a certain point
              p = .25;
            pout += log(p) ;
            nout++ ;
            
          }
        pin /= nin ; pout /= nout ;
        if (FZERO(pout))
          pout = 1e-5;
        MRIsetVoxVal(mri_ll, x, y, z, 0, 2*pin-pout) ;
        MRIsetVoxVal(mri_pin, x, y, z, 0, pin) ;
        MRIsetVoxVal(mri_pout, x, y, z, 0, pout) ;
      }
    }
  }
  MRIwrite(mri_ll, "ll.mgz") ;
  MRIwrite(mri_pin, "pin.mgz") ;
  MRIwrite(mri_pout, "pout.mgz") ;
  MRIfree(&mri_mask) ; MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ; MRIfree(&mri_dist) ;
  MatrixFree(&m_vox2vox_vector) ; MRIfree(&mri_pin) ; MRIfree(&mri_pout) ;

  return(mri_ll) ;
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

    if (target_label > 0)
    {
      vi[vno].h_in = histograms[target_label] ;
      vi[vno].l_in = target_label ;
    }
  }

  return(NO_ERROR) ;
}

#if 1
#define MAX_LL 10

static double
externalSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int         vno ;
  VERTEX      *v ;
  double      sse ;
  Real        val, xv, yv, zv ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, parms->mri_ll, &xv, &yv, &zv) ;
    MRIsampleVolume(parms->mri_ll, xv, yv, zv, &val);
    if (vno == Gdiag_no) 
      printf("l_extern: v %d, ll vox (%2.1f, %2.1f, %2.1f) = %2.2f\n", 
             vno, xv, yv, zv, val) ;

    if (val > MAX_LL)
      val = MAX_LL ;
    sse += (MAX_LL-val)*(MAX_LL-val) ;
  }

  return(sse * parms->l_external) ;
}

static double
externalGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int         vno ;
  VERTEX      *v ;
  double      sse, mag, error ;
  Real        val, xv, yv, zv, nx, ny, nz ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, parms->mri_ll, &xv, &yv, &zv) ;
    MRIsampleVolume(parms->mri_ll, xv, yv, zv, &val);
    if (vno == Gdiag_no) 
      printf("l_extern: v %d, ll vox (%2.1f, %2.1f, %2.1f) = %2.2f\n", 
             vno, xv, yv, zv, val) ;

    if (val > MAX_LL)
      val = MAX_LL ;
    error = MAX_LL-val ;
    sse += error*error ;

    MRISvertexNormalInVoxelCoords(mris, parms->mri_ll, vno, &nx, &ny, &nz) ;
    MRIsampleVolumeDerivativeScale(parms->mri_ll, xv, yv, zv, nx, ny, nz, &mag,v->val2) ;

    if (vno == Gdiag_no)
      printf("l_extern: v %d, voxel (%2.1f, %2.1f, %2.1f), N=(%2.2f, %2.2f, %2.2f), mag = %2.2f\n",
             vno, xv, yv, zv, nx, ny, nz, mag) ;

    if (FZERO(mag))
      mag = 1 ;
    mag /= fabs(mag) ; mag *= parms->l_external ;
    if (vno == Gdiag_no)
      printf("l_extern:   dx = (%2.2f, %2.2f, %2.2f)\n", mag*nx, mag*ny, mag*nz) ;

    if (error > v->d)
      DiagBreak() ;
    v->d = error ;
    v->dx += mag*nx ; 
    v->dy += mag*ny ; 
    v->dz += mag*nz ; 
    if (!devFinite(v->dx) || !devFinite(v->dy) || !devFinite(v->dz))
      DiagBreak() ;
  }

  return(sse * parms->l_external) ;
}

#else
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
    sse += dx*dx + dy*dy + dz*dz ;
  }

  return(sse * parms->l_external) ;
}
static double
externalSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int         vno ;
  VERTEX      *v ;
  double      dx, dy, dz, sse, dist ;
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

    dist = dx*dx + dy*dy + dz*dz;
    
    if (dist > v->d)
      DiagBreak() ;
    v->d = dist ;
    sse += dist;
  }

  return(sse * parms->l_external) ;
}

#endif

static double
externalRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) {
  return(sqrt(externalSSE(mris, parms)/mris->nvertices)) ;
}

