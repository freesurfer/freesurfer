/**
 * @brief generic surface deformation binary.
 *
 * binary for deforming a surface of a voxel segmentation to more smoothly and
 * accurately represent the border.
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
#include "gcaboundary.h"

static int compute_gradient_target_positions(MRI_SURFACE *mris, MRI *mri_intensities, VERTEX_INFO *vi, 
                                             float current_sigma);
static GCA *make_gca(char *label_vol_name, char *intensity_vol_name,
										 TRANSFORM *transform) ;
static MRI *compute_pmap(GCA *gca, TRANSFORM *transform, MRI *mri_intensities, int target_label, MRI *mri_dist, double pad, MRI *mri_pmap) ;
static MRI *compute_pmap_with_gcab(GCA *gca, TRANSFORM *transform, MRI *mri_intensities, MRI *mri_labels, int target_label, MRI *mri_dist, double pad, GCAB *gcab, MRI *mri_pmap) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int build_label_histograms(MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histos) ;
static MRI *compute_target_intensities_with_gcab(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, 
                                                 HISTOGRAM **histograms, VERTEX_INFO *vi,float sigma, 
                                                 TRANSFORM *transform, GCA *gca, int target_label, float resolution,
                                                 GCAB *gcab) ;
static MRI *compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, 
                                       HISTOGRAM **histograms,
                                      VERTEX_INFO *vi, float sigma,
                                      TRANSFORM *transform, GCA *gca, int label, float resolution) ;
static int compute_target_labels(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms,
                                 VERTEX_INFO *vi) ;


const char *Progname ;

static float resolution = 8.0 ;
static char *label_vol_name = NULL ;
static char *intensity_vol_name = NULL ; // for building GCA from hires data
static char *gca_fname = NULL ;
static int min_averages = 0 ;
static int max_averages = 4 ;
static float sigma = 1.0 ;
static int vavgs = 0 ;
static int nbrs = 2 ;
static int target_label = -1 ;
static int renormalize_gca = 1 ;
static char *renorm_seg_fname = NULL ;
static int read_ll = 0 ;

static int use_grad = 0 ;
static double max_grad_dist= 1.5 ; // mm away from current boundary position

static void usage_exit(int code) ;

static char *gca_write_fname = NULL ;
static INTEGRATION_PARMS parms ;

static double externalLLGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalLLSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalLLRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) ;

static double externalGradGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalGradSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) ;
static double externalGradRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) ;

static GCAB *gcab ;

int
main(int argc, char *argv[]) {
  char        *out_fname, **av ;
  int         ac, nargs, i ;
  MRI         *mri_intensities, *mri_labels/*, *mri_kernel, *mri_smooth=NULL*/,
              *mri_ll = NULL ;
  MRI_SURFACE *mris ;
  int         msec, minutes, seconds, n_averages ;
  float        current_sigma ;
  Timer start ;
  char        *cp ;
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

  std::string cmdline = getAllInfo(argc, argv, "mris_ca_deform");

  nargs = handleVersionOption(argc, argv, "mris_ca_deform");
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

  if (argc < 6)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;
  printf("reading input surface from %s\n", argv[1]) ;
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

  mri_intensities = MRIread(argv[4]) ;
  if (mri_intensities == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not read intensity volume from %s\n", Progname,argv[4]);

  if (label_vol_name)   // build gca from hires label/intensity vols
	gca = make_gca(label_vol_name, intensity_vol_name, transform) ;
  else
	{
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
	  if (renormalize_gca == 2)
		{
		  MRI *mri_seg = MRIread(renorm_seg_fname) ;
		  if (mri_seg == NULL)
			exit(Gerror) ;
		  GCArenormalizeToExample(gca, mri_seg, mri_intensities) ;
		  MRIfree(&mri_seg) ;
		}
	  else 
		{
		  GCAmapRenormalizeWithAlignment(gca, mri_intensities, transform, NULL, "", NULL, 0) ;
		  if (gca_write_fname)
			GCAwrite(gca, gca_write_fname) ;
		}
	}
  GCAregularizeConditionalDensities(gca, 0.5) ;
  GCAregularizeCovariance(gca, 0.8) ;
  {
	MRI *mri_tmp = GCAbuildMostLikelyVolume(gca, NULL) ;
	TransformRas2Vox(transform, mri_intensities, mri_tmp) ;
	MRIfree(&mri_tmp) ;
  }
  TransformInvert(transform, mri_intensities) ;
  
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  if (nbrs > 1)
    MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
  build_label_histograms(mri_labels, mri_intensities, histos) ;
  
  vi = (VERTEX_INFO *)calloc(mris->nvertices, sizeof(VERTEX_INFO)) ;
  if (vi == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d vertex info table", Progname, mris->nvertices) ;
  parms.user_parms = (void *)vi ;
  compute_target_labels(mris, mri_labels, mri_intensities, histos, vi) ;
  current_sigma = sigma ;
  if (use_grad)
  {
    gMRISexternalGradient = externalGradGradient ;
    gMRISexternalSSE = externalGradSSE ;
    gMRISexternalRMS = externalGradRMS ;
  }
  else
  {
    gMRISexternalGradient = externalLLGradient ;
    gMRISexternalSSE = externalLLSSE ;
    gMRISexternalRMS = externalLLRMS ;
  }

  MRISaverageVertexPositions(mris, 1) ;

  for (n_averages = max_averages, i = 0 ;
       n_averages >= min_averages ;
       n_averages /= 2, current_sigma /= 2, i++) {
    MRISsetVal2(mris, current_sigma) ;
    parms.sigma = current_sigma ;
#if 0
    fprintf(stderr, "smoothing T1 volume with sigma = %2.3f\n",
            current_sigma) ;
    mri_kernel = MRIgaussian1d(current_sigma, 100) ;
    if (!mri_smooth)
      mri_smooth = MRIclone(mri_intensities, NULL) ;
    MRIconvolveGaussian(mri_intensities, mri_smooth, mri_kernel) ;
    MRIfree(&mri_kernel) ;
#endif

    //    if (i == 0)
    {
      if (read_ll)
        mri_ll = MRIread("ll.mgz") ;
      else if (gcab != NULL)
        mri_ll = compute_target_intensities_with_gcab(mris, mri_labels, mri_intensities, histos, vi, current_sigma, 
                                                      transform, gca, target_label,resolution, gcab) ;
      else if (use_grad)
        compute_gradient_target_positions(mris, mri_intensities, vi, current_sigma);
      else
        mri_ll = compute_target_intensities(mris, mri_labels, mri_intensities, histos, vi, current_sigma, 
                                        transform, gca, target_label,resolution) ;
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
#if 0
    //    parms.l_external /= 2;
    parms.rmin = 2 ;
    parms.rmax = 10 ;
    parms.l_nlspring *= 10 ;
    MRISpositionSurface(mris, mri_intensities, mri_intensities, &parms) ;
    parms.l_nlspring /= 10 ;
#endif
    if (!n_averages)
      break ;
  }

  MRISaddCommandLine(mris, cmdline) ;
  printf("writing output surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = start.milliseconds() ;
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
  if (!stricmp(option, "dt")) {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting time step to %2.3f\n", parms.dt) ;
  }
  else if (!stricmp(option, "tspring")) {
    parms.l_tspring = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_tspring = %2.3f\n", parms.l_tspring) ;
  } else if (!stricmp(option, "resolution")) {
    resolution = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting resolution to 1/%2.1f of the voxel size\n", resolution) ;
  } else if (!stricmp(option, "sigma")) {
    sigma = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting max sigma to %2.3f\n", sigma) ;
  } else if (!stricmp(option, "grad")) {
    use_grad = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "%susing only intensity gradient information\n", use_grad ? "" : "not ") ;
  } else if (!stricmp(option, "grad_dist")) {
    max_grad_dist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting max grad dist to %2.2f mm\n", max_grad_dist) ;
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
    fprintf(stderr, "writing renormalized GCA to %s\n", gca_write_fname) ;
  } else if (!stricmp(option, "rmax")) {
    parms.rmax = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting rmax to %2.2f\n", parms.rmax) ;
  } else if (!stricmp(option, "gcab")) {
    gcab = GCABread(argv[2], NULL) ;
    if (gcab == NULL)
      ErrorExit(ERROR_BADFILE, "%s: could not read boundary atlas from %s", argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using boundary atlas in %s\n", argv[2]) ;
  } else if (!stricmp(option, "rmin")) {
    parms.rmin = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting rmin to %2.2f\n", parms.rmin) ;
  } else if (!stricmp(option, "renorm")) {
    renormalize_gca = atoi(argv[2]) ;
    if (renormalize_gca > 1)
    {
      nargs = 2 ;
      renorm_seg_fname = argv[3] ;
      fprintf(stderr, "%srenormalizing GCA using segmentation %s\n", 
              renormalize_gca ? "" : "not ", renorm_seg_fname) ;
    }
    else
    {
      nargs = 1 ;
      fprintf(stderr, "%srenormalizing GCA\n", renormalize_gca ? "" : "not ") ;
    }
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
  } else if (!stricmp(option, "make_gca")) {
		label_vol_name = argv[2] ;
		intensity_vol_name = argv[3] ;
    nargs = 2 ;
    fprintf(stderr, "building gca from label vol %s and intensity vol %s\n",
						label_vol_name, intensity_vol_name) ;
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
  printf("usage: %s [options] <input surface> <label vol> <transform> <intensity vol> <output surface>\n", Progname) ;
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

          h = MRIhistogramLabel(mri_intensities, mri_labels, l, 50) ;
          histos[l] = HISTOsmooth(h, NULL, 2) ;
          HISTOmakePDF(histos[l], histos[l]) ;
          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
          {
            char fname[STRLEN] ;
            sprintf(fname, "label%d.plt", l) ;
            HISTOplot(histos[l], fname) ;
          }
        }

        labels[l] = 1 ;
      }

  return(NO_ERROR) ;
}

#define PTHRESH  0.005

static MRI *
compute_target_intensities(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, 
                           HISTOGRAM **histograms, VERTEX_INFO *vi,float sigma, 
                           TRANSFORM *transform, GCA *gca, int target_label, float resolution) {
  int      nin, nout ;
  double   val, nx, ny, nz, xi, yi, zi, d, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2, mag ;
  MRI      *mri_ll ;
  MRI      *mri_grad, *mri_dx, *mri_dy, *mri_dz, *mri_mask, 
           *mri_tmp, *mri_dist, *mri_pin, *mri_pout, *mri_pmap ;
  int       x, y, z, pad ;
  MATRIX    *m_vox2vox, *m_vox2vox_vector ;
  VECTOR    *v1, *v2 ;
  double    vsize ;

#define WSIZE 3

  mri_tmp = MRISfillInterior(mris, mri_labels->xsize/resolution, NULL) ;
  mri_tmp->c_r += mri_labels->c_r ;
  mri_tmp->c_a += mri_labels->c_a ;
  mri_tmp->c_s += mri_labels->c_s ;
  pad = ceil(3/mri_tmp->xsize) ;
  mri_mask = MRIextractRegionAndPad(mri_tmp, NULL, NULL, pad) ;
  mri_dist = MRIdistanceTransform(mri_mask, NULL, 1, nint(10/mri_mask->xsize), 
                                  DTRANS_MODE_SIGNED, NULL) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_dist, "d.mgz") ;

  mri_grad = MRIsobel(mri_dist, NULL, NULL) ;
  mri_dx = MRIxDerivative(mri_dist, NULL) ;
  mri_dy = MRIyDerivative(mri_dist, NULL) ;
  mri_dz = MRIzDerivative(mri_dist, NULL) ;
  MRIfree(&mri_grad) ; MRIfree(&mri_tmp) ; 
  mri_ll = MRIclone(mri_dist, NULL) ;
  mri_pin = MRIclone(mri_ll, NULL) ;
  mri_pout = MRIclone(mri_ll, NULL) ;

  mri_ll->outside_val = -100 ;
  pad *= mri_dist->xsize ;
  mri_pmap = compute_pmap(gca, transform, mri_intensities, target_label, mri_dist,pad,NULL) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_pmap, "p.mgz") ;
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

  vsize = mri_pmap->xsize ;
  for (x = 0 ; x < mri_ll->width; x++)
  {
    if (((x+1)%10) == 0)
      printf("x=%d\n", x) ;
    for (y = 0 ; y < mri_ll->height ; y++)
    {
      for (z = 0 ; z < mri_ll->depth ; z++)
      {
        double xw, yw, zw, pout, pin, p ;

        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
        if (val > pad)
        {
          MRIsetVoxVal(mri_ll, x, y, z, 0, -1000) ;
          continue ;
        }
        nx = MRIgetVoxVal(mri_dx, x, y, z, 0) ;
        ny = MRIgetVoxVal(mri_dy, x, y, z, 0) ;
        nz = MRIgetVoxVal(mri_dz, x, y, z, 0) ;
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
        {
          MRIsetVoxVal(mri_ll, x, y, z, 0, -10) ;
          continue ;
        }
#if 0
        V3_X(v1) = nx ; V3_Y(v1) = ny ; V3_Z(v1) = nz ;
        MatrixMultiply(m_vox2vox_vector, v1, v2) ; // to aseg/norm coords
        nx = V3_X(v2) ; ny = V3_Y(v2) ; nz = V3_Z(v2) ; 
#endif
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
#if 0
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xi = V3_X(v2) ; yi = V3_Y(v2) ; zi = V3_Z(v2) ; 
#else
        xi = x ; yi = y ; zi = z ;
#endif

        for (nin=0, pin = 0.0, d = -WSIZE/vsize ; d < 0 ; 
             d += vsize)
        {
          if (DZERO(d))
            continue ;
          xw = xi+d*nx ; yw = yi+d*ny ; zw = zi+d*nz ;
#if 0          
          MRIsampleVolume(mri_intensities, xw, yw, zw, &val);  // for debugging
          p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
#else
          MRIsampleVolume(mri_pmap, xw, yw, zw, &p) ;
#endif
          
          if (DZERO(p))
            p = 1e-5;
          nin++ ;
          pin += log(p) ;
        }
        
        d = 0.5 ;  // distance to sample outwards - d1 and d2 will be in tangent directions
        for (nout=0, pout = 0.0, d1 = -1/vsize ; d1 <= 1.0/vsize ; d1 += vsize)
          for (d2 = -1.0/vsize ; d2 <= 1.0/vsize ; d2 += vsize)
          {
            xw = xi+ d*nx + d1*e1x + d2*e2x ; yw = yi + d*ny + d1*e1y + d2*e2y ; zw = zi + d*nz + d1*e1z + d2*e2z ;
#if 0
            MRIsampleVolume(mri_intensities, xw, yw, zw, &val);
            p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
#else
            MRIsampleVolume(mri_pmap, xw, yw, zw, &p) ;
#endif
            
#define PSAT 0.1
            if (p < PSAT)  // don't care how unlikely it is if its past a certain point
              p = PSAT;
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
  MRIwrite(mri_pmap, "p.mgz") ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_ll, "ll.mgz") ;
    MRIwrite(mri_pin, "pin.mgz") ;
    MRIwrite(mri_pout, "pout.mgz") ;
  }
  MRIfree(&mri_mask) ; MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ; MRIfree(&mri_dist) ;
  MatrixFree(&m_vox2vox_vector) ; MRIfree(&mri_pin) ; MRIfree(&mri_pout) ;

  return(mri_ll) ;
}

static MRI *
compute_target_intensities_with_gcab(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, 
                                     HISTOGRAM **histograms, VERTEX_INFO *vi,float sigma, 
                                     TRANSFORM *transform, GCA *gca, int target_label, float resolution,
                                     GCAB *gcab) {
  int      nin, nout ;
  double   val, nx, ny, nz, xi, yi, zi, d, e1x, e2x, e1y, e2y, e1z, e2z, d1, d2, mag ;
  MRI      *mri_ll ;
  MRI      *mri_grad, *mri_dx, *mri_dy, *mri_dz, *mri_mask, 
           *mri_tmp, *mri_dist, *mri_pin, *mri_pout, *mri_pmap ;
  int       x, y, z, pad ;
  MATRIX    *m_vox2vox, *m_vox2vox_vector ;
  VECTOR    *v1, *v2 ;
  double    vsize ;

#define WSIZE 3

  mri_tmp = MRISfillInterior(mris, mri_labels->xsize/resolution, NULL) ;
  mri_tmp->c_r += mri_labels->c_r ;
  mri_tmp->c_a += mri_labels->c_a ;
  mri_tmp->c_s += mri_labels->c_s ;
  pad = ceil(3/mri_tmp->xsize) ;
  mri_mask = MRIextractRegionAndPad(mri_tmp, NULL, NULL, pad+1) ;
  mri_dist = MRIdistanceTransform(mri_mask, NULL, 1, nint(10/mri_mask->xsize), 
                                  DTRANS_MODE_SIGNED, NULL) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_dist, "d.mgz") ;

  mri_grad = MRIsobel(mri_dist, NULL, NULL) ;
  mri_dx = MRIxDerivative(mri_dist, NULL) ;
  mri_dy = MRIyDerivative(mri_dist, NULL) ;
  mri_dz = MRIzDerivative(mri_dist, NULL) ;
  MRIfree(&mri_grad) ; MRIfree(&mri_tmp) ; 
  mri_ll = MRIclone(mri_dist, NULL) ;
  mri_pin = MRIclone(mri_ll, NULL) ;
  mri_pout = MRIclone(mri_ll, NULL) ;

  mri_ll->outside_val = -100 ;
  pad *= mri_dist->xsize ;
  mri_pmap = compute_pmap_with_gcab(gca, transform, mri_intensities, mri_labels, target_label, mri_dist,pad, gcab, NULL) ;
  MRIwrite(mri_pmap, "p.mgz") ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_pmap, "p.mgz") ;
    MRIwrite(mri_dist, "d.mgz") ; MRIwrite(mri_mask, "m.mgz") ;
    MRIwrite(mri_dx, "dx.mgz") ; MRIwrite(mri_dy, "dy.mgz") ; 
    MRIwrite(mri_dz, "dz.mgz") ; 
  }
  MRIcopy(mri_pmap, mri_ll) ;
  return(mri_ll) ;

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_mask, mri_intensities) ;
  m_vox2vox_vector = MatrixCopy(m_vox2vox, NULL) ;
  *MATRIX_RELT(m_vox2vox_vector, 1, 4) = 0.0 ;
  *MATRIX_RELT(m_vox2vox_vector, 2, 4) = 0.0 ;
  *MATRIX_RELT(m_vox2vox_vector, 3, 4) = 0.0 ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;

  vsize = mri_pmap->xsize ;
  for (x = 0 ; x < mri_ll->width; x++)
  {
    if (((x+1)%10) == 0)
      printf("x=%d\n", x) ;
    for (y = 0 ; y < mri_ll->height ; y++)
    {
      for (z = 0 ; z < mri_ll->depth ; z++)
      {
        double xw, yw, zw, pout, pin, p ;

        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        val = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
        if (val > pad)
        {
          MRIsetVoxVal(mri_ll, x, y, z, 0, -1000) ;
          continue ;
        }
        nx = MRIgetVoxVal(mri_dx, x, y, z, 0) ;
        ny = MRIgetVoxVal(mri_dy, x, y, z, 0) ;
        nz = MRIgetVoxVal(mri_dz, x, y, z, 0) ;
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
        {
          MRIsetVoxVal(mri_ll, x, y, z, 0, -10) ;
          continue ;
        }
#if 0
        V3_X(v1) = nx ; V3_Y(v1) = ny ; V3_Z(v1) = nz ;
        MatrixMultiply(m_vox2vox_vector, v1, v2) ; // to aseg/norm coords
        nx = V3_X(v2) ; ny = V3_Y(v2) ; nz = V3_Z(v2) ; 
#endif
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
#if 0
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xi = V3_X(v2) ; yi = V3_Y(v2) ; zi = V3_Z(v2) ; 
#else
        xi = x ; yi = y ; zi = z ;
#endif

        for (nin=0, pin = 0.0, d = -WSIZE/vsize ; d < 0 ; 
             d += vsize)
        {
          if (DZERO(d))
            continue ;
          xw = xi+d*nx ; yw = yi+d*ny ; zw = zi+d*nz ;
#if 0          
          MRIsampleVolume(mri_intensities, xw, yw, zw, &val);  // for debugging
          p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
#else
          MRIsampleVolume(mri_pmap, xw, yw, zw, &p) ;
#endif
          
          if (DZERO(p))
            p = 1e-5;
          nin++ ;
          pin += log(p) ;
        }
        
        d = 0.5 ;  // distance to sample outwards - d1 and d2 will be in tangent directions
        for (nout=0, pout = 0.0, d1 = -1/vsize ; d1 <= 1.0/vsize ; d1 += vsize)
          for (d2 = -1.0/vsize ; d2 <= 1.0/vsize ; d2 += vsize)
          {
            xw = xi+ d*nx + d1*e1x + d2*e2x ; yw = yi + d*ny + d1*e1y + d2*e2y ; zw = zi + d*nz + d1*e1z + d2*e2z ;
#if 0
            MRIsampleVolume(mri_intensities, xw, yw, zw, &val);
            p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
#else
            MRIsampleVolume(mri_pmap, xw, yw, zw, &p) ;
#endif
            
#define PSAT 0.1
            if (p < PSAT)  // don't care how unlikely it is if its past a certain point
              p = PSAT;
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
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_ll, "ll.mgz") ;
    MRIwrite(mri_pin, "pin.mgz") ;
    MRIwrite(mri_pout, "pout.mgz") ;
  }
  MRIfree(&mri_mask) ; MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ; MRIfree(&mri_dist) ;
  MatrixFree(&m_vox2vox_vector) ; MRIfree(&mri_pin) ; MRIfree(&mri_pout) ;

  return(mri_ll) ;
}

static int
compute_target_labels(MRI_SURFACE *mris, MRI *mri_labels, MRI *mri_intensities, HISTOGRAM **histograms, VERTEX_INFO *vi) {
  VERTEX   *v ;
  int      vno, l_out, l_in ;
  double   xv, yv, zv, val, nx, ny, nz ;
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

#define MAX_LL 1  // was 10
#define MAX_DIST 5

static double
externalLLSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int         vno ;
  VERTEX      *v ;
  double      sse/*, error*/, d ;
  double      val, xv, yv, zv, nx, ny, nz, max_d, max_ll, x2, y2, z2, dist, 
              prev_dist, xvd, yvd, zvd, nxd, nyd, nzd, mag, sigma=4.0 ;
  MRI         *mri_ll = parms->mri_ll, *mri_dist ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  MRISsmoothSurfaceNormals(mris, 10) ;
  //  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  //  MRISaverageVertexPositions(mris, 10) ;
  mri_dist = MRIScomputeDistanceToSurface(mris, NULL, mri_ll->xsize)  ;
#if 0
  mri_kernel = MRIgaussian1d(2, 100) ;
  mri_smooth = MRIclone(mri_dist, NULL) ;
  MRIconvolveGaussian(mri_dist, mri_smooth, mri_kernel) ;
  MRIcopy(mri_smooth, mri_dist) ;
  MRIfree(&mri_kernel) ; MRIfree(&mri_smooth) ;
#endif

  //  MRISrestoreVertexPositions(mris, TMP_VERTICES) ; 
  //  MRIScomputeMetricProperties(mris) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_dist, "d.mgz") ;
  MRIreInitCache(mri_dist) ;
  MRIreInitCache(parms->mri_ll) ;
  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, mri_dist, &xvd, &yvd, &zvd) ;
    MRISvertexToVoxel(mris, v, parms->mri_ll, &xv, &yv, &zv) ;
    MRIsampleVolume(mri_ll, xv, yv, zv, &val);
    max_ll = val; max_d = 0 ;

    MRISvertexNormalInVoxelCoords(mris, mri_ll, vno, &nx, &ny, &nz) ;
    MRISvertexNormalInVoxelCoords(mris, mri_dist, vno, &nxd, &nyd, &nzd) ;
    prev_dist = 0 ;
    for (d = -mri_ll->xsize/2 ; d >= -MAX_DIST/mri_ll->xsize ; d -= mri_ll->xsize/2)
    {
      x2 = xvd + d*nxd ; y2 = yvd + d*nyd ; z2 = zvd + d*nzd ;
      MRIsampleVolume(mri_dist, x2, y2, z2, &dist) ; dist = fabs(dist) ;
      MRIsampleVolumeDerivativeScale(mri_dist, x2, y2, z2, nxd, nyd, nzd, &mag, sigma) ;
#if 1
      if (mag < 0)
        break ;
#else
      if (dist < prev_dist) // closer to the opposite side of the surface
        break ;
#endif
      prev_dist = dist ;
      x2 = xv + d*nx ; y2 = yv + d*ny ; z2 = zv + d*nz ;
      MRIsampleVolume(mri_ll, x2, y2, z2, &val);
      if (val > max_ll)
      {
        max_ll = val ;
        max_d = d ;
      }
    }
    prev_dist = 0 ;
    for (d = 0 ; d <= MAX_DIST/mri_ll->xsize ; d += mri_ll->xsize/2)
    {
      x2 = xvd + d*nxd ; y2 = yvd + d*nyd ; z2 = zvd + d*nzd ;
      MRIsampleVolume(mri_dist, x2, y2, z2, &dist) ; dist = fabs(dist) ;
      MRIsampleVolumeDerivativeScale(mri_dist, x2, y2, z2, nxd, nyd, nzd, &mag, sigma) ;
#if 1
      if (mag < 0)
        break ;
#else
      if (dist < prev_dist) // closer to the opposite side of the surface
        break ;
#endif
      prev_dist = dist ;

      x2 = xv + d*nx ; y2 = yv + d*ny ; z2 = zv + d*nz ;
      MRIsampleVolume(mri_ll, x2, y2, z2, &val);
      if (val > max_ll)
      {
        max_ll = val ;
        max_d = d ;
      }
    }
    if (vno == Gdiag_no) 
      printf("E_extern: max_ll = %2.5f, max_d = %2.2f\n", max_ll, max_d) ;
    sse += max_d*max_d ;
  }

  MRIfree(&mri_dist) ;
  return(sse * parms->l_external) ;
}

static double
externalLLGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) 
{
  int         vno ;
  VERTEX      *v ;
  double      sse, error, d ;
  double      val, xv, yv, zv, nx, ny, nz, max_d, max_ll, x2, y2, z2, dist, 
              prev_dist, xvd, yvd, zvd, nxd, nyd, nzd, mag, sigma = 4.0 ;
  MRI         *mri_ll = parms->mri_ll, *mri_dist ;

  if (FZERO(parms->l_external))
    return(0.0) ;

  MRISsmoothSurfaceNormals(mris, 10) ;
  //  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  //  MRISaverageVertexPositions(mris, 10) ;
  mri_dist = MRIScomputeDistanceToSurface(mris, NULL, mri_ll->xsize)  ;
  //  MRISrestoreVertexPositions(mris, TMP_VERTICES) ; 
  //  MRIScomputeMetricProperties(mris) ;
#if 0
  mri_kernel = MRIgaussian1d(2, 100) ;
  mri_smooth = MRIclone(mri_dist, NULL) ;
  MRIconvolveGaussian(mri_dist, mri_smooth, mri_kernel) ;
  MRIcopy(mri_smooth, mri_dist) ;
  MRIfree(&mri_kernel) ; MRIfree(&mri_smooth) ;
#endif
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_dist, "d.mgz") ;
  for (sse = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    MRISvertexToVoxel(mris, v, mri_dist, &xvd, &yvd, &zvd) ;
    MRISvertexToVoxel(mris, v, mri_ll, &xv, &yv, &zv) ;
    MRIsampleVolume(mri_ll, xv, yv, zv, &val);
    if (vno == Gdiag_no) 
      printf("l_extern: v %d, ll vox (%2.1f, %2.1f, %2.1f) = %2.6f\n", 
             vno, xv, yv, zv, val) ;

    if (val > MAX_LL)
      val = MAX_LL ;
    error = MAX_LL-val ;
    sse += error*error ;

    MRISvertexNormalInVoxelCoords(mris, mri_ll, vno, &nx, &ny, &nz) ;
    MRISvertexNormalInVoxelCoords(mris, mri_dist, vno, &nxd, &nyd, &nzd) ;
    v->d = error ;
    max_ll = val; max_d = 0 ;

    prev_dist = 0 ;
    for (d = -mri_ll->xsize/2 ; d >= -MAX_DIST/mri_ll->xsize ; d -= mri_ll->xsize/2)
    {
      x2 = xvd + d*nxd ; y2 = yvd + d*nyd ; z2 = zvd + d*nzd ;
      MRIsampleVolume(mri_dist, x2, y2, z2, &dist) ; dist = fabs(dist) ;
      MRIsampleVolumeDerivativeScale(mri_dist, x2, y2, z2, nxd, nyd,nzd, &mag, sigma) ;
#if 1
      if (mag < 0)
        break ;
#else
      if (dist < prev_dist) // closer to the opposite side of the surface
        break ;
#endif
      prev_dist = dist ;

      x2 = xv + d*nx ; y2 = yv + d*ny ; z2 = zv + d*nz ;
      MRIsampleVolume(mri_ll, x2, y2, z2, &val);
      if (val > max_ll)
      {
        max_ll = val ;
        max_d = d ;
      }
    }
    prev_dist = 0 ;
    for (d = 0 ; d <= MAX_DIST/mri_ll->xsize ; d += mri_ll->xsize/2)
    {
      x2 = xvd + d*nxd ; y2 = yvd + d*nyd ; z2 = zvd + d*nzd ;
      MRIsampleVolume(mri_dist, x2, y2, z2, &dist) ; dist = fabs(dist) ;
      MRIsampleVolumeDerivativeScale(mri_dist, x2, y2, z2, nxd, nyd, nzd, &mag, sigma) ;
#if 1
      if (mag < 0)
        break ;
#else
      if (dist < prev_dist) // closer to the opposite side of the surface
        break ;
#endif
      prev_dist = dist ;

      x2 = xv + d*nx ; y2 = yv + d*ny ; z2 = zv + d*nz ;
      MRIsampleVolume(mri_ll, x2, y2, z2, &val);
      if (val > max_ll)
      {
        max_ll = val ;
        max_d = d ;
      }
    }
    if (vno == Gdiag_no)
      printf("l_extern: v %d, voxel (%2.1f, %2.1f, %2.1f), N=(%2.2f, %2.2f, %2.2f), max ll = %2.4f at distance %2.2f\n",
             vno, xv, yv, zv, nx, ny, nz, max_ll, max_d) ;
    mag = max_d ;
    if (fabs(mag) > 1)
      mag /= fabs(mag) ; 
    mag *= parms->l_external ;
    if (vno == Gdiag_no)
      printf("l_extern:   dx_surf = (%2.2f, %2.2f, %2.2f)\n", mag*v->nx, mag*v->ny, mag*v->nz) ;

    if (error > v->d)
      DiagBreak() ;
    v->dx += mag*v->nx ; 
    v->dy += mag*v->ny ; 
    v->dz += mag*v->nz ; 
    if (!devFinite(v->dx) || !devFinite(v->dy) || !devFinite(v->dz))
      DiagBreak() ;
  }

  MRIfree(&mri_dist) ;
  return(sse * parms->l_external) ;
}


static double
externalLLRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms) {
  return(sqrt(externalLLSSE(mris, parms)/mris->nvertices)) ;
}

static GCA *
make_gca(char *label_vol_name, char *intensity_vol_name, TRANSFORM *transform)
{
	MRI        *mri_labels, *mri_intensity, *mri_tmp ;
	GCA        *gca ;
	MRI_REGION box ;

	mri_labels = MRIread(label_vol_name) ;
	if (mri_labels == NULL)
		ErrorExit(ERROR_NOFILE, "%s: could not read label vol from %s",
							Progname, label_vol_name) ;

	mri_intensity = MRIread(label_vol_name) ;
	if (mri_intensity == NULL)
		ErrorExit(ERROR_NOFILE, "%s: could not read intensity vol from %s",
							Progname, intensity_vol_name) ;

  MRIboundingBox(mri_labels, 0, &box) ;
	mri_tmp = MRIextractRegionAndPad(mri_labels, NULL, &box, 10) ;
	MRIfree(&mri_labels) ; mri_labels = mri_tmp ;
	mri_tmp = MRIextractRegionAndPad(mri_intensity, NULL, &box, 10) ;
	MRIfree(&mri_intensity) ; mri_intensity = mri_tmp ;

	gca = GCAalloc(1, 8, 8, mri_labels->width, mri_labels->height, 
								 mri_labels->depth, 0) ;
	GCAtrain(gca, mri_intensity, mri_labels, transform, NULL, 0) ;
	GCAcompleteMeanTraining(gca) ;
	GCAtrainCovariances(gca, mri_intensity, mri_labels, transform) ;
	GCAcompleteCovarianceTraining(gca) ;
	MRIfree(&mri_labels) ;
	MRIfree(&mri_intensity) ;
	return(gca) ;
}

static MRI *
compute_pmap(GCA *gca, TRANSFORM *transform, MRI *mri_intensities, int target_label, 
             MRI *mri_dist, double pad, MRI *mri_pmap)
{
  int       x, y, z ;
  double    p, dist, xw, yw, zw ;
  VECTOR    *v1, *v2 ;
  MATRIX    *m_vox2vox ;

  if (mri_pmap == NULL)
  {
    mri_pmap = MRIalloc(mri_dist->width,mri_dist->height,mri_dist->depth,MRI_FLOAT);
    if (mri_pmap == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate (%d, %d, %d) pmap",
                Progname, mri_dist->width,mri_dist->height,mri_dist->depth) ;
    MRIcopyHeader(mri_dist, mri_pmap) ;
  }
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_dist, mri_intensities) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;

  for (x = 0 ; x < mri_pmap->width; x++)
  {
    if (((x+1)%10) == 0)
      printf("x=%d\n", x) ;
    for (y = 0 ; y < mri_pmap->height ; y++)
    {
      for (z = 0 ; z < mri_pmap->depth ; z++)
      {
        dist = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
        if (dist > pad)
          continue ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xw = V3_X(v2) ; yw = V3_Y(v2) ; zw = V3_Z(v2) ; 
        p = GCAcomputeLabelLikelihood(gca, transform, mri_intensities, xw, yw, zw, target_label);
        MRIsetVoxVal(mri_pmap, x, y, z, 0, p) ;
      }
    }
  }

  MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
  return(mri_pmap) ;
}

static MRI *
compute_pmap_with_gcab(GCA *gca, TRANSFORM *transform, MRI *mri_intensities, MRI *mri_labels, int target_label, 
                       MRI *mri_dist, double pad, GCAB *gcab, MRI *mri_pmap)
{
  int       x, y, z ;
  double    p, dist, xw, yw, zw, nx, ny, nz, mag, xw1, yw1, zw1,pout,pin,pgrad;
  VECTOR    *v1, *v2 ;
  MATRIX    *m_vox2vox ;
  MRI       *mri_dist1mm, *mri_pout, *mri_pin, *mri_pgrad, *mri_tmp ;


  if (mri_pmap == NULL)
  {
    mri_pmap = MRIalloc(mri_dist->width,mri_dist->height,mri_dist->depth,MRI_FLOAT);
    if (mri_pmap == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate (%d, %d, %d) pmap",
                Progname, mri_dist->width,mri_dist->height,mri_dist->depth) ;
    MRIcopyHeader(mri_dist, mri_pmap) ;
  }
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_dist, mri_intensities) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = VECTOR_ELT(v2, 4) = 1.0 ;

  mri_tmp = MRIclone(mri_labels, NULL) ;
  MRIcopyLabel(mri_labels, mri_tmp, target_label) ;
  mri_dist1mm = MRIdistanceTransform(mri_tmp, NULL, target_label, 10, DTRANS_MODE_SIGNED, NULL) ;
  MRIfree(&mri_tmp) ;

  mri_pgrad = MRIclone(mri_pmap, NULL) ;
  mri_pin = MRIclone(mri_pmap, NULL) ;
  mri_pout = MRIclone(mri_pmap, NULL) ;
  for (x = 0 ; x < mri_pmap->width; x++)
  {
    if (((x+1)%10) == 0)
      printf("x=%d\n", x) ;
    for (y = 0 ; y < mri_pmap->height ; y++)
    {
      for (z = 0 ; z < mri_pmap->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        dist = MRIgetVoxVal(mri_dist, x, y, z, 0) ;
        if (dist > pad)
          continue ;
        MRIsampleVolumeGradient(mri_dist, x, y, z, &nx, &ny, &nz) ;
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
          continue ;
        nx /= mag ; ny /= mag ; nz /= mag ;

        V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xw = V3_X(v2) ; yw = V3_Y(v2) ; zw = V3_Z(v2) ; 

        V3_X(v1) = x+2*nx ; V3_Y(v1) = y+2*ny ; V3_Z(v1) = z+2*nz ;
        MatrixMultiply(m_vox2vox, v1, v2) ; // to aseg/norm coords
        xw1 = V3_X(v2) ; yw1 = V3_Y(v2) ; zw1 = V3_Z(v2) ; 
        nx = xw1-xw ;  ny = yw1-yw ;  nz = zw1-zw ; 
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
          continue ;
        nx /= mag ; ny /= mag ; nz /= mag ;

        p = GCABgetProbability(gcab, mri_intensities, mri_dist1mm, transform, xw, yw, zw, nx, ny, nz) ;
        pout = GCABgetPout(gcab, mri_intensities, mri_dist1mm, transform, xw, yw, zw, nx, ny, nz) ;
        pin = GCABgetPin(gcab, mri_intensities, mri_dist1mm, transform, xw, yw, zw, nx, ny, nz) ;
        pgrad = GCABgetPgrad(gcab, mri_intensities, mri_dist1mm, transform, xw, yw, zw, nx, ny, nz) ;
        MRIsetVoxVal(mri_pgrad, x, y, z, 0, pgrad) ;
        MRIsetVoxVal(mri_pin, x, y, z, 0, pin) ;
        MRIsetVoxVal(mri_pout, x, y, z, 0, pout) ;
        p = pow(pin * pout * pgrad, 0.3333) ;
        MRIsetVoxVal(mri_pmap, x, y, z, 0, p) ;
      }
    }
  }

  MRIwrite(mri_pin, "pin.mgz") ;
  MRIwrite(mri_pout, "pout.mgz") ;
  MRIwrite(mri_pgrad, "pgrad.mgz") ;
  MRIfree(&mri_dist1mm) ;
  MRIfree(&mri_pout) ;
  MRIfree(&mri_pin) ;
  MRIfree(&mri_pgrad) ;
  MatrixFree(&m_vox2vox) ; VectorFree(&v1) ; VectorFree(&v2) ;
  return(mri_pmap) ;
}
static int
compute_gradient_target_positions(MRI_SURFACE *mris, MRI *mri_intensities, VERTEX_INFO *vi, 
                                  float current_sigma)
{
  return(NO_ERROR) ;
}

static double
externalGradGradient(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  return(0.0) ;
}
static double
externalGradSSE(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  return(0.0) ;
}
static double
externalGradRMS(MRI_SURFACE *mris,INTEGRATION_PARMS *parms)
{
  return(0.0) ;
}

