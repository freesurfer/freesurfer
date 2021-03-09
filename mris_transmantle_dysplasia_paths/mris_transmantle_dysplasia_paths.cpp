/**
 * @brief MCMC for computing posterior of splines connecting cortex with ventricle
 *
 * Fit a Catmull Rom spline to each point in the cortex, initializing it with a connection along the
 * shortest path to the lateral ventricles, then use MCMC to estimate the posterior distribution.
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
#include "const.h"
#include "timer.h"
#include "version.h"
#include "transform.h"
#include "mrisurf.h"
#include "cma.h"
#include "icosahedron.h"
#include "voxlist.h"
#include "pdf.h"
#include "tritri.h"

#include "romp_support.h"


int main(int argc, char *argv[]) ;
static VOXEL_LIST *compute_path_to_ventricles(MRI_SURFACE *mris, int vno, MRI *mri_ventricle_dist_grad, MRI *mri_aseg) ;
static int get_option(int argc, char *argv[]) ;
static MRI *compute_migration_probabilities(MRI_SURFACE *mris, MRI *mri, MRI *mri_aseg, TRANSFORM *xform, MRI *mri_pvals, int spline_control_points, int mcmc_samples, int read_flag) ;
static VOXEL_LIST *find_optimal_spline(VOXEL_LIST *vl, MRI *mri_intensity, MRI *mri_aseg, MRI *mri_wm_dist, double gm_mean, int nsamples,MRI_SURFACE *mris, int spline_control_points, int mcmc_samples, double spline_length_penalty, double spline_nonwm_penalty, double spline_interior_penalty, MRI *mri_posterior, VOXEL_LIST *vl_posterior,
				       double *pentropy, MRI *mri_total_posterior) ;

static MRI *compute_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg, MRI *mri_posterior_on_spline) ;

static MRI *compute_filtered_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg, MRI *mri_posterior_on_spline) ;


#define SPLINE_WM_DIST    0x0001
#define SPLINE_ABOVE      0x0002
#define SPLINE_BELOW      0x0004
#define SPLINE_LENGTH     0x0008
#define SPLINE_SIGNED     0x0010

//static int energy_flags = SPLINE_WM_DIST | SPLINE_LENGTH | SPLINE_ABOVE | SPLINE_BELOW ;
static int energy_flags = SPLINE_WM_DIST | SPLINE_LENGTH | SPLINE_BELOW | SPLINE_ABOVE ;

static double compute_spline_energy(VOXEL_LIST *vl, MRI *mri_intensity, MRI *mri_aseg, MRI *mri_wm_dist, double mean, int which, double spline_length_penalty, double spline_nonwm_penalty, double spline_interior_penalty) ;

const char *Progname ;
static void usage_exit(int code) ;

static double noise = .1 ;
static int spline_control_points = 5 ;
static double spline_length_penalty = 5 ;
static double spline_nonwm_penalty = 200 ;
static double spline_interior_penalty = 1000 ;
static double max_wm_dist = -2.5 ;
static int read_flag = 0 ;

static char *label_name = NULL ;
static double proposal_sigma = 5.0 ; // stddev of noise distribution
static double acceptance_sigma = .5 ;
static int mcmc_samples = 10000 ;

static int randomize_data = 0 ;

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI_SURFACE  *mris ;
  MRI          *mri, *mri_aseg, *mri_pvals ;
  TRANSFORM    *xform ;

  nargs = handleVersionOption(argc, argv, "mris_transmantle_dysplasia_paths");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  setRandomSeed(0L) ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;
  
  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load surface from %s", Progname, argv[1]) ;
  
  if (label_name)
  {
    LABEL *area ;
    area = LabelRead(NULL, label_name) ;
    if (area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label from %s", Progname, label_name) ;
    LabelRipRestOfSurface(area,mris) ;
    LabelFree(&area) ;
  }

  mri = MRIread(argv[2]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load intensity volume from %s", Progname, argv[2]) ;
  mri_aseg = MRIread(argv[3]) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load aseg volume from %s", Progname, argv[3]) ;
  if (!stricmp(argv[4], "identity.nofile"))
    xform = TransformAlloc(LINEAR_VOX_TO_VOX, mri) ;
  else
    xform = TransformRead(argv[4]) ;
  if (xform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load transform from %s", Progname, argv[4]) ;

  mri_pvals = compute_migration_probabilities(mris,mri,mri_aseg,xform,NULL,spline_control_points,mcmc_samples,read_flag);

  printf("writing path log probabilities to %s\n", argv[5]) ;
  MRIwrite(mri_pvals, argv[5]) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "simulated neural migration took %d minutes and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "RAND"))
  {
    randomize_data = 1 ;
    printf("randomizing input intensities\n") ;
  }
  else switch (toupper(*option)) {
  case 'R':
    read_flag = 1 ;
    break ;
  case 'L':
    label_name = argv[2] ;
    nargs = 1 ;
    printf("reading cortex label from %s\n", label_name) ;
    break ;
  case 'P':
    mcmc_samples = atoi(argv[2]) ;
    printf("sampling from %d paths in MCMC\n", mcmc_samples) ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
    break ;
  case 'N':
    noise = atof(argv[2]) ;
    printf("using noise level = %2.2f to deflect ventricle distance transform negative gradient vectors\n",noise) ;
    nargs = 1 ;
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
  printf("usage: %s [options] <surface> <aseg volume> <intensity volume> <xform> <output volume>\n", 
         Progname) ;
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  printf("\tn - noise-sensitivity normalize inverse (default=1)") ;
  exit(code) ;
}



static MRI *
MRIbuildPossibleMigrationPaths(MRI *mri_aseg, MRI *mri_wm_interior, MRI *mri_paths)
{
  int  x, y,z, nadded, label, n ;
  MRI  *mri_dilated = NULL ;

  if (mri_paths == NULL)
    mri_paths = MRIclone(mri_aseg, NULL);

  MRIcopyLabel(mri_aseg, mri_paths, Left_Lateral_Ventricle) ;
  MRIcopyLabel(mri_aseg, mri_paths, Right_Lateral_Ventricle) ;
  MRIbinarize(mri_paths, mri_paths, 1, 0, 1) ;

  n = 1 ;
  do
  {
    nadded = 0 ;
    mri_dilated = MRIdilate(mri_paths, mri_dilated) ;
    for (x = 0 ; x < mri_aseg->width ; x++)
      for (y = 0 ; y < mri_aseg->height ; y++)
	for (z = 0 ; z < mri_aseg->depth ; z++)
	{
	  if (MRIgetVoxVal(mri_wm_interior, x, y, z, 0) == 0)
	    continue ;
	  label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	  if (IS_WMH(label) == 0)
	    continue ;
	  if (MRIgetVoxVal(mri_dilated, x, y, z, 0) == 0)  // not next to a voxel in the current set
	    continue ;
	  if (MRIgetVoxVal(mri_paths, x, y, z, 0))        // already in the current set
	    continue ;
	  MRIsetVoxVal(mri_paths, x, y, z, 0, n+1) ;
	  nadded++ ;
	}
//    printf("round %d: %d voxels added\n", n, nadded) ;
    n++ ;
  } while (nadded > 0); 

  n = 0 ;
  do
  {
    nadded = 0 ;
    mri_dilated = MRIdilate(mri_paths, mri_dilated) ;
    for (x = 0 ; x < mri_aseg->width ; x++)
      for (y = 0 ; y < mri_aseg->height ; y++)
	for (z = 0 ; z < mri_aseg->depth ; z++)
	{
	  if (MRIgetVoxVal(mri_dilated, x, y, z, 0) == 0)  // not next to a voxel in the current set
	    continue ;
	  if (MRIgetVoxVal(mri_paths, x, y, z, 0))        // already in the current set
	    continue ;
	  MRIsetVoxVal(mri_paths, x, y, z, 0, MRIgetVoxVal(mri_dilated, x, y, z, 0)+1) ;
	  nadded++ ;
	}
    printf("round %d: %d voxels added\n", n, nadded) ;
    n++ ;
  } while (nadded > 0 && n < 5); 

  MRIfree(&mri_dilated) ;
  return(mri_paths) ;
}


static MRI *
compute_migration_probabilities(MRI_SURFACE *mris, MRI *mri_intensity, MRI *mri_aseg,
				TRANSFORM *xform, MRI *mri_pvals, int spline_control_points, int mcmc_samples, int read_flag)
{
  int    vno, x, y, z, nvox, label, nvox_wm, n ;
  MATRIX *m_intensity2aseg, *v_aseg, *v_intensity ;
  double gm_mean, gm_var, val, wm_mean, wm_var, entropy=1.0 ;
  VERTEX  *v;
  MRI    *mri_wm_interior, *mri_wm_dist, *mri_kernel, *mri_tmp, *mri_posterior, 
    *mri_entropy, *mri_total_posterior,
    *mri_eroded_aseg, *mri_possible_migration_paths,*mri_path_grad, *mri_splines;
  VOXEL_LIST *vl, *vl_spline, *vl_posterior ;
  char fname[STRLEN], path[STRLEN] ;
  MRI  *mri_posterior_on_splines, *mri_filtered_posterior_on_spline ;


  wm_mean = gm_mean = wm_var = gm_var = 0 ;

  if (read_flag)
  {
    char fname[STRLEN], path[STRLEN] ;

    sprintf(fname, "%s/%s.splines.mgz", 
	    FileNamePath(mris->fname, path), 
	    mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    printf("reading optimal splines from %s\n", fname) ;
    mri_splines = MRIread(fname) ;
    if (mri_splines == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read splines from %s", Progname, fname) ;
  }
  else
    mri_splines = MRIalloc(mris->nvertices, spline_control_points, 3, MRI_FLOAT) ;
  mri_entropy = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;
  mri_posterior = MRIcloneDifferentType(mri_aseg, MRI_FLOAT) ;
  mri_total_posterior = MRIcloneDifferentType(mri_aseg, MRI_FLOAT) ;
  mri_eroded_aseg  = MRIclone(mri_aseg, NULL) ;
  MRIcopyLabel(mri_aseg, mri_eroded_aseg, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_eroded_aseg, Right_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_eroded_aseg, WM_hypointensities) ;
  MRIbinarize(mri_eroded_aseg, mri_eroded_aseg, 1, 0, 1) ;
  MRIerode(mri_eroded_aseg, mri_eroded_aseg) ;
  MRIreplaceValuesOnly(mri_eroded_aseg, mri_eroded_aseg, 1, Left_Cerebral_White_Matter) ;
  mri_wm_interior = MRIcloneDifferentType(mri_intensity, MRI_FLOAT) ;
  MRISfillInterior(mris, mri_intensity->xsize, mri_wm_interior) ;
  mri_possible_migration_paths = MRIbuildPossibleMigrationPaths(mri_aseg, mri_wm_interior, NULL);

  mri_wm_dist = MRIdistanceTransform(mri_wm_interior, NULL, 1, 5, DTRANS_MODE_SIGNED, NULL);
  MRIfree(&mri_wm_interior) ;
  if (mri_pvals == NULL)
    mri_pvals = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;

  v_intensity = VectorAlloc(4, MATRIX_REAL) ;
  v_aseg = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_aseg, 4) = 1.0 ;
  VECTOR_ELT(v_intensity, 4) = 1.0 ;
  m_intensity2aseg = MRIgetVoxelToVoxelXform(mri_intensity, mri_aseg) ;
  
  for (nvox_wm = nvox = x = 0 ; x < mri_intensity->width ; x++)
    for (y = 0 ; y < mri_intensity->height ; y++)
      for (z = 0 ; z < mri_intensity->depth ; z++)
      {
	V3_X(v_intensity) = x ;
	V3_Y(v_intensity) = y ;
	V3_Z(v_intensity) = z ;
	MatrixMultiply(m_intensity2aseg, v_intensity, v_aseg) ;
	label = MRIgetVoxVal(mri_aseg, nint(V3_X(v_aseg)), nint(V3_Y(v_aseg)), nint(V3_Z(v_aseg)), 0) ;
	if (IS_CORTEX(label))
	{
	  nvox++ ;
	  val = MRIgetVoxVal(mri_intensity, x, y, z, 0) ;
	  gm_mean += val ; gm_var += val*val ;
	}
	if (IS_WM(label) && MRIgetVoxVal(mri_wm_dist, x, y, z, 0) <= -1)
	{
	  nvox_wm++ ;
	  val = MRIgetVoxVal(mri_intensity, x, y, z, 0) ;
	  wm_mean += val ; wm_var += val*val ;
	}
      }

  gm_mean /= nvox ; gm_var = gm_var / nvox - gm_mean*gm_mean ;
  wm_mean /= nvox_wm ; wm_var = wm_var / nvox_wm - wm_mean*wm_mean ;
  printf("%d gm voxels found: %2.1f +- %2.1f\n", nvox, gm_mean, sqrt(gm_var)) ;
  printf("%d wm voxels found: %2.1f +- %2.1f\n", nvox_wm, wm_mean, sqrt(wm_var)) ;

  if (randomize_data)
  {
    printf("randomizing input intensities\n") ;
    for (nvox_wm = nvox = x = 0 ; x < mri_intensity->width ; x++)
      for (y = 0 ; y < mri_intensity->height ; y++)
	for (z = 0 ; z < mri_intensity->depth ; z++)
	{
	  V3_X(v_intensity) = x ;
	  V3_Y(v_intensity) = y ;
	  V3_Z(v_intensity) = z ;
	  MatrixMultiply(m_intensity2aseg, v_intensity, v_aseg) ;
	  label = MRIgetVoxVal(mri_aseg, nint(V3_X(v_aseg)), nint(V3_Y(v_aseg)), nint(V3_Z(v_aseg)), 0) ;
	  if (IS_CORTEX(label))
	    val = sqrt(gm_var) * PDFgaussian() + gm_mean ;
	  else if (IS_WMH(label))
	    val = sqrt(wm_var) * PDFgaussian() + wm_mean ;
	  else
	    continue ;
	  MRIsetVoxVal(mri_intensity, x, y, z, 0, val) ;
	}
  }

  mri_kernel = MRIgaussian1d(1, 100) ;
  mri_tmp = MRIconvolveGaussian(mri_possible_migration_paths, NULL, mri_kernel) ;
  MRIfree(&mri_possible_migration_paths) ; mri_possible_migration_paths = mri_tmp ;
  MRIfree(&mri_kernel) ;
  mri_path_grad = MRIsobel(mri_possible_migration_paths, NULL, NULL) ;
  MRInormalizeSequence(mri_path_grad, 1.0) ;
  vl_posterior = VLSTalloc(mri_aseg->width*mri_aseg->height*mri_aseg->depth) ;
  if (Gdiag_no >= 0)
  {
    printf("writing paths and gradients\n") ;
    MRIwrite(mri_path_grad, "pg.mgz") ;
    MRIwrite(mri_possible_migration_paths, "p.mgz") ;
  }
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  v = NULL ; n = 0 ;
  vl = vl_spline = NULL ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(n, v, vl, vl_spline, entropy, gm_mean, mcmc_samples) shared(mri_intensity, mri_aseg, mri_path_grad,mri_splines) schedule(static,1)
#endif
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    ROMP_PFLB_begin
    
    if ((vno % (mris->nvertices/200)) == 0)
      printf("processed %d of %d: %2.1f%%\n", vno, mris->nvertices, 100.0*vno/(float)mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (read_flag == 0)
    {
      vl = compute_path_to_ventricles(mris, vno, mri_path_grad, mri_aseg) ;
      if (vl == NULL)
	ROMP_PFLB_continue ;
      vl_spline = find_optimal_spline(vl, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, mcmc_samples,mris,spline_control_points, mcmc_samples, spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty,
				      mri_posterior, vl_posterior, &entropy, mri_total_posterior);
      MRIsetVoxVal(mri_entropy, vno, 0, 0, 0, -log(entropy)) ;
      if (vl_spline == NULL)
      {
	VLSTfree(&vl) ;
	ROMP_PFLB_continue ;
      }
      if (vno == Gdiag_no)
      {
	double     init_energy, energy ;
	VOXEL_LIST *vl_init_spline, *vl_interp ;
	
	energy = compute_spline_energy(vl_spline, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, energy_flags,
				       spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty) ;
	printf("vno %d: final energy = %2.2f\n", vno, energy) ;
	VLSTwriteLabel(vl_spline, "spline.optimal.cpts.label", mris, mri_aseg) ;
	vl_interp = VLSTinterpolate(vl_spline, .5) ;
	VLSTwriteLabel(vl_interp, "spline.optimal.label", mris, mri_aseg) ;
	val = compute_spline_energy(vl_spline, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, energy_flags,
				    spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty) ;
	
	vl_init_spline = VLSTsplineFit(vl, spline_control_points) ;
	VLSTwriteLabel(vl, "spline.init.label", mris, mri_aseg) ;
	VLSTwriteLabel(vl_init_spline, "spline.init.cpts.label", mris, mri_aseg) ;
	init_energy = compute_spline_energy(vl_init_spline, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, energy_flags,
					    spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty) ;
	printf("v %d: initial energy = %2.2f\n", vno, init_energy) ;
	DiagBreak() ;
	VLSTfree(&vl_interp) ;
	VLSTfree(&vl_init_spline) ;
      }
      VLSTfree(&vl) ;
    }
    else
      vl_spline = VLSTfromMRI(mri_splines, vno) ;
    val = compute_spline_energy(vl_spline, mri_intensity, mri_eroded_aseg, mri_wm_dist, wm_mean, 
				SPLINE_ABOVE, spline_length_penalty, 0, spline_interior_penalty) ;


    for (n = 0 ; n < spline_control_points ; n++)
    {
      MRIsetVoxVal(mri_splines, vno, n, 0, 0, vl_spline->xd[n]) ;
      MRIsetVoxVal(mri_splines, vno, n, 1, 0, vl_spline->yd[n]) ;
      MRIsetVoxVal(mri_splines, vno, n, 2, 0, vl_spline->zd[n]) ;
    }
    VLSTfree(&vl_spline) ;
//    MRIsetVoxVal(mri_pvals, vno, 0, 0, 0, exp(val/100.0)) ;
    MRIsetVoxVal(mri_pvals, vno, 0, 0, 0, val) ;
  
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  if (read_flag == 0)
  {
    if (randomize_data)
      sprintf(fname, "%s/%s.splines.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.splines.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    printf("writing optimal splines to %s\n", fname) ;
    MRIwrite(mri_splines, fname) ;

    if (randomize_data)
      sprintf(fname, "%s/%s.entropy.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.entropy.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 

    printf("writing entropies to to %s\n", fname) ;
    MRIwrite(mri_entropy, fname) ;
    if (randomize_data)
      sprintf(fname, "%s/%s.posterior.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.posterior.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    printf("writing integrated posterior to %s\n", fname) ;
    MRIwrite(mri_total_posterior, fname) ;
  }
  else    // read previously computed data in
  {
    if (randomize_data)
      sprintf(fname, "%s/%s.splines.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.splines.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    printf("reading optimal splines from %s\n", fname) ;
    mri_splines = MRIread(fname) ;

    if (randomize_data)
      sprintf(fname, "%s/%s.entropy.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.entropy.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 

    printf("reading entropies from to %s\n", fname) ;
    mri_entropy = MRIread(fname) ;
    if (randomize_data)
      sprintf(fname, "%s/%s.posterior.rand.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    else
      sprintf(fname, "%s/%s.posterior.mgz", 
	      FileNamePath(mris->fname, path), 
	      mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
    printf("reading integrated posterior from %s\n", fname) ;
    mri_total_posterior = MRIread(fname) ;

  }
  mri_posterior_on_splines = compute_posterior_on_paths(mris, mri_splines, mri_total_posterior, mri_aseg,NULL) ;
  if (randomize_data)
    sprintf(fname, "%s/%s.spline_posterior.rand.mgz", 
	    FileNamePath(mris->fname, path), 
	    mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
  else
    sprintf(fname, "%s/%s.spline_posterior.mgz", 
	    FileNamePath(mris->fname, path), 
	    mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
  printf("writing posterior projected onto optimal splines to %s\n", fname) ;
  MRIwrite(mri_posterior_on_splines, fname) ;

  mri_filtered_posterior_on_spline = compute_filtered_posterior_on_paths(mris, mri_splines, mri_total_posterior, mri_aseg,NULL) ;
  if (randomize_data)
    sprintf(fname, "%s/%s.spline_posterior.filtered.rand.mgz", 
	    FileNamePath(mris->fname, path), 
	    mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
  else
    sprintf(fname, "%s/%s.spline_posterior.filtered.mgz", 
	    FileNamePath(mris->fname, path), 
	    mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh"); 
  printf("writing posterior projected and filtered onto optimal splines  to %s\n", fname) ;
  MRIwrite(mri_filtered_posterior_on_spline, fname) ;

  MRIfree(&mri_splines) ; MRIfree(&mri_posterior); VLSTfree(&vl_posterior) ; MRIfree(&mri_entropy) ;
  MRIfree(&mri_total_posterior) ;
  MRIfree(&mri_possible_migration_paths) ; MRIfree(&mri_path_grad) ;
  VectorFree(&v_intensity) ; VectorFree(&v_aseg) ; MatrixFree(&m_intensity2aseg) ; MRIfree(&mri_wm_dist) ;
  MRIfree(&mri_eroded_aseg) ;
  return(mri_pvals) ;
}

#define MAX_VENTRICLE_DIST 500
static VOXEL_LIST *
compute_path_to_ventricles(MRI_SURFACE *mris, int vno, MRI *mri_ventricle_dist_grad, MRI *mri_aseg)
{
  double     xv, yv, zv, step, dx, dy, dz, norm ;
  VOXEL_LIST *vl ;
  int        label, max_steps ;
  VERTEX     *v ;

  if (vno == Gdiag_no)
    DiagBreak() ;
  step = mri_ventricle_dist_grad->xsize/2 ;

  v = &mris->vertices[vno] ;
  max_steps = nint(ceil(MAX_VENTRICLE_DIST/step)) ;
  vl = VLSTalloc(max_steps+1) ; vl->nvox = 0 ;

  // once to compute length, and again to add voxels to list
  MRISsurfaceRASToVoxelCached(mris, mri_aseg, v->x, v->y, v->z, &xv, &yv, &zv) ;
  do
  {
    MRIsampleVolumeFrame(mri_ventricle_dist_grad, xv, yv, zv, 0, &dx) ;
    MRIsampleVolumeFrame(mri_ventricle_dist_grad, xv, yv, zv, 1, &dy) ;
    MRIsampleVolumeFrame(mri_ventricle_dist_grad, xv, yv, zv, 2, &dz) ;
    norm = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (!FEQUAL(norm,1))
      DiagBreak() ;
    if (FZERO(norm))
    {
      dy = randomNumber(0.0, 1.0) ;
      dy = randomNumber(0.0, 1.0) ;
      dz = randomNumber(0.0, 1.0) ;
      norm = sqrt(dx*dx + dy*dy + dz*dz) ;
    }
    dx /= norm ; dy /= norm ; dz /= norm ;
    xv -= step*dx ; yv -= step*dy ; zv -= step*dz ;
    label = MRIgetVoxVal(mri_aseg, nint(xv), nint(yv), nint(zv), 0) ;
    if (vl->nvox >= max_steps)
      ErrorReturn(NULL, (ERROR_BADPARM, "could not find path to ventricles for vno %d", vno)) ;
    VLSTadd(vl, nint(xv), nint(yv), nint(zv), xv, yv, zv);
  } while (IS_LAT_VENT(label) == 0) ;
  return(vl) ;
}

static int burnin = 1000 ;
static int jump = 5 ;

static VOXEL_LIST *
find_optimal_spline(VOXEL_LIST *vl, MRI *mri_intensity, MRI *mri_aseg, MRI *mri_wm_dist, double gm_mean, int nsamples, MRI_SURFACE *mris, int spline_control_points, int mcmc_samples, double spline_length_penalty, double spline_nonwm_penalty, double spline_interior_penalty, MRI *mri_posterior, VOXEL_LIST *vl_posterior,
		    double *pentropy, MRI *mri_total_posterior)
{
  VOXEL_LIST    *vl_spline, *vl_spline_optimal, *vl_spline_current ;
  double        energy, xn, yn, zn, best_energy, acceptance_val, current_energy, entropy ;
  int           n, cnum, label, nposterior = 0, lastn = 0 ;

  if (vl_posterior)
    vl_posterior->nvox = 0 ;  // reset for new list

  if (vl->nvox <= 2)
    ErrorReturn(NULL, (ERROR_NOFILE, "find_optimal_spline: input voxlist is too short (%d)", vl->nvox)) ;

  if (spline_control_points > vl->nvox)
    spline_control_points = vl->nvox ;
  vl_spline_optimal = VLSTsplineFit(vl, spline_control_points) ;
  best_energy = compute_spline_energy(vl_spline_optimal, mri_intensity, mri_aseg, mri_wm_dist, gm_mean,energy_flags,
				     spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty) ;
  if (DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    VOXEL_LIST *vl_interp = VLSTinterpolate(vl_spline_optimal, .1) ;

    sprintf(fname, "spline.000.label") ;
    VLSTwriteLabel(vl_interp, fname, mris, mri_intensity) ;
    sprintf(fname, "spline.000.cpts.label") ;
    VLSTwriteLabel(vl_spline_optimal, fname, mris, mri_intensity) ;
    VLSTfree(&vl_interp) ;
  }

  vl_spline_current = VLSTcopy(vl_spline_optimal, NULL, 0, vl_spline_optimal->nvox) ;
  current_energy = best_energy ;
  for (n = 0 ; n < nsamples ; n++)
  {
    cnum = (int)randomNumber(1.0, spline_control_points-.001) ;
    vl_spline = VLSTcopy(vl_spline_current, NULL, 0, vl_spline_current->nvox) ;
    xn = proposal_sigma * PDFgaussian() ;
    yn = proposal_sigma * PDFgaussian() ;
    zn = proposal_sigma * PDFgaussian() ;

    // project out tangential component
    if (cnum != spline_control_points-1)   // don't project out tangential component for end point
    {
      double tx, ty, tz, dot, norm ;
      int    km1, kp1 ;

      km1 = cnum-1 ; kp1 = cnum+1 ;
      if (cnum == 0)
	km1 = 0 ;
      else if (cnum == spline_control_points-1)
	kp1 = spline_control_points-1 ;
      tx = vl_spline->xd[kp1] - vl_spline->xd[km1] ;
      ty = vl_spline->yd[kp1] - vl_spline->yd[km1] ;
      tz = vl_spline->zd[kp1] - vl_spline->zd[km1] ;
      norm = sqrt(tx*tx + ty*ty + tz*tz) ;
      if (DZERO(norm))
	norm = 1.0 ;
      tx /= norm ; ty /= norm ;  tz /= norm ;
      dot = tx * xn + ty*yn + tz*zn ;
      // remove tangential component
      xn = xn - dot*tx ; yn = yn - dot*ty ; zn = zn - dot*tz ;
    }
    vl_spline->xd[cnum] += xn ; vl_spline->yd[cnum] += yn ; vl_spline->zd[cnum] += zn ;
    vl_spline->xi[cnum] = nint(vl_spline->xd[cnum]) ;
    vl_spline->yi[cnum] = nint(vl_spline->yd[cnum]) ;
    vl_spline->zi[cnum] = nint(vl_spline->zd[cnum]) ;
    if (cnum == spline_control_points-1)  // last one must be in ventricles
    {
      label = MRIgetVoxVal(mri_aseg, vl_spline->xi[cnum], vl_spline->yi[cnum], vl_spline->zi[cnum],0) ;
      if (IS_LAT_VENT(label) == 0)  // last control point must be in ventricles
      {
	VLSTfree(&vl_spline) ;
	continue ;
      }
    }
    energy = compute_spline_energy(vl_spline, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, energy_flags,
				     spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty) ;
    acceptance_val = exp((current_energy-energy)/acceptance_sigma) ;
    if (randomNumber(0.0, 1.0) < acceptance_val)
    {
      VLSTfree(&vl_spline_current) ;
      vl_spline_current = vl_spline ;
      current_energy = energy ;
      if (n > burnin && mri_posterior && n > (lastn+jump))
      {
	
	VOXEL_LIST *vl_interp = VLSTinterpolate(vl_spline_current, 1.0) ;
	VLSTinterpolateSplineIntoVolume(vl_spline_current, mri_total_posterior, 
					mri_posterior->xsize, vl_posterior, 1.0/vl_interp->nvox) ;
	VLSTinterpolateSplineIntoVolume(vl_spline_current, mri_posterior, mri_posterior->xsize, vl_posterior,
	  1.0/vl_interp->nvox) ;
	VLSTfree(&vl_interp) ;
	nposterior++ ;
	lastn = n ;
      }

      if (best_energy > energy)
      {
	best_energy = energy ;
	VLSTfree(&vl_spline_optimal) ;
	vl_spline_optimal = VLSTcopy(vl_spline_current, NULL, 0, vl_spline_current->nvox) ;
	if (DIAG_VERBOSE_ON)
	{
	  char fname[STRLEN] ;
	  VOXEL_LIST *vl_interp = VLSTinterpolate(vl_spline_optimal, 1) ;
	  
	  sprintf(fname, "spline.%3.3d.label", n+1) ;
	  printf("%d: new optimal energy %2.2f, perturbation %d: (%2.1f, %2.1f, %2.1f)\n", 
		 n, best_energy, cnum, xn, yn, zn) ;
	  VLSTwriteLabel(vl_interp, fname, mris, mri_intensity) ;
	  sprintf(fname, "spline.%3.3d.cpts.label", n+1) ;
	  VLSTwriteLabel(vl_spline_optimal, fname, mris, mri_intensity) ;
	  VLSTfree(&vl_interp) ;
	}
      }
    }
    else
      VLSTfree(&vl_spline) ;
  }

  entropy = VLSTcomputeEntropy(vl_posterior, mri_posterior, nposterior) ;
  *pentropy = entropy ;
  return(vl_spline_optimal) ;
}

static double
compute_spline_energy(VOXEL_LIST *vl, MRI *mri_intensity, MRI *mri_aseg, MRI *mri_wm_dist, double mean, int which, 
		      double spline_length_penalty, double spline_nonwm_penalty, double spline_interior_penalty) 
{
  int         n, label, num ;
  VOXEL_LIST *vl_interp ;
  double     val, wm_dist ;
  double     neg_log_p ;

  neg_log_p = 0 ;
  vl_interp = VLSTinterpolate(vl, 1) ;

  for (n = num = 0 ; n < vl_interp->nvox ; n++)
  {
    if (vl_interp->xi[n] == Gx &&  vl_interp->yi[n] == Gy &&  vl_interp->zi[n] == Gz)
      DiagBreak() ;
    MRIsampleVolume(mri_wm_dist, vl_interp->xd[n], vl_interp->yd[n], vl_interp->zd[n], &wm_dist) ;
    label = MRIgetVoxVal(mri_aseg, vl_interp->xi[n], vl_interp->yi[n], vl_interp->zi[n], 0) ;
    if ((which & SPLINE_WM_DIST) && (wm_dist > max_wm_dist))
      neg_log_p += spline_interior_penalty*(wm_dist-max_wm_dist);
    if (IS_WMH(label) == 0)
    {
      neg_log_p += spline_nonwm_penalty ;
      num++ ;
      continue ;
    }
    MRIsampleVolume(mri_intensity, vl_interp->xd[n], vl_interp->yd[n], vl_interp->zd[n], &val) ;
#if 0
    if (which & SPLINE_ABOVE && val > mean)// don't compute intensities of non-wm tissue. Length penalty will handle this
    {
      neg_log_p += spline_nonwm_penalty ;
      continue ;
    }
#endif
    num++ ;
    if (val < mean)  // value is below the mean
    {
      if (which & SPLINE_BELOW)
	neg_log_p += SQR(val-mean) ;
      if (which & SPLINE_SIGNED)
	neg_log_p -= SQR(val-mean) ;
    }
    else   // value is above the mean
    {
      if (which & SPLINE_ABOVE)
	neg_log_p += (val-mean) ;
    }
  }

  // add in length penalty
  if (num == 0)
    num = 1 ;
  if (neg_log_p < 0)
    neg_log_p = 0 ;
  if (which & SPLINE_LENGTH)
    neg_log_p = spline_length_penalty*sqrt(vl_interp->nvox) + sqrt(neg_log_p/num) ;
  else
    neg_log_p= sqrt(neg_log_p/num) ;

  VLSTfree(&vl_interp) ;
  return(neg_log_p) ;
}

static MRI *
compute_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg,
			   MRI *mri_posterior_on_spline)
{
  int         vno ;
  VOXEL_LIST  *vl_spline ;
  double      mean ;
  MRI         *mri_mask, *mri_masked_posterior ;

  if (0)
  {
    double sigma = 2.0 ;
    MRI *mri_kernel = MRIgaussian1d(sigma, 11) ;
    mri_posterior_on_spline = MRIlaplacian(mri_total_posterior, NULL) ;
    MRIscalarMul(mri_posterior_on_spline, mri_posterior_on_spline, -1) ;
    MRIwrite(mri_posterior_on_spline, "lap.mgz") ;
    MRIconvolveGaussian(mri_posterior_on_spline, mri_posterior_on_spline, mri_kernel) ;
    MRIwrite(mri_posterior_on_spline, "slap.mgz") ;
    MRIcopy(mri_posterior_on_spline, mri_total_posterior) ;
    MRIfree(&mri_kernel) ;
    MRIfree(&mri_posterior_on_spline) ;
  }

  if (mri_posterior_on_spline == NULL)
    mri_posterior_on_spline = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;

  mri_mask = MRIcloneDifferentType(mri_total_posterior, MRI_UCHAR) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, WM_hypointensities) ;
  MRIbinarize(mri_mask, mri_mask, 1, 0, 1) ;
  mri_masked_posterior = MRImask(mri_total_posterior, mri_mask, NULL, 0, 0) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;

    vl_spline = VLSTfromMRI(mri_splines, vno) ;

    mean = VLSTcomputeSplineMedian(vl_spline, mri_masked_posterior, mri_total_posterior->xsize/2) ;
    mean = VLSTcomputeSplineSegmentMean(vl_spline, mri_masked_posterior, mri_total_posterior->xsize/2,
					8, 12) ;
    VLSTfree(&vl_spline) ;
    MRIsetVoxVal(mri_posterior_on_spline, vno, 0, 0, 0, mean) ;
  }

  MRIfree(&mri_masked_posterior) ; MRIfree(&mri_mask) ;

  return(mri_posterior_on_spline) ;
}

#define ANGLE_STEP RADIANS(10)

static MRI *
compute_filtered_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg, MRI *mri_filtered_posterior_on_spline) 
{
  int         vno, nm1, np1, nsamples, n, nspline, outside_bad ;
  VOXEL_LIST  *vl_spline, *vl ;
  double      tangent[3], normal1[3], normal2[3], norm, theta, radius = 4.0, p[3], x1[3], y1[3],
    xscale, yscale, outside, total, val ;
  MRI         *mri_mask, *mri_masked_posterior ;

  if (mri_filtered_posterior_on_spline == NULL)
    mri_filtered_posterior_on_spline = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT) ;

  mri_mask = MRIcloneDifferentType(mri_total_posterior, MRI_UCHAR) ;
  MRIcopyLabel(mri_aseg, mri_mask, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, Right_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_mask, WM_hypointensities) ;
  MRIbinarize(mri_mask, mri_mask, 1, 0, 1) ;
  mri_masked_posterior = MRImask(mri_total_posterior, mri_mask, NULL, 0, 0) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;

    vl_spline = VLSTfromMRI(mri_splines, vno) ;
    vl = VLSTinterpolate(vl_spline, mri_total_posterior->xsize) ;

    for (total = 0.0, nspline = n = 0 ; n < vl->nvox ; n++)
    {
      if (MRIgetVoxVal(mri_mask, nint(vl->xi[n]), nint(vl->yi[n]), nint(vl->zi[n]), 0) == 0)
	continue ; // not in wm
      MRIsampleVolume(mri_total_posterior, vl->xd[n], vl->yd[n], vl->zd[n], &val) ;
      nm1 = n-1 ; np1 = n+1 ;
      if (n == 0)
	nm1 = n ;
      else if (n == vl->nvox-1)
	np1 = n ;

      tangent[0] = vl->xd[np1]-vl->xd[nm1] ;
      tangent[1] = vl->yd[np1]-vl->yd[nm1] ;
      tangent[2] = vl->zd[np1]-vl->zd[nm1] ;
      norm = VLEN(tangent) ; tangent[0] /= norm ; tangent[1] /= norm ; tangent[2] /= norm ;
      // pick some random, no colinear vector by permuting tangent
      if (!FEQUAL(tangent[0], tangent[1]))
      {
	normal1[0] = tangent[1] ;
	normal1[1] = tangent[0] ;
	normal1[2] = tangent[2] ;
      }
      else  // if [0] and [1] are equal, 1 and 2 can't be for a unit vector
      {
	normal1[0] = tangent[0] ;
	normal1[1] = tangent[2] ;
	normal1[2] = tangent[1] ;
      }
      CROSS(normal2, tangent, normal1) ;
      norm = VLEN(normal2) ; normal2[0] /= norm ; normal2[1] /= norm ; normal2[2] /= norm ;
      CROSS(normal1, tangent, normal2) ;
      norm = VLEN(normal1) ; normal1[0] /= norm ; normal1[1] /= norm ; normal1[2] /= norm ;

      outside_bad = 0 ;
      for (nsamples = 0, outside = theta = 0.0 ; theta <= 2*M_PI ; theta += ANGLE_STEP)
      {
	xscale = radius * cos(theta) ;
	yscale = radius * sin(theta) ;
	SCALAR_MUL(x1, xscale, normal1) ;
	SCALAR_MUL(y1, yscale, normal2) ;
	ADD(p, x1, y1) ;
	if (MRIgetVoxVal(mri_mask, nint(vl->xd[n]+p[0]), nint(vl->yd[n]+p[1]), nint(vl->zd[n]+p[2]), 0) == 0)
	{
	  outside_bad = 1 ;
	  break ; // not in wm
	}
	MRIsampleVolume(mri_total_posterior, vl->xd[n]+p[0], vl->yd[n]+p[1], vl->zd[n]+p[2], &val) ;
	nsamples++ ;
	outside += val ;
      }
      if (outside_bad)
	continue ;
      nspline++ ;
      outside /= (double)nsamples ;
      MRIsampleVolume(mri_total_posterior, vl->xd[n], vl->yd[n], vl->zd[n], &val) ;
      total += val / outside ;
    }

    VLSTfree(&vl_spline) ;
    VLSTfree(&vl) ;
    if (nspline > 0)
      MRIsetVoxVal(mri_filtered_posterior_on_spline, vno, 0, 0, 0, total/(double)nspline) ;
  }

  MRIfree(&mri_masked_posterior) ; MRIfree(&mri_mask) ;

  return(mri_filtered_posterior_on_spline) ;
}
