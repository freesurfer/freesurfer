/**
 * @brief MCMC for computing posterior of splines connecting cortical parcellation with itself
 *
 * Fit a Catmull Rom spline to each pair of points in the cortex, initializing it with a connection along the
 * shortest interior path  between them, then use MCMC to estimate the posterior distribution.
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
#include "cmat.h"
#include "fsinit.h"

#include "romp_support.h"



#define MIN_SPLINE_CONTROL_POINTS  3

#define SPLINE_WM_DIST    0x0001
#define SPLINE_LENGTH     0x0008
#define SPLINE_SIGNED     0x0010
#define MIN_LH_CORTEX 1000
#define MAX_LH_CORTEX 1999
#define MIN_RH_CORTEX 2000
#define MAX_RH_CORTEX 2999
#define MIN_CORTEX    MIN_LH_CORTEX
#define MAX_CORTEX    MAX_RH_CORTEX
#define MAX_LABELS  (MAX_CORTEX-MIN_CORTEX+1)
#define XHEMI_OFFSET  1000

#define LAPLACE_SOURCE  -1
#define LAPLACE_TARGET  1
#define LAPLACE_OUTSIDE 2
#if 0
static VOXEL_LIST *compute_path_to_ventricles(MRI_SURFACE *mris, int vno, MRI *mri_ventricle_dist_grad, MRI *mri_aseg) ;
static MRI *compute_migration_probabilities(MRI_SURFACE *mris, MRI *mri, MRI *mri_aseg, TRANSFORM *xform, MRI *mri_pvals, int spline_control_points, int mcmc_samples, int read_flag) ;
static MRI *compute_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg, MRI *mri_posterior_on_spline) ;

static MRI *compute_filtered_posterior_on_paths(MRI_SURFACE *mris, MRI *mri_splines, MRI *mri_total_posterior, MRI *mri_aseg, MRI *mri_posterior_on_spline) ;



static double spline_nonwm_penalty = 200 ;
static double max_wm_dist = -2.5 ;


#endif

//static VOXEL_LIST *MRIfindBestSpline(MRI *mri_aseg, MRI *mri_wm_dist, int label1_target, int label2_target, int ncontrol) ;

static int use_laplace = 0 ;
static int hemi = 0 ;
static int label1_target = -1 ;
static int label2_target = -1 ;
static int ncontrol = 5 ;

static int energy_flags = SPLINE_WM_DIST | SPLINE_LENGTH ;
static double spline_length_penalty = 2 ;
static double proposal_sigma = .5 ; // stddev of noise distribution
static double acceptance_sigma = 6 ;
static double compute_spline_energy(VOXEL_LIST *vl, MRI *mri_wm_dist, int flags, double spline_length_penalty, double spline_interior_penalty) ;
static VOXEL_LIST *find_optimal_spline(VOXEL_LIST *vl, MRI *mri_aseg, MRI *mri_wm_dist, int mcmc_samples, double spline_length_penalty, double spline_interior_penalty) ;

static MRI *MRIinteriorDistanceTransform(MRI *mri_aseg, MRI *mri_wm_interior, MRI *mri_paths, int label, int hemi) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static VOXEL_LIST *compute_spline_initialization(MRI *mri_aseg, MRI *mri_wm, MRI *mri_wm_dist, MRI *mri_label1_dist, MRI *mri_dist_grad, int label1, int label2, int min_spline_control_points) ;
static VOXEL_LIST *compute_path_to_label(MRI *mri_aseg, MRI *mri_wm, MRI *mri_dist_grad, MRI *mri_dist, int label1, int label2) ;

static int min_spline_control_points = MIN_SPLINE_CONTROL_POINTS ;
static int write_diags = 0 ;
static int randomize_data = 0 ;
static int mcmc_samples = 10000 ;
static double noise = .1 ;

static double spline_interior_penalty = 100 ;

static int vol_thresh = 75 ;

extern VOXLIST *MRIcomputeLaplaceStreamline(MRI *mri_laplace, int max_steps, float x0, float y0, float z0,
					    float source_val,float target_val, float outside_val);


const char *Progname ;
static void usage_exit(int code) ;

static char sdir[STRLEN] = "" ;

static int xhemi = 0 ;  // if 1, only do homologous ROIs across the hemis

int
main(int argc, char *argv[]) {
  char         fname[STRLEN], *subject, base_name[STRLEN] ;
  int          nargs, x, y, z, labels[MAX_LABELS], label_found[MAX_LABELS], nlabels ;
  int          msec, minutes, seconds, label, label2 ;
  Timer start ;
  MRI          *mri_aseg, *mri_wm, *mri_label1_dist, *mri_dist_grad, *mri_smooth, *mri_wm_dist ;
  VOXEL_LIST   **vl_splines[MAX_LABELS], *vl ;
  CMAT         *cmat ;
  MRI          *mri_tmp, *mri_wm_only ;

  nargs = handleVersionOption(argc, argv, "mris_init_global_tractography");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  setRandomSeed(-1L) ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  FSinit() ;
  start.reset() ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;
  if (sdir[0] == 0)  // not specified on command line
  {
    char *cp = getenv("SUBJECTS_DIR") ;
    if (cp == NULL)
      ErrorExit(ERROR_UNSUPPORTED, "%s: SUBJECTS_DIR must be specified on the command line or in the env",
		Progname) ;
    strcpy(sdir, cp) ;
  }
  subject = argv[1] ;
  
  int req = snprintf(fname, STRLEN, "%s/%s/mri/%s.mgz", sdir, subject, argv[2]) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load aparc+aseg from %s", Progname, fname) ;

  FileNameRemoveExtension(argv[3], base_name) ;

  memset(label_found, 0, sizeof(labels)) ;
  for (nlabels = x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
	label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	if ((label >= MIN_CORTEX && label <= MAX_CORTEX) && (label_found[label-MIN_CORTEX] == 0))
	{
	  if (MRIvoxelsInLabel(mri_aseg, label) < vol_thresh && !xhemi)
	  {
	    label_found[label-MIN_CORTEX] = -1 ;
	    printf("ignoring small label %s (%d)\n", cma_label_to_name(label), label) ;
	    continue ;
	  }

	  label_found[label-MIN_CORTEX] = 1 ;
	  labels[nlabels++] = label ;
	}
      }
  
  printf("%d unique cortical labels found\n", nlabels) ;

  mri_wm = MRIclone(mri_aseg, NULL) ;
//  MRIcopyLabel(mri_aseg, mri_wm, Right_VentralDC) ;
//  MRIcopyLabel(mri_aseg, mri_wm, Left_VentralDC) ;
//  MRIcopyLabel(mri_aseg, mri_wm, Brain_Stem) ;
  MRIcopyLabel(mri_aseg, mri_wm, Left_Cerebellum_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_wm, Right_Cerebellum_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_wm, Left_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_wm, Right_Cerebral_White_Matter) ;
  MRIcopyLabel(mri_aseg, mri_wm, CC_Posterior) ;
  MRIcopyLabel(mri_aseg, mri_wm, CC_Mid_Posterior) ;
  MRIcopyLabel(mri_aseg, mri_wm, CC_Central) ;
  MRIcopyLabel(mri_aseg, mri_wm, CC_Mid_Anterior) ;
  MRIcopyLabel(mri_aseg, mri_wm, CC_Anterior) ;
  MRIcopyLabel(mri_aseg, mri_wm, WM_hypointensities) ;
  MRIcopyLabel(mri_aseg, mri_wm, Left_VentralDC) ;
  MRIcopyLabel(mri_aseg, mri_wm, Right_VentralDC) ;
  MRIcopyLabel(mri_aseg, mri_wm, Brain_Stem) ;
  if (label1_target > 0)  // operate in two-label mode
    MRIcopyLabel(mri_aseg, mri_wm, label1_target) ;
  if (label2_target > 0)  // operate in two-label mode
    MRIcopyLabel(mri_aseg, mri_wm, label2_target) ;
  MRIbinarize(mri_wm, mri_wm, 1, 0, 1) ;
  mri_wm_dist = MRIdistanceTransform(mri_wm, NULL, 1, 25, DTRANS_MODE_SIGNED, NULL);
  mri_wm_only = MRIcopy(mri_wm, NULL) ; // target label will be added to mri_wm later
//  labels[0] = 1024 ; 
//  labels[1] = 1035 ;
  if (label1_target > 0)  // operate in two-label mode
  {
    VOXEL_LIST *vl_spline ;

    if (use_laplace)
    {
      MRI        *mri_laplace ;
      int xm, ym, zm, xv, yv, zv, xk, yk, zk, xi, yi, zi ;
      double dist, min_dist ;
      VOXLIST *vl ;

      mri_laplace = MRIsolveLaplaceEquation(mri_wm_only, mri_aseg, label2_target,label1_target,LAPLACE_SOURCE,LAPLACE_TARGET,LAPLACE_OUTSIDE) ;
      
      min_dist = 1e10;
      xm = ym = zm = 0 ;
      for (xv = 0 ; xv < mri_aseg->width ; xv++)
	for (yv = 0 ; yv < mri_aseg->height ; yv++) 
	  for (zv = 0 ; zv < mri_aseg->depth ; zv++) 
	  {
	    if (MRIgetVoxVal(mri_aseg, xv, yv, zv,0) == label1_target)
	      DiagBreak() ;
	    if (MRIgetVoxVal(mri_aseg, xv, yv, zv,0) == label1_target && MRIlabelsInNbhd6(mri_wm,xv,yv,zv,1))
	    {
	      for (xk = -1 ; xk <= 1 ; xk++)
		for (yk = -1 ; yk <= 1 ; yk++)
		  for (zk = -1 ; zk <= 1 ; zk++)
		  {
		    if (abs(xk) + abs(yk) + abs(zk) != 1)
		      continue ;  // enforce 6-connectivity
		    xi = mri_laplace->xi[(xv+xk)] ; yi = mri_laplace->yi[(yv+yk)] ; zi = mri_laplace->zi[(zv+zk)] ;
		    if (MRIgetVoxVal(mri_wm,xi,yi,zi,0) == 0)
		      continue ;
		    dist = MRIgetVoxVal(mri_laplace, xi, yi, zi, 0) ;
		    if (dist < min_dist)
		    {
		      xm = xi ; ym = yi ; zm = zi ; min_dist = dist ;
		    }
		  }
	    }
	  }
      vl = MRIcomputeLaplaceStreamline(mri_laplace, 1500, xm, ym, zm, LAPLACE_SOURCE,LAPLACE_TARGET,LAPLACE_OUTSIDE) ;
      if (vl == NULL)
	ErrorExit(Gerror, "") ;
      
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
	MRIwrite(mri_laplace, "lap.mgz");
      printf("writing label to %s\n", argv[3]);
      VLSTwriteLabel(vl, argv[3], NULL, mri_aseg) ;
      exit(0) ;
    }

    mri_tmp = MRIclone(mri_aseg, NULL) ;
    MRIcopyLabel(mri_aseg, mri_tmp, label2_target) ;
    MRIbinarize(mri_tmp, mri_tmp, 1, 0, 2) ;
    MRIcopyLabel(mri_tmp, mri_wm, 2) ;  // this label will be a fixed point in the smoothing
    MRIfree(&mri_tmp) ;

    MRIcopy(mri_wm_only, mri_wm) ;
    mri_tmp = MRInbrThresholdLabel(mri_aseg, NULL, label2_target, 0, 1, 5) ;  // remove isolated voxels in label
    mri_label1_dist = MRIinteriorDistanceTransform(mri_tmp, mri_wm, NULL, label2_target, hemi) ;
    MRIfree(&mri_tmp) ;

    mri_smooth = MRIsmoothLabel6Connected(mri_label1_dist, mri_wm, NULL, 500, 1, 2, .5) ;
    MRIcopy(mri_wm_only, mri_wm) ;
    {
      int x, y, z, i, l ;
      float val ;
      
      mri_tmp = MRIcopy(mri_smooth, NULL) ;
      for (i = 0 ; i < 10 ; i++)
      {
#if 0
#ifdef HAVE_OPENMP
        val = 0 ; l = 0 ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(val, l) shared(mri_aseg, mri_tmp, labels, y, z) schedule(static,1)
#endif
#endif
	for (x = 0 ; x < mri_aseg->width ; x++) {
	  for (y = 0 ; y < mri_aseg->height ; y++) {
	    for (z = 0 ; z < mri_aseg->depth ; z++) {
	      l = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	      if (IS_WMH(l) || l == label2_target)
		continue ;
	      val = MRImaxInRegion(mri_smooth, x, y, z, 1)  + 1 ;  // one greater than max
	      MRIsetVoxVal(mri_tmp, x, y, z, 0, val) ;
	    }
	  }
	}
	MRIcopy(mri_tmp, mri_smooth) ;
      }
      MRIfree(&mri_tmp) ;
    }

    MRIfree(&mri_label1_dist) ; mri_label1_dist = mri_smooth ;
//    mri_label1_dist = mri_laplace ;
    mri_dist_grad = MRIsobel(mri_label1_dist, NULL, NULL) ;
    MRInormalizeSequence(mri_dist_grad, 1.0) ;
    vl_spline = compute_spline_initialization(mri_aseg, mri_wm, 
					      mri_wm_dist, mri_label1_dist,
					      mri_dist_grad,
					      label2_target, label1_target,
					      min_spline_control_points) ;
    if (vl_spline == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not find path between labels", Progname) ;

    printf("writing label to %s\n", argv[3]);
    VLSTwriteLabel(vl_spline, argv[3], NULL, mri_aseg) ;
    
    exit(0) ;
  }

  if (Gdiag_no > 0)
    nlabels = 3;
  for (label = 0 ; label < nlabels ; label++)
    vl_splines[label] = (VOXEL_LIST **)calloc(nlabels, sizeof(VOXEL_LIST *)) ;

  for (label = 0 ; label < nlabels-1 ; label++)
  {
    printf("************************* processing label %d of %d --> %s (%d) **************************\n", 
	   label+1, nlabels, cma_label_to_name(labels[label]), labels[label]) ;
    MRIcopy(mri_wm_only, mri_wm) ;
    if (label == 16)
      DiagBreak() ;
    mri_tmp = MRInbrThresholdLabel(mri_aseg, NULL, labels[label], 0, 1, 5) ;  // remove isolated voxels in label
    mri_label1_dist = MRIinteriorDistanceTransform(mri_tmp, mri_wm, NULL, labels[label], hemi) ;
    MRIfree(&mri_tmp) ;

    // temporarily add target label to wm so it affects smoothing
    mri_tmp = MRIclone(mri_aseg, NULL) ;
    MRIcopyLabel(mri_aseg, mri_tmp, labels[label]) ;
    MRIbinarize(mri_tmp, mri_tmp, 1, 0, 2) ;
    MRIcopyLabel(mri_tmp, mri_wm, 2) ;  // this label will be a fixed point in the smoothing
    MRIfree(&mri_tmp) ;
    mri_smooth = MRIsmoothLabel6Connected(mri_label1_dist, mri_wm, NULL, 500, 1, 2, .5) ;
    MRIcopy(mri_wm_only, mri_wm) ;
    {
      int x, y, z, i, l ;
      float val ;
      
      mri_tmp = MRIcopy(mri_smooth, NULL) ;
      for (i = 0 ; i < 10 ; i++)
      {
	ROMP_PF_begin
#if 0
#ifdef HAVE_OPENMP
	val = 0 ; l = 0 ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(val, l) shared(mri_aseg, mri_tmp, labels, y, z) schedule(static,1)
#endif
#endif
	for (x = 0 ; x < mri_aseg->width ; x++) {
          ROMP_PFLB_begin
          
	  for (y = 0 ; y < mri_aseg->height ; y++) {
	    for (z = 0 ; z < mri_aseg->depth ; z++)
	    {
	      l = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	      if (IS_WMH(l) || l == labels[label])
		continue ;
	      val = MRImaxInRegion(mri_smooth, x, y, z, 1)  + 1 ;  // one greater than max
	      MRIsetVoxVal(mri_tmp, x, y, z, 0, val) ;
	    }
	  }
          
          ROMP_PFLB_end
	}
        ROMP_PF_end

	MRIcopy(mri_tmp, mri_smooth) ;
      }
      MRIfree(&mri_tmp) ;
    }

    MRIfree(&mri_label1_dist) ; mri_label1_dist = mri_smooth ;
    mri_dist_grad = MRIsobel(mri_label1_dist, NULL, NULL) ;
    MRInormalizeSequence(mri_dist_grad, 1.0) ;
    if (0 && write_diags) {
      sprintf(fname, "dist.%d.mgz", labels[label]) ;
      MRIwrite(mri_label1_dist, fname) ;
      sprintf(fname, "grad.%d.mgz", labels[label]) ;
      MRIwrite(mri_dist_grad, fname) ;
    }

    for (label2 = 0 ; label2 < nlabels ; label2++)
    {
      if (xhemi &&   // only do corresponding ROIs
	  ((labels[label]+XHEMI_OFFSET != labels[label2]) &&
	   (labels[label2]+XHEMI_OFFSET != labels[label])))
	continue ;
      if (label2 == label || (xhemi == 0 && vl_splines[label2][label]))
	continue ;  // if this spline has already been successfully computed from the other direction
      vl = compute_spline_initialization(mri_aseg, mri_wm, mri_wm_dist, mri_label1_dist,
					 mri_dist_grad,
					 labels[label], labels[label2], 
					 min_spline_control_points) ;
      if (vl == NULL)
	continue ;
      if (label < label2)
	vl_splines[label][label2] = vl ;
      else
	vl_splines[label2][label] = vl ;

      if (write_diags)
      {
	sprintf(fname, "%s.%d.%d.label", base_name, labels[label], labels[label2]) ;
	printf("writing %d control point spline to %s for %s (%d) --> %s (%d)\n",
	       vl_splines[label][label2]->nvox, fname, cma_label_to_name(labels[label2]), labels[label2], 
	       cma_label_to_name(labels[label]), labels[label]);
	printf("\t(%d: [%d %d %d]) --> (%d: [%d %d %d])\n", 
	       labels[label2], vl->xi[0], vl->yi[0], vl->zi[0], 
	       labels[label], vl->xi[vl->nvox-1], vl->yi[vl->nvox-1], vl->zi[vl->nvox-1]) ; 
	VLSTwriteLabel(vl_splines[label][label2], fname, NULL, mri_aseg) ;
      }
      else
	printf("%d control point spline computed for %s (%d) --> %s (%d)\n",vl->nvox, 
	       cma_label_to_name(labels[label2]), labels[label2], cma_label_to_name(labels[label]), labels[label]);
      if (xhemi)
	break ;
    }
    MRIfree(&mri_label1_dist) ; MRIfree(&mri_dist_grad) ;
  }

  printf("writing initial spline estimates to %s\n", argv[3]) ;
  cmat = CMATalloc(nlabels, labels) ;
  for (label = 0 ; label < nlabels-1 ; label++)
    for (label2 = label+1 ; label2 < nlabels ; label2++)
    {
      if (vl_splines[label][label2] == NULL)
	continue ;
      cmat->splines[label][label2] = VLSTtoLabel(vl_splines[label][label2],  NULL, mri_aseg) ;
      cmat->splines[label2][label] = VLSTtoLabel(vl_splines[label][label2],  NULL, mri_aseg) ;
    }

  CMATwrite(cmat, argv[3]) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "minimal WM path-length calculation took %d hours, %d minutes and %d seconds.\n", minutes/60, minutes%60, seconds) ;
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
  else if (!stricmp(option, "LAPLACE"))
  {
    use_laplace = 1 ;
    printf("using Laplace equation to compute initialization\n") ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "LH"))
  {
    hemi = LEFT_HEMISPHERE ;
    printf("assuming label is left hemisphere\n") ;
  }
  else if (!stricmp(option, "RH"))
  {
    hemi = RIGHT_HEMISPHERE ;
    printf("assuming label is right hemisphere\n") ;
  }
  else if (!stricmp(option, "XHEMI"))
  {
    xhemi = 1 ;
    printf("only computing inter-hemispheric splines\n") ;
  }
  else if (!stricmp(option, "NCONTROL"))
  {
    ncontrol = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d control points in spline fit\n", ncontrol) ;
  }
  else if (!stricmp(option, "openmp")) {
    char str[STRLEN] ;
    sprintf(str, "OMP_NUM_THREADS=%d", atoi(argv[2]));
    putenv(str) ;
#ifdef HAVE_OPENMP
    omp_set_num_threads(atoi(argv[2]));
#else
    fprintf(stderr, "Warning - built without openmp support\n");
#endif
    nargs = 1 ;
    fprintf(stderr, "Setting %s\n", str) ;
  }
  else if (!stricmp(option, "LABELS"))
  {
    label1_target = atof(argv[2]) ;
    label2_target = atof(argv[3]) ;
    nargs = 2 ;
    printf("computing optimal initial spline connecting %s (%d) and %s (%d)\n",
	   cma_label_to_name(label1_target), label1_target,
	   cma_label_to_name(label2_target), label2_target) ;
  }
  else switch (toupper(*option)) {
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
  printf("usage: %s [options] <subject> <parcellation> <output volume>\n", 
         Progname) ;
  printf(
    "\t\n"
  );
  exit(code) ;
}

#define MAX_THICKNESS 5

static MRI *
MRIinteriorDistanceTransform(MRI *mri_aseg, MRI *mri_wm_interior, MRI *mri_paths, int target_label, int hemi)
{
  int         x, y,z, nadded, n, max_thick_vox, i, xi, yi, zi, xk, yk, zk, label, nbr_label ;
  VOXEL_LIST  *vl_current, *vl_next ;
  MRI         *mri_dilated = NULL ;
  float        fmin, fmax ;

  if (mri_paths == NULL)
    mri_paths = MRIcloneDifferentType(mri_aseg, MRI_FLOAT);

  MRIcopyLabel(mri_aseg, mri_paths, target_label) ;
  MRIbinarize(mri_paths, mri_paths, 1, 0, 1) ;

  n = 1 ;
  vl_next = VLSTalloc(mri_aseg->depth*mri_aseg->width*mri_aseg->height) ;
  vl_current = VLSTalloc(mri_aseg->depth*mri_aseg->width*mri_aseg->height) ;
  vl_next->nvox = vl_current->nvox = 0 ;
  VLSTcreate(mri_paths, 1, 1, vl_current, 0, 0) ;
  do
  {
    nadded = vl_next->nvox = 0 ;
    for (i = 0 ; i < vl_current->nvox ; i++)
    {
      x = vl_current->xi[i] ; y = vl_current->yi[i] ; z = vl_current->zi[i] ;
      label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
      for (xk = -1 ; xk <= 1 ; xk++)
	for (yk = -1 ; yk <= 1 ; yk++)
	  for (zk = -1 ; zk <= 1 ; zk++)
	  {
	    if (abs(xk) + fabs(yk) + abs(zk) != 1)
	      continue ;
	    xi = mri_aseg->xi[xk+x] ; yi = mri_aseg->yi[yk+y] ; zi = mri_aseg->zi[zk+z] ;
	    if (MRIgetVoxVal(mri_wm_interior, xi, yi, zi, 0) == 0)
	      continue ;
	    nbr_label = MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
	    if (IS_VENTRAL_DC(label) && (label != nbr_label))
	      continue ;   // no crossing in DC
	    if (MRIgetVoxVal(mri_paths, xi, yi, zi, 0))        // already in the current set
	      continue ;
	    if (label >= MIN_CORTEX)  //  make first step from gm into wm of the correct hemi
	    {
	      // user might use a label that needs to be explicitly identified
	      if ((label >= MIN_LH_CORTEX && label <= MAX_LH_CORTEX) ||
		  ((label == target_label) && (hemi == LEFT_HEMISPHERE))) // in the left hemi
	      {
		if (MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) == Right_Cerebral_White_Matter)
		  continue ;
	      }
	      else   // in right hemi gray matter
	      {
		if (MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) == Left_Cerebral_White_Matter)
		  continue ;
	      }
	    }
	    if (xi == Gx && yi == Gy && zi == Gz)
	      DiagBreak() ;
	    MRIsetVoxVal(mri_paths, xi, yi, zi, 0, n+1) ;
	    VLSTadd(vl_next, xi, yi, zi, xi, yi, zi) ;
	    nadded++ ;
	  }
	}
//    printf("round %d: %d voxels added\n", n, nadded) ;
    n++ ;
    VLSTcopy(vl_next, vl_current, 0, vl_next->nvox) ;
  } while (nadded > 0); 

  
  // extend distance transform into gray matter
  n = 0 ;  max_thick_vox = nint(MAX_THICKNESS / mri_aseg->xsize) ;
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
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  MRIsetVoxVal(mri_paths, x, y, z, 0, MRIgetVoxVal(mri_dilated, x, y, z, 0)+1) ;
	  nadded++ ;
	}
//    printf("round %d: %d voxels added\n", n, nadded) ;
    n++ ;
  } while (nadded > 0 && n < max_thick_vox); 
  
  
  // replace all exterior zeroes with valuse one greater than the max adjacent value
  MRInonzeroValRange(mri_paths, &fmin, &fmax) ;
  VLSTcreate(mri_paths, fmin, fmax, vl_current, 0, 1) ;  // create from only the border voxels
  n = 0 ;
  do
  {
    int val ;

    nadded = vl_next->nvox = 0 ;
    for (i = 0 ; i < vl_current->nvox ; i++)
    {
      x = vl_current->xi[i] ; y = vl_current->yi[i] ; z = vl_current->zi[i] ;
      val = MRImaxInRegion(mri_paths, x, y, z, 1)  + 1 ;  // one greater than max
      for (xk = -1 ; xk <= 1 ; xk++)
	for (yk = -1 ; yk <= 1 ; yk++)
	  for (zk = -1 ; zk <= 1 ; zk++)
	  {
	    if (abs(xk) + fabs(yk) + abs(zk) != 1)
	      continue ;
	    xi = mri_aseg->xi[xk+x] ; yi = mri_aseg->yi[yk+y] ; zi = mri_aseg->zi[zk+z] ;
	    if (MRIgetVoxVal(mri_paths, xi, yi, zi, 0) != 0)
	      continue ;
	    if (MRIgetVoxVal(mri_paths, xi, yi, zi, 0))        // already in the current set
	      continue ;
	    if (xi == Gx && yi == Gy && zi == Gz)
	      DiagBreak() ;
	    MRIsetVoxVal(mri_paths, xi, yi, zi, 0, val) ;
	    VLSTadd(vl_next, xi, yi, zi, xi, yi, zi) ;
	    nadded++ ;
	  }
    }
      
//    printf("round %d: %d voxels added\n", n, nadded) ;
    n++ ;
    VLSTcopy(vl_next, vl_current, 0, vl_next->nvox) ;
  } while (nadded > 0 && n < 20) ;


  VLSTfree(&vl_next) ; VLSTfree(&vl_current) ; MRIfree(&mri_dilated) ;
  return(mri_paths) ;
}

#if 0
static MRI *
compute_migration_probabilities(MRI_SURFACE *mris, MRI *mri_intensity, MRI *mri_aseg,
				TRANSFORM *xform, MRI *mri_pvals, int min_spline_control_points, int mcmc_samples, int read_flag)
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
    mri_splines = MRIalloc(mris->nvertices, min_spline_control_points, 3, MRI_FLOAT) ;
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
  v = NULL ; n = 0 ; vl = vl_spline = NULL ;
#pragma omp parallel for if_ROMP(experimental) firstprivate(n, v, vl, vl_spline, entropy, gm_mean, mcmc_samples) shared(mri_intensity, mri_aseg, mri_path_grad,mri_splines) schedule(static,1)
#endif
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    ROMP_PFLB_begin
    if ((vno % (mris->nvertices/200)) == 0)
      printf("processed %d of %d: %2.1f%%\n", vno, mris->nvertices, 100.0*vno/(float)mris->nvertices) ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      ROMP_PFLB_continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (read_flag == 0)
    {
      vl = compute_path_to_ventricles(mris, vno, mri_path_grad, mri_aseg) ;
      if (vl == NULL)
	ROMP_PFLB_continue ;
      vl_spline = find_optimal_spline(vl, mri_intensity, mri_aseg, mri_wm_dist, gm_mean, mcmc_samples,mris,min_spline_control_points, mcmc_samples, spline_length_penalty, spline_nonwm_penalty, spline_interior_penalty,
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
	
	vl_init_spline = VLSTsplineFit(vl, min_spline_control_points) ;
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


    for (n = 0 ; n < min_spline_control_points ; n++)
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
#endif

#define ALLOWABLE_SPLINE_DIST 1.0

static VOXEL_LIST *compute_optimal_number_of_control_points(MRI *mri_aseg, VOXEL_LIST *vl, int min_points, double allowable_dist) ;

static VOXEL_LIST *
compute_spline_initialization(MRI *mri_aseg, MRI *mri_wm, MRI *mri_wm_dist, MRI *mri_label1_dist, MRI *mri_dist_grad, int label1, int label2, int min_spline_control_points)
{
  VOXEL_LIST  *vl_spline, *vl, *vl_optimal, *vl_init_spline, *vl_optimal_final, *vl_interp ;

  if (label1 == 2010 && label2 == 1008)
    DiagBreak() ;
  vl = compute_path_to_label(mri_aseg, mri_wm, mri_dist_grad, mri_label1_dist, label1, label2) ;
  if (vl == NULL || vl->nvox < MIN_SPLINE_CONTROL_POINTS)
  {
    if (vl)
      VLSTfree(&vl) ;
    return(NULL) ;
  }
  vl_spline = compute_optimal_number_of_control_points(mri_aseg, vl, min_spline_control_points, ALLOWABLE_SPLINE_DIST) ;

  // use more control points to give spline flexibility to minimize energy functional
  vl_interp = VLSTinterpolate(vl_spline, 1) ;
  vl_init_spline = VLSTsplineFit(vl_interp, MIN(2*vl_spline->nvox, vl_interp->nvox)) ;
  vl_optimal = find_optimal_spline(vl_init_spline, mri_aseg, mri_wm_dist, mcmc_samples, 
				   spline_length_penalty, spline_interior_penalty);
  VLSTfree(&vl_interp) ;
  vl_interp = VLSTinterpolate(vl_optimal, 1) ;
  vl_optimal_final = compute_optimal_number_of_control_points(mri_aseg, vl_interp, min_spline_control_points, ALLOWABLE_SPLINE_DIST) ;


  if (write_diags ||  (label1 == 2010 && label2 == 1008)) // YYYY
  {
    printf("writing v[123].label\n") ;
    VLSTwriteLabel(vl, "vl1.label", 0L, mri_aseg)  ;
    VLSTwriteLabel(vl_spline, "vl2.label", 0L, mri_aseg)  ;
    VLSTwriteLabel(vl_optimal, "vl3.label", 0L, mri_aseg)  ;
    VLSTwriteLabel(vl_optimal_final, "vl4.label", 0L, mri_aseg)  ;
    DiagBreak() ;
  }

  VLSTfree(&vl_optimal) ; VLSTfree(&vl_init_spline) ; VLSTfree(&vl_spline) ;
  VLSTfree(&vl) ;VLSTfree(&vl_interp) ;

  return(vl_optimal_final) ;
}


//static double momentum = .5 ;
/*
  starting for the centroid of label2, march down gradient of label1 distance transform, staying
  in the interior of the WM, until we reach label1
*/
#define MAX_PARCEL_DIST 200
static VOXEL_LIST *
compute_path_to_label(MRI *mri_aseg, MRI *mri_wm, MRI *mri_dist_grad, MRI *mri_dist, int label1, int label2)
{
  double     step, xc, yc, zc, dist, min_dist ;
  VOXEL_LIST *vl ;
  int        label = 0, max_steps, nvox, xk, yk, zk, xi, yi, zi, xm, ym, zm, odx, ody, odz, dx, dy, dz, xv, yv, zv, found,
    backtrack_steps = 10, current_label, restarted = 0 ;
  MRI        *mri_path ;
  char        fname[STRLEN] ;

  mri_path = MRIclone(mri_aseg, NULL) ;

  step = mri_dist_grad->xsize/2 ;

  max_steps = nint(ceil(MAX_PARCEL_DIST/step)) ;
  vl = VLSTalloc(max_steps+1) ; vl->nvox = 0 ;

  // compute voxel in label2 that is closest to centroid
  for (xc = yc = zc = 0.0, nvox = xv = 0 ; xv < mri_aseg->width ; xv++)
    for (yv = 0 ; yv < mri_aseg->height ; yv++) 
      for (zv = 0 ; zv < mri_aseg->depth ; zv++) 
      {
	if (MRIgetVoxVal(mri_aseg, xv, yv, zv,0) == label2)
	{
	  xc += xv ; yc += yv ; zc += zv ; nvox++ ;
	}
      }
  xc /= nvox; yc /= nvox ; zc /= nvox ; 
//  printf("centroid found at (%2.0f, %2.0f, %2.0f)\n", xc, yc, zc) ;

  // find point in label that is closest to target label and adjacent to wm
  // find point in label that is closest to centroid and adjacent to wm
  min_dist = mri_aseg->width + mri_aseg->height + mri_aseg->depth ;
  xm = ym = zm = 0 ;
  for (xv = 0 ; xv < mri_aseg->width ; xv++)
    for (yv = 0 ; yv < mri_aseg->height ; yv++) 
      for (zv = 0 ; zv < mri_aseg->depth ; zv++) 
      {
	if (MRIgetVoxVal(mri_aseg, xv, yv, zv,0) == label2 && MRIlabelsInNbhd6(mri_wm,xv,yv,zv,1))
	{

	  for (xk = -1 ; xk <= 1 ; xk++)
	    for (yk = -1 ; yk <= 1 ; yk++)
	      for (zk = -1 ; zk <= 1 ; zk++)
	      {
		if (abs(xk) + abs(yk) + abs(zk) != 1)
		  continue ;  // enforce 6-connectivity
		xi = mri_dist->xi[(xv+xk)] ; yi = mri_dist->yi[(yv+yk)] ; zi = mri_dist->zi[(zv+zk)] ;
		if (MRIgetVoxVal(mri_wm,xi,yi,zi,0) == 0)
		  continue ;
		dist = MRIgetVoxVal(mri_dist, xi, yi, zi, 0) ;
		if (dist < min_dist)
		{
		  xm = xi ; ym = yi ; zm = zi ; min_dist = dist ;
		}
	      }
	}
      }
//  printf("initial spline anchor set to (%d, %d, %d)\n", xm, ym, zm) ;

  // move downwards along distance transform direction moving through appropriate labels until we get to dest
  xv = xm ; yv = ym ; zv = zm ;
  MRIsetVoxVal(mri_path, xv, yv, zv, 0, vl->nvox+1) ;
  dx = dy = dz = odx = ody = odz = 0 ;
  VLSTaddWithValue(vl, nint(xv), nint(yv), nint(zv), xv, yv, zv, MRIgetVoxVal(mri_dist,xv,yv,zv,0), 
		   MRIgetVoxVal(mri_aseg,xv,yv,zv,0)); 
  do
  {
    if (vl->nvox == Gdiag_no)
      DiagBreak() ;

    // examine all 6-connected voxels except where we came from to see what direction to move in
    min_dist = mri_aseg->width + mri_aseg->height + mri_aseg->depth ;
    current_label = MRIgetVoxVal(mri_aseg, xv, yv, zv, 0) ;
    for (found = 0, xk = -1 ; xk <= 1 ; xk++)
      for (yk = -1 ; yk <= 1 ; yk++)
	for (zk = -1 ; zk <= 1 ; zk++)
	{
	  if (abs(xk) + abs(yk) + abs(zk) != 1)
	    continue ;  // enforce 6-connectivity
	  xi = mri_dist->xi[(xv+xk)] ;
	  yi = mri_dist->yi[(yv+yk)] ;
	  zi = mri_dist->zi[(zv+zk)] ;
	  dist = MRIgetVoxVal(mri_dist, xi, yi, zi, 0) ;
	  label = MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
	  // find the smallest distance that we haven't visited with a feasible label
	  if ( (label == label1 && current_label != label2) ||     // found target label after leaving source label
	       ((dist < min_dist) &&                               // or find min_dist voxel that hasn't been visited
		IS_WMH(label) &&                                   // if we are moving into wm
		nint(MRIgetVoxVal(mri_path, xi, yi, zi, 0)) == 0)) // and we have never visited this voxel before
	  {
	    dx = xi-xv ;  dy = yi-yv ; dz = zi-zv ; 
	    found = 1 ;
	    min_dist = dist ;
	  }
	}
#define PATH_START 10
    if (found == 0)
    {
      if (vl->nvox-1 > (backtrack_steps-PATH_START))   // near the start of the path, start somewhere random
	vl->nvox = MAX(1, vl->nvox-backtrack_steps) ;
      else
	vl->nvox = vl->nvox-(int)randomNumber(1, (vl->nvox-1)) ;
      printf("backing up %d steps to %d and restarting minimization\n", backtrack_steps, vl->nvox) ;
      xv = vl->xi[vl->nvox-1] ; yv = vl->yi[vl->nvox-1] ; zv = vl->zi[vl->nvox-1] ;
      if (vl->nvox <= 1)
      {
	if (restarted)
	  MRIclear(mri_path) ;
	restarted++ ;
	if (restarted > 3)
	  return(NULL) ;    // fail
      }
      continue ;
      
//      printf("!!!!!!!!!!!!!!  path to label %d from label %d step could not be found !!!!!!!!!!!\n", label1, label2) ;
      if (write_diags)
      {
	MRIwrite(mri_path, "p.mgz") ;
	MRIwrite(mri_dist, "d.mgz") ;
	VLSTwriteLabel(vl, "vl.label", NULL, mri_aseg) ;
      }
      
      printf("vl: (%d, %d, %d) --> (%d, %d, %d)\n",
	     vl->xi[0], vl->yi[0], vl->zi[0], vl->xi[vl->nvox-1], vl->yi[vl->nvox-1], vl->zi[vl->nvox-1]) ;
      DiagBreak() ;   // should never happen
    }
    xv += dx ; yv += dy ; zv += dz ;
    if (MRIindexNotInVolume(mri_aseg, xv, yv, zv) || (vl->nvox >= max_steps))  // zzz
    {
      printf("!!!!!!!!!!!!!!  path to label %d from label %d left volume or too long !!!!!!!!!!!\n", label1, label2) ;
      if (write_diags)
      {
	MRIwrite(mri_path, "p.mgz") ;
	MRIwrite(mri_dist, "d.mgz") ;
	VLSTwriteLabel(vl, "vl.label", NULL, mri_aseg) ;
      }
      
      printf("vl: (%d, %d, %d) --> (%d, %d, %d)\n",
	     vl->xi[0], vl->yi[0], vl->zi[0], vl->xi[vl->nvox-1], vl->yi[vl->nvox-1], vl->zi[vl->nvox-1]) ;
      DiagBreak() ;
      backtrack_steps += 10 ;
      vl->nvox = 0 ;
      MRIsetValues(mri_path, 0) ;
      printf("restarting  path finding with backtracking length = %d\n", backtrack_steps) ;
      xv = xm ; yv = ym ; zv = zm ; odx = ody = odz = 0 ;
      continue ;
      VLSTfree(&vl) ; MRIfree(&mri_path) ;
      ErrorReturn(NULL, (ERROR_BADPARM, "!!!!!!!!!!!!!! path to label %d from label %d left volume !!!!!!!!!!!", label1, label2)) ;
    }
    odx = dx ; ody = dy ; odz = dz ;
    VLSTaddWithValue(vl, nint(xv), nint(yv), nint(zv), xv, yv, zv, MRIgetVoxVal(mri_dist,xv,yv,zv,0), 
		     MRIgetVoxVal(mri_aseg,xv,yv,zv,0)); 
    MRIsetVoxVal(mri_path, xv, yv, zv, 0, vl->nvox+1) ;
    label = MRIgetVoxVal(mri_aseg, xv, yv, zv, 0) ;
  } while (label != label1) ;

  if (write_diags)   // XXXX
  {
    sprintf(fname, "path.%d.%d.mgz", label1, label2) ;
    printf("writing path to %s\n", fname) ;
    MRIwrite(mri_path, fname) ;
    sprintf(fname, "path.%d.%d.label", label1, label2) ;
    printf("writing path to %s\n", fname) ;
    VLSTwriteLabel(vl, fname, NULL, mri_aseg) ;
  }
  MRIfree(&mri_path) ;
  return(vl) ;
}

static VOXEL_LIST *
compute_optimal_number_of_control_points(MRI *mri_aseg, VOXEL_LIST *vl, int min_points, double allowable_dist)
{
  int          ncontrol, max_points ;
  double       dist ;
  MRI          *mri_dist = NULL ;
  VOXEL_LIST  *vl_spline = NULL, *vl_spline_list ;

  vl->mri = mri_aseg ;
  max_points = MIN(10*min_points, vl->nvox) ;
  for (ncontrol =  min_points ; ncontrol <= max_points ; ncontrol++)
  {
    if (vl_spline)
      VLSTfree(&vl_spline) ;
    vl_spline = VLSTsplineFit(vl, ncontrol) ;
    vl_spline_list = VLSTinterpolate(vl_spline, 1) ;
    dist = VLSThausdorffDistance(vl, vl_spline_list, 10*allowable_dist, &mri_dist) ;
    VLSTfree(&vl_spline_list) ;
    if (dist < allowable_dist)
      break ;
  }
  if (mri_dist)
    MRIfree(&mri_dist) ;
  return(vl_spline) ;
}

#if 0
// used for computing posterior
static int burnin = 1000 ;
static int jump = 5 ;
#endif

static VOXEL_LIST *
find_optimal_spline(VOXEL_LIST *vl, MRI *mri_aseg, MRI *mri_wm_dist, int mcmc_samples, double spline_length_penalty, double spline_interior_penalty)
{
  VOXEL_LIST    *vl_spline, *vl_spline_optimal, *vl_spline_current ;
  double        energy, xn, yn, zn, best_energy, acceptance_val, current_energy ;
  int           n, cnum, label, label1, label2 ;

  label1 = MRIgetVoxVal(mri_aseg, vl->xi[0], vl->yi[0], vl->zi[0],0) ;
  label2 = MRIgetVoxVal(mri_aseg, vl->xi[vl->nvox-1], vl->yi[vl->nvox-1], vl->zi[vl->nvox-1],0) ;
  if (vl->nvox <= 2)
    ErrorReturn(NULL, (ERROR_NOFILE, "find_optimal_spline: input voxlist is too short (%d)", vl->nvox)) ;

  vl_spline_optimal = VLSTcopy(vl, NULL, 0, vl->nvox) ;
  best_energy = compute_spline_energy(vl_spline_optimal, mri_wm_dist,energy_flags, spline_length_penalty, spline_interior_penalty) ;
  if (DIAG_VERBOSE_ON)
  {
    char fname[STRLEN] ;
    VOXEL_LIST *vl_interp = VLSTinterpolate(vl_spline_optimal, .1) ;

    sprintf(fname, "spline.000.label") ;
    VLSTwriteLabel(vl_interp, fname, NULL, mri_aseg) ;
    sprintf(fname, "spline.000.cpts.label") ;
    VLSTwriteLabel(vl_spline_optimal, fname, NULL, mri_aseg) ;
    VLSTfree(&vl_interp) ;
  }

  vl_spline_current = VLSTcopy(vl_spline_optimal, NULL, 0, vl_spline_optimal->nvox) ;
  current_energy = best_energy ;
  for (n = 0 ; n < mcmc_samples ; n++)
  {
    cnum = (int)randomNumber(1.0, vl_spline_current->nvox-.001) ;
    vl_spline = VLSTcopy(vl_spline_current, NULL, 0, vl_spline_current->nvox) ;
    xn = proposal_sigma * PDFgaussian() ;
    yn = proposal_sigma * PDFgaussian() ;
    zn = proposal_sigma * PDFgaussian() ;

#if 0
    // project out tangential component
    if (cnum != vl_spline_current->nvox-1)   // don't project out tangential component for end point
    {
      double tx, ty, tz, dot, norm ;
      int    km1, kp1 ;

      km1 = cnum-1 ; kp1 = cnum+1 ;
      if (cnum == 0)
	km1 = 0 ;
      else if (cnum == vl_spline_current->nvox-1)
	kp1 = vl_spline_current->nvox-1 ;
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
#endif
    vl_spline->xd[cnum] += xn ; vl_spline->yd[cnum] += yn ; vl_spline->zd[cnum] += zn ;
    vl_spline->xi[cnum] = nint(vl_spline->xd[cnum]) ;
    vl_spline->yi[cnum] = nint(vl_spline->yd[cnum]) ;
    vl_spline->zi[cnum] = nint(vl_spline->zd[cnum]) ;
    if (cnum == vl_spline_current->nvox-1)  // last one must be in label2
    {
      label = MRIgetVoxVal(mri_aseg, vl_spline->xi[cnum], vl_spline->yi[cnum], vl_spline->zi[cnum],0) ;
      if (label != label2)  // last control point must be in label2
      {
	VLSTfree(&vl_spline) ;
	continue ;
      }
    }
    if (cnum == 0)  // last one must be in label1
    {
      label = MRIgetVoxVal(mri_aseg, vl_spline->xi[cnum], vl_spline->yi[cnum], vl_spline->zi[cnum],0) ;
      if (label != label1)  // first control point must be in label1
      {
	VLSTfree(&vl_spline) ;
	continue ;
      }
    }
    energy = compute_spline_energy(vl_spline, mri_wm_dist,energy_flags, spline_length_penalty, spline_interior_penalty) ;
    acceptance_val = exp((current_energy-energy)/acceptance_sigma) ;
    if (randomNumber(0.0, 1.0) < acceptance_val)
    {
      VLSTfree(&vl_spline_current) ;
      vl_spline_current = vl_spline ;
      current_energy = energy ;
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
	  VLSTwriteLabel(vl_interp, fname, NULL, mri_aseg) ;
	  sprintf(fname, "spline.%3.3d.cpts.label", n+1) ;
	  VLSTwriteLabel(vl_spline_optimal, fname, NULL, mri_aseg) ;
	  VLSTfree(&vl_interp) ;
	}
      }
    }
    else
      VLSTfree(&vl_spline) ;
  }

  if (DIAG_VERBOSE_ON)
    best_energy = compute_spline_energy(vl_spline_optimal, mri_wm_dist,energy_flags, spline_length_penalty, spline_interior_penalty) ;
  return(vl_spline_optimal) ;
}

static double
compute_spline_energy(VOXEL_LIST *vl, MRI *mri_wm_dist, int which, 
		      double spline_length_penalty, double spline_interior_penalty) 
{
  int         n, num ;
  VOXEL_LIST *vl_interp ;
  double     wm_dist ;
  double     interior_energy, length_energy, energy ;

  interior_energy = 0 ;
  vl_interp = VLSTinterpolate(vl, 1) ;

  for (n = num = 0 ; n < vl_interp->nvox ; n++)
  {
    if (vl_interp->xi[n] == Gx &&  vl_interp->yi[n] == Gy &&  vl_interp->zi[n] == Gz)
      DiagBreak() ;
    MRIsampleVolume(mri_wm_dist, vl_interp->xd[n], vl_interp->yd[n], vl_interp->zd[n], &wm_dist) ;
    if (wm_dist > 0)
      interior_energy = MAX(interior_energy, wm_dist);
    num++ ;
  }

  // add in specified penalty
  length_energy = vl_interp->nvox ;

  energy = 0 ;
  if (which & SPLINE_LENGTH)
    energy += spline_length_penalty*length_energy ;   // used to be sqrt(length_energy)
  if (which & SPLINE_WM_DIST)
    energy += spline_interior_penalty*interior_energy ;  // used to be sqrt(interior_energy) ;
  
  VLSTfree(&vl_interp) ;
  return(energy) ;
}
