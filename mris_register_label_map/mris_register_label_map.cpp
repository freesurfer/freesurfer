/**
 * @brief register the resting state patterns from one label map to another based on priors
 *
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



/*!
\file mris_register_label_map.c
\brief program to register an individual subject's resting state map to a group average
\author Bruce Fischl

*/



/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"

#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "mrisurf_project.h"

#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "cma.h"
#include "romp_support.h"
#include "voxlist.h"
#include "matrix.h"

#define MAX_SUBJECTS 1000
#define MAX_RUNS      100


static int powell_minimize(VECTOR **v_D, MATRIX **m_I, int nsubjects, MRI *mri_mask, VECTOR *v_weights) ;
static MRI *compute_voxlist_surface_correlations(VOXEL_LIST *vl, int num_maps, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat) ; 
static MRI *compute_voxlist_surface_correlations_across_runs(VOXEL_LIST *vl, int num_maps, MRI **mri_fvol, MRI **mri_fsurf, int runs, MRI *mri_dst) ;
VOXEL_LIST *rank_voxels(MRI *mri_stats, int num_maps) ;
//static MRI *compute_volume_correlations_in_surface_label(LABEL *area, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat, MRI *mri_mask)  ;
//static int compute_volume_correlation_at_vertex(MRI *mri_fvol, MRI *mri_fsurf, int vno, int frame, MRI *mri_cmat, MRI *mri_mask)  ;
//static MRI *accumulate_subcortical_map(MRI **mri_fvol, MRI **mri_fsurf, LABEL *area, int runs, MRI *mri_stats, MRI *mri_mask) ;
static MATRIX *compute_subcortical_map_weights(MRI_SURFACE *mris, MRI *mri_fvol[MAX_SUBJECTS][MAX_RUNS], MRI *mri_fsurf[MAX_SUBJECTS][MAX_RUNS], LABEL **labels, int runs, int nsubjects, MRI *mri_mask, int ndilate, char *out_prefix) ;
static int compute_surface_correlation_map_at_voxel(MRI *mri_fvol, MRI *mri_fsurf, int x, int y, int z, int frame, MRI *mri_cmat)  ;
#if 0
static MRI *compute_surface_correlations_in_aseg_label(MRI *mri_aseg, int aseg_label, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat)  ;
static MRI *average_aseg_label(MRI *mri_aseg, int aseg_label, MRI **mri_fvol, MRI **mri_fsurf, int runs, MRI *mri_dst) ;
#endif

static MRI *sample_fixed(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats) ;
static double compute_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_mov_avg, MRI *mri_stats) ;
static double sample_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, double *pvar) ;
static double sample_hemi_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, int frame, double *pvar) ;
static double sample_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MHT *mht, MRI *mri_stats, MRI *mri_corrmat_mov, MRI *mri_label_avg_mov, int vno, double x, double y, double z);
static double sample_hemi_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, MRI *mri_target_label, int vno, double x, double y, double z);
static int compute_warp_gradient(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed,  MRI_SURFACE *mris_rh_fixed,  MHT *mht_lh, MHT *mht_rh, MRI *mri_cmat_permuted,  MRI *mri_stats, MRI *mri_label_avg, MRI *mri_corrmat_mov, double dt) ;
static int compute_hemi_warp_gradient(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed,  MHT *mht,  MRI *mri_stats, MRI *mri_label_avg, double dt) ;
static int warp_hemi(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MRI *mri_target_label, MRI *mri_stats, int nsubjects, double tol, double dt)  ;
static int  warp_surface(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MRI *mri_cmat_mov, MRI *mri_stats, double tol, LABEL *area, int offset, double dt) ;
static int compute_vertex_permutation(MRI_SURFACE *mris_mov_mov, MRI_SURFACE *mris_fixed, MHT *mht,  int *vertices) ;
static MRI *average_within_label(MRI *mri_cmat, LABEL *area, int offset, MRI *mri_avg) ;
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static int write_snapshot(VECTOR *v_weights, MATRIX **m_I, MRI *mri_mask, const char *prefix, int iter, int nsubjects) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
static MRI *compute_mean_and_variance_across_frames(MRI **mri_label_avg, int nsubjects)  ;
static MRI *compute_mean_and_variance(MRI **mri_label_avg, int nsubjects)  ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

#define MAX_ITERS 500

static int create_only = 1 ;
static int ndilate = 0 ;
static int use_powell = 0 ;    
static int max_iters = MAX_ITERS ;
static int max_grad_averages = 16 ;
static int min_grad_averages = 16 ;
static MRI *mri_aseg = NULL ;
static int aseg_label = -1 ;
static int runs = 1 ;
static int downsample = 0 ;
static char *fmri_vol = NULL ;
static char *fmri_surf = NULL ;
static int num_maps = 100 ;

static char *cmat_name = NULL ;
static char *trgsubject = NULL ;
static char  *label_name = NULL ;
static char *prior_name = NULL ;
static char *TempVolFile=NULL;
static char *SUBJECTS_DIR = NULL;
static char **subjects = NULL ;
static int nsubjects = 0 ;
static char *hemi = NULL ;
static const char *ohemi = NULL ;
static int offset = 0 ;
static char *output_name = NULL ;
static int write_diags = 1 ;
static double tol = 0.01 ;   // terminate when error func % change is smaller than this
static double dt = .1 ;

static int subcortical_labels[] = 
{
  Left_Thalamus_Proper,
  Right_Thalamus_Proper,
  Left_Putamen,
  Right_Putamen,
  Left_Pallidum,
  Right_Pallidum,
  Left_Caudate,
  Right_Caudate,
  Brain_Stem,
  Left_VentralDC,
  Right_VentralDC,
  Left_Accumbens_area,
  Right_Accumbens_area,
  Left_Cerebellum_Cortex,
  Right_Cerebellum_Cortex
} ;
#define SUBCORTICAL_LABEL_COUNT (sizeof(subcortical_labels) / sizeof(subcortical_labels[0]))



/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int         nargs, n, ic, nvert;
  char        *subject, fname[STRLEN] ;
  MRI         *mri_cmat, *mri_prior, *mri_stats = NULL, *mri_label_avg[MAX_SUBJECTS], *mri_cmat_mov,
    *mri_fvol[MAX_SUBJECTS][MAX_RUNS], *mri_fsurf[MAX_SUBJECTS][MAX_RUNS], *mri_target_label = NULL ;
  LABEL       *labels[MAX_SUBJECTS], *target_area ;
  MRI_SURFACE *mris_lh_mov, *mris_rh_mov, *mris_lh_fixed, *mris_rh_fixed ;
  MATRIX      *m_map_weights ;

  nargs = handleVersionOption(argc, argv, "mris_register_label_map");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  if (SUBJECTS_DIR == NULL)
  {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR == NULL) {
      printf("ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  if (downsample > 0 && mri_aseg)
  {
    MRI *mri_tmp2 = MRIdownsampleN(mri_aseg, NULL, downsample, downsample, downsample, 0), *mri_tmp ;
    mri_tmp = MRIresample(mri_aseg, mri_tmp2, SAMPLE_NEAREST) ;
    MRIfree(&mri_aseg) ; mri_aseg = mri_tmp ; MRIfree(&mri_tmp2) ;
  }
  for (n = 0 ; n < nsubjects ; n++)
  {
    subject = subjects[n] ;
    printf("processing subject %s\n", subject) ;

    if (fmri_vol)   // using correlation with subcortical  voxels to drive warp
    {
      char name[STRLEN] ;
      int  run ;

      for (run = 0 ; run < runs ; run++)
      {
	MRI *mri_tmp ;
	sprintf(name, fmri_surf, hemi, run+1) ;
	sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, subject, name) ;
	mri_fsurf[n][run] = MRIread(fname) ;
	if (mri_fsurf[n][run] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not load surface time series from %s", Progname, fname) ;
	if (mri_fsurf[n][run]->height != 1 || mri_fsurf[n][run]->depth != 1)
	{
	  mri_tmp = mri_reshape(mri_fsurf[n][run], 
				mri_fsurf[n][run]->width*mri_fsurf[n][run]->height*mri_fsurf[n][run]->depth,
				1, 1, mri_fsurf[n][run]->nframes) ;
	  MRIfree(&mri_fsurf[n][run]) ;
	  mri_fsurf[n][run] = mri_tmp ;
	}

	sprintf(name, fmri_vol, run+1) ;
	sprintf(fname, "%s/%s/mri/%s", SUBJECTS_DIR, subject, name) ;
	mri_fvol[n][run] = MRIread(fname) ;
	if (mri_fvol[n][run] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not load volume time series from %s", Progname, fname) ;
      }
      if (label_name)
      {
        sprintf(fname, "%s/%s/label/%s.%s", SUBJECTS_DIR, subject, hemi, label_name) ;
        labels[n] = LabelRead("", fname) ;
        if (!labels[n])
	  ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", Progname, fname) ;
      }
      else
      {
	for (run = 0 ; run < runs ; run++)
        {
	  MRIfree(&mri_fvol[n][run]) ;
	  MRIfree(&mri_fsurf[n][run]) ;
        }
      }
    }
    else   // using surface-based correlation to drive warp
    {
      sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, subject, cmat_name) ;
      mri_cmat = MRIread(fname) ;
      if (!mri_cmat)
	ErrorExit(ERROR_NOFILE, "%s: could not read cmat from %s", Progname, fname) ;
      
      sprintf(fname, "%s/%s/label/%s.%s", SUBJECTS_DIR, subject, hemi, label_name) ;
      labels[n] = LabelRead("", fname) ;
      if (!labels[n])
	ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", Progname, fname) ;
      if (stricmp(hemi, "rh") == 0)
	offset = mri_cmat->width/2 ;  // lower right block is for rh
      mri_label_avg[n] = average_within_label(mri_cmat, labels[n], offset, NULL) ;
      MRIfree(&mri_cmat) ;

      if (write_diags)
      {
	sprintf(fname, "%s.%s.%s.label_avg.mgz", hemi, subject, output_name) ;
	mri_label_avg[n]->width /= 2 ;
	MRIwrite(mri_label_avg[n], fname) ;
	mri_label_avg[n]->width *= 2 ;
      }
    }
  }
  if (mri_aseg == NULL)  // using surface-based label to drive warp
  {
    mri_stats = compute_mean_and_variance_across_frames(mri_label_avg, nsubjects) ;
    if (write_diags)
    {
      sprintf(fname, "%s.%s.label_avg.mgz", hemi, output_name) ;
      mri_stats->width /= 2 ;
      MRIwriteFrame(mri_stats, fname, 0) ;
      mri_stats->width *= 2 ;
    }
  }
  else   // using subcortical maps to drive warp
  {
    char       name[STRLEN] ;
    int        run ;
//    VOXEL_LIST *vl ;

    if (label_name)  // compute subcortical map that has highest correlation with surface labels
    {
      MRI        *mri_mask ;
      int        n ;
      MRI_SURFACE *mris ;
      char        fname[STRLEN] ;

      mri_mask = MRIclone(mri_aseg, NULL) ;
      for (n = 0 ; n < SUBCORTICAL_LABEL_COUNT ; n++)
	MRIcopyLabel(mri_aseg, mri_mask, subcortical_labels[n]) ;

#if 0
      mri_stats = NULL ;
      for (n = 0 ; n < nsubjects ; n++)
	mri_stats = accumulate_subcortical_map(mri_fvol[n], mri_fsurf[n], labels[n], runs, mri_stats, mri_mask) ;
      MRIcomputeMeanAndStandardDeviation(mri_stats, mri_stats, nsubjects) ;
      vl = rank_voxels(mri_stats, num_maps) ;
      for (n = 0 ; n < nsubjects ; n++)
	mri_label_avg[n] = compute_voxlist_surface_correlations_across_runs(vl, num_maps, mri_fvol[n], mri_fsurf[n], runs, NULL) ;
      MRIfree(&mri_stats) ;
      mri_stats = compute_mean_and_variance(mri_label_avg, nsubjects) ;
#endif
      sprintf(fname, "%s/fsaverage%d/surf/%s.inflated", SUBJECTS_DIR, 5, hemi) ;
      mris = MRISread(fname) ;
      if (!mris)
	ErrorExit(ERROR_NOFILE, "%s: could not read *surface from %s", Progname, fname) ;
      //  load the data for the target subject into the array spot [nsubjects]
      subjects[nsubjects] = trgsubject ;
      for (run = 0 ; run < runs ; run++)
      {
	MRI *mri_tmp ;
	sprintf(name, fmri_surf, hemi, run+1) ;
	sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, trgsubject, name) ;
	mri_fsurf[nsubjects][run] = MRIread(fname) ;
	if (mri_fsurf[nsubjects][run] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not load surface time series from %s", Progname, fname) ;
	if (mri_fsurf[nsubjects][run]->height != 1 || mri_fsurf[nsubjects][run]->depth != 1)
	{
	  mri_tmp = mri_reshape(mri_fsurf[nsubjects][run], 
				mri_fsurf[nsubjects][run]->width*mri_fsurf[nsubjects][run]->height*mri_fsurf[nsubjects][run]->depth,
				1, 1, mri_fsurf[nsubjects][run]->nframes) ;
	  MRIfree(&mri_fsurf[nsubjects][run]) ;
	  mri_fsurf[nsubjects][run] = mri_tmp ;
	}
      
	sprintf(name, fmri_vol, run+1) ;
	sprintf(fname, "%s/%s/mri/%s", SUBJECTS_DIR, trgsubject, name) ;
	mri_fvol[nsubjects][run] = MRIread(fname) ;
	if (mri_fvol[nsubjects][run] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not load volume time series from %s", Progname, fname) ;
      }
      if (label_name)
      {
        sprintf(fname, "%s/%s/label/%s.%s", SUBJECTS_DIR, trgsubject, hemi, label_name) ;
        labels[nsubjects] = LabelRead("", fname) ;
        if (!labels[nsubjects])
	  ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", Progname, fname) ;
      }
      m_map_weights = compute_subcortical_map_weights(mris, mri_fvol, mri_fsurf, labels, runs, nsubjects, mri_mask, ndilate, output_name) ;
      mri_target_label = NULL ; // FIX THIS!
      MRIfree(&mri_mask) ;
    }    
    else
      mri_stats = compute_mean_and_variance(mri_label_avg, nsubjects) ;
    if (write_diags)
    {
      sprintf(fname, "%s.%s.label_avg.mgz", hemi, output_name) ;
      printf("writing mean and variance volume to %s\n", fname) ;
      MRIwrite(mri_stats, fname) ;
    }

    //  load the data for the target subject
    for (run = 0 ; run < runs ; run++)
    {
      MRI *mri_tmp ;
      sprintf(name, fmri_surf, hemi, run+1) ;
      sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, trgsubject, name) ;
      mri_fsurf[n][run] = MRIread(fname) ;
      if (mri_fsurf[n][run] == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load surface time series from %s", Progname, fname) ;
      if (mri_fsurf[n][run]->height != 1 || mri_fsurf[n][run]->depth != 1)
      {
	mri_tmp = mri_reshape(mri_fsurf[n][run], 
			      mri_fsurf[n][run]->width*mri_fsurf[n][run]->height*mri_fsurf[n][run]->depth,
			      1, 1, mri_fsurf[n][run]->nframes) ;
	MRIfree(&mri_fsurf[n][run]) ;
	mri_fsurf[n][run] = mri_tmp ;
      }
      
      sprintf(name, fmri_vol, run+1) ;
      sprintf(fname, "%s/%s/mri/%s", SUBJECTS_DIR, trgsubject, name) ;
      mri_fvol[n][run] = MRIread(fname) ;
      if (mri_fvol[n][run] == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load volume time series from %s", Progname, fname) ;
    }
    
//   mri_target_label = compute_voxlist_surface_correlations_across_runs(vl, num_maps, mri_fvol[n], mri_fsurf[n], runs, NULL) ;
    
    for (run = 0 ; run < runs ; run++)
    {
      MRIfree(&mri_fvol[n][run]) ;
      MRIfree(&mri_fsurf[n][run]) ;
    }
  }

  ic  = 5 ;
  nvert = mri_aseg ? MRInvox(mri_label_avg[0])/mri_label_avg[0]->nframes : MRInvox(mri_label_avg[0])/2 ;
  switch (nvert)
  {
  case 40962: ic = 6 ; break ;
  case 10242: ic = 5 ; break ;
  default: ErrorExit(ERROR_UNSUPPORTED, "only IC 5 or 6 supported (nvertices = %d detected)\n", MRInvox(mri_label_avg[0])/2);
  }

  sprintf(fname, "%s/fsaverage%d/surf/lh.sphere", SUBJECTS_DIR, ic) ;
  mris_lh_mov = MRISread(fname) ;
  if (!mris_lh_mov)
    ErrorExit(ERROR_NOFILE, "%s: could not read *surface from %s", Progname, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere", SUBJECTS_DIR, ic) ;
  mris_rh_mov = MRISread(fname) ;
  if (!mris_rh_mov)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, fname) ;

  mris_lh_fixed = MRISread(fname) ;
  if (!mris_lh_fixed)
    ErrorExit(ERROR_NOFILE, "%s: could not read *surface from %s", Progname, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere", SUBJECTS_DIR, ic) ;
  mris_rh_fixed = MRISread(fname) ;
  if (!mris_rh_fixed)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, fname) ;

  MRISsetNeighborhoodSizeAndDist(mris_lh_mov, 2) ;
  MRISsetNeighborhoodSizeAndDist(mris_rh_mov, 2) ;
  MRISsetNeighborhoodSizeAndDist(mris_lh_fixed, 2) ;
  MRISsetNeighborhoodSizeAndDist(mris_rh_fixed, 2) ;
  MRIScomputeSecondFundamentalForm(mris_lh_mov) ;
  MRIScomputeSecondFundamentalForm(mris_rh_mov) ;
  MRIScomputeSecondFundamentalForm(mris_lh_fixed) ;
  MRIScomputeSecondFundamentalForm(mris_rh_fixed) ;

  if (mri_aseg == NULL)   // using surface label and priors to drive warp
  {
    sprintf(fname, "%s/fsaverage%d/label/%s.%s", SUBJECTS_DIR, ic, hemi, prior_name) ;
    mri_prior = MRIread(fname) ;
    if (!mri_prior)
      ErrorExit(ERROR_NOFILE, "%s: could not prior overlay from %s", Progname, fname) ;
    
    sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, trgsubject, cmat_name) ;
    mri_cmat_mov = MRIread(fname) ;
    if (!mri_cmat_mov)
      ErrorExit(ERROR_NOFILE, "%s: could not read movable cmat from %s", Progname, fname) ;
    
    sprintf(fname, "%s/fsaverage%d/label/%s.MT.thresh.label", SUBJECTS_DIR, ic, hemi) ;
    target_area = LabelRead("", fname) ;
    if (!target_area)
      ErrorExit(ERROR_NOFILE, "%s: could not read movable cmat from %s", Progname, fname) ;

    
    warp_surface(mris_lh_mov, mris_rh_mov, mris_lh_fixed, mris_rh_fixed, mri_cmat_mov, mri_stats, tol, target_area, offset, dt) ;
  }
  else  // warping based on correlations with subcortical label(s)
  {
    warp_hemi(mris_lh_mov, mris_lh_fixed, mri_target_label, mri_stats, nsubjects, tol, dt) ;
  }
  sprintf(fname, "%s/fsaverage%d/surf/lh.sphere.%s", SUBJECTS_DIR, ic, output_name) ;
  printf("writing output surface to %s\n", fname) ;
  MRISwrite(mris_lh_mov, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere.%s", SUBJECTS_DIR, ic, output_name) ;
  printf("writing output surface to %s\n", fname) ;
  MRISwrite(mris_rh_mov, fname) ;
  return 0;
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--sdir"))
    {
      SUBJECTS_DIR = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--prior"))
    {
      prior_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--create_only"))
    {
      create_only = 1 ;
      nargsused = 0 ;
    }
    else if (!strcasecmp(option, "--dilate"))
    {
      ndilate = atoi(pargv[0]) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--debug_voxel"))
    {
      Gx = atoi(pargv[0]) ;
      Gy = atoi(pargv[1]) ;
      Gz = atoi(pargv[2]) ;
      nargsused = 3 ;
      printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    }
    else if (!strcasecmp(option, "--ds"))
    {
      downsample = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("downsampling aseg volume %d times (time series will be downsampled to match resolution)\n", 
	     downsample) ;
    }
    else if (!strcasecmp(option, "--averages"))
    {
      max_grad_averages = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("averaging gradient %d times before applying\n", max_grad_averages) ;
      if (max_grad_averages < min_grad_averages)
	min_grad_averages = max_grad_averages ;
    }
    else if (!strcasecmp(option, "--min_averages"))
    {
      min_grad_averages = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("minimum averaging of gradients set to %d times\n", min_grad_averages) ;
    }
    else if (!strcasecmp(option, "--max_iters"))
    {
      max_iters = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("limiting minimization to %d iterations\n", max_iters) ;
    }
    else if (!strcasecmp(option, "--maps"))
    {
      num_maps = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("using top %d voxel maps for registration\n", num_maps) ;
    }
    else if (!strcasecmp(option, "--dt"))
    {
      dt = atof(pargv[0]) ;
      nargsused = 1 ;
      printf("using timestep dt = %2.3f\n", dt) ;
    }
    else if (!strcasecmp(option, "--runs"))
    {
      runs = atoi(pargv[0]) ;
      nargsused = 1 ;
      printf("loading %d runs\n", runs) ;
    }
    else if (!strcasecmp(option, "--fmri"))
    {
      fmri_vol = pargv[0] ;
      fmri_surf = pargv[1] ;
      nargsused = 2 ;
      printf("loading volume (%s) and surface (%s) rfSMRI time courses for each subject\n",
	     fmri_vol, fmri_surf) ;
    }
    else if (!strcasecmp(option, "--aseg"))
    {
      mri_aseg = MRIread(pargv[0]) ;
      if (mri_aseg == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load aseg from %s", Progname, pargv[0]) ;
      aseg_label = atoi(pargv[1]) ;
      printf("using maps from label %s (%d) to align cortices\n", cma_label_to_name(aseg_label), aseg_label) ;
      nargsused = 2 ;
    }
    else if (!strcasecmp(option, "--v"))
    {
      Gdiag_no = atoi(pargv[0]) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--tol"))
    {
      tol = atof(pargv[0]) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--trgsubject"))
    {
      trgsubject = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--output"))
    {
      output_name = pargv[0] ;
      printf("setting output name = %s\n", output_name) ;
      fflush(stdout) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--hemi"))
    {
      hemi = pargv[0] ;
      nargsused = 1 ;
      if (strcmp(hemi, "lh") == 0)
	ohemi = "rh" ;
      else
	ohemi = "lh" ;
    }
    else if (!strcasecmp(option, "--subjects"))
    {
      int nextarg, i ;

      for (nextarg = 0 ; nextarg < nargc ; nextarg++)
	if (pargv[nextarg][0] == '-' && pargv[nextarg][1] == '-')
	  break ;
      printf("parsing --subjects with argc = %d\n", nextarg) ;
      fflush(stdout) ;
      nsubjects = nargsused = nextarg ;
//      subjects = pargv ;
      subjects = (char **)calloc(nsubjects+1, sizeof(char *)) ;
      for (i = 0 ; i < nsubjects ; i++)
	subjects[i] = pargv[i] ;
    }
    else if (!strcasecmp(option, "--label"))
    {
      label_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--cmat"))
    {
      cmat_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--temp-vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      TempVolFile = pargv[0];
      nargsused = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void usage_exit(void)
\brief Prints usage and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_usage(void)
\brief Prints usage and returns (does not exit)
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --subjects  SUBJECT_LIST: list of training subjects \n");
  printf("   --trgsubject SUBJECT: name of target subject \n");
  printf("   --prior      PRIOR: name of prior surface overlay \n");
  printf("   --label      LABEL: name of label for each subject \n");
  printf("   --temp-vol volfile : template volume \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --sdir      SUBJECTS_DIR\n");
  printf("   --version   print out version and exit\n");
  printf("   --v         VNO  debug this vertex\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  return;
}

/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void dump_options(FILE *fp)
\brief Prints command-line options to the given file pointer
\param FILE *fp - file pointer
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}

static MRI *
compute_mean_and_variance(MRI **mri_label_avg, int nsubjects) 
{
  MRI    *mri_stats, *mri ;
  int    n, x, y, z, f, nframes ;
  float  val, mean, var ;

  mri = mri_label_avg[0] ;
  // the first N frames are means, the 2nd N are variances
  nframes = mri_label_avg[0]->nframes ;
  mri_stats = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 2*nframes) ;
  MRIcopyHeader(mri, mri_stats) ;
  for (f = 0 ; f < mri->nframes ; f++)
  {
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
	for (z = 0 ; z < mri->depth ; z++)
	  for (n = 0 ; n < nsubjects ; n++)
	  {
	    mri = mri_label_avg[n] ;
	    val = MRIgetVoxVal(mri, x, y, z, f) ;
	    mean = MRIgetVoxVal(mri_stats, x, y, z, f) ;
	    var = MRIgetVoxVal(mri_stats, x, y, z, f+nframes) ;
	    mean += val ; var += (val*val) ;
	    MRIsetVoxVal(mri_stats, x, y, z, f, mean) ;
	    MRIsetVoxVal(mri_stats, x, y, z, f+nframes, var) ;
	  }
  }

  for (f = 0 ; f < mri->nframes ; f++)
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
	for (z = 0 ; z < mri->depth ; z++)
	{
	  mean = MRIgetVoxVal(mri_stats, x, y, z, f) ;
	  var = MRIgetVoxVal(mri_stats, x, y, z, f+nframes) ;
	  mean /= nsubjects ;
	  var = var / nsubjects - mean*mean ;
	  if (var < .0001)
	    var = .01 ; 
	  MRIsetVoxVal(mri_stats, x, y, z, f, mean) ;
	  MRIsetVoxVal(mri_stats, x, y, z, f+nframes, var) ;
	}
  return(mri_stats) ;
}

static MRI *
compute_mean_and_variance_across_frames(MRI **mri_label_avg, int nsubjects) 
{
  MRI    *mri_stats, *mri ;
  int    n, x, y, z, f ;
  float  val, mean, var ;

  mri = mri_label_avg[0] ;
  mri_stats = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 2) ;
  MRIcopyHeader(mri, mri_stats) ;
  for (n = 0 ; n < nsubjects ; n++)
  {
    mri = mri_label_avg[n] ;
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
	for (z = 0 ; z < mri->depth ; z++)
	  for (f = 0 ; f < mri->nframes ; f++)
	  {
	    val = MRIgetVoxVal(mri, x, y, z, f) ;
	    mean = MRIgetVoxVal(mri_stats, x, y, z, 0) ;
	    var = MRIgetVoxVal(mri_stats, x, y, z, 1) ;
	    mean += val ; var += (val*val) ;
	    MRIsetVoxVal(mri_stats, x, y, z, 0, mean) ;
	    MRIsetVoxVal(mri_stats, x, y, z, 1, var) ;
	  }
  }

  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
	mean = MRIgetVoxVal(mri_stats, x, y, z, 0) ;
	var = MRIgetVoxVal(mri_stats, x, y, z, 1) ;
	mean /= nsubjects ;
	var = var / nsubjects - mean*mean ;
	if (var < .0001)
	  var = .01 ; 
	MRIsetVoxVal(mri_stats, x, y, z, 0, mean) ;
	MRIsetVoxVal(mri_stats, x, y, z, 1, var) ;
      }
  return(mri_stats) ;
}

static MRI *
average_within_label(MRI *mri_cmat, LABEL *area, int offet, MRI *mri_avg) 
{
  int  n, vno, x, y, z ;
  float avg, val ;

  if (mri_avg == NULL)
  {
    mri_avg = MRIallocSequence(mri_cmat->width, mri_cmat->height, mri_cmat->depth, MRI_FLOAT, 1) ;
    MRIcopyHeader(mri_cmat, mri_avg) ;
  }

  for (x = 0 ; x < mri_cmat->width ; x++)
    for (y = 0 ; y < mri_cmat->height ; y++)
      for (z = 0 ; z < mri_cmat->depth ; z++)
      {
	for (avg = 0.0, n = 0 ; n < area->n_points ; n++)
	{
	  vno = area->lv[n].vno + offset ;
	  if (vno == x)
	    continue ;  // don't include diagonal in estimate
	  
	  val = MRIgetVoxVal(mri_cmat, x, y, z, vno) ;
	  if (x == Gdiag_no)
	  {
	    printf("vno %d, map %d, adding %2.2f\n", x, vno, val) ;
	    DiagBreak() ;
	  }
	  avg += val ;
	}
	avg /= (area->n_points-1) ;  // didn't  include diagonal
	if (fabs(avg) > 1000)
	  DiagBreak() ;
	MRIsetVoxVal(mri_avg, x, y, z, 0, avg) ;
      }

  return(mri_avg) ;
}


static int
compute_vertex_permutation(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht,  int *vertices)
{
  int    vno ;
  VERTEX *v, *vfixed ;

  for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    float min_dist;
    vfixed = MHTfindClosestVertex2(mht, mris_fixed, mris_mov, v, &min_dist) ;
    vertices[vno] = vfixed - mris_fixed->vertices ;
    
  }    
  
  return(NO_ERROR) ;
}

static MRI *
permute_corrmat(MRI *mri_corrmat_src, int *lh_vertices, int *rh_vertices, MRI *mri_corrmat_dst)
{
  int vno_src, vno_dst, r, offset ;
  float val ;

  if (mri_corrmat_dst == NULL)
    mri_corrmat_dst = MRIclone(mri_corrmat_src, NULL) ;

  offset = mri_corrmat_src->width/2 ;
  for (vno_src = 0 ; vno_src < mri_corrmat_src->width ; vno_src++)
  {
    if (vno_src < offset)   // lh
    {
      vno_dst = lh_vertices[vno_src] ;
      if (vno_src != vno_dst)
	DiagBreak() ;
    }
    else
    {
      vno_dst = rh_vertices[vno_src-offset] ;
      if (vno_src-offset != vno_dst)
	DiagBreak() ;
    }

    // swap rows
    for (r = 0 ; r <  mri_corrmat_src->width ; r++)
    {
      val = MRIgetVoxVal(mri_corrmat_src, vno_src, 0,0, r) ;
      MRIsetVoxVal(mri_corrmat_dst, vno_dst,0, 0, r, val) ;
    }
    // swap cols
    for (r = 0 ; r <  mri_corrmat_src->width ; r++)
    {
      val = MRIgetVoxVal(mri_corrmat_src, r, 0,0, vno_src) ;
      MRIsetVoxVal(mri_corrmat_dst, r,0, 0, vno_dst, val) ;
    }
  }

  return(mri_corrmat_dst) ;
}


static double
compute_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_mov_avg, MRI *mri_stats)
{
  int     vno ;
  double  sse ;
  double  val, mean, var ;
  VERTEX  *v ;

  for (sse = 0.0, vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    mean = sample_stats(mris_mov, mris_fixed, mht, mri_stats, v->x, v->y, v->z, &var) ;
    val = MRIgetVoxVal(mri_mov_avg, vno, 0,0, 0) ;
    if (FZERO(var))
      var = 1.0 ;
    sse += SQR(mean-val)/var ;
  }
  return(sqrt(sse/mri_mov_avg->width)) ;
}

static double
sample_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MHT *mht, MRI *mri_stats, MRI *mri_corrmat_mov, MRI *mri_label_avg, int vno, double x, double y, double z)
{
  double mean, var, val, error ;

  mean = sample_stats(mris_mov, mris_lh_fixed, mht, mri_stats, x, y, z, &var) ;
  val = MRIgetVoxVal(mri_label_avg, vno, 0, 0, 0) ;
  error = SQR(val-mean)/var ;
  
  return(error) ;
}
static MRI *
sample_fixed(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats) 
{
  MRI    *mri ;
  double mean, var ;
  int    vno ;
  VERTEX *v ;

  mri = MRIalloc(mris_mov->nvertices, 1, 1, MRI_FLOAT) ;
  MRIcopyHeader(mri_stats, mri) ;
  for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    mean = sample_stats(mris_mov, mris_fixed, mht, mri_stats, v->x, v->y, v->z, &var) ;
    MRIsetVoxVal(mri, vno, 0, 0, 0, mean) ;
  }
  return(mri) ;
}
static double
sample_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, double *pvar)
{
  int   fno ;
  FACE  *face ;
  double fdist, val0, val1, val2, mean, var ;

  MHTfindClosestFaceGeneric(mht, mris_fixed, x, y, z, 8000, -1, 0, &face, &fno, &fdist);

  // interpolate mean
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,0) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,0) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,0) ;
  mean = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(mean)>100)
    DiagBreak() ;

  // interpolate variance
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,1) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,1) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,1) ;
  var = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(var)>100)
    DiagBreak() ;

  *pvar = var ;
  return(mean) ;
}


static int
compute_warp_gradient(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed,  MRI_SURFACE *mris_rh_fixed,  MHT *mht_lh, MHT *mht_rh, MRI *mri_cmat_permuted,  MRI *mri_stats, MRI *mri_label_avg, MRI *mri_corrmat_mov, double dt)
{
  int    vno ;
  VERTEX *v ;
  double  e1p, e1m, e2p, e2m, x, y, z, dx, dy, dz ;

  for (vno = 0 ; vno < mris_lh_mov->nvertices ; vno++)
  {
    v = &mris_lh_mov->vertices[vno] ;

    // sample vertex energy in plus and minus two principal directions
    x = v->x + v->e1x*dt ; y = v->y + v->e1y*dt ; z = v->z + v->e1z*dt ; 
    e1p  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x - v->e1x*dt ; y = v->y - v->e1y*dt ; z = v->z - v->e1z*dt ; 
    e1m  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x + v->e2x*dt ; y = v->y + v->e2y*dt ; z = v->z + v->e2z*dt ; 
    e2p  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x - v->e2x*dt ; y = v->y - v->e2y*dt ; z = v->z - v->e2z*dt ; 
    e2m  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    if (e1p < e1m && e1p < e2p && e1p < e2m)   // move in e1p dir
    {
      dx = v->e1x*dt ; dy = v->e1y*dt ; dz = v->e1z*dt ;
    }
    else if (e1m < e2p && e1m < e2m)
    {
      dx = -v->e1x*dt ; dy = -v->e1y*dt ; dz = -v->e1z*dt ;
    }
    else if (e2m < e2p)
    {
      dx = -v->e2x*dt ; dy = -v->e2y*dt ; dz = -v->e2z*dt ;
    }
    else
    {
      dx = v->e2x*dt ; dy = v->e2y*dt ; dz = v->e2z*dt ;
    }
    v->dx += dx ; v->dy += dy ; v->dz += dz ;
  }

  return(NO_ERROR) ;
}

static int
warp_surface(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MRI *mri_cmat_mov, MRI *mri_stats, double tol, LABEL *area, int offset, double dt) 
{
  int  *lh_vertices, *rh_vertices, iter = 0 ;
  MHT  *mht_lh, *mht_rh, *mht_lh_faces, *mht_rh_faces ;
  MRI  *mri_cmat_permuted = NULL, *mri_mov_avg ;
  double last_error, error, pct_change ;

  lh_vertices = (int *)calloc(mris_lh_fixed->nvertices, sizeof(int)) ;
  rh_vertices = (int *)calloc(mris_rh_fixed->nvertices, sizeof(int)) ;
  mht_lh = MHTcreateVertexTable_Resolution(mris_lh_fixed, CURRENT_VERTICES, ceil(mris_lh_fixed->avg_vertex_dist));
  mht_rh = MHTcreateVertexTable_Resolution(mris_rh_fixed, CURRENT_VERTICES, ceil(mris_rh_fixed->avg_vertex_dist));
  mht_lh_faces = MHTcreateFaceTable_Resolution(mris_lh_fixed, CURRENT_VERTICES, ceil(mris_lh_fixed->avg_vertex_dist)) ;
  mht_rh_faces = MHTcreateFaceTable_Resolution(mris_rh_fixed, CURRENT_VERTICES, ceil(mris_rh_fixed->avg_vertex_dist)) ;
  compute_vertex_permutation(mris_lh_mov, mris_lh_fixed, mht_lh, lh_vertices) ;
  compute_vertex_permutation(mris_rh_mov, mris_rh_fixed, mht_rh, rh_vertices) ;
  mri_cmat_permuted = permute_corrmat(mri_cmat_mov, lh_vertices, rh_vertices, mri_cmat_permuted) ;
  mri_mov_avg = average_within_label(mri_cmat_permuted, area, offset, NULL) ;
  last_error = compute_error(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_mov_avg, mri_stats) ;
  if (write_diags)
  {
    char fname[STRLEN] ;
    MRI  *mri_fixed ; 

    mri_fixed = sample_fixed(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_stats) ;
    sprintf(fname, "lh.label_avg.%3.3d.mgz",  iter) ;
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_fixed, fname) ;
    MRIfree(&mri_fixed) ;
  }

  printf("%3.3d: error = %2.3f\n", iter, last_error) ;
  do
  {
    compute_warp_gradient(mris_lh_mov, mris_rh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh_faces, mht_rh_faces, mri_cmat_permuted, mri_stats, mri_mov_avg, mri_cmat_permuted, dt) ;
    MRISapplyGradient(mris_lh_mov, 1) ; MRISprojectOntoSphere(mris_lh_mov, mris_lh_mov, mris_lh_mov->radius) ;
    MRISapplyGradient(mris_rh_mov, 1) ; MRISprojectOntoSphere(mris_rh_mov, mris_rh_mov, mris_rh_mov->radius) ;

    mri_cmat_permuted = permute_corrmat(mri_cmat_mov, lh_vertices, rh_vertices, mri_cmat_permuted) ;
    average_within_label(mri_cmat_permuted, area, offset, mri_mov_avg) ;
    if (write_diags)
    {
      char fname[STRLEN] ;
      MRI  *mri_fixed ; 

      sprintf(fname, "lh.label_avg.%3.3d.mgz",  iter+1) ;
      printf("writing snapshot to %s\n", fname) ;
      mri_fixed = sample_fixed(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_stats) ;
      MRIwrite(mri_fixed, fname) ;
      MRIfree(&mri_fixed) ;
    }
    compute_vertex_permutation(mris_lh_mov, mris_lh_fixed, mht_lh, lh_vertices) ;
    compute_vertex_permutation(mris_rh_mov, mris_rh_fixed, mht_rh, rh_vertices) ;
    error = compute_error(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_mov_avg, mri_stats) ;
    pct_change = 100*(last_error - error) / last_error ;
    printf("%3.3d: error = %2.3f (%2.3f%%)\n", iter+1, error, pct_change) ;
    last_error = error ;
  } while (pct_change > tol && iter++ < max_iters) ;
    

  free(lh_vertices) ; free(rh_vertices) ; MHTfree(&mht_lh) ; MHTfree(&mht_rh) ;
  return(NO_ERROR) ;
}
static int
compute_surface_correlation_map_at_voxel(MRI *mri_fvol, MRI *mri_fsurf, int x, int y, int z, int frame, MRI *mri_cmat) 
{
  int   vno, t ;
  float corr, sval, vval, stotal, stotalsq, vtotal, vtotalsq, smean, vmean, sstd, vstd, denom ;
  
  for ( vtotal = vtotalsq = 0.0, t = 0 ;  t < mri_fvol->nframes ; t++)
  {
    vval = MRIgetVoxVal(mri_fvol, x, y, z, t) ;
    vtotal += vval ; vtotalsq += vval*vval ;
  }
  vmean = vtotal / (float)mri_fvol->nframes ;
  vstd = sqrt(vtotalsq / (float)mri_fvol->nframes  - vmean*vmean) ;
  
  if (FZERO(vstd))   // volume time series is empty
  {
    for (vno = 0 ; vno < mri_cmat->width ; vno++)
      MRIsetVoxVal(mri_cmat, vno, 0, 0, frame, -10) ;
    return(NO_ERROR) ;
  }

  for (vno = 0 ; vno < mri_cmat->width ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (stotal = stotalsq = 0.0, t = 0 ;  t < mri_fsurf->nframes ; t++)
    {
      sval = MRIgetVoxVal(mri_fsurf, vno, 0, 0, t) ;
      stotal += sval ; stotalsq += sval*sval ;
    }
    smean = stotal / (float)mri_fsurf->nframes ;
    sstd = sqrt(stotalsq / (float)mri_fsurf->nframes  - smean*smean) ;
    if (FZERO(sstd))
    {
      DiagBreak() ;
      corr = -10;
    }
    else
    {
      denom = sstd*vstd ;
      for (corr = 0.0, t = 0 ;  t < mri_fsurf->nframes ; t++)
      {
	sval = MRIgetVoxVal(mri_fsurf, vno, 0, 0, t) ;
	vval = MRIgetVoxVal(mri_fvol, x, y, z, t) ;
	corr += (sval-smean)*(vval-vmean) ;
      }
      corr = corr / ((float)mri_fsurf->nframes * denom) ;
      if (fabs(corr) > 1)
	DiagBreak() ;
    }
    MRIsetVoxVal(mri_cmat, vno, 0, 0, frame, corr) ;
  }
  return(NO_ERROR) ;
}

#if 0

static MRI *
compute_surface_correlations_in_aseg_label(MRI *mri_aseg, int aseg_label, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat) 
{
  int  nvox, x, y, z, n ;

  nvox = MRIvoxelsInLabel(mri_aseg, aseg_label) ;
  mri_cmat = MRIallocSequence(mri_fsurf->width*mri_fsurf->height*mri_fsurf->depth,1,1,MRI_FLOAT,nvox) ;

  // compute correlation maps for every point in the target label
  for (n = x = 0 ; x < mri_aseg->width ; x++)
  {
    for (y = 0 ; y < mri_aseg->height ; y++)
    {
      for (z = 0 ; z < mri_aseg->depth  ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	if (MRIgetVoxVal(mri_aseg, x, y, z, 0) == aseg_label)
	{
	  compute_surface_correlation_map_at_voxel(mri_fvol, mri_fsurf, x, y, z, n, mri_cmat) ;
	  n++ ;
	}
      }
    }
  }

  return(mri_cmat) ;
}

static MRI *
average_aseg_label(MRI *mri_aseg, int aseg_label, MRI **mri_fvol, MRI **mri_fsurf, int runs, MRI *mri_dst)
{
  int  r ;
  MRI  *mri_cmat, *mri_run, *mri_template ;

  mri_template = MRIallocHeader(mri_aseg->width, mri_aseg->height,mri_aseg->depth,MRI_FLOAT,mri_fvol[0]->nframes);
  MRIcopyHeader(mri_aseg, mri_template) ;
  for (r = 0 ; r < runs ; r++)
  {
    if (MRImatchDimensions(mri_template, mri_fvol[r]) == 0)  // resample rsFMRI to match aseg
      mri_run = MRIresample(mri_fvol[r], mri_template, SAMPLE_TRILINEAR) ;
    else
      mri_run = mri_fvol[r] ;
    mri_cmat = compute_surface_correlations_in_aseg_label(mri_aseg, aseg_label, mri_run, mri_fsurf[r], NULL) ;
    if (mri_dst == NULL)
      mri_dst = mri_cmat ;
    else
    {
      mri_dst = MRIadd(mri_cmat, mri_dst, mri_dst) ;
      MRIfree(&mri_cmat) ;
    }
  }

  MRIscalarMul(mri_dst, mri_dst, 1.0/(float)runs) ; // make it an average
  MRIfree(&mri_template) ;
  if (write_diags && 0)
  {
    static int cno = 1 ;
    char fname[STRLEN] ;
    sprintf(fname, "c%d.mgz", cno++) ;
    printf("writing overlap maps to %s\n", fname) ;
    MRIwrite(mri_dst, fname) ;
  }

  return(mri_dst) ;
}
#endif


#if 0
static double
sample_label_avg(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_label_avg, double x, double y, double z, double *pvar)
{
  int   fno ;
  FACE  *face ;
  double fdist, val0, val1, val2, mean, var ;

  MHTfindClosestFaceGeneric(mht, mris_fixed, x, y, z, 8000, -1, 0, &face, &fno, &fdist);

  // interpolate mean
  val0 = MRIgetVoxVal(mri_label_avg, face->v[0], 0, 0,0) ;
  val1 = MRIgetVoxVal(mri_label_avg, face->v[1], 0, 0,0) ;
  val2 = MRIgetVoxVal(mri_label_avg, face->v[2], 0, 0,0) ;
  mean = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(mean)>100)
    DiagBreak() ;

  // interpolate variance
  val0 = MRIgetVoxVal(mri_label_avg, face->v[0], 0, 0,1) ;
  val1 = MRIgetVoxVal(mri_label_avg, face->v[1], 0, 0,1) ;
  val2 = MRIgetVoxVal(mri_label_avg, face->v[2], 0, 0,1) ;
  var = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(var)>100)
    DiagBreak() ;

  *pvar = var ;
  return(mean) ;
}
#endif


static double
sample_hemi_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, MRI *mri_label_avg, int vno, double x, double y, double z)
{
  double mean, var, val, error, total_error ;
  int    frame ;

  for (total_error = 0.0, frame = 0 ; frame < mri_label_avg->nframes ; frame++)
  {
    mean = sample_hemi_stats(mris_mov, mris_fixed, mht, mri_stats, x, y, z, frame, &var) ;
    val = MRIgetVoxVal(mri_label_avg, vno, 0, 0, frame) ;
    error = SQR(val-mean)/var ;
    total_error += error ;
  }
  
  return(total_error) ;
}

static int
compute_hemi_warp_gradient(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, MRI *mri_target_label, double dt)
{
  int    vno ;
  VERTEX *v ;
  double  e1p, e1m, e2p, e2m, x, y, z, dx, dy, dz ;

  for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;

    // sample vertex energy in plus and minus two principal directions
    x = v->x + v->e1x*dt ; y = v->y + v->e1y*dt ; z = v->z + v->e1z*dt ; 
    e1p  = sample_hemi_vertex_error(mris_mov, mris_fixed, mht, mri_stats, mri_target_label, vno, x, y, z) ;
    x = v->x - v->e1x*dt ; y = v->y - v->e1y*dt ; z = v->z - v->e1z*dt ; 
    e1m  = sample_hemi_vertex_error(mris_mov, mris_fixed, mht, mri_stats, mri_target_label, vno, x, y, z) ;
    x = v->x + v->e2x*dt ; y = v->y + v->e2y*dt ; z = v->z + v->e2z*dt ; 
    e2p  = sample_hemi_vertex_error(mris_mov, mris_fixed, mht, mri_stats, mri_target_label, vno, x, y, z) ;
    x = v->x - v->e2x*dt ; y = v->y - v->e2y*dt ; z = v->z - v->e2z*dt ; 
    e2m  = sample_hemi_vertex_error(mris_mov, mris_fixed, mht, mri_stats, mri_target_label, vno, x, y, z) ;
    if (e1p < e1m && e1p < e2p && e1p < e2m)   // move in e1p dir
    {
      dx = v->e1x*dt ; dy = v->e1y*dt ; dz = v->e1z*dt ;
    }
    else if (e1m < e2p && e1m < e2m)
    {
      dx = -v->e1x*dt ; dy = -v->e1y*dt ; dz = -v->e1z*dt ;
    }
    else if (e2m < e2p)
    {
      dx = -v->e2x*dt ; dy = -v->e2y*dt ; dz = -v->e2z*dt ;
    }
    else
    {
      dx = v->e2x*dt ; dy = v->e2y*dt ; dz = v->e2z*dt ;
    }
    v->dx += dx ; v->dy += dy ; v->dz += dz ;
  }

  return(NO_ERROR) ;
}

static double
sample_hemi_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, int frame, double *pvar)
{
  int   fno, var_offset = mri_stats->nframes/2 ;
  FACE  *face ;
  double fdist, val0, val1, val2, mean, var ;

  MHTfindClosestFaceGeneric(mht, mris_fixed, x, y, z, 8000, -1, 0, &face, &fno, &fdist);

  // interpolate mean
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,frame) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,frame) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,frame) ;
  mean = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(mean)>100)
    DiagBreak() ;

  // interpolate variance
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,frame+var_offset) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,frame+var_offset) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,frame+var_offset) ;
  var = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(var)>100)
    DiagBreak() ;

  *pvar = var ;
  return(mean) ;
}
static double
compute_hemi_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht_faces, MRI *mri_stats, MRI *mri_target_label, int nsubjects)
{
  int     vno, frame ;
  double  sse ;
  double  val, mean, var ;
  VERTEX  *v ;

  for (sse = 0.0, vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    for (frame = 0 ; frame < mri_target_label->nframes ; frame++)
    {
      mean = sample_hemi_stats(mris_mov, mris_fixed, mht_faces, mri_stats, v->x, v->y, v->z, frame, &var) ;
      val = MRIgetVoxVal(mri_target_label, vno, 0,0, frame) ;
      if (FZERO(var))
	var = 1.0 ;
      sse += SQR(mean-val)/var ;
    }
  }
  return(sqrt(sse/(mris_mov->nvertices*mri_target_label->nframes))) ;
}
static int
warp_hemi(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MRI *mri_target_label, MRI *mri_stats, int nsubjects, double tol, double dt) 
{
  MHT  *mht_vertices, *mht_faces ;
  int  iter, grad_averages, iters_so_far = 0 ;
  double error, last_error, pct_change, orig_dt ;

  mht_vertices = MHTcreateVertexTable_Resolution(mris_fixed, CURRENT_VERTICES, ceil(mris_fixed->avg_vertex_dist));
  mht_faces    = MHTcreateFaceTable_Resolution  (mris_fixed, CURRENT_VERTICES, ceil(mris_fixed->avg_vertex_dist));
#if 0
  if (write_diags)
  {
    char fname[STRLEN] ;
    MRI  *mri_fixed ; 

    mri_fixed = sample_fixed(mris_mov, mris_fixed, mht_faces, mri_stats) ;
    sprintf(fname, "lh.label_avg.%3.3d.mgz",  iter) ;
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_fixed, fname) ;
    MRIfree(&mri_fixed) ;
  }
#endif

  if (write_diags)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.%s.%3.3d", hemi, output_name, 0) ;
    printf("writing snapshot to %s\n", fname) ;
    MRISwrite(mris_mov, fname) ;
  }
  last_error = compute_hemi_error(mris_mov, mris_fixed, mht_faces, mri_stats, mri_target_label, nsubjects) ;
  orig_dt = dt ;
  printf("%3.3d: dt = %2.3f, error = %2.4f\n", iter=0, orig_dt, last_error) ;
  for (grad_averages = max_grad_averages ; grad_averages >= min_grad_averages ;grad_averages /= 4)
  {
    iter = 0 ;
    do
    {
      compute_hemi_warp_gradient(mris_mov, mris_fixed, mht_faces, mri_stats, mri_target_label, dt) ;
      MRISaverageGradients(mris_mov, grad_averages) ;
      MRISapplyGradient(mris_mov, 1) ; MRISprojectOntoSphere(mris_mov, mris_mov, mris_mov->radius) ;
      
      if (write_diags && (((iter+1)%10) == 0))
      {
	char fname[STRLEN] ;
	sprintf(fname, "%s.%s.%3.3d", hemi, output_name, iter+1+iters_so_far) ;
	printf("writing snapshot to %s\n", fname) ;
	MRISwrite(mris_mov, fname) ;
      }
      error = compute_hemi_error(mris_mov, mris_fixed, mht_faces, mri_stats, mri_target_label, nsubjects) ;
      pct_change = 100*(last_error - error) / last_error ;
      printf("%3.3d: dt = %2.3f, avgs = %d, error = %2.4f (%2.4f%%)\n", iter+1+iters_so_far, dt, grad_averages, error, pct_change) ;
      last_error = error ;
    } while (pct_change > tol && iter++ < max_iters) ;
    iters_so_far += iter ;
    dt /= 2 ;
    if (grad_averages == 0)
      break ;
  }
  
  
  MHTfree(&mht_faces) ; MHTfree(&mht_vertices) ;
  return(NO_ERROR) ;
}

#if  0
static MRI *
accumulate_subcortical_map(MRI **mri_fvol, MRI **mri_fsurf, LABEL *area, int runs, MRI *mri_stats, MRI *mri_mask)
{
  int  r ;
  MRI  *mri_cmat, *mri_avg = NULL, *mri_tmp ;

  // compute the average correlation within the label, then average that across runs
  for (r = 0 ; r < runs ; r++)
  {
    mri_cmat = compute_volume_correlations_in_surface_label(area, mri_fvol[r], mri_fsurf[r], NULL, mri_mask) ;
    mri_tmp = MRIaverageFrames(mri_cmat, NULL, 0, mri_cmat->nframes-1) ;
    if (Gx >= 0)
    {
      float val = MRIgetVoxVal(mri_tmp, Gx, Gy, Gz, 0) ;
      printf("frame avg (%d, %d, %d) = %2.3f\n", Gx, Gy, Gz, val) ;
    }
    MRIfree(&mri_cmat) ;
    if (mri_avg == NULL)
      mri_avg = mri_tmp ;
    else
    {
      MRIadd(mri_tmp, mri_avg, mri_avg) ;
      MRIfree(&mri_tmp) ;
    }
  }

  MRIscalarMul(mri_avg, mri_avg, 1.0/(float)runs) ; // make it an average
  if (Gx >= 0)
  {
    float val = MRIgetVoxVal(mri_avg, Gx, Gy, Gz, 0) ;
    printf("run avg (%d, %d, %d) = %2.3f\n", Gx, Gy, Gz, val) ;
  }

  // accumulate means and standard deviations into mri_stats (passed from caller)
  if (mri_stats == NULL)
  {
    mri_stats = MRIallocSequence(mri_avg->width, mri_avg->height, mri_avg->depth, MRI_FLOAT, 2) ;
    MRIcopyHeader(mri_avg, mri_stats) ;
  }
  MRIaddToFrame(mri_stats, mri_avg, mri_stats, 0, 0) ;
  MRIsqr(mri_avg, mri_avg) ;
  MRIaddToFrame(mri_stats, mri_avg, mri_stats, 1, 1) ;
  if (Gx >= 0)
  {
    float val1 = MRIgetVoxVal(mri_stats, Gx, Gy, Gz, 0) ;
    float val2 = MRIgetVoxVal(mri_stats, Gx, Gy, Gz, 1) ;
    printf("stats (%d, %d, %d) = %2.3f, %2.3f\n", Gx, Gy, Gz, val1, val2) ;
  }
  if (write_diags && 0)
  {
    static int cno = 1 ;
    char fname[STRLEN] ;
    sprintf(fname, "c%d.mgz", cno++) ;
    printf("writing overlap maps to %s\n", fname) ;
    MRIwrite(mri_avg, fname) ;
  }

  MRIfree(&mri_avg) ;
  return(mri_stats) ;

}
static int
compute_volume_correlation_at_vertex(MRI *mri_fvol, MRI *mri_fsurf, int vno, int frame, MRI *mri_cmat, MRI *mri_mask) 
{
  int   t, x, y, z, xmax, ymax, zmax ;
  float corr, sval, vval, norma, normb, max_corr ;
  static MRI   *mri_surf = NULL; 

  xmax = ymax = zmax = max_corr = 0.0 ;
  if (vno == Gdiag_no)
    DiagBreak() ;
  for (x = 0 ; x < mri_fvol->width ; x++)  
    for (y = 0 ; y < mri_fvol->height ; y++)
      for (z = 0 ; z < mri_fvol->depth ; z++)
      {
	if (mri_mask && MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
	  continue ;
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	for (norma = normb = corr = 0.0, t = 0 ;  t < mri_fsurf->nframes ; t++)
	{
	  sval = MRIgetVoxVal(mri_fsurf, vno, 0, 0, t) ;
	  vval = MRIgetVoxVal(mri_fvol, x, y, z, t) ;
	  corr += sval * vval ;
	  norma += sval*sval ;
	  normb += vval*vval ;
	}
	if (FZERO(norma) || FZERO(normb))
	  corr = 0 ;
	else
	  corr /= (sqrt(norma) * sqrt(normb)) ;
	if (corr > max_corr)
	{
	  max_corr = corr ; xmax = x ; ymax = y ; zmax = z ;
	}
	MRIsetVoxVal(mri_cmat, x, y, z, frame, corr) ;
      }

  if (write_diags > 1)
  {
    if (mri_surf == NULL)
      mri_surf = MRIallocSequence(mri_fsurf->width, 1, 1, MRI_FLOAT, 100) ;
    compute_surface_correlation_map_at_voxel(mri_fvol, mri_fsurf, xmax, ymax, zmax, frame, mri_surf) ;
    MRIwrite(mri_surf, "s.mgz") ;
  }

  return(NO_ERROR) ;
}
static MRI *
compute_volume_correlations_in_surface_label(LABEL *area, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat, MRI *mri_mask) 
{
  int  n ;

  mri_cmat = MRIallocSequence(mri_fvol->width, mri_fvol->height, mri_fvol->depth,MRI_FLOAT,area->n_points) ;
  MRIcopyHeader(mri_fvol, mri_cmat) ;

  // compute correlation maps for every point in the target label
  for (n = 0 ; n < area->n_points ; n++)
    compute_volume_correlation_at_vertex(mri_fvol, mri_fsurf, area->lv[n].vno, n, mri_cmat, mri_mask) ;

  return(mri_cmat) ;
}
#endif

VOXEL_LIST *
rank_voxels(MRI *mri_stats, int num_maps) 
{
  VOXEL_LIST *vl ;
  MRI        *mri_corr_coef ;

  mri_corr_coef = MRIdivideFrames(mri_stats, mri_stats, 0, 1, NULL) ;

  vl = VLSTcreate(mri_corr_coef, .00001, 10 , NULL, 0, 0) ;
  VLSTsort(vl, vl) ;
  MRIfree(&mri_corr_coef) ;
  return(vl) ;
}

static MRI *
compute_voxlist_surface_correlations(VOXEL_LIST *vl, int num_maps, MRI *mri_fvol, MRI *mri_fsurf, MRI *mri_cmat) 
{
  int  x, y, z, n ;

  mri_cmat = MRIallocSequence(mri_fsurf->width*mri_fsurf->height*mri_fsurf->depth,1,1,MRI_FLOAT,num_maps) ;

  // compute correlation maps for every point in the target label
  for (n = 0 ; n < num_maps ; n++)
  {
    x = vl->xi[n] ; y = vl->yi[n] ;z = vl->zi[n] ;
    if (x == Gx && y == Gy && z == Gz)
      DiagBreak() ;
    compute_surface_correlation_map_at_voxel(mri_fvol, mri_fsurf, x, y, z, n, mri_cmat) ;
  }

  return(mri_cmat) ;
}
static MRI *
compute_voxlist_surface_correlations_across_runs(VOXEL_LIST *vl, int num_maps, MRI **mri_fvol, MRI **mri_fsurf, int runs, MRI *mri_dst)
{
  int  r ;
  MRI  *mri_cmat ;

  for (r = 0 ; r < runs ; r++)
  {
    mri_cmat = compute_voxlist_surface_correlations(vl,  num_maps, mri_fvol[r], mri_fsurf[r], NULL) ;
    if (mri_dst == NULL)
      mri_dst = mri_cmat ;
    else
    {
      mri_dst = MRIadd(mri_cmat, mri_dst, mri_dst) ;
      MRIfree(&mri_cmat) ;
    }
  }

  MRIscalarMul(mri_dst, mri_dst, 1.0/(float)runs) ; // make it an average
  if (0 && write_diags) {
    static int cno = 1 ;
    char fname[STRLEN] ;
    sprintf(fname, "c%d.mgz", cno++) ;
    printf("writing overlap maps to %s\n", fname) ;
    MRIwrite(mri_dst, fname) ;
  }

  return(mri_dst) ;
}
static VECTOR *
compute_desired_map(MRI_SURFACE *mris, int cols,  LABEL *area, int ndilate) 
{
  VECTOR  *v ;
  int     n ;

  area = LabelCopy(area, NULL) ;
  v = VectorAlloc(cols, MATRIX_REAL) ;
  for (n = 0 ; n < area->n_points ; n++)
    VECTOR_ELT(v, area->lv[n].vno+1) = 1 ;

  if (ndilate > 0)
  {
    MRISclearMarks(mris) ; 
#if 1
    LabelDilate(area, mris, 1, CURRENT_VERTICES) ; // have a ring outside of label that isn't included in negative area
#endif
    LabelMark(area, mris) ;   
    MRIScopyMarkedToMarked2(mris) ;
    
    LabelDilate(area, mris, ndilate, CURRENT_VERTICES) ;
    LabelMark(area, mris) ;
    for (n = 0 ; n < mris->nvertices ; n++)
      if (mris->vertices[n].marked && mris->vertices[n].marked2 == 0)  // in ring around label
	VECTOR_ELT(v, n+1) = -.1/(float)ndilate ;
  }

  LabelFree(&area) ;  //  a local copy - not what was passed by caller
  return(v) ;
}

static double
compute_weight_functional(VECTOR **v_D, MATRIX **m_I, VECTOR *v_weights, int nsubjects)
{
  double total_error ;
  int    n ;
  VECTOR *v_weight_T = VectorTranspose(v_weights, NULL) ;
  
  total_error = 0.0 ;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate (v_D, v_weight_T) reduction(+:total_error) schedule(static,1)
#endif
  for (n = 0 ; n < nsubjects ; n++)
  {
    ROMP_PFLB_begin
    VECTOR *I_weighted, *I_weighted_T ;
    double error ;

    I_weighted = MatrixMultiply(v_weight_T, m_I[n], NULL) ;
    I_weighted_T = MatrixTranspose(I_weighted, NULL) ;

    error = VectorSSE(v_D[n], I_weighted_T) ;
    total_error += error ;
    VectorFree(&I_weighted) ;  VectorFree(&I_weighted_T) ; 
    ROMP_PFLB_end
  }
  ROMP_PF_end
  VectorFree(&v_weight_T) ;
  return(total_error) ;
}
static VECTOR *
compute_gradient_wrt_weights(VECTOR **v_D, MATRIX **m_I, VECTOR *v_weights, int nsubjects, VECTOR *v_gradient)
{
  int    n, l, nvertices, nvox ;
  VECTOR *v_weight_T = VectorTranspose(v_weights, NULL), *v_diff[MAX_SUBJECTS] ;
  
  if (v_gradient == NULL)
    v_gradient = VectorClone(v_weights) ;
  else
    VectorClear(v_gradient) ;
  nvertices = v_D[0]->rows ; nvox = v_weights->rows ;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) firstprivate (v_D) shared(v_diff)  schedule(static,1)
#endif
  for (n = 0 ; n < nsubjects ; n++)
  {
    ROMP_PFLB_begin
    VECTOR *I_weighted, *I_weighted_T ;
    I_weighted = MatrixMultiply(v_weight_T, m_I[n], NULL) ;
    I_weighted_T = VectorTranspose(I_weighted, NULL) ;
    v_diff[n] = MatrixSubtract(v_D[n], I_weighted_T, NULL) ;
    VectorFree(&I_weighted) ; VectorFree(&I_weighted_T) ; 
    ROMP_PFLB_end
  }
  ROMP_PF_end
  VectorFree(&v_weight_T) ;

  for (l = 0 ; l < nvox ; l++)
  {
    ROMP_PF_begin
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(experimental) firstprivate (l, v_gradient, nvertices) schedule(static,1)
#endif
    for (n = 0 ; n < nsubjects ; n++)
    {
      ROMP_PFLB_begin
      int i ;
      for (i = 0 ; i < nvertices ; i++)
      {
	double  vd, In;
	vd = VECTOR_ELT(v_diff[n], i+1) ;
	In = *MATRIX_RELT(m_I[n], l+1, i+1) ;
	VECTOR_ELT(v_gradient, l+1) += vd * (-In) ;
      }
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }
  for (n = 0 ; n < nsubjects ; n++)
    VectorFree(&v_diff[n]) ;
  VectorScalarMul(v_gradient, 1.0/((float)nvertices*(float)nvox), v_gradient) ;
  return(v_gradient) ;
}
static int max_wt_iters = 100000 ;
static double wt_dt = .9 ;
static double wt_mom = .75 ;

static VECTOR *
optimize_weight_functional(VECTOR **v_D, MATRIX **m_I, int nsubjects, MRI *mri_mask, char *out_prefix)
{
  MATRIX *v_weights ;
  double  last_error, error, pct_error ;
  int     iter ;
  VECTOR  *v_d_weights = NULL ;
  VECTOR  *v_new_weights = NULL ;
  FILE    *log_fp ;
  char    fname[STRLEN] ;
  VECTOR  *v_last_delta ;

  v_weights = VectorAlloc(m_I[0]->rows, MATRIX_REAL) ;
  v_d_weights = VectorClone(v_weights) ;
  v_last_delta = VectorClone(v_weights) ;
  last_error = error = compute_weight_functional(v_D, m_I, v_weights, nsubjects) ;

  sprintf(fname, "%s.log", out_prefix) ;
  log_fp = fopen(fname, "w") ;
  printf("%4.4d: error = %2.3f\n", 0, error) ;
  fprintf(log_fp, "%4.4d: error = %2.3f\n", 0, error) ; fflush(log_fp) ;
  for (iter = 0 ; iter < max_wt_iters ; iter++)
  {
    last_error = error ;
    v_d_weights = compute_gradient_wrt_weights(v_D, m_I, v_weights, nsubjects, v_d_weights) ;
    VectorScalarMul(v_d_weights, -wt_dt, v_d_weights) ;
    VectorScalarMul(v_last_delta, wt_mom, v_last_delta) ;
    VectorAdd(v_d_weights, v_last_delta, v_d_weights) ;
    VectorCopy(v_d_weights, v_last_delta) ;
    v_new_weights = VectorAdd(v_d_weights, v_weights, v_new_weights) ;
    error = compute_weight_functional(v_D, m_I, v_new_weights, nsubjects) ;
    pct_error = 100*(last_error - error) / (last_error) ;
    if ((iter+1)%100 == 0)
    {
      printf("%4.4d: error = %2.3f (%2.3f%%)\n", iter+1, error, pct_error) ;
      fprintf(log_fp, "%4.4d: error = %2.3f (%2.3f%%)\n", iter+1, error, pct_error) ; fflush(log_fp) ;
    }
    if (write_diags && ((iter+1)%100 == 0))
      write_snapshot(v_weights, m_I, mri_mask, out_prefix, iter+1, nsubjects) ;
    MatrixCopy(v_new_weights, v_weights) ;
  }

  VectorFree(&v_d_weights) ;
  return(v_weights) ;
}
static MATRIX *
compute_subcortical_map_weights(MRI_SURFACE *mris, MRI *mri_fvol[MAX_SUBJECTS][MAX_RUNS], MRI *mri_fsurf[MAX_SUBJECTS][MAX_RUNS], LABEL **labels, int runs, int nsubjects, MRI *mri_mask, int ndilate, char *out_prefix) 
{
  int        n ;
  MATRIX     *m_I[MAX_SUBJECTS] ;
  VECTOR     *v_weights, *v_D[MAX_SUBJECTS] ;
  VOXEL_LIST *vl ;
  MRI        *mri_cmat ;

  /* compute  array of surface maps.
     m_I = nvox x nvertices arrays of surface maps
     v_D = desried map (1 in label, 0 outside, or maybe -1 at in dilated region
  */
  vl = VLSTcreate(mri_mask, .001, 256 , NULL, 0, 0) ;
  if (write_diags)
  {
    char fname[STRLEN] ;
    MRI *mri_order = VLSTwriteOrderToMRI(vl, NULL) ;
    sprintf(fname, "%s.order.mgz", out_prefix) ;
    printf("writing voxel list order to %s\n", fname) ;
    MRIwrite(mri_order, fname) ;
    MRIfree(&mri_order) ;
  }

//  vl->nvox = 1000;

  printf("creating target maps...\n") ;
  for (n = 0 ; n <= nsubjects ; n++)
  {
    v_D[n] = compute_desired_map(mris, mri_fsurf[0][0]->width,  labels[n], ndilate) ;
    if (write_diags)
    {
      char fname[STRLEN] ;
      int  vno ;
      float m ;
      MRISclearCurvature(mris) ;  
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
	if (vno == Gdiag_no)
	  DiagBreak() ;
	m = VECTOR_ELT(v_D[n], vno+1) ;
	mris->vertices[vno].curv = m ;
      }
      sprintf(fname, "%s.%s.desired.%s.mgz", hemi, out_prefix, subjects[n]) ;
      printf("writing target map to %s\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
    }
  }
  printf("reading timeseries and creating correlations\n") ;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  mri_cmat = NULL ;
  #pragma omp parallel for if_ROMP(experimental) firstprivate (mri_cmat) shared(m_I, vl, mri_fvol, mri_fsurf, runs) schedule(static,1)
#endif
  for (n = 0 ; n <= nsubjects ; n++)
  {
    ROMP_PFLB_begin
    mri_cmat = compute_voxlist_surface_correlations_across_runs(vl, vl->nvox, mri_fvol[n], mri_fsurf[n], runs, NULL) ;
    m_I[n] = MRIcopyFramesToMatrixRows(mri_cmat, NULL, 0, vl->nvox, 1) ;
    if (write_diags)
    {
      char fname[STRLEN] ;

      sprintf(fname, "%s.%s.cmat.%s.mgz", hemi, out_prefix, subjects[n]) ;
      printf("writing subject cmat to %s\n", fname) ;
      MRIwrite(mri_cmat, fname) ;
    }
    MRIfree(&mri_cmat) ;
    ROMP_PFLB_end
  }
  ROMP_PF_end

  if (create_only)
    exit(0) ;
  printf("creating linear combination of %d subcortical maps\n", vl->nvox) ;
  if (use_powell)
  {
    v_weights = VectorAlloc(m_I[0]->rows, MATRIX_REAL) ;
    powell_minimize(v_D, m_I, nsubjects, mri_mask, v_weights) ;
  }
  else
    v_weights = optimize_weight_functional(v_D, m_I, nsubjects, mri_mask, out_prefix) ;
  if (write_diags)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.weights.txt", out_prefix) ;
    MatrixWriteTxt(fname, v_weights) ;
  }
  VLSTfree(&vl) ;
  for (n = 0 ; n <= nsubjects ; n++)
  {
    MatrixFree(&m_I[n]) ;
    VectorFree(&v_D[n]) ;
  }

  if (write_diags)
  {
    MRI  *mri_wts, *mri_surf ;
    VECTOR *I_weighted, *v_weights_T, *I_map ;

    mri_wts = MRIclone(mri_mask, NULL);
    mri_surf = MRIalloc(mri_fsurf[0][0]->width,  mri_fsurf[0][0]->height,  mri_fsurf[0][0]->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_fsurf[0][0], mri_surf) ;

    v_weights_T = VectorTranspose(v_weights, NULL) ;
    for (n = 0 ; n < vl->nvox ; n++)
    {
      MRIsetVoxVal(mri_wts, vl->xi[n], vl->yi[n], vl->zi[n], 0, VECTOR_ELT(v_weights, n+1)) ;
      I_weighted = MatrixMultiply(v_weights_T, m_I[n], NULL) ;
      if (I_map == NULL)
	I_map = VectorClone(I_weighted) ;
      VectorAdd(I_weighted, I_map, I_map) ;
      VectorFree(&I_weighted) ;
    }

    for (n = 0 ; n < I_map->cols ; n++)
      MRIsetVoxVal(mri_surf, n, 0, 0, 0, VECTOR_ELT(I_map, n+1)) ;
    MRIwrite(mri_wts, "w.mgz") ;
    MRIwrite(mri_surf, "s.mgz") ;

    VectorFree(&v_weights_T) ; VectorFree(&I_map) ;
    MRIfree(&mri_wts) ; MRIfree(&mri_surf) ;
  }

  return(v_weights) ;
}

static int
write_snapshot(VECTOR *v_weights, MATRIX **m_I, MRI *mri_mask, const char *prefix, int iter, int nsubjects)
{
  MRI         *mri_wts, *mri_surf ;
  VECTOR      *v_weights_T, *I_map ;
  char        fname[STRLEN] ;
  int         n, i ;
  VOXEL_LIST  *vl ;

  vl = VLSTcreate(mri_mask, .001, 256 , NULL, 0, 0) ;

  mri_wts = MRIclone(mri_mask, NULL);
  mri_surf = MRIalloc(m_I[0]->cols,  1,  1, MRI_FLOAT) ;
  
  v_weights_T = VectorTranspose(v_weights, NULL) ;
  for (n = 0 ; n <= nsubjects ; n++)
  {
    I_map = MatrixMultiply(v_weights_T, m_I[n], NULL) ;
    for (i = 0 ; i < I_map->cols ; i++)
      MRIsetVoxVal(mri_surf, i, 0, 0, 0, RVECTOR_ELT(I_map, i+1)) ;
    sprintf(fname, "%s.%s.%3.3d.map.mgz", prefix, subjects[n], iter) ; printf("writing map to %s\n", fname) ;
    MRIwrite(mri_surf, fname) ;
    VectorFree(&I_map) ;
  }

  for (n = 0 ; n < m_I[0]->rows ; n++)
    MRIsetVoxVal(mri_wts, vl->xi[n], vl->yi[n], vl->zi[n], 0, VECTOR_ELT(v_weights, n+1)) ;
  
  sprintf(fname, "%s.%3.3d.wts.mgz", prefix, iter) ;   printf("writing wt image to %s\n", fname) ;
  MRIwrite(mri_wts, fname) ;
  sprintf(fname, "%s.%3.3d.wts.txt", prefix, iter) ;   printf("writing wt vector to %s\n", fname) ;
  MatrixWriteTxt(fname, v_weights) ;
    
  VectorFree(&v_weights_T) ; MRIfree(&mri_wts) ; MRIfree(&mri_surf) ; VLSTfree(&vl) ;
  return(NO_ERROR) ;
}
#include "numerics.h"
static float compute_powell_sse(float *p) ;

#ifdef TOL
#undef TOL
#endif
#define TOL 1e-5

static VECTOR *G_v_weights ;
static int G_nsubjects ;
static MATRIX **G_m_I ;
static VECTOR **G_v_D ;
static MRI *G_mri_mask ;


static float
compute_powell_sse(float *p)
{
  int   i ;

  for (i = 1 ; i <= G_v_weights->rows ;  i++)
    VECTOR_ELT(G_v_weights, i) = p[i] ;
  return(compute_weight_functional(G_v_D, G_m_I, G_v_weights, G_nsubjects)) ;
}


static int max_powell_iter = 10 ;
static double powell_ftol = .001 ;
static double powell_linmin_tol = .001 ;

static int
powell_minimize(VECTOR **v_D, MATRIX **m_I, int nsubjects, MRI *mri_mask, VECTOR *v_weights)
{
  float *p, **xi, fret, fstart, min_sse;
  int   i, r, c, iter ;

  G_v_D = v_D ;
  G_m_I = m_I ;
  G_nsubjects = nsubjects ;
  G_v_weights = v_weights ;
  G_mri_mask = mri_mask ;

  p = vector(1, v_weights->rows) ;
  xi = matrix(1, v_weights->rows, 1, v_weights->rows) ;
  for (i = 1 ; i <= v_weights->rows ; i++)
    p[i] = VECTOR_ELT(v_weights, i) ;

  for (r = 1 ; r <= v_weights->rows ; r++)
  {
    for (c = 1 ; c <= v_weights->rows ; c++)
    {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  min_sse = compute_powell_sse(p) ;
//  Gdiag_no = 1 ;
  OpenPowell2(p, xi, v_weights->rows, powell_ftol, powell_linmin_tol, max_powell_iter, &iter, &fret, compute_powell_sse);
  for (i = 1 ; i <= v_weights->rows ; i++)
    VECTOR_ELT(v_weights, i)  = p[i] ;
  if (write_diags)
    write_snapshot(v_weights, m_I, mri_mask, "after_powell", iter+1, nsubjects) ;
  do
  {
    // reinitialize powell directions
    for (r = 1 ; r <= v_weights->rows ; r++)
    {
      for (c = 1 ; c <= v_weights->rows ; c++)
      {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell2(p, xi, v_weights->rows, powell_ftol, powell_linmin_tol, max_powell_iter, &iter, &fret, compute_powell_sse);
    for (i = 1 ; i <= v_weights->rows ; i++)
      VECTOR_ELT(v_weights, i) = p[i] ;
    if (write_diags)
      write_snapshot(v_weights, m_I, mri_mask, "after_second_powell", iter+1, nsubjects) ;

    if (write_diags)
      printf("best alignment after powell: %2.3f (%d steps)\n", fret, iter) ;
    if ((fstart-fret)/fstart < TOL)
      break ;
  }
  while (fret < fstart) ;

  free_matrix(xi, 1, v_weights->rows, 1, v_weights->rows) ;
  free_vector(p, 1, v_weights->rows) ;
  return(NO_ERROR) ;
}
