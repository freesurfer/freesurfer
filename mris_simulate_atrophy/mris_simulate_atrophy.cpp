/**
 * @brief simulate cortical atrophy
  *
 * program to simulate atrophic changes in the cortical by
 * darkening the T1 intensities according to decrease volume fractions
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "transform.h"
#include "label.h"

#define WM_VAL         1
#define GM_VAL         2
#define CSF_VAL        3
#define SUBCORT_GM_VAL 4

static MRI *MRISsimulateAtrophy(MRI *mri_norm, MRI *mri_unpv_intensities, MRI  *mri_wm, MRI *mri_subcort_gm, MRI *mri_cortex, MRI *mri_csf,
				LABEL *area, double atrophy_frac, MRI *mri_norm_atrophy, MRI **pmri_cortex_out, MRI **pmri_csf_out) ;
static void patch_csf_vol(MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort,  MRI *mri_vfrac_csf) ;
MRI *add_aseg_structures_outside_ribbon(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst,
                                       int wm_val, int gm_val, int csf_val) ;
int MRIcomputePartialVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox, 
                                     MRI *mri_seg, MRI *mri_wm, MRI *mri_subcort_gm, MRI *mri_cortex, 
				     MRI *mri_csf,
                                     int wm_val, int subcort_gm_val, int cortex_val, int csf_val) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
MRI *MRIsimulateAtrophy(MRI *mri_norm, MRI *mri_aseg,  int target_label, 
                        int *border_labels, int nlabels, float atrophy_pct, 
                        MRI *mri_norm_atrophy)  ;

static MRI *compute_unpartial_volumed_intensities(MRI *mri_src, MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort, 
						 MRI *mri_vfrac_csf, 
						  int whalf0,  double sigma, MRI *mri_dst, int separate_frames) ;

static float noise_sigma = 4 ;
static char *sdir = NULL ;

static const char *T1_name = "brain.finalsurfs.mgz" ;
static const char *white_name = "white" ;
static const char *pial_name = "pial" ;
static const char *aseg_name = "aseg.mgz" ;
static float resolution = .5 ;

static float noise_min = -1 ;
static float noise_max = 0 ;
static float noise_step = 0 ;

static float atrophy_min = -1 ;
static float atrophy_max = 0 ;
static float atrophy_step = 0 ;

static int whalf = 1 ;
static double sigma = 1 ;

int
main(int argc, char *argv[])
{
  char        **av, *out_fname, *subject, *hemi, buf[STRLEN], fname[STRLEN] ;
  int          ac, nargs, msec, minutes, seconds, nvox ;
  Timer start ;
  MRI         *mri_norm, *mri_norm_atrophy, *mri_noise, *mri_pial, *mri_seg, *mri_aseg, *mri_cortex,*mri_wm,
    *mri_csf, *mri_subcort_gm, *mri_ribbon, *mri_unpv_intensities, *mri_tmp, *mri_noisy_atrophy, *mri_cortex_atrophy, *mri_csf_atrophy ;
  LABEL       *area ;
  float        atrophy_frac ;
  MRI_SURFACE  *mris_white_lh, *mris_pial_lh, *mris_white_rh, *mris_pial_rh ;
  MATRIX       *m_vox2vox ;
  char         extension[STRLEN] ;

  nargs = handleVersionOption(argc, argv, "mris_simulate_atrophy");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit(1) ;

  if (sdir == NULL)
  {
    char *cp = getenv("SUBJECTS_DIR") ;
    if (cp == NULL)
      ErrorExit(ERROR_UNSUPPORTED, "%s: SUBJECTS_DIR must be specified in env or on cmdline with -sdir", Progname) ;
    sdir = buf ;
    strcpy(buf, cp) ;
  }
  subject = argv[1] ; hemi = argv[2] ;
  sprintf(fname, "%s/%s/label/%s.%s", sdir, subject, hemi, argv[3]) ;
  area = LabelRead(NULL, fname) ;
  if (area == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not load label %s", Progname, fname) ;
  atrophy_frac = atof(argv[4]) ;
  if (atrophy_frac > 1)
  {
    if (atrophy_frac > 100)
      ErrorExit(ERROR_BADPARM, "%s: must specify atrophy in [0 1]",
		Progname) ;
    atrophy_frac /= 100 ; // assume it was specified in pct
  }
  out_fname = argv[5] ;
  FileNameExtension(out_fname, extension) ;
  FileNameRemoveExtension(out_fname, out_fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, T1_name) ;
  mri_norm = MRIread(fname) ;
  if (mri_norm == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/mri/%s", sdir, subject, aseg_name) ;
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/rh.%s", sdir, subject, white_name) ;
  mris_white_rh = MRISread(fname) ;
  if (mris_white_rh == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh white surface from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/lh.%s", sdir, subject, white_name) ;
  mris_white_lh = MRISread(fname) ;
  if (mris_white_lh == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh white surface from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/rh.%s", sdir, subject, pial_name) ;
  mris_pial_rh = MRISread(fname) ;
  if (mris_pial_rh == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read rh pial surface from %s", Progname, fname) ;


  sprintf(fname, "%s/%s/surf/lh.%s", sdir, subject, pial_name) ;
  mris_pial_lh = MRISread(fname) ;
  if (mris_pial_lh == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read lh pial surface from %s", Progname, fname) ;

  if (strcmp(hemi, "rh") == 0)
    LabelToCurrent(area, mris_pial_rh) ;
  else
    LabelToCurrent(area, mris_pial_lh) ;

  if (noise_min < 0)  // only one step
  {
    noise_min = noise_max = noise_sigma ; noise_step = 1 ;
  }
  if (atrophy_min < 0)  // only one step
  {
    atrophy_min = atrophy_max = atrophy_frac ; atrophy_step = 1 ;
  }

  nvox = (int)ceil(mri_norm->width/resolution); 
  mri_pial = MRIalloc(nvox, nvox, nvox, MRI_UCHAR) ;
  MRIsetResolution(mri_pial, resolution, resolution, resolution) ;

  mri_pial->xstart = -resolution*mri_pial->width/2.0 ;
  mri_pial->xend = resolution*mri_pial->width/2.0 ;
  mri_pial->ystart = -resolution*mri_pial->height/2.0 ;
  mri_pial->yend = resolution*mri_pial->height/2.0 ;
  mri_pial->zstart = -resolution*mri_pial->depth/2.0 ;
  mri_pial->zend = resolution*mri_pial->depth/2 ;
  mri_pial->c_r = mri_aseg->c_r ; mri_pial->c_a = mri_aseg->c_a ; mri_pial->c_s = mri_aseg->c_s ;
  MRIreInitCache(mri_pial) ; 

  printf("filling interior of lh pial surface...\n") ;
  MRISfillInterior(mris_pial_lh, resolution, mri_pial) ;
  mri_seg = MRIclone(mri_pial, NULL) ;
  mri_tmp = MRIclone(mri_pial, NULL) ;
  printf("filling interior of rh pial surface...\n") ;
  MRISfillInterior(mris_pial_rh, resolution, mri_tmp) ;
  MRIcopyLabel(mri_tmp, mri_pial, 1) ;
  MRIclear(mri_tmp) ;
  printf("filling interior of lh white matter surface...\n") ;
  MRISfillWhiteMatterInterior(mris_white_lh, mri_aseg, mri_seg, resolution,
                              WM_VAL, SUBCORT_GM_VAL, CSF_VAL);
  printf("filling interior of rh white matter surface...\n") ;
  MRISfillWhiteMatterInterior(mris_white_rh, mri_aseg, mri_tmp, resolution,
                              WM_VAL, SUBCORT_GM_VAL, CSF_VAL);
  MRIcopyLabel(mri_tmp, mri_seg, WM_VAL) ;
  MRIcopyLabel(mri_tmp, mri_seg, SUBCORT_GM_VAL) ;
  MRIcopyLabel(mri_tmp, mri_seg, CSF_VAL) ;
  MRIfree(&mri_tmp) ;

  mri_ribbon = MRInot(mri_seg, NULL) ;
  MRIcopyLabel(mri_seg, mri_pial, CSF_VAL) ;
  MRIreplaceValuesOnly(mri_pial, mri_pial, CSF_VAL, 0) ;
  MRIand(mri_ribbon, mri_pial, mri_ribbon, 1) ;
  MRIbinarize(mri_ribbon, mri_ribbon, 1, 0, GM_VAL) ;
  MRIcopyLabel(mri_ribbon, mri_seg, GM_VAL) ;
  MRIreplaceValuesOnly(mri_seg, mri_seg, CSF_VAL, 0) ;
  add_aseg_structures_outside_ribbon(mri_seg, mri_aseg, mri_seg, WM_VAL, SUBCORT_GM_VAL, CSF_VAL) ;

  m_vox2vox = MRIgetVoxelToVoxelXform(mri_seg, mri_norm) ;
  mri_cortex = MRIalloc(mri_norm->width, mri_norm->height, mri_norm->depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_norm, mri_cortex) ;
  mri_subcort_gm = MRIclone(mri_cortex, NULL) ;
  mri_wm = MRIclone(mri_cortex, NULL) ;
  mri_csf = MRIclone(mri_cortex, NULL) ;
  printf("computing partial volume fractions...\n") ;
  MRIcomputePartialVolumeFractions(mri_norm, m_vox2vox, mri_seg, mri_wm, mri_subcort_gm, mri_cortex, mri_csf,
				   WM_VAL, SUBCORT_GM_VAL, GM_VAL, 0) ;
  if (Gdiag & DIAG_WRITE)
  {
    printf("writing volume fractions...\n") ;
    MRIwrite(mri_wm, "wm.vfrac.mgz") ;
    MRIwrite(mri_csf, "csf.vfrac.mgz") ;
    MRIwrite(mri_cortex, "cortex.vfrac.mgz") ;
    MRIwrite(mri_subcort_gm, "gm.vfrac.mgz") ;
  }
  patch_csf_vol(mri_wm, mri_cortex, mri_subcort_gm,  mri_csf) ;
  mri_unpv_intensities =   
    compute_unpartial_volumed_intensities(mri_norm, mri_wm,  mri_cortex, mri_subcort_gm, 
					  mri_csf, whalf,  sigma, NULL, 1) ;
  if (Gdiag & DIAG_WRITE)
    MRIwrite(mri_unpv_intensities, "pvi.mgz") ;

  for (atrophy_frac = atrophy_min ; atrophy_frac <= atrophy_max ; atrophy_frac += atrophy_step)
  {
    mri_norm_atrophy =  MRISsimulateAtrophy(mri_norm, mri_unpv_intensities, mri_wm, mri_subcort_gm, mri_cortex, mri_csf,
					    area, atrophy_frac, NULL, &mri_cortex_atrophy, &mri_csf_atrophy) ;
    
    sprintf(fname, "%s.gm.atrophy%2.1f.%s", out_fname, atrophy_frac, extension) ;
    printf("writing atrophic gm vfracs to %s\n", fname) ;
    MRIwrite(mri_cortex_atrophy, fname) ;
    sprintf(fname, "%s.csf.atrophy%2.1f.%s", out_fname, atrophy_frac, extension) ;
    printf("writing atrophic csf vfracs to %s\n", fname) ;
    MRIwrite(mri_csf_atrophy, fname) ;
    sprintf(fname, "%s.atrophy.0.0.noise.0.0.%s", out_fname, extension) ;
    printf("writing simulated atrophy with noise sigma = 0 image to %s\n", fname) ;
    MRIwrite(mri_norm_atrophy, fname) ;

    for (mri_noisy_atrophy = NULL, noise_sigma = noise_min ; noise_sigma <= noise_max ; noise_sigma += noise_step)
    {
      mri_noise = MRIrandn(mri_norm->width, mri_norm->height, mri_norm->depth, 1, 0, noise_sigma, NULL) ;
      MRImaskZero(mri_noise, mri_norm, mri_noise) ;
      mri_noisy_atrophy = MRIadd(mri_noise, mri_norm_atrophy, mri_noisy_atrophy) ;
      MRIfree(&mri_noise) ;
      
      sprintf(fname, "%s.atrophy.%2.2f.noise.%2.1f.%s", out_fname, atrophy_frac, noise_sigma, extension) ;
      printf("writing simulated atrophy (%2.2f) with noise sigma = %2.1f image to %s\n", atrophy_frac, noise_sigma, fname) ;
      MRIwrite(mri_noisy_atrophy, fname) ;
    }
#if 0
    {
      char extension[STRLEN], fname2[STRLEN] ;
      FileNameExtension(out_fname, extension) ;
      FileNameRemoveExtension(out_fname, out_fname) ;
      sprintf(fname2, "%s.synth.%s", out_fname, extension) ;
      printf("writing synthesized volume with no atrophy to %s\n", fname2) ;
      MRIwrite(mri_no_atrophy, fname2) ;
    }
#endif
  }
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr,
          "simulation took %d minutes and %d seconds.\n", minutes, seconds) ;

  exit(0) ;
  return(0) ;
}



static int get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    sdir = argv[2] ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "nrange") || !stricmp(option, "noise_range"))
  {
    if (sscanf(argv[2], "%f:%f:%f", &noise_min, &noise_step, &noise_max) != 3)
      ErrorExit(ERROR_BADPARM, "%s: couldn't parse min:step:max from '%s'", Progname, argv[2]) ;
    printf("stepping through noise values from %2.1f to %2.1f in steps of %2.2f\n", 
	   noise_min, noise_max, noise_step) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "arange") || !stricmp(option, "atrophy_range"))
  {
    if (sscanf(argv[2], "%f:%f:%f", &atrophy_min, &atrophy_step, &atrophy_max) != 3)
      ErrorExit(ERROR_BADPARM, "%s: couldn't parse min:step:max from '%s'", Progname, argv[2]) ;
    printf("stepping through atrophy values from %2.1f to %2.1f in steps of %2.2f\n", 
	   atrophy_min, atrophy_max, atrophy_step) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "T1") || !stricmp(option, "input"))
  {
    T1_name = argv[2] ;
    printf("using %s as input volume\n", T1_name) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
  case 'W':
    whalf = atoi(argv[2]) ;
    printf("using half window size = %d\n", whalf) ;
    if (whalf < 1)
      ErrorExit(ERROR_BADPARM, "half window size must be >= 1 to estimate unpartial-volumed intensities") ;
    nargs = 1 ;
    break ;
  case 'S':
    sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting Gaussian smoothing sigma to %2.3fmm\n", sigma) ;
    break ;
  case 'N':
    noise_sigma = atof(argv[2]) ;
    printf("using noise level %f\n", noise_sigma) ;
    nargs = 1 ;
    break ;
  case 'R':
    resolution = atof(argv[2]) ;
    printf("using resolution = %2.3f\n", resolution) ;
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

static void usage_exit(int code)
{
  printf("usage: %s [options] <subject> <hemi> <label> <atrophy fraction> <output volume>\n",
         Progname) ;
  printf("\noptions:\n");
  printf("  -a <atrophy %%> - %% atrophy to simulate of structure\n");
  printf("  -n <sigma>     - gaussian noise level to add\n") ;
  exit(code) ;
}
int
MRIcomputePartialVolumeFractions(MRI *mri_src, MATRIX *m_vox2vox, 
                                 MRI *mri_seg, MRI *mri_wm, MRI *mri_subcort_gm, MRI *mri_cortex,
				 MRI *mri_csf,
                                 int wm_val, int subcort_gm_val, int cortex_val, int csf_val)
{
  int    x, y, z, xs, ys, zs, label ;
  VECTOR *v1, *v2 ;
  MRI    *mri_counts ;
  float  val, count ;
  MATRIX *m_inv ;

  m_inv = MatrixInverse(m_vox2vox, NULL) ;
  if (m_inv == NULL)
  {
    MatrixPrint(stdout, m_vox2vox) ;
    ErrorExit(ERROR_BADPARM, "MRIcomputePartialVolumeFractions: non-invertible vox2vox matrix");
  }
  mri_counts = MRIcloneDifferentType(mri_src, MRI_INT) ;

  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  for (x = 0 ; x < mri_seg->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_seg->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
        if (xs >= 0 && ys >= 0 && zs >= 0 &&
            xs < mri_src->width && ys < mri_src->height && zs < mri_src->depth)
        {
          val = MRIgetVoxVal(mri_counts, xs, ys, zs, 0) ;
          MRIsetVoxVal(mri_counts, xs, ys, zs, 0, val+1) ;

          label = MRIgetVoxVal(mri_seg, x, y, z, 0) ;
          if (label == csf_val)
          {
            val = MRIgetVoxVal(mri_csf, xs, ys, zs, 0) ;
            MRIsetVoxVal(mri_csf, xs, ys, zs, 0, val+1) ;
          }
          else if (label == wm_val)
          {
            val = MRIgetVoxVal(mri_wm, xs, ys, zs, 0) ;
            MRIsetVoxVal(mri_wm, xs, ys, zs, 0, val+1) ;
          }
          else if (label == subcort_gm_val)
          {
            val = MRIgetVoxVal(mri_subcort_gm, xs, ys, zs, 0) ;
            MRIsetVoxVal(mri_subcort_gm, xs, ys, zs, 0, val+1) ;
          }
          else if (label == cortex_val)
          {
            val = MRIgetVoxVal(mri_cortex, xs, ys, zs, 0) ;
            MRIsetVoxVal(mri_cortex, xs, ys, zs, 0, val+1) ;
          }
          else
            DiagBreak() ;
        }
      }
    }
  }

  for (x = 0 ; x < mri_src->width ; x++)
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        count = MRIgetVoxVal(mri_counts, x, y, z, 0) ;
        if (count >= 1)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          val = MRIgetVoxVal(mri_wm, x, y, z, 0) ;
          MRIsetVoxVal(mri_wm, x, y, z, 0, val/count) ;
          val = MRIgetVoxVal(mri_subcort_gm, x, y, z, 0) ;
          MRIsetVoxVal(mri_subcort_gm, x, y, z, 0, val/count) ;
          val = MRIgetVoxVal(mri_cortex, x, y, z, 0) ;
          MRIsetVoxVal(mri_cortex, x, y, z, 0, val/count) ;
          val = MRIgetVoxVal(mri_csf, x, y, z, 0) ;
          MRIsetVoxVal(mri_csf, x, y, z, 0, val/count) ;
        }
        else  // sample in other direction
        {
          V3_X(v1) = x ; V3_Y(v1) = y ; V3_Z(v1) = z ;
          MatrixMultiply(m_inv, v1, v2) ;
          MatrixMultiply(m_inv, v1, v2) ;
          xs = nint(V3_X(v2)) ; ys = nint(V3_Y(v2)) ; zs = nint(V3_Z(v2)) ;
          if (xs >= 0 && ys >= 0 && zs >= 0 &&
              xs < mri_seg->width && ys < mri_seg->height && zs < mri_seg->depth)
          {
            label = MRIgetVoxVal(mri_seg, xs, ys, zs, 0) ;
            if (label == csf_val)
              MRIsetVoxVal(mri_csf, x, y, z, 0, 1) ;
            else if (label == wm_val)
              MRIsetVoxVal(mri_wm, x, y, z, 0, 1) ;
            else if (label == subcort_gm_val)
              MRIsetVoxVal(mri_subcort_gm, x, y, z, 0, 1) ;
            else if (cortex_val)
              MRIsetVoxVal(mri_cortex, x, y, z, 0, 1) ;
            else
              DiagBreak() ;
          }
        }
      }
  VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_inv) ;
  MRIfree(&mri_counts) ;

  return(NO_ERROR) ;
}
MRI *
add_aseg_structures_outside_ribbon(MRI *mri_src, MRI *mri_aseg, MRI *mri_dst,
                                   int wm_val, int gm_val, int csf_val)
{
  VECTOR *v1, *v2 ;
  MATRIX *m_vox2vox ;
  int    x, y, z, xa, ya, za, seg_label, aseg_label ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_src, mri_aseg) ;

  for (x = 0 ; x < mri_dst->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_dst->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        seg_label = nint(MRIgetVoxVal(mri_dst, x, y, z, 0)) ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xa = (int)(nint(V3_X(v2))) ;
        ya = (int)(nint(V3_Y(v2))) ;
        za = (int)(nint(V3_Z(v2))) ;
        if (xa < 0 || ya < 0 || za < 0 ||
            xa >= mri_aseg->width || ya >= mri_aseg->height || za >= mri_aseg->depth)
          continue ;
        if (xa == Gx && ya == Gy && za == Gz)
          DiagBreak() ;
        aseg_label = nint(MRIgetVoxVal(mri_aseg, xa, ya, za, 0)) ;
        if (seg_label != 0 && !IS_MTL(aseg_label))  // already labeled and not amyg/hippo, skip it
          continue ;
        switch (aseg_label)
        {
        case Left_Cerebellum_White_Matter:
        case Right_Cerebellum_White_Matter:
        case Brain_Stem:
          MRIsetVoxVal(mri_dst, x, y, z, 0, wm_val) ;
          break ;
	case Left_Hippocampus:
	case Right_Hippocampus:
	case Left_Amygdala:
	case Right_Amygdala:
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_Cortex:
        case Left_Pallidum:
        case Right_Pallidum:
        case Left_Thalamus_Proper:
        case Right_Thalamus_Proper:
        case Right_Putamen:
        case Left_Putamen:
        case Right_Caudate:
        case Left_Caudate:
        case Left_Accumbens_area:
        case Right_Accumbens_area:  // remove them from cortex
          MRIsetVoxVal(mri_dst, x, y, z, 0, gm_val) ;
          break ;
        default:
          break ;
        }
      }
    }
  }

  VectorFree(&v1) ; VectorFree(&v2) ; MatrixFree(&m_vox2vox) ;
  return(mri_dst) ;
}
static void
patch_csf_vol(MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort,  MRI *mri_vfrac_csf)
{
  int x, y, z ;
  double v ;

  for (x = 0 ; x < mri_vfrac_wm->width ; x++)
    for (y = 0 ; y < mri_vfrac_wm->height ; y++)
      for (z = 0 ; z < mri_vfrac_wm->depth ; z++)
      {
	v = MRIgetVoxVal(mri_vfrac_wm, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_cortex, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_subcort, x, y, z, 0) ;
	v += MRIgetVoxVal(mri_vfrac_csf, x, y, z, 0) ;
	if (FZERO(v))
	  MRIsetVoxVal(mri_vfrac_csf, x, y, z, 0, 1.0) ;
      }
}
static MRI *
compute_unpartial_volumed_intensities(MRI *mri_src, MRI *mri_vfrac_wm, MRI *mri_vfrac_cortex, MRI *mri_vfrac_subcort, MRI *mri_vfrac_csf, 
				      int whalf0,  double sigma, MRI *mri_dst, int separate_frames) 
{
  int     xi, yi, zi, xk, yk, zk, x, y, z, nvox, whalfx, whalfy, whalfz, num, row ;
  double  w, distsq, dx, dy, dz, norm, total_gm, total_csf, total_wm, vwm, vgm, vcsf ;
  MATRIX   *m_A_pinv, *m_A3, *m_A2, *m_A1, *m_A ;
  VECTOR  *v_I, *v_s3, *v_s2, *v_s1, *v_s ;
  float   wm, gm, csf ;

  whalfx = (int)ceil(whalf0 / mri_src->xsize) ;
  whalfy = (int)ceil(whalf0 / mri_src->ysize) ;
  whalfz = (int)ceil(whalf0 / mri_src->zsize) ;
  nvox = (2*whalfx + 1) * (2*whalfy + 1) * (2*whalfz + 1) ;
  m_A3 = MatrixAlloc(nvox, 3, MATRIX_REAL) ;
  m_A2 = MatrixAlloc(nvox, 2, MATRIX_REAL) ;
  m_A1 = MatrixAlloc(nvox, 1, MATRIX_REAL) ;
  v_s3 = VectorAlloc(3, MATRIX_REAL) ;
  v_s2 = VectorAlloc(2, MATRIX_REAL) ;
  v_s1 = VectorAlloc(1, MATRIX_REAL) ;
  v_I = VectorAlloc(nvox, MATRIX_REAL) ;


  if (mri_dst == NULL)
  {
    if (separate_frames)
    {
      mri_dst = MRIallocSequence(mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT, 3) ;
      MRIcopyHeader(mri_src, mri_dst) ;
    }
    else
      mri_dst = MRIcloneDifferentType(mri_src, MRI_FLOAT) ;
  }

  for (x = 0 ; x < mri_src->width ; x++)
  {
    if (!(x % 32))
      printf("%d of %d\n", x, mri_src->width) ;
    for (y = 0 ; y < mri_src->height ; y++)
      for (z = 0 ; z < mri_src->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;

	total_gm = total_csf = total_wm = 0.0 ;
	for (norm = 0.0, num = 0, xk = -whalfx ; xk <= whalfx ; xk++)
	{
	  dx = xk*mri_src->xsize ;
	  xi = mri_src->xi[x+xk] ;
	  for (yk = -whalfy ; yk <= whalfy ; yk++)
	  {
	    dy = yk*mri_src->ysize ;
	    yi = mri_src->yi[y+yk] ;
	    for (zk = -whalfz ; zk <= whalfz ; zk++)
	    {
	      dz = zk*mri_src->zsize ;
	      zi = mri_src->zi[z+zk] ;
	      distsq = dx*dx + dy*dy + dz*dz ;
	      w = exp(-.5*distsq/(sigma*sigma)) ;
	      norm += w ;
	      VECTOR_ELT(v_I, num+1) = w*MRIgetVoxVal(mri_src, xi, yi, zi, 0) ;
	      vwm = MRIgetVoxVal(mri_vfrac_wm, xi, yi, zi, 0) ;
	      vgm = MRIgetVoxVal(mri_vfrac_cortex, xi, yi, zi, 0) + MRIgetVoxVal(mri_vfrac_subcort, xi, yi, zi, 0) ;
	      vcsf = MRIgetVoxVal(mri_vfrac_csf, xi, yi, zi, 0)  ; ;
	      *MATRIX_RELT(m_A3, num+1, 1) = w*vwm ; *MATRIX_RELT(m_A3, num+1, 2) = w*vgm ; *MATRIX_RELT(m_A3, num+1, 3) = w*vcsf ;
	      total_gm += w*vgm ; total_csf += w*vcsf ; total_wm += w*vwm ;
	      num++ ;
	      if (num  == Gdiag_no)
		DiagBreak() ;
	    }
	  }
	}

	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	if (total_csf < total_gm/100 || total_csf < total_wm/100)
	  total_csf = 0 ;
	if (total_gm < total_csf/100 || total_gm < total_wm/100)
	  total_gm = 0 ;
	if (total_wm < total_gm/100 || total_wm < total_csf/100)
	  total_wm = 0 ;

	if (!FZERO(total_gm))   // some gm in this voxel
	{
	  if (!FZERO(total_wm) && !FZERO(total_csf))
	  {
	    v_s = v_s3 ;
	    m_A = m_A3 ;
	  }

	  else if (!FZERO(total_wm))   // estimate wm and gm, but not csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 2) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else if (!FZERO(total_csf))  // estimate gm and csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 2) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 3) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else   // only gm in this voxel
	  {
	    for (row = 1 ; row < m_A3->rows ; row++)
	      *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 2) ;

	    v_s = v_s1 ;
	    m_A = m_A1 ;
	  }
	}
	else if (!FZERO(total_wm))  // some wm in voxel, but no gm
	{
	  if (!FZERO(total_csf))  // estimate wm and csf
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	    {
	      *MATRIX_RELT(m_A2, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;
	      *MATRIX_RELT(m_A2, row, 2) = *MATRIX_RELT(m_A3, row, 3) ;
	    }
	    v_s = v_s2 ;
	    m_A = m_A2 ;
	  }
	  else   // only wm in this voxel
	  {
	    for (row = 1 ; row <= m_A3->rows ; row++)
	      *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 1) ;

	    v_s = v_s1 ;
	    m_A = m_A1 ;
	  }
	}
	else  // only csf in this region
	{
	  for (row = 1 ; row <= m_A3->rows ; row++)
	    *MATRIX_RELT(m_A1, row, 1) = *MATRIX_RELT(m_A3, row, 3) ;

	    v_s = v_s1 ;
	    m_A = m_A1 ;
	}

	m_A_pinv = MatrixPseudoInverse(m_A, NULL) ;
	if (m_A_pinv == NULL)
	  continue ;

	MatrixMultiply(m_A_pinv, v_I, v_s) ;
	vwm = MRIgetVoxVal(mri_vfrac_wm, x, y, z, 0) ;
	vcsf = MRIgetVoxVal(mri_vfrac_csf, x, y, z, 0) ;
	vgm = MRIgetVoxVal(mri_vfrac_cortex, x, y, z, 0) + MRIgetVoxVal(mri_vfrac_subcort, x, y, z, 0) ;
	wm = gm = csf = 0.0 ;
	if (!FZERO(total_wm))  // wm in 1st col
	{
	  wm = (float)VECTOR_ELT(v_s, 1) ; 
	  if (!FZERO(total_gm))
	  {
	    gm = (float)VECTOR_ELT(v_s, 2) ; 
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 3) ; 
	  }
	  else   // wm but no gm, csf in col 2
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 2) ; 
	}
	else  // no wm in this region
	{
	  if (!FZERO(total_gm))  // some gm in this voxel
	  {
	    gm = (float)VECTOR_ELT(v_s, 1) ; 
	    if (!FZERO(total_csf))
	      csf = (float)VECTOR_ELT(v_s, 2) ; 
	  }
	  else  // only csf in this voxel
	    csf = (float)VECTOR_ELT(v_s, 1) ; 
	}

	if (separate_frames)
	{
	  MRIsetVoxVal(mri_dst, x, y, z, 0, wm) ;
	  MRIsetVoxVal(mri_dst, x, y, z, 1, gm) ;
	  MRIsetVoxVal(mri_dst, x, y, z, 2, csf) ;
	}
	else  // store the value for the class with the largest volume fraction
	{
	  if (vwm > vgm && vwm > vcsf) // wm always in first col
	    MRIsetVoxVal(mri_dst, x, y, z, 0, wm) ;  // mostly wm
	  else if (vgm > vcsf)   // mostly gm
	    MRIsetVoxVal(mri_dst, x, y, z, 0, gm) ;  
	  else   // wm not estimated, in col 1
	    MRIsetVoxVal(mri_dst, x, y, z, 0, csf) ;  

	  if (!devFinite(MRIgetVoxVal(mri_dst, x, y, z, 0)))
	  {
	    DiagBreak() ;
	    MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
	  }
	}
	MatrixFree(&m_A_pinv) ;
      }
  }


  MatrixFree(&m_A3) ; VectorFree(&v_s3) ; VectorFree(&v_I) ;
  MatrixFree(&m_A2) ; VectorFree(&v_s2) ; 
  MatrixFree(&m_A1) ; VectorFree(&v_s1) ; 
  return(mri_dst) ;
}
static MRI *
MRISsimulateAtrophy(MRI *mri_norm, MRI *mri_unpv_intensities, MRI  *mri_wm, MRI *mri_subcort_gm, MRI *mri_cortex, MRI *mri_csf,
		    LABEL *area, double atrophy_frac, MRI *mri_norm_atrophy, MRI **pmri_cortex_out, MRI **pmri_csf_out) 
{
  MRI *mri_filled, *mri_csf_out, *mri_cortex_out, *mri_unpv_intensities_out ;
  int x, y, z, n, xi, yi, zi, xk, yk, zk, xv, yv, zv ;
  float gm_frac, csf_frac, wm_frac, wm_intensity, gm_intensity, csf_intensity, out_val, gm_reduced, max_csf ;
  LABEL  *lvox ;

  mri_unpv_intensities_out = MRIcopy(mri_unpv_intensities, NULL) ;
  mri_csf_out = MRIcopy(mri_csf, NULL) ;
  mri_cortex_out = MRIcopy(mri_cortex, NULL) ;

  mri_norm_atrophy = MRIcopy(mri_norm, mri_norm_atrophy) ;
  mri_filled = MRIclone(mri_norm, NULL) ;

  lvox = LabelToVoxel(area, mri_norm, NULL) ;  // convert label to voxel coords
  for (n = 0 ; n < lvox->n_points ; n++)
  {
    xv = nint(lvox->lv[n].x) ; yv = nint(lvox->lv[n].y) ; zv = nint(lvox->lv[n].z) ;
    for (xk = -1 ; xk <= 1 ; xk++)
      for (yk = -1 ; yk <= 1 ; yk++)
	for (zk = -1 ; zk <= 1 ; zk++)
	{
	  xi = mri_norm->xi[xv+xk] ;  yi = mri_norm->yi[yv+yk] ;  zi = mri_norm->zi[zv+zk] ;
	  if (MRIgetVoxVal(mri_filled, xi, yi, zi, 0)) // already added atrophy here
	    continue ;
	  if (xi == Gx && yi == Gy && zi == Gz)
	    DiagBreak() ;
	  csf_frac = MRIgetVoxVal(mri_csf, xi, yi, zi, 0) ;
	  gm_frac = MRIgetVoxVal(mri_cortex, xi, yi, zi, 0) ;
	  if (csf_frac > 0)
	    max_csf = 1.0 ;
	  else
	    max_csf = MRImaxInNbhd6Connected(mri_csf,  xi,  yi, zi, 0) ;
	  if (max_csf > 0 && gm_frac > 0)  // can only simulate in voxels where there is some of each
	  {
	    gm_reduced = atrophy_frac * gm_frac * max_csf ;
	    if (xi == Gx && yi == Gy && zi == Gz)
	      DiagBreak() ;
	    gm_frac -= gm_reduced ; csf_frac += gm_reduced ;
	    MRIsetVoxVal(mri_csf_out, xi, yi, zi, 0, csf_frac) ;
	    MRIsetVoxVal(mri_cortex_out, xi, yi, zi, 0, gm_frac) ;
	    MRIsetVoxVal(mri_filled, xi, yi, zi, 0, 1) ;
	  }
	}
  }

  LabelFree(&lvox) ; 
  for (x = 0 ; x < mri_norm->width ; x++)
    for (y = 0 ; y < mri_norm->height ; y++)
      for (z = 0 ; z < mri_norm->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	wm_frac = MRIgetVoxVal(mri_wm, x, y, z, 0) ;
	csf_frac = MRIgetVoxVal(mri_csf_out, x, y, z, 0) ;
	gm_frac = MRIgetVoxVal(mri_cortex_out, x, y, z, 0) ;
	gm_frac += MRIgetVoxVal(mri_subcort_gm, x, y, z, 0) ;

#if 0
	if (FEQUAL(csf_frac,1))   // retain contralateral hemi
	  continue ;
#endif
	wm_intensity = MRIgetVoxVal(mri_unpv_intensities_out, x, y, z, 0) ;
	gm_intensity = MRIgetVoxVal(mri_unpv_intensities_out, x, y, z, 1) ;
	csf_intensity = MRIgetVoxVal(mri_unpv_intensities_out, x, y, z, 2) ;
	out_val = wm_frac * wm_intensity + gm_frac * gm_intensity + csf_frac*csf_intensity ;
	MRIsetVoxVal(mri_norm_atrophy, x, y, z, 0, out_val) ;
      }

  MRIfree(&mri_filled) ; 
  // these are copies of those supplied by the caller - not the original volumes (alloced at start)
  MRIfree(&mri_unpv_intensities_out) ; 
  if (pmri_csf_out)
    *pmri_csf_out = mri_csf_out ; 
  else
    MRIfree(&mri_csf_out) ; 
  if (pmri_cortex_out)
    *pmri_cortex_out = mri_cortex_out ;
  else
    MRIfree(&mri_cortex_out) ;
  return(mri_norm_atrophy) ;
}

