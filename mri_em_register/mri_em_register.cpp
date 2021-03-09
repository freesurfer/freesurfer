/**
 * @brief linear registration to a gca atlas
 *
 * pick a bunch of samples and use them to find the transform that
 * maximizes the likelihood of the image given the atlas and the transform
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef HAVE_OPENMP
#include "romp_support.h"
#endif
#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"
#include "mri2.h"
#include "connectcomp.h"

#include "emregisterutils.h"
#include "findtranslation.h"
#include "fsinit.h"

#define DEFAULT_MIN_SCALES 3
#define MM_FROM_EXTERIOR  5  // distance into brain mask to go when erasing super bright CSF voxels

static int MIN_SCALES = DEFAULT_MIN_SCALES ;
static int clamp_set = 0 ;
static double Gclamp = 6 ;   // robust threshold - everything less likely than -Gclamp will be set to -Gclamp
static int remove_cerebellum = 0 ;
static int mark_gcas_classes(GCA_SAMPLE *gcas, int nsamples) ;
static void printUsage(void);

static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;

static int bigvent = 0 ;
static int skull = 0 ;  /* if 1, aligning to image with skull */
static int rigid = 0 ;

GCA_SAMPLE *Ggcas ;

#define MAX_SCALE_PCT 0.15
static float max_scale_pct = MAX_SCALE_PCT ;
static int Gscale_samples = 0 ;
int robust = 0 ;
/*
  allowable distance from an unknown sample to one in brain. Default
  is 1 implying just a ring of unknowns alloowed outside brain. If aligning
  to skull images this gets bigger to allow the full csf/skull context to
  be used in the alignment.
*/
static int unknown_nbr_spacing = 1 ;
static int vent_spacing = -1 ;

int use_variance = 0 ;

static float label_scales[MAX_CMA_LABELS] ;
static float label_offsets[MAX_CMA_LABELS] ;

static double find_optimal_linear_xform
(GCA *gca, GCA_SAMPLE *gcas,
 MRI *mri,
 int nsamples, MATRIX *m_L, MATRIX *m_origin,
 float min_angle, float max_angle,
 float min_scale, float max_scale,
 float min_trans, float max_trans,
 float angle_steps, float scale_steps, float trans_steps,
 int nreductions);

const char         *Progname ;
static MORPH_PARMS  parms ;

static char *T2_mask_fname = NULL ;
static double T2_thresh = 0 ;
static char *aparc_aseg_fname = NULL ;
static char *mask_fname = NULL ;
static char *norm_fname = NULL ;

static char *xform_name = NULL ;
static char *long_reg_fname = NULL ;

static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;

static int map_to_flash = 0 ;
static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static int baby = 0 ;

float G_wm_mean, G_gm_mean, G_fluid_mean ;
static int nomap = 0 ;

static char *sample_fname = NULL ;
static char *transformed_sample_fname = NULL ;
static char *normalized_transformed_sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 0 ;

#define MAX_SPACING  8
static int max_spacing = MAX_SPACING ;
static int nscales = 1 ;

static int use_contrast = 0 ;
static float min_prior = MIN_PRIOR ;
static double tol = 0.001 ;
static double tx = 0.0 ;
static double ty = 0.0 ;
static double tz = 0.0 ;
static double rzrot = 0.0 ;
static double rxrot = 0.0 ;
static double ryrot = 0.0 ;

int exvivo = 0 ;
static int remove_lh = 0 ;
static int remove_rh = 0 ;

static FILE *diag_fp = NULL ;

static LTA *Glta = NULL ;
static TRANSFORM *transform = NULL ;

static int translation_only = 0 ;
static int get_option(int argc, char *argv[]) ;
static int register_mri
(MRI *mri_in, GCA *gca, MP *parms, int passno, int spacing) ;

static char *renormalization_fname = NULL ;
static char *tissue_parms_fname = NULL ;
static int center = 1 ;
static int nreductions = 1 ;
static int noscale = 0 ;
static int noiscale = 0 ;
static double prior_iscale = -1 ;
static int num_xforms = 1 ;
static int transform_loaded = 0 ;
static char *gca_mean_fname = NULL ;

static int ninsertions = 0 ;
static int insert_labels[MAX_INSERTIONS] ;
static int insert_intensities[MAX_INSERTIONS] ;
static int insert_coords[MAX_INSERTIONS][3] ;
static int insert_whalf[MAX_INSERTIONS] ;

static MATRIX *find_optimal_transform
(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas,
 int nsamples, MATRIX *m_L, int passno,
 int write_iterations, int spacing) ;

static double blur_sigma = 0.0f ;

/*
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

#define NPARMS           12
#define NSAMPLES        (NPARMS*20)
#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;
static int nsamples = NSAMPLES ;

static MRI *
apply_transform(MRI *mri, GCA *gca, MATRIX *m_L) 
{
  MRI *mri_aligned ;

  mri->nframes = 1 ;
  Glta->xforms[0].m_L = m_L ;
  mri_aligned = MRIlinearTransform(mri, NULL, m_L) ;
  GCAcopyDCToMRI(gca, mri_aligned) ;
  mri->nframes = gca->ninputs ;

  if (mri->xsize < gca->xsize || mri->ysize < gca->ysize || mri->zsize < gca->zsize)
  {
    MRI *mri_tmp ;
    mri_tmp = MRIextract(mri_aligned, NULL, 0, 0, 0, gca->width, gca->height, gca->depth) ;
    GCAcopyDCToMRI(gca, mri_tmp) ;
    MRIfree(&mri_aligned) ;
    mri_aligned = mri_tmp ;
    mri_aligned->xsize = gca->xsize ; mri_aligned->ysize = gca->ysize ; mri_aligned->zsize = gca->zsize ;
  }
  return(mri_aligned) ;
}

char *rusage_file=NULL;
int n_omp_threads = 1;

int
main(int argc, char *argv[])
{  
  if (0) 
  {
      int i;
      const char* sep = "";
      for (i = 0; i < argc; i++) 
      {
	fputs(sep, stderr);
        fputs(argv[i], stderr);
	sep = " ";
      }
      fputs("\n", stderr);
  }

  char         *gca_fname, *in_fname, *out_fname, fname[STRLEN], **av ;
  MRI          *mri_in, *mri_tmp, *mri_dst ;
  GCA          *gca /*, *gca_tmp, *gca_reduced*/ ;
  int          ac, nargs, i, ninputs, scale, spacing;
  int          exclude_list[MAX_CMA_LABEL+1] ;
  int          msec, minutes, seconds, min_left_cbm, min_right_cbm ;
  Timer start ;
  float        old_log_p, log_p ;

  FSinit() ;

  nargs = handleVersionOption(argc, argv, "mri_em_register");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  memset(exclude_list, 0, sizeof(exclude_list))  ;
  /*    exclude_list[0] = 1 ;*/
  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.max_levels = 0 ;
  parms.dt = 5e-6 ;  /* was 5e-6 */
  parms.tol = 1e-5 ;
  parms.momentum = 0.8 ;
  parms.niterations = 25 ;
  Progname = argv[0] ;

  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
  {
    printUsage();
    exit(1);
  }

#ifdef HAVE_OPENMP
  n_omp_threads = omp_get_max_threads(); 
  printf("\n== Number of threads available to %s for OpenMP = %d == \n", Progname, n_omp_threads);
#endif

  ninputs = argc-3 ;
  printf("reading %d input volumes...\n", ninputs) ;
  gca_fname = argv[ninputs+1] ;
  out_fname = argv[ninputs+2] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  printf("logging results to %s.log\n", parms.base_name) ;
  parms.rigid = rigid ;

  if (skull)
  {
    int l ;

    for (l = 0 ; l < MAX_CMA_LABELS ; l++)
      if (IS_BRAIN(l))
      {
        exclude_list[l] = 1 ;
      }
  }
  start.reset() ;
  ///////////  read GCA //////////////////////////////////////////////////
  printf("reading '%s'...\n", gca_fname) ;
  fflush(stdout) ;

  {
      Timer start ;
      start.reset() ;
      gca = GCAread(gca_fname) ;
      int msec = start.milliseconds() ;
      int seconds = nint((float)msec/1000.0f) ;
      printf("GCAread took %d minutes and %d seconds.\n",
         seconds / 60, seconds % 60) ;    
  }
  
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;
  if (gca->ninputs > 1 && clamp_set == 0)
  {
    printf("disabling robust clamp for multispectral data\n") ;
    Gclamp = 2000 ;
  }

  if (baby)
  {
    GCArenormalizeClass(gca, CSF_CLASS, 240.0/190.0) ;  // for T2 space
  }
  if (exvivo)
  {
    GCArenormalizeClass(gca, GM_CLASS, 180.0/140.0) ; // for exvivo contrast
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      printf("writing test.gca\n") ;
      GCAwrite(gca, "test.gca") ;
    }
  }
  if (remove_lh)
  {
    Gvx = nint(gca->width*.4) ; // only one hemi - assume it is left/right centered in FOV
    GCAremoveHemi(gca, 1) ;  // for exvivo contrast
  }
  if (remove_rh)
  {
    GCAremoveHemi(gca, 0) ;  // for exvivo contrast
  }
  if (remove_cerebellum)
  {
    GCAremoveLabel(gca, Brain_Stem) ;
    GCAremoveLabel(gca, Left_Cerebellum_Cortex) ;
    GCAremoveLabel(gca, Left_Cerebellum_White_Matter) ;
    GCAremoveLabel(gca, Right_Cerebellum_White_Matter) ;
    GCAremoveLabel(gca, Right_Cerebellum_Cortex) ;
  }

  /////////  -novar option //////////////////////////////////////////////
  if (novar)
  {
    GCAregularizeCovariance(gca,1.0);
  }
  //    GCAunifyVariance(gca) ;

  GCAfixSingularCovarianceMatrices(gca) ;
  GCAapplyRenormalization(gca, label_scales, label_offsets, 0) ;

  ////////// -renorm fname ////////////////////////////////////////////
  if (renormalization_fname)
  {
    FILE   *fp ;
    int    *labels, nlines, i ;
    float  *intensities, f1, f2 ;
    char   *cp, line[STRLEN] ;

    fp = fopen(renormalization_fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                Progname, renormalization_fname) ;

    cp = fgetl(line, 199, fp) ;
    nlines = 0 ;
    while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ;
      intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
      {
        DiagBreak() ;
      }
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ;
    free(intensities) ;
  }

  //////////////////////////////////////////////////////////////
  // create a list of MRI volumes
  for (i = 0 ; i < ninputs ; i++)
  {
    in_fname = argv[i+1] ;
    printf("reading '%s'...\n", in_fname) ;
    fflush(stdout) ;
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
                Progname, in_fname) ;

    TRs[i] = mri_tmp->tr ;
    fas[i] = mri_tmp->flip_angle ;
    TEs[i] = mri_tmp->te ;
    if (mask_fname)
    {
      MRI *mri_mask ;
      int val ;

      mri_mask = MRIread(mask_fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, mask_fname) ;
      for (val = 0 ; val < MIN_WM_VAL ; val++)
	MRImask(mri_tmp, mri_mask, mri_tmp, val, 0) ;

      MRIfree(&mri_mask) ;
      if (parms.write_iterations != 0 && Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s_masked", parms.base_name) ;
        printf("writing masked volume to %s...\n", fname) ;
        MRIwriteImageViews(mri_tmp, fname, IMAGE_SIZE) ;
        sprintf(fname, "%s_masked.mgz", parms.base_name) ;
        fflush(stdout);
        MRIwrite(mri_tmp, fname) ;
      }
    }
    if (T2_mask_fname)
    {
      MRI *mri_T2, *mri_aparc_aseg = nullptr;

      mri_T2 = MRIread(T2_mask_fname) ;
      if (!mri_T2)
	ErrorExit(ERROR_NOFILE, "%s: could not open T2 mask volume %s.\n",
		  Progname, mask_fname) ;
      if (aparc_aseg_fname)   // use T2 and aparc+aseg to remove non-brain stuff
      {
 	mri_aparc_aseg = MRIread(aparc_aseg_fname) ;
	if (mri_aparc_aseg == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not open aparc+aseg volume %s.\n",
		    Progname, aparc_aseg_fname) ;
      }

      MRImask_with_T2_and_aparc_aseg(mri_tmp, mri_tmp, mri_T2, mri_aparc_aseg, T2_thresh, MM_FROM_EXTERIOR) ;
      MRIfree(&mri_T2) ; MRIfree(&mri_aparc_aseg) ;
    }
    if (i == 0)
    {
      mri_in = MRIallocSequence(mri_tmp->width,
                                mri_tmp->height,
                                mri_tmp->depth,
                                mri_tmp->type,
                                ninputs) ;
      MRIcopyHeader(mri_tmp, mri_in) ;
    }
    MRIcopyFrame(mri_tmp, mri_in, 0, i) ;
    MRIfree(&mri_tmp) ;
  }
  //////////////////////////////////////////////////////////////
  if (alpha > 0)
  {
    mri_in->flip_angle = alpha ;
  }
  if (TR > 0)
  {
    mri_in->tr = TR ;
  }
  if (TE > 0)
  {
    mri_in->te = TE ;
  }

  /////////////////////  -flash ////////////////////////////////
  if (map_to_flash || gca->type==GCA_PARAM)
  {
    GCA *gca_tmp ;

    printf("mapping GCA into %d-dimensional FLASH space...\n",
           mri_in->nframes) ;
    // that means gca->ninputs = nframes
    gca_tmp = GCAcreateFlashGCAfromParameterGCA
              (gca, TRs, fas, TEs,
               mri_in->nframes, GCA_DEFAULT_NOISE_PARAMETER) ;
    // now the type is set gca->type = GCA_FLASH
    GCAfree(&gca) ;
    gca = gca_tmp ;
    if (gca->ninputs > 1)
      /* multispectral - normalize to remove bias field */
    {
      GCAnormalizeMeans(gca, 100) ;
      MRInormalizeSequence(mri_in, 100) ;
    }
    else
    {
      GCAhistoScaleImageIntensities(gca, mri_in, skull==0) ;
    }
    if (ninputs != gca->ninputs)
      ErrorExit
      (ERROR_BADPARM,
       "%s: must specify %d inputs, not %d for this atlas\n",
       Progname, gca->ninputs, ninputs) ;
    /*                GCAhistoScaleImageIntensities(gca, mri_in, skull==0) ;*/
    if (novar)
    {
      GCAregularizeCovariance(gca,1.0);
    }
    //      GCAunifyVariance(gca) ;
  }
  ///////////////////////////////////////////////////////////////
  else if (gca->type == GCA_FLASH)
  {
    GCA *gca_tmp ;
    int need_map_flag = 0;
    int n;

    if (gca->ninputs != mri_in->nframes)
    {
      need_map_flag = 1;
    }
    else
    {
      for (n = 0 ; n < mri_in->nframes; n++)
      {
        if (!FZERO(gca->TRs[n] - TRs[n]))
        {
          need_map_flag = 1;
        }
        if (!FZERO(gca->FAs[n] - fas[n]))
        {
          need_map_flag = 1;
        }
        if (!FZERO(gca->TEs[n] - TEs[n]))
        {
          need_map_flag = 1;
        }
      }
    }

    if (nomap == 0 && need_map_flag == 1)
    {
      printf("GCAcreateFLASHGCAfromFlashGCA...\n");
      gca_tmp = GCAcreateFlashGCAfromFlashGCA
                (gca, TRs, fas, TEs, mri_in->nframes) ;
      GCAfree(&gca) ;
      gca = gca_tmp ;
    }
#if 0
    if (gca->ninputs > 1)
      /* multispectral - normalize to remove bias field */
    {
      GCAnormalizeMeans(gca, 100) ;
      MRInormalizeSequence(mri_in, 100) ;
    }
    else
#endif
      GCAhistoScaleImageIntensities(gca, mri_in, skull == 0) ;
    if (novar)
    {
      GCAregularizeCovariance(gca,1.0);
    }
    //      GCAunifyVariance(gca) ;
  }
  if (ninputs != gca->ninputs)
    ErrorExit
    (ERROR_BADPARM,
     "%s: must specify %d input volumes, not %d for this atlas\n",
     Progname, gca->ninputs, ninputs) ;

  parms.vgca = (void *)gca ;
  if (exvivo == 0)
  {
    printf("freeing gibbs priors...") ;
    fflush(stdout);
    GCAfreeGibbs(gca) ;
    printf("done.\n") ;
    fflush(stdout);
  }

  //////////////////////////// -example option ////////////////////////
  if (example_T1)
  {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not read example segmentation from %s",
       Progname, example_segmentation) ;
    mri_T1 = MRIread(example_T1) ;
    if (!mri_T1)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not read example T1 from %s",
       Progname, example_T1) ;
    printf("scaling atlas intensities using specified examples...\n") ;
    fflush(stdout);
    MRIeraseBorderPlanes(mri_seg, 1) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ;
    MRIfree(&mri_T1) ;
  }

  if (tissue_parms_fname)   /* use FLASH forward model */
  {
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;
  }

  ////////////////////////// -fsample option //////////////////////////
  if (sample_fname)
  {
    GCAwriteSamples(gca, mri_in, parms.gcas, nsamples, sample_fname) ;
    printf("samples written to %s\n", sample_fname) ;
  }
#if 0
  if (gca_reduced != gca)
  {
    GCAfree(&gca_reduced) ;
  }
#endif
  /////////////////////////  -d tx ty tz ////////////////////////////////
  if (!FZERO(tx) || !FZERO(ty) || !FZERO(tz))
  {
    MRI *mri_tmp ;

    printf("translating second volume by (%2.1f, %2.1f, %2.1f)\n",
           tx, ty, tz) ;
    mri_tmp = MRItranslate(mri_in, NULL, tx, ty, tz) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  ////////////////////////  -r rxrot ryrot rzrot ///////////////////////
  // note the rotation is only around each axis
  // not around the center of the volume
  // Therefore, the volume may get rotated out.....
  // Rotation is yRot*xRot*zRot composition applying on the right
  if (!FZERO(rzrot))
  {
    MRI *mri_tmp ;

    printf(
      "rotating second volume by %2.1f degrees around Z axis\n",
      (float)DEGREES(rzrot)) ;
    mri_tmp = MRIrotateZ_I(mri_in, NULL, rzrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(rxrot))
  {
    MRI *mri_tmp ;

    printf(
      "rotating second volume by %2.1f degrees around X axis\n",
      (float)DEGREES(rxrot)) ;
    mri_tmp = MRIrotateX_I(mri_in, NULL, rxrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (!FZERO(ryrot))
  {
    MRI *mri_tmp ;

    printf(
      "rotating second volume by %2.1f degrees around Y axis\n",
      (float)DEGREES(ryrot)) ;
    mri_tmp = MRIrotateY_I(mri_in, NULL, ryrot) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  //////////////////// if -t is not used
  if (!transform_loaded)   /* wasn't preloaded */
  {
    MATRIX *m_voxsize, *m_tmp ;
    // allocate only one transform
    // mri_in is used only to set x0, y0, z0
    // Note that vox-to-vox transform
    parms.transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, mri_in) ;
    Glta = parms.lta = (LTA *)transform->xform ;
    ////////////////////////////////////////////////////
    // now start working (remember this is vox-to-vox transform)
    parms.lta->xforms[0].m_L = MatrixIdentity(4, NULL) ;
    printf("accounting for voxel sizes in initial transform\n") ;
    m_voxsize = MatrixIdentity(4, NULL) ;
    *MATRIX_RELT(m_voxsize, 1,1) = mri_in->xsize ;
    *MATRIX_RELT(m_voxsize, 2,2) = mri_in->ysize ;
    *MATRIX_RELT(m_voxsize, 3,3) = mri_in->zsize ;
    m_tmp = MatrixMultiply(m_voxsize, parms.lta->xforms[0].m_L,NULL) ;
    MatrixCopy(m_tmp, parms.lta->xforms[0].m_L) ;
    MatrixFree(&m_voxsize) ; MatrixFree(&m_tmp) ;
  }

  //////////  -b option ///////////////////////////////////////////////////
  if (!FZERO(blur_sigma))
  {
    MRI *mri_tmp, *mri_kernel ;

    mri_kernel = MRIgaussian1d(blur_sigma, 100) ;
    mri_tmp = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_in) ;
    mri_in = mri_tmp ;
  }
  if (ninsertions > 0)
    GCAinsertLabels(gca, mri_in, transform, ninsertions,
                    insert_labels, insert_intensities, insert_coords,
                    insert_whalf) ;

  if (parms.write_iterations != 0)
  {
    char fname[STRLEN] ;
    MRI  *mri_gca ;
    mri_gca = MRIclone(mri_in, NULL) ;
    GCAbuildMostLikelyVolume(gca, mri_gca) ;
    sprintf(fname, "%s_target", parms.base_name) ;
    MRIwriteImageViews(mri_gca, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s_target.mgz", parms.base_name) ;
    printf("writing target volume to %s...\n", fname) ;
    fflush(stdout);
    MRIwrite(mri_gca, fname) ;
    MRIfree(&mri_gca) ;
  }
  i = 0 ;

  ////////////////////////////////////////////////////////////////
  // change scale up to nscales (spacing is halved)
  // default nscales = 1
  for (spacing = max_spacing, scale = 0 ;
       scale < nscales ;
       scale++, spacing /= 2)
  {
    if (0 && skull)
      parms.gcas = GCAfindExteriorSamples(gca, &nsamples,spacing,
                                          min_prior,unknown_nbr_spacing, 0) ;
    else if (use_contrast) // -contrast option
    {
      parms.gcas = GCAfindContrastSamples(gca,&nsamples, spacing,min_prior);
    }
    else
      parms.gcas = GCAfindStableSamples
                   (gca, &nsamples,spacing,
                    min_prior,exclude_list,unknown_nbr_spacing, vent_spacing) ;
    Ggcas = parms.gcas ;  // for diags
    mark_gcas_classes(parms.gcas, nsamples) ;
    printf("************************************************\n");
    printf("spacing=%d, using %d sample points, tol=%2.2e...\n",
           spacing, nsamples, parms.tol) ;
    printf("************************************************\n");
    fflush(stdout);
    parms.nsamples = nsamples ;
    if (sample_fname)
    {
      GCAwriteSamples(gca, mri_in, parms.gcas, nsamples, sample_fname) ;
      printf("samples written\n") ;
    }
    old_log_p = local_GCAcomputeLogSampleProbability
      (gca, parms.gcas, mri_in, ((LTA *)(transform->xform))->xforms[0].m_L, nsamples, exvivo, Gclamp) ;
    // real work done here
    register_mri(mri_in, gca, &parms, i, spacing) ;
    // calculate log_p
    log_p = local_GCAcomputeLogSampleProbability
      (gca, parms.gcas, mri_in, ((LTA *)(transform->xform))->xforms[0].m_L, nsamples, exvivo, Gclamp) ;
    printf("pass %d, spacing %d: log(p) = %2.3f (old=%2.3f)\n",
           i+1, spacing, log_p, old_log_p) ;
    GCAfreeSamples(&parms.gcas, nsamples) ;
    parms.tol *= 10 ;
    i++ ;
  }

  // change nsamples to all samples
  if (0 && skull)
    parms.gcas = GCAfindExteriorSamples(gca, &nsamples,spacing/2,
                                        min_prior,unknown_nbr_spacing, 0) ;
  else
    parms.gcas = GCAfindAllSamples(gca, &nsamples,
                                   exclude_list, unknown_nbr_spacing);//HACK, bigvent) ;
  mark_gcas_classes(parms.gcas, nsamples) ;
  parms.nsamples = nsamples ;
  parms.tol = 1e-7 ;
  parms.start_t++ ;

  //////////////////// diagnostics //////////////////////////////////////////
  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0))
  {
    MRI *mri_aligned ;

    mri_aligned = apply_transform(mri_in, gca, parms.lta->xforms[0].m_L) ;
    sprintf(fname, "%s%03d", parms.base_name, parms.start_t) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s%03d.mgz", parms.base_name, parms.start_t) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;

    /*                Glta->xforms[0].m_L = m_L ;*/
    sprintf(fname, "%s%3.3d_fsamples.mgz",
            parms.base_name, parms.start_t) ;
    GCAtransformAndWriteSamples
    (gca, mri_in, parms.gcas, nsamples, fname, transform) ;
    sprintf(fname, "%s%3.3d_pvals.mgz",
            parms.base_name, parms.start_t) ;
    GCAtransformAndWriteSamplePvals
    (gca, mri_in, parms.gcas, nsamples, fname, transform) ;
    sprintf(fname, "%s%3.3d_means.mgz",
            parms.base_name, parms.start_t) ;
    GCAtransformAndWriteSampleMeans
    (gca, mri_in, parms.gcas, nsamples, fname, transform) ;
  }
  /////////////////////////////////////////////////////////////////////////
  if (exvivo)
  {
    GCAupdateDistributions(gca, mri_in, transform);
    GCAwrite(gca, "exvivo.gca") ;

  }

  parms.start_t++ ;
  printf("transform before final EM align:\n") ;
  MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
  printf("\n") ;
  printf("**************************************************\n");
  printf(" EM alignment process ...\n");
  printf(" Computing final MAP estimate using %d samples. \n", nsamples) ;
  printf("**************************************************\n");
  fflush(stdout);
  parms.mri_in = mri_in ;  /* for diagnostics */
  parms.clamp = Gclamp ;
  MRIemAlign(mri_in, gca, &parms, parms.lta->xforms[0].m_L) ;


  printf("final transform:\n") ;
  MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
  printf("\n") ;

  /////////////////////////diagnostics/////////////////////////////////
  if ((Gdiag & DIAG_WRITE) && (parms.write_iterations != 0))
  {
    MRI *mri_aligned ;

    mri_aligned = apply_transform(mri_in, gca, parms.lta->xforms[0].m_L) ;
    sprintf(fname, "%s%03d", parms.base_name, parms.start_t) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    sprintf(fname, "%s%03d.mgz", parms.base_name, parms.start_t) ;
    MRIwrite(mri_aligned, fname) ;
    MRIfree(&mri_aligned) ;
  }
  /////////////////////////////////////////////////////////////////////
  printf("writing output transformation to %s...\n", out_fname) ;
  fflush(stdout);
  // writing transform section here
  // create gca volume for outputting dirction cosines and c_(ras)
  mri_dst = MRIallocHeader(gca->width,
                           gca->height,
                           gca->depth,
                           mri_in->type,
                           mri_in->nframes);
  GCAcopyDCToMRI(gca, mri_dst);
  strcpy(mri_dst->fname,gca_fname); // copy gca name
  if (!stricmp(out_fname+strlen(out_fname)-3, "XFM"))
  {
    printf("converting xform to RAS...\n") ;
    printf("initial:\n") ;
    fflush(stdout);
    MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
    MRIvoxelXformToRasXform
    (mri_in, mri_dst, parms.lta->xforms[0].m_L, parms.lta->xforms[0].m_L) ;
    printf("final:\n") ;
    MatrixPrint(stdout, parms.lta->xforms[0].m_L) ;
    parms.lta->type = LINEAR_RAS_TO_RAS ;
    fflush(stdout);
  }
  else
  {
    parms.lta->type = LINEAR_VOX_TO_VOX ;
  }

  // add src and dst info
  getVolGeom(mri_in, &parms.lta->xforms[0].src);
  getVolGeom(mri_dst, &parms.lta->xforms[0].dst);

  LTAwriteEx(parms.lta, out_fname) ;

  ///////////////////////////////////////////// end of writing transform
  if (parms.lta->type == LINEAR_RAS_TO_RAS)  /* convert back to voxel */
  {
    printf("converting xform back to voxel...\n") ;
    fflush(stdout);
    MRIrasXformToVoxelXform
    (mri_in, mri_dst, parms.lta->xforms[0].m_L, parms.lta->xforms[0].m_L) ;
    parms.lta->type = LINEAR_VOX_TO_VOX ;
  }
  MRIfree(&mri_dst);
  //////////////////////////////////////////////////////////////////////
  if (transformed_sample_fname)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.pvals.mgz",parms.base_name) ;
    printf("writing transformed samples to %s and pvals to %s...\n",
           transformed_sample_fname, fname) ;
    fflush(stdout);
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples,
                                transformed_sample_fname, transform) ;
    printf("samples written\n") ;
    fflush(stdout);
    GCAtransformAndWriteSamplePvals
    (gca, mri_in, parms.gcas, nsamples, fname, transform) ;
  }
  ///////////////////////////////////////////////////////////////////////
  if (norm_fname)
  {
    int   *ordered_indices, i, label, nused, nleft_cbm, nright_cbm;
    int   nleft_used, nright_used ;
    MRI   *mri_norm ;

    local_GCAcomputeLogSampleProbability(gca, parms.gcas, mri_in,
                                         ((LTA *)(transform->xform))->xforms[0].m_L, nsamples, exvivo, Gclamp) ;
#if 0
    GCAnormalizedLogSampleProbability(gca, parms.gcas, mri_in,
                                      transform, nsamples, exvivo) ;
#endif
    /* make "unknowns" the bottom of the list */
    for (nleft_cbm = nright_cbm = nused = i = 0 ; i < nsamples ; i++)
    {
      label = parms.gcas[i].label ;
      if (((label != Left_Cerebral_White_Matter ) &&
           (label != Right_Cerebral_White_Matter ) &&
           (label != Left_Cerebellum_White_Matter ) &&
           (label != Right_Cerebellum_White_Matter )))
      {
        parms.gcas[i].log_p = -100000 ;
        parms.gcas[i].label = -1 ;
      }
      else
      {
        nused++ ;
        if (label == Left_Cerebellum_White_Matter )
        {
          nleft_cbm++ ;
        }
        else if (label == Right_Cerebellum_White_Matter)
        {
          nright_cbm++ ;
        }
      }
    }
    GCAremoveOutlyingSamples(gca, parms.gcas, mri_in,
                             transform, nsamples, 2.0) ;

    /* rank samples by log probability */
    ordered_indices = (int *)calloc(nsamples, sizeof(int)) ;
    GCArankSamples(gca, parms.gcas, nsamples, ordered_indices) ;

#if 0
    if (DIAG_VERBOSE_ON)
    {
      for (nleft_used = nright_used = i = 0 ; i < nsamples ; i++)
      {
        if (parms.gcas[ordered_indices[i]].label ==
            Left_Cerebellum_White_Matter)
        {
          nleft_used++ ;
        }
        else if (parms.gcas[ordered_indices[i]].label ==
                 Right_Cerebellum_White_Matter)
        {
          nright_used++ ;
        }
        if (parms.gcas[ordered_indices[i]].label ==
            Right_Cerebral_White_Matter)
        {
          DiagBreak() ;
        }
        if (parms.gcas[ordered_indices[i]].x == Gx &&
            parms.gcas[ordered_indices[i]].y == Gy &&
            parms.gcas[ordered_indices[i]].z == Gz)
        {
          DiagBreak() ;
        }
        if ((parms.gcas[ordered_indices[i]].label ==
             Right_Cerebral_White_Matter) &&
            (parms.gcas[ordered_indices[i]].prior > 0.95) &&
            fabs(MRIvox(mri_in, parms.gcas[ordered_indices[i]].x,
                        parms.gcas[ordered_indices[i]].y,
                        parms.gcas[ordered_indices[i]].z)-\
                 parms.gcas[ordered_indices[i]].mean) < 5)
        {
          DiagBreak() ;
        }
      }
    }
#endif

    /* remove the least likely samples */
    printf("sorting %d (%d/%d l/r cerebellum) white matter points by "
           "likelihood\n", nused, nleft_cbm, nright_cbm) ;
    for (nleft_used = nright_used = i = 0 ;
         i < nint(nused*ctl_point_pct);
         i++)
    {
      if (parms.gcas[ordered_indices[i]].label ==
          Left_Cerebellum_White_Matter)
      {
        nleft_used++ ;
      }
      else if (parms.gcas[ordered_indices[i]].label ==
               Right_Cerebellum_White_Matter)
      {
        nright_used++ ;
      }
#if 0
      if (diag_fp)
        fprintf(diag_fp, "%d %d %d %2.1f %d\n",
                parms.gcas[ordered_indices[i]].x,
                parms.gcas[ordered_indices[i]].y,
                parms.gcas[ordered_indices[i]].z,
                parms.gcas[ordered_indices[i]].mean,
                parms.gcas[ordered_indices[i]].label) ;
#endif
      if (parms.gcas[ordered_indices[i]].label ==
          Right_Cerebral_White_Matter)
      {
        DiagBreak() ;
      }
      if (parms.gcas[ordered_indices[i]].x == Gx &&
          parms.gcas[ordered_indices[i]].y == Gy &&
          parms.gcas[ordered_indices[i]].z == Gz)
      {
        DiagBreak() ;
      }

    }

    min_left_cbm = nint(nleft_cbm*ctl_point_pct+.9) ;
    min_right_cbm = nint(nright_cbm*ctl_point_pct+.9) ;
    printf("%d/%d (l/r) cerebellar points initially in "
           "top %d%%, min (%d,%d)\n"
           , nleft_used,nright_used, (int)(ctl_point_pct*100.0f),
           min_left_cbm, min_right_cbm) ;

    for (i = nint(nused*ctl_point_pct) ; i < nsamples ; i++)
    {
      if ((parms.gcas[ordered_indices[i]].label ==
           Left_Cerebellum_White_Matter) && nleft_used < min_left_cbm)
      {
        nleft_used++ ;
      }
      else  if ((parms.gcas[ordered_indices[i]].label ==
                 Right_Cerebellum_White_Matter) &&
                nright_used< min_right_cbm)
      {
        nright_used++ ;
      }
      else
      {
        parms.gcas[ordered_indices[i]].label = -1 ;
      }
#if 0
      if (diag_fp && parms.gcas[ordered_indices[i]].label > 0)
        fprintf(diag_fp, "%d %d %d %d %d %2.1f %d\n",
                i, ordered_indices[i],
                parms.gcas[ordered_indices[i]].x,
                parms.gcas[ordered_indices[i]].y,
                parms.gcas[ordered_indices[i]].z,
                parms.gcas[ordered_indices[i]].mean,
                parms.gcas[ordered_indices[i]].label) ;
#endif

      if (parms.gcas[ordered_indices[i]].label > 0 &&
          MRIvox(mri_in, parms.gcas[ordered_indices[i]].x,
                 parms.gcas[ordered_indices[i]].y,
                 parms.gcas[ordered_indices[i]].z) < 50)
      {
        DiagBreak() ;
      }
    }

#if 0
    /* replace sample label with pct rank so that we can write it out easily
       for diagnostic purposes (a bit of a hack) */
    for (i = 0 ; i < nsamples ; i++)
    {
      float pct ;

      pct = 100.0f*(float)(nsamples-i)/nsamples;
      parms.gcas[ordered_indices[i]].label = pct ;
    }
    GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples,
                                norm_fname, transform) ;
#endif

    if (nint(nused*ctl_point_pct) == 0)
      ErrorPrintf(ERROR_BADPARM,
                  "%s: too few control points (%d) for normalization",
                  Progname, nused) ;
    else
    {
      if (normalized_transformed_sample_fname)
        GCAtransformAndWriteSamples(gca, mri_in, parms.gcas, nsamples,
                                    normalized_transformed_sample_fname,
                                    transform) ;

      mri_norm = GCAnormalizeSamples(mri_in, gca, parms.gcas, nsamples,
                                     transform, ctl_point_fname) ;
      printf("writing normalized volume to %s...\n", norm_fname) ;
      fflush(stdout);
      if (MRIwrite(mri_norm, norm_fname)  != NO_ERROR)
        ErrorExit
        (ERROR_BADFILE,
         "%s: could not write normalized volume to %s",
         Progname, norm_fname);

      MRIfree(&mri_norm) ;
    }
  }

#if 0
  if (gca)
  {
    GCAfree(&gca) ;
  }
#endif
  if (mri_in)
  {
    MRIfree(&mri_in) ;
  }

  // Print usage stats to the terminal (and a file is specified)
  //PrintRUsage(RUSAGE_SELF, "mri_em_register ", stdout);
  //if(rusage_file) WriteRUsage(RUSAGE_SELF, "", rusage_file);
  printf("#VMPC# mri_em_register VmPeak  %d\n",GetVmPeak());

  ///////////////////////////////////////////////////////////////
  msec = start.milliseconds() ;
  printf("FSRUNTIME@ mri_em_register %7.4f hours %d threads\n", msec/(1000.0*60.0*60.0), n_omp_threads);
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("registration took %d minutes and %d seconds.\n",minutes, seconds) ;
  if (diag_fp)
  {
    fclose(diag_fp) ;
  }
  exit(0) ;
  return(0) ;
}


static int
register_mri
(MRI *mri_in, GCA *gca, MORPH_PARMS *parms, int passno, int spacing)
{
  MATRIX  *m_L ;
  // get the stored transform (vox-to-vox transform)
  m_L = MatrixCopy(parms->lta->xforms[0].m_L, NULL) ;

  ////////////////////////////////////////////////////////////////////////
  fprintf(stdout, "register_mri: find_optimal_transform\n");
  find_optimal_transform(mri_in, gca, parms->gcas, parms->nsamples,m_L,passno,
                         parms->write_iterations, spacing);

  /* make sure transform and lta are the same (sorry - retrofitting!) */
  if (!parms->lta)
  {
    parms->transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    Glta = parms->lta = (LTA *)transform->xform ;
  }
  else
  {
#if 0
    parms->transform = transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
    transform->xform = (void *)parms->lta ;
#endif
  }

  MatrixCopy(m_L, parms->lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_SHOW)
  {
    printf("global search transform:\n") ;
    MatrixPrint(stdout, m_L) ;
  }

  parms->start_t++ ;
  printf("***********************************************\n");
  printf("Computing MAP estimate using %d samples...\n", parms->nsamples) ;
  printf("***********************************************\n");
  fflush(stdout);
  parms->mri_in = mri_in ;  /* for diagnostics */

  //////////////// calling MRIemAlign() //////////////////////////
  if (exvivo == 0)
  {
    MRIemAlign(mri_in, gca, parms, m_L) ;
  }
  MatrixCopy(m_L, parms->lta->xforms[0].m_L) ;

  printf("Resulting transform:\n") ;
  MatrixPrint(stdout, parms->lta->xforms[0].m_L) ;
  printf("\n") ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRI *mri_aligned ;
    char fname[STRLEN] ;

    mri_aligned = apply_transform(mri_in, gca, parms->lta->xforms[0].m_L) ;
    sprintf(fname, "%s_after_alignment", parms->base_name) ;
    MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
    MRIfree(&mri_aligned) ;
  }

  return(NO_ERROR) ;
}

#define DEFAULT_MAX_STEPS 5
static double MAX_ANGLES = DEFAULT_MAX_STEPS ;
static double MAX_ANGLE  = RADIANS(30);
//static double MIN_ANGLE  = RADIANS(2) ;

static int max_angles = DEFAULT_MAX_STEPS ;
static int max_scales = DEFAULT_MAX_STEPS ;
static int MAX_TRANS_STEPS = DEFAULT_MAX_STEPS ;
static double MAX_TRANS = 30 ;

static MATRIX *
find_optimal_transform
(MRI *mri, GCA *gca, GCA_SAMPLE *gcas, int nsamples,
 MATRIX *m_L, int passno, int write_iterations, int spacing)
{
  MATRIX   *m_origin ;
  MRI      *mri_gca  ;
  double   gca_means[3], /*in_means[3], dx, dy, dz,*/ max_log_p, old_max,
           max_angle, angle_steps, min_scale, max_scale, scale_steps, scale,
           delta, mean ;
  int      niter, good_step, done, nscales, scale_samples ;
  float      min_search_scale ;
#if 0
  int        min_real_bin, mri_peak ;
  float      min_real_val, fmax, fmin ;
  MRI_REGION box, gca_box ;
  HISTOGRAM *h_mri, *h_smooth ;
#endif

  fprintf(stdout,
          "find_optimal_transform: nsamples %d, passno %d, spacing %d\n",
          nsamples, passno, spacing);

#define MIN_SEARCH_SCALE 0.1
  min_search_scale = MIN_SEARCH_SCALE ;

  // used by various spacing value

  if (Gscale_samples)
  {
    scale_samples = Gscale_samples ;
  }
  else
  {
    if (spacing >= 16)
    {
      scale_samples = 5 ;
    }
    else
    {
      scale_samples = 3 ;
    }
  }
  if (spacing >= 8)
  {
    min_search_scale /= 4;
  }

  /////////////////////////////////////////////////////////////////////////////
  max_log_p = local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L,nsamples, exvivo, Gclamp) ;

  // create volume from gca with the size of input
  mri_gca = MRIclone(mri, NULL) ;
  GCAmri(gca, mri_gca) ;// set the values mri_gca has the same DC as gca has
  //
  MRIcenterOfMass(mri_gca, gca_means, 0) ;
  // unit matrix
  m_origin = MatrixIdentity(4, NULL) ;
  // set the translation  (center = 1)
  *MATRIX_RELT(m_origin, 1, 4) = gca_means[0]*(float)center ;
  *MATRIX_RELT(m_origin, 2, 4) = gca_means[1]*(float)center ;
  *MATRIX_RELT(m_origin, 3, 4) = gca_means[2]*(float)center ;
  *MATRIX_RELT(m_origin, 4, 4) = 1 ;

  if (passno == 0)
  {
    if (noiscale)     // default is noiscale = 0 and thus perform
    {
      if (prior_iscale > 0)
	MRIscalarMulFrame(mri, mri, prior_iscale, 0);
    }
    else   // perform histogram scaling
    {
      if (Gdiag & DIAG_WRITE && write_iterations > 0)
      {
	char fname[STRLEN] ;
	MRI  *mri_aligned ;
	
	Glta->xforms[0].m_L = m_L ;
	mri_aligned = apply_transform(mri, gca, m_L) ;
	sprintf(fname, "%s_before_intensity.mgz", parms.base_name) ;
	printf("writing snapshot to %s...\n", fname) ;
	fflush(stdout);
	MRIwrite(mri_aligned, fname) ;
	MRIfree(&mri_aligned) ;
      }
      GCAhistoScaleImageIntensities(gca, mri, skull == 0) ;
    }

    //////////////// diagnostics ////////////////////////////////
    if (Gdiag & DIAG_WRITE && write_iterations > 0)
    {
      char fname[STRLEN] ;
      MRI  *mri_aligned ;
      
      mri_aligned = apply_transform(mri, gca, m_L) ;
      sprintf(fname, "%s000.mgz", parms.base_name) ;

      printf("writing snapshot to %s...\n", fname) ;
      fflush(stdout);
      MRIwrite(mri_aligned, fname) ;
      sprintf(fname, "%s000", parms.base_name) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;

      sprintf(fname, "%s000_fsamples.mgz", parms.base_name) ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
				  fname, transform) ;
      sprintf(fname, "%s000_pvals.mgz", parms.base_name) ;
      GCAtransformAndWriteSamplePvals(gca, mri, gcas, nsamples,
				      fname, transform) ;
      sprintf(fname, "%s000_means.mgz", parms.base_name) ;
      GCAtransformAndWriteSampleMeans
	(gca, mri, parms.gcas, nsamples, fname, transform) ;
      MRIfree(&mri_aligned) ;
    }
  }

  if (passno == 0 && !transform_loaded)   /* only first time*/
  {
    /////////////////////////////////////////////////////////////
    /* first align centroids */
    if (gca_mean_fname)
    {
      printf("writing gca volume to %s...\n", gca_mean_fname) ;
      MRIwrite(mri_gca, gca_mean_fname) ;
      printf("done\n") ;
      fflush(stdout);
    }

    printf("initial log_p = %2.3f\n", max_log_p) ;
    ///////////////////////////////////////////////////////////////////////
#if 0
#if 1
    MRIvalRange(mri, &fmin, &fmax) ;
    h_mri = MRIhistogram(mri, nint(fmax-fmin+1)) ;
    h_mri->counts[0] = 0 ; /* ignore background */
    h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
    mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins/3) ;
    min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25) ;
    min_real_val = h_smooth->bins[min_real_bin] ;
    printf("using real data threshold=%2.1f\n", min_real_val) ;

    MRIfindApproximateSkullBoundingBox(mri, min_real_val, &box) ;
    HISTOfree(&h_mri) ;
    HISTOfree(&h_smooth) ;
    printf("input bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
           box.x, box.y, box.z, box.x+box.dx, box.y+box.dy, box.z+box.dz) ;

    MRIvalRange(mri_gca, &fmin, &fmax) ;
    h_mri = MRIhistogram(mri_gca, nint(fmax-fmin+1)) ;
    h_mri->counts[0] = 0 ; /* ignore background */
    h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
    mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins/3) ;
    min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25) ;
    min_real_val = h_smooth->bins[min_real_bin] ;
    printf("using real data threshold=%2.1f\n", min_real_val) ;
    MRIfindApproximateSkullBoundingBox(mri_gca, min_real_val, &gca_box) ;
    HISTOfree(&h_mri) ;
    HISTOfree(&h_smooth) ;
    printf("gca bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
           gca_box.x, gca_box.y, gca_box.z,
           gca_box.x+gca_box.dx, gca_box.y+gca_box.dy,gca_box.z+gca_box.dz) ;
    in_means[0] = box.x + 0.5*box.dx ;
    in_means[1] = box.y + 0.5*box.dy ;
    in_means[2] = box.z + 0.5*box.dz ;
    gca_means[0] = gca_box.x + 0.5*gca_box.dx ;
    gca_means[1] = gca_box.y + 0.5*gca_box.dy ;
    gca_means[2] = gca_box.z + 0.5*gca_box.dz ;
#else
    MRIcenterOfMass(mri, in_means, 0) ;
    printf("input centroid (%2.1f, %2.1f, %2.1f), "
           "gca centroid (%2.1f, %2.1f, %2.1f)\n",
           in_means[0], in_means[1], in_means[2],
           gca_means[0], gca_means[1], gca_means[2]) ;
#endif

#if 0
    in_means[1] = box.y+box.dx*0.55 ;
    printf("resetting superior/inferior centroid to %2.1f\n", in_means[1]) ;
#endif

    /* now apply translation to take in centroid to ref centroid */
    dx = gca_means[0] - in_means[0] ;
    dy = gca_means[1] - in_means[1] ;
    dz = gca_means[2] - in_means[2] ;
    *MATRIX_RELT(m_L, 1, 4) = dx ;
    *MATRIX_RELT(m_L, 2, 4) = dy ;
    *MATRIX_RELT(m_L, 3, 4) = dz ;
    max_log_p = local_GCAcomputeLogSampleProbability
      (gca, gcas, mri, m_L,nsamples, exvivo, Gclamp) ;
    printf("initial translation: (%2.1f, %2.1f, %2.1f): log p = %2.3f\n",
           dx,dy,dz, max_log_p) ;
#else ///////////////this is executed  ////////////////////////////////////
    fprintf(stdout, "************************************************\n");
    fprintf(stdout, "First Search limited to translation only.\n");
    fprintf(stdout, "************************************************\n");
    if (skull)  // remove excess neck
    {
      MRI_REGION box ;
      int        min_real_bin, mri_peak ;
      float      min_real_val, fmax, fmin ;
      HISTOGRAM  *h_mri, *h_smooth ;

      MRIvalRange(mri, &fmin, &fmax) ;
      h_mri = MRIhistogram(mri, nint(fmax-fmin+1)) ;
      h_mri->counts[0] = 0 ; /* ignore background */
      h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
      mri_peak = HISTOfindHighestPeakInRegion(h_smooth, 0, h_smooth->nbins/3) ;
      min_real_bin = HISTOfindEndOfPeak(h_smooth, mri_peak, .25) ;
      min_real_val = 2*h_smooth->bins[min_real_bin] ;

      if (mriConformed(mri) && min_real_val > fmax/2)
      {
        min_real_val = fmax/2 ;
      }
      printf("using real data threshold=%2.1f\n", min_real_val) ;
      MRIfindApproximateSkullBoundingBox(mri, min_real_val, &box) ;
      box.y = (box.y+1.1*box.dx) ;  // use saggital to estimate s/i extent
      box.dy = (mri->height-1)-box.y ;
      box.x = box.z = 0 ;
      box.dx = mri->width ;
      box.dz = mri->depth;
      if (mriConformed(mri))
      {
        if (box.y < 190)
        {
          box.y = 190 ;
          box.dy = (mri->height-1)-box.y ;
        }
      }
      MRIfillBox(mri, &box, 0) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRIwrite(mri, "i.mgz") ;
      }
      HISTOfree(&h_mri) ;
      HISTOfree(&h_smooth) ;
    }
    max_log_p = find_optimal_translation(gca, gcas, mri, nsamples, m_L,
                                         -200, 200, 19, 7, Gclamp) ;
    max_log_p = local_GCAcomputeLogSampleProbability
      (gca, gcas, mri, m_L,nsamples, exvivo, Gclamp) ;
    fprintf(stdout,
            "Found translation: (%2.1f, %2.1f, %2.1f): log p = %4.3f\n",
            *MATRIX_RELT(m_L, 1, 4),
            *MATRIX_RELT(m_L, 2, 4),
            *MATRIX_RELT(m_L, 3, 4),
            max_log_p) ;

    /*
    printf( "%s: Debugging exit\n", __FUNCTION__ );
    exit( EXIT_SUCCESS );
    */

#endif
    /////////////////////////////////////////////////////////////////////////

    parms.start_t++ ;   // translation counts as one timestep
    if (write_iterations != 0)
    {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = apply_transform(mri, gca, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, parms.start_t) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgz", parms.base_name, parms.start_t) ;
      printf("writing image after centering to %s...\n", fname) ;
      fflush(stdout);
      MRIwrite(mri_aligned, fname) ;
      Glta->xforms[0].m_L = m_L ;
      sprintf(fname, "%s%03d_fsamples.mgz", parms.base_name, parms.start_t) ;
      GCAtransformAndWriteSamples(gca, mri, gcas, nsamples,
                                  fname, transform) ;
      sprintf(fname, "%s%03d_pvals.mgz", parms.base_name, parms.start_t) ;
      GCAtransformAndWriteSamplePvals(gca, mri, gcas, nsamples,
                                      fname, transform) ;
      sprintf(fname, "%s%03d_means.mgz", parms.base_name, parms.start_t) ;
      GCAtransformAndWriteSampleMeans(gca, mri, gcas, nsamples,
                                      fname, transform) ;
      MRIfree(&mri_aligned) ;
      mri_aligned = MRIinverseLinearTransform(mri, NULL, m_L) ;
      sprintf(fname, "%s%03d.inv.mgz", parms.base_name, parms.start_t) ;
      MRIwrite(mri_aligned, fname) ;
      MRIfree(&mri_aligned) ;
    }
    MRIfree(&mri_gca) ;
  }

  max_angle = MAX_ANGLE ;
  angle_steps = max_angles ;
  max_scale = 1+max_scale_pct ;
  min_scale = 1-max_scale_pct ;
  scale_steps = max_scales ;

  //////////////// loop here ////////////////////////////////////////////
  niter = nscales = 0 ; 
  scale = 1.0 ;
  good_step = 0 ;
  done = 0 ;
  do
  {
    Timer start ;
    start.reset() ;

    old_max = max_log_p ;
    printf("****************************************\n");
    printf("Nine parameter search.  iteration %d nscales = %d ...\n",
           niter, nscales);
    printf("****************************************\n");
    fflush(stdout);
    max_log_p = find_optimal_linear_xform
                (gca, gcas, mri, nsamples,
                 m_L, m_origin,
                 -(max_angle*scale),
                 (max_angle*scale),
                 min_scale,
                 max_scale,
                 -scale*(spacing/16.0)*MAX_TRANS,
                 scale*(spacing/16.0)*MAX_TRANS,
                 max_angles, scale_samples, /*MAX_TRANS_STEPS*/3, 2);
    fflush(stdout);


    if (write_iterations != 0)
    {
      char fname[STRLEN] ;
      MRI *mri_aligned ;

      mri_aligned = apply_transform(mri, gca, m_L) ;
      sprintf(fname, "%s%03d", parms.base_name, parms.start_t+niter+1) ;
      MRIwriteImageViews(mri_aligned, fname, IMAGE_SIZE) ;
      sprintf(fname, "%s%03d.mgz",
              parms.base_name, parms.start_t+niter+1) ;
      printf("writing %s\n", fname) ;
      MRIwrite(mri_aligned, fname) ;
#if 0
      MRIwrite(mri_aligned, fname) ;
#else
      Glta->xforms[0].m_L = m_L ;
      sprintf(fname, "%s%3.3d_fsamples.mgz",
              parms.base_name, parms.start_t+niter+1) ;
      GCAtransformAndWriteSamples
      (gca, mri, gcas, nsamples, fname, transform) ;
      sprintf(fname, "%s%3.3d_means.mgz",
              parms.base_name, parms.start_t+niter+1) ;
      GCAtransformAndWriteSampleMeans
      (gca, mri, parms.gcas, nsamples, fname, transform) ;
#endif
      MRIfree(&mri_aligned) ;
    }
    fprintf(stdout,
            "Result so far: scale %2.3f: max_log_p=%2.3f, "
            "old_max_log_p =%2.3f (thresh=%2.1f)\n",
            scale,max_log_p, old_max, old_max+fabs(tol*old_max)) ;
    MatrixPrint(stdout, m_L);
    fflush(stdout);
    /* search a finer nbhd (if do-while continues) */
    if ((max_log_p < old_max+fabs(tol*old_max))) /* couldn't take a step */
    {
      scale *= 0.25 ;
      if (scale < min_search_scale)
      {
        printf("min search scale %f reached\n", min_search_scale);
        break ;
      }
      mean = (max_scale + min_scale)/2 ;
      delta = (max_scale - min_scale)/2 ;
      max_scale = 1.0 + delta*scale ;
      min_scale = 1.0 - delta*scale ;
      done = (good_step == 0) ;
      good_step = 0 ;
      printf("reducing scale to %2.4f\n", scale) ;
      nscales++ ;
    }
    else
    {
      good_step = 1 ;  /* took at least one good step at this scale */
    }

    int msec = start.milliseconds();
    int seconds = nint((float)msec/1000.0f) ;
    printf("iteration took %d minutes and %d seconds.\n", seconds / 60, seconds % 60) ;
	    
    niter++ ;
  }
  while (nscales < MIN_SCALES || (done == FALSE)) ;

  parms.start_t += niter ;
  MatrixFree(&m_origin) ;
  return(m_L) ;
}




/*----------------------------------------------------------------------
  ----------------------------------------------------------------------*/
#include "mri_em_register.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_em_register_help_xml,
                mri_em_register_help_xml_len);
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE"))
  {
    // seems like not used.
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_dist = %2.2f\n", parms.l_dist) ;
  }
  else if (!strcmp(option, "-HELP")||!strcmp(option, "-USAGE"))
  {
    printUsage();
    exit(1) ;
  }
  else if (!strcmp(option, "NOMAP"))
  {
    // seems not used
    nomap = 1 ;
  }
  else if (!strcmp(option, "BIGVENT"))
  {
    bigvent = 1 ;
    vent_spacing = 30 ;
    printf("handling expanded ventricles with vent_spacing = %d\n", vent_spacing) ;
  }
  else if (!strcmp(option, "VAR"))
  {
    use_variance = 1 ;
    printf("minimizing within-label intensity variance\n") ;
  }
  else if (!strcmp(option, "LSCALE"))
  {
    int l ;
    l = atoi(argv[2]) ;
    label_scales[l] = atof(argv[3]) ;
    nargs = 2 ;
    printf("scaling label %s by %2.2f\n", cma_label_to_name(l), label_scales[l]) ;
    for (l = 0 ; l < MAX_CMA_LABELS ; l++)
      if (FZERO(label_scales[l]))
	label_scales[l] = 1.0 ;
  }
  else if (!strcmp(option, "CLAMP"))
  {
    clamp_set = 1 ;
    Gclamp = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting robust clamp to %2.3f\n", Gclamp) ;
  }
  else if (!strcmp(option, "ROBUST"))
  {
    robust = 1 ;
  }
  else if (!stricmp(option, "FLASH"))
  {
    map_to_flash = 1 ;
    printf("using FLASH forward model to predict intensity values...\n") ;
  }
  else if (!stricmp(option, "MAX_ANGLE"))
  {
    MAX_ANGLE = RADIANS(atof(argv[2])) ;
    printf("using %2.2f deg as max angle for rotational search\n",
           DEGREES(MAX_ANGLE)) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "MAX_SCALE"))
  {
    max_scale_pct = atof(argv[2]) ;
    printf("using %2.2f%% as max scaling delta ([%2.3f %2.3f])\n",
           100*max_scale_pct, 1-max_scale_pct, 1+max_scale_pct) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "LH"))
  {
    remove_rh = 1  ;
    printf("removing right hemisphere labels\n") ;
  }
  else if (!stricmp(option, "RH"))
  {
    remove_lh = 1  ;
    printf("removing left hemisphere labels\n") ;
  }
  else if (!stricmp(option, "EXVIVO"))
  {
    exvivo = 1 ;
    printf("assuming input is ex vivo image with unknown conditional distributions\n") ;
  }
  else if (!stricmp(option, "BABY"))
  {
    baby = 1 ;
    printf("using baby brain intensity model\n") ;
  }
  else if (!stricmp(option, "INSERT"))
  {
    if (ninsertions >= MAX_INSERTIONS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many insertions (%d) specified\n",
                Progname, ninsertions) ;

    insert_labels[ninsertions] = atoi(argv[2]) ;
    insert_intensities[ninsertions] = atoi(argv[3]) ;
    insert_coords[ninsertions][0] = atoi(argv[4]) ;
    insert_coords[ninsertions][1] = atoi(argv[5]) ;
    insert_coords[ninsertions][2] = atoi(argv[6]) ;
    insert_whalf[ninsertions] = atoi(argv[7]) ;
    printf("inserting label %d (%s) at (%d, %d, %d) with intensity = %d\n",
           insert_labels[ninsertions],
           cma_label_to_name(insert_labels[ninsertions]),
           insert_coords[ninsertions][0],
           insert_coords[ninsertions][1],
           insert_coords[ninsertions][2],
           insert_intensities[ninsertions]) ;

    ninsertions++ ;
    nargs = 6 ;
  }
  else if (!stricmp(option, "T2MASK"))
  {
    T2_mask_fname = argv[2] ;
    T2_thresh = atof(argv[3]) ;
    nargs = 2 ;
    printf("using T2 volume %s thresholded at %f to mask input volume...\n", 
	   T2_mask_fname, T2_thresh) ;
  }
  else if (!stricmp(option, "AMASK"))
  {
    aparc_aseg_fname = argv[2] ;
    T2_mask_fname = argv[3] ;
    T2_thresh = atof(argv[4]) ;
    nargs = 3 ;
    printf("using aparc+aseg vol %s and T2 volume %s thresholded at %f to mask input volume...\n", 
	   aparc_aseg_fname, T2_mask_fname, T2_thresh) ;
  }
  else if (!stricmp(option, "THREADS"))
  {
#ifdef HAVE_OPENMP
    sscanf(argv[2], "%d", &n_omp_threads);
    omp_set_num_threads(n_omp_threads);
    printf("threads %d\n", n_omp_threads);
#else
    printf("multithreading is not available\n");
#endif
    nargs = 1 ;
  }
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "SKULL"))
  {
    unknown_nbr_spacing = 5 ;
    printf("aligning to atlas containing skull, "
           "setting unknown_nbr_spacing = %d\n",
           unknown_nbr_spacing) ;
  }
  else if (!stricmp(option, "rusage"))
  {
    // resource usage
    rusage_file = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "NOCEREBELLUM"))
  {
    remove_cerebellum = 1 ;
    printf("removing cerebellum from atlas\n") ;
  }
  else if (!strcmp(option, "RIGID"))
  {
    rigid = 1 ;
    printf("constraining transform to be rigid\n") ;
  }
  else if (!strcmp(option, "UNS"))
  {
    unknown_nbr_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting unknown_nbr_spacing = %d\n", unknown_nbr_spacing) ;
  }
  else if (!strcmp(option, "VENT_SPACING"))
  {
    vent_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("disallowing WM voxels within %d mm of ventricles in atlas (useful for big ventricles)\n", vent_spacing) ;
  }
  /////// debug options //////////////////////////////////
  else if (!strcmp(option, "DIAG"))
  {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not open diag file %s for writing",
       Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!strcmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %s (%d)\n",
           cma_label_to_name(Ggca_label), Ggca_label) ;
  }
  ////////// TR, TE, Alpha ////////////////////////////////
  else if (!strcmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!strcmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!strcmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!strcmp(option, "EXAMPLE"))
  {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  /////////////// writing out various samples /////////////////
  else if (!strcmp(option, "SAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "FSAMPLES") || !strcmp(option, "ISAMPLES"))
  {
    transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing transformed control points to %s...\n",
           transformed_sample_fname) ;
  }
  else if (!strcmp(option, "NSAMPLES"))
  {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n",
           normalized_transformed_sample_fname) ;
  }
  ///////////////////
  else if (!strcmp(option, "CONTRAST"))
  {
    use_contrast = 1 ;
    printf("using contrast to find labels...\n") ;
  }
  else if (!stricmp(option, "RENORM") || !stricmp(option, "RENORMALIZE"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!strcmp(option, "FLASH_PARMS"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!strcmp(option, "TRANSONLY"))
  {
    translation_only = 1 ;
    printf("only computing translation parameters...\n") ;
  }
  else if (!strcmp(option, "WRITE_MEAN"))
  {
    gca_mean_fname = argv[2] ;
    nargs = 1 ;
    printf("writing gca means to %s...\n", gca_mean_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!strcmp(option, "SPACING"))
  {
    max_spacing = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using max GCA spacing %d...\n", max_spacing) ;
  }
  else if (!stricmp(option, "SCALES") || !stricmp(option, "SCALES"))
  {
    nscales = MIN_SCALES = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding optimal linear transform over %d scales...\n", MIN_SCALES);
  }
  else if (!stricmp(option, "NSCALES"))
  {
    Gscale_samples = atoi(argv[2]) ;
    nargs = 1 ;
    printf("sampling scaling %d times at each scale\n", Gscale_samples) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    printf("dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "TOL"))
  {
    tol = parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "CENTER"))
  {
    center = 1 ;
    printf("using GCA centroid as origin of transform\n") ;
  }
  else if (!strcmp(option, "NOSCALE"))
  {
    noscale = 1 ;
    printf("disabling scaling...\n") ;
  }
  else if (!strcmp(option, "NOISCALE"))
  {
    noiscale = 1 ;
    printf("disabling intensity scaling...\n") ;
  }
  else if (!stricmp(option, "ISCALE"))
  {
    noiscale = 1 ;
    prior_iscale = atof(argv[2]) ;
    nargs = 1 ;
    printf("scaling image intensities by %2.3f\n", prior_iscale) ;
  }
  else if (!strcmp(option, "NUM"))
  {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    printf("finding a total of %d linear transforms\n", num_xforms) ;
  }
  else if (!strcmp(option, "AREA"))
  {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_area = %2.2f\n", parms.l_area) ;
  }
  else if (!strcmp(option, "NLAREA"))
  {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_nlarea = %2.2f\n", parms.l_nlarea) ;
  }
  else if (!strcmp(option, "LEVELS"))
  {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("levels = %d\n", parms.levels) ;
  }
  else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    printf("l_intensity = %2.2f\n", parms.l_intensity) ;
  }
  else if (!stricmp(option, "reduce"))
  {
    nreductions = atoi(argv[2]) ;
    nargs = 1 ;
    printf("reducing input images %d times before aligning...\n",
           nreductions) ;
  }
  else if (!stricmp(option, "nsamples"))
  {
    nsamples = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using %d samples of GCA...\n", nsamples) ;
  }
  else if (!stricmp(option, "norm"))
  {
    norm_fname = argv[2] ;
    nargs = 1 ;
    printf("intensity normalizing and writing to %s...\n",norm_fname);
  }
  else if (!stricmp(option, "trans"))
  {
    MAX_TRANS = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max translation search range to be %2.1f\n", MAX_TRANS) ;
  }
  else if (!stricmp(option, "steps"))
  {
    max_angles = atoi(argv[2]) ;
    nargs = 1 ;
    printf("taking %d angular steps...\n", max_angles) ;
  }
  else switch (*option)
    {
    case 'L':   /* for longitudinal analysis */
    {
      TRANSFORM *reg_transform ;

      xform_name = argv[2] ;
      long_reg_fname = argv[3] ;
      nargs = 2 ;
      printf("reading previously computed atlas xform %s "
             "and applying registration %s\n",
             xform_name, long_reg_fname) ;
      parms.transform = transform = TransformRead(argv[2]) ;
      if (transform == NULL)
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read transform from %s",
         Progname, argv[2]) ;
      Glta = parms.lta = (LTA *)transform->xform ;
      reg_transform = TransformRead(argv[3]) ;
      if (reg_transform == NULL)
        ErrorExit
        (ERROR_NOFILE,
         "%s: could not read registration from %s",
         Progname, argv[3]) ;
      transform_loaded = 1 ;
      TransformInvert(reg_transform, NULL) ;
      MatrixMultiply(((LTA *)(transform->xform))->xforms[0].m_L,
                     ((LTA *)(reg_transform->xform))->inv_xforms[0].m_L,
                     ((LTA *)(transform->xform))->xforms[0].m_L) ;
      TransformFree(&reg_transform) ;
    }
    break ;
    case 'F':
      ctl_point_fname = argv[2] ;
      nargs = 1 ;
      printf("reading manually defined control points from %s\n",
             ctl_point_fname) ;
      break ;
    case 'D':
      tx = atof(argv[2]) ;
      ty = atof(argv[3]) ;
      tz = atof(argv[4]) ;
      nargs = 3 ;
      break ;
    case 'R':
      rxrot = RADIANS(atof(argv[2])) ;
      ryrot = RADIANS(atof(argv[3])) ;
      rzrot = RADIANS(atof(argv[4])) ;
      nargs = 3 ;
      break ;
    case 'T':
      parms.transform = transform = TransformRead(argv[2]) ;
      if (transform == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read xform from %s",
                  Progname, argv[2]) ;
      Glta = parms.lta = (LTA *)transform->xform ;
#if 0
      parms.lta = LTAreadEx(argv[2]) ; // used to be LTAread()
#endif
      if (!parms.lta)
        ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                  Progname, argv[2]) ;
      if (parms.lta->type!=LINEAR_VOX_TO_VOX)
        ErrorExit(ERROR_BADFILE, "%s: must be LINEAR_VOX_TO_VOX (=0), but %d",
                  Progname, argv[2], parms.lta->type) ;
      nargs = 1 ;
      printf("using previously computed transform %s\n", argv[2]) ;
      if (parms.lta->type != LINEAR_VOX_TO_VOX)
      {
        fprintf(stdout,
                "ERROR: must use LINEAR_VOX_TO_VOX (=0) transform. "
                "The type was %d.\n",
                parms.lta->type);
        exit(1);
      }
      transform_loaded = 1 ;
      break ;
    case 'B':
      blur_sigma = atof(argv[2]) ;
      nargs = 1 ;
      printf("blurring input image with sigma=%2.3f\n", blur_sigma);
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'S':
#if 0
      parms.sigma = atof(argv[2]) ;
      printf("using sigma=%2.3f as upper bound on blurring.\n",
             parms.sigma) ;
      nargs = 1 ;
#else
      Gscale_samples = max_scales = MAX_ANGLES = MAX_TRANS_STEPS = max_angles = (float)atoi(argv[2]) ;
      nargs = 1 ;
      printf("examining %2.0f different trans/rot/scale values...\n",
             MAX_ANGLES);
#endif
      break ;
    case '?':
    case 'H':
    case 'U':
      printUsage();
      exit(1) ;
      break ;
    case 'N':
      parms.niterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("niterations = %d\n", parms.niterations) ;
      break ;
    case 'W':
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("write iterations = %d\n", parms.write_iterations) ;
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'P':
      ctl_point_pct = atof(argv[2]) ;
      nargs = 1 ;
      printf("using top %2.1f%% wm points as control points....\n",
             100.0*ctl_point_pct) ;
      break ;
    case 'M':
      parms.momentum = atof(argv[2]) ;
      nargs = 1 ;
      printf("momentum = %2.2f\n", parms.momentum) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }
  fflush(stdout);

  return(nargs) ;
}





/*/////////////////////////////////////////////////////////////
  search 9-dimensional parameter space
*/
static double
find_optimal_linear_xform
(GCA *gca, GCA_SAMPLE *gcas,
 MRI *mri,
 int nsamples, MATRIX *m_L, MATRIX *m_origin,
 float min_angle, float max_angle,
 float min_scale, float max_scale,
 float min_trans, float max_trans,
 float angle_steps, float scale_steps, float trans_steps,
 int nreductions )
{

  MATRIX   *m_rot, *m_x_rot, *m_y_rot, *m_z_rot, *m_tmp,*m_L_tmp,*m_origin_inv,
           *m_tmp2, *m_scale, *m_trans, *m_tmp3 = NULL ;

  double x_max_rot;
  double y_max_rot, z_max_rot, delta_rot;
  double x_max_scale, y_max_scale, z_max_scale;
  double delta_scale, delta_trans;
  double max_log_p, mean_angle;
  double mean_scale, x_max_trans, y_max_trans, z_max_trans, mean_trans ;
  double x_trans, y_trans, z_trans;
  double x_scale, y_scale, z_scale;
  double x_angle, y_angle, z_angle;
  double log_p;
  int i;

  if (rigid)
  {
    min_scale = max_scale = 1.0 ;
  }

  m_trans = MatrixIdentity(4, NULL) ;
  m_origin_inv = MatrixCopy(m_origin, NULL) ;
  *MATRIX_RELT(m_origin_inv, 1, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 2, 4) *= -1 ;
  *MATRIX_RELT(m_origin_inv, 3, 4) *= -1 ;
  m_L_tmp = m_x_rot = m_y_rot = m_z_rot = m_rot = m_tmp = m_tmp2 = NULL ;
  x_max_trans = y_max_trans = z_max_trans \
                              = x_max_rot = y_max_rot = z_max_rot = 0.0 ;
  x_max_scale = y_max_scale = z_max_scale = 1.0f ;
  m_scale = MatrixIdentity(4, NULL) ;
  max_log_p = local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L, nsamples, exvivo, Gclamp) ;

  // Loop a set number of times to polish transform

  for (i = 0 ; i < nreductions ; i++)
  {
    delta_trans = (max_trans-min_trans) / (trans_steps-1) ;
    delta_scale = (max_scale-min_scale) / (scale_steps-1) ;
    if (FZERO(delta_scale) || rigid)
    {
      delta_scale = max_scale ;
    }

    if (angle_steps == 1)
    {
      min_angle = max_angle = 0.0 ;
      delta_rot = 1 ;
    }
    else
    {
      delta_rot = (max_angle-min_angle) / (angle_steps-1) ;
    }
    if (Gdiag & DIAG_SHOW)
    {
      printf("  scanning %2.2f degree nbhd (%2.1f)\n"
             "  scale %2.3f->%2.3f (step %2.3f), "
             "trans %2.2f->%2.2f (step %2.2f)\n",
             (float)DEGREES(max_angle), (float)DEGREES(delta_rot),
             min_scale,max_scale, delta_scale,
             min_trans, max_trans, delta_trans);
      fflush(stdout) ;
    }

    // scale /////////////////////////////////////////////////////////////
    for (x_scale = min_scale ; x_scale <= max_scale ; x_scale += delta_scale)
    {
      /*      printf("x_scale = %2.3f\n", x_scale) ;*/
      *MATRIX_RELT(m_scale, 1, 1) = x_scale ;
      for (y_scale = min_scale ;
           y_scale <= max_scale ;
           y_scale += delta_scale)
      {
        *MATRIX_RELT(m_scale, 2, 2) = y_scale ;
        for (z_scale= min_scale ;
             z_scale <= max_scale;
             z_scale += delta_scale)
        {
          *MATRIX_RELT(m_scale, 3, 3) = z_scale ;

          /* reset translation values */
          *MATRIX_RELT(m_scale, 1, 4) =
            *MATRIX_RELT(m_scale, 2, 4) =
              *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
          m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m_scale) ;

          // angle //////////////////////////////
          for (x_angle = min_angle ;
               x_angle <= max_angle ;
               x_angle += delta_rot)
          {
            m_x_rot = MatrixReallocRotation
                      (4, x_angle, X_ROTATION, m_x_rot) ;
            for (y_angle = min_angle ;
                 y_angle <= max_angle ;
                 y_angle += delta_rot)
            {
              m_y_rot = MatrixReallocRotation
                        (4, y_angle, Y_ROTATION, m_y_rot);
              m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
              for (z_angle= min_angle;
                   z_angle <= max_angle;
                   z_angle += delta_rot)
              {
                m_z_rot = MatrixReallocRotation
                          (4, z_angle,Z_ROTATION,m_z_rot);
                m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
                m_tmp2 = MatrixMultiply
                         (m_rot, m_origin_inv, m_tmp2) ;
                MatrixMultiply(m_origin, m_tmp2, m_rot) ;

                m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
                m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

                // translation //////////
                for (x_trans = min_trans ;
                     x_trans <= max_trans ;
                     x_trans += delta_trans)
                {
                  *MATRIX_RELT(m_trans, 1, 4) = x_trans ;
                  for (y_trans = min_trans ;
                       y_trans <= max_trans ;
                       y_trans += delta_trans)
                  {
                    *MATRIX_RELT(m_trans, 2, 4) = y_trans ;
                    for (z_trans= min_trans ;
                         z_trans <= max_trans ;
                         z_trans += delta_trans)
                    {
                      *MATRIX_RELT(m_trans, 3, 4) =
                        z_trans ;

                      m_L_tmp = MatrixMultiply
                                (m_trans, m_tmp3, m_L_tmp) ;

                      log_p = local_GCAcomputeLogSampleProbability(gca, gcas, mri, m_L_tmp, nsamples, exvivo, Gclamp);
                      if (log_p > max_log_p)
                      {
                        if (exvivo)
                          printf("current estimates G=%d, W=%d, F=%d\n",
                                 (int)G_gm_mean, (int)G_wm_mean, (int)G_fluid_mean) ;
                        max_log_p = log_p ;
                        x_max_scale = x_scale ;
                        y_max_scale = y_scale ;
                        z_max_scale = z_scale ;
                        x_max_rot = x_angle ;
                        y_max_rot = y_angle ;
                        z_max_rot = z_angle ;
                        x_max_trans = x_trans ;
                        y_max_trans = y_trans ;
                        z_max_trans = z_trans ;
                      }
#if 0
                      printf( "%s: log_p = %f\n", __FUNCTION__, log_p );
                      printf( "%s: Translation (%4.2f, %4.2f, %4.2f)\n",
                              __FUNCTION__, x_trans, y_trans, z_trans );
                      printf( "%s: Rotation (%4.2f, %4.2f, %4.2f)\n",
                              __FUNCTION__, x_angle, y_angle, z_angle );
                      printf( "%s: Scale (%4.2f, %4.2f, %4.2f)\n",
                              __FUNCTION__, x_scale, y_scale, z_scale );
                      exit( 0 );
#endif
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (Gdiag & DIAG_SHOW)
    {
      printf("  max log p = %2.3f @ R=(%2.3f,%2.3f,%2.3f),"
             "S=(%2.3f,%2.3f,%2.3f), T=(%2.1f,%2.1f,%2.1f)\n",
             max_log_p, DEGREES(x_max_rot), DEGREES(y_max_rot),
             DEGREES(z_max_rot),x_max_scale, y_max_scale, z_max_scale,
             x_max_trans, y_max_trans,z_max_trans) ;
    }


    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_scale, 1, 4) =
      *MATRIX_RELT(m_scale, 2, 4) = *MATRIX_RELT(m_scale, 3, 4) = 0.0f ;
    *MATRIX_RELT(m_scale,1,1) = x_max_scale ;
    *MATRIX_RELT(m_scale,2,2) = y_max_scale ;
    *MATRIX_RELT(m_scale,3,3) = z_max_scale ;
    m_tmp = MatrixMultiply(m_scale, m_origin_inv, m_tmp) ;
    m_scale = MatrixMultiply(m_origin, m_tmp, m_scale) ;

    x_max_scale = y_max_scale = z_max_scale = 1.0 ;

    mean_scale = (max_scale + min_scale) / 2 ;
    delta_scale = (max_scale-min_scale)/4 ;
    min_scale = mean_scale - delta_scale ;
    max_scale = mean_scale + delta_scale ;

    /* update L to reflect new maximum and search around it */
    m_x_rot = MatrixReallocRotation(4, x_max_rot, X_ROTATION, m_x_rot) ;
    m_y_rot = MatrixReallocRotation(4, y_max_rot, Y_ROTATION, m_y_rot) ;
    m_z_rot = MatrixReallocRotation(4, z_max_rot, Z_ROTATION, m_z_rot) ;
    m_tmp = MatrixMultiply(m_y_rot, m_x_rot, m_tmp) ;
    m_rot = MatrixMultiply(m_z_rot, m_tmp, m_rot) ;
    m_tmp2 = MatrixMultiply(m_rot, m_origin_inv, m_tmp2) ;
    m_rot = MatrixMultiply(m_origin, m_tmp2, m_rot) ;


    m_tmp2 = MatrixMultiply(m_scale, m_rot, m_tmp2) ;
    m_tmp3 = MatrixMultiply(m_tmp2, m_L, m_tmp3) ;

    /* update L to reflect new maximum and search around it */
    *MATRIX_RELT(m_trans, 1, 4) = x_max_trans ;
    *MATRIX_RELT(m_trans, 2, 4) = y_max_trans ;
    *MATRIX_RELT(m_trans, 3, 4) = z_max_trans ;
    m_L_tmp = MatrixMultiply(m_trans, m_tmp3, m_L_tmp) ;



    MatrixCopy(m_L_tmp, m_L) ;

    x_max_trans = y_max_trans = z_max_trans = 0.0 ;
    /* we've translated transform by old maxs */
    mean_trans = (max_trans + min_trans) / 2 ;
    delta_trans = (max_trans-min_trans)/4 ;
    min_trans = mean_trans - delta_trans ;
    max_trans = mean_trans + delta_trans ;

    /* we've rotated transform to old max */
    x_max_rot = y_max_rot = z_max_rot = 0.0 ;

    mean_angle = (max_angle + min_angle) / 2 ;
    delta_rot = (max_angle-min_angle)/4 ;
    min_angle = mean_angle - delta_rot ;
    max_angle = mean_angle + delta_rot ;
  }

  MatrixFree(&m_x_rot) ;
  MatrixFree(&m_y_rot) ;
  MatrixFree(&m_z_rot) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_tmp) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_tmp2) ;
  MatrixFree(&m_trans) ;
  MatrixFree(&m_tmp3) ;

  return(max_log_p) ;
}

static int
mark_gcas_classes(GCA_SAMPLE *gcas, int nsamples)
{
  int i ;

  for (i = 0 ; i < nsamples ; i++)
  {
    if (IS_GRAY_MATTER(gcas[i].label))
    {
      gcas[i].tissue_class = GM_CLASS ;
    }
    else if (IS_WHITE_MATTER(gcas[i].label))
    {
      gcas[i].tissue_class = WM_CLASS ;
    }
    else if (IS_FLUID(gcas[i].label))
    {
      gcas[i].tissue_class = FLUID_CLASS ;
    }
    else
    {
      gcas[i].tissue_class = OTHER_CLASS ;
    }
  }
  return(NO_ERROR) ;
}

