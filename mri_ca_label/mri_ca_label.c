/**
 * @file  mri_ca_label.c
 * @brief anisotropic nonstationary markov random field labeling
 *
 * Program for computing the MAP segmentation modeling the labeling as an
 * anisotropic nonstationary markov random field (based on manually labeled
 * data compiled into an atlas and stored in a .gca file)
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/05/21 15:15:06 $
 *    $Revision: 1.99 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <sys/utsname.h>
#include <unistd.h>


#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "gca.h"
#include "mri_conform.h"
#include "transform.h"
#include "gcamorph.h"
#include "cma.h"
#include "histo.h"
#include "tags.h"
#include "mrinorm.h"
#include "version.h"
#include "fsinit.h"

static int remove_cerebellum = 0 ;
static int remove_lh = 0 ;
static int remove_rh = 0 ;

static int GCAremoveWMSA(GCA *gca) ;
static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static char *save_gca_fname = NULL ;

static float Glabel_scales[MAX_CMA_LABELS] ;
static float Glabel_offsets[MAX_CMA_LABELS] ;

#define MAX_READS 100
static int nreads = 0 ;
static char *read_intensity_fname[MAX_READS] ;
static float regularize = 0 ;
static float regularize_mean = 0 ;
static int avgs = 0 ;
static int norm_PD = 0;
static int map_to_flash = 0 ;

static int reclassify_unlikely = 0 ;
static int unlikely_wsize = 9 ;
static float unlikely_sigma = .5 ;
static float unlikely_mah_dist_thresh = 4 ;  // more than this many sigmas from mean means unlikely

static int wmsa = 0 ;   // apply wmsa postprocessing (using T2/PD data)
static int nowmsa = 0 ; // remove all wmsa labels from the atlas

static int handle_expanded_ventricles = 0;

static int renorm_with_histos = 0 ;

static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;

#if 0
static int load_val_vector(VECTOR *v_means,
                           MRI *mri_inputs,
                           int x, int y, int z) ;
static double compute_conditional_density(MATRIX *m_inv_cov,
                                          VECTOR *v_means,
                                          VECTOR *v_vals) ;
#endif

int is_possible( GCA *gca,
                 MRI *mri,
                 TRANSFORM *transform,
                 int x, int y, int z, int label );

int insert_thin_temporal_white_matter( MRI *mri_in,
                                       MRI *mri_labeled,
                                       GCA *gca,
                                       TRANSFORM *transform,
                                       GCA *gca_all );

int distance_to_label( MRI *mri_labeled, int label, int x,
                       int y, int z, int dx, int dy,
                       int dz, int max_dist );

int preprocess( MRI *mri_in,
                MRI *mri_labeled,
                GCA *gca,
                TRANSFORM *transform,
                MRI *mri_fixed );

int edit_hippocampus( MRI *mri_in,
                      MRI *mri_labeled,
                      GCA *gca,
                      TRANSFORM *transform,
                      MRI *mri_fixed );

int edit_amygdala( MRI *mri_in,
                   MRI *mri_labeled,
                   GCA *gca,
                   TRANSFORM *transform,
                   MRI *mri_fixed );

int MRIcountNbhdLabels( MRI *mri,
                        int x, int y, int z, int label );

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


static char *tl_gca_fname = NULL ;
static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static char *tissue_parms_fname = NULL ;
static char *mask_volume_fname = NULL ;
static int histo_norm = 0 ;

char *Progname ;
static void usage_exit(int code) ;

static int novar = 0 ;
static char *renormalization_fname = NULL ;
static int renormalize_wsize = 0 ;
static int renormalize_iter = 0 ;
static int renormalize_new = 0 ;
static int renormalize_align = 0 ;
static int no_old_renormalize = 1;
static int filter = 0 ;
static float pthresh = .7 ;
#if 0
static float thresh = 0.5 ;
#endif
static char *reg_fname = NULL ;
static char *read_fname =  NULL ;

static char *wm_fname = NULL ;
static int fixed_flag = 0 ;
static char *heq_fname = NULL ;
static int max_iter = 200 ;
static int mle_niter = 2 ;
static int no_gibbs = 0 ;
static int anneal = 0 ;
static char *mri_fname = NULL ;
static int hippocampus_flag = 1 ;

#define CMA_PARCELLATION  0
static int parcellation_type = CMA_PARCELLATION ;

MRI *insert_wm_segmentation( MRI *mri_labeled, MRI *mri_wm,
                             int parcellation_type, int fixed_flag,
                             GCA *gca, TRANSFORM *transform );

extern char *gca_write_fname ;
extern int gca_write_iterations ;

//static int expand_flag = TRUE ;
static int expand_flag = FALSE ;
static int expand_ventricle_flag = FALSE ;
static int conform_flag = FALSE ;
struct utsname uts;
char *cmdline2, cwd[2000];

int main(int argc, char *argv[])
{
  char         **av ;
  int          ac, nargs, extra ;
  char         *in_fname, *out_fname,  *gca_fname, *xform_fname ;
  MRI          *mri_inputs, *mri_labeled, *mri_fixed = NULL, *mri_tmp ;
  int          msec, minutes, seconds, ninputs, input ;
  struct timeb start ;
  GCA          *gca ;
  TRANSFORM     *transform ;

  char cmdline[CMD_LINE_LEN] ;

  FSinit() ;
  make_cmd_version_string
  (argc, argv,
   "$Id: mri_ca_label.c,v 1.99 2012/05/21 15:15:06 fischl Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mri_ca_label.c,v 1.99 2012/05/21 15:15:06 fischl Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;
  cmdline2 = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  printf("sysname  %s\n",uts.sysname);
  printf("hostname %s\n",uts.nodename);
  printf("machine  %s\n",uts.machine);
  printf("\n");
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("cd %s\n",cwd);
  printf("%s\n",cmdline2);
  printf("\n");
  fflush(stdout);  fflush(stderr);

  setRandomSeed(-1L) ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (getenv("BUILD_GCA_HISTO") != NULL)
  {
    int *counts, i, max_i ;

    fprintf(stderr, "reading gca from %s...\n", argv[1]) ;
    gca = GCAread(argv[1]) ;
    if (!gca)
      ErrorExit(ERROR_NOFILE, "%s: could not read classifier array from %s",
                Progname, argv[1]) ;


    if (gca->flags & GCA_NO_MRF)
      ErrorExit(ERROR_BADPARM, "%s: gca %s built without markov priors",
                Progname, argv[1]) ;

    GCAhisto(gca, 100, &counts) ;

    max_i = 0 ;
    for (i = 1 ; i < 100 ; i++)
      if (counts[i] > 0)
      {
        max_i = i ;
      }

    for (i = 1 ; i <= max_i ; i++)
    {
      printf("%d %d\n", i, counts[i]) ;
    }
    exit(1) ;
  }

  if (argc < 5)
  {
    usage_exit(1) ;
  }

  in_fname = argv[1] ;
  xform_fname = argv[argc-3];
  gca_fname = argv[argc-2] ;
  out_fname = argv[argc-1] ;
  ninputs = argc-4 ;

  printf("reading %d input volumes...\n", ninputs) ;

  /*  fprintf(stderr,
      "mri_inputs read: xform %s\n", mri_inputs->transform_fname) ;*/
  printf("reading classifier array from %s...\n", gca_fname) ;
  fflush(stdout);  fflush(stderr);
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier array from %s",
              Progname, gca_fname) ;

  if (remove_lh)
  {
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
  if (nowmsa)
  {
    GCAremoveWMSA(gca) ;
  }
  extra = 0 ;
  if (gca->flags & GCA_XGRAD)
  {
    extra += ninputs ;
  }
  if (gca->flags & GCA_YGRAD)
  {
    extra += ninputs ;
  }
  if (gca->flags & GCA_ZGRAD)
  {
    extra += ninputs ;
  }

  if (gca->ninputs != (ninputs+extra) &&
      !map_to_flash &&
      gca->type != GCA_FLASH &&
      gca->type != GCA_PARAM)
    ErrorExit
    (ERROR_BADPARM,
     "%s: gca requires %d inputs, %d specified on command line",
     Progname, gca->ninputs, ninputs) ;

  if (avgs)
  {
    GCAmeanFilterConditionalDensities(gca, avgs) ;
  }

  // gathering inputs
  for (input = 0 ; input < ninputs ; input++)
  {
    in_fname = argv[1+input] ;
    printf("reading input volume from %s...\n", in_fname) ;
    fflush(stdout);  fflush(stderr);
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, in_fname) ;

    if (alpha > 0)
    {
      mri_tmp->flip_angle = alpha ;
    }
    if (TR > 0)
    {
      mri_tmp->tr = TR ;
    }
    if (TE > 0)
    {
      mri_tmp->te = TE ;
    }

    TRs[input] = mri_tmp->tr ;
    fas[input] = mri_tmp->flip_angle ;
    TEs[input] = mri_tmp->te ;
    if (conform_flag)
    {
      MRI *mri_tmp2 ;

      mri_tmp2 = MRIconform(mri_tmp) ;
      mri_tmp = mri_tmp2 ;
    }

    if (input == 0)
    {
      mri_inputs =
        MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                         mri_tmp->type, ninputs+extra) ;
      if (!mri_inputs)
        ErrorExit
        (ERROR_NOMEMORY,
         "%s: could not allocate input volume %dx%dx%dx%d",
         mri_tmp->width, mri_tmp->height, mri_tmp->depth,ninputs) ;
      MRIcopyHeader(mri_tmp, mri_inputs) ;
    }

    if (filter)
    {
      MRI *mri_dir, /**mri_grad,*/ *mri_kernel, *mri_smooth ;

      mri_kernel = MRIgaussian1d(1, 15) ;
      mri_smooth = MRIconvolveGaussian(mri_tmp, NULL, mri_kernel) ;
      MRIfree(&mri_kernel) ;
#if 0
      mri_grad = MRIsobel(mri_smooth, NULL, NULL) ;
      mri_dir = MRIdirectionMapUchar(mri_grad, NULL, 5) ;
      MRIfree(&mri_grad) ;
      MRIwrite(mri_dir, "dir.mgz") ;
#else
      mri_dir = MRIgradientDir2ndDerivative(mri_tmp, NULL, 5) ;
      MRIscaleAndMultiply(mri_tmp, 128.0, mri_dir, mri_tmp) ;
      MRIwrite(mri_dir, "lap.mgz") ;
      MRIwrite(mri_tmp, "filtered.mgz") ;
#endif
      MRIfree(&mri_dir) ;
      MRIfree(&mri_smooth) ;
      exit(1) ;
    }
    MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  MRIaddCommandLine(mri_inputs, cmdline) ;
  fflush(stdout);  fflush(stderr);

  // -example fname  option
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
    MRIeraseBorderPlanes(mri_seg, 1) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ;
    MRIfree(&mri_T1) ;
  }
  // -flash_parms fname option
  if (tissue_parms_fname)   /* use FLASH forward model */
  {
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_inputs) ;
  }
  //  -mri fname option
  if (mri_fname)
  {
    GCAbuildMostLikelyVolume(gca, mri_inputs) ;
    MRIwrite(mri_inputs, mri_fname) ;
    exit(0) ;
  }
  // -renorm fname option
  GCAapplyRenormalization(gca, Glabel_scales, Glabel_offsets, 0) ;
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
  fflush(stdout);  fflush(stderr);
  //
  if (gca->type == GCA_FLASH)
  {
    GCA *gca_tmp ;
    int need_map_flag = 0;
    int n;

    if (gca->ninputs != mri_inputs->nframes)
    {
      need_map_flag = 1;
    }
    else
    {
      for (n = 0 ; n < mri_inputs->nframes; n++)
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

    if (need_map_flag)
    {
      printf("mapping %d-dimensional flash atlas "
             "into %d-dimensional input space\n", gca->ninputs, ninputs) ;

      gca_tmp = GCAcreateFlashGCAfromFlashGCA
                (gca, TRs, fas, TEs, mri_inputs->nframes) ;
      GCAfree(&gca) ;
      gca = gca_tmp ;
    }

    if (novar)
    {
      GCAregularizeCovariance(gca,1.0);
    }
    //      GCAunifyVariance(gca) ;

    GCAhistoScaleImageIntensities(gca, mri_inputs, 1) ;
  }
  // -flash option
  if (map_to_flash || gca->type == GCA_PARAM)
  {
    GCA *gca_tmp ;

    if (novar)
    {
      GCAregularizeCovariance(gca,1.0);
    }
    //      GCAunifyVariance(gca) ;

    gca_tmp = GCAcreateFlashGCAfromParameterGCA
              (gca, TRs, fas, TEs,
               mri_inputs->nframes, GCA_DEFAULT_NOISE_PARAMETER) ;
    GCAfree(&gca) ;
    gca = gca_tmp ;
    GCAhistoScaleImageIntensities(gca, mri_inputs, 1) ;
  }
  else if (histo_norm)
  {
    GCAhistoScaleImageIntensities(gca, mri_inputs, 1) ;
  }


  if (gca->flags & GCA_GRAD)
  {
    int i, start = ninputs ;
    MRI *mri_kernel, *mri_smooth, *mri_grad, *mri_tmp ;

    mri_kernel = MRIgaussian1d(1.0, 30) ;
    mri_smooth = MRIconvolveGaussian(mri_inputs, NULL, mri_kernel) ;

    if (mri_inputs->type != MRI_FLOAT)
    {
      mri_tmp = MRISeqchangeType(mri_inputs, MRI_FLOAT, 0, 0, 1) ;
      MRIfree(&mri_inputs) ;
      mri_inputs = mri_tmp ;
    }
    start = ninputs ;
    if (gca->flags & GCA_XGRAD)
    {
      for (i = 0 ; i < ninputs ; i++)
      {
        mri_grad = MRIxSobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }
    if (gca->flags & GCA_YGRAD)
    {
      for (i = 0 ; i < ninputs ; i++)
      {
        mri_grad = MRIySobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }
    if (gca->flags & GCA_ZGRAD)
    {
      for (i = 0 ; i < ninputs ; i++)
      {
        mri_grad = MRIzSobel(mri_smooth, NULL, i) ;
        MRIcopyFrame(mri_grad, mri_inputs, 0, start+i) ;
        MRIfree(&mri_grad) ;
      }
      start += ninputs ;
    }

    MRIfree(&mri_kernel) ;
    MRIfree(&mri_smooth) ;
  }

  if (novar)
  {
    GCAregularizeCovariance(gca,1.0);
  }
  //    GCAunifyVariance(gca) ;
  if (regularize > 0)
  {
    GCAregularizeCovariance(gca, regularize) ;
  }

  if (stricmp(xform_fname, "none"))
  {
    GCA_MORPH *gcam;
    printf("reading transform from %s...\n", xform_fname) ;
    transform = TransformRead(xform_fname) ;
    if (!transform)
    {
      ErrorExit(ERROR_NOFILE, "%s: could not open transform", xform_fname) ;
    }

    if (TransformFileNameType(xform_fname) == MORPH_3D_TYPE)
    {
      gcam = (GCA_MORPH *)(transform->xform);
      printf("Atlas used for the 3D morph was %s\n", gcam->atlas.fname);
    }

    TransformInvert(transform, mri_inputs) ;
  }
  else
  {
    transform = TransformAlloc(LINEAR_VOX_TO_VOX, NULL) ;
  }

  if (norm_PD)
  {
    GCAnormalizePD(gca, mri_inputs, transform) ;
  }

  if (Ggca_x >= 0 && Gx < 0)
  {
    GCAsourceVoxelToNode(gca, mri_inputs, transform, Ggca_x, Ggca_y, Ggca_z,&Gx, &Gy, &Gz) ;
    printf("source voxel (%d, %d, %d) maps to node (%d, %d, %d)\n", Ggca_x, Ggca_y, Ggca_z, Gx, Gy, Gz) ;
    GCAdump(gca, mri_inputs, Ggca_x, Ggca_y, Ggca_z, transform, stdout, 0) ;
  }
  if (nreads > 0)
  {
    float label_scales[MAX_CMA_LABELS], label_offsets[MAX_CMA_LABELS] ;
    float label_scales_total[MAX_CMA_LABELS],
          label_offsets_total[MAX_CMA_LABELS];
    char  *fname ;
    int   i, l ;

    memset(label_scales_total, 0, sizeof(label_scales_total)) ;
    memset(label_offsets_total, 0, sizeof(label_offsets_total)) ;

    for (i = 0 ; i < nreads ; i++)
    {
      fname = read_intensity_fname[i] ;
      printf("reading label scales and offsets from %s\n", fname) ;
      GCAreadLabelIntensities(fname, label_scales, label_offsets) ;
      for (l = 0; l < MAX_CMA_LABELS ; l++)
      {
        label_scales_total[l] += label_scales[l] ;
        label_offsets_total[l] += label_offsets[l] ;
      }
    }
    for (l = 0; l < MAX_CMA_LABELS ; l++)
    {
      label_scales_total[l] /= (float)nreads ;
      label_offsets_total[l] /= (float)nreads ;
    }

    GCAapplyRenormalization(gca, label_scales_total, label_offsets_total, 0) ;
  }
  if (heq_fname)
  {
    MRI *mri_eq ;

    mri_eq = MRIread(heq_fname) ;
    if (!mri_eq)
      ErrorExit(ERROR_NOFILE,
                "%s: could not read histogram equalization volume %s",
                Progname, heq_fname) ;
    MRIhistoEqualize(mri_inputs, mri_eq, mri_inputs, 30, 170) ;
    MRIfree(&mri_eq) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      fprintf(stderr, "writing equalized volume to %s...\n", "heq.mgz") ;
      MRIwrite(mri_inputs, "heq.mgz") ;
    }
  }
  fflush(stdout);  fflush(stderr);

  GCAfixSingularCovarianceMatrices(gca) ;
  if (renorm_with_histos)
  {
    GCAmapRenormalizeWithHistograms
      (gca, mri_inputs, transform, NULL, "ca_label",NULL,NULL,NULL,NULL);
  }
  if (read_fname != NULL && reg_fname == NULL)  /* use given segmentation */
  {
    //read in initial segmentation from file read_fname
    mri_labeled = MRIread(read_fname) ;
    if (!mri_labeled)
      ErrorExit(ERROR_NOFILE, "%s: could not read segmentation from %s",
                Progname, read_fname) ;
    if (Ggca_x >= 0)
      printf("label(%d, %d, %d) = %s (%d), norm=%2.0f\n",
             Ggca_x, Ggca_y, Ggca_z,
             cma_label_to_name(MRIvox(mri_labeled, Ggca_x, Ggca_y, Ggca_z)),
             MRIvox(mri_labeled, Ggca_x, Ggca_y, Ggca_z),
             MRIgetVoxVal(mri_inputs, Ggca_x, Ggca_y, Ggca_z, 0)) ;
    if (regularize_mean > 0)
    {
      GCAregularizeConditionalDensities(gca, regularize_mean) ;
    }
    if (Ggca_x >= 0)
    {
      GCAdump(gca, mri_inputs, Ggca_x, Ggca_y, Ggca_z, transform, stdout, 0) ;
    }
    GCAreclassifyUsingGibbsPriors
    (mri_inputs, gca, mri_labeled, transform, max_iter,
     mri_fixed, 0, NULL);
  }
  else    // don't read old one in - create new labeling
  {
    if (reg_fname == NULL)   //so read_fname must be NULL too
    {
      printf("labeling volume...\n") ;
      // create labeled volume by MAP rule
      mri_labeled = GCAlabel(mri_inputs, gca, NULL, transform) ;
      // -wm fname option
      if (wm_fname)
      {
        MRI *mri_wm ;

        mri_wm = MRIread(wm_fname) ;
        if (!mri_wm)
          ErrorExit
          (ERROR_NOFILE,
           "%s: could not read wm segmentation from %s",
           Progname, wm_fname) ;
        // put wm into fixed
        mri_fixed = insert_wm_segmentation
                    (mri_labeled,mri_wm,parcellation_type,
                     fixed_flag, gca, transform);
        if (DIAG_VERBOSE_ON)
        {
          fprintf(stderr,
                  "writing patched labeling to %s...\n", out_fname) ;
          MRIwrite(mri_labeled, out_fname) ;
        }
        MRIfree(&mri_wm) ;
      } //if (wm_fname)
      else
        // just clone the labeled one
      {
        mri_fixed = MRIclone(mri_labeled, NULL) ;
      }

      if (gca_write_iterations != 0)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s_pre.mgz", gca_write_fname) ;
        printf("writing snapshot to %s...\n", fname) ;
        MRIwrite(mri_labeled, fname) ;
      }
      // renormalize iteration
      if (renormalize_align)
      {
        FILE *logfp ;
        char base_name[STRLEN] ;

        FileNameOnly(out_fname, base_name) ;
        FileNameRemoveExtension(base_name, base_name) ;
        if (Gdiag & DIAG_WRITE)
        {
          char fname[STRLEN] ;
          strcpy(fname, base_name) ;
          strcat(fname, ".log") ;
          logfp = fopen(fname, "w") ;
        }
        else
        {
          logfp = NULL ;
        }
        if (!no_old_renormalize)
        {
          GCAmapRenormalize(gca, mri_inputs, transform) ;
        }

        // GCA Renormalization with Alignment:
        // run it twice so that the histograms overlap on the 2nd run
        float label_scales[MAX_CMA_LABELS],
          label_offsets[MAX_CMA_LABELS],
          label_peaks[MAX_CMA_LABELS];
        int label_computed[MAX_CMA_LABELS];
        // initial call (returning the label_* infos)
        GCAcomputeRenormalizationWithAlignment
        (gca, mri_inputs, transform,
         logfp, base_name, NULL, handle_expanded_ventricles,
         label_scales,label_offsets,label_peaks,label_computed) ;
        // sequential call gets passed the results from first call
        // will overwrite the intensity.txt file with combinded results
        GCAseqRenormalizeWithAlignment
        (gca, mri_inputs, transform,
         logfp, base_name, NULL, handle_expanded_ventricles,
         label_scales,label_offsets,label_peaks,label_computed) ;
        if (regularize_mean > 0)
        {
          GCAregularizeConditionalDensities(gca, regularize_mean) ;
        }
        if (save_gca_fname)
        {
          GCAwrite(gca, save_gca_fname) ;
        }
        GCAlabel(mri_inputs, gca, mri_labeled, transform) ;
        {
          MRI *mri_imp ;
          mri_imp = GCAmarkImpossible(gca, mri_labeled, NULL, transform) ;
          if (Gdiag & DIAG_WRITE)
          {
            MRIwrite(mri_imp, "gca_imp.mgz") ;
          }
          MRIfree(&mri_imp) ;
        }
      }
      else if (renormalize_iter > 0)
      {
        if (renormalize_new)
        {
          GCAmapRenormalizeByClass(gca, mri_inputs, transform) ;
        }
        else
        {
          GCAmapRenormalize(gca, mri_inputs, transform) ;
        }
        if (save_gca_fname)
        {
          GCAwrite(gca, save_gca_fname) ;
        }
        printf("relabeling volume...\n") ;
        if (regularize_mean > 0)
        {
          GCAregularizeConditionalDensities(gca, regularize_mean) ;
        }
        GCAlabel(mri_inputs, gca, mri_labeled, transform) ;
      }
      else if (regularize_mean > 0)
      {
        GCAregularizeConditionalDensities(gca, regularize_mean) ;
      }
      preprocess(mri_inputs, mri_labeled, gca, transform, mri_fixed) ;
      if (fixed_flag == 0)
      {
        MRIfree(&mri_fixed) ;
      }
    } //if(reg_fname == NULL)
    else  /* processing longitudinal data (old version) */
    {
      ErrorExit(ERROR_BADPARM, "%s ERROR: the -l option is currently not supported. Debugging and testing needed.\n",
                Progname) ;
      // Now, both an intial seg and a transformation are given
      // suppose the seg is from tp1 and
      // transform is to align tp1 to current tp
      TRANSFORM *transform_long ;
      MRI       *mri_tmp ;
      VOL_GEOM vgm_in;
      mri_labeled = MRIread(read_fname) ;
      if (!mri_labeled)
        ErrorExit(ERROR_NOFILE, "%s: could not read segmentation from %s",
                  Progname, read_fname) ;
      printf("applying transform %s to previously "
             "computed segmentation %s\n",
             reg_fname, read_fname) ;
      transform_long = TransformRead(reg_fname) ;
      if (transform_long == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not open registration file %s",
                  Progname, reg_fname) ;

      {
        LTA *lta = (LTA *)transform_long->xform ;
        if (lta->xforms[0].src.valid == 0)
        {
          LTAmodifySrcDstGeom(lta, mri_labeled, NULL); //add src info
        }
        if (lta->xforms[0].dst.valid == 0)
        {
          LTAmodifySrcDstGeom(lta, NULL, mri_inputs); //add dst information
        }
        LTAchangeType(lta, LINEAR_VOX_TO_VOX) ;

        getVolGeom(mri_inputs, &vgm_in);
        if (vg_isEqual(&lta->xforms[0].dst, &vgm_in) == 0)
        {
          printf("%s: WARNING: dst volume of lta "
                 "doesn't match that of input volume\n",Progname);
          printf("Volume geometry for lta-dst:\n");
          vg_print(&lta->xforms[0].dst);
          printf("Volume geometry for input volume is:\n");
          vg_print(&vgm_in);
        }

        getVolGeom(mri_labeled, &vgm_in);
        if (vg_isEqual(&lta->xforms[0].src, &vgm_in) == 0)
        {
          printf("%s: WARNING: src volume of lta "
                 "doesn't match that of tp1 label volume\n",Progname);
          printf("Volume geometry for lta-src:\n");
          vg_print(&lta->xforms[0].src);
          printf("Volume geometry for tp1 label volume is:\n");
          vg_print(&vgm_in);

        }

        transform_long->type = LINEAR_VOX_TO_VOX ;
      }
      if (transform_long->type != LINEAR_VOX_TO_VOX)
        ErrorExit
        (ERROR_BADPARM,
         "%s: transform type (%d) must be LINEAR_VOX_TO_VOX",
         Progname, transform_long->type) ;
      mri_tmp = MRIalloc
                (mri_inputs->width,
                 mri_inputs->height,
                 mri_inputs->depth,
                 mri_labeled->type);
      MRIcopyHeader(mri_inputs, mri_tmp);

      mri_tmp =
        MRIlinearTransformInterp
        (mri_labeled,
         mri_tmp,
         ((LTA *)(transform_long->xform))->xforms[0].m_L,
         SAMPLE_NEAREST) ;

      MRIfree(&mri_labeled) ;
      mri_labeled = mri_tmp ;

      if (Ggca_x >= 0)
        printf("label(%d, %d, %d) = %s (%d)\n",
               Ggca_x, Ggca_y, Ggca_z,
               cma_label_to_name(MRIvox
                                 (mri_labeled, Ggca_x, Ggca_y, Ggca_z)),
               MRIvox(mri_labeled, Ggca_x, Ggca_y, Ggca_z)) ;
      TransformFree(&transform_long) ;
      if (gca_write_iterations != 0)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s_pre.mgz", gca_write_fname) ;
        printf("writing snapshot to %s...\n", fname) ;
        MRIwrite(mri_labeled, fname) ;
      }

      if (renormalize_align)
      {
        FILE *logfp ;
        char base_name[STRLEN] ;

        FileNameOnly(out_fname, base_name) ;
        FileNameRemoveExtension(base_name, base_name) ;
        if (Gdiag & DIAG_WRITE)
        {
          char fname[STRLEN] ;
          strcpy(fname, base_name) ;
          strcat(fname, ".log") ;
          logfp = fopen(fname, "w") ;
        }
        else
        {
          logfp = NULL ;
        }
        GCAmapRenormalizeWithAlignment
        (gca, mri_inputs, transform,
         logfp, base_name, NULL, handle_expanded_ventricles) ;
      }
      if (save_gca_fname)
      {
        GCAwrite(gca, save_gca_fname) ;
      }
      //This normalization seems to bias WM to gets larger
#if 0
      if (ninputs == 1)
      {
        GCArenormalizeToExample(gca, mri_labeled, mri_inputs);
      }
      else
      {
        printf("Warning: should renormalize GCA, but current "
               "code only support this for single-channel atlas\n");
      }
#else
      // renormalize iteration
      if (renormalize_iter > 0)
      {
        if (renormalize_new)
        {
          GCAmapRenormalizeByClass(gca, mri_inputs, transform) ;
        }
        else
        {
          GCAmapRenormalize(gca, mri_inputs, transform) ;
        }
      }

#endif

    } //else /* processing long data */

    if (!no_gibbs)
    {
      if (anneal)
      {
        GCAanneal(mri_inputs, gca, mri_labeled, transform, max_iter) ;
      }
      else
      {
        if (Ggca_x >= 0)
        {
          GCAdump(gca, mri_inputs, Ggca_x, Ggca_y, Ggca_z, transform, stdout, 0) ;
        }
        GCAreclassifyUsingGibbsPriors
        (mri_inputs, gca, mri_labeled, transform, max_iter,
         mri_fixed, 0, NULL);
      }
    }
  }

  if (read_fname == NULL && 0)
  {
    GCAmaxLikelihoodBorders(gca, mri_inputs, mri_labeled, mri_labeled,transform,mle_niter, 5.0);
  }
  if (expand_ventricle_flag)
  {
    GCAexpandVentricle(gca, mri_inputs, mri_labeled, mri_labeled, transform,
                       Left_Lateral_Ventricle) ;
    GCAexpandVentricle(gca, mri_inputs, mri_labeled, mri_labeled, transform,
                       Right_Lateral_Ventricle) ;
#if 0
    GCAexpandVentricle(gca, mri_inputs, mri_labeled, mri_labeled, transform,
                       Left_Inf_Lat_Vent) ;
    GCAexpandVentricle(gca, mri_inputs, mri_labeled, mri_labeled, transform,
                       Right_Inf_Lat_Vent) ;
#endif
  }
  if (expand_flag)
  {
    GCAexpandCortex(gca, mri_inputs, mri_labeled, mri_labeled, transform) ;
  }

  if (mask_volume_fname)
  {
    MRI *mri_mask ;
    printf("reading volume %s for masking...\n", mask_volume_fname) ;
    mri_mask = MRIread(mask_volume_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not read mask volume from %s",
                Progname, mask_volume_fname) ;
    // if mask has some value > 1, then keep the original
    // if      has some value < 1, then the value is set to 0
    MRIthresholdMask(mri_labeled, mri_mask, mri_labeled, 1, 0) ;
    MRIfree(&mri_mask) ;
  }

  if (gca_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_post.mgz", gca_write_fname) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
  }

  if (reclassify_unlikely)
  {
    MRI *mri_sigma ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_labeled, "aseg_before.mgz") ;
    }
    mri_sigma = GCAcomputeOptimalScale(gca, transform, mri_inputs, mri_labeled, NULL,
                                       2*(gca->node_spacing+1), 0.5, 2, .5) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_sigma, "sigma.mgz") ;
    }

#if 0
    GCAreclassifyUnlikelyVoxels(gca, transform, mri_inputs, mri_labeled,
                                mri_labeled, unlikely_mah_dist_thresh,
                                unlikely_wsize, unlikely_sigma) ;
#else
    GCAreclassifyVoxelsAtOptimalScale(gca, transform, mri_inputs, mri_labeled,
                                      mri_labeled,
                                      mri_sigma, unlikely_wsize) ;
#endif
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_labeled, "aseg_after.mgz") ;
    }
    MRIfree(&mri_sigma) ;
  }
  if (read_fname == NULL)
  {
    GCAconstrainLabelTopology(gca, mri_inputs, mri_labeled, mri_labeled, transform) ;
  }
  if (wmsa)
  {
    MRI *mri_tmp ;
    mri_tmp = GCAlabelWMandWMSAs(gca, mri_inputs, mri_labeled, NULL, transform) ;
    MRIfree(&mri_labeled) ;
    mri_labeled = mri_tmp ;
  }

  mri_labeled->ct = gca->ct ;  // embed color table in output volume
  /*  GCAfree(&gca) ; */
  MRIfree(&mri_inputs) ;
#if 0
  if (filter)
  {
    MRI *mri_tmp ;

    printf("filtering labeled volume...\n") ;
    mri_tmp = MRIthreshModeFilter(mri_labeled, NULL, filter, thresh) ;
    MRIfree(&mri_labeled) ;
    mri_labeled = mri_tmp ;
  }
#endif

  printf("writing labeled volume to %s...\n", out_fname) ;
  if (MRIwrite(mri_labeled, out_fname) != NO_ERROR)
  {
    ErrorExit(Gerror, "%s: MRIwrite(%s) failed", Progname, out_fname) ;
  }

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("auto-labeling took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  return(0) ;
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

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "NOGIBBS"))
  {
    no_gibbs = 1 ;
    printf("disabling gibbs priors...\n") ;
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
  else if (!strcmp(option, "NOCEREBELLUM"))
  {
    remove_cerebellum = 1 ;
    printf("removing cerebellum from atlas\n") ;
  }
  else if (!stricmp(option, "nowmsa"))
  {
    nowmsa = 1 ;
    printf("disabling WMSA labels\n") ;
  }
  else if (!stricmp(option, "read_intensities") || !stricmp(option, "ri"))
  {
    read_intensity_fname[nreads] = argv[2] ;
    nargs = 1 ;
    printf("reading intensity scaling from %s...\n", read_intensity_fname[nreads]) ;
    nreads++ ;
    if (nreads > MAX_READS)
    {
      ErrorExit(ERROR_UNSUPPORTED, "%s: too many intensity files specified (max %d)", Progname, MAX_READS);
    }
  }
  else if (!stricmp(option, "WM"))
  {
    wm_fname = argv[2] ;
    nargs = 1 ;
    printf("inserting white matter segmentation from %s...\n", wm_fname) ;
  }
  else if (!stricmp(option, "SD"))
  {
    setenv("SUBJECTS_DIR",argv[2],1) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-HELP")||!stricmp(option, "-USAGE"))
  {
    usage_exit(0) ;
  }
  else if (!stricmp(option, "histo-norm"))
  {
    printf("using prior subject histograms for initial GCA renormalization\n") ;
    renorm_with_histos = 1 ;
  }
  else if (!stricmp(option, "SAVE_GCA"))
  {
    save_gca_fname = argv[2] ;
    nargs = 1 ;
    printf("saving renormalized gca to %s\n", save_gca_fname) ;
  }
  else if (!stricmp(option, "RELABEL_UNLIKELY"))
  {
    reclassify_unlikely = atoi(argv[2]) ;
    unlikely_wsize = atoi(argv[3]) ;
    unlikely_sigma = atof(argv[4]) ;
    unlikely_mah_dist_thresh = atof(argv[5]) ;
    if (reclassify_unlikely)
      printf("relabeling unlikely voxels more than %2.1f sigmas from label mean, with sigma=%2.1fmm, wsize=%dmm\n",
             unlikely_mah_dist_thresh, unlikely_sigma, unlikely_wsize) ;
    else
    {
      printf("not relabeling unlikely voxels\n") ;
    }
    nargs = 4 ;
  }
  else if (!stricmp(option, "CONFORM"))
  {
    conform_flag = TRUE ;
    printf("resampling input volume(s) to be 256^3 and 1mm^3\n") ;
  }
  else if (!stricmp(option, "WMSA"))
  {
    wmsa = 1 ;
    printf("relabeling wm and wmsa in postprocessing\n") ;
  }
  else if (!stricmp(option, "NORMPD"))
  {
    norm_PD = TRUE ;
    printf("normalizing PD image (2nd input) to GCA means[1]\n") ;
  }
  else if (!stricmp(option, "bigventricles"))
  {
    handle_expanded_ventricles = 1 ;
    printf("handling expanded ventricles...\n") ;
  }
  else if (!stricmp(option, "nobigventricles"))
  {
    handle_expanded_ventricles = 0 ;
    printf("not handling expanded ventricles...\n") ;
  }
  else if (!stricmp(option, "write_probs"))
  {
    G_write_probs = argv[2] ;
    nargs = 1 ;
    printf("writing label probabilities to %s...\n", G_write_probs) ;
  }
  else if (!stricmp(option, "TL"))
  {
    tl_gca_fname = argv[2] ;
    nargs = 1 ;
    printf("using gca %s to label thin temporal lobe...\n", tl_gca_fname) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Ggca_x,Ggca_y,Ggca_z) ;
  }
  else if (!stricmp(option, "DEBUG_NODE"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "DEBUG_PRIOR"))
  {
    Gxp = atoi(argv[2]) ;
    Gyp = atoi(argv[3]) ;
    Gzp = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %d\n", Ggca_label) ;
  }
  else if (!stricmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!stricmp(option, "expand"))
  {
    expand_ventricle_flag = TRUE ;
    printf("expanding ventricles in postprocessing...\n") ;
  }
  else if (!stricmp(option, "EXAMPLE"))
  {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!stricmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!stricmp(option, "PTHRESH"))
  {
    nargs = 1 ;
    pthresh = atof(argv[2]) ;
    printf("using p threshold %2.2f for adaptive renormalization\n",
           pthresh) ;
  }
  else if (!stricmp(option, "NITER"))
  {
    mle_niter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("applying max likelihood for %d iterations...\n", mle_niter) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance in classification\n") ;
  }
  else if (!stricmp(option, "LSCALE"))
  {
    int l ;
    l = atoi(argv[2]) ;
    Glabel_scales[l] = atof(argv[3]) ;
    nargs = 2 ;
    printf("scaling label %s by %2.2f\n", cma_label_to_name(l), Glabel_scales[l]) ;
    for (l = 0 ; l < MAX_CMA_LABELS ; l++)
      if (FZERO(Glabel_scales[l]))
	Glabel_scales[l] = 1.0 ;
  }
  else if (!stricmp(option, "REGULARIZE"))
  {
    regularize = atof(argv[2]) ;
    printf("regularizing covariance to be %2.2f Cclass + %2.2f Cpooled\n",
           (1-regularize),regularize) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "REGULARIZE_MEAN"))
  {
    regularize_mean = atof(argv[2]) ;
    printf("regularizing means to be %2.2f u(global) + %2.2f u(r)\n",
           regularize_mean, 1-regularize_mean) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "align-cross-sequence") ||
           !stricmp(option, "align") ||
           !stricmp(option, "align-after"))
  {
    regularize = .5 ;
    regularize_mean = 0.5 ;
    avgs = 2 ;
    renormalize_align = 1 ;
    printf("renormalizing sequences with structure alignment, "
           "equivalent to:\n") ;
    printf("\t-renormalize\n\t-renormalize_mean "
           "%2.3f\n\t-regularize %2.3f\n",regularize_mean, regularize) ;
  }
  else if (!stricmp(option, "no_old_renormalize"))
  {
    no_old_renormalize  = 1; //this option turns off
    // the initial GCAmapRenormalize() when -align option is used
  }
  else if (!stricmp(option, "cross-sequence") ||
           !stricmp(option, "cross_sequence") ||
           !stricmp(option, "cross-sequence-new") ||
           !stricmp(option, "cross_sequence-new"))
  {
    if (!stricmp(option, "cross-sequence-new") ||
        !stricmp(option, "cross_sequence-new"))
    {
      renormalize_new = 1 ;
    }

    regularize = .5 ;
    renormalize_iter = 1 ;
    renormalize_wsize = 9 ;
    avgs = 2 ;
    printf("labeling across sequences, equivalent to:\n") ;
    printf("\t-renormalize %d %d\n\t-a %d\n\t-regularize %2.3f\n",
           renormalize_iter, renormalize_wsize,
           avgs, regularize) ;
  }
  else if (!stricmp(option, "nohippo"))
  {
    hippocampus_flag = 0 ;
    printf("disabling auto-editting of hippocampus\n") ;
  }
  else if (!stricmp(option, "FWM"))
  {
    fixed_flag = 1 ;
    wm_fname = argv[2] ;
    nargs = 1 ;
    printf("inserting fixed white matter segmentation from %s...\n",
           wm_fname);
  }
  else if (!stricmp(option, "MRI"))
  {
    mri_fname = argv[2] ;
    nargs = 1 ;
    printf("building most likely MR volume and writing to %s...\n",
           mri_fname);
  }
  else if (!stricmp(option, "HEQ"))
  {
    heq_fname = argv[2] ;
    nargs = 1 ;
    printf("reading template for histogram equalization from %s...\n",
           heq_fname) ;
  }
  else if (!stricmp(option, "RENORM"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!stricmp(option, "FLASH"))
  {
    map_to_flash = 1 ;
    printf("using FLASH forward model to predict intensity values...\n") ;
  }
  else if (!stricmp(option, "FLASH_PARMS"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!stricmp(option, "renormalize"))
  {
    renormalize_wsize = atoi(argv[2]) ;
    renormalize_iter = atoi(argv[3]) ;
    nargs = 2 ;
    printf("renormalizing class means %d times after "
           "initial labeling using %d window size\n",
           renormalize_iter, renormalize_wsize) ;
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'L': //longitudinal option
      //option for longitudinal processing
      //read_fname is the label for tp1
      //reg_fname is the linear transform from tp1 to current study
      read_fname = argv[2] ;
      reg_fname = argv[3] ;
      nargs = 2;
      break ;
    case 'R':
      read_fname = argv[2] ;
      printf("reading previously labeled volume from %s...\n", read_fname) ;
      nargs = 1 ;
      break ;
    case 'H':
      histo_norm = 1 ;
      printf("using GCA to histogram normalize input image...\n") ;
      break ;
    case 'A':
      avgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf
      (stderr,
       "applying mean filter %d times to conditional densities...\n", avgs) ;
      break ;
    case 'W':
      gca_write_iterations = atoi(argv[2]) ;
      gca_write_fname = argv[3] ;
      nargs = 2 ;
      printf("writing out snapshots of gibbs process "
             "every %d iterations to %s\n",
             gca_write_iterations, gca_write_fname) ;
      break ;
    case 'M':
      mask_volume_fname = argv[2] ;
      nargs = 1 ;
      printf("using %s to mask final labeling...\n", mask_volume_fname) ;
      break ;
    case 'E':
      expand_flag = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'N':
      max_iter = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting max iterations to %d...\n", max_iter) ;
      break ;
    case 'F':
#if 0
      filter = atoi(argv[2]) ;
      thresh = atof(argv[3]) ;
      nargs = 2 ;
      printf("applying thresholded (%2.2f) mode filter %d times to output of "
             "labelling\n",thresh,filter);
#else
      filter = 1 ;
#endif
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
#include "mri_ca_label.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_ca_label_help_xml, mri_ca_label_help_xml_len);
  exit(code);
}


static int cma_editable_labels[] =
{
  Left_Hippocampus,
  Right_Hippocampus,
  Left_Cerebral_Cortex,
  Left_Cerebral_White_Matter,
  Right_Cerebral_Cortex,
  Right_Cerebral_White_Matter,
#if 1
  Left_Amygdala,
  Right_Amygdala
#endif
} ;
#define NEDITABLE_LABELS sizeof(cma_editable_labels) / sizeof(cma_editable_labels[0])

MRI *
insert_wm_segmentation( MRI *mri_labeled, MRI *mri_wm,
                        int parcellation_type,
                        int fixed_flag, GCA *gca,TRANSFORM *transform )
{
  int      x, y, z, width, depth, height, change_label[1000], n, label,
           nchanged, rh, lh, xn, yn, zn ;
  MRI      *mri_fixed ;
  double   wm_prior ;
  GCA_PRIOR *gcap ;

  mri_fixed = MRIclone(mri_wm, NULL) ;

  memset(change_label, 0, sizeof(change_label)) ;
  for (n = 0 ; n < NEDITABLE_LABELS ; n++)
  {
    change_label[cma_editable_labels[n]] = 1 ;
  }
  width = mri_wm->width ;
  height = mri_wm->height ;
  depth = mri_wm->depth ;

  for (nchanged = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 103 && y == 121 && z == 122)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL ||
            MRIvox(mri_wm,x,y,z)>=200)
        {
          continue ;
        }
        label = MRIvox(mri_labeled, x, y, z) ;
        if (!change_label[label])
        {
          continue ;
        }

        if (!GCAsourceVoxelToNode
            (gca, mri_labeled, transform, x, y, z, &xn, &yn, &zn))
        {
          gcap = getGCAP(gca, mri_labeled, transform, x, y, z) ;
          wm_prior = 0.0 ;
          for (n = 0 ; n < gcap->nlabels ; n++)
            if (gcap->labels[n] == Right_Cerebral_White_Matter ||
                gcap->labels[n] == Left_Cerebral_White_Matter)
            {
              wm_prior = gcap->priors[n] ;
            }
#define PRIOR_THRESH 0.1
#define FIXED_PRIOR  0.5
          if (wm_prior < PRIOR_THRESH)
          {
            continue ;
          }
          lh = MRIcountNbhdLabels(mri_labeled, x, y, z,
                                  Left_Cerebral_White_Matter) ;
          lh += MRIcountNbhdLabels(mri_labeled, x, y, z,
                                   Left_Cerebral_Cortex) ;
          rh = MRIcountNbhdLabels(mri_labeled, x, y, z,
                                  Right_Cerebral_White_Matter) ;
          rh += MRIcountNbhdLabels(mri_labeled, x, y, z,
                                   Right_Cerebral_Cortex) ;
          if (rh > lh)
          {
            label = Right_Cerebral_White_Matter ;
          }
          else
          {
            label = Left_Cerebral_White_Matter ;
          }
          if (label != MRIvox(mri_labeled, x, y, z))
          {
            nchanged++ ;
          }
          MRIvox(mri_labeled, x, y, z) = label ;
#if 0
          if (MRIvox(mri_wm, x, y, z) < 140)
#endif
            MRIvox(mri_fixed, x, y, z) = fixed_flag ;
        }// !GCA
      }
    }
  }

  printf("%d labels changed to wm\n", nchanged) ;
  return(mri_fixed) ;
}


int
MRIcountNbhdLabels( MRI *mri,
                    int x, int y, int z, int label )
{
  int     total, xi, yi, zi, xk, yk, zk ;

  for (total = 0, zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIvox(mri, xi, yi, zi) == label)
        {
          total++ ;
        }
      }
    }
  }

  return(total) ;
}


static int cma_expandable_labels[] =
{
#if 0
  Left_Putamen,
  Right_Putamen,
  Left_Amygdala,
  Right_Amygdala,
  Left_Caudate,
  Right_Caudate,
#endif
  Left_Pallidum,
  Right_Pallidum
#if 0
  Left_Lateral_Ventricle,
  Left_Inf_Lat_Vent,
  Right_Lateral_Ventricle,
  Right_Inf_Lat_Vent
#endif
} ;
#define NEXPANDABLE_LABELS sizeof(cma_expandable_labels) / sizeof(cma_expandable_labels[0])

int
preprocess( MRI *mri_inputs, MRI *mri_labeled,
            GCA *gca, TRANSFORM *transform,
            MRI *mri_fixed )
{
  int i ;

  GCArelabel_cortical_gray_and_white
  (gca, mri_inputs, mri_labeled,mri_labeled,transform);
  if (gca_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_cortex.mgz", gca_write_fname) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
  }
  if (tl_gca_fname)
  {
    GCA  *gca_tl ;

    gca_tl = GCAread(tl_gca_fname) ;
    if (!gca_tl)
      ErrorExit(ERROR_NOFILE, "%s: could not read temporal lobe GCA file %s",
                Progname, tl_gca_fname) ;

    if (renormalize_wsize > 0)  /* gca was renormalized -
                                                             update gca_tl to new means */
    {
      GCArenormalizeFromAtlas(gca_tl, gca) ;
    }

    insert_thin_temporal_white_matter
    (mri_inputs, mri_labeled, gca_tl, transform,gca) ;
    GCAfree(&gca_tl) ;
    if (gca_write_iterations != 0)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s_temporal.mgz", gca_write_fname) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri_labeled, fname) ;
    }
  }

  if (hippocampus_flag)
  {
    edit_hippocampus(mri_inputs, mri_labeled, gca, transform, mri_fixed) ;
    edit_amygdala(mri_inputs, mri_labeled, gca, transform, mri_fixed) ;
  }
  if (expand_flag)
  {
    for (i = 0 ; i < NEXPANDABLE_LABELS ; i++)
      GCAexpandLabelIntoWM
      (gca, mri_inputs, mri_labeled, mri_labeled, transform,
       mri_fixed, cma_expandable_labels[i]) ;
#if 0
    GCAexpandVentricleIntoWM
    (gca, mri_inputs, mri_labeled, mri_labeled, transform,
     mri_fixed,Left_Lateral_Ventricle) ;
    GCAexpandVentricleIntoWM
    (gca, mri_inputs, mri_labeled, mri_labeled, transform,
     mri_fixed,Left_Inf_Lat_Vent) ;
    GCAexpandVentricleIntoWM
    (gca, mri_inputs, mri_labeled, mri_labeled, transform,
     mri_fixed,Right_Lateral_Ventricle) ;
    GCAexpandVentricleIntoWM
    (gca, mri_inputs, mri_labeled, mri_labeled, transform,
     mri_fixed, Right_Inf_Lat_Vent) ;
#endif
  }
#if 0
  GCAmaxLikelihoodBorders
  (gca, mri_inputs, mri_labeled, mri_labeled, transform,2,2.0) ;
#endif
  return(NO_ERROR) ;
}

int
distance_to_label( MRI *mri_labeled, int label,
                   int x, int y, int z,
                   int dx, int dy, int dz, int max_dist )
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ;
    yi = y + d * dy ;
    zi = z + d * dz ;
    xi = mri_labeled->xi[xi] ;
    yi = mri_labeled->yi[yi] ;
    zi = mri_labeled->zi[zi];
    if (MRIvox(mri_labeled, xi, yi, zi) == label)
    {
      break ;
    }
  }

  return(d) ;
}

int
edit_hippocampus( MRI *mri_inputs, MRI *mri_labeled,
                  GCA *gca,
                  TRANSFORM *transform,
                  MRI *mri_fixed )
{
  int   width, height, depth, x, y, z, label, nchanged, dleft, olabel,
        dright, dpos, dant, dup, ddown, dup_hippo, ddown_gray, i, left,
        ddown_ventricle, yi, found_wm, found_hippo, was_wm, was_hippo,
        was_unknown, total_changed, xn, yn, zn ;
  MRI   *mri_tmp ;
  float phippo, pwm ;
  GC1D  *gc ;

  mri_tmp = MRIcopy(mri_labeled, NULL) ;

  width = mri_inputs->width ;
  height = mri_inputs->height ;
  depth = mri_inputs->depth ;

  total_changed = 0 ;
  for (nchanged = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIvox(mri_labeled, x, y, z) ;
        if (!GCAsourceVoxelToNode
            (gca, mri_labeled, transform, x, y, z, &xn, &yn, &zn))
        {
          gc = GCAfindGC(gca, xn, yn, zn, Left_Hippocampus) ;
          if (!gc)
          {
            gc = GCAfindGC(gca, xn, yn, zn, Right_Hippocampus) ;
          }
          if (!gc)
          {
            continue ;
          }

          left = 0 ;
          switch (label)
          {
            /*
               see if there is cortical gray matter
               medial to us. If so, change it to hippocampus.
            */
          case Left_Lateral_Ventricle:
          case Left_Inf_Lat_Vent:    /* medial is in negative
                                        x direction */
          {
            int xi ;

            xi = mri_inputs->xi[x-1] ;
            if (MRIvox(mri_labeled, xi, y, z) ==
                Left_Cerebral_Cortex)
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   xi, y, z, Left_Hippocampus) == 0)
              {
                continue ;
              }
              if (xi == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                printf("label at (%d, %d, %d) changed "
                       "from %s to %s\n",
                       xi, y, z,
                       cma_label_to_name(Left_Cerebral_Cortex),
                       cma_label_to_name(Left_Hippocampus)) ;


              }
              MRIvox(mri_tmp, xi, y, z) = Left_Hippocampus ;
              nchanged++ ;
              total_changed++ ;
            }
            break ;
          }
          case Right_Lateral_Ventricle:
          case Right_Inf_Lat_Vent:   /* medial is in positive
                                        x direction */
          {
            int xi ;

            xi = mri_inputs->xi[x+1] ;
            if (MRIvox(mri_labeled, xi, y, z) ==
                Right_Cerebral_Cortex)
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   xi, y, z, Right_Hippocampus) == 0)
              {
                continue ;
              }
              if (xi == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                printf
                ("label at (%d, %d, %d) changed from "
                 "%s to %s\n",
                 xi, y, z,
                 cma_label_to_name(Right_Cerebral_Cortex),
                 cma_label_to_name(Right_Hippocampus)) ;
              }
              MRIvox(mri_tmp, xi, y, z) = Right_Hippocampus ;
              nchanged++ ;
              total_changed++ ;
            }
            break ;
          }
          case Left_Hippocampus:
            pwm = GCAlabelProbability
                  (mri_inputs, gca, transform, (float)x, (float)y, (float)z,
                   Left_Cerebral_White_Matter) ;
            phippo = GCAlabelProbability
                     (mri_inputs, gca, transform, (float)x, (float)y, (float)z,
                      Left_Hippocampus) ;
            if (phippo > 2*pwm)
            {
              continue ;  /* don't let it change */
            }

            dleft = distance_to_label
                    (mri_labeled, Left_Cerebral_White_Matter,
                     x,y,z,-1,0,0,3);
            dright = distance_to_label
                     (mri_labeled, Left_Cerebral_White_Matter,
                      x,y,z,1,0,0,3);
            ddown = distance_to_label
                    (mri_labeled,
                     Left_Cerebral_White_Matter,x,y,z,0,1,0,2);
            dup = distance_to_label
                  (mri_labeled,
                   Left_Cerebral_White_Matter,x,y,z,0,-1,0,2);
            dpos = distance_to_label
                   (mri_labeled,
                    Left_Cerebral_White_Matter,x,y,z,0,0,-1,3);
            dant = distance_to_label
                   (mri_labeled,
                    Left_Cerebral_White_Matter,x,y,z,0,0,1,3);

            if ((dpos <= 2 && dant <= 2) ||
                (dleft <= 2 && dright <= 2) ||
                (dup <= 1 && (dleft == 1 || dright == 1)))
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   x, y, z, Left_Cerebral_White_Matter) == 0)
              {
                continue ;
              }
              nchanged++ ;
              total_changed++ ;
              MRIvox(mri_tmp, x, y, z) =
                Left_Cerebral_White_Matter ;
            }
            break ;
          case Right_Hippocampus:
            pwm = GCAlabelProbability
                  (mri_inputs, gca, transform, (float)x, (float)y, (float)z,
                   Right_Cerebral_White_Matter) ;
            phippo = GCAlabelProbability
                     (mri_inputs, gca, transform, (float)x, (float)y, (float)z,
                      Right_Hippocampus) ;
            if (phippo > 2*pwm)
            {
              continue ;  /* don't let it change */
            }

            dleft = distance_to_label
                    (mri_labeled, Right_Cerebral_White_Matter,
                     x,y,z,-1,0,0,3);
            dright = distance_to_label
                     (mri_labeled, Right_Cerebral_White_Matter,
                      x,y,z,1,0,0,3);
            ddown = distance_to_label
                    (mri_labeled,
                     Right_Cerebral_White_Matter,x,y,z,0,1,0,2);
            dup = distance_to_label
                  (mri_labeled,
                   Right_Cerebral_White_Matter,x,y,z,0,-1,0,2);
            dpos = distance_to_label
                   (mri_labeled,
                    Right_Cerebral_White_Matter,x,y,z,0,0,-1,3);
            dant = distance_to_label
                   (mri_labeled,
                    Right_Cerebral_White_Matter,x,y,z,0,0,1,3);
            if ((dpos <= 2 && dant <= 2) ||
                (dleft <= 2 && dright <= 2) ||
                (dup <= 1 && (dleft == 1 || dright == 1)))
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   x, y, z, Right_Cerebral_White_Matter) == 0)
              {
                continue ;
              }
              nchanged++ ;
              total_changed++ ;
              MRIvox(mri_tmp, x, y, z) =
                Right_Cerebral_White_Matter ;
            }
            break ;
          default:
            break ;
          }
        } //!GCA
      }
    }
  }

  MRIcopy(mri_tmp, mri_labeled) ;
  for (i = 0 ; i < 2 ; i++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Hippocampus:
            ddown = distance_to_label
                    (mri_tmp, Left_Cerebral_White_Matter,
                     x,y,z,0,1,0,2);
            dup = distance_to_label
                  (mri_tmp,
                   Left_Cerebral_White_Matter,x,y,z,0,-1,0,3);

            ddown_gray = distance_to_label
                         (mri_tmp,Left_Cerebral_Cortex,
                          x,y,z,0,1,0,3);
            dup_hippo = distance_to_label
                        (mri_tmp,
                         Left_Hippocampus,x,y,z,0,-1,0,3);

            ddown_ventricle =
              MIN(distance_to_label(mri_tmp, Left_Lateral_Ventricle,
                                    x,y,z,0,1,0,3),
                  distance_to_label(mri_tmp, Left_Inf_Lat_Vent,
                                    x,y,z,0,1,0,3)) ;

            if (ddown_ventricle < ddown_gray)
            {
              continue ;
            }  /* gray should never be
                  superior to ventricle */

            /* closer to superior white matter
               than hp and gray below */
            if (((dup < dup_hippo) && ddown > ddown_gray))
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   x, y, z, Left_Cerebral_Cortex) == 0)
              {
                continue ;
              }
              nchanged++ ;
              total_changed++ ;
              MRIvox(mri_tmp, x, y, z) = Left_Cerebral_Cortex ;
            }
            break ;
          case Right_Hippocampus:
            ddown = distance_to_label
                    (mri_tmp,Right_Cerebral_White_Matter,
                     x,y,z,0,1,0,3);
            dup = distance_to_label
                  (mri_tmp,Right_Cerebral_White_Matter,
                   x,y,z,0,-1,0,3);

            ddown_gray = distance_to_label
                         (mri_tmp,
                          Right_Cerebral_Cortex,
                          x,y,z,0,1,0,3);
            dup_hippo = distance_to_label
                        (mri_tmp,
                         Right_Hippocampus,x,y,z,0,-1,0,3);

            ddown_ventricle =
              MIN(distance_to_label
                  (mri_tmp, Right_Lateral_Ventricle,
                   x,y,z,0,1,0,3),
                  distance_to_label
                  (mri_tmp, Right_Inf_Lat_Vent,
                   x,y,z,0,1,0,3)) ;

            if (ddown_ventricle < ddown_gray)
            {
              continue ;
            }  /* gray should never be
                  superior to ventricle */

            /* closer to superior white matter
               than hp and gray below */
            if (((dup < dup_hippo) && ddown > ddown_gray))
            {
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   x, y, z, Right_Cerebral_Cortex) == 0)
              {
                continue ;
              }
              nchanged++ ;
              total_changed++ ;
              MRIvox(mri_tmp, x, y, z) = Right_Cerebral_Cortex ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_labeled) ;
  }

  /* change gray to hippocampus based on wm */
  for (i = 0 ; i < 3 ; i++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            left = 1 ;
          case Right_Cerebral_Cortex:
            /*
              if the current label is gray,
              and there is white matter below
              and hippocampus above, change to hippocampus.
            */
            ddown = distance_to_label
                    (mri_labeled,
                     left ? Left_Cerebral_White_Matter :
                     Right_Cerebral_White_Matter,
                     x,y,z,0,1,0,3);
            dup = distance_to_label
                  (mri_labeled,
                   left ?  Left_Hippocampus :
                   Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 3)
            {
              label = left ? Left_Hippocampus : Right_Hippocampus ;
              if (GCAlabelExists
                  (gca, mri_labeled, transform,
                   x, y, z, label) == 0)
              {
                continue ;
              }

              nchanged++ ;
              total_changed++ ;
              MRIvox(mri_tmp, x, y, z) = label ;
            }
            break ;
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_labeled) ;
  }

  /* change hippocampal voxels that have a run of wm superior to them
     followed by a run of hippo.
  */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIvox(mri_tmp, x, y, z) ;

        left = 0 ;
        was_unknown = 0 ;
        switch (label)
        {
        case Left_Hippocampus:
          left = 1 ;
        case Right_Hippocampus:
#define MIN_UNKNOWN 5
          for (i = 1 ; i <= MIN_UNKNOWN+3 ; i++)
          {
            yi = mri_tmp->yi[y+i] ;
            olabel = MRIvox(mri_tmp, x, yi, z) ;
            if (olabel == Left_Hippocampus ||
                olabel == Right_Hippocampus)
            {
              was_unknown = 0 ;
              continue ;
            }
            if (!IS_UNKNOWN(olabel))  /* don't change it */
            {
              break ;
            }
            if (++was_unknown >= MIN_UNKNOWN)
            {
              break ;
            }
          }
          if (was_unknown >= MIN_UNKNOWN)
          {
            if (GCAlabelExists
                (gca, mri_labeled, transform,
                 x, y, z, Left_Cerebral_Cortex) == 0 &&
                GCAlabelExists
                (gca, mri_labeled, transform,
                 x, y, z, Right_Cerebral_Cortex) == 0)
            {
              continue ;
            }
            nchanged++ ;
            total_changed++ ;
            MRIvox(mri_tmp, x, y, z) =
              left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf
              ("(%d, %d, %d) %s changed to %s "
               "in edit_hippocampus\n",
               x, y, z, cma_label_to_name(label),
               cma_label_to_name(MRIvox(mri_tmp, x, y, z))) ;
            }
          }

          break ;
        default:
          break ;
        }
      }
    }
  }

  /* if voxel is hippocampus with only hippocampus and a bunch of unknowns
     inferior to it, change it to cortical gray.
  */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIvox(mri_tmp, x, y, z) ;

        left = 0 ;
        found_wm = found_hippo = was_wm = was_hippo = 0 ;
        switch (label)
        {
        case Left_Hippocampus:
          left = 1 ;
        case Right_Hippocampus:
          for (i = 1 ; i <= 10 ; i++)
          {
            yi = mri_tmp->yi[y-i] ;
            olabel = MRIvox(mri_tmp, x, yi, z) ;
            if (found_wm)  /* check for hippo */
            {
              if (olabel == Left_Hippocampus ||
                  olabel == Right_Hippocampus)
              {
                was_hippo++ ;
                if (was_hippo >= 2)
                {
                  found_hippo = 1 ;
                }
              }
              else
              {
                was_hippo = 0 ;
              }
            }
            else if (olabel == Right_Cerebral_White_Matter ||
                     olabel == Left_Cerebral_White_Matter)
            {
              was_wm++ ;
              if (was_wm >= 3)
              {
                found_wm = 1 ;
              }
            }
            else
            {
              was_wm = 0 ;
            }
          }
          if (found_wm && found_hippo)
          {
            if (GCAlabelExists
                (gca, mri_labeled, transform,
                 x, y, z, Left_Cerebral_Cortex) == 0 &&
                GCAlabelExists
                (gca, mri_labeled, transform,
                 x, y, z, Right_Cerebral_Cortex) == 0)
            {
              continue ;
            }
            nchanged++ ;
            total_changed++ ;
            MRIvox(mri_tmp, x, y, z) =
              left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf("(%d, %d, %d) %s changed to "
                     "%s in edit_hippocampus\n",
                     x, y, z, cma_label_to_name(label),
                     cma_label_to_name(MRIvox(mri_tmp, x, y, z))) ;
            }
          }

          break ;
        default:
          break ;
        }
      }
    }
  }

  /* do a close on the hippocampus (only on cortical gray matter voxels) */
  do
  {
    for (nchanged = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (GCAlabelExists
              (gca, mri_labeled, transform,
               x, y, z, Left_Hippocampus) == 0 &&
              GCAlabelExists
              (gca, mri_labeled, transform,
               x, y, z, Right_Hippocampus) == 0)
          {
            continue ;
          }

          label = MRIvox(mri_labeled, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            if (MRIcountNbhdLabels
                (mri_labeled, x, y,z,Left_Hippocampus) >= 20)
            {
              nchanged++ ;
              total_changed++ ;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                printf("label at (%d, %d, %d) changed from "
                       "%s to %s\n",
                       x, y, z,
                       cma_label_to_name(Left_Cerebral_Cortex),
                       cma_label_to_name(Left_Hippocampus)) ;
              }
              MRIvox(mri_tmp, x, y, z) = Left_Hippocampus ;
            }
            break ;
          case Right_Cerebral_Cortex:
            if (MRIcountNbhdLabels
                (mri_labeled, x,y,z,Right_Hippocampus) >= 20)
            {
              nchanged++ ;
              total_changed++ ;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                printf("label at (%d, %d, %d) changed from "
                       "%s to %s\n",
                       x, y, z,
                       cma_label_to_name(Right_Cerebral_Cortex),
                       cma_label_to_name(Right_Hippocampus)) ;
              }
              MRIvox(mri_tmp, x, y, z) = Right_Hippocampus ;
            }
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_labeled) ;
  }
  while (nchanged > 0) ;

  /* check for long  runs of hippocampus next to gm */
  do
  {
    int xi ;
    for (nchanged = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }

          label = MRIvox(mri_labeled, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Hippocampus:
            xi = mri_labeled->xi[x+1] ;
            if (MRIvox(mri_labeled, xi, y, z) ==
                Left_Cerebral_Cortex)
            {
              found_hippo = 1 ;
              for (i = 0 ; found_hippo && i < 8 ; i++)
              {
                xi = mri_labeled->xi[x-i] ;
                if (MRIvox(mri_labeled, xi, y, z) !=
                    Left_Hippocampus)
                {
                  found_hippo = 0 ;
                }
              }
              if (found_hippo)
              {
                xi = mri_labeled->xi[x+1] ;
                if (GCAlabelExists
                    (gca, mri_labeled, transform,
                     xi, y, z, Left_Hippocampus) == 0)
                {
                  continue ;
                }
                nchanged++ ;
                total_changed++ ;
                if (xi == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  printf
                  ("label at (%d, %d, %d) changed from "
                   "%s to %s\n",
                   xi, y, z,
                   cma_label_to_name(Left_Cerebral_Cortex),
                   cma_label_to_name(Left_Hippocampus)) ;
                }
                MRIvox(mri_labeled, xi, y, z) =
                  Left_Hippocampus ;
              }
            }
            break ;
          case Right_Hippocampus:
            xi = mri_labeled->xi[x+1] ;
            if (MRIvox(mri_tmp, xi, y, z) == Right_Cerebral_Cortex)
            {
              found_hippo = 1 ;
              for (i = 0 ; found_hippo && i < 8 ; i++)
              {
                xi = mri_labeled->xi[x+i] ;
                if (MRIvox(mri_labeled, xi, y, z) !=
                    Right_Hippocampus)
                {
                  found_hippo = 0 ;
                }
              }
              if (found_hippo)
              {
                xi = mri_labeled->xi[x+1] ;
                if (GCAlabelExists
                    (gca, mri_labeled, transform,
                     xi, y, z, Right_Hippocampus) == 0)
                {
                  continue ;
                }
                nchanged++ ;
                total_changed++ ;
                if (xi == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  printf
                  ("label at (%d, %d, %d) changed from "
                   "%s to %s\n",
                   xi, y, z,
                   cma_label_to_name(Right_Cerebral_Cortex),
                   cma_label_to_name(Right_Hippocampus)) ;
                }
                MRIvox(mri_labeled, xi, y, z) =
                  Right_Hippocampus ;
              }
            }
            break ;
          }
        }
      }
    }
  }
  while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d hippocampal voxels changed.\n", total_changed) ;
  return(NO_ERROR) ;
}


int
edit_amygdala( MRI *mri_inputs,
               MRI *mri_labeled,
               GCA *gca,
               TRANSFORM *transform,
               MRI *mri_fixed )
{
  int   width, height, depth, x, y, z, label, nchanged, olabel,
        i, left, yi, found_wm, found_amygdala, was_wm, was_amygdala ;
  MRI   *mri_tmp ;

  mri_tmp = MRIcopy(mri_labeled, NULL) ;

  width = mri_inputs->width ;
  height = mri_inputs->height ;
  depth = mri_inputs->depth ;

  nchanged = 0 ;
  /* change amygdala voxels that have a run of wm superior to them
     followed by a run of amygdala.
  */
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIvox(mri_tmp, x, y, z) ;

        left = 0 ;
        found_wm = found_amygdala = was_wm = was_amygdala = 0 ;
        switch (label)
        {
        case Left_Amygdala:
          left = 1 ;
        case Right_Amygdala:
          for (i = 1 ; i <= 10 ; i++)
          {
            yi = mri_tmp->yi[y-i] ;
            olabel = MRIvox(mri_tmp, x, yi, z) ;
            if (found_wm)  /* check for amygdala */
            {
              if (olabel == Left_Amygdala ||
                  olabel == Right_Amygdala)
              {
                was_amygdala++ ;
                if (was_amygdala >= 2)
                {
                  found_amygdala = 1 ;
                }
              }
              else
              {
                was_amygdala = 0 ;
              }
            }
            else if (olabel == Right_Cerebral_White_Matter ||
                     olabel == Left_Cerebral_White_Matter)
            {
              was_wm++ ;
              if (was_wm >= 3)
              {
                found_wm = 1 ;
              }
            }
            else
            {
              was_wm = 0 ;
            }
          }
          label = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;

          if (found_wm &&
              found_amygdala &&
              GCAisPossible
              (gca, mri_tmp, label, transform, x, y, z,0))
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = label ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf("(%d, %d, %d) %s changed to "
                     "%s in edit_amygdala\n",
                     x, y, z, cma_label_to_name(label),
                     cma_label_to_name(MRIvox(mri_tmp, x, y, z))) ;
            }
          }

          break ;
        default:
          break ;
        }
      }
    }
  }

  MRIcopy(mri_tmp, mri_labeled) ;
  MRIfree(&mri_tmp) ;
  printf("%d amygdala voxels changed.\n", nchanged) ;
  return(NO_ERROR) ;
}

int
insert_thin_temporal_white_matter( MRI *mri_inputs, MRI *mri_labeled,
                                   GCA *gca, TRANSFORM *transform,
                                   GCA *gca_all )
{
  MRI       *mri_tmp, *mri_probs, *mri_tmp_labels ;
  int       width, height, depth, x, y, z, n, nsamples, i, xp, yp, zp ;
  int       xmin, xmax, ymin, ymax, zmin, zmax,yi,zi, yimax, ximax, zimax,
            **added, nchanged, label, r, c, v ;
  GCA_PRIOR *gcap ;
  GC1D      *gc ;
  GCA_SAMPLE *gcas ;
  double     p, pmax ;

  mri_tmp = MRIclone(mri_inputs, NULL) ;
  mri_tmp_labels = MRIclone(mri_labeled, NULL) ;

  width = mri_tmp->width ;
  height = mri_tmp->height ;
  depth = mri_tmp->depth ;
  mri_probs = MRIalloc(width, height, depth, MRI_UCHAR) ;

  xmin = width ;
  ymin = height ;
  zmin = depth ;
  xmax = ymax = zmax = 0 ;
  for (nsamples = x = 0 ; x < width ; x++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        if (!GCAsourceVoxelToPrior
            (gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp))
        {
          gcap = &gca->priors[xp][yp][zp] ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            if (!IS_UNKNOWN(gcap->labels[n]))
            {
              if (x < xmin)
              {
                xmin = x ;
              }
              if (x > xmax)
              {
                xmax = x ;
              }
              if (y < ymin)
              {
                ymin = y ;
              }
              if (y > ymax)
              {
                ymax = y ;
              }
              if (z < zmin)
              {
                zmin = z ;
              }
              if (z > zmax)
              {
                zmax = z ;
              }

              nsamples++ ;
              break ;
            }
          }
        }// !GCA
      }
    }
  }

  printf("allocating %d TL samples, box [%d, %d, %d] -> [%d, %d, %d]...\n",
         nsamples, xmin, ymin, zmin, xmax, ymax, zmax) ;
  gcas = calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  if (!gcas)
  {
    ErrorExit(ERROR_NOMEMORY, "could not allocate gcas for TL insertion") ;
  }

  /* put in the coordinates of the samples */
  for (i = 0, x = xmin ; x <= xmax ; x++)
  {
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        if (!GCAsourceVoxelToPrior
            (gca, mri_inputs, transform, x, y, z, &xp, &yp, &zp))
        {
          gcap = &gca->priors[xp][yp][zp] ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            if (!IS_UNKNOWN(gcap->labels[n]))
            {
              gc =
                GCAfindPriorGC(gca, xp, yp, zp, gcap->labels[n]) ;
              gcas[i].x = x ;
              gcas[i].y = y ;
              gcas[i].z = z ;
              gcas[i].xp = xp ;
              gcas[i].yp = yp ;
              gcas[i].zp = zp ;
              gcas[i].label = gcap->labels[n] ;
              gcas[i].prior = getPrior(gcap, gcap->labels[n]) ;
              for (r = v = 0 ; r < gca->ninputs ; r++)
              {
                gcas[i].means[r] = gc->means[r] ;
                for (c = r ; c < gca->ninputs ; c++, v++)
                {
                  gcas[i].covars[v] = gc->covars[v] ;
                }
              }
              i++ ;
              break ;
            }
          }
        }// !GCA
      }
    }
  }
  GCAcomputeLogSampleProbabilityUsingCoords(gca, gcas, mri_inputs,
      transform, nsamples, DEFAULT_CLAMP) ;
  for (i = 0 ; i < nsamples ; i++)
  {
    double p ;
    p = exp(gcas[i].log_p) * 255 * 20 ;
    if (p > 255.0)
    {
      p = 255.0 ;
    }
    MRIvox(mri_probs, gcas[i].x, gcas[i].y, gcas[i].z) = (char)p ;
    if (gcas[i].x == Ggca_x && gcas[i].y == Ggca_y && gcas[i].z == Ggca_z)
    {
      DiagBreak() ;
    }
    MRIvox(mri_tmp_labels, gcas[i].x, gcas[i].y, gcas[i].z) = gcas[i].label ;
  }

  added = (int **)calloc(width, sizeof(int *)) ;
  if (!added)
  {
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate added array", Progname);
  }
  for (x = 0 ; x < width ; x++)
  {
    added[x] = (int *)calloc(depth, sizeof(int)) ;
    if (!added[x])
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate added array[%d]",
                Progname, x);
  }

  /******** Left Hemisphere ***********/

  /* find highest probability pair over whole temporal lobe,
     and use them to bootstrap process.
  */
  pmax = 0.0 ;
  ximax = yimax = zimax = 0 ;
  for (i = 0, x = xmin ; x <= xmax ; x++)
  {
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        if (MRIvox(mri_tmp_labels,x,y,z) != Left_Cerebral_White_Matter)
        {
          continue ;
        }
        yi = mri_probs->yi[y+1] ;
        p = MRIvox(mri_probs, x,  y, z) + MRIvox(mri_probs, x, yi, z) ;
        if (p > pmax)
        {
          pmax = p ;
          ximax = x ;
          yimax = y ;
          zimax = z ;
        }
      }
    }
  }


  printf("using (%d, %d, %d) as seed point for LH\n", ximax, yimax, zimax) ;

  /* put the seed point back in */
  MRIvox(mri_tmp, ximax, yimax, zimax) = 255 ;
  MRIvox(mri_tmp, ximax, yimax+1, zimax) = 255 ;
  xmin = MAX(0,ximax-1) ;
  xmax = MIN(ximax+1,width-1) ;
  ymin = MAX(0,yimax-1) ;
  ymax = MIN(yimax+2, height-1) ;
  zmin = MAX(0,zimax-1) ;
  zmax = MIN(zimax+1, depth-1) ;
  added[ximax][zimax] = 1 ;  /* already added to this plane */
  do
  {
    nchanged = 0 ;

    for (x = xmin ; x <= xmax ; x++)
    {
      for (z = zmin ; z <= zmax ; z++)
      {
        if (added[x][z])
        {
          continue ;
        }
        pmax = 0.0 ;
        yimax = -1 ;
        for (y = ymin ; y <= ymax ; y++)
        {
          if (MRIvox(mri_tmp, x, y, z))
          {
            continue ;
          }
#if 0
          /* enforce connectivity of temporal wm */
          if (!MRIvox(mri_tmp, x+1, y, z)  &&
              !MRIvox(mri_tmp, x-1, y, z) &&
              (!MRIvox(mri_tmp, x, y, z+1) &&
               !MRIvox(mri_tmp, x, y, z-1)))
          {
            continue ;
          }
#endif
          yi = mri_probs->yi[y+1] ;
          p =
            MRIvox(mri_probs, x,  y, z) +
            MRIvox(mri_probs, x, yi, z) ;
          if (p > pmax)
          {
            pmax = p ;
            yimax = y ;
          }
        }

        if (yimax < 0)
        {
          continue ;  /* couldn't find one */
        }

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        if (x == Ggca_x && z == Ggca_z)
        {
          DiagBreak() ;
        }
        added[x][z] = 1 ;
        if (x == xmin)
        {
          xmin-- ;
        }
        if (x == xmax)
        {
          xmax++ ;
        }
        if (yimax == ymin)
        {
          ymin-- ;
        }
        if (yimax == ymax)
        {
          ymax++ ;
        }
        if (z == zmin)
        {
          zmin-- ;
        }
        if (z == zmax)
        {
          zmax++ ;
        }

        MRIvox(mri_tmp, x, yimax, z) = 255 ;
        MRIvox(mri_tmp, x, yimax+1, z) = 255 ;
        nchanged++ ;
      }
    }
  }
  while (nchanged > 0) ;


  /******** Right Hemisphere ***********/

  /* recompute the bounding box */
  xmin = width ;
  ymin = height ;
  zmin = depth ;
  xmax = ymax = zmax = 0 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    x = gcas[i].x ;
    y = gcas[i].y ;
    z = gcas[i].z ;
    if (x < xmin)
    {
      xmin = x ;
    }
    if (x > xmax)
    {
      xmax = x ;
    }
    if (y < ymin)
    {
      ymin = y ;
    }
    if (y > ymax)
    {
      ymax = y ;
    }
    if (z < zmin)
    {
      zmin = z ;
    }
    if (z > zmax)
    {
      zmax = z ;
    }
  }

  /* find highest probability pair over whole temporal lobe,
     and use them to bootstrap process.
  */
  pmax = 0.0 ;
  ximax = yimax = zimax = 0 ;
  for (i = 0, x = xmin ; x <= xmax ; x++)
  {
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        if (MRIvox(mri_tmp_labels,x,y,z) != Right_Cerebral_White_Matter)
        {
          continue ;
        }
        yi = mri_probs->yi[y+1] ;
        p = MRIvox(mri_probs, x,  y, z) + MRIvox(mri_probs, x, yi, z) ;
        if (p > pmax)
        {
          pmax = p ;
          ximax = x ;
          yimax = y ;
          zimax = z ;
        }
      }
    }
  }


  printf("using (%d, %d, %d) as seed point for RH\n", ximax, yimax, zimax) ;

  /* put the seed point back in */
  MRIvox(mri_tmp, ximax, yimax, zimax) = 255 ;
  MRIvox(mri_tmp, ximax, yimax+1, zimax) = 255 ;
  xmin = MAX(0,ximax-1) ;
  xmax = MIN(ximax+1,width-1) ;
  ymin = MAX(0,yimax-1) ;
  ymax = MIN(yimax+2, height-1) ;
  zmin = MAX(0,zimax-1) ;
  zmax = MIN(zimax+1, depth-1) ;
  added[ximax][zimax] = 1 ;  /* already added to this plane */
  do
  {
    nchanged = 0 ;

    for (x = xmin ; x <= xmax ; x++)
    {
      for (z = zmin ; z <= zmax ; z++)
      {
        if (added[x][z])
        {
          continue ;
        }
        pmax = 0.0 ;
        yimax = -1 ;
        for (y = ymin ; y <= ymax ; y++)
        {
          if (MRIvox(mri_tmp, x, y, z))
          {
            continue ;
          }
#if 0
          /* enforce connectivity of temporal wm */
          if (!MRIvox(mri_tmp, x+1, y, z)  &&
              !MRIvox(mri_tmp, x-1, y, z) &&
              (!MRIvox(mri_tmp, x, y, z+1) &&
               !MRIvox(mri_tmp, x, y, z-1)))
          {
            continue ;
          }
#endif
          yi = mri_probs->yi[y+1] ;
          p =
            MRIvox(mri_probs, x,  y, z) +
            MRIvox(mri_probs, x, yi, z) ;
          if (p > pmax)
          {
            pmax = p ;
            yimax = y ;
          }
        }

        if (yimax < 0)
        {
          continue ;  /* couldn't find one */
        }

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        if (x == Ggca_x && z == Ggca_z)
        {
          DiagBreak() ;
        }
        added[x][z] = 1 ;
        if (x == xmin)
        {
          xmin-- ;
        }
        if (x == xmax)
        {
          xmax++ ;
        }
        if (yimax == ymin)
        {
          ymin-- ;
        }
        if (yimax == ymax)
        {
          ymax++ ;
        }
        if (z == zmin)
        {
          zmin-- ;
        }
        if (z == zmax)
        {
          zmax++ ;
        }

        MRIvox(mri_tmp, x, yimax, z) = 255 ;
        MRIvox(mri_tmp, x, yimax+1, z) = 255 ;
        nchanged++ ;
      }
    }
  }
  while (nchanged > 0) ;

  /* now do some spackling and flossing */
  /* recompute the bounding box */
  xmin = width ;
  ymin = height ;
  zmin = depth ;
  xmax = ymax = zmax = 0 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    x = gcas[i].x ;
    y = gcas[i].y ;
    z = gcas[i].z ;
    if (x < xmin)
    {
      xmin = x ;
    }
    if (x > xmax)
    {
      xmax = x ;
    }
    if (y < ymin)
    {
      ymin = y ;
    }
    if (y > ymax)
    {
      ymax = y ;
    }
    if (z < zmin)
    {
      zmin = z ;
    }
    if (z > zmax)
    {
      zmax = z ;
    }
  }
  do
  {
    nchanged = 0 ;

    for (x = xmin ; x <= xmax ; x++)
    {
      for (z = zmin ; z <= zmax ; z++)
      {
#if 0
        if (added[x][z])
        {
          continue ;
        }
#endif
        for (y = ymin ; y <= ymax ; y++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (MRIvox(mri_tmp, x, y, z))
          {
            continue ;
          }
          if (z == 0 || z>= depth)
          {
            DiagBreak() ;
          }
          if (x <= 0 || x >= width || z <= 0 || z >= depth)
          {
            continue ;
          }
          if ((MRIvox(mri_tmp, x+1, y, z) &&
               MRIvox(mri_tmp, x-1,y,z)) ||
              (MRIvox(mri_tmp, x, y, z-1) &&
               MRIvox(mri_tmp, x,y,z+1)))
          {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              DiagBreak() ;
            }

            nchanged++ ;
            MRIvox(mri_tmp, x,y,z) =
              MAX(MAX(MAX(MRIvox(mri_tmp, x+1, y, z),
                          MRIvox(mri_tmp, x-1,y,z)),
                      MRIvox(mri_tmp, x, y, z-1)),
                  MRIvox(mri_tmp, x,y,z+1)) ;
          }
        }
      }
    }
  }
  while (nchanged > 0) ;

  do
  {
    nchanged = 0 ;
    for (x = xmin ; x <= xmax ; x++)
    {
      for (z = zmin ; z <= zmax ; z++)
      {
#if 0
        if (added[x][z])
        {
          continue ;
        }
#endif
        for (y = ymin ; y <= ymax ; y++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (MRIvox(mri_tmp, x, y, z))
          {
            continue ;
          }
          yi = mri_tmp->yi[y+1] ;
          if (MRIvox(mri_tmp,x,yi,z)) /* check inferior */
          {
            yi = mri_tmp->yi[y-1] ;
            zi = mri_tmp->zi[z+1] ;
            if (MRIvox(mri_tmp, x, yi, zi)) /* inferior voxel on */
            {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                DiagBreak() ;
              }
              nchanged++ ;
              MRIvox(mri_tmp, x,y,z) = 255 ;
            }
            else
            {
              zi = mri_tmp->zi[z-1] ;
              if (MRIvox(mri_tmp, x, yi,  zi))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  DiagBreak() ;
                }
                nchanged++ ;
                MRIvox(mri_tmp, x,y,z) = 255 ;
              }
            }
          }

          yi = mri_tmp->yi[y-1] ;
          if (MRIvox(mri_tmp,x,yi,z)) /* check suprior */
          {
            yi = mri_tmp->yi[y+1] ;
            zi = mri_tmp->zi[z+1] ;
            if (MRIvox(mri_tmp, x, yi, zi)) /* inferior voxel on */
            {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                DiagBreak() ;
              }
              nchanged++ ;
              MRIvox(mri_tmp, x,y,z) = 255 ;
            }
            else
            {
              zi = mri_tmp->zi[z-1] ;
              if (MRIvox(mri_tmp, x, yi,  zi))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  DiagBreak() ;
                }
                nchanged++ ;
                MRIvox(mri_tmp, x,y,z) = 255 ;
              }
            }
          }
        }
      }
    }
  }
  while (nchanged > 0) ;


  /* do a couple of iterations adding maxiumum likelihood wm voxels
     mark these separately as it will grow into the main body of the
     wm and hence not bound the hippocampus.
  */
  for (i = 0 ; i < 3 ; i++)
  {
    int labels[MAX_CMA_LABEL+1], nlabels, max_label, n ;
    double probs[MAX_CMA_LABEL+1], max_p ;

    for (x = xmin ; x <= xmax ; x++)
    {
      for (z = zmin ; z <= zmax ; z++)
      {
#if 0
        if (added[x][z])
        {
          continue ;
        }
#endif
        for (y = ymin ; y <= ymax ; y++)
        {
          if (MRIvox(mri_tmp, x, y, z))
          {
            continue ;
          }
          if (MRIneighborsOn(mri_tmp, x, y, z, 1) == 0)
          {
            continue ;
          }
          nlabels =
            GCAcomputeVoxelLikelihoods
            (gca_all,mri_inputs,x,y,z,transform,labels,probs);
          for (max_p = 0, max_label = -1, n = 0 ; n < nlabels ; n++)
          {
            if (probs[n] > max_p)
            {
              max_p = probs[n] ;
              max_label = labels[n] ;
            }
          }
          if (IS_WM(max_label))
          {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              DiagBreak() ;
            }
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = 128 ;
          }
        }
      }
    }
  }

  /*
     change hippocampal labels inferior to TL to GM, and GM superior
     to TL to hippocampus.
  */
  for (x = xmin ; x <= xmax ; x++)
  {
    for (z = zmin ; z <= zmax ; z++)
    {
#if 0
      if (added[x][z])
      {
        continue ;
      }
#endif
      for (y = ymin-1 ; y <= ymax+1 ; y++)
      {
        if (y < 0 || y >= height)
        {
          continue ;
        }

        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }

        if (MRIvox(mri_tmp, x, y, z) != 255) /* not temporal wm */
        {
          continue ;
        }

        /* examine label inferior to this one.
           If it's hippo or ventricle,
           change it to cortex.
        */
        for (i = 1 ; i < 10 ; i++)
        {
          /* check for hippocampus inferior */
          yi = mri_tmp->yi[y-i] ;
          label = MRIvox(mri_labeled,x,yi,z) ;
#if 0
          if (label == Left_Cerebral_Cortex ||
              label == Right_Cerebral_Cortex)
          {
            int left ;

            left = (label == Left_Cerebral_Cortex) ;
            label = left ? Left_Hippocampus : Right_Hippocampus ;
            if (x == Ggca_x && yi == Ggca_y && z == Ggca_z)
            {
              printf
              ("(%d, %d, %d) %s changed to %s in insert tl\n",
               x, yi, z,
               cma_label_to_name(MRIvox(mri_labeled,x,yi, z)),
               cma_label_to_name(label)) ;
            }
            MRIvox(mri_labeled, x, yi, z) = label ;
          }
#else
          if (!GCAsourceVoxelToPrior
              (gca_all, mri_inputs, transform,
               x, yi, z, &xp, &yp, &zp))
          {
            gcap = &gca_all->priors[xp][yp][zp] ;
            if (label == Left_Cerebral_Cortex ||
                label == Right_Cerebral_Cortex)
            {
              int left ;

              left = (label == Left_Cerebral_Cortex) ;
              label = left ? Left_Hippocampus : Right_Hippocampus ;
              for (n = 0 ; n < gcap->nlabels ; n++)
              {
                if (gcap->labels[n] ==
                    label)  /* it's possible */
                {
                  if (x == Ggca_x &&
                      yi == Ggca_y &&
                      z == Ggca_z)
                  {
                    printf
                    ("(%d, %d, %d) %s changed to %s "
                     "in insert tl\n",
                     x, yi, z,
                     cma_label_to_name
                     (MRIvox(mri_labeled,x,yi, z)),
                     cma_label_to_name(label)) ;
                  }
                  MRIvox(mri_labeled, x, yi, z) = label ;
                }
              }
            }
          }//!GCA
#endif
          /* check for hippocampus inferior */
          yi = mri_tmp->yi[y+i] ;
          label = MRIvox(mri_labeled,x,yi,z) ;
          if (IS_HIPPO(label) || IS_LAT_VENT(label))
          {
            int left ;

            left =
              (label == Left_Hippocampus) ||
              (label == Left_Lateral_Ventricle) ||
              (label == Left_Inf_Lat_Vent) ;
            label =
              left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
            if (is_possible
                (gca, mri_labeled, transform, x, yi, z, label))
            {
              if (x == Ggca_x && yi == Ggca_y && z == Ggca_z)
              {
                printf
                ("(%d, %d, %d) %s changed to %s "
                 "in insert tl\n",
                 x, yi, z,
                 cma_label_to_name
                 (MRIvox(mri_labeled,x,yi, z)),
                 cma_label_to_name(label)) ;
              }
              MRIvox(mri_labeled, x, yi, z) = label ;
            }
          }
        }
      }
    }
  }

#if 0
  MRIclose(mri_tmp, mri_tmp) ;
  MRIdilate6(mri_tmp, mri_tmp) ;
#endif

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    printf("writing temporal lobe volume...") ;
    MRIwrite(mri_tmp, "temp.mgz") ;
    MRIwrite(mri_probs, "probs.mgz") ;
  }

  for (i = 0 ; i < nsamples ; i++)
  {
    x = gcas[i].x ;
    y = gcas[i].y ;
    z = gcas[i].z ;
    if (MRIvox(mri_tmp, x, y, z))
    {
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
      {
        printf("changing voxel (%d, %d, %d) from %s to %s\n",
               x, y, z, cma_label_to_name(MRIvox(mri_labeled,x, y,z)),
               cma_label_to_name(gcas[i].label)) ;
        DiagBreak() ;
      }
      MRIvox(mri_labeled, x, y, z) = gcas[i].label ;
    }
  }

  free(gcas) ;

  MRIfree(&mri_tmp) ;
  MRIfree(&mri_probs) ;
  MRIfree(&mri_tmp_labels) ;
  for (x = 0 ; x < width ; x++)
  {
    free(added[x]) ;
  }
  free(added) ;

  return(NO_ERROR) ;
}

int
is_possible( GCA *gca,
             MRI *mri,
             TRANSFORM *transform,
             int x, int y, int z, int label )
{
  GCA_PRIOR  *gcap ;
  int        n ;

  gcap = getGCAP(gca, mri, transform, x, y, z) ;
  for (n = 0 ; n < gcap->nlabels ; n++)
    if (gcap->labels[n] == label)
    {
      return(1) ;
    }
  return(0) ;
}

#if 0
#ifdef WSIZE
#undef WSIZE
#endif
#ifdef WHALF
#undef WHALF
#endif

#define WSIZE  5
#define WHALF  ((WSIZE-1)/2)

MRI *
GCAlabelWMandWMSAs( GCA *gca,
                    MRI *mri_inputs,
                    MRI *mri_src_labels,
                    MRI *mri_dst_labels,
                    TRANSFORM *transform )
{
  int    h, wm_label, wmsa_label, x, y, z, label, nwm, nwmsa, nunknown, ngm,
         ncaudate, caudate_label, gm_label, n, found, i;
  MATRIX *m_cov_wm, *m_cov_wmsa, *m_inv_cov_wmsa, *m_inv_cov_wm, *m_I,
         *m_cov_un, *m_inv_cov_un;
  VECTOR *v_mean_wm, *v_mean_wmsa, *v_vals, *v_dif_label, *v_dif_wmsa,
         *v_mean_caudate, *v_mean_un ;
  double pwm, pwmsa, wmsa_dist, wm_dist, wm_mdist, wmsa_mdist ;
  GCA_PRIOR *gcap ;
  MRI       *mri_tmp = NULL ;

  mri_dst_labels = MRIcopy(mri_src_labels, mri_dst_labels) ;

  v_vals = VectorAlloc(mri_inputs->nframes, MATRIX_REAL) ;
  v_dif_label = VectorAlloc(mri_inputs->nframes, MATRIX_REAL) ;
  v_dif_wmsa = VectorAlloc(mri_inputs->nframes, MATRIX_REAL) ;
  m_I = MatrixIdentity(mri_inputs->nframes, NULL) ;
  for (h = 0 ; h <= 1 ; h++)
  {
    if (h == 0) // lh
    {
      wm_label = Left_Cerebral_White_Matter ;
      wmsa_label = Left_WM_hypointensities ;
      caudate_label = Left_Caudate ;
      gm_label = Left_Cerebral_Cortex ;
    }
    else
    {
      wm_label = Right_Cerebral_White_Matter ;
      wmsa_label = Right_WM_hypointensities ;
      caudate_label = Right_Caudate ;
      gm_label = Right_Cerebral_Cortex ;
    }

    GCAcomputeLabelMeansAndCovariances(gca, Unknown, &m_cov_un, &v_mean_un) ;
    GCAcomputeLabelMeansAndCovariances(gca, wm_label, &m_cov_wm, &v_mean_wm) ;
    GCAcomputeLabelMeansAndCovariances(gca, caudate_label, &m_cov_wm, &v_mean_caudate) ;
    GCAcomputeLabelMeansAndCovariances(gca, wmsa_label, &m_cov_wmsa, &v_mean_wmsa) ;
    m_inv_cov_wm = MatrixInverse(m_cov_wm, NULL) ;
    if (m_inv_cov_wm == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not compute inverse covariance for %s (%d)",
                Progname, cma_label_to_name(wm_label), wm_label) ;
    m_inv_cov_un = MatrixInverse(m_cov_un, NULL) ;
    if (m_inv_cov_un == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not compute inverse covariance for %s (%d)",
                Progname, cma_label_to_name(Unknown), Unknown) ;
    m_inv_cov_wmsa = MatrixInverse(m_cov_wmsa, NULL) ;
    if (m_inv_cov_wmsa == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not compute inverse covariance for %s (%d)",
                Progname, cma_label_to_name(wmsa_label), wmsa_label) ;

    // do max likelihood reclassification of possible wmsa voxels
    // if they are in a nbhd with likely labels
    for (x = 0 ; x < mri_inputs->width ; x++)
    {
      for (y = 0 ; y < mri_inputs->height ; y++)
      {
        for (z = 0 ; z < mri_inputs->depth ; z++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          label = MRIgetVoxVal(mri_src_labels, x, y, z, 0) ;
          if (label != wm_label && label != wmsa_label && label != Unknown)
          {
            continue ;
          }
          // only process it if it's in the body of the wm
          nwm = MRIlabelsInNbhd(mri_src_labels, x, y, z, WHALF, wm_label) ;
          nwmsa = MRIlabelsInNbhd(mri_src_labels, x, y, z,WHALF, wmsa_label) ;

          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("(%d, %d, %d) - %s (nbrs = %d + %d = %2.2f%%)\n",
                   x, y, z, cma_label_to_name(label),
                   nwm, nwmsa, (double)(nwm+nwmsa)*100.0/(WSIZE*WSIZE*WSIZE));
          if (label == Unknown)
          {
            // only unknowns that are close to wm
            if (nwm+nwmsa < 0.5*WSIZE*WSIZE*WSIZE)
            {
              continue ;
            }
            nunknown = MRIlabelsInNbhd(mri_src_labels, x, y, z,WHALF, Unknown) ;
            if (nwm+nwmsa+nunknown < 0.9*WSIZE*WSIZE*WSIZE)
            {
              continue ;
            }
          }
          else if (nwm+nwmsa < .9*WSIZE*WSIZE*WSIZE)   // somewhat arbitrary - the bulk of the nbhd
          {
            continue ;
          }

          gcap = getGCAP(gca, mri_dst_labels, transform, x, y, z) ;
          for (found = n = 0 ; n < gcap->nlabels ; n++)
            if ((IS_WHITE_CLASS(gcap->labels[n]) && gcap->priors[n] > 0.1) ||
                IS_HYPO(gcap->labels[n]))
            {
              found = 1 ;
            }
          if (found == 0)  // no chance of wm or wmsa here
          {
            continue ;
          }

          if (label == Unknown)
          {
            DiagBreak() ;
          }
          load_val_vector(v_vals, mri_inputs, x, y, z) ;
          pwm = compute_conditional_density(m_inv_cov_wm, v_mean_wm, v_vals) ;
          pwmsa = compute_conditional_density(m_inv_cov_wmsa, v_mean_wmsa, v_vals) ;
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("         - pwm = %2.3e, pwmsa = %2.3e\n",
                   pwm, pwmsa) ;
          if (label == wm_label && pwmsa > pwm)
          {
            wm_dist = VectorDistance(v_mean_wm, v_vals) ;
            wmsa_dist = VectorDistance(v_mean_wmsa, v_vals) ;
            wm_mdist = MatrixMahalanobisDistance(v_mean_wm, m_inv_cov_wm, v_vals) ;
            wmsa_mdist = MatrixMahalanobisDistance(v_mean_wmsa, m_inv_cov_wmsa, v_vals) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              printf("         - wm_dist = %2.0f, wmsa_dist = %2.0f, mdists = (%2.0f, %2.0f)\n",
                     wm_dist, wmsa_dist, wm_mdist, wmsa_mdist) ;
            if ((wm_dist > wmsa_dist) && (wm_mdist > wmsa_mdist))
            {
              VectorSubtract(v_vals, v_mean_wm, v_dif_label) ;
              VectorSubtract(v_vals, v_mean_wmsa, v_dif_wmsa) ;
              if (
                ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < fabs(VECTOR_ELT(v_dif_label,1))) &&
                 (fabs(VECTOR_ELT(v_dif_wmsa,2)) < fabs(VECTOR_ELT(v_dif_label,2))) &&
                 (fabs(VECTOR_ELT(v_dif_wmsa,3)) < fabs(VECTOR_ELT(v_dif_label,3)))) ||
                ((2*wmsa_dist < wm_dist) && (2*wmsa_mdist < wm_mdist)))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("changing label from %s to %s\n",
                         cma_label_to_name(label),
                         cma_label_to_name(wmsa_label)) ;
                if (label == Unknown)
                {
                  DiagBreak() ;
                }
                label = wmsa_label ;
              }
            }
          }
          MRIsetVoxVal(mri_dst_labels, x, y, z, 0, label) ;
        }
      }
    }

    // now do 3 iterations of region growing
    for (i = 0 ; i < 3 ; i++)
    {
      mri_tmp = MRIcopy(mri_dst_labels, mri_tmp) ;
      for (x = 0 ; x < mri_inputs->width ; x++)
      {
        for (y = 0 ; y < mri_inputs->height ; y++)
        {
          for (z = 0 ; z < mri_inputs->depth ; z++)
          {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              DiagBreak() ;
            }
            label = MRIgetVoxVal(mri_dst_labels, x, y, z, 0) ;
            if (label != wm_label && label != Unknown && label != caudate_label)
            {
              continue ;
            }
            load_val_vector(v_vals, mri_inputs, x, y, z) ;
            nwmsa = MRIlabelsInNbhd(mri_dst_labels, x, y, z, 1, wmsa_label) ;
            if (nwmsa < 1)
            {
              continue ;
            }
            gcap = getGCAP(gca, mri_dst_labels, transform, x, y, z) ;
            for (found = n = 0 ; n < gcap->nlabels ; n++)
              if ((IS_WHITE_CLASS(gcap->labels[n]) && gcap->priors[n] > 0.1) ||
                  IS_HYPO(gcap->labels[n]))
              {
                found = 1 ;
              }
            if (found == 0)  // no chance of wm or wmsa here
            {
              continue ;
            }

            // only process it if it's in the body of the wm
#undef WSIZE
#define WSIZE 5
#define WHALF ((WSIZE-1)/2)

            nwm = MRIlabelsInNbhd(mri_tmp, x, y, z, WHALF, wm_label) ;
            nwmsa = MRIlabelsInNbhd(mri_tmp, x, y, z,WHALF, wmsa_label) ;
            nunknown = MRIlabelsInNbhd(mri_tmp, x, y, z,WHALF, Unknown) ;
            ncaudate = MRIlabelsInNbhd(mri_tmp, x, y, z,WHALF, caudate_label) ;
            ngm = MRIlabelsInNbhd(mri_tmp, x, y, z,WHALF, gm_label) ;

            if (ngm+ncaudate+nwm+nwmsa+nunknown < .9*WSIZE*WSIZE*WSIZE)  // somewhat arbitrary - the bulk of the nbhd
            {
              continue ;
            }
            ngm = MRIlabelsInNbhd(mri_tmp, x, y, z,1, gm_label) ;
            if (ngm > 0)  // not if there are any nearest nbrs that are gm
            {
              continue ;
            }
            if (nwm + nwmsa == 0)
            {
              continue ;
            }
            wm_dist = VectorDistance(v_mean_wm, v_vals) ;
            wmsa_dist = VectorDistance(v_mean_wmsa, v_vals) ;
            VectorSubtract(v_vals, v_mean_wmsa, v_dif_wmsa) ;
            if (label == caudate_label)
            {
              VectorSubtract(v_vals, v_mean_caudate, v_dif_label) ;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf("         - wm_dist = %2.0f, wmsa_dist = %2.0f\n",
                       wm_dist, wmsa_dist) ;
              if ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < fabs(VECTOR_ELT(v_dif_label,1))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,2)) < fabs(VECTOR_ELT(v_dif_label,2))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,3)) < fabs(VECTOR_ELT(v_dif_label,3))))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("changing label from %s to %s\n",
                         cma_label_to_name(label),
                         cma_label_to_name(wmsa_label)) ;
                label = wmsa_label ;
              }
            }
            else if (label == wm_label)
            {
              VectorSubtract(v_vals, v_mean_wm, v_dif_label) ;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf("         - wm_dist = %2.0f, wmsa_dist = %2.0f\n",
                       wm_dist, wmsa_dist) ;
              if (((fabs(VECTOR_ELT(v_dif_wmsa,1)) < fabs(VECTOR_ELT(v_dif_label,1))) &&
                   (fabs(VECTOR_ELT(v_dif_wmsa,2)) < fabs(VECTOR_ELT(v_dif_label,2))) &&
                   (fabs(VECTOR_ELT(v_dif_wmsa,3)) < fabs(VECTOR_ELT(v_dif_label,3)))) ||
                  (wmsa_dist*3 < wm_dist))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("changing label from %s to %s\n",
                         cma_label_to_name(label),
                         cma_label_to_name(wmsa_label)) ;
                if (label == Unknown)
                {
                  DiagBreak() ;
                }
                label = wmsa_label ;
              }
            }
            else if (label == Unknown)
            {
              VectorSubtract(v_vals, v_mean_un, v_dif_label) ;
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                printf("         - wm_dist = %2.0f, wmsa_dist = %2.0f\n",
                       wm_dist, wmsa_dist) ;
              if ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < fabs(VECTOR_ELT(v_dif_label,1))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,2)) < fabs(VECTOR_ELT(v_dif_label,2))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,3)) < fabs(VECTOR_ELT(v_dif_label,3))))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("changing label from %s to %s\n",
                         cma_label_to_name(label),
                         cma_label_to_name(wmsa_label)) ;
                if (label == Unknown)
                {
                  DiagBreak() ;
                }
                label = wmsa_label ;
              }
            }
            MRIsetVoxVal(mri_dst_labels, x, y, z, 0, label) ;
          }
        }
      }
    }
    MatrixFree(&m_cov_un) ;
    MatrixFree(&m_inv_cov_un) ;
    MatrixFree(&m_cov_wm) ;
    MatrixFree(&m_inv_cov_wm) ;
    MatrixFree(&m_cov_wmsa) ;
    MatrixFree(&m_inv_cov_wmsa) ;
    VectorFree(&v_mean_wm) ;
    VectorFree(&v_mean_wmsa) ;
    VectorFree(&v_mean_caudate) ;
    VectorFree(&v_mean_un) ;
  }
  VectorFree(&v_vals) ;
  VectorFree(&v_dif_label) ;
  VectorFree(&v_dif_wmsa) ;
  MRIfree(&mri_tmp) ;
  return(mri_dst_labels) ;
}

double
compute_conditional_density( MATRIX *m_inv_cov,
                             VECTOR *v_means,
                             VECTOR *v_vals )
{
  double  p, dist, det ;
  int     ninputs ;

  ninputs = m_inv_cov->rows ;

  det = MatrixDeterminant(m_inv_cov) ;
  dist = MatrixMahalanobisDistance(v_means, m_inv_cov, v_vals) ;
  p = (1.0 / (pow(2*M_PI,ninputs/2.0)*sqrt(1.0/det))) * exp(-0.5*dist) ;
  return(p) ;
}


int
load_val_vector( VECTOR *v_means,
                 MRI *mri_inputs,
                 int x, int y, int z )
{
  int  n ;

  for (n = 0 ; n < mri_inputs->nframes ; n++)
  {
    VECTOR_ELT(v_means, n+1) = MRIgetVoxVal(mri_inputs, x, y, z, n) ;
  }
  return(NO_ERROR) ;
}

#endif

int
GCAremoveWMSA( GCA *gca )
{
  int        x, y, z, n, found, i ;
  GCA_PRIOR  *gcap ;
  GCA_NODE   *gcan ;
  double     ptotal ;

  for (x = 0 ; x < gca->prior_width ; x++)
  {
    for (y = 0 ; y < gca->prior_height ; y++)
    {
      for (z = 0 ; z < gca->prior_depth ; z++)
      {
        gcap = &gca->priors[x][y][z] ;
        if (gcap==NULL)
        {
          continue;
        }
        found = 0 ;
        for (n = 0 ; n < gcap->nlabels ; n++)
        {
          if (IS_HYPO(gcap->labels[n]))
          {
            found = 1 ;
            gcap->priors[n] = 1e-10 ;
            break ;
          }
        }
        if (found) // renormalize priors
        {
          ptotal = 0 ;
          for (n = 0 ; n < gcap->nlabels ; n++)
          {
            ptotal += gcap->priors[n] ;
          }
          if (!FZERO(ptotal))
            for (n = 0 ; n < gcap->nlabels ; n++)
            {
              gcap->priors[n] /= ptotal ;
            }
        }
      }
    }
  }

  for (x = 0 ; x < gca->node_width ; x++)
  {
    for (y = 0 ; y < gca->node_height ; y++)
    {
      for (z = 0 ; z < gca->node_depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        if (gcan==NULL)
        {
          continue;
        }
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (IS_HYPO(gcan->labels[n]))
          {
            for (i = 0 ; i < gca->ninputs ; i++)
            {
              gcan->gcs[n].means[i] = -1e4 ;
            }
            break ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}

