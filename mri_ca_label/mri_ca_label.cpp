/**
 * @brief anisotropic nonstationary markov random field labeling
 *
 * Program for computing the MAP segmentation modeling the labeling as an
 * anisotropic nonstationary markov random field (based on manually labeled
 * data compiled into an atlas and stored in a .gca file)
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
#include <sys/utsname.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#ifdef HAVE_OPENMP // mrisurf.c has numerous parallelized functions
#include "romp_support.h"
#endif

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
#include "ventfix.h"

static char *write_likelihood = NULL ;
static double PRIOR_FACTOR = 1.0 ;

static char *read_renorm_fname = NULL ;
static char *write_renorm_fname = NULL ;

static double Gvent_topo_dist = 3 ;
static double Gvent_topo_volume_thresh1 = 50 ;
static double Gvent_topo_volume_thresh2 = 100 ;

static int remove_cerebellum = 0 ;
static int remove_lh = 0 ;
static int remove_rh = 0 ;

static int GCAremoveWMSA(GCA *gca) ;
static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static char *save_gca_fname = NULL ;

static char *surf_name = NULL ;
static char *surf_dir = NULL ;
static int rescale = 0 ;
static float Glabel_scales[MAX_CMA_LABELS] ;
static float Glabel_offsets[MAX_CMA_LABELS] ;

static float gsmooth_sigma = -1 ;
#define MAX_READS 100
static int nreads = 0 ;
static char *read_intensity_fname[MAX_READS] ;
static float regularize = 0 ;
static float regularize_mean = 0 ;
static int avgs = 0 ;
static int norm_PD = 0;
static int map_to_flash = 0 ;

static int reclassify_unlikely = 0 ;
static double unlikely_prior_thresh = 0.3 ;
static int unlikely_wsize = 9 ;
#if 0
static float unlikely_sigma = .5 ;
// more than this many sigmas from mean means unlikely:
static float unlikely_mah_dist_thresh = 4 ;  
#endif

static int wmsa = 0 ;   // apply wmsa postprocessing (using T2/PD data)
static int nowmsa = 0 ; // remove all wmsa labels from the atlas
static int fcd = 0 ;  // don't check for focal cortical dysplasias

static int handle_expanded_ventricles = 0;

static int renorm_with_histos = 0 ;

static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;

static MRI *GCArelabelUnlikely(GCA *gca,
                               MRI *mri_inputs,
                               TRANSFORM *transform,
                               MRI *mri_src_labeled,
                               MRI *mri_dst_labeled,
			       MRI *mri_independent_posterior,
                               double prior_thresh,
                               int whalf)  ;
static MRI *fix_putamen(GCA *gca,
                        MRI *mri_inputs,
                        MRI *mri_imp,
                        TRANSFORM *transform,
                        MRI *mri_src_labeled,
                        MRI *mri_dst_labeled,
                        double prior_thresh)  ;
MRI *insert_wm_bet_putctx(MRI *seg, int topo, const char *psfile, MRI *out);
int insert_wm_bet_putctx_topo = 0;
static MRI *gcaCheckForFCDs(MRI *mri_src, MRI *mri_dst, GCA *gca, TRANSFORM *transform, MRI *mri_inputs) ;
#if 0
static int load_val_vector(VECTOR *v_means,
                           MRI *mri_inputs,
                           int x, int y, int z) ;
static double compute_conditional_density(MATRIX *m_inv_cov,
    VECTOR *v_means,
    VECTOR *v_vals) ;
#endif

static MRI *replace_cortex_far_from_surface_with_wmsa(MRI *mri_inputs,
    MRI *mri_src,
    MRI *mri_dst,
    MRI_SURFACE *mris_lh,
    MRI_SURFACE *mris_rh) ;
static MRI *MRIexpandVentricles(MRI *mri_labeled_src,
                                MRI *mri_labeled,
                                GCA *gca,
                                MRI *mri_inputs,
                                int num_expansions) ;
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

const char *Progname ;
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
static char *wmsa_probs = NULL;

int FixVents=0;
int nitersFixVents = -1;
int nmaxFixVents = 5000;
int topoFixVents = 1;

#define CMA_PARCELLATION  0
static int parcellation_type = CMA_PARCELLATION ;

MRI *insert_wm_segmentation( MRI *mri_labeled, MRI *mri_wm,
                             int parcellation_type, int fixed_flag,
                             GCA *gca, TRANSFORM *transform );
int MRItoUCHAR(MRI **pmri);
extern char *gca_write_fname ;
extern int gca_write_iterations ;

//static int expand_flag = TRUE ;
static int expand_flag = FALSE ;
static int expand_ventricle_flag = FALSE ;
static int conform_flag = FALSE ;
struct utsname uts;
char *cmdline2, cwd[2000];
char *rusage_file=NULL;
char *PreGibbsFile=NULL;
int n_omp_threads;
MRI *InsertFromSeg=NULL;
std::vector<int> InsertFromSegIndices;
int MRIinsertFromSeg(MRI *mri_labeled, MRI *InsertFromSeg, std::vector<int> InsertFromSegIndices, int ZeroInsertIndex);
int InsertCblumFromSeg = -1;
int main(int argc, char *argv[])
{
  char         **av ;
  int          ac, nargs, extra ;
  char         *in_fname, *out_fname,  *gca_fname, *xform_fname;
  MRI          *mri_inputs, *mri_labeled, *mri_fixed = NULL, *mri_tmp, *mri_independent_posterior = NULL ;
  int          msec, minutes, seconds, ninputs, input ;
  Timer start ;
  GCA          *gca ;
  TRANSFORM     *transform ;


  FSinit() ;
  std::string cmdline = getAllInfo(argc, argv, "mri_ca_label");

  nargs = handleVersionOption(argc, argv, "mri_ca_label");
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
  fflush(stdout);
  fflush(stderr);

  Progname = argv[0];


  setRandomSeed(-1L) ;
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

#ifdef HAVE_OPENMP
  n_omp_threads = omp_get_max_threads();
  printf("\n== Number of threads available to for OpenMP = %d == \n",n_omp_threads);
#else
  printf("Do not have OpenMP\n");
  n_omp_threads = 1;
#endif

  if (getenv("BUILD_GCA_HISTO") != NULL)
  {
    int *counts, i, max_i ;

    fprintf(stderr, "reading gca from %s\n", argv[1]) ;
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

  printf("reading %d input volumes\n", ninputs) ;

  /*  fprintf(stderr,
      "mri_inputs read: xform %s\n", mri_inputs->transform_fname) ;*/
  printf("reading classifier array from %s\n", gca_fname) ;
  fflush(stdout);
  fflush(stderr);
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

  if (gsmooth_sigma > 0)
  {
    GCA *gca_smooth ;
    gca_smooth = GCAsmooth(gca, gsmooth_sigma) ;
    GCAfree(&gca) ;
    gca = gca_smooth ;
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
    printf("reading input volume from %s\n", in_fname) ;
    fflush(stdout);
    fflush(stderr);
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
  fflush(stdout);
  fflush(stderr);

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
  if (rescale)
  {
    GCAapplyRenormalization(gca, Glabel_scales, Glabel_offsets, 0) ;
  }
  fflush(stdout);
  fflush(stderr);
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
    printf("reading transform from %s\n", xform_fname) ;
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

  if (Ggca_x >= 0 && Gx < 0)
  {
    GCAsourceVoxelToNode(gca, mri_inputs, transform, 
                         Ggca_x, Ggca_y, Ggca_z,
                         &Gx, &Gy, &Gz) ;
    printf("source voxel (%d, %d, %d) maps to node (%d, %d, %d)\n",
           Ggca_x, Ggca_y, Ggca_z, Gx, Gy, Gz) ;
    GCAdump(gca, mri_inputs, Ggca_x, Ggca_y, Ggca_z, transform, stdout, 0) ;
  }
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
    printf("reading %d labels from %s\n", nlines,renormalization_fname) ;
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
  if (norm_PD)
  {
    GCAnormalizePD(gca, mri_inputs, transform) ;
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
      fprintf(stderr, "writing equalized volume to %s\n", "heq.mgz") ;
      MRIwrite(mri_inputs, "heq.mgz") ;
    }
  }
  fflush(stdout);
  fflush(stderr);

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
             cma_label_to_name(MRIgetVoxVal(mri_labeled, Ggca_x, Ggca_y, Ggca_z,0)),
             (int)MRIgetVoxVal(mri_labeled, Ggca_x, Ggca_y, Ggca_z,0),
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
     mri_fixed, 0, NULL, PRIOR_FACTOR, PRIOR_FACTOR);
  }
  else    // don't read old one in - create new labeling
  {
    if (reg_fname == NULL)   //so read_fname must be NULL too
    {
      printf("labeling volume...\n") ;
      // create labeled volume by MAP rule
      mri_labeled = GCAlabel(mri_inputs, gca, NULL, transform) ;
      mri_independent_posterior = MRIcopy(mri_labeled, NULL) ;
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
                  "writing patched labeling to %s\n", out_fname) ;
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
        printf("writing snapshot to %s\n", fname) ;
        MRIwrite(mri_labeled, fname) ;
      }
      // renormalize iteration
      if (renormalize_align)
      {
        FILE *logfp ;
        char base_name[STRLEN] ;

        //FileNameOnly(out_fname, base_name) ;
        FileNameRemoveExtension(out_fname, base_name) ;
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

        if (read_renorm_fname)
        {
          GCAfree(&gca) ;
          gca = GCAread(read_renorm_fname) ;
        }
        else
        {
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
	  if (Ggca_x >= 0)
	    printf("label(%d, %d, %d) = %s (%d), norm=%2.0f\n",
		   Ggca_x, Ggca_y, Ggca_z,
		   cma_label_to_name(MRIgetVoxVal(mri_labeled, Ggca_x, Ggca_y, Ggca_z, 0)),
		   (int)MRIgetVoxVal(mri_labeled, Ggca_x, Ggca_y, Ggca_z,0),
		   MRIgetVoxVal(mri_inputs, Ggca_x, Ggca_y, Ggca_z, 0)) ;
          if (write_renorm_fname)
          {
            GCAwrite(gca, write_renorm_fname) ;
          }
        }

        if (regularize_mean > 0)
        {
          GCAregularizeConditionalDensities(gca, regularize_mean) ;
        }
        if (save_gca_fname)
        {
          GCAwrite(gca, save_gca_fname) ;
        }
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
          printf("reading %d labels from %s\n",
                 nlines,renormalization_fname) ;
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
        fflush(stdout);
        fflush(stderr);
        GCAlabel(mri_inputs, gca, mri_labeled, transform) ;

        /*
          try to create a better starting labeling by relabeling regions 
          of a label that disagree with the atlas label all at once 
          (the MRF stuff is prone to local minima where flipping one
          label increases the energy, but flipping a bunch all at 
          once lowers it)
        */
        if (reclassify_unlikely)
        {
          GCArelabelUnlikely(gca, mri_inputs, transform, 
                             mri_labeled, mri_labeled, mri_independent_posterior,
                             unlikely_prior_thresh, (unlikely_wsize-1)/2) ;
        }
        {
          MRI *mri_imp ;
          mri_imp = GCAmarkImpossible(gca, mri_labeled, NULL, transform) ;
	  
	  // This is never run because of the -0.1
          fix_putamen(gca, mri_inputs, mri_imp, transform,  mri_labeled, mri_labeled,-.1) ;
          if (Gdiag & DIAG_WRITE)
          {
            MRIwrite(mri_imp, "gca_imp.mgz") ;
          }
	  // Do this after Gibbs
	  //insert_wm_bet_putctx(mri_labeled, insert_wm_bet_putctx_topo, "nofile", mri_labeled);
          MRIfree(&mri_imp) ;
        }
      } // renormalize_align
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
      ErrorExit(ERROR_BADPARM,
                "%s ERROR: the -l option is currently not supported."
                " Debugging and testing needed.\n",
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
               cma_label_to_name(MRIgetVoxVal
                                 (mri_labeled, Ggca_x, Ggca_y, Ggca_z,0)),
               (int)MRIgetVoxVal(mri_labeled, Ggca_x, Ggca_y, Ggca_z,0)) ;
      TransformFree(&transform_long) ;
      if (gca_write_iterations != 0)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s_pre.mgz", gca_write_fname) ;
        printf("writing snapshot to %s\n", fname) ;
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

    if(PreGibbsFile){
      printf("Saving pre-gibbs to %s\n",PreGibbsFile);
      int err;
      err = MRIwrite(mri_labeled,PreGibbsFile);
      if(err) exit(1);
    }

    if (!no_gibbs)
    {
      if (anneal)
      {
        GCAanneal(mri_inputs, gca, mri_labeled, 
                  transform, max_iter, PRIOR_FACTOR) ;
      }
      else
      {
        if (Ggca_x >= 0)
        {
          GCAdump(gca, mri_inputs, Ggca_x, Ggca_y, Ggca_z,
                  transform, stdout, 0) ;
        }
	printf("Reclassifying using Gibbs Priors\n");
        GCAreclassifyUsingGibbsPriors
        (mri_inputs, gca, mri_labeled, transform, max_iter,
         mri_fixed, 0, NULL, PRIOR_FACTOR, PRIOR_FACTOR);
        if (reclassify_unlikely)
        {
	  int w ;

	  for (w = (unlikely_wsize-1)/2 ; w >= 1 ; w--)
	  {
	    GCArelabelUnlikely(gca, mri_inputs, transform, 
			       mri_labeled, mri_labeled, mri_independent_posterior,
			       unlikely_prior_thresh, w) ;
	    if (gca_write_iterations != 0)
	    {
	      char fname[STRLEN] ;
	      sprintf(fname, "%s_un.w%d.mgz", gca_write_fname, w) ;
	      printf("writing snapshot to %s\n", fname) ;
	      MRIwrite(mri_labeled, fname) ;
	    }
	  }
        }
      }
    }
  }

  if (read_fname == NULL && 0)
  {
    GCAmaxLikelihoodBorders(gca, mri_inputs, mri_labeled, mri_labeled,
                            transform,mle_niter, 5.0);
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

  if (handle_expanded_ventricles)
  {
    MRIexpandVentricles(mri_labeled, mri_labeled, gca, mri_inputs, 10) ;
  }
  if (gca_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_post.mgz", gca_write_fname) ;
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
  }

  if (0 && reclassify_unlikely) {
    MRI *mri_sigma ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_labeled, "aseg_before.mgz") ;
    }
    mri_sigma = GCAcomputeOptimalScale(gca, transform, mri_inputs,
                                       mri_labeled, NULL,
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
    GCAconstrainLabelTopology(gca, mri_inputs, 
                              mri_labeled, mri_labeled, transform,
			      Gvent_topo_dist, Gvent_topo_volume_thresh1,
			      Gvent_topo_volume_thresh2) ;
  }
  if (wmsa)
  {
    MRI *mri_tmp ;
    if (surf_name)
    {
      char fname[STRLEN] ;
      MRI_SURFACE *mris_lh, *mris_rh ;
      sprintf(fname, "%s/lh.%s", surf_dir, surf_name) ;
      mris_lh = MRISread(fname) ;
      if (mris_lh == NULL)
      {
        ErrorExit(ERROR_NOFILE, 
                  "%s: could not read lh from %s",
                  Progname, fname) ;
      }
      sprintf(fname, "%s/rh.%s", surf_dir, surf_name) ;
      mris_rh = MRISread(fname) ;
      if (mris_rh == NULL)
      {
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read rh from %s",
                  Progname, fname) ;
      }
      replace_cortex_far_from_surface_with_wmsa(mri_inputs,
                                                mri_labeled,
                                                mri_labeled,
                                                mris_lh,
                                                mris_rh) ;

    }
    mri_tmp = GCAlabelWMandWMSAs(gca, mri_inputs, mri_labeled, NULL, transform);
    MRIfree(&mri_labeled) ;
    mri_labeled = mri_tmp ;
  }
  if (wmsa_probs != NULL)
//if (wmsa_probs)
  {
    MRI *probs_tmp ;
    char fname[STRLEN] ;
    printf("working on WMSA probs now...\n") ;
    sprintf(fname, "%s.mgz", wmsa_probs) ;
    probs_tmp=GCAsampleToVolWMSAprob(mri_inputs, gca, transform, NULL);
    printf("writing WMSA probability volume to %s....\n", fname) ;
    //MRIwrite(probs_tmp,fname);
    if (MRIwrite(probs_tmp, fname) != NO_ERROR)
    {
      ErrorExit(Gerror, "%s: MRIwrite(%s) failed", Progname, fname);
    }
  }

  // embed color lookup table from GCA
  mri_labeled->ct = gca->ct;

  // if GCA has no ctab, embed the default instead
  if (!mri_labeled->ct) mri_labeled->ct = CTABreadDefault();

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
  if (fcd)
    gcaCheckForFCDs(mri_labeled, mri_labeled, gca, transform, mri_inputs) ;

  if(insert_wm_bet_putctx_topo>0){
    printf("Inserting WM between putamen and cortex topo=%d\n",insert_wm_bet_putctx_topo);
    insert_wm_bet_putctx(mri_labeled, insert_wm_bet_putctx_topo, "nofile", mri_labeled);
  }

  if(FixVents){
    printf("Fixing vents  niters = %d  nmax = %d topo = %d\n", nitersFixVents, nmaxFixVents, topoFixVents);
    char segids[5] = "4,43";
    MRI *newseg = VentFix::fixasegps(mri_labeled, mri_inputs, &(segids[0]), 0.5, nitersFixVents, nmaxFixVents, topoFixVents);
    MRIfree(&mri_labeled);
    mri_labeled = newseg;
  }

  if(InsertFromSeg){
    printf("Inserting from seg\n");
    MRIinsertFromSeg(mri_labeled, InsertFromSeg, InsertFromSegIndices, InsertCblumFromSeg);
  }

  printf("writing labeled volume to %s\n", out_fname) ;
  if (MRIwrite(mri_labeled, out_fname) != NO_ERROR)
  {
    ErrorExit(Gerror, "%s: MRIwrite(%s) failed", Progname, out_fname) ;
  }

  if (G_write_probs != NULL)
  {
    MRI *mri_probs ;
    char fname[STRLEN] ;

    mri_probs = GCAcomputeTopNProbabilities(mri_inputs,gca,mri_labeled,NULL,transform,3);
    sprintf(fname, "%s.mgz", G_write_probs) ;
    printf("writing ordered labels and  probabilities to %s\n", fname) ;
    MRIwrite(mri_probs, fname) ;
    MRIfree(&mri_probs) ;
  }
  if (write_likelihood != NULL)
  {
    MRI *mri_probs ;
    char fname[STRLEN] ;

    mri_probs = GCAcomputeLikelihoodImage(gca, mri_inputs, mri_labeled, transform);
    sprintf(fname, "%s", write_likelihood) ;
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
      MRIthresholdMask(mri_probs, mri_mask, mri_probs, 1, 0) ; 
    }
    printf("writing likelihood image to to %s\n", fname) ;
    MRIwrite(mri_probs, fname) ;
    MRIfree(&mri_probs) ;
  }
  MRIfree(&mri_inputs) ;
  GCAfree(&gca);

  // Print usage stats to the terminal (and a file is specified)
  PrintRUsage(RUSAGE_SELF, "mri_ca_label ", stdout);
  if(rusage_file) WriteRUsage(RUSAGE_SELF, "", rusage_file);

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("mri_ca_label took %d minutes and %d seconds.\n", minutes, seconds) ;
  printf("mri_ca_label done\n");
  exit(0);
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
  else if (!stricmp(option, "THREADS"))
  {
    sscanf(argv[2],"%d",&n_omp_threads);
    #ifdef HAVE_OPENMP
    omp_set_num_threads(n_omp_threads);
    printf("Setting threads to %d\n",n_omp_threads);
    #else
    printf("dont have openmp \n");
    #endif
    nargs = 1 ;
  }
  else if (!stricmp(option, "PREGIBBS"))
  {
    PreGibbsFile = argv[2];
    printf("Saving pre-gibbs to %s\n",PreGibbsFile);
    nargs = 1 ;
  }
  else if (!stricmp(option, "LH"))
  {
    remove_rh = 1  ;
    printf("removing right hemisphere labels\n") ;
  }
  else if (!stricmp(option, "vent_topo_dist"))
  {
    Gvent_topo_dist = atof(argv[2]) ;
    printf("setting ventricle topology distance threshold to %2.1fmm (default=3)\n", Gvent_topo_dist) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "vent_topo_volume_thresh1"))
  {
    Gvent_topo_volume_thresh1 = atof(argv[2]) ;
    printf("setting ventricle topology volume1 threshold to %2.1fmm^3 (default=50)\n", Gvent_topo_volume_thresh1) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "vent_topo_volume_thresh2"))
  {
    Gvent_topo_volume_thresh2 = atof(argv[2]) ;
    printf("setting ventricle topology volume2 threshold to %2.1fmm^3 (default=100)\n", Gvent_topo_volume_thresh2) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "RH"))
  {
    remove_lh = 1  ;
    printf("removing left hemisphere labels\n") ;
  }
  else if (!stricmp(option, "PRIOR"))
  {
    PRIOR_FACTOR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using Gibbs prior factor = %2.3f\n", PRIOR_FACTOR) ;
  }
  else if (!stricmp(option, "NOCEREBELLUM"))
  {
    remove_cerebellum = 1 ;
    printf("removing cerebellum from atlas\n") ;
  }
  else if (!stricmp(option, "nowmsa"))
  {
    nowmsa = 1 ;
    printf("disabling WMSA labels\n") ;
  }
  else if(!stricmp(option, "insert-from-seg") || !stricmp(option, "sa-insert-from-seg"))
  {
    // -insert-from-seg    InserFromSeg.mgz index1 <index2 ...>
    // -sa-insert-from-seg InserFromSeg.mgz index1 <index2 ...> InputSeg OutputSeg
    InsertFromSeg = MRIread(argv[2]);
    if(!InsertFromSeg) exit(1);
    printf("Inserting from seg %s  ",argv[2]);
    nargs = 1;
    int k=3;
    while(1){
      if(!argv[k]) break;
      if(!isdigit(argv[k][0])) break;
      InsertFromSegIndices.push_back(atoi(argv[k]));
      printf("%s ",argv[k]);
      k++;
      nargs++;
    }
    printf("\n");
    if(InsertFromSegIndices.size()==0) {
      printf("ERROR: -insert-from-seg needs at least one index to insert\n");
      exit(1);
    }
    MRI *InputSeg=NULL;
    if(!stricmp(option, "sa-insert-from-seg")){
      if(argv[k]==NULL || argv[k+1]==NULL){
	printf("ERROR: -sa-insert-from-seg needs an input and output seg\n");
	exit(1);
      }
      InputSeg = MRIread(argv[k]);
      if(!InputSeg) exit(1);
      printf("Inserting from seg\n");
      MRIinsertFromSeg(InputSeg, InsertFromSeg, InsertFromSegIndices,-1);
      int err = MRIwrite(InputSeg,argv[k+1]);
      exit(err);
    }
    printf("\n") ;
  }
  else if(!stricmp(option, "cblum-from-seg") || !stricmp(option, "sa-cblum-from-seg"))
  {
    // same as insert-from-seg, but uses cblum cortex and wm and also zeros voxels
    // that are cblum in the input but CSF (24) in the output. This is a hack to 
    // preventa bunch of random CSF voxels in the output
    // -cblum-from-seg    InserFromSeg.mgz 
    // -sa-cblum-from-seg InserFromSeg.mgz InputSeg OutputSeg
    InsertCblumFromSeg = 24;
    InsertFromSeg = MRIread(argv[2]);
    if(!InsertFromSeg) exit(1);
    printf("Inserting cblum from seg %s  ",argv[2]);
    InsertFromSegIndices.push_back(7);
    InsertFromSegIndices.push_back(8);
    InsertFromSegIndices.push_back(46);
    InsertFromSegIndices.push_back(47);
    nargs = 1;
    MRI *InputSeg=NULL;
    if(!stricmp(option, "sa-cblum-from-seg")){
      if(argv[2]==NULL || argv[3]==NULL){
	printf("ERROR: -sa-cblum-from-seg needs an input and output seg\n");
	exit(1);
      }
      InputSeg = MRIread(argv[3]);
      if(!InputSeg) exit(1);
      printf("Inserting from seg\n");
      MRIinsertFromSeg(InputSeg, InsertFromSeg, InsertFromSegIndices,24);
      printf("Writing to %s\n",argv[4]);
      int err = MRIwrite(InputSeg,argv[4]);
      printf("mi_ca_label done\n") ;
      exit(err);
    }
    printf("\n") ;
  }


  else if (!stricmp(option, "insert-wm-bet-putctx")){
    sscanf(argv[2],"%d",&insert_wm_bet_putctx_topo);
    nargs = 1 ;
  }
  else if (!stricmp(option, "sa-insert-wm-bet-putctx")){
    // stand-alone option to apply insert_wm_bet_putctx()
    // 2=segvol 3=topo 4=outputsegvol 5=psfile
    MRI *seg = MRIread(argv[2]);
    if(seg==NULL) exit(1);
    int topo;
    sscanf(argv[3],"%d",&topo);
    insert_wm_bet_putctx(seg,topo,argv[5],seg);
    printf("Writing to %s\n",argv[4]);fflush(stdout);
    int err = MRIwrite(seg,argv[4]);
    exit(err);
  }
  else if (!stricmp(option, "fcd"))
  {
    fcd = 1 ;
//    reclassify_unlikely = 0 ;
    printf("checking for focal cortical dysplasias\n") ;
  }
  else if (!stricmp(option, "read_intensities") || !stricmp(option, "ri"))
  {
    read_intensity_fname[nreads] = argv[2] ;
    nargs = 1 ;
    printf("reading intensity scaling from %s\n",
           read_intensity_fname[nreads]) ;
    nreads++ ;
    if (nreads > MAX_READS)
    {
      ErrorExit(ERROR_UNSUPPORTED, 
                "%s: too many intensity files specified (max %d)",
                Progname, MAX_READS);
    }
  }
  else if (!stricmp(option, "WM"))
  {
    wm_fname = argv[2] ;
    nargs = 1 ;
    printf("inserting white matter segmentation from %s\n", wm_fname) ;
  }
  else if (!stricmp(option, "SURF"))
  {
    surf_dir = argv[2] ;
    surf_name = argv[3] ;
    nargs = 2 ;
    printf("using surfaces ?h.%s in directory %s\n", surf_name, surf_dir) ;
  }
  else if (!stricmp(option, "SD"))
  {
    setenv("SUBJECTS_DIR",argv[2],1) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "rusage"))
  {
    // resource usage
    rusage_file = argv[2] ;
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
#if 0
    reclassify_unlikely = atoi(argv[2]) ;
    unlikely_wsize = atoi(argv[3]) ;
    unlikely_sigma = atof(argv[4]) ;
    unlikely_mah_dist_thresh = atof(argv[5]) ;
    if (reclassify_unlikely)
      printf("relabeling unlikely voxels more than "
             "%2.1f sigmas from label mean, with sigma=%2.1fmm, wsize=%dmm\n",
             unlikely_mah_dist_thresh, unlikely_sigma, unlikely_wsize) ;
    else
    {
      printf("not relabeling unlikely voxels\n") ;
    }
    nargs = 4 ;
#else
    reclassify_unlikely = 1 ;
    unlikely_wsize = atoi(argv[2]) ;
    unlikely_prior_thresh = atof(argv[3]) ;
    printf("relabeling unlikely voxels with "
           "window_size = %d and prior threshold %2.2f\n",
           unlikely_wsize, unlikely_prior_thresh) ;
    nargs = 2 ;
#endif
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
    printf("writing label probabilities to %s\n", G_write_probs) ;
  }
  else if (!stricmp(option, "write_likelihood"))
  {
    write_likelihood = argv[2] ;
    nargs = 1 ;
    printf("writing image likelihoods unders labeling to %s\n", 
	   write_likelihood) ;
  }
  else if (!stricmp(option, "wmsa_probs"))
  {
    wmsa_probs = argv[2] ;
    //wmsa_probs = 1;
    printf("Writing WMSA probabilities to file %s\n", wmsa_probs) ;
    nargs =1 ;
    //sprintf(wmsaprob_fname,"%s.mgz",wmsa_probs);
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
    printf("debugging prior (%d, %d, %d)\n", Gxp,Gyp,Gzp) ;
  }
  else if (!stricmp(option, "VENT-FIX"))
  {
    /* This labels underlabeled vertices in ventricles. It iteratively
       grows the ventricle into 0-valued voxels until one of three
       stopping criteria is reached: (1) number of iters exceeds
       niters (if niters < 0, then no max iters). (2) the number of
       changed voxels exceeds nmax. or (3) there are no more 0-valued
       voxels that neighbor a ventricle voxel. A neighbor is defined
       by the topo: 1=face, 2=face+edge, 3=face+edge+corner. Typical
       -1 7000 1. This may fail if the unlabeled ventricle is not
       completely surrounded by non-zero segments. */
    FixVents = 1;
    nitersFixVents = atoi(argv[2]) ; // -1
    nmaxFixVents = atoi(argv[3]) ; // 7000
    topoFixVents = atoi(argv[4]) ; // 1
    nargs = 3 ;
  }
  else if (!stricmp(option, "SA-VENT-FIX"))
  {
    // -sa-vent-fix niters nmax topo inseg brainmask outseg
    nitersFixVents = atoi(argv[2]) ; // -1
    nmaxFixVents = atoi(argv[3]) ; // 7000
    topoFixVents = atoi(argv[4]) ; // 1
    MRI *inseg = MRIread(argv[5]); 
    if(inseg == NULL) exit(1);
    MRI *brainmask = MRIread(argv[6]); 
    if(brainmask == NULL) exit(1);
    printf("Fixing vents  niters = %d  nmax = %d topo = %d\n", nitersFixVents, nmaxFixVents, topoFixVents);
    char segids[5] = "4,43";
    MRI *newseg = VentFix::fixasegps(inseg, brainmask, &(segids[0]), 0.5, nitersFixVents, nmaxFixVents, topoFixVents);
    int err = MRIwrite(newseg,argv[7]);
    exit(err);
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %d\n", Ggca_label) ;
  }
  else if (!stricmp(option, "DEBUG"))
  {
    Gdiag = DIAG_WRITE | DIAG_VERBOSE_ON;
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
  else if (!stricmp(option, "GSMOOTH"))
  {
    gsmooth_sigma = atof(argv[2]) ;
    printf("smoothing atlas with a Gaussian with sigma = %2.2f mm\n",
           gsmooth_sigma) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "LSCALE"))
  {
    int l ;
    l = atoi(argv[2]) ;
    rescale = 1 ;
    Glabel_scales[l] = atof(argv[3]) ;
    nargs = 2 ;
    printf("scaling label %s by %2.2f\n",
           cma_label_to_name(l), Glabel_scales[l]) ;
    for (l = 0 ; l < MAX_CMA_LABELS ; l++)
      if (FZERO(Glabel_scales[l]))
      {
        Glabel_scales[l] = 1.0 ;
      }
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
    printf("inserting fixed white matter segmentation from %s\n",
           wm_fname);
  }
  else if (!stricmp(option, "MRI"))
  {
    mri_fname = argv[2] ;
    nargs = 1 ;
    printf("building most likely MR volume and writing to %s\n",
           mri_fname);
  }
  else if (!stricmp(option, "HEQ"))
  {
    heq_fname = argv[2] ;
    nargs = 1 ;
    printf("reading template for histogram equalization from %s\n",
           heq_fname) ;
  }
  else if (!strcmp(option, "RENORM") || !stricmp(option, "RENORMALIZE"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s\n",
           renormalization_fname) ;
  }
  else if (!stricmp(option, "WRITE_RENORM"))
  {
    write_renorm_fname = argv[2] ;
    nargs = 1 ;
    printf("writing renormalized GCA to in %s\n", write_renorm_fname) ;
  }
  else if (!stricmp(option, "READ_RENORM"))
  {
    read_renorm_fname = argv[2] ;
    nargs = 1 ;
    printf("reading renormalized GCA from %s\n", read_renorm_fname) ;
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
  else if (!stricmp(option, "renormalize_iter"))
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
      printf("reading previously labeled volume from %s\n", read_fname) ;
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
        label = MRIgetVoxVal(mri_labeled, x, y, z,0) ;
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
          if (label != MRIgetVoxVal(mri_labeled, x, y, z,0))
          {
            nchanged++ ;
          }
          MRIsetVoxVal(mri_labeled, x, y, z, 0, label) ;
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
        if (MRIgetVoxVal(mri, xi, yi, zi,0) == label)
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
    printf("writing snapshot to %s\n", fname) ;
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
      printf("writing snapshot to %s\n", fname) ;
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
    if (MRIgetVoxVal(mri_labeled, xi, yi, zi,0) == label)
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
        label = MRIgetVoxVal(mri_labeled, x, y, z,0) ;
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
            if (MRIgetVoxVal(mri_labeled, xi, y, z,0) ==
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
              MRIsetVoxVal(mri_tmp, xi, y, z,0, Left_Hippocampus) ;
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
            if (MRIgetVoxVal(mri_labeled, xi, y, z, 0) ==
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
              MRIsetVoxVal(mri_tmp, xi, y, z, 0,  Right_Hippocampus) ;
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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Left_Cerebral_White_Matter) ;
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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Right_Cerebral_White_Matter) ;
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
          label = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;

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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Left_Cerebral_Cortex) ;
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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Right_Cerebral_Cortex) ;
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
          label = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;

          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            left = 1 ;
#if __GNUC__  >= 8
	    [[gnu::fallthrough]];
#endif
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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, label) ;
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
        label = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;

        left = 0 ;
        was_unknown = 0 ;
        switch (label)
        {
        case Left_Hippocampus:
          left = 1 ;
#if __GNUC__  >= 8
      [[gnu::fallthrough]];
#endif
        case Right_Hippocampus:
#define MIN_UNKNOWN 5
          for (i = 1 ; i <= MIN_UNKNOWN+3 ; i++)
          {
            yi = mri_tmp->yi[y+i] ;
            olabel = MRIgetVoxVal(mri_tmp, x, yi, z, 0) ;
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
            MRIsetVoxVal(mri_tmp, x, y, z, 0, left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf
              ("(%d, %d, %d) %s changed to %s "
               "in edit_hippocampus\n",
               x, y, z, cma_label_to_name(label),
               cma_label_to_name(MRIgetVoxVal(mri_tmp, x, y, z,0))) ;
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
        label = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;

        left = 0 ;
        found_wm = found_hippo = was_wm = was_hippo = 0 ;
        switch (label)
        {
        case Left_Hippocampus:
          left = 1 ;
#if __GNUC__  >= 8
      [[gnu::fallthrough]];
#endif
        case Right_Hippocampus:
          for (i = 1 ; i <= 10 ; i++)
          {
            yi = mri_tmp->yi[y-i] ;
            olabel = MRIgetVoxVal(mri_tmp, x, yi, z, 0) ;
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
            MRIsetVoxVal(mri_tmp, x, y, z, 0, 
                         left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf("(%d, %d, %d) %s changed to "
                     "%s in edit_hippocampus\n",
                     x, y, z, cma_label_to_name(label),
                     cma_label_to_name(MRIgetVoxVal(mri_tmp, x, y, z,0))) ;
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

          label = MRIgetVoxVal(mri_labeled, x, y, z, 0) ;

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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Left_Hippocampus) ;
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
              MRIsetVoxVal(mri_tmp, x, y, z, 0, Right_Hippocampus) ;
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

          label = MRIgetVoxVal(mri_labeled, x, y, z, 0) ;

          left = 0 ;
          switch (label)
          {
          case Left_Hippocampus:
            xi = mri_labeled->xi[x+1] ;
            if (MRIgetVoxVal(mri_labeled, xi, y, z, 0) ==
                Left_Cerebral_Cortex)
            {
              found_hippo = 1 ;
              for (i = 0 ; found_hippo && i < 8 ; i++)
              {
                xi = mri_labeled->xi[x-i] ;
                if (MRIgetVoxVal(mri_labeled, xi, y, z, 0) !=
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
                MRIsetVoxVal(mri_labeled, xi, y, z,0, Left_Hippocampus) ;
              }
            }
            break ;
          case Right_Hippocampus:
            xi = mri_labeled->xi[x+1] ;
            if (MRIgetVoxVal(mri_tmp, xi, y, z, 0) == Right_Cerebral_Cortex)
            {
              found_hippo = 1 ;
              for (i = 0 ; found_hippo && i < 8 ; i++)
              {
                xi = mri_labeled->xi[x+i] ;
                if (MRIgetVoxVal(mri_labeled, xi, y, z, 0) !=
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
                MRIsetVoxVal(mri_labeled, xi, y, z, 0, Right_Hippocampus) ;
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
        label = MRIgetVoxVal(mri_tmp, x, y, z, 0) ;

        left = 0 ;
        found_wm = found_amygdala = was_wm = was_amygdala = 0 ;
        switch (label)
        {
        case Left_Amygdala:
          left = 1 ;
#if __GNUC__  >= 8
      [[gnu::fallthrough]];
#endif
        case Right_Amygdala:
          for (i = 1 ; i <= 10 ; i++)
          {
            yi = mri_tmp->yi[y-i] ;
            olabel = MRIgetVoxVal(mri_tmp, x, yi, z, 0) ;
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
            MRIsetVoxVal(mri_tmp, x, y, z, 0, label) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              printf("(%d, %d, %d) %s changed to "
                     "%s in edit_amygdala\n",
                     x, y, z, cma_label_to_name(label),
                     cma_label_to_name(MRIgetVoxVal(mri_tmp, x, y, z, 0))) ;
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
  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
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
              gcas_setPrior(gcas[i], getPrior(gcap, gcap->labels[n]));
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
    MRIsetVoxVal(mri_probs, gcas[i].x, gcas[i].y, gcas[i].z, 0, p) ;
    if (gcas[i].x == Ggca_x && gcas[i].y == Ggca_y && gcas[i].z == Ggca_z)
    {
      DiagBreak() ;
    }
    MRIsetVoxVal(mri_tmp_labels, 
                 gcas[i].x, gcas[i].y, gcas[i].z,
                 0,
                 gcas[i].label) ;
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
        if (MRIgetVoxVal(mri_tmp_labels,x,y,z, 0) != Left_Cerebral_White_Matter)
        {
          continue ;
        }
        yi = mri_probs->yi[y+1] ;
        p = MRIgetVoxVal(mri_probs, x,  y, z, 0) + 
          MRIgetVoxVal(mri_probs, x, yi, z,0) ;
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
  MRIsetVoxVal(mri_tmp, ximax, yimax, zimax, 0, 255) ;
  MRIsetVoxVal(mri_tmp, ximax, yimax+1, zimax, 0, 255) ;
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
          if (MRIgetVoxVal(mri_tmp, x, y, z, 0))
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
            MRIgetVoxVal(mri_probs, x,  y, z, 0) +
            MRIgetVoxVal(mri_probs, x, yi, z, 0) ;
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

        MRIsetVoxVal(mri_tmp, x, yimax, z, 0, 255) ;
        MRIsetVoxVal(mri_tmp, x, yimax+1, z, 0, 255) ;
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
        if (MRIgetVoxVal(mri_tmp_labels,x,y,z,0) != Right_Cerebral_White_Matter)
        {
          continue ;
        }
        yi = mri_probs->yi[y+1] ;
        p = MRIgetVoxVal(mri_probs, x,  y, z, 0) + 
          MRIgetVoxVal(mri_probs, x, yi, z, 0) ;
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
  MRIsetVoxVal(mri_tmp, ximax, yimax, zimax, 0, 255) ;
  MRIsetVoxVal(mri_tmp, ximax, yimax+1, zimax, 0, 255) ;
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
          if (MRIgetVoxVal(mri_tmp, x, y, z, 0))
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
            MRIgetVoxVal(mri_probs, x,  y, z, 0) +
            MRIgetVoxVal(mri_probs, x, yi, z,0) ;
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

        MRIsetVoxVal(mri_tmp, x, yimax, z, 0, 255) ;
        MRIsetVoxVal(mri_tmp, x, yimax+1, z, 0,  255) ;
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
          if (MRIgetVoxVal(mri_tmp, x, y, z, 0))
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
          if ((MRIgetVoxVal(mri_tmp, x+1, y, z, 0) &&
               MRIgetVoxVal(mri_tmp, x-1,y,z, 0)) ||
              (MRIgetVoxVal(mri_tmp, x, y, z-1, 0) &&
               MRIgetVoxVal(mri_tmp, x,y,z+1,0)))
          {
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            {
              DiagBreak() ;
            }

            nchanged++ ;
            MRIsetVoxVal(mri_tmp, x,y,z, 0,
                         MAX(MAX(MAX(MRIgetVoxVal(mri_tmp, x+1, y, z, 0),
                                     MRIgetVoxVal(mri_tmp, x-1,y,z, 0)),
                                 MRIgetVoxVal(mri_tmp, x, y, z-1, 0)),
                             MRIgetVoxVal(mri_tmp, x,y,z+1, 0))) ;
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
          if (MRIgetVoxVal(mri_tmp, x, y, z,0))
          {
            continue ;
          }
          yi = mri_tmp->yi[y+1] ;
          if (MRIgetVoxVal(mri_tmp,x,yi,z,0)) /* check inferior */
          {
            yi = mri_tmp->yi[y-1] ;
            zi = mri_tmp->zi[z+1] ;
            if (MRIgetVoxVal(mri_tmp, x, yi, zi,0)) /* inferior voxel on */
            {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                DiagBreak() ;
              }
              nchanged++ ;
              MRIsetVoxVal(mri_tmp, x,y,z, 0, 255) ;
            }
            else
            {
              zi = mri_tmp->zi[z-1] ;
              if (MRIgetVoxVal(mri_tmp, x, yi,  zi, 0))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  DiagBreak() ;
                }
                nchanged++ ;
                MRIsetVoxVal(mri_tmp, x,y,z, 0, 255) ;
              }
            }
          }

          yi = mri_tmp->yi[y-1] ;
          if (MRIgetVoxVal(mri_tmp,x,yi,z,0)) /* check suprior */
          {
            yi = mri_tmp->yi[y+1] ;
            zi = mri_tmp->zi[z+1] ;
            if (MRIgetVoxVal(mri_tmp, x, yi, zi,0)) /* inferior voxel on */
            {
              if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              {
                DiagBreak() ;
              }
              nchanged++ ;
              MRIsetVoxVal(mri_tmp, x,y,z,0,255) ;
            }
            else
            {
              zi = mri_tmp->zi[z-1] ;
              if (MRIgetVoxVal(mri_tmp, x, yi,  zi,0))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                {
                  DiagBreak() ;
                }
                nchanged++ ;
                MRIsetVoxVal(mri_tmp, x,y,z, 0,  255) ;
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
          if (MRIgetVoxVal(mri_tmp, x, y, z,0))
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
            MRIsetVoxVal(mri_tmp, x, y, z, 0, 128) ;
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

        if (MRIgetVoxVal(mri_tmp, x, y, z,0) != 255) /* not temporal wm */
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
          label = MRIgetVoxVal(mri_labeled,x,yi,z,0) ;
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
               cma_label_to_name(MRIgetVoxVal(mri_labeled,x,yi, z,0)),
               cma_label_to_name(label)) ;
            }
            MRIsetVoxVal(mri_labeled, x, yi, z, 0,  label) ;
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
                     (MRIgetVoxVal(mri_labeled,x,yi, z,0)),
                     cma_label_to_name(label)) ;
                  }
                  MRIsetVoxVal(mri_labeled, x, yi, z, 0, label) ;
                }
              }
            }
          }//!GCA
#endif
          /* check for hippocampus inferior */
          yi = mri_tmp->yi[y+i] ;
          label = MRIgetVoxVal(mri_labeled,x,yi,z,0) ;
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
                 (MRIgetVoxVal(mri_labeled,x,yi, z,0)),
                 cma_label_to_name(label)) ;
              }
              MRIsetVoxVal(mri_labeled, x, yi, z, 0, label) ;
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
    if (MRIgetVoxVal(mri_tmp, x, y, z,0))
    {
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
      {
        printf("changing voxel (%d, %d, %d) from %s to %s\n",
               x, y, z, cma_label_to_name(MRIgetVoxVal(mri_labeled,x, y,z,0)),
               cma_label_to_name(gcas[i].label)) ;
        DiagBreak() ;
      }
      MRIsetVoxVal(mri_labeled, x, y, z, 0, gcas[i].label) ;
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

    GCAcomputeLabelMeansAndCovariances
      (gca, Unknown, &m_cov_un, &v_mean_un) ;
    GCAcomputeLabelMeansAndCovariances
      (gca, wm_label, &m_cov_wm, &v_mean_wm) ;
    GCAcomputeLabelMeansAndCovariances
      (gca, caudate_label, &m_cov_wm, &v_mean_caudate) ;
    GCAcomputeLabelMeansAndCovariances
      (gca, wmsa_label, &m_cov_wmsa, &v_mean_wmsa) ;
    m_inv_cov_wm = MatrixInverse(m_cov_wm, NULL) ;
    if (m_inv_cov_wm == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
                Progname, cma_label_to_name(wm_label), wm_label) ;
    m_inv_cov_un = MatrixInverse(m_cov_un, NULL) ;
    if (m_inv_cov_un == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
                Progname, cma_label_to_name(Unknown), Unknown) ;
    m_inv_cov_wmsa = MatrixInverse(m_cov_wmsa, NULL) ;
    if (m_inv_cov_wmsa == NULL)
      ErrorExit(ERROR_BADPARM,
                "%s: could not compute inverse covariance for %s (%d)",
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
          else if (nwm+nwmsa < .9*WSIZE*WSIZE*WSIZE) 
          {// somewhat arbitrary - the bulk of the nbhd
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
          pwmsa = compute_conditional_density
            (m_inv_cov_wmsa, v_mean_wmsa, v_vals) ;
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("         - pwm = %2.3e, pwmsa = %2.3e\n",
                   pwm, pwmsa) ;
          if (label == wm_label && pwmsa > pwm)
          {
            wm_dist = VectorDistance(v_mean_wm, v_vals) ;
            wmsa_dist = VectorDistance(v_mean_wmsa, v_vals) ;
            wm_mdist = MatrixMahalanobisDistance
              (v_mean_wm, m_inv_cov_wm, v_vals) ;
            wmsa_mdist = MatrixMahalanobisDistance
              (v_mean_wmsa, m_inv_cov_wmsa, v_vals) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              printf("         - wm_dist = %2.0f, wmsa_dist = %2.0f,"
                     " mdists = (%2.0f, %2.0f)\n",
                     wm_dist, wmsa_dist, wm_mdist, wmsa_mdist) ;
            if ((wm_dist > wmsa_dist) && (wm_mdist > wmsa_mdist))
            {
              VectorSubtract(v_vals, v_mean_wm, v_dif_label) ;
              VectorSubtract(v_vals, v_mean_wmsa, v_dif_wmsa) ;
              if (
                ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < 
                  fabs(VECTOR_ELT(v_dif_label,1))) &&
                 (fabs(VECTOR_ELT(v_dif_wmsa,2)) < 
                  fabs(VECTOR_ELT(v_dif_label,2))) &&
                 (fabs(VECTOR_ELT(v_dif_wmsa,3)) < 
                  fabs(VECTOR_ELT(v_dif_label,3)))) ||
                ((2*wmsa_dist < wm_dist) && (2*wmsa_mdist < wm_mdist)))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("GCAlabelWMandWMSAs 1: changing label from %s to %s\n",
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

            if (ngm+ncaudate+nwm+nwmsa+nunknown < .9*WSIZE*WSIZE*WSIZE)  
            {// somewhat arbitrary - the bulk of the nbhd
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
              if ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < 
                   fabs(VECTOR_ELT(v_dif_label,1))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,2)) < 
                   fabs(VECTOR_ELT(v_dif_label,2))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,3)) < 
                   fabs(VECTOR_ELT(v_dif_label,3))))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("GCAlabelWMandWMSAs 2: changing label from %s to %s\n",
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
              if (((fabs(VECTOR_ELT(v_dif_wmsa,1)) < 
                    fabs(VECTOR_ELT(v_dif_label,1))) &&
                   (fabs(VECTOR_ELT(v_dif_wmsa,2)) < 
                    fabs(VECTOR_ELT(v_dif_label,2))) &&
                   (fabs(VECTOR_ELT(v_dif_wmsa,3)) < 
                    fabs(VECTOR_ELT(v_dif_label,3)))) ||
                  (wmsa_dist*3 < wm_dist))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("GCAlabelWMandWMSAs 3: changing label from %s to %s\n",
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
              if ((fabs(VECTOR_ELT(v_dif_wmsa,1)) < 
                   fabs(VECTOR_ELT(v_dif_label,1))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,2)) < 
                   fabs(VECTOR_ELT(v_dif_label,2))) &&
                  (fabs(VECTOR_ELT(v_dif_wmsa,3)) < 
                   fabs(VECTOR_ELT(v_dif_label,3))))
              {
                if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  printf("GCAlabelWMandWMSAs 4: changing label from %s to %s\n",
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

#include "voxlist.h"

#define EXCLUDED(l)   (((l) == Unknown) || ((l) == CSF) || ((l) == Third_Ventricle))
static MRI *
MRIexpandVentricles(MRI *mri_labeled_src, 
                    MRI *mri_labeled,
                    GCA *gca, 
                    MRI *mri_inputs,
                    int num_expansions)
{
  MRI       *mri_vent, *mri_vent_orig ;
  double    mean, std, thresh ;
  VOXLIST   *vl = NULL ;
  int       nadded, x, y, z, label, i, n ;
  MRI       *mri_dilated = NULL, *mri_border = NULL ;
  float     label_means[MAX_CMA_LABELS], label_stds[MAX_CMA_LABELS], vdist, ldist, val ;

  if (mri_labeled == NULL)
  {
    mri_labeled = MRIcopy(mri_labeled_src, NULL) ;
  }

  for (label = 0 ; label <= MAX_CMA_LABEL ; label++)
  {
    GCAlabelMean(gca, label, &label_means[label]) ;
    GCAlabelVar(gca, label, &label_stds[label]) ;
    label_stds[label] = sqrt(label_stds[label]) ;
  }

  mri_vent = MRIclone(mri_labeled_src, NULL) ;
  MRIcopyLabel(mri_labeled_src, mri_vent, Left_Lateral_Ventricle) ;
  MRIcopyLabel(mri_labeled_src, mri_vent, Right_Lateral_Ventricle) ;
  MRIcopyLabel(mri_labeled_src, mri_vent, Left_Inf_Lat_Vent) ;
  MRIcopyLabel(mri_labeled_src, mri_vent, Right_Inf_Lat_Vent) ;
//  MRIcopyLabel(mri_labeled_src, mri_vent, Third_Ventricle) ;
  MRIbinarize(mri_vent, mri_vent, 1, 0, 1) ;
  mean = MRImeanAndStdInLabel(mri_inputs, mri_vent, 1, &std) ;
  mri_vent_orig = MRIcopy(mri_vent, NULL) ;
  thresh = mean+.5*std ;
  printf("expanding ventricles using %2.1f + 0.5 * %2.1f = %2.1f as threshold\n", mean, std, thresh) ;

  for (n = 0 ; n < num_expansions ; n++)
  {
    mri_dilated = MRIdilate(mri_vent, mri_dilated) ;
    mri_border = MRIsubtract(mri_dilated, mri_vent, mri_border) ;
    vl = VLSTcreate(mri_border, 1, 1, NULL, 0, 0) ;
    nadded = 0 ;
    for (i = 0 ; i < vl->nvox ; i++)
    {
      x = vl->xi[i] ;
      y = vl->yi[i] ;
      z = vl->zi[i] ;
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        DiagBreak() ;

      label = MRIgetVoxVal(mri_labeled, x, y, z, 0)  ;
      if (EXCLUDED(label))
	continue ;

      val =  MRIgetVoxVal(mri_inputs, x, y, z, 0) ;
      vdist = sqrt(SQR((val - mean) / std)) ;
      ldist = sqrt(SQR((val - label_means[label]) / label_stds[label])) ;
      if (vdist < 0.5*ldist)
      {
        label = MRIgetVoxVal(mri_labeled_src, x, y, z, 0) ;
        if (!EXCLUDED(label))
        {
          nadded++ ;
          MRIsetVoxVal(mri_vent, x, y, z, 0, 1) ;
        }
      }
    }
    MRIclear(mri_dilated) ;
    MRIclear(mri_border) ;
  }
  printf("adding %d ventricular labels\n", nadded) ;

  for (x = 0 ; x < mri_vent->width;  x++)
    for (y = 0 ; y < mri_vent->height;  y++)
      for (z = 0 ; z < mri_vent->depth;  z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          DiagBreak() ;

        if (MRIgetVoxVal(mri_vent, x, y, z, 0) > 0 && 
            MRIgetVoxVal(mri_vent_orig, x, y, z, 0)  == 0) 
        {// new ventricle label
          int lh, rh, third, lhi, rhi ;
          lh = MRIlabelsInNbhd(mri_labeled_src, x, y, z, 
                               num_expansions, Left_Lateral_Ventricle) ;
          lhi = MRIlabelsInNbhd(mri_labeled_src, x, y, z, 
                               num_expansions, Left_Inf_Lat_Vent) ;
          rh = MRIlabelsInNbhd(mri_labeled_src, x, y, z, 
                               num_expansions, Right_Lateral_Ventricle) ;
          rhi = MRIlabelsInNbhd(mri_labeled_src, x, y, z, 
                               num_expansions, Right_Inf_Lat_Vent) ;
          third = MRIlabelsInNbhd(mri_labeled_src, x, y, z, 
                                  num_expansions, Third_Ventricle) ;
          if (lh > rh && lh > third && lh > lhi && lh > rhi)
          {
            label = Left_Lateral_Ventricle ;
          }
          else if (rh > third && rh > lhi && rh > rhi)
          {
            label = Right_Lateral_Ventricle ;
          }
          else if (third > lhi && third > rhi)
          {
            label = Third_Ventricle ;
          }
	  else if (lhi > rhi)
          {
            label = Left_Inf_Lat_Vent ;
          }
	  else
          {
            label = Right_Inf_Lat_Vent ;
          }
          MRIsetVoxVal(mri_labeled, x, y, z, 0, label) ;
        }
      }

  MRIfree(&mri_vent) ;
  MRIfree(&mri_vent_orig) ;
  return(mri_labeled) ;
}

static MRI *
replace_cortex_far_from_surface_with_wmsa(MRI *mri_inputs,
                                          MRI *mri_src,
                                          MRI *mri_dst,
                                          MRI_SURFACE *mris_lh,
                                          MRI_SURFACE *mris_rh)
{
  MRI     *mri_lh_dist, *mri_rh_dist ;
  int     x, y, z, label, f, nchanged = 0 ;
  VECTOR *v_wm_mean, *v_wmsa_mean, *v_gm_mean, *v_tmp ;
  float  wm_dist, wmsa_dist, gm_dist ;


  mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_lh_dist = MRIcloneDifferentType(mri_dst, MRI_FLOAT) ;
  MRIScomputeDistanceToSurface(mris_lh, mri_lh_dist, mri_dst->xsize) ;
  mri_rh_dist = MRIcloneDifferentType(mri_dst, MRI_FLOAT) ;
  MRIScomputeDistanceToSurface(mris_rh, mri_rh_dist, mri_dst->xsize) ;

  v_tmp = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Left_Cerebral_Cortex)  ;
  v_gm_mean = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Right_Cerebral_Cortex)  ;
  VectorAdd(v_tmp, v_gm_mean, v_gm_mean) ;
  VectorScalarMul(v_gm_mean, 0.5, v_gm_mean) ;
  VectorFree(&v_tmp) ;

  v_tmp = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Left_Cerebral_White_Matter)  ;
  v_wm_mean = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Right_Cerebral_White_Matter)  ;
  VectorAdd(v_tmp, v_wm_mean, v_wm_mean) ;
  VectorScalarMul(v_wm_mean, 0.5, v_wm_mean) ;
  VectorFree(&v_tmp) ;

  v_tmp = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Left_WM_hypointensities)  ;
  v_wmsa_mean = MRImeanInLabelMultispectral
    (mri_inputs, mri_dst, Right_WM_hypointensities)  ;
  VectorAdd(v_tmp, v_wmsa_mean, v_wmsa_mean) ;
  VectorScalarMul(v_wmsa_mean, 0.5, v_wmsa_mean) ;
//  VectorFree(&v_tmp) ; don't free it - will use it later

  for (x = 0 ; x < mri_dst->width ; x++)
    for (y = 0 ; y < mri_dst->height ; y++)
      for (z = 0 ; z < mri_dst->depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
        if (IS_CORTEX(label) == 0)
        {
          continue ;
        }
        if (MRIgetVoxVal(mri_lh_dist, x, y, z, 0) < -1.5 || 
            MRIgetVoxVal(mri_rh_dist, x, y, z, 0) < -1.5)
        {
          for (f = 0 ; f < mri_inputs->nframes; f++)
          {
            VECTOR_ELT(v_tmp, f+1) = MRIgetVoxVal(mri_inputs, x, y, z, f) ;
          }

          wm_dist = VectorDistance(v_wm_mean, v_tmp) ;
          wmsa_dist = VectorDistance(v_wmsa_mean, v_tmp) ;
          gm_dist = VectorDistance(v_gm_mean, v_tmp) ;
          if (wmsa_dist < gm_dist && wmsa_dist < wm_dist)
          {
            nchanged++ ;
            if (label == Left_Cerebral_Cortex)
            {
              label = Left_WM_hypointensities ;
            }
            else
            {
              label = Right_WM_hypointensities ;
            }
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
              printf("voxel (%d, %d, %d) is interior to surface "
                     "and changed to %s (int dists %2.0f, %2.0f, %2.0f,"
                     " surf dists %2.1f, %2.1f)\n",
                     x, y, z, cma_label_to_name(label),
                     wm_dist, gm_dist, wmsa_dist,
                     MRIgetVoxVal(mri_lh_dist, x, y, z, 0),
                     MRIgetVoxVal(mri_rh_dist, x, y, z, 0)) ;

            MRIsetVoxVal(mri_dst, x, y, z, 0, label) ;
          }
          else if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
            printf("voxel (%d, %d, %d) is interior to surface "
                   "but not changed (int dists %2.0f, %2.0f, %2.0f, "
                   "surf dists %2.1f, %2.1f)\n",
                   x, y, z,
                   wm_dist, gm_dist, wmsa_dist,
                   MRIgetVoxVal(mri_lh_dist, x, y, z, 0),
                   MRIgetVoxVal(mri_rh_dist, x, y, z, 0)) ;
        }
        else if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          printf("voxel (%d, %d, %d) is not interior to surface "
                 "and not changed (surf dists %2.1f, %2.1f)\n",
                 x, y, z,
                 MRIgetVoxVal(mri_lh_dist, x, y, z, 0),
                 MRIgetVoxVal(mri_rh_dist, x, y, z, 0)) ;
      }

  printf("%d voxels interior to surface changed from cortex to WMSA\n",
         nchanged) ;
  MRIfree(&mri_lh_dist) ;
  MRIfree(&mri_rh_dist) ;
  return(mri_dst) ;
}

/*
  \fn int MRItoUCHAR(MRI **pmri)
  \brief Changes the mri to uchar if it is not uchar already and the
  range is between 0 and 255 inclusive. This is mainly used by the
  output of mri_ca_label to change the output to uchar since the inner
  workings of the labeling were changed to use int. With this, the
  aseg output will be uchar and so not break anything downstream.
 */
int MRItoUCHAR(MRI **pmri)
{
  MRI *mri, *mri2;
  int c,r,s,f;
  double v, vmin, vmax;

  mri = *pmri;

  if(mri->type == MRI_UCHAR)
  {
    return(0);
  }

  vmax = MRIgetVoxVal(mri,0,0,0,0);
  vmin = vmax;
  for(c=0; c < mri->width; c++)
  {
    for(r=0; r < mri->height; r++)
    {
      for(s=0; s < mri->depth; s++)
      {
        for(f=0; f < mri->nframes; f++)
        {
          v = MRIgetVoxVal(mri,c,r,s,f);
          if(v < vmin)
          {
            vmin = v;
          }
          if(v > vmax)
          {
            vmax = v;
          }
        }
      }
    }
  }

  printf("MRItoUCHAR: min=%g, max=%g\n",vmin,vmax);

  if(vmin < 0 || vmax > 255)
  {
    printf("MRItoUCHAR: range too large, not changing type\n");
    return(0);
  }

  printf("MRItoUCHAR: converting to UCHAR\n");

  mri2 = MRIcloneBySpace(mri, MRI_UCHAR, mri->nframes);

  for(c=0; c < mri->width; c++)
  {
    for(r=0; r < mri->height; r++)
    {
      for(s=0; s < mri->depth; s++)
      {
        for(f=0; f < mri->nframes; f++)
        {
          v = MRIgetVoxVal(mri,c,r,s,f);
          MRIsetVoxVal(mri2,c,r,s,f,nint(v));
        }
      }
    }
  }

  MRIfree(&mri);
  *pmri = mri2;
  return(1);
}

static MRI *
fix_putamen(GCA *gca,
            MRI *mri_inputs,
            MRI *mri_imp,
            TRANSFORM *transform,
            MRI *mri_src_labeled,
            MRI *mri_dst_labeled,
            double prior_thresh)
{
  int x, y, z, nchanged, label, n, left, right, above, below;
  int gm, wm, iter, total_changed ;
  double     pwm, pgm ;
  GCA_PRIOR *gcap ;
  // This was never tested
  // This function is called above with prior_thresh < 0, so never run

  if (mri_dst_labeled == NULL)
  {
    mri_dst_labeled = MRIcopy(mri_src_labeled, NULL) ;
  }
  if (prior_thresh < 0)
  {
    return(mri_dst_labeled) ;
  }

  printf("fix_putamen()\n");
  for (iter = total_changed = 0 ; iter < 5 ; iter++)
  {
    for (nchanged = x = 0 ; x < mri_imp->width ; x++)
      for (y = 0 ; y < mri_imp->height ; y++)
        for (z = 0 ; z < mri_imp->depth ; z++)
        {
          if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_imp, x, y, z, 0) > 0)
          {
            gcap = getGCAP(gca, mri_inputs, transform, x, y, z) ;
            for (n = 0 ; n < gcap->nlabels ; n++)
              if (IS_PUTAMEN(gcap->labels[n]))
              {
                break ;
              }
            if (n < gcap->nlabels)   // found putamen
            {
              label = MRIgetVoxVal(mri_src_labeled, x, y, z, 0) ;
              if (IS_PUTAMEN(label) && gcap->priors[n] < prior_thresh)
              {
                left = MRIgetVoxVal(mri_src_labeled, x+1, y, z, 0) ;
                right = MRIgetVoxVal(mri_src_labeled, x-1, y, z, 0) ;
                above = MRIgetVoxVal(mri_src_labeled, x, y-1, z, 0) ;
                below = MRIgetVoxVal(mri_src_labeled, x, y+1, z, 0) ;
                if (((IS_CORTEX(left) || IS_WM(left)) &&
                     (IS_CORTEX(right) || IS_WM(right))) ||
                    (((IS_CORTEX(above) || IS_WM(above)) &&
                      (IS_CORTEX(below) || IS_WM(below)))))
                {
                  nchanged++ ;
                  if  (label == Left_Putamen)
                  {
                    gm = Left_Cerebral_Cortex ;
                    wm = Left_Cerebral_White_Matter ;
                  }
                  else
                  {
                    gm = Right_Cerebral_Cortex ;
                    wm = Right_Cerebral_White_Matter ;
                  }
                  pwm = GCAlabelProbability
                    (mri_inputs, gca, transform, 
                     (float)x, (float)y, (float)z,wm) ;
                  pgm = GCAlabelProbability
                    (mri_inputs, gca, transform,
                     (float)x, (float)y, (float)z,gm) ;
                  if (pwm > pgm) // change it to white matter
                  {
                    label = wm ;
                  }
                  else           // change it to gray matter
                  {
                    label = gm ;
                  }
                  MRIsetVoxVal(mri_dst_labeled, x, y, z, 0, label) ;
                  if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
                  {
                    printf("changing putamen label at (%d, %d, %d) "
                           "to %s (%d)\n",
                           x, y, z, cma_label_to_name(label), label) ;
                  }
                }
              }
              else   // consider changing it to putamen from something else
              {
              }
            }
          }
        }
    total_changed += nchanged ;
    if (nchanged == 0)
    {
      break ;
    }
  }

  printf("%d putamen labels changed to cortex/white matter\n", total_changed) ;
  return(mri_dst_labeled) ;
}


VOXEL_LIST *
find_unlikely_voxels_in_region(GCA *gca,
                               TRANSFORM *transform,
                               MRI *mri_labeled,
                               double prior_thresh,
                               MRI *mri_prior_labels,
                               MRI *mri_priors,
                               int x, int y, int z,
                               int whalf)
{
  int         xi, yi, zi, xk, yk, zk, nvox, label, prior_label ;
  double      prior ;
  GCA_PRIOR   *gcap ;
  VOXEL_LIST  *vl ;

  for (nvox = 0, xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri_labeled->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri_labeled->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri_labeled->zi[z+zk] ;
        label = MRIgetVoxVal(mri_labeled, xi, yi, zi, 0) ;
	if (IS_WMSA(label))  // applying spatial priors to WMSAs doesn't work
	  continue ;
        prior_label = MRIgetVoxVal(mri_prior_labels, xi, yi, zi, 0) ;
        if (label == prior_label)
        {
          continue ;
        }
        gcap = getGCAP(gca, mri_labeled, transform, xi, yi, zi) ;
        prior = getPrior(gcap, label) ;
        if (prior < prior_thresh)
        {
          nvox++ ;
        }
      }
    }
  }
  if (nvox == 0)
  {
    return(NULL) ;
  }
  vl = VLSTalloc(nvox) ;
  vl->nvox = 0 ;
  vl->mri = mri_labeled ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri_labeled->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri_labeled->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri_labeled->zi[z+zk] ;
        if (xi == Gx && yi == Gy && zi == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_labeled, xi, yi, zi, 0) ;
	if (IS_WMSA(label))  // applying spatial priors to WMSAs doesn't work
	  continue ;
        prior_label = MRIgetVoxVal(mri_prior_labels, xi, yi, zi, 0) ;
        if (label == prior_label)
        {
          continue ;
        }
        gcap = getGCAP(gca, mri_labeled, transform, xi, yi, zi) ;
        prior = getPrior(gcap, label) ;
        if (prior < prior_thresh)
        {
          VLSTadd(vl, xi, yi, zi, xi, yi, zi) ;
        }
      }
    }
  }
  VLSTsample(vl, mri_labeled) ;

  return(vl) ;
}

static int
change_unlikely_voxels(GCA *gca,
                       MRI *mri_dst_label, 
                       MRI *mri_inputs,
                       TRANSFORM *transform,
                       VOXEL_LIST *vl,
                       int label_to_change, 
                       MRI *mri_prior_labels)
{
  int       i, max_label, x, y, z, nchanged ;
#if 0
  int       n ;
  double    p, max_p ;
  GCA_PRIOR *gcap ;
#endif

  for (nchanged = i = 0 ; i < vl->nvox ; i++)
  {
    if (nint(vl->vsrc[i]) != label_to_change)
    {
      continue ;
    }
    nchanged++ ;
    x = vl->xi[i] ;
    y = vl->yi[i] ;
    z = vl->zi[i] ;
    max_label = MRIgetVoxVal(mri_prior_labels, x, y, z, 0) ;
    MRIsetVoxVal(mri_dst_label, x, y, z, 0, max_label) ; // use prior label
    if (0 && x == Ggca_x && y == Ggca_y && z == Ggca_z)
    {
      printf("change_unlikely_voxels 1: changing label at (%d, %d, %d) from %s to %s\n",
             Ggca_x, Ggca_y, Ggca_z, 
             cma_label_to_name(label_to_change),
             cma_label_to_name(max_label));
    }
#if 0
    gcap = getGCAP(gca, mri_inputs, transform, 
                   vl->xi[i], vl->yi[i], vl->zi[i]) ;
    max_p = 0.0 ;
    max_label = 0 ;
    for (n = 0 ; n < gcap->nlabels ; n++)
    {
      if (x == Gx && y == Gy && z == Gz)
      {
        DiagBreak() ;
      }
      if (gcap->labels[n] == label_to_change)   
      {// don't let it be the one we are changing
        continue ;
      }
      p = GCAlabelProbability(mri_inputs, gca, transform, 
                              (float)x, (float)y, (float)z, gcap->labels[n]) ;
      if (p > max_p)
      {
        max_label = gcap->labels[n] ;
        max_p = p ;
      }
    }
    if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
    {
      printf("change_unlikely_voxels 2: changing label at (%d, %d, %d) from %s to %s\n",
             Ggca_x, Ggca_y, Ggca_z,
             cma_label_to_name(label_to_change),
             cma_label_to_name(max_label));
    }
    MRIsetVoxVal(mri_dst_label, x, y, z, 0, max_label) ;
#endif
  }

  return(nchanged) ;
}

static MRI *
GCArelabelUnlikely(GCA *gca,
                   MRI *mri_inputs,
                   TRANSFORM *transform, 
                   MRI *mri_src_labeled, 
                   MRI *mri_dst_labeled,
                   MRI *mri_independent_posterior,
                   double prior_thresh,
                   int whalf)
{
  int         x, y, z, nchanged, total_changed, i, nindices, 
              label, index, w  ;
  short       *x_indices, *y_indices, *z_indices ;
  MRI        *mri_prior_labels, *mri_priors, *mri_unchanged, *mri_tmp ;
  double      prior  ;

  mri_prior_labels = MRIclone(mri_dst_labeled, NULL) ;
  mri_priors = MRIcloneDifferentType(mri_dst_labeled, MRI_FLOAT) ;
  mri_unchanged = MRIcloneDifferentType(mri_dst_labeled, MRI_UCHAR) ;
  mri_tmp = MRIclone(mri_unchanged, NULL) ;

  if (mri_dst_labeled == NULL)
  {
    mri_dst_labeled = MRIcopy(mri_src_labeled, NULL) ;
  }

  for (x = 0 ; x < mri_inputs->width ; x++)
    for (y = 0 ; y < mri_inputs->height ; y++)
      for (z = 0 ; z < mri_inputs->depth ; z++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          DiagBreak() ;
        }
        label = GCAgetMaxPriorLabelAtVoxel(gca, mri_dst_labeled, 
                                           x, y, z, transform, &prior) ;
        MRIsetVoxVal(mri_prior_labels, x, y, z, 0, label) ;
        MRIsetVoxVal(mri_priors, x, y, z, 0, prior) ;
      }

  nindices = mri_inputs->width * mri_inputs->height * mri_inputs->depth ;
  x_indices = (short *)calloc(nindices, sizeof(short)) ;
  y_indices = (short *)calloc(nindices, sizeof(short)) ;
  z_indices = (short *)calloc(nindices, sizeof(short)) ;
  for (i = total_changed = 0 ; i < 5 ; i++)
  {
    nchanged = 0 ;
    MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices, z_indices) ;
#if 0 //def HAVE_OPENMP   doesn't work
#pragma omp parallel for if_ROMP(experimental) firstprivate(mri_independent_posterior, mri_inputs, gca, mri_prior_labels, whalf, prior_thresh, Ggca_x, Ggca_y, Ggca_z) shared(mri_unchanged, mri_priors,transform, x_indices,y_indices,z_indices, mri_dst_labeled) reduction(+:nchanged)
#endif
    for (index = 0 ; index < nindices ; index++)
    {
      double      posterior_before, posterior_after  ;
      int         label, x, y, z, nchanged_tmp, old_label = 0 ;
      VOXEL_LIST  *vl ;

      x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
      if ((int)MRIgetVoxVal(mri_unchanged, x, y, z, 0) == 1)
	continue ;  // this nbhd not changed from last call

      MRIsetVoxVal(mri_unchanged, x, y, z, 0, 1) ;
      label = MRIgetVoxVal(mri_dst_labeled, x, y, z, 0) ;
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        DiagBreak() ;

      if (IS_WMSA(label))
	continue ;

      // don't process it if it the highest prior or the highest posterior
      if ((label == (int)MRIgetVoxVal(mri_prior_labels,x,y,z,0) &&
	   (mri_independent_posterior && 
	    label == (int)MRIgetVoxVal(mri_independent_posterior,x,y,z,0))))
      {
        continue ;
      }

#if 0
      if ((GCAisPossible(gca, mri_inputs, Left_Putamen,
                         transform, x, y, z, 0) == 0) &&
          (GCAisPossible(gca, mri_inputs, Left_Putamen,
                         transform, x, y, z, 0) == 0))
      {
        continue ;
      }
      if (!IS_PUTAMEN(label))   // disable everything but putamen for now
      {
        continue ;
      }
#endif
      posterior_before = 
        GCAwindowPosteriorLogProbability
        (gca, mri_dst_labeled, mri_inputs, transform, x, y, z,  whalf)  ;

      vl = find_unlikely_voxels_in_region
        (gca, transform, mri_dst_labeled, prior_thresh, 
         mri_prior_labels, mri_priors, x, y, z, whalf) ;

      if (vl == NULL)
        continue;
      if (Ggca_x >= 0)
	old_label = MRIgetVoxVal(mri_src_labeled, Ggca_x, Ggca_y, Ggca_z, 0) ;
	
      if (Ggca_x >= 0 && VLSTinList(vl, Ggca_x, Ggca_y, Ggca_z))
      {
        MRI *mri ;
        mri = VLSTtoMri(vl, NULL) ;
        MRIwrite(mri, "vl.mgz") ;
        MRIfree(&mri) ;
        DiagBreak() ;
      }

      nchanged_tmp = 
        change_unlikely_voxels
        (gca, mri_dst_labeled, mri_inputs, transform, 
         vl, label,mri_prior_labels) ;
      posterior_after = 
        GCAwindowPosteriorLogProbability
        (gca, mri_dst_labeled, mri_inputs, transform, x, y, z,  whalf)  ;
      if (posterior_after <= posterior_before)
      {
        VLSTvsrcToMri(vl, mri_dst_labeled) ;  // undo label change
      }
      else
      {
	if (Ggca_x >= 0)
	{
	  int max_label = MRIgetVoxVal(mri_dst_labeled, Ggca_x, Ggca_y, Ggca_z, 0) ;
	  if  (old_label != max_label && VLSTinList(vl, Ggca_x, Ggca_y, Ggca_z))
	  {
	    
	    printf("iter = %d: change_unlikely_voxels: changing label at (%d, %d, %d) from %s to %s\n",
		   i, Ggca_x, Ggca_y, Ggca_z, 
		   cma_label_to_name(old_label),
		   cma_label_to_name(max_label));
	  }
	}
	MRIsetVoxVal(mri_unchanged, x, y, z, 0, 0) ;
        nchanged += nchanged_tmp ;
        DiagBreak() ;
      }
      VLSTfree(&vl) ;
    }
    total_changed += nchanged ;
    printf("%d voxels changed in iteration %d of "
           "unlikely voxel relabeling\n", nchanged, i) ;
    if (!nchanged)
      break ;
    for (w = 0 ; w < whalf ; w++)
    {
      MRIerode(mri_unchanged, mri_tmp) ;
      MRIcopy(mri_tmp, mri_unchanged) ;
    }
  }

  MRIfree(&mri_prior_labels) ; MRIfree(&mri_priors) ; MRIfree(&mri_unchanged) ; MRIfree(&mri_tmp) ;
  return(mri_dst_labeled) ;
}

static MRI *
gcaCheckForFCDs(MRI *mri_src, MRI *mri_dst, GCA *gca, TRANSFORM *transform, MRI *mri_inputs) 
{
  int        x, y, z, label, out_label, nchanged, iter = 0 ;
  float      left_means[2], right_means[2], global_gm_mean, global_wm_mean ;

  GCAlabelMean(gca, Left_Cerebral_Cortex, left_means) ;
  GCAlabelMean(gca, Right_Cerebral_Cortex, right_means) ;
  global_gm_mean = 0.5 * (left_means[0] + right_means[0]) ;
  GCAlabelMean(gca, Left_Cerebral_White_Matter, left_means) ;
  GCAlabelMean(gca, Right_Cerebral_White_Matter, right_means) ;
  global_wm_mean = 0.5 * (left_means[0] + right_means[0]) ;

  printf("checking for FCDS - global means: WM=%2.0f, GM=%2.0f\n", global_wm_mean, global_gm_mean);
  mri_dst = MRIcopy(mri_src, mri_dst) ;

  do
  {
    nchanged = 0 ;
    for (x = 0 ; x < mri_src->width ; x++)
      for (y = 0 ; y < mri_src->height ; y++)
	for (z = 0 ; z < mri_src->depth ; z++)
	{
	  out_label = label = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
	  if (Ggca_x == x && Ggca_y == y && Ggca_z == z)
	    printf("label at (%d, %d, %d) = %s\n", x, y, z, cma_label_to_name(label)) ;
	  if (IS_HYPO(label))
	  {
	    out_label = Left_WM_hypointensities ? Left_non_WM_hypointensities : Right_non_WM_hypointensities;
	  }
	  else if (IS_WM(label))
	  {
	    GCA_PRIOR *gcap ;
	    GC1D      *gm_gc, *wm_gc ;
	    float     wm_prior, gm_prior, gm_dist, wm_dist, val, thresh ;
	    
	    if (MRIcountNbhdLabels(mri_src, x, y, z, Left_WM_hypointensities) > 0 ||
		MRIcountNbhdLabels(mri_src, x, y, z, Right_WM_hypointensities) > 0 ||
		MRIcountNbhdLabels(mri_src, x, y, z, WM_hypointensities) > 0 ||
		MRIcountNbhdLabels(mri_dst, x, y, z, Left_non_WM_hypointensities) > 0 ||
		MRIcountNbhdLabels(mri_dst, x, y, z, Right_non_WM_hypointensities) > 0)
	      
	      thresh = 1.3 ;
	    
	    val = MRIgetVoxVal(mri_inputs, x, y, z, 0) ;
	    gcap = getGCAP(gca, mri_src, transform,  x,  y, z) ;
	    wm_prior = getPrior(gcap, label) ;
	    if (label == Left_Cerebral_White_Matter)
	    {
	      gm_prior = getPrior(gcap, Left_Cerebral_Cortex) ;
	      gm_gc = GCAfindSourceGC(gca, mri_src, transform, x, y, z, Left_Cerebral_Cortex) ;
	    }
	    else
	    {
	      gm_prior = getPrior(gcap, Right_Cerebral_Cortex) ;
	      gm_gc = GCAfindSourceGC(gca, mri_src, transform, x, y, z, Left_Cerebral_Cortex) ;
	    }
	    wm_gc = GCAfindSourceGC(gca, mri_src, transform, x, y, z, Left_Cerebral_Cortex) ;
	    gm_dist = gm_gc ? fabs(val - gm_gc->means[0]) : fabs(val-global_gm_mean);
	    wm_dist = wm_gc ? fabs(val - wm_gc->means[0]) : fabs(val-global_wm_mean);
	    if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
	      printf("%d, %d, %d) = %2.0f, WM dist = %2.1f, prior = %2.2f, GM: dist = %2.1f,  prior = %2.2f, thresh %2.1f\n",
		     x, y, z, val, wm_dist, wm_prior, gm_dist, gm_prior, thresh) ;
	    if (wm_prior > gm_prior && thresh*gm_dist < wm_dist)
	    {
	      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
		printf(" CHANGING LABEL\n") ;
	      out_label = (label == Left_Cerebral_White_Matter) 
		? Left_non_WM_hypointensities : Right_non_WM_hypointensities ;
	    }
	  }
	  if (label != out_label && (x == Ggca_x && y == Ggca_y && z == Ggca_z))
	    printf("(%d %d %d): changing from %s to %s\n", x,y,z,cma_label_to_name(label),
		   cma_label_to_name(out_label)) ;
	  if (label != out_label)
	  {
	    nchanged++ ;
	    MRIsetVoxVal(mri_dst, x, y, z, 0, out_label);
	  }
	}
    printf("iter %d: %d changed\n", iter, nchanged) ;
  } while (nchanged > 0) ;
  return(mri_dst) ;
}

/*!
  \fn MRI *insert_wm_bet_putctx(MRI *seg, int topo, char *psfile, MRI *out)
  \brief Finds all putamen voxels with a bordering cortical voxel and
  changes them to WM. The rational is that anatomically putamen and
  cortex always have WM between them. Putamen is replaced with WM
  because it is the one usually overlabed (usually claustrum gets
  labeled as putamen). This is a quick fix as we should really be
  labeling claustrum.
*/
MRI *insert_wm_bet_putctx(MRI *seg, int topo, const char *psfile, MRI *out)
{
  printf("insert_wm_bet_putctx() topo=%d\n",topo);  fflush(stdout);

  if(out == NULL){
    out = MRIalloc(seg->width,seg->height,seg->depth,seg->type);
    if(out==NULL) return(NULL);
  }
  MRIcopyHeader(seg,out);
  MRIcopyPulseParameters(seg,out);

  MATRIX *vox2ras, *crs, *ras;
  fsPointSet ps;
  if(psfile != NULL && strcmp(psfile,"nofile")!=0){
    vox2ras = MRIxfmCRS2XYZ(seg,0);
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
    ras = NULL;
  }

  int c, r, s, nchanged=0;
  for(c = 0 ; c < seg->width ; c++) {
    for(r = 0 ; r < seg->height ; r++){
      for(s = 0 ; s < seg->depth ; s++) {
	int label = MRIgetVoxVal(seg, c, r, s, 0) ;
	if(!IS_PUTAMEN(label)) continue;
	int hit = 0;
	for(int dc = -1; dc < 2; dc++){
	  for(int dr = -1; dr < 2; dr++){
	    for(int ds = -1; ds < 2; ds++){
	      int dsum = abs(dc)+abs(dr)+abs(ds);
	      if(dsum > topo) continue;
	      int adjlabel = MRIgetVoxVal(seg,c+dc,r+dr,s+ds,0);
	      if(IS_CORTEX(adjlabel)) {
		hit++;
		if(c == Gx && r == Gy && s == Gz){
		  printf("%3d %3d %3d %d %d  %d %d %d\n",c,r,s,label,hit,c+dc,r+dr,s+ds); 
		  fflush(stdout);
		}
	      }
	    }
	  }
	}
	if(hit==0) continue; // cortex vox does not border putamen
	int wmlabel;
	if(label == Left_Putamen) wmlabel  = Left_Cerebral_White_Matter ;
	else                      wmlabel  = Right_Cerebral_White_Matter ;
	MRIsetVoxVal(out,c,r,s,0,wmlabel) ;
	nchanged++;
	if(psfile != NULL && strcmp(psfile,"nofile")!=0){
	  // Make a point set of the putamen voxels that were changed
	  crs->rptr[1][1] = c;
	  crs->rptr[2][1] = r;
	  crs->rptr[3][1] = s;
	  ras = MatrixMultiplyD(vox2ras,crs,ras);
	  ps.add((double)ras->rptr[1][1],(double)ras->rptr[2][1],(double)ras->rptr[3][1],0);
	}
      }//s
    }//r
  }//c

  if(psfile != NULL && strcmp(psfile,"nofile")!=0){
    ps.save(psfile);
    MatrixFree(&vox2ras);
    MatrixFree(&crs);
    MatrixFree(&ras);
  }

  printf("#INSERT_WM_BET_PUTCTX %d putamen voxels changed to WM\n", nchanged);
  fflush(stdout);
  return(out);
}


int MRIinsertFromSeg(MRI *mri_labeled, MRI *InsertFromSeg, std::vector<int> InsertFromSegIndices, int ZeroInsertIndex)
{
  int nchanged=0;
  for(int c=0; c < mri_labeled->width; c++){
    for(int r=0; r < mri_labeled->height; r++){
      for(int s=0; s < mri_labeled->depth; s++){
	// Get the indices of this voxel for both segs
	int i1 = MRIgetVoxVal(mri_labeled,c,r,s,0);
	int i2 = MRIgetVoxVal(InsertFromSeg,c,r,s,0);

	if(i1 == i2) continue; // both the same so keep as is

	// Check whether i1, i2 are in the list
	int i1InList=0, i2InList=0;
	for(int n=0; n < InsertFromSegIndices.size(); n++){
	  if(i1==InsertFromSegIndices[n]){
	    i1InList = 1;
	    break;
	  }
	}
	for(int n=0; n < InsertFromSegIndices.size(); n++){
	  if(i2==InsertFromSegIndices[n]){
	    i2InList = 1;
	    break;
	  }
	}
	// Neither are in the list, so keep as is
	if(!i1InList && !i2InList) continue; 

	// At this point, at least one of them is in the list, so
	// replace it. This will change things on the boundary. Eg, if
	// i1 is cblum wm and i2 is brainstem, then i1 will become
	// brainstem eventhough brainstem is not in the list. This is
	// probably ok for the main application (replacing aseg cblum
	// with synthseg cblum), but I worry a little bit about
	// causing some artifacts near the boundaries. zeros voxels
	// that are cblum in the input but CSF (24) in the output. The
	// ZeroInsertIndex is a hack; when applying this to cblum,
	// setting ZeroInsertIndex=24 will prevent a bunch of random
	// CSF (24) voxels in the output
	if(i2 == ZeroInsertIndex) i2 = 0;
	MRIsetVoxVal(mri_labeled,c,r,s,0,i2);
	nchanged++;
	
      }
    }
  }
  printf("MRIinsertFromSeg() changed %d\n",nchanged);
  return(nchanged);
}
