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
#include "gca.h"
#include "transform.h"
#include "cma.h"
#include "histo.h"
#include "mrinorm.h"

static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;

static int distance_to_label(MRI *mri_labeled, int label, int x, 
                             int y, int z, int dx, int dy, 
                             int dz, int max_dist) ;
static int preprocess(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta, 
                      MRI *mri_fixed) ;
static int edit_hippocampus(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta, 
                      MRI *mri_fixed) ;
static int edit_amygdala(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta, 
                         MRI *mri_fixed) ;
static int MRIcountNbhdLabels(MRI *mri, int x, int y, int z, int label) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static char *tissue_parms_fname = NULL ;
static char *mask_volume_fname = NULL ;

char *Progname ;
static void usage_exit(int code) ;

static int novar = 0 ;
static char *renormalization_fname = NULL ;
static int renormalize = 0 ;
static int filter = 0 ;
#if 0
static float thresh = 0.5 ;
#endif
static int read_flag = 0 ;
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
static MRI *insert_wm_segmentation(MRI *mri_labeled, MRI *mri_wm, 
                                  int parcellation_type, int fixed_flag,
                                   GCA *gca, LTA *lta) ;

extern char *gca_write_fname ;
extern int gca_write_iterations ;

static int expand_flag = TRUE ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  char   *in_fname, *out_fname,  *gca_fname, *xform_fname ;
  MRI    *mri_in, *mri_labeled, *mri_fixed = NULL ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  GCA     *gca ;
  LTA     *lta ;

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
      
    GCAhisto(gca, 100, &counts) ;

    max_i = 0 ;
    for (i = 1 ; i < 100 ; i++)
      if (counts[i] > 0)
        max_i = i ;

    for (i = 1 ; i <= max_i ; i++)
      printf("%d %d\n", i, counts[i]) ;
    exit(1) ;
  }

  if (argc < 5)
    usage_exit(1) ;

  in_fname = argv[1] ;
  xform_fname = argv[2];
  gca_fname = argv[3] ;
  out_fname = argv[4] ;

  
  printf("reading input volume from %s...\n", in_fname) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
              Progname, in_fname) ;

  if (filter)
  {
    MRI *mri_dir, /**mri_grad,*/ *mri_kernel, *mri_smooth ;

    mri_kernel = MRIgaussian1d(1, 15) ;
    mri_smooth = MRIconvolveGaussian(mri_in, NULL, mri_kernel) ;
    MRIfree(&mri_kernel) ;
#if 0
    mri_grad = MRIsobel(mri_smooth, NULL, NULL) ;
    mri_dir = MRIdirectionMapUchar(mri_grad, NULL, 5) ;
    MRIfree(&mri_grad) ; 
    MRIwrite(mri_dir, "dir.mgh") ;
#else
    mri_dir = MRIgradientDir2ndDerivative(mri_in, NULL, 5) ;
    MRIscaleAndMultiply(mri_in, 128.0, mri_dir, mri_in) ;
    MRIwrite(mri_dir, "lap.mgh") ;
    MRIwrite(mri_in, "filtered.mgh") ;
#endif
    MRIfree(&mri_dir) ; MRIfree(&mri_smooth) ;
    exit(1) ;
  }

  /*  fprintf(stderr, "mri_in read: xform %s\n", mri_in->transform_fname) ;*/
  printf("reading classifier array from %s...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier array from %s",
              Progname, gca_fname) ;

  if (alpha > 0)
    mri_in->flip_angle = alpha ;
  if (TR > 0)
    mri_in->tr = TR ;
  if (TE > 0)
    mri_in->te = TE ;

  if (example_T1)
  {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE,"%s: could not read example segmentation from %s",
                Progname, example_segmentation) ;
    mri_T1 = MRIread(example_T1) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE,"%s: could not read example T1 from %s",
                Progname, example_T1) ;
    printf("scaling atlas intensities using specified examples...\n") ;
    MRIeraseBorderPlanes(mri_seg) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ; MRIfree(&mri_T1) ;
  }

  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;

  if (mri_fname)
  {
    GCAbuildMostLikelyVolume(gca, mri_in) ;
    MRIwrite(mri_in, mri_fname) ;
    exit(0) ;
  }
  if (novar)
    GCAunifyVariance(gca) ;

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
      labels[i] = (int)f1 ; intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
        DiagBreak() ;
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ; free(intensities) ;
  }
  GCAhistoScaleImageIntensities(gca, mri_in) ;
  if (stricmp(xform_fname, "none"))
  {
    lta = LTAread(xform_fname) ;
    if (!lta)
      ErrorExit(ERROR_NOFILE, "%s: could not open transform", xform_fname) ;
  }
  else
    lta = LTAalloc(1, NULL) ;

  if (heq_fname)
  {
    MRI *mri_eq ;

    mri_eq = MRIread(heq_fname) ;
    if (!mri_eq)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read histogram equalization volume %s", 
                Progname, heq_fname) ;
    MRIhistoEqualize(mri_in, mri_eq, mri_in, 30, 170) ;
    MRIfree(&mri_eq) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      fprintf(stderr, "writing equalized volume to %s...\n", "heq.mgh") ;
      MRIwrite(mri_in, "heq.mgh") ;
    }
  }

  if (read_flag)
  {
    mri_labeled = MRIread(read_fname) ;
    if (!mri_labeled)
      ErrorExit(ERROR_NOFILE, "%s: could not read parcellation from %s",
                Progname, out_fname) ;
  }
  else
  {
    printf("labeling volume...\n") ;
    mri_labeled = GCAlabel(mri_in, gca, NULL, lta) ;
    if (wm_fname)
    {
      MRI *mri_wm ;
      
      mri_wm = MRIread(wm_fname) ;
      if (!mri_wm)
        ErrorExit(ERROR_NOFILE, "%s: could not read wm segmentation from %s",
                  Progname, wm_fname) ;
      mri_fixed = insert_wm_segmentation(mri_labeled,mri_wm,parcellation_type,
                                         fixed_flag, gca, lta);
      if (DIAG_VERBOSE_ON)
      {
        fprintf(stderr, "writing patched labeling to %s...\n", out_fname) ;
        MRIwrite(mri_labeled, out_fname) ;
      }
      MRIfree(&mri_wm) ;
    }
    else
      mri_fixed = MRIclone(mri_labeled, NULL) ;
    
    if (gca_write_iterations != 0)
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s_pre.mgh", gca_write_fname) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri_labeled, fname) ;
    }
    preprocess(mri_in, mri_labeled, gca, lta, mri_fixed) ;
    if (renormalize)
      GCArenormalizeLabels(mri_in, mri_labeled, gca, lta) ;
    if (fixed_flag == 0)
      MRIfree(&mri_fixed) ;
    if (!no_gibbs)
    {
      if (anneal)
        GCAanneal(mri_in, gca, mri_labeled, lta, max_iter) ;
      else
        GCAreclassifyUsingGibbsPriors(mri_in, gca, mri_labeled, lta, max_iter,
                                      mri_fixed, 0);
    }
  }
  GCAmaxLikelihoodBorders(gca, mri_in, mri_labeled, mri_labeled,lta,mle_niter,
                          5.0);
  if (expand_flag)
  {
    GCAexpandVentricle(gca, mri_in, mri_labeled, mri_labeled, lta,   
                       Left_Lateral_Ventricle) ;
    GCAexpandVentricle(gca, mri_in, mri_labeled, mri_labeled, lta,   
                       Right_Lateral_Ventricle) ;
#if 0
    GCAexpandVentricle(gca, mri_in, mri_labeled, mri_labeled, lta,   
                       Left_Inf_Lat_Vent) ;
    GCAexpandVentricle(gca, mri_in, mri_labeled, mri_labeled, lta,   
                       Right_Inf_Lat_Vent) ;
#endif
    GCAexpandCortex(gca, mri_in, mri_labeled, mri_labeled, lta) ;
  }

  if (mask_volume_fname)
  {
    MRI *mri_mask ;
    printf("reading volume %s for masking...\n", mask_volume_fname) ;
    mri_mask = MRIread(mask_volume_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not read mask volume from %s",
                Progname, mask_volume_fname) ;
    MRIthresholdMask(mri_labeled, mri_mask, mri_labeled, 1, 0) ;
    MRIfree(&mri_mask) ;
  }
    
  if (gca_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_post.mgh", gca_write_fname) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
  }

  GCAconstrainLabelTopology(gca, mri_in, mri_labeled, mri_labeled, lta) ;
  GCAfree(&gca) ; MRIfree(&mri_in) ;
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
  MRIwrite(mri_labeled, out_fname) ;

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("auto-labeling took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  exit(0) ;
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
  else if (!stricmp(option, "WM"))
  {
    wm_fname = argv[2] ;
    nargs = 1 ;
    printf("inserting white matter segmentation from %s...\n", wm_fname) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Ggca_x,Ggca_y,Ggca_z) ;
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
    printf("inserting fixed white matter segmentation from %s...\n",wm_fname);
  }
  else if (!stricmp(option, "MRI"))
  {
    mri_fname = argv[2] ;
    nargs = 1 ;
    printf("building most likely MR volume and writing to %s...\n", mri_fname);
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
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!stricmp(option, "renormalize"))
  {
    renormalize = atoi(argv[2]) ;
    nargs = 1 ;
    printf("renormalizing class means after initial labeling\n") ;
  }
  else switch (toupper(*option))
  {
  case 'R':
    read_flag = 1 ;
    read_fname = argv[2] ;
    nargs = 1 ;
    break ;
  case 'A':
    anneal = 1 ;
    fprintf(stderr, "using simulated annealing to find optimum\n") ;
    break ;
  case 'W':
    gca_write_iterations = atoi(argv[2]) ;
    gca_write_fname = argv[3] ;
    nargs = 2 ;
    printf("writing out snapshots of gibbs process every %d iterations to %s\n",
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
static void
usage_exit(int code)
{
  printf("usage: %s [options] <input volume> <xform> <gca file>"
         " <output volume>\n", Progname) ;
  exit(code) ;
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

static MRI *
insert_wm_segmentation(MRI *mri_labeled, MRI *mri_wm, int parcellation_type, 
                       int fixed_flag, GCA *gca,LTA *lta)
{
  int      x, y, z, width, depth, height, change_label[1000], n, label,
           nchanged, rh, lh, xn, yn, zn ;
  MRI      *mri_fixed ;
  GCA_NODE *gcan ;
  double   wm_prior ;

  mri_fixed = MRIclone(mri_wm, NULL) ;

  memset(change_label, 0, sizeof(change_label)) ;
  for (n = 0 ; n < NEDITABLE_LABELS ; n++)
    change_label[cma_editable_labels[n]] = 1 ;
  width = mri_wm->width ; height = mri_wm->height ; depth = mri_wm->depth ;

  for (nchanged = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 103 && y == 121 && z == 122)
          DiagBreak() ;
        if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL || MRIvox(mri_wm,x,y,z)>=200)
          continue ;
        label = MRIvox(mri_labeled, x, y, z) ;
        if (!change_label[label])
          continue ;
        
        GCAsourceVoxelToNode(gca, mri_labeled, lta, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        wm_prior = 0.0 ;
        for (n = 0 ; n < gcan->nlabels ; n++)
          if (gcan->labels[n] == Right_Cerebral_White_Matter ||
              gcan->labels[n] == Left_Cerebral_White_Matter)
            wm_prior = gcan->gcs[n].prior ;
#define PRIOR_THRESH 0.1
#define FIXED_PRIOR  0.5
        if (wm_prior < PRIOR_THRESH)
          continue ;
        lh = MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                Left_Cerebral_White_Matter) ;
        lh += MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                Left_Cerebral_Cortex) ;
        rh = MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                Right_Cerebral_White_Matter) ;
        rh += MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                Right_Cerebral_Cortex) ;
        if (rh > lh)
          label = Right_Cerebral_White_Matter ;
        else
          label = Left_Cerebral_White_Matter ;
        if (label != MRIvox(mri_labeled, x, y, z))
          nchanged++ ;
        MRIvox(mri_labeled, x, y, z) = label ;
#if 0
        if (MRIvox(mri_wm, x, y, z) < 140)
#endif
          MRIvox(mri_fixed, x, y, z) = fixed_flag ;
      }
    }
  }

  printf("%d labels changed to wm\n", nchanged) ;
  return(mri_fixed) ;
}

static int
MRIcountNbhdLabels(MRI *mri, int x, int y, int z, int label)
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
          total++ ;
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

static int
preprocess(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta,
           MRI *mri_fixed)
{
  int i ;

  GCArelabel_cortical_gray_and_white(gca, mri_in, mri_labeled,mri_labeled,lta);
  if (gca_write_iterations != 0)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s_cortex.mgh", gca_write_fname) ;
    printf("writing snapshot to %s...\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
  }
  if (hippocampus_flag)
  {
    edit_hippocampus(mri_in, mri_labeled, gca, lta, mri_fixed) ;
    edit_amygdala(mri_in, mri_labeled, gca, lta, mri_fixed) ;
  }
  if (expand_flag) 
  {
    for (i = 0 ; i < NEXPANDABLE_LABELS ; i++)
      GCAexpandLabelIntoWM(gca, mri_in, mri_labeled, mri_labeled, lta,
                           mri_fixed, cma_expandable_labels[i]) ;
#if 0
    GCAexpandVentricleIntoWM(gca, mri_in, mri_labeled, mri_labeled, lta,
                             mri_fixed,Left_Lateral_Ventricle) ;
    GCAexpandVentricleIntoWM(gca, mri_in, mri_labeled, mri_labeled, lta,
                             mri_fixed,Left_Inf_Lat_Vent) ;
    GCAexpandVentricleIntoWM(gca, mri_in, mri_labeled, mri_labeled, lta,
                             mri_fixed,Right_Lateral_Ventricle) ;
    GCAexpandVentricleIntoWM(gca, mri_in, mri_labeled, mri_labeled, lta,
                             mri_fixed, Right_Inf_Lat_Vent) ;
#endif
  }
#if 0
  GCAmaxLikelihoodBorders(gca, mri_in, mri_labeled, mri_labeled, lta,2,2.0) ;
#endif
  return(NO_ERROR) ;
}

static int
distance_to_label(MRI *mri_labeled, int label, int x, int y, int z, int dx, 
                  int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ; yi = y + d * dy ; zi = z + d * dz ;
    xi = mri_labeled->xi[xi] ; 
    yi = mri_labeled->yi[yi] ; 
    zi = mri_labeled->zi[zi];
    if (MRIvox(mri_labeled, xi, yi, zi) == label)
      break ;
  }

  return(d) ;
}
static int
edit_hippocampus(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta, 
                 MRI *mri_fixed)
{
  int   width, height, depth, x, y, z, label, nchanged, dleft, olabel,
        dright, dpos, dant, dup, ddown, dup_hippo, ddown_gray, i, left,
        ddown_ventricle, yi, found_wm, found_hippo, was_wm, was_hippo,
        was_unknown ;
  MRI   *mri_tmp ;
  float phippo, pwm ;

  mri_tmp = MRIcopy(mri_labeled, NULL) ;

  width = mri_in->width ; height = mri_in->height ; depth = mri_in->depth ;

  for (nchanged = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
          DiagBreak() ;
        label = MRIvox(mri_labeled, x, y, z) ;

        left = 0 ;
        switch (label)
        {
        case Left_Hippocampus:
          pwm = GCAlabelProbability(mri_in, gca, lta, x, y, z, 
                                    Left_Cerebral_White_Matter) ;
          phippo = GCAlabelProbability(mri_in, gca, lta, x, y, z, 
                                       Left_Hippocampus) ;
          if (phippo > 2*pwm)
            continue ;  /* don't let it change */

          dleft = distance_to_label(mri_labeled, Left_Cerebral_White_Matter,
                                    x,y,z,-1,0,0,3);
          dright = distance_to_label(mri_labeled, Left_Cerebral_White_Matter,
                                     x,y,z,1,0,0,3);
          ddown = distance_to_label(mri_labeled,
                                  Left_Cerebral_White_Matter,x,y,z,0,1,0,2);
          dup = distance_to_label(mri_labeled,
                                   Left_Cerebral_White_Matter,x,y,z,0,-1,0,2);
          dpos = distance_to_label(mri_labeled,
                                   Left_Cerebral_White_Matter,x,y,z,0,0,-1,3);
          dant = distance_to_label(mri_labeled,
                                   Left_Cerebral_White_Matter,x,y,z,0,0,1,3);

          if ((dpos <= 2 && dant <= 2) ||
              (dleft <= 2 && dright <= 2) ||
              (dup <= 1 && (dleft == 1 || dright == 1)))
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = Left_Cerebral_White_Matter ;
          }
          break ;
        case Right_Hippocampus:
          pwm = GCAlabelProbability(mri_in, gca, lta, x, y, z, 
                                    Right_Cerebral_White_Matter) ;
          phippo = GCAlabelProbability(mri_in, gca, lta, x, y, z, 
                                       Right_Hippocampus) ;
          if (phippo > 2*pwm)
            continue ;  /* don't let it change */

          dleft = distance_to_label(mri_labeled, Right_Cerebral_White_Matter,
                                    x,y,z,-1,0,0,3);
          dright = distance_to_label(mri_labeled, Right_Cerebral_White_Matter,
                                     x,y,z,1,0,0,3);
          ddown = distance_to_label(mri_labeled,
                                  Right_Cerebral_White_Matter,x,y,z,0,1,0,2);
          dup = distance_to_label(mri_labeled,
                                   Right_Cerebral_White_Matter,x,y,z,0,-1,0,2);
          dpos = distance_to_label(mri_labeled,
                                   Right_Cerebral_White_Matter,x,y,z,0,0,-1,3);
          dant = distance_to_label(mri_labeled,
                                   Right_Cerebral_White_Matter,x,y,z,0,0,1,3);
          if ((dpos <= 2 && dant <= 2) ||
              (dleft <= 2 && dright <= 2) ||
              (dup <= 1 && (dleft == 1 || dright == 1)))
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = Right_Cerebral_White_Matter ;
          }
          break ;
        default:
          break ;
        }
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
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Hippocampus:
            ddown = distance_to_label(mri_tmp, Left_Cerebral_White_Matter,
                                      x,y,z,0,1,0,2);
            dup = distance_to_label(mri_tmp,
                                    Left_Cerebral_White_Matter,x,y,z,0,-1,0,3);
            
            ddown_gray = distance_to_label(mri_tmp,Left_Cerebral_Cortex,
                                           x,y,z,0,1,0,3);
            dup_hippo = distance_to_label(mri_tmp,
                                          Left_Hippocampus,x,y,z,0,-1,0,3);
            
            ddown_ventricle = 
              MIN(distance_to_label(mri_tmp, Left_Lateral_Ventricle,
                                    x,y,z,0,1,0,3),
                  distance_to_label(mri_tmp, Left_Inf_Lat_Vent,
                                    x,y,z,0,1,0,3)) ;

            if (ddown_ventricle < ddown_gray)
              continue ;  /* gray should never be superior to ventricle */

            /* closer to superior white matter than hp and gray below */
            if (((dup < dup_hippo) && ddown > ddown_gray))
            {
              nchanged++ ;
              MRIvox(mri_tmp, x, y, z) = Left_Cerebral_Cortex ;
            }
            break ;
          case Right_Hippocampus:
            ddown = distance_to_label(mri_tmp,Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_tmp,Right_Cerebral_White_Matter,
                                    x,y,z,0,-1,0,3);
            
            ddown_gray = distance_to_label(mri_tmp,
                                           Right_Cerebral_Cortex,
                                           x,y,z,0,1,0,3);
            dup_hippo = distance_to_label(mri_tmp,
                                          Right_Hippocampus,x,y,z,0,-1,0,3);
            
            ddown_ventricle = 
              MIN(distance_to_label(mri_tmp, Right_Lateral_Ventricle,
                                    x,y,z,0,1,0,3),
                  distance_to_label(mri_tmp, Right_Inf_Lat_Vent,
                                    x,y,z,0,1,0,3)) ;

            if (ddown_ventricle < ddown_gray)
              continue ;  /* gray should never be superior to ventricle */

            /* closer to superior white matter than hp and gray below */
            if (((dup < dup_hippo) && ddown > ddown_gray))
            {
              nchanged++ ;
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
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            left = 1 ;
          case Right_Cerebral_Cortex:
            /*
              if the current label is gray, and there is white matter below
              and hippocampus above, change to hippocampus.
            */
            ddown = distance_to_label(mri_labeled,
                                      left ? Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 3)
            {
              nchanged++ ;
              MRIvox(mri_tmp, x, y, z) = left ? Left_Hippocampus : 
                Right_Hippocampus ;
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
          DiagBreak() ;
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
            if (olabel == Left_Hippocampus || olabel == Right_Hippocampus)
            {
              was_unknown = 0 ;
              continue ;
            }
            if (!IS_UNKNOWN(olabel))  /* don't change it */
              break ;
            if (++was_unknown >= MIN_UNKNOWN)
              break ;
          }
          if (was_unknown >= MIN_UNKNOWN)
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = left ? Left_Cerebral_Cortex : 
              Right_Cerebral_Cortex ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)  
            {
              printf("(%d, %d, %d) %s changed to %s in edit_hippocampus\n",
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
          DiagBreak() ;
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
              if (olabel == Left_Hippocampus || olabel == Right_Hippocampus)
              {
                was_hippo++ ;
                if (was_hippo >= 2)
                  found_hippo = 1 ;
              }
              else
                was_hippo = 0 ;
            }
            else if (olabel == Right_Cerebral_White_Matter ||
                     olabel == Left_Cerebral_White_Matter)
            {
              was_wm++ ;
              if (was_wm >= 3)
                found_wm = 1 ;
            }
            else
              was_wm = 0 ;
          }
          if (found_wm && found_hippo)
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = left ? Left_Cerebral_Cortex : 
              Right_Cerebral_Cortex ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)  
            {
              printf("(%d, %d, %d) %s changed to %s in edit_hippocampus\n",
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
  printf("%d hippocampal voxels changed.\n", nchanged) ;
  return(NO_ERROR) ;
}
static int
edit_amygdala(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta, 
                 MRI *mri_fixed)
{
  int   width, height, depth, x, y, z, label, nchanged, olabel,
        i, left, yi, found_wm, found_amygdala, was_wm, was_amygdala ;
  MRI   *mri_tmp ;

  mri_tmp = MRIcopy(mri_labeled, NULL) ;

  width = mri_in->width ; height = mri_in->height ; depth = mri_in->depth ;

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
          DiagBreak() ;
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
              if (olabel == Left_Amygdala || olabel == Right_Amygdala)
              {
                was_amygdala++ ;
                if (was_amygdala >= 2)
                  found_amygdala = 1 ;
              }
              else
                was_amygdala = 0 ;
            }
            else if (olabel == Right_Cerebral_White_Matter ||
                     olabel == Left_Cerebral_White_Matter)
            {
              was_wm++ ;
              if (was_wm >= 3)
                found_wm = 1 ;
            }
            else
              was_wm = 0 ;
          }
          if (found_wm && found_amygdala)
          {
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = left ? Left_Cerebral_Cortex : 
              Right_Cerebral_Cortex ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z)  
            {
              printf("(%d, %d, %d) %s changed to %s in edit_amygdala\n",
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
