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

static int MRIcountNbhdLabels(MRI *mri, int x, int y, int z, int label) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

static int filter = 0 ;
static float thresh = 0.5 ;
static int read_flag = 0 ;

static char *wm_fname = NULL ;
static char *heq_fname = NULL ;
static int max_iter = 25 ;
static int no_gibbs = 0 ;
static int anneal = 0 ;
static char *mri_fname = NULL ;

#define CMA_PARCELLATION  0
static int parcellation_type = CMA_PARCELLATION ;
static MRI *insert_wm_segmentation(MRI *mri_labeled, MRI *mri_wm, 
                                  int parcellation_type) ;

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

  /*  fprintf(stderr, "mri_in read: xform %s\n", mri_in->transform_fname) ;*/
  printf("reading classifier array from %s...\n", gca_fname) ;
  gca = GCAread(gca_fname) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read classifier array from %s",
              Progname, gca_fname) ;

  if (mri_fname)
  {
    GCAbuildMostLikelyVolume(gca, mri_in) ;
    MRIwrite(mri_in, mri_fname) ;
    exit(0) ;
  }
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
    fprintf(stderr, "writing equalized volume to %s...\n", out_fname) ;
    MRIfree(&mri_eq) ;
    MRIwrite(mri_in, out_fname) ;
  }

  if (read_flag)
  {
    mri_labeled = MRIread(out_fname) ;
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
      mri_fixed = insert_wm_segmentation(mri_labeled,mri_wm,parcellation_type);
      if (DIAG_VERBOSE_ON)
      {
        fprintf(stderr, "writing patched labeling to %s...\n", out_fname) ;
        MRIwrite(mri_labeled, out_fname) ;
      }
      MRIfree(&mri_wm) ;
    }
    if (!no_gibbs)
    {
      if (anneal)
        GCAanneal(mri_in, gca, mri_labeled, lta, max_iter) ;
      else
        GCAreclassifyUsingGibbsPriors(mri_in, gca, mri_labeled, lta, max_iter,
                                      mri_fixed);
    }
  }
  GCAconstrainLabelTopology(gca, mri_in, mri_labeled, mri_labeled, lta) ;
  GCAfree(&gca) ; MRIfree(&mri_in) ;
  if (filter)
  {
    MRI *mri_tmp ;

    printf("filtering labeled volume...\n") ;
    mri_tmp = MRIthreshModeFilter(mri_labeled, NULL, filter, thresh) ;
    MRIfree(&mri_labeled) ;
    mri_labeled = mri_tmp ;
  }

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
    printf("reading white matter segmentation from %s...\n", wm_fname) ;
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
  else switch (toupper(*option))
  {
  case 'R':
    read_flag = 1 ;
    break ;
  case 'A':
    anneal = 1 ;
    fprintf(stderr, "using simulated annealing to find optimum\n") ;
    break ;
  case 'N':
    max_iter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting max iterations to %d...\n", max_iter) ;
    break ;
  case 'F':
    filter = atoi(argv[2]) ;
    thresh = atof(argv[3]) ;
    nargs = 2 ;
    printf("applying thresholded (%2.2f) mode filter %d times to output of "
           "labelling\n",thresh,filter);
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

#define LEFT_HIPPOCAMPUS               17
#define RIGHT_HIPPOCAMPUS              53
#define LEFT_CEREBRAL_CORTEX           3
#define LEFT_CEREBRAL_WHITE_MATTTER    2
#define RIGHT_CEREBRAL_CORTEX          42
#define RIGHT_CEREBRAL_WHITE_MATTTER   41
#define LEFT_AMYGDALA                  18
#define RIGHT_AMYGDALA                 54


static int cma_editable_labels[] = 
{
  LEFT_HIPPOCAMPUS,
  RIGHT_HIPPOCAMPUS,
  LEFT_CEREBRAL_CORTEX,
  LEFT_CEREBRAL_WHITE_MATTTER,
  RIGHT_CEREBRAL_CORTEX,
  RIGHT_CEREBRAL_WHITE_MATTTER,
  LEFT_AMYGDALA,
  RIGHT_AMYGDALA
} ;
#define NEDITABLE_LABELS sizeof(cma_editable_labels) / sizeof(cma_editable_labels[0])

static MRI *
insert_wm_segmentation(MRI *mri_labeled, MRI *mri_wm, 
                                  int parcellation_type)
{
  int      x, y, z, width, depth, height, change_label[1000], n, label,
           nchanged, rh, lh ;
  MRI      *mri_fixed ;

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
        if (x == 153 && y == 127 && z == 128)
          DiagBreak() ;
        if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL || MRIvox(mri_wm,x,y,z)>=200)
          continue ;
        label = MRIvox(mri_labeled, x, y, z) ;
        if (!change_label[label])
          continue ;
        
        lh = MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                LEFT_CEREBRAL_WHITE_MATTTER) ;
        lh += MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                LEFT_CEREBRAL_CORTEX) ;
        rh = MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                RIGHT_CEREBRAL_WHITE_MATTTER) ;
        rh += MRIcountNbhdLabels(mri_labeled, x, y, z, 
                                RIGHT_CEREBRAL_CORTEX) ;
        if (rh > lh)
          label = RIGHT_CEREBRAL_WHITE_MATTTER ;
        else
          label = LEFT_CEREBRAL_WHITE_MATTTER ;
        if (label != MRIvox(mri_labeled, x, y, z))
          nchanged++ ;
        MRIvox(mri_labeled, x, y, z) = label ;
        if (MRIvox(mri_wm, x, y, z) < 140)
          MRIvox(mri_fixed, x, y, z) = 1 ;
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

