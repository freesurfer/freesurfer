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

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int replaceLabels(MRI *mri_seg) ;

static int gca_flags = GCA_NO_FLAGS ;

char *Progname ;
static void usage_exit(int code) ;
static char *mask_fname = NULL ;
static char *insert_fname = NULL ;
static int  insert_label = 0 ;
static char *histo_fname = NULL ;

static GCA_PARMS parms ;
static char *seg_dir = "seg" ;
static char *orig_dir = "orig" ;
static char *xform_name = "talairach.xfm" ;
static int prune = 0 ;
static float smooth = -1 ;

static int ninputs = 1 ;  /* T1 intensity */
static int navgs = 0 ;

static char subjects_dir[STRLEN] ;
static char *heq_fname = NULL ;

int
main(int argc, char *argv[])
{
  char         **av, fname[STRLEN], *out_fname, *subject_name, *cp ;
  int          ac, nargs, i, n, noint = 0, options ;
  int          msec, minutes, seconds, nsubjects ;
  struct timeb start ;
  GCA          *gca, *gca_prune = NULL ;
  MRI          *mri_seg, *mri_T1, *mri_eq = NULL ;
  TRANSFORM    *transform ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  parms.use_gradient = 0 ;
  parms.spacing = 4.0f ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment", 
                Progname);
    strcpy(subjects_dir, cp) ;
    if (argc < 3)
      usage_exit(1) ;
  }


  if (heq_fname)
  {
    mri_eq = MRIread(heq_fname) ;
    if (!mri_eq)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not read histogram equalization volume %s", 
                Progname, heq_fname) ;
  }

  out_fname = argv[argc-1] ;
  nsubjects = argc-2 ;
  for (options = i = 0 ; i < nsubjects ; i++)
  {
    if (argv[i+1][0] == '-')
    {
      nsubjects-- ; options++ ;
    }
  }

  printf("training on %d subject and writing results to %s\n",
          nsubjects, out_fname) ;

  n = 0 ;
  do
  {
    gca = GCAalloc(ninputs, parms.spacing, DEFAULT_VOLUME_SIZE, 
                   DEFAULT_VOLUME_SIZE,DEFAULT_VOLUME_SIZE, gca_flags);

    for (nargs = i = 0 ; i < nsubjects+options ; i++)
    {
      subject_name = argv[i+1] ;
      if (stricmp(subject_name, "-NOINT") == 0)
      {
        printf("not using intensity information for subsequent subjects...\n");
        noint = 1 ; nargs++ ;
        continue ;
      }
      else if (stricmp(subject_name, "-INT") == 0)
      {
        printf("using intensity information for subsequent subjects...\n");
        noint = 0 ; nargs++ ;
        continue ;
      }
      printf("processing subject %s, %d of %d...\n", subject_name,i+1-nargs,
             nsubjects);
      sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, seg_dir) ;
      if (DIAG_VERBOSE_ON)
        printf("reading segmentation from %s...\n", fname) ;
      mri_seg = MRIread(fname) ;
      if (!mri_seg)
        ErrorExit(ERROR_NOFILE, "%s: could not read segmentation file %s",
                  Progname, fname) ;
      if (insert_fname)
      {
        MRI *mri_insert ;

        sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, insert_fname) ;
        mri_insert = MRIread(fname) ;
        if (mri_insert == NULL)
          ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s for insertion",
                    Progname, insert_fname) ;

        MRIbinarize(mri_insert, mri_insert, 1, 0, insert_label) ;
        MRIcopyLabel(mri_insert, mri_seg, insert_label) ;
        MRIfree(&mri_insert) ;
      }
        
      replaceLabels(mri_seg) ;
      MRIeraseBorderPlanes(mri_seg) ;
      
      sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, orig_dir) ;
      if (DIAG_VERBOSE_ON)
        printf("reading co-registered T1 from %s...\n", fname) ;
      mri_T1 = MRIread(fname) ;
      if (!mri_T1)
        ErrorExit(ERROR_NOFILE, "%s: could not read T1 data from file %s",
                  Progname, fname) ;
      
      if (mask_fname)
      {
        MRI *mri_mask ;
        
        sprintf(fname, "%s/%s/mri/%s", subjects_dir, subject_name, mask_fname);
        printf("reading volume %s for masking...\n", fname) ;
        mri_mask = MRIread(fname) ;
        if (!mri_mask)
          ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                    Progname, fname) ;
        
        MRImask(mri_T1, mri_mask, mri_T1, 0, 0) ;
        MRIfree(&mri_mask) ;
      }
      if (mri_eq && !noint)
      {
        printf("histogram equalizing input image...\n") ;
        MRIhistoEqualize(mri_T1, mri_eq, mri_T1, 30, 170) ;
      }
      
      if (xform_name)
      {
        sprintf(fname, "%s/%s/mri/transforms/%s", 
                subjects_dir, subject_name, xform_name) ;
        if (DIAG_VERBOSE_ON)
          printf("reading transform from %s...\n", fname) ;
        transform = TransformRead(fname) ;
        if (!transform)
          ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
                    Progname, fname) ;
        TransformInvert(transform, mri_T1) ;
      }
      else
        transform = TransformAlloc(LINEAR_VOXEL_TO_VOXEL, NULL) ;
      
      GCAtrain(gca, mri_T1, mri_seg, transform, gca_prune, noint) ;
      MRIfree(&mri_seg) ; MRIfree(&mri_T1) ; TransformFree(&transform) ;
    }
    GCAcompleteTraining(gca) ;
    if (gca_prune)
      GCAfree(&gca_prune) ;
    gca_prune = gca ;
  } while (n++ < prune) ;

  if (smooth > 0)
  {
    printf("regularizing conditional densities with smooth=%2.2f\n", smooth) ;
    GCAregularizeConditionalDensities(gca, smooth) ;
  }
  if (navgs)
  {
    printf("applying mean filter %d times to conditional densities\n", navgs) ;
    GCAmeanFilterConditionalDensities(gca, navgs) ;
  }

  printf("writing trained GCA to %s...\n", out_fname) ;
  GCAwrite(gca, out_fname) ;

  if (histo_fname)
  {
    FILE *fp ;
    int   histo_counts[10000], xn, yn, zn, max_count ;
    GCA_NODE  *gcan ;

    memset(histo_counts, 0, sizeof(histo_counts)) ;
    fp = fopen(histo_fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open histo file %s",
                Progname, histo_fname) ;

    max_count = 0 ;
    for (xn = 0 ; xn < gca->width;  xn++)
    {
      for (yn = 0 ; yn < gca->height ; yn++)
      {
        for (zn = 0 ; zn < gca->depth ; zn++)
        {
          gcan = &gca->nodes[xn][yn][zn] ;
          if (gcan->nlabels < 1)
            continue ;
          if (gcan->nlabels == 1 && IS_UNKNOWN(gcan->labels[0]))
            continue ;
          histo_counts[gcan->nlabels]++ ;
          if (gcan->nlabels > max_count)
            max_count = gcan->nlabels ;
        }
      }
    }
    max_count = 20 ;
    for (xn = 1 ; xn < max_count ;  xn++)
      fprintf(fp, "%d %d\n", xn, histo_counts[xn]) ;
    fclose(fp) ;
  }

  GCAfree(&gca) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("classifier array training took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
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
  if (!stricmp(option, "GRADIENT"))
  {
    parms.use_gradient = 1 ;
    ninputs += 3 ;  /* components of the gradient */
  }
  else if (!stricmp(option, "SPACING"))
  {
    parms.spacing = atof(argv[2]) ;
    nargs = 1 ;
    printf("spacing nodes every %2.1f mm\n", parms.spacing) ;
  }
  else if (!stricmp(option, "NOMRF"))
  {
    gca_flags |= GCA_NO_MRF ;
    printf("not computing MRF statistics...\n") ;
  }
  else if (!stricmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!stricmp(option, "DEBUG_NODE"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Ggca_x,Ggca_y,Ggca_z) ;
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging label %d\n", Ggca_label) ;
  }
  else if (!stricmp(option, "INSERT"))
  {
    insert_fname = argv[2] ;
    insert_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("inserting non-zero vals from %s as label %d...\n", insert_fname,insert_label);
  }
  else if (!stricmp(option, "PRUNE"))
  {
    prune = atoi(argv[2]) ;
    nargs = 1 ;
    printf("pruning classifier %d times after initial training\n", prune) ;
  }
  else if (!stricmp(option, "HEQ"))
  {
    heq_fname = argv[2] ;
    nargs = 1 ;
    printf("reading template for histogram equalization from %s...\n", 
           heq_fname) ;
  }
  else if (!stricmp(option, "T1"))
  {
    orig_dir = argv[2] ;
    nargs = 1 ;
    printf("reading T1 data from subject's mri/%s directory\n",
            orig_dir) ;
  }
  else if (!stricmp(option, "PARC_DIR") || !stricmp(option, "SEG_DIR"))
  {
    seg_dir = argv[2] ;
    nargs = 1 ;
    printf("reading segmentation from subject's mri/%s directory\n",
            seg_dir) ;
  }
  else if (!stricmp(option, "XFORM"))
  {
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("reading xform from %s\n", xform_name) ;
  }
  else if (!stricmp(option, "NOXFORM"))
  {
    xform_name = NULL ;
    printf("disabling application of xform...\n") ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "SMOOTH"))
  {
    smooth = atof(argv[2]) ;
    if (smooth <= 0 || smooth > 1)
      ErrorExit(ERROR_BADPARM, 
                "%s: smoothing parameter %2.1f must be in [0,1]\n",
                Progname, smooth) ;
    nargs = 1 ;
    printf("imposing %2.1f smoothing on conditional statistics\n", smooth) ;
  }
  else switch (toupper(*option))
  {
  case 'A':
    navgs = atoi(argv[2]) ;
    printf("applying %d mean filters to classifiers after training\n",navgs);
    nargs = 1 ;
    break ;
  case 'H':
    histo_fname = argv[2] ;
    nargs = 1 ;
    printf("writing histogram of classes/voxel to %s\n", histo_fname) ;
    break; 
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
  printf("usage: %s [options] <subject 1> <subject 2> ... <output file>\n",
         Progname) ;
  printf(
         "\t-spacing  - spacing of classifiers in canonical space\n");
  printf("\t-gradient - use intensity gradient as input to classifier.\n") ;
  exit(code) ;
}
static int input_labels[] = {  
  Left_Cerebral_Exterior,
  Right_Cerebral_Exterior,
  Left_Cerebellum_Exterior,
  Right_Cerebellum_Exterior
} ;
static int output_labels[] = {  
  Left_Cerebral_Cortex,
  Right_Cerebral_Cortex,
  Left_Cerebellum_Cortex,
  Right_Cerebellum_Cortex
} ;
  
static int
replaceLabels(MRI *mri_seg)
{
  int    i ;

  for (i = 0 ; i < sizeof(output_labels)/sizeof(output_labels[0]) ; i++)
    MRIreplaceValues(mri_seg, mri_seg, input_labels[i], output_labels[i]) ;
  return(NO_ERROR) ;
}

