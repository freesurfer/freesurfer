#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "gca.h"
#include "transform.h"
#include "version.h"

static char vcid[] = "$Id: mri_extract_label.c,v 1.6 2004/06/30 14:59:27 vicka Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  extract_labeled_image(MRI *mri_in, MATRIX *m, int label, MRI *mri_out) ;

char *Progname ;
static char *out_like_fname = NULL ;

static char *xform_fname = NULL ;
static LTA  *lta = NULL;
static float sigma = 0 ;

int
main(int argc, char *argv[])
{
  char        **av, *in_vol, *out_vol, out_fname[STRLEN] ;
  int         ac, nargs, i, invert_flag = 0, ras_flag = 0, label ;
  MRI         *mri_in, *mri_out, *mri_kernel, *mri_smoothed ;
  MATRIX      *m ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_extract_label.c,v 1.6 2004/06/30 14:59:27 vicka Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  in_vol = argv[1] ;
  out_vol = argv[argc-1] ;

  printf("reading volume from %s...\n", in_vol) ;
  mri_in = MRIread(in_vol) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, 
              in_vol) ;
  if (out_like_fname)
  {
    MRI *mri_tmp = MRIread(out_like_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",out_like_fname) ;
    mri_out = MRIalloc(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type) ;
    /*    MRIcopyHeader(mri_tmp, mri_out) ;*/
    MRIfree(&mri_tmp) ;
  }
  else
    mri_out = MRIalloc(256, 256, 256, mri_in->type) ;

  if (lta)
  {
    m = MatrixCopy(lta->xforms[0].m_L, NULL) ;
    
    if (lta->type == LINEAR_RAS_TO_RAS)  /* convert it to a voxel transform */
    {
      ras_flag = 1 ;
    }
    else if (ras_flag)
      ErrorExit(ERROR_UNSUPPORTED, "%s: transforms must be all RAS or all voxel",Progname) ;
    
    if (invert_flag)
    {
      MATRIX *m_tmp ;
      printf("inverting transform...\n") ;
      m_tmp = MatrixInverse(m, NULL) ;
      if (!m_tmp)
        ErrorExit(ERROR_BADPARM, "%s: transform is singular!") ;
      MatrixFree(&m) ; m = m_tmp ;
      invert_flag = 0 ;
    }
    LTAfree(&lta) ;
    if (ras_flag)  /* convert it to a voxel transform */
    {
      MATRIX *m_tmp ;
      printf("converting RAS xform to voxel xform...\n") ;
      m_tmp = MRIrasXformToVoxelXform(mri_in, mri_out, m, NULL) ;
      MatrixFree(&m) ; m = m_tmp ;
    }
  }
  else
    m = NULL ;

  for (i = 2 ; i < argc-1 ; i++)
  {
    label = atoi(argv[i]) ;
    printf("extracting label %d (%s)\n", label, cma_label_to_name(label)) ;
    extract_labeled_image(mri_in, m, label, mri_out) ;
    if (!FZERO(sigma))
    {
      printf("smoothing extracted volume...\n") ;
      mri_kernel = MRIgaussian1d(sigma, 30) ;
      mri_smoothed = MRIconvolveGaussian(mri_out, NULL, mri_kernel) ;
      MRIfree(&mri_out) ; mri_out = mri_smoothed ;
    }
    /* removed for gcc3.3
     * vsprintf(out_fname, out_vol, (va_list) &label) ;
    */
    sprintf(out_fname, out_vol, label) ;
    printf("writing output to %s.\n", out_fname) ;
    MRIwrite(mri_out, out_fname) ;
  }
  

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "out_like") || !stricmp(option, "ol"))
  {
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  }
  else switch (toupper(*option))
  {
  case 'S':
    sigma = atof(argv[2]) ;
    printf("applying sigma=%2.1f smoothing kernel after extraction...\n",sigma) ;
    nargs = 1 ;
    break ;
  case 'T':
    xform_fname = argv[2] ;
    printf("reading and applying transform %s...\n", xform_fname) ;
    nargs = 1 ;
    lta = LTAread(xform_fname) ;
    if (!lta)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", 
                Progname, xform_fname) ;

    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input volume> <label 1> <label 2> ... <output name>\n",
          Progname) ;
  fprintf(stderr, "where options are:\n") ;
  fprintf(stderr, 
          "\t-s <sigma>\tapply a Gaussian smoothing kernel\n"
          "\t-t <xform file>\tapply the transform in <xform file> to extracted volume\n");
}

static void
print_help(void)
{
  fprintf(stderr, 
          "\nThis program will extract a set of labeled voxels from an image\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static int
extract_labeled_image(MRI *mri_src, MATRIX *mA, int label, MRI *mri_dst)
{
  MRI  *mri_binarized ;

  mri_binarized = MRIclone(mri_src, NULL) ;
  MRIcopyLabel(mri_src, mri_binarized, label) ;
  MRIbinarize(mri_binarized, mri_binarized, 1, 0, 255) ;
  if (mA)
    MRIlinearTransform(mri_binarized, mri_dst, mA) ;
  else
    MRIcopy(mri_binarized, mri_dst) ;
  MRIfree(&mri_binarized) ;
  return(NO_ERROR) ;
}

