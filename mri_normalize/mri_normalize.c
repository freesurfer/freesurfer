#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "timer.h"
#include "proto.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(int code) ;

static int conform = 0 ;
static int gentle_flag = 0 ;

char *Progname ;

static MRI_NORM_INFO  mni ;
static int verbose = 1 ;
static int num_3d_iter = 2 ;

static float intensity_above = 20 ;
static float intensity_below = 10 ;

static char *control_point_fname ;

static char *control_volume_fname = NULL ;
static char *bias_volume_fname = NULL ;

static int no1d = 0 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, n ;
  MRI    *mri_src, *mri_dst = NULL ;
  char   *in_fname, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_normalize.c,v 1.22 2003/06/16 18:12:21 fischl Exp $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  mni.max_gradient = MAX_GRADIENT ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(0) ;

  if (argc < 1)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[1] ;

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no output name specified", Progname) ;
  out_fname = argv[2] ;

  if (verbose)
    fprintf(stderr, "reading from %s...\n", in_fname) ;
  mri_src = MRIread(in_fname) ;
  if (!mri_src)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s", 
              Progname, in_fname) ;
#if 0
  if ((mri_src->type != MRI_UCHAR) ||
      (!(mri_src->xsize == 1 && mri_src->ysize == 1 && mri_src->zsize == 1)))
#else
    if (conform || mri_src->type != MRI_UCHAR)
#endif
  {
    MRI  *mri_tmp ;

    fprintf(stderr, 
            "downsampling to 8 bits and scaling to isotropic voxels...\n") ;
    mri_tmp = MRIconform(mri_src) ;
    mri_src = mri_tmp ;
  }

  if (verbose)
    fprintf(stderr, "normalizing image...\n") ;
  TimerStart(&start) ;
  if (!no1d)
  {
    MRInormInit(mri_src, &mni, 0, 0, 0, 0, 0.0f) ;
    mri_dst = MRInormalize(mri_src, NULL, &mni) ;
    if (!mri_dst)
      ErrorExit(ERROR_BADPARM, "%s: normalization failed", Progname) ;
  }
  else
  {
    mri_dst = MRIcopy(mri_src, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_BADPARM, "%s: could not allocate volume", Progname) ;
  }

  if (control_point_fname)
    MRI3dUseFileControlPoints(mri_dst, control_point_fname) ;
  if (control_volume_fname)
    MRI3dWriteControlPoints(control_volume_fname) ;
  if (bias_volume_fname)
    MRI3dWriteBias(bias_volume_fname) ;
  for (n = 0 ; n < num_3d_iter ; n++)
  {
    fprintf(stderr, "3d normalization pass %d of %d\n", n+1, num_3d_iter) ;
    if (gentle_flag)
      MRI3dGentleNormalize(mri_dst, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE,
                           mri_dst,
                           intensity_above/2, intensity_below/2,
                           control_point_fname != NULL && !n && no1d);
    else
      MRI3dNormalize(mri_dst, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE,mri_dst,
                     intensity_above, intensity_below,
                     control_point_fname != NULL && !n && no1d);
  }

  if (verbose)
    fprintf(stderr, "writing output to %s\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "3D bias adjustment took %d minutes and %d seconds.\n", 
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
  if (!stricmp(option, "no1d"))
  {
    no1d = 1 ;
    fprintf(stderr, "disabling 1d normalization...\n") ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = 1 ;
    fprintf(stderr, "interpolating and embedding volume to be 256^3...\n") ;
  }
  else if (!stricmp(option, "gentle"))
  {
    gentle_flag = 1 ;
    fprintf(stderr, "performing kinder gentler normalization...\n") ;
  }
  else switch (toupper(*option))
  {
  case 'W':
    control_volume_fname = argv[2] ;
    bias_volume_fname = argv[3] ;
    nargs = 2 ;
    fprintf(stderr, "writing ctrl pts to   %s\n", control_volume_fname) ;
    fprintf(stderr, "writing bias field to %s\n", bias_volume_fname) ;
    break ;
  case 'F':
    control_point_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using control points from file %s...\n", 
           control_point_fname) ;
    break ;
  case 'A':
    intensity_above = atof(argv[2]) ;
    fprintf(stderr, "using control point with intensity %2.1f above target.\n",
            intensity_above) ;
    nargs = 1 ;
    break ;
  case 'B':
    intensity_below = atof(argv[2]) ;
    fprintf(stderr, "using control point with intensity %2.1f below target.\n",
            intensity_below) ;
    nargs = 1 ;
    break ;
  case 'G':
    mni.max_gradient = atof(argv[2]) ;
    fprintf(stderr, "using max gradient = %2.3f\n", mni.max_gradient) ;
    nargs = 1 ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
    num_3d_iter = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "performing 3d normalization %d times\n", num_3d_iter) ;
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
static void
usage_exit(int code)
{
    printf("usage: %s <input volume> <output volume>\n\n", Progname) ;
    printf("\t-slope <float s>  set the curvature slope (both n and p)\n");
    printf("\t-pslope <float p> set the curvature pslope (default=%2.1f)\n", pslope);
    printf("\t-nslope <float n> set the curvature nslope (default=%2.1f)\n", nslope);
    printf("\t-debug_voxel <int x y z> set voxel for debugging\n");
    printf("\t-auto             automatically detect class statistics (default)\n");
    printf("\t-noauto           don't automatically detect class statistics\n");
    printf("\t-log              log to ./segment.dat\n");
    printf("\t-keep             keep wm edits. maintains all values of 0 and 255\n");
    printf("\t-ghi, -gray_hi <int h> set the gray matter high limit (default=%d)\n", gray_hi);
    printf("\t-wlo, -wm_low  <int l> set the white matter low limit (default=%d)\n", wm_low);
    printf("\t-whi, -wm_hi <int h>   set the white matter high limit (default=%d)\n", wm_hi);
    printf("\t-nseg <int n>      thicken the n largest thin strands (default=%d)\n", nsegments);
    printf("\t-thicken           toggle thickening step (default=ON)\n");
    printf("\t-fillbg            toggle filling of the basal ganglia (default=OFF)\n");
    printf("\t-fillv             toggle filling of the ventricles (default=OFF)\n");
    printf("\t-b <float s>       set blur sigma (default=%2.2f)\n", blur_sigma);
    printf("\t-n <int i>         set # iterations of border classification (default=%d)\n", niter);
    printf("\t-t <int t>         set limit to thin strands in mm (default=%d)\n", thickness);
    printf("\t-v                 verbose\n");
    printf("\t-p <float p>       set % threshold (default=%2.2f)\n", pct);
    printf("\t-x <filename>      extract options from filename\n");
    printf("\t-w <int w>         set wsize (default=%d)\n", wsize);
    printf("\t-u                 usage\n");
    exit(code);
}

