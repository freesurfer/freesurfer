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
#include "tags.h"
#include "version.h"
#include "cma.h"

static MRI *MRIremoveWMOutliers(MRI *mri_src, 
				MRI *mri_src_ctrl, 
				MRI *mri_dst_ctrl) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(int code) ;
static MRI *compute_bias(MRI *mri_src, MRI *mri_dst, MRI *mri_bias) ;
static int conform = 1 ;
static int gentle_flag = 0 ;
static int nosnr = 1 ;

static float bias_sigma = 8.0 ;

static char *mask_fname ;
char *Progname ;

static int prune = 0 ;  /* off by default */
static MRI_NORM_INFO  mni = {} ;
static int verbose = 1 ;
static int num_3d_iter = 2 ;

static float intensity_above = 25 ;
static float intensity_below = 10 ;

static char *control_point_fname ;

static char *aseg_fname = NULL ;
//static int aseg_wm_labels[] =
// { Left_Cerebral_White_Matter, Right_Cerebral_White_Matter, Brain_Stem} ;
static int aseg_wm_labels[] =
  { Left_Cerebral_White_Matter, Right_Cerebral_White_Matter} ;
#define NWM_LABELS (sizeof(aseg_wm_labels) / sizeof(aseg_wm_labels[0]))

static char *control_volume_fname = NULL ;
static char *bias_volume_fname = NULL ;
static int read_flag = 0 ;

static int no1d = 0 ;
static int file_only = 0 ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, n ;
  MRI    *mri_src, *mri_dst = NULL, *mri_bias, *mri_orig, *mri_aseg = NULL ;
  char   *in_fname, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
    (argc, argv,
     "$Id: mri_normalize.c,v 1.48 2006/11/18 02:11:55 fischl Exp $",
     "$Name:  $",
     cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mri_normalize.c,v 1.48 2006/11/18 02:11:55 fischl Exp $",
     "$Name:  $");
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
  if (!mriConformed(mri_src) && conform > 0)
  {
    printf("unconformed source detected - conforming...\n") ;
    mri_src = MRIconform(mri_src) ;
  }

  MRIaddCommandLine(mri_src, cmdline) ;
  if (mask_fname)
  {
    MRI *mri_mask ;

    mri_mask = MRIread(mask_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, mask_fname) ;
    MRImask(mri_src, mri_mask, mri_src, 0, 0) ;
    MRIfree(&mri_mask) ;
  }
  if (read_flag)
  {
    MRI *mri_ctrl ;
    double scale ;

    mri_bias = MRIread(bias_volume_fname) ;
    if (!mri_bias)
      ErrorExit
        (ERROR_BADPARM,
         "%s: could not read bias volume %s", Progname, bias_volume_fname) ;
    mri_ctrl = MRIread(control_volume_fname) ;
    if (!mri_ctrl)
      ErrorExit
        (ERROR_BADPARM,
         "%s: could not read control volume %s",
         Progname, control_volume_fname) ;
    MRIbinarize(mri_ctrl, mri_ctrl, 1, 0, 128) ;
    mri_dst = MRImultiply(mri_bias, mri_src, NULL) ;
    scale = MRImeanInLabel(mri_dst, mri_ctrl, 128) ;
    printf("mean in wm is %2.0f, scaling by %2.2f\n", scale, 110/scale) ;
    scale = 110/scale ;
    MRIscalarMul(mri_dst, mri_dst, scale) ;
    MRIwrite(mri_dst, out_fname) ;
    exit(0) ;
  }

#if 0
#if 0
  if ((mri_src->type != MRI_UCHAR) ||
      (!(mri_src->xsize == 1 && mri_src->ysize == 1 && mri_src->zsize == 1)))
#else
    if (conform || (mri_src->type != MRI_UCHAR && conform > 0))
#endif
    {
      MRI  *mri_tmp ;

      fprintf
        (stderr,
         "downsampling to 8 bits and scaling to isotropic voxels...\n") ;
      mri_tmp = MRIconform(mri_src) ;
      mri_src = mri_tmp ;
    }
#endif

  if (aseg_fname)
  {
    mri_aseg = MRIread(aseg_fname) ;
    if (mri_aseg == NULL)
      ErrorExit
        (ERROR_NOFILE,
         "%s: could not read aseg from file %s", Progname, aseg_fname) ;
  }
  else
    mri_aseg = NULL ;

  if (verbose)
    fprintf(stderr, "normalizing image...\n") ;
  TimerStart(&start) ;

  if (control_point_fname)
    MRI3dUseFileControlPoints(mri_src, control_point_fname) ;
  if (control_volume_fname)
    // this just setup writing control-point volume saving
    MRI3dWriteControlPoints(control_volume_fname) ;

  /* first do a gentle normalization to get
     things in the right intensity range */
  if (control_point_fname != NULL)  /* do one pass with only
                                       file control points first */
    mri_dst =
      MRI3dGentleNormalize(mri_src, 
                           NULL, 
                           DEFAULT_DESIRED_WHITE_MATTER_VALUE,
                           NULL,
                           intensity_above, 
                           intensity_below/2,1, 
                           bias_sigma);
  else
    mri_dst = MRIcopy(mri_src, NULL) ;
	
  if (mri_aseg)
  {
    MRI *mri_ctrl, *mri_bias ;
    int  i ;

    mri_ctrl = MRIclone(mri_aseg, NULL) ;
    for (i = 0 ; i < NWM_LABELS ; i++)
      MRIcopyLabel(mri_aseg, mri_ctrl, aseg_wm_labels[i]) ;
    fprintf(stderr,"removing outliers in the aseg WM...\n") ;
    MRIremoveWMOutliers(mri_dst, mri_ctrl, mri_ctrl) ;
    MRIbinarize(mri_ctrl, mri_ctrl, 1, CONTROL_NONE, CONTROL_MARKED) ;
    MRInormAddFileControlPoints(mri_ctrl, CONTROL_MARKED) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_ctrl, "norm_ctrl.mgz") ;
    mri_bias = MRIbuildBiasImage(mri_dst, mri_ctrl, NULL, bias_sigma) ;
    MRIfree(&mri_ctrl) ; MRIfree(&mri_aseg) ;
    mri_dst = MRIapplyBiasCorrectionSameGeometry
      (mri_dst, mri_bias, mri_dst, DEFAULT_DESIRED_WHITE_MATTER_VALUE) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_dst, "norm_1.mgz") ;
  }
  else
  {
    if (!no1d)
    {
      MRInormInit(mri_src, &mni, 0, 0, 0, 0, 0.0f) ;
      mri_dst = MRInormalize(mri_src, NULL, &mni) ;
      if (!mri_dst)
        ErrorExit(ERROR_BADPARM, "%s: normalization failed", Progname) ;
    }
    else
    {
      if ((file_only && nosnr) ||
          ((gentle_flag != 0) && (control_point_fname != NULL)))
	    {
	      if (mri_dst == NULL)
          mri_dst = MRIcopy(mri_src, NULL) ;
	    }
      else
	    {
	      printf("computing initial normalization using SNR...\n") ;
	      mri_dst = MRInormalizeHighSignalLowStd
          (mri_src, mri_dst, bias_sigma,
           DEFAULT_DESIRED_WHITE_MATTER_VALUE) ;
	    }
      if (!mri_dst)
        ErrorExit
          (ERROR_BADPARM, "%s: could not allocate volume", Progname) ;
    }
  }

  if (file_only == 0)
    MRI3dGentleNormalize(mri_dst, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE,
                         mri_dst,
                         intensity_above, intensity_below/2,
                         file_only, bias_sigma);
  mri_orig = MRIcopy(mri_dst, NULL) ;
  for (n = 0 ; n < num_3d_iter ; n++)
  {
    if (file_only)
      break ;
    fprintf(stderr, "3d normalization pass %d of %d\n", n+1, num_3d_iter) ;
    if (gentle_flag)
      MRI3dGentleNormalize(mri_dst, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE,
                           mri_dst,
                           intensity_above/2, intensity_below/2,
                           file_only, bias_sigma);
    else
      MRI3dNormalize(mri_orig, mri_dst, DEFAULT_DESIRED_WHITE_MATTER_VALUE,
                     mri_dst,
                     intensity_above, intensity_below,
                     file_only, prune, bias_sigma);
  }

  if (bias_volume_fname)
  {
    mri_bias = compute_bias(mri_src, mri_dst, NULL) ;
    printf("writing bias field to %s....\n", bias_volume_fname) ;
    MRIwrite(mri_bias, bias_volume_fname) ;
    MRIfree(&mri_bias) ;
  }

  if (verbose)
    fprintf(stderr, "writing output to %s\n", out_fname) ;
  MRIwrite(mri_dst, out_fname) ;
  msec = TimerStop(&start) ;

  MRIfree(&mri_src);
  MRIfree(&mri_dst);

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
  else if (!stricmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!stricmp(option, "monkey"))
  {
    no1d = 1 ;
    num_3d_iter = 1 ;
    printf("disabling 1D normalization and "
           "setting niter=1, make sure to use "
           "-f to specify control points\n") ;
  }
  else if (!stricmp(option, "nosnr"))
  {
    nosnr = 1 ;
    printf("disabling SNR normalization\n") ;
  }
  else if (!stricmp(option, "snr"))
  {
    nosnr = 0 ;
    printf("enabling SNR normalization\n") ;
  }
  else if (!stricmp(option, "sigma"))
  {
    bias_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("using Gaussian smoothing of bias field, sigma=%2.3f\n",
           bias_sigma) ;
  }
  else if (!stricmp(option, "conform"))
  {
    conform = atoi(argv[1]) ;
    nargs = 1 ;
    fprintf(stderr, "%sinterpolating and embedding volume to be 256^3...\n", conform ? "": "not ") ;
  }
  else if (!stricmp(option, "noconform"))
  {
    conform = 0 ;
    fprintf(stderr, "%sinterpolating and embedding volume to be 256^3...\n", conform ? "": "not ") ;
  }
  else if (!stricmp(option, "aseg") || !stricmp(option, "segmentation"))
  {
    aseg_fname = argv[2] ;
    nargs = 1  ;
    fprintf(stderr,
            "using segmentation for initial intensity normalization\n") ;
  }
  else if (!stricmp(option, "gentle"))
  {
    gentle_flag = 1 ;
    fprintf(stderr, "performing kinder gentler normalization...\n") ;
  }
  else if (!stricmp(option, "file_only") ||
           !stricmp(option, "fonly") ||
           !stricmp(option, "fileonly"))
  {
    file_only = 1 ;
    control_point_fname = argv[2] ;
    no1d = 1 ;
    nargs = 1 ;
    fprintf(stderr, "using control points from file %s...\n",
            control_point_fname) ;
    fprintf(stderr, "only using file control points...\n") ;
  }
  else switch (toupper(*option))
  {
  case 'D':
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    break ;
  case 'V':
    Gvx = atoi(argv[2]) ;
    Gvy = atoi(argv[3]) ;
    Gvz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging alternative voxel (%d, %d, %d)\n", Gvx, Gvy, Gvz) ;
    break ;
  case 'P':
    prune = atoi(argv[2]) ;
    nargs = 1 ;
    printf("turning control point pruning %s\n", prune > 0 ? "on" : "off") ;
    if (prune == 0)
      prune = -1 ;
    break ;
  case 'R':
    read_flag = 1 ;
    nargs = 2 ;
    control_volume_fname = argv[2] ;
    bias_volume_fname = argv[3] ;
    printf("reading bias field from %s and ctrl points from %s\n",
           bias_volume_fname, control_volume_fname) ;
    break ;
  case 'W':
    control_volume_fname = argv[2] ;
    bias_volume_fname = argv[3] ;
    nargs = 2 ;
    printf("writing ctrl pts to   %s\n", control_volume_fname) ;
    printf("writing bias field to %s\n", bias_volume_fname) ;
    break ;
  case 'F':
    control_point_fname = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using control points from file %s...\n",
            control_point_fname) ;
    break ;
  case 'A':
    intensity_above = atof(argv[2]) ;
    fprintf(stderr,
            "using control point with intensity %2.1f above target.\n",
            intensity_above) ;
    nargs = 1 ;
    break ;
  case 'B':
    intensity_below = atof(argv[2]) ;
    fprintf(stderr,
            "using control point with intensity %2.1f below target.\n",
            intensity_below) ;
    nargs = 1 ;
    break ;
  case 'G':
    mni.max_gradient = atof(argv[2]) ;
    fprintf(stderr, "using max gradient = %2.3f\n", mni.max_gradient) ;
    nargs = 1 ;
    break ;
#if 0
  case 'V':
    verbose = !verbose ;
    break ;
#endif
  case 'N':
    num_3d_iter = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "performing 3d normalization %d times\n", num_3d_iter) ;
    break ;
  case '?':
  case 'U':
  case 'H':
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
  printf("%s input output\n\n", Progname) ;
  printf("  -no1d              disable 1d normalization\n");
  printf("  -conform           interpolate and embed volume to be 256^3\n");
  printf("  -gentle            perform kinder gentler normalization\n");
  printf("  -f <path to file>  use control points "
         "file (usually control.dat)\n");
  printf("  -w <mri_vol c> <mri_vol b> write ctrl point(c) "
         "and bias field(b) volumes\n");
  printf("  -a <float a>       use control point with "
         "intensity a above target (default=%2.1f)\n", intensity_above);
  printf("  -b <float b>       use control point with "
         "intensity b below target (default=%2.1f)\n", intensity_below);
  printf("  -g <float g>       use max intensity/mm "
         "gradient g (default=%2.3f)\n", mni.max_gradient);
  printf("  -v                 verbose\n");
  printf("  -v Gvx Gvy Gvz : for debuggin\n");
  printf("  -n <int n>         use n 3d normalization "
         "iterations (default=%d)\n", num_3d_iter);
  printf("  -u or -h            print usage\n");
  printf("  -prune <boolean>     turn pruning of control points "
         "on/off (default=off). Useful if white is expanding into gm\n") ;
  printf("  -MASK maskfile \n");
  printf("  -monkey : turns of 1d, sets num_3d_iter=1\n");
  printf("  -nosnr : disable snr normalization\n");
  printf("  -sigma sigma : smooth bias field\n");
  printf("  -aseg aseg\n");
  printf("  -d Gx Gy Gz : for debuggin\n");
  printf("  -r controlpoints biasfield : for reading\n");
  printf("  \n");

  exit(code);
}

static MRI *
compute_bias(MRI *mri_src, MRI *mri_dst, MRI *mri_bias)
{
  int x, y, z ;
  float bias, src, dst ;

  if (!mri_bias)
    mri_bias = MRIalloc
      (mri_src->width, mri_src->height, mri_src->depth, MRI_FLOAT) ;

  MRIcopyHeader(mri_src, mri_bias) ;
  for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0; y < mri_src->height ; y++)
        {
          for (z = 0 ; z < mri_src->depth ; z++)
            {
              src = MRIgetVoxVal(mri_src, x, y, z, 0) ;
              dst = MRIgetVoxVal(mri_dst, x, y, z, 0) ;
              if (FZERO(src))
                bias = 1 ;
              else
                bias = dst/src ;
              MRIsetVoxVal(mri_bias, x, y, z, 0, bias) ;
            }
        }
    }

  return(mri_bias) ;
}

static MRI *
MRIremoveWMOutliers(MRI *mri_src, MRI *mri_src_ctrl, MRI *mri_dst_ctrl)
{
  MRI  *mri_inside, *mri_bin ;
  HISTOGRAM *histo, *hsmooth ;
  int   wm_peak, x, y, z, nremoved ;
  float thresh ;
  double val, lmean, max ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_src_ctrl, "sc.mgz") ;
  mri_dst_ctrl = MRIerode(mri_src_ctrl, mri_dst_ctrl) ;
  mri_inside = MRIerode(mri_dst_ctrl, NULL) ;
  MRIbinarize(mri_inside, mri_inside, 1, 0, 1) ;

  histo = MRIhistogramLabel(mri_src, mri_inside, 1, 256) ;
  hsmooth = HISTOcopy(histo, NULL) ;
  HISTOsmooth(histo, hsmooth, 2) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      HISTOplot(histo, "h.plt") ;
      HISTOplot(hsmooth, "hs.plt") ;
    }
  wm_peak = HISTOfindHighestPeakInRegion(hsmooth, 1, hsmooth->nbins-1) ;
  wm_peak = hsmooth->bins[wm_peak] ;
  thresh = wm_peak-10 ;
  printf("using wm threshold %2.1f for removing exterior voxels\n", thresh) ;

  // now remove stuff that's on the border and is pretty dark
  for (nremoved = x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
	{
	  for (z = 0 ; z < mri_src->depth ; z++)
	    {
	      if (x == Gx && y == Gy && z == Gz)
		DiagBreak() ;
	      /* if it's a control point, 
		 it's not in the interior of the wm, 
		 and it's T1 val is too low */
	      if (MRIgetVoxVal(mri_dst_ctrl, x, y, z, 0) == 0)
		continue ;  // not a  control point

	      /* if it's way far from the wm mode 
		 then remove it even if it's in the interior */
	      val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
	      if (val < thresh-5)
		{
		  MRIsetVoxVal(mri_dst_ctrl, x, y, z, 0, 0) ;
		  nremoved++ ;
		}

	      if (nint(MRIgetVoxVal(mri_inside, x, y, z, 0)) > 0)  
		// don't process interior voxels further
		continue ;   // in the interior
	      if (val < thresh)
		{
		  MRIsetVoxVal(mri_dst_ctrl, x, y, z, 0, 0) ;
		  nremoved++ ;
		}
	      else
		{
		  lmean = 
		    MRImeanInLabelInRegion(mri_src, mri_inside, 1, x, y, z, 7);
		  if (val < lmean-10)
		    {
		      MRIsetVoxVal(mri_dst_ctrl, x, y, z, 0, 0) ;
		      nremoved++ ;
		    }
		}
	    }
	}
    }

#if 0
  for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
	{
	  for (z = 0 ; z < mri_src->depth ; z++)
	    {
	      if (x == Gx && y == Gy && z == Gz)
		DiagBreak() ;
	      /* if it's a control point, 
		 it's not in the interior of the wm, 
		 and it's T1 val is too low */
	      if (MRIgetVoxVal(mri_dst_ctrl, x, y, z, 0) == 0)
		continue ;  // not a  control point
	      if (MRIcountNonzeroInNbhd(mri_dst_ctrl,3, x, y, z)<=2)
		{
		  MRIsetVoxVal(mri_dst_ctrl, x, y, z, 0, 0) ;
		  nremoved++ ;
		}
	    }
	}
    }
#endif

  /* now take out voxels that have too big an intensity diff 
     with surrounding ones */
  mri_bin = MRIbinarize(mri_dst_ctrl, NULL, 1, 0, 1) ;
  for (x = 0 ; x < mri_src->width ; x++)
    {
      for (y = 0 ; y < mri_src->height ; y++)
	{
	  for (z = 0 ; z < mri_src->depth ; z++)
	    {
	      if (x == Gx && y == Gy && z == Gz)
		DiagBreak() ;
	      /* if it's a control point, 
		 it's not in the interior of the wm, 
		 and it's T1 val is too low */
	      if (MRIgetVoxVal(mri_dst_ctrl, x, y, z, 0) == 0)
		continue ;  // not a  control point
	      val = MRIgetVoxVal(mri_src, x, y, z, 0) ;
	      max = MRImaxInLabelInRegion(mri_src, mri_bin, 1, x, y, z, 3);
	      if (val+7 < max)
		{
		  MRIsetVoxVal(mri_dst_ctrl, x, y, z, 0, 0) ;
		  nremoved++ ;
		}
	    }
	}
    }

  fprintf(stderr, "%d control points removed\n", nremoved) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_dst_ctrl, "dc.mgz") ;
  HISTOfree(&histo) ; HISTOfree(&hsmooth) ; 
  MRIfree(&mri_inside) ; MRIfree(&mri_bin) ;
  return(mri_dst_ctrl) ;
}

