/**
 * @brief Normalize the white-matter, based on control points.
 *
 * The variation in intensity due to the B1 bias field is corrected.
 *
 * Reference:
 * "Cortical Surface-Based Analysis I: Segmentation and Surface
 * Reconstruction", Dale, A.M., Fischl, B., Sereno, M.I.
 * (1999) NeuroImage 9(2):179-194
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
#include "timer.h"
#include "proto.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "tags.h"
#include "version.h"
#include "cma.h"


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(int code) ;
static MRI *find_control_points_always_in_WM(MRI *mri_aseg, MRI *mri_ctrl) ;
static MRI *remove_points_not_in_range(MRI *mri_brain, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, float min_target, float max_target) ;

#if 0
static MRI *remove_gradient_outliers(MRI *mri_in, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, double max_grad) ;
static MRI *remove_absolute_outliers(MRI *mri_in, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, double min_intensity) ;
static float max_grad = 5 ;  // ctrl points that are this much below any nbr will be removed
static float min_intensity = 95 ;

#endif

static int intensity_pad = 2;

static char *mask_fname = NULL ;
static const char *aseg_name = "aseg.mgz" ;
static const char *brain_name = "brain.mgz" ;

static char *control_volume_fname = NULL ; ;
static char *bias_volume_fname = NULL ;
static char *control_point_fname = NULL ;
static float bias_sigma = 1.0 ;
static float cross_time_sigma = -1 ;

const char *Progname ;

#define MAX_TPS 1000

int
main(int argc, char *argv[])
{
  char   *tp_fname, *in_fname, *out_fname, *tp_names[MAX_TPS] ;
  int    nargs, t, ntps ;
  char line[STRLEN], *cp, fname[STRLEN], bdir[STRLEN], *bname, sdir[STRLEN] ;
  MRI   *mri_norm = NULL, *mri_aseg = NULL, *mri_ctrl, *mri_bias, *mri_dst, *mri_brain = NULL ;
  FILE  *fp ;
  LTA   *ltas[MAX_TPS] ;

  nargs = handleVersionOption(argc, argv, "mri_long_normalize");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  tp_fname = argv[1] ;
  in_fname = argv[2] ;
  out_fname = argv[3] ;

  ntps = FileNumberOfEntries(tp_fname) ;
  printf("reading time point file %s with %d timepoints\n", tp_fname, ntps) ;
  
  fp = fopen(tp_fname, "r") ;
  if (ntps <= 0 || fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read base tp file %s", Progname, tp_fname) ;

  FileNamePath(tp_fname, bdir) ;
  bname = strrchr(bdir, '/')+1 ;
  if (bname == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not parse base name from %s", Progname, bdir) ;
  strcpy(sdir, bdir) ;
  cp = strrchr(sdir, '/') ;
  if (cp == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not parse SUBJECTS_DIR from %s", Progname, sdir) ;
  *cp = 0 ;

  for (t = 0 ; t < ntps ; t++)
  {
    MRI *mri_tmp, *mri_tmp2 ;

     cp = fgetl(line, 199, fp) ;
     tp_names[t] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
     strcpy(tp_names[t], cp) ;

     int req = snprintf(fname, STRLEN, "%s/%s.long.%s/mri/%s", sdir, cp, bname, in_fname);
     if (req >= STRLEN) {
       std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
     }
     printf("reading input volume %s\n", fname) ;
     mri_tmp = MRIread(fname) ;
     if (mri_tmp == NULL)
       ErrorExit(Gerror, NULL) ;
     if (t == 0)
     {
       mri_norm = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type, ntps);
       MRIcopyHeader(mri_tmp, mri_norm) ;
     }
     MRIcopyFrame(mri_tmp, mri_norm, 0, t) ;
     MRIfree(&mri_tmp) ;

     req = snprintf(fname, STRLEN, "%s/%s/mri/%s", sdir, cp, brain_name);
     if (req >= STRLEN) {
       std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
     }
     printf("reading input volume %s\n", fname) ;
     mri_tmp = MRIread(fname) ;
     if (mri_tmp == NULL)
       ErrorExit(Gerror, NULL) ;
     if (t == 0)
     {
       mri_brain = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type, ntps);
       MRIcopyHeader(mri_tmp, mri_brain) ;
     }

     req = snprintf(fname, STRLEN, "%s/mri/transforms/%s_to_%s.lta", bdir, bname, tp_names[t]);
     if (req >= STRLEN) {
       std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
     }

     printf("reading input transform %s\n", fname) ;
     ltas[t] = LTAread(fname) ;
     if (ltas[t] == NULL)
       ErrorExit(Gerror, NULL) ;
     if (ltas[t]->type != LINEAR_VOX_TO_VOX)
       ErrorExit(ERROR_UNSUPPORTED, "%s: transforms must be linear vox to vox", Progname) ;
     mri_tmp2 = MRIinverseLinearTransform(mri_tmp, NULL, ltas[t]->xforms[0].m_L) ;
     MRIcopyFrame(mri_tmp2, mri_brain, 0, t) ;
     MRIfree(&mri_tmp) ; MRIfree(&mri_tmp2) ;

     req = snprintf(fname, STRLEN, "%s/%s.long.%s/mri/%s", sdir, cp, bname, aseg_name);
     if (req >= STRLEN) {
       std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
     }
     printf("reading input volume %s\n", fname) ;
     mri_tmp = MRIread(fname) ;
     if (mri_tmp == NULL)
       ErrorExit(Gerror, NULL) ;
     if (t == 0)
     {
       mri_aseg = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type, ntps);
       MRIcopyHeader(mri_tmp, mri_aseg) ;
     }
     MRIcopyFrame(mri_tmp, mri_aseg, 0, t) ;
     MRIfree(&mri_tmp) ;
  }
  fclose(fp) ;
//  MRIwrite(mri_aseg, "a.mgz");
//  MRIwrite(mri_brain, "b.mgz");
  mri_ctrl = find_control_points_always_in_WM(mri_aseg, NULL) ;
//  remove_gradient_outliers(mri_norm, mri_ctrl, mri_ctrl, max_grad) ;
  for (t = 0 ; t < ntps ; t++)
  {
    MRI *mri_tmp = NULL ;
    float scale ;
    
    mri_tmp = MRIcopyFrame(mri_norm, mri_tmp, t, 0) ;
    scale = MRImeanInLabel(mri_tmp, mri_ctrl, CONTROL_MARKED) ;
    printf("mean in wm is %2.0f, scaling by %2.2f\n", scale, 110/scale) ;
    scale = 110/scale ;
    MRIscalarMul(mri_tmp, mri_tmp, scale) ;
    MRIcopyFrame(mri_tmp, mri_norm, 0, t) ;
  }
  remove_points_not_in_range(mri_brain, mri_ctrl, mri_ctrl, DEFAULT_DESIRED_WHITE_MATTER_VALUE-intensity_pad, DEFAULT_DESIRED_WHITE_MATTER_VALUE+intensity_pad) ;

//  MRIwrite(mri_ctrl, "c.mgz") ;
//  remove_absolute_outliers(mri_norm, mri_ctrl, mri_ctrl, min_intensity) ;
  printf("using %d final control points\n", MRIcountNonzero(mri_ctrl)) ;
  for (t = 0 ; t < ntps ; t++)
  {
    MRI *mri_tmp = NULL ;
    
    mri_tmp = MRIcopyFrame(mri_norm, mri_tmp, t, 0) ;
    mri_bias = MRIbuildBiasImage(mri_tmp, mri_ctrl, NULL, bias_sigma) ;
    mri_dst = MRIapplyBiasCorrectionSameGeometry(mri_tmp, mri_bias, NULL, DEFAULT_DESIRED_WHITE_MATTER_VALUE) ;
    int req = snprintf(fname, STRLEN, "%s/%s.long.%s/mri/%s", sdir, tp_names[t], bname, out_fname);
    if (req >= STRLEN) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing output to %s\n", fname) ;
    MRIwrite(mri_dst, fname) ;
  }

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
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage_exit(0);
  }
  else if (!stricmp(option, "cross_time_sigma"))
  {
    cross_time_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing images with Parzen window with sigma=  = %2.2f\n", cross_time_sigma) ;

  }
  else if (!stricmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else switch (toupper(*option))
  {
  case 'P':
    intensity_pad = atoi(argv[2]) ;
    nargs = 1;
    printf("using intensity pad %d (default=2)\n", intensity_pad) ;
    break ;
  case 'S':
    bias_sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("smoothing bias field with sigma = %2.2f\n", bias_sigma) ;
    break ;
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
  case 'A':
    aseg_name = argv[1] ;
    nargs = 1 ;
    printf("using aseg %s for normalization\n", aseg_name) ;
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
    printf( "using control points from file %s...\n",
	    control_point_fname) ;
    break ;
  case '?':
  case 'U':
  case 'H':
    usage_exit(0) ;
  break ;
  default:
    printf( "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

#include "mri_long_normalize.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_long_normalize_help_xml,mri_long_normalize_help_xml_len);
  exit(code);
}


static MRI *
find_control_points_always_in_WM(MRI *mri_aseg, MRI *mri_ctrl)
{
  int  x, y, z, f, always_wm, label ;

  if (mri_ctrl == NULL)
  {
    mri_ctrl = MRIalloc(mri_aseg->width, mri_aseg->height, mri_aseg->depth,MRI_UCHAR) ;
    MRIcopyHeader(mri_aseg, mri_ctrl) ;
  }


  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
	always_wm = 1 ;
	for (f = 0 ; f < mri_aseg->nframes ; f++)
	{
	  label = MRIgetVoxVal(mri_aseg, x, y, z, f) ;
	  if (IS_WM(label) == 0)
	  {
	    always_wm = 0 ;
	    break ;
	  }
	}
	if (always_wm)
	  MRIsetVoxVal(mri_ctrl, x, y, z, 0, CONTROL_MARKED) ;
      }

  return(mri_ctrl) ;
}

static MRI *
remove_points_not_in_range(MRI *mri_brain, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, float min_target, float max_target) 
{
  int  x, y, z, f, removed, nremoved = 0 ;
  float val ;

  if (mri_ctrl_dst == NULL)
    mri_ctrl_dst = MRIclone(mri_ctrl_src, NULL) ;

  for (x = 0 ; x < mri_brain->width ; x++)
    for (y = 0 ; y < mri_brain->height ; y++)
      for (z = 0 ; z < mri_brain->depth ; z++)
      {
	if (MRIgetVoxVal(mri_ctrl_src, x, y, z, 0) != CONTROL_MARKED)
	  continue ;
	MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_MARKED) ;  // tentative
	removed = 0 ;
	for (f = 0 ; f < mri_brain->nframes ; f++)
	{
	  val = MRIgetVoxVal(mri_brain, x, y, z, f) ;
	  if ((val < min_target) || (val > max_target))
	  {
	    removed = 1 ;
	    break ;
	  }
	}
	if (removed)
	{
	  MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_NONE) ;  // erase it
	  nremoved++ ;
	}
      }
  printf("%d intensity outlier control points removed\n", nremoved) ;

  return(mri_ctrl_dst) ;
}
#if 0
static MRI *
remove_gradient_outliers(MRI *mri_norm, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, double max_grad) 
{
  int  x, y, z, f, xk, yk, zk, xi, yi, zi, removed, nremoved = 0 ;
  float val0, val ;

  if (mri_ctrl_dst == NULL)
    mri_ctrl_dst = MRIclone(mri_ctrl_src, NULL) ;


  for (x = 0 ; x < mri_norm->width ; x++)
    for (y = 0 ; y < mri_norm->height ; y++)
      for (z = 0 ; z < mri_norm->depth ; z++)
      {
	if (MRIgetVoxVal(mri_ctrl_src, x, y, z, 0) != CONTROL_MARKED)
	  continue ;
	MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_MARKED) ;  // tentative
	for (f = 0 ; f < mri_norm->nframes ; f++)
	{
	  val0 = MRIgetVoxVal(mri_norm, x, y, z, f) ;
	  removed = 0 ;
	  for (xk = -1 ; xk <= 1 ; xk++)
	  {
	    xi = mri_norm->xi[x+xk] ;
	    for (yk = -1 ; yk <= 1 ; yk++)
	    {
	      yi = mri_norm->yi[y+yk] ;
	      for (zk = -1 ; zk <= 1 ; zk++)
	      {
		zi = mri_norm->zi[z+zk] ;
		val = MRIgetVoxVal(mri_norm, xi, yi, zi, f) ;
		if (val0-val > max_grad)
		{
		  MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_NONE); // reverse decision
		  removed = 1 ;
		  break ;
		}
	      }
	      if (removed)
		break ;
	    }
	    if (removed)
	      break ;
	  }
	  nremoved += removed ;
	}
      }
  printf("%d gradient outlier control points removed\n", nremoved) ;

  return(mri_ctrl_dst) ;
}
static MRI *
remove_absolute_outliers(MRI *mri_norm, MRI *mri_ctrl_src, MRI *mri_ctrl_dst, double min_intensity) 
{
  int  x, y, z, f, removed, nremoved = 0 ;
  float val ;

  if (mri_ctrl_dst == NULL)
    mri_ctrl_dst = MRIclone(mri_ctrl_src, NULL) ;

  for (x = 0 ; x < mri_norm->width ; x++)
    for (y = 0 ; y < mri_norm->height ; y++)
      for (z = 0 ; z < mri_norm->depth ; z++)
      {
	if (MRIgetVoxVal(mri_ctrl_src, x, y, z, 0) != CONTROL_MARKED)
	  continue ;
	MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_MARKED) ;  // tentative
	removed = 0 ;
	for (f = 0 ; f < mri_norm->nframes ; f++)
	{
	  val = MRIgetVoxVal(mri_norm, x, y, z, f) ;
	  if (val < min_intensity)
	    removed = 1 ;
	}
	if (removed)
	{
	  MRIsetVoxVal(mri_ctrl_dst, x, y, z, 0, CONTROL_NONE) ; // erase it
	  nremoved++ ;
	}
      }
  printf("%d intensity outlier control points removed\n", nremoved) ;

  return(mri_ctrl_dst) ;
}

#endif
