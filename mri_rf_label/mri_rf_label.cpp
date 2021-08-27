/**
 * @brief random forest classification
 *
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
#include "rfa.h"
#include "mrisegment.h"

static int remove_wmsas_close_to_surface(char **surf_names, int nsurfs, MRI *mri_labeled, double surface_dist) ;
static int postprocess_segmentation_with_aseg(MRI *mri_labeled, MRI *mri_aseg, int min_voxels)  ;
static int postprocess_grow_wmsas(MRI *mri_labeled, MRI *mri_pvals, float pthresh, MRI *mri_aseg)  ;
#define MAX_READS 100
static int nreads = 0 ;
static char *read_intensity_fname[MAX_READS] ;
static int avgs = 0 ;

static int wmsa = 0 ;   // apply wmsa postprocessing (using T2/PD data)
static int nowmsa = 0 ; // remove all wmsa labels from the atlas

static double TRs[MAX_GCA_INPUTS] ;
static double fas[MAX_GCA_INPUTS] ;
static double TEs[MAX_GCA_INPUTS] ;
static int nsurfs = 0 ;
static char *surface_names[MAX_SURFACES] ;
static double surface_dist = 1 ;



int distance_to_label( MRI *mri_labeled, int label, int x,
                       int y, int z, int dx, int dy,
                       int dz, int max_dist );


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;


static double TR = -1 ;
static double alpha = -1 ;
static double TE = -1 ;
static char *mask_volume_fname = NULL ;

const char *Progname ;
static void usage_exit(int code) ;

static int filter = 0 ;

static int conform_flag = FALSE ;
static int single_classifier_flag = 0 ;
static char *only_nbrs_rf_fname = 0 ;   // only pick voxels that are on borders of a wmsa to train
static GCA *gca ;
static float  wmsa_thresh = 0.0 ;
static float pthresh = -1 ;
static float wm_thresh = .8 ;
static int wmsa_whalf = 3 ;
static MRI *relabel_wmsa_nbrs_with_random_forest(RANDOM_FOREST *rf, TRANSFORM *transform, GCA *gca, 
						 MRI *mri_inputs, MRI *mri_labeled) ;
static MRI *label_with_random_forest(RANDOM_FOREST *rf, TRANSFORM *transform, GCA *gca,
				     float wm_thresh, MRI *mri_in, MRI *mri_labeled, int wmsa_whalf,
				     MRI *mri_aseg, MRI **pmri_pvals) ;

static MRI *mri_aseg = NULL; 
static int min_voxels = 6 ; // remove wmsa segments smaller than this

int
main(int argc, char *argv[])
{
  char         **av ;
  int          ac, nargs, extra = 0 ;
  char         *in_fname, *out_fname,  *rfa_fname, *xform_fname ;
  MRI          *mri_inputs, *mri_labeled, *mri_tmp, *mri_pvals ;
  int          msec, minutes, seconds, ninputs, input ;
  Timer start ;
  TRANSFORM     *transform ;
  RFA          *rfa = NULL ;
  RANDOM_FOREST *rf = NULL ;  // if single_classifier_flag is true


  std::string cmdline = getAllInfo(argc, argv, "mri_rf_label");

  nargs = handleVersionOption(argc, argv, "mri_rf_label");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

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

  if (argc < 5)
  {
    usage_exit(1) ;
  }

  in_fname = argv[1] ;
  xform_fname = argv[argc-3];
  rfa_fname = argv[argc-2] ;
  out_fname = argv[argc-1] ;
  ninputs = argc-4 ;

  printf("reading %d input volumes...\n", ninputs) ;

  /*  fprintf(stderr,
      "mri_inputs read: xform %s\n", mri_inputs->transform_fname) ;*/
  if (single_classifier_flag)
  {
    rf = RFread(rfa_fname) ;
  }
  else
  {
    printf("reading random forest classifier array from %s...\n", rfa_fname) ;
    rfa = RFAread(rfa_fname) ;
    if (!rfa)
      ErrorExit(ERROR_NOFILE, "%s: could not read random forest classifier array from %s",
		Progname, rfa_fname) ;
    
    if (rfa->ninputs != (ninputs+extra+3))  // +3 is for spatial coords or priors
      ErrorExit
	(ERROR_BADPARM,
	 "%s: rfa requires %d inputs, %d specified on command line",
	 Progname, rfa->ninputs, ninputs) ;
  }

  // gathering inputs
  for (input = 0 ; input < ninputs ; input++)
  {
    in_fname = argv[1+input] ;
    printf("reading input volume from %s...\n", in_fname) ;
    mri_tmp = MRIread(in_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, in_fname) ;

    if (alpha > 0)
      mri_tmp->flip_angle = alpha ;
    if (TR > 0)
      mri_tmp->tr = TR ;
    if (TE > 0)
      mri_tmp->te = TE ;

    TRs[input] = mri_tmp->tr ; fas[input] = mri_tmp->flip_angle ; TEs[input] = mri_tmp->te ;
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
      mri_dir = MRIgradientDir2ndDerivative(mri_tmp, NULL, 5) ;
      MRIscaleAndMultiply(mri_tmp, 128.0, mri_dir, mri_tmp) ;
      MRIwrite(mri_dir, "lap.mgz") ;
      MRIwrite(mri_tmp, "filtered.mgz") ;
      MRIfree(&mri_dir) ;
      MRIfree(&mri_smooth) ;
      exit(1) ;
    }
    MRIcopyFrame(mri_tmp, mri_inputs, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  MRIaddCommandLine(mri_inputs, cmdline) ;


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

#if 0
  if (Ggca_x >= 0 && Gx < 0)
  {
    RFAsourceVoxelToNode(rfa, mri_inputs, transform, Ggca_x, Ggca_y, Ggca_z,&Gx, &Gy, &Gz) ;
    printf("source voxel (%d, %d, %d) maps to node (%d, %d, %d)\n", Ggca_x, Ggca_y, Ggca_z, Gx, Gy, Gz) ;
  }
#endif

  printf("labeling volume...\n") ;
  // create labeled volume using array of random forests
  if (single_classifier_flag)
  {
    mri_labeled = label_with_random_forest(rf, transform, gca, wm_thresh, mri_inputs, NULL, wmsa_whalf,
					   mri_aseg, &mri_pvals) ;
    if (only_nbrs_rf_fname)
    {
      RANDOM_FOREST *rf_nbrs = RFread(only_nbrs_rf_fname) ;
      if (rf_nbrs == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read nbr random forest from %s",
		  Progname, only_nbrs_rf_fname) ;
      relabel_wmsa_nbrs_with_random_forest(rf_nbrs, transform, gca, mri_inputs, mri_labeled) ;
    }
    if (mri_aseg)
      postprocess_segmentation_with_aseg(mri_labeled, mri_aseg, min_voxels) ;
    if (nsurfs > 0)
	remove_wmsas_close_to_surface(surface_names, nsurfs, mri_labeled, surface_dist) ;
    if (pthresh >= 0)
      postprocess_grow_wmsas(mri_labeled, mri_pvals, pthresh, mri_aseg) ;
  }
  else
    mri_labeled = RFAlabel(mri_inputs, rfa, NULL, transform) ;

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


#if 0  
  mri_labeled->ct = rfa->ct ;  // embed color table in output volume
#endif

  MRIfree(&mri_inputs) ;

  printf("writing labeled volume to %s...\n", out_fname) ;
  if (MRIwrite(mri_labeled, out_fname) != NO_ERROR)
  {
    ErrorExit(Gerror, "%s: MRIwrite(%s) failed", Progname, out_fname) ;
  }

  msec = start.milliseconds() ;
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
  if (!stricmp(option, "aseg"))
  {
    printf("reading aseg  from %s\n", argv[2]) ;
    mri_aseg = MRIread(argv[2]) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load aseg from %s", Progname, argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "NBRS"))
  {
    only_nbrs_rf_fname = argv[2] ;
    printf("reading random forest from %s and applying to WMSA nbrs\n", only_nbrs_rf_fname) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "MIN_VOXELS"))
  {
    min_voxels = atoi(argv[2]) ;
    printf("removing WMSA segments that have fewer than %d voxels\n", min_voxels) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "surface"))
  {
    surface_names[nsurfs++] = argv[2] ;
    printf("removing WMSA segments that are more than %2.1f mm interior to the ?h.%s surfaces\n", 
	   surface_dist, argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "surface_dist"))
  {
    surface_dist = atof(argv[2]) ;
    printf("removing WMSA segments that are more than %2.1f mm interior surfaces\n", surface_dist) ;
    nargs = 1 ;
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
  else if (!stricmp(option, "-HELP")||!stricmp(option, "-USAGE"))
  {
    usage_exit(0) ;
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
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
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
  else if (!stricmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
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
  else if (!stricmp(option, "WMSA_WHALF"))
  {
    nargs = 1 ;
    wmsa_whalf = atoi(argv[2]) ;
    printf("only examing voxels that wmsa occurred within %d voxels of\n", wmsa_whalf) ;
  }
  else if (!stricmp(option, "THRESH"))
  {
    nargs = 1 ;    
    wmsa_thresh = atof(argv[2]) ;
    printf("only labeling WMSAs that exceed threshold %2.2f\n", wmsa_thresh) ;
  }
  else if (!stricmp(option, "PTHRESH"))
  {
    nargs = 1 ;    
    pthresh = atof(argv[2]) ;
    printf("relabeling voxels adjacent to WMSA that have p-vals < %2.3f\n", pthresh) ;
  }
  else switch (toupper(*option))
    {
    case '1':
      single_classifier_flag = 1 ;
      gca = GCAread(argv[2]) ;
      if (gca == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read gca from %s", argv[2]) ;
      nargs = 1 ;
      printf("training a single classifier instead of an array using gca %s\n", argv[2]) ;
      break ;
    case 'T':
      wm_thresh = atof(argv[2]) ;
      nargs = 1 ;
      printf("thresholding wm priors at %f to build training set\n", wm_thresh) ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'A':
      avgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf
      (stderr,
       "applying mean filter %d times to conditional densities...\n", avgs) ;
      break ;
    case 'M':
      mask_volume_fname = argv[2] ;
      nargs = 1 ;
      printf("using %s to mask final labeling...\n", mask_volume_fname) ;
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
#include "mri_rf_label.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_rf_label_help_xml, mri_rf_label_help_xml_len);
  exit(code);
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


static MRI *
label_with_random_forest(RANDOM_FOREST *rf, TRANSFORM *transform, GCA *gca,
			 float wm_thresh, MRI *mri_in, MRI *mri_labeled, int wmsa_whalf, MRI *mri_aseg, 
			 MRI **pmri_pvals)
{
  int     x, y, z, wsize, label ;
  double  *feature, xatlas, yatlas, zatlas, pval ;
  MRI     *mri_wmsa_possible, *mri_pvals ;

  wsize = nint(pow((rf->nfeatures-3)/mri_in->nframes, 1.0/3)) ;
  if (mri_labeled == NULL)
  {
    mri_labeled = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT) ;
    MRIcopyHeader(mri_in, mri_labeled) ;
  }
  mri_pvals = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 3) ;
  MRIcopyHeader(mri_in, mri_pvals) ;
  mri_wmsa_possible = MRIclone(mri_labeled, NULL) ;
  for (x = 0 ; x < mri_in->width ; x++)
    for (y = 0 ; y < mri_in->height ; y++)
      for (z = 0 ; z < mri_in->depth ; z++)
	if (is_possible_wmsa(gca, mri_in, transform, x, y, z, 0))
	  MRIsetVoxVal(mri_wmsa_possible, x, y, z, 0, 1) ;
  for ( ; wmsa_whalf > 0 ; wmsa_whalf--)
    MRIdilate(mri_wmsa_possible, mri_wmsa_possible) ;

  feature = (double *)calloc(rf->nfeatures, sizeof(double)) ;
  if (feature == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d-len feature vector", rf->nfeatures) ;

  if (Gx >= 0)    // diagnostics
  {
    int whalf = (wsize-1)/2, n, i ;
    char buf[STRLEN] ;

    rf->feature_names = (char **)calloc(rf->nfeatures, sizeof(char *)) ;
    for (i = 0, x = -whalf ; x <= whalf ; x++)
      for (y = -whalf ; y <= whalf ; y++)
	for (z = -whalf ; z <= whalf ; z++)
	  for (n = 0 ; n < mri_in->nframes ; n++, i++)
	  {
	    switch (n)
	    {
	    default:
	    case 0: sprintf(buf, "T1(%d, %d, %d)", x, y, z) ; break ;
	    case 1: sprintf(buf, "T2(%d, %d, %d)", x, y, z) ; break ;
	    case 2: sprintf(buf, "PD(%d, %d, %d)", x, y, z) ; break ;
	    }
	    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
	    strcpy(rf->feature_names[i], buf) ;
	  }
    rf->feature_names[i++] = const_cast<char*>("nonzero count");
    rf->feature_names[i++] = const_cast<char*>("gm prior");
    rf->feature_names[i++] = const_cast<char*>("wm prior");
    rf->feature_names[i++] = const_cast<char*>("csf prior");
  }

  for (x = 0 ; x < mri_in->width ; x++)
    for (y = 0 ; y < mri_in->height ; y++)
      for (z = 0 ; z < mri_in->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	{
	  printf("voxel (%d, %d, %d); WM prior = %2.2f, cortex prior = %f, wmsa possible = %d\n",
		 x,y,z, 
		 wm_prior(gca, mri_in, transform, x, y, z),
		 cortex_prior(gca, mri_in, transform, x, y, z),
		 (int)MRIgetVoxVal(mri_wmsa_possible,x,y,z,0));
	  DiagBreak() ;
	}
	MRIsetVoxVal(mri_pvals, x, y, z, 0, wm_prior(gca, mri_in, transform, x, y, z)) ;
	MRIsetVoxVal(mri_pvals, x, y, z, 1, cortex_prior(gca, mri_in, transform, x, y, z)) ;
	MRIsetVoxVal(mri_pvals, x, y, z, 2, 1) ; // will be changed later if admissable voxel
	if (((int)MRIgetVoxVal(mri_wmsa_possible,x,y,z,0) == 0) &&
	    ((wm_prior(gca, mri_in, transform, x, y, z) < wm_thresh) ||
	     (cortex_prior(gca, mri_in, transform, x, y, z) > .5)))
	{
	  if (x == Gx && y == Gy && z == Gz)
	    printf("voxel (%d, %d, %d); WM prior = %2.2f, cortex prior = %f\n",
		   x,y,z, 
		   wm_prior(gca, mri_in, transform, x, y, z),
		   cortex_prior(gca, mri_in, transform, x, y, z)) ;
	  continue ;
	}
	TransformSourceVoxelToAtlas(transform, mri_in, x, y, z, &xatlas, &yatlas, &zatlas) ;
	extract_feature(mri_in, wsize, x, y, z, feature, xatlas, yatlas, zatlas) ;
	if (mri_aseg)
	  feature[rf->nfeatures-4] = MRIcountCSFInNbhd(mri_aseg, 5, x, y, z) ;
	else
	  feature[rf->nfeatures-4] = 0 ;
	feature[rf->nfeatures-3] = 100*gm_prior(gca, mri_in, transform, x, y, z) ;
	feature[rf->nfeatures-2] = 100*wm_prior(gca, mri_in, transform, x, y, z) ;
	feature[rf->nfeatures-1] = 100*csf_prior(gca, mri_in, transform, x, y, z) ;
	if (x == Gx && y == Gy && z == Gz)
	{
	  int j ;
	  printf("voxel (%d, %d, %d) ", x, y, z) ;
	  printf("\nfeature = ") ;
	  for (j = 0 ; j < rf->nfeatures ; j++)
	    printf("%.0f ", feature[j]) ;
	  printf("\n") ;
	  Gdiag |= DIAG_VERBOSE ;
	  DiagBreak() ;
	}
	label = RFclassify(rf, feature, &pval, -1);
	if (label > 0)
	  DiagBreak() ;
	if (label > 0 && pval > wmsa_thresh)
	  MRIsetVoxVal(mri_labeled, x, y, z, 0, label) ;
	MRIsetVoxVal(mri_pvals, x, y, z, 2, pval) ;
	if (x == Gx && y == Gy && z == Gz)
	  Gdiag &= ~DIAG_VERBOSE ;
      }

  if (Gx >= 0)
  {
    printf("writing pvals.mgz...\n") ;
    MRIwrite(mri_pvals, "pvals.mgz") ;
  }
  if (pmri_pvals)
    *pmri_pvals = mri_pvals ;
  else
    MRIfree(&mri_pvals) ;
  free(feature) ; MRIfree(&mri_wmsa_possible) ;
  return(mri_labeled) ;
}
static MRI *
relabel_wmsa_nbrs_with_random_forest(RANDOM_FOREST *rf, TRANSFORM *transform, GCA *gca, 
				     MRI *mri_in, MRI *mri_labeled)
{
  int     x, y, z, wsize, label, non, noff, nunchanged = 0, new_label, total ;
  double  *feature, pval ;
  MRI     *mri_mask, *mri_orig_labeled ;

  wsize = nint(pow((rf->nfeatures-3)/mri_in->nframes, 1.0/3)) ;
  if (mri_labeled == NULL)
  {
    mri_labeled = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_SHORT) ;
    MRIcopyHeader(mri_in, mri_labeled) ;
  }
  mri_orig_labeled = MRIcopy(mri_labeled, NULL) ;
  mri_mask = MRIcopy(mri_labeled, NULL) ;
  MRIdilate(mri_mask, mri_mask) ;  // nbrs

  feature = (double *)calloc(rf->nfeatures, sizeof(double)) ;
  if (feature == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d-len feature vector", rf->nfeatures) ;

  if (Gx >= 0)    // diagnostics
  {
    int whalf = (wsize-1)/2, n, i ;
    char buf[STRLEN] ;

    rf->feature_names = (char **)calloc(rf->nfeatures, sizeof(char *)) ;
    for (i = 0, x = -whalf ; x <= whalf ; x++)
      for (y = -whalf ; y <= whalf ; y++)
	for (z = -whalf ; z <= whalf ; z++)
	  for (n = 0 ; n < mri_in->nframes ; n++, i++)
	  {
	    switch (n)
	    {
	    default:
	    case 0: sprintf(buf, "T1(%d, %d, %d)", x, y, z) ; break ;
	    case 1: sprintf(buf, "T2(%d, %d, %d)", x, y, z) ; break ;
	    case 2: sprintf(buf, "PD(%d, %d, %d)", x, y, z) ; break ;
	    }
	    rf->feature_names[i] = (char *)calloc(strlen(buf)+1, sizeof(char)) ;
	    strcpy(rf->feature_names[i], buf) ;
	  }
    rf->feature_names[i++] = const_cast<char*>("gm prior");
    rf->feature_names[i++] = const_cast<char*>("wm prior");
    rf->feature_names[i++] = const_cast<char*>("csf prior");
  }

  for (non = noff = x = 0 ; x < mri_in->width ; x++)
    for (y = 0 ; y < mri_in->height ; y++)
      for (z = 0 ; z < mri_in->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	{
	  printf("voxel (%d, %d, %d); WM prior = %2.2f, cortex prior = %f, mask = %d\n",
		 x,y,z, 
		 wm_prior(gca, mri_in, transform, x, y, z),
		 cortex_prior(gca, mri_in, transform, x, y, z),
		 (int)MRIgetVoxVal(mri_mask,x,y,z,0));
	  DiagBreak() ;
	}
	label = MRIgetVoxVal(mri_orig_labeled, x, y, z, 0) ;
	if ((IS_WMSA(label) && MRIneighborsOff(mri_orig_labeled,x,y,z,1) == 0) ||
	    MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)
	  continue ;
	
	extract_feature(mri_in, wsize, x, y, z, feature, 0, 0, 0) ;
	feature[rf->nfeatures-3] = 100*gm_prior(gca, mri_in, transform, x, y, z) ;
	feature[rf->nfeatures-2] = 100*wm_prior(gca, mri_in, transform, x, y, z) ;
	feature[rf->nfeatures-1] = 100*csf_prior(gca, mri_in, transform, x, y, z) ;
	if (x == Gx && y == Gy && z == Gz)
	{
	  int j ;
	  printf("voxel (%d, %d, %d) ", x, y, z) ;
	  printf("\nfeature = ") ;
	  for (j = 0 ; j < rf->nfeatures ; j++)
	    printf("%.0f ", feature[j]) ;
	  printf("\n") ;
	  Gdiag |= DIAG_VERBOSE ;
	  DiagBreak() ;
	}
	new_label = RFclassify(rf, feature, &pval, -1);
	if (x == Gx && y == Gy && z == Gz)
	  printf("\nlabel = %d, pval = %2.2f\n", new_label, pval) ;
	if (new_label > 0)
	  DiagBreak() ;
	MRIsetVoxVal(mri_labeled, x, y, z, 0, new_label) ;
	if (new_label == label)
	  nunchanged++ ;
	else if (new_label > 0)
	  non++ ;
	else
	  noff++ ;
	if (x == Gx && y == Gy && z == Gz)
	  Gdiag &= ~DIAG_VERBOSE ;
      }

  total = non+noff+nunchanged ;
  printf("%2.0f%% (%d) WMSA added and %2.0f%% (%d) removed, %2.0f%% (%d)\n", 
	 100.0*(float)non/total, non, 
	 100.0*(float)noff/total, noff, 
	 100.0*(float)nunchanged/total, nunchanged) ; 
  free(feature) ; MRIfree(&mri_mask) ; MRIfree(&mri_orig_labeled) ;
  return(mri_labeled) ;
}

static int close_labels[] = 
{ 
  CSF,
  Left_Lateral_Ventricle,
  Right_Lateral_Ventricle,
  Third_Ventricle,
  Fourth_Ventricle,
  CC_Posterior,
  CC_Mid_Posterior,
  CC_Central,
  CC_Mid_Anterior,
  CC_Anterior,
  Left_Inf_Lat_Vent, 
  Right_Inf_Lat_Vent, 
  Left_choroid_plexus,
  Right_choroid_plexus,
  Unknown, 


  Left_Amygdala,
  Right_Amygdala,
  Left_Putamen,
  Right_Putamen,
  Left_Pallidum,
  Right_Pallidum
} ;

#define NCLOSE_LABELS (sizeof(close_labels) / sizeof(close_labels[0]))

static int
postprocess_segmentation_with_aseg(MRI *mri_labeled, MRI *mri_aseg, int min_voxels) 
{
  int  x, y, z, aseg_label, in_wmsa_label, out_wmsa_label, l, nchanged = 0 ;
  MRI_SEGMENTATION *mseg ;
  MRI              *mri_closed ;

  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	in_wmsa_label = MRIgetVoxVal(mri_labeled, x, y, z, 0);
	if (in_wmsa_label == 0) // not labeled WMSA
	  continue ;
	aseg_label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	out_wmsa_label = in_wmsa_label ;
	if (NOT_TRAINING_LABEL(aseg_label))
	  out_wmsa_label = 0 ;
	switch (aseg_label)
	{
	case Unknown:
	case Right_Putamen:
	case Left_Putamen:
	case Right_Thalamus:
	case Left_Thalamus:
	case Left_Lateral_Ventricle:
	case Right_Lateral_Ventricle:
	case Left_Inf_Lat_Vent:
	case Right_Inf_Lat_Vent:
	  out_wmsa_label = 0 ;
	  break ;
	default:
	  break ;
	}
	MRIsetVoxVal(mri_labeled, x, y, z, 0, out_wmsa_label) ;
	if (aseg_label != WM_hypointensities)
	  for (l = 0 ; l < NCLOSE_LABELS ; l++)
	  {
	    if (MRIlabelsInNbhd(mri_aseg, x, y, z, 1, close_labels[l]) > 0)
	      out_wmsa_label = 0 ;
	  }
	if (out_wmsa_label != in_wmsa_label)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    printf("voxel (%d, %d, %d): WMSA label removed during postprocessing\n", x,y,z) ;
	  nchanged++ ;
	}
	MRIsetVoxVal(mri_labeled, x, y, z, 0, out_wmsa_label) ;
      }
  mri_closed = MRIclose(mri_labeled, NULL) ;
  mseg = MRIsegment(mri_labeled, .5, 1000) ;
  MRIfree(&mri_closed) ;
  nchanged += MRIeraseSmallSegments(mseg, mri_labeled, min_voxels) ;
  MRIsegmentFree(&mseg) ;
  printf("%d labels wmsa removed during aseg postprocessing\n", nchanged) ;
  return(NO_ERROR) ;
}

static int
remove_wmsas_close_to_surface(char **surface_names, int nsurfs, MRI *mri_labeled, double surface_dist)
{
  MRI         *mri_dist, *mri_dist_total ;
  int         x, y, z, nchanged, i ;
  MRI_SURFACE *mris ;
  
  for (i = 0 ; i < nsurfs ; i++)
  {
    mris = MRISread(surface_names[i]) ;
    if (mris == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, surface_names[i]) ;

    mri_dist = MRIcloneDifferentType(mri_labeled, MRI_FLOAT) ;
    MRIScomputeDistanceToSurface(mris, mri_dist, mri_dist->xsize) ;

    if (i == 0)
      mri_dist_total = mri_dist ;
    else
    {
      MRImin(mri_dist_total, mri_dist, mri_dist_total) ;
      MRIfree(&mri_dist) ;
    }
    MRISfree(&mris) ;
  }

  for (nchanged = x = 0 ; x < mri_labeled->width ;  x++)
    for (y = 0 ; y < mri_labeled->height ;  y++)
      for (z = 0 ; z < mri_labeled->depth ;  z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	if (MRIgetVoxVal(mri_labeled, x, y, z, 0) == 0)
	  continue ;  // not a wmsa
	if (MRIgetVoxVal(mri_dist_total, x, y, z, 0) > surface_dist)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    printf("removing voxel (%d, %d, %d): surface_dist = %2.3f\n",x,y,z,MRIgetVoxVal(mri_dist, x, y, z, 0));
	  nchanged++ ;
	  MRIsetVoxVal(mri_labeled, x, y, z, 0, 0) ;
	}
      }

  printf("%d voxels removed due to surface proximity\n", nchanged) ;
  MRIfree(&mri_dist_total) ;
  return(nchanged) ;
}

static int
postprocess_grow_wmsas(MRI *mri_labeled, MRI *mri_pvals, float pthresh, MRI *mri_aseg) 
{
  int   x, y, z, nchanged, total_changed, label ;
  float pval ;
  MRI   *mri_tmp ;

  total_changed = 0 ;
  mri_tmp = MRIcopy(mri_labeled, NULL) ;
  do
  {
    nchanged = 0 ;
    for (x = 0 ;  x < mri_pvals->width; x++)
      for (y = 0 ;  y < mri_pvals->height; y++)
	for (z = 0 ;  z < mri_pvals->depth; z++)
	{
	  if (Gx == x && Gy == y && Gz == z)
	    DiagBreak() ;
	  label = MRIgetVoxVal(mri_labeled, x, y, z, 0) ;
	  if (label)  
	    continue ;   // already a wmsa
	  pval = MRIgetVoxVal(mri_pvals, x, y, z, 2) ;
	  if (pval > pthresh)
	    continue ;    // too confident in the non-wmsa labeling
	  label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
	  if (NOT_TRAINING_LABEL(label))
	    continue ;

	  if (MRIcountNonzeroInNbhd(mri_labeled, 3, x, y, z) > 0)
	  {
	    if (Gx == x && Gy == y && Gz == z)
	      printf("voxel (%d, %d, %d): changed to WMSA\n", x,y,z) ;
	    MRIsetVoxVal(mri_tmp, x, y, z, 0, 1) ;
	    nchanged++ ;
	  }
	}
    MRIcopy(mri_tmp, mri_labeled) ;
    total_changed += nchanged ;
    printf("%d voxels changed to wmsa in region growing\n", nchanged) ;
  } while (nchanged > 0) ;
  MRIfree(&mri_tmp) ;
  return(total_changed) ;
}

