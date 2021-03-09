/**
 * @brief program for training voxel-based classifiers
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
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "mriclass.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "const.h"
#include "classify.h"
#include "version.h"
#include "rforest.h"
#include "cma.h"
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

static int features = FEATURE_INTENSITY | FEATURE_MEAN3 | FEATURE_DIRECTION |
                      FEATURE_CPOLV_MEDIAN5 ;

static int extract = 0 ;
static int classifier = CLASSIFIER_RFOREST ;
static char priors_fname[100] = "none" ;
static int  verbose = 0 ;

const char *Progname ;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

#define NCLUSTERS  6

static int nclusters = 0 ;
static int train_cpolv = 0 ;

static RBF_PARMS rbf_parms = {
                               {
                                 NCLUSTERS/2, NCLUSTERS, NCLUSTERS, NCLUSTERS, NCLUSTERS, NCLUSTERS/2
                               }
                             } ;

#define MAX_LONG 1000   // max # of timepoints
#define MAX_VOLS 10     // max # of input volumes

static int wsize = 1 ;  // use intensities in a 3x3x3 window around each point
static int nlong = 0 ;
static char *long_names[MAX_LONG] ;
static double long_times[MAX_LONG] ;

static int nvols = 0 ;
static const char *vol_names[MAX_VOLS] = 
{
  "norm.mgz"
} ;

static const char *seg_name = "wmsa/wmsa.mgz" ;

static char *classify_name = NULL ;
static const char *aseg_name = "aseg.mgz" ;

static char sdir[STRLEN] = "" ;
static RF_PARMS rf_parms ;
static const char *wmsa_class_names[] =
{
  "NOT WM",
  "WMSA 0",
  "WMSA 1",
  "WMSA 2",
  "WMSA 3",
  "WMSA 4",
  "WMSA 5",
  "WMSA 6",
  "WMSA 7",
} ;
#if 0
static int nwmsa_classes = sizeof(wmsa_class_names) / sizeof(wmsa_class_names[0]);
#endif

static RANDOM_FOREST *train_rforest(const char *training_file_name, RF_PARMS *parms, int wsize, int nlong, char **long_names, int nvols, const char **vol_names, const char *sdir, const char *seg_name, double *paccuracy) ;

static int classify_subjects(RANDOM_FOREST *rf, const char *training_file_name, int wsize, int nlong, 
				char **long_names, int nvols, const char **vol_names, const char *sdir, 
				const char *seg_name, const char *classify_name);


int
main(int argc, char *argv[]) {
  MRIC    *mric ;
  RANDOM_FOREST *rf ;
  char    *training_file_name, *output_file_name, *cp ;
  int     nargs, error, i ;
  double    accuracy =0 ;

  nargs = handleVersionOption(argc, argv, "mri_train");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  rf_parms.ntrees = 500 ;
  rf_parms.nsteps = 10 ;
  rf_parms.training_fraction = .001 ;
  rf_parms.feature_fraction = 1 ;
  rf_parms.max_depth = 10 ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,"usage: %s <training file name>  <output file>",
              Progname);

  if (strlen(sdir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (cp == NULL)
      ErrorExit(ERROR_UNSUPPORTED, "SUBJECTS_DIR must be in env or on cmdline with -sdir") ;
    exit(1) ;
    strcpy(sdir, cp) ;
  }

  training_file_name = argv[1] ;
  output_file_name = argv[2] ;

  switch (classifier)
  {
  case CLASSIFIER_RFOREST:
    rf = train_rforest(training_file_name, &rf_parms, wsize, nlong, long_names, nvols, vol_names, sdir, seg_name, &accuracy) ;
    {
      struct flock fl;
      int    fd;
      char   line[MAX_LINE_LEN] ;

      printf("writing results to train.log file\n") ;
      fd = open("train.log", O_WRONLY|O_APPEND|O_CREAT, S_IRWXU|S_IRWXG);
      if (fd < 0)
	ErrorExit(ERROR_NOFILE, "%s: could not open test log file", Progname);
      
      fcntl(fd, F_SETLKW, &fl);  /* F_GETLK, F_SETLK, F_SETLKW */
      sprintf(line, "%f %d %d %f\n", 
	      rf->training_fraction,
	      rf->max_depth,
	      rf->ntrees,
	      accuracy) ;
      if (write(fd, line, (strlen(line))*sizeof(char)) < (strlen(line))*sizeof(char))
        ErrorPrintf(ERROR_NOFILE, "%s: could not write %d bytes to %s", Progname, (strlen(line))*sizeof(char));
      fl.l_type   = F_UNLCK;  /* tell it to unlock the region */
      fcntl(fd, F_SETLK, &fl); /* set the region to unlocked */
      close(fd) ;
      if (classify_name)
	classify_subjects(rf, training_file_name, wsize, nlong, long_names, nvols, vol_names, sdir, seg_name, classify_name) ;
    }
    break ;
  default:
    if (nclusters > 0) {
      for (i = 0 ; i < NCLASSES ; i++) {
	if (ISWHITE(i) || i == GRAY_MATTER)
	  rbf_parms.max_clusters[i] = nclusters ;
	else
        rbf_parms.max_clusters[i] = nclusters/2 ;
      }
    } else
    nclusters = NCLUSTERS ;
    
    if (train_cpolv) {
      for (i = 0 ; i < NCLASSES ; i++) {
	if (ISWHITE(i))
	  rbf_parms.max_clusters[i] = nclusters/3 ;
      else
        if (i == CSF)
          rbf_parms.max_clusters[i] = nclusters ;
        else
          rbf_parms.max_clusters[i] = 0 ;
      }
    }
    
    mric = MRICalloc(1, &classifier, &features, (void *)&rbf_parms) ;
    if ((strlen(priors_fname) > 1) && stricmp(priors_fname, "none"))
      error = MRICtrain(mric, training_file_name, priors_fname) ;
    else
      error = MRICtrain(mric, training_file_name, NULL) ;
    
    if (error != NO_ERROR)
      ErrorExit(error, "training failed.\n") ;
    
    MRICwrite(mric, output_file_name) ;
    MRICfree(&mric) ;
    break ;
  }
  
  exit(0) ;
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "cpolv"))
    train_cpolv = 1 ;
  else if (!stricmp(option, "long")) {
    if (nlong >= MAX_LONG)
      ErrorExit(ERROR_NOMEMORY, "%s: too many long names specified (max=%d)\n",
		nlong) ;
    long_names[nlong] = argv[2] ;
    long_times[nlong] = atof(argv[3]) ;
    nargs = 2 ;
    printf("using longitudinal timepoint prefix %s at time %2.1f\n", 
	   long_names[nlong], long_times[nlong]) ;
    nlong++ ;
  } else if (!stricmp(option, "vol")) {
    if (nlong >= MAX_VOLS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many vol names specified (max=%d)\n",
		nvols) ;
    vol_names[nvols] = argv[2] ;
    nargs = 1 ;
    printf("using input volume named %s\n", vol_names[nvols]) ;
    nvols++ ;
  } else if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    nargs = 1 ;
    printf("using SUBJECTS_DIR = %s\n", sdir);
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "classify")) {
    classify_name = argv[2] ;
    nargs = 1 ;
    printf("writing out of bag classifications to %s\n", classify_name);
  } else if (!stricmp(option, "ntrees")) {
    rf_parms.ntrees = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using ntrees = %d\n", rf_parms.ntrees) ;
  } else if (!stricmp(option, "seg")) {
    seg_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as training segmentation volume\n", seg_name) ;
  } else if (!stricmp(option, "training_fraction")) {
    rf_parms.training_fraction = atof(argv[2]) ;
    nargs = 1 ;
    printf("using training = %lf\n", rf_parms.training_fraction) ;
  } else if (!stricmp(option, "feature_fraction")) {
    rf_parms.feature_fraction = atof(argv[2]) ;
    nargs = 1 ;
    printf("using feature fraction = %lf\n", rf_parms.feature_fraction) ;
  } else if (!stricmp(option, "max_depth")) {
    rf_parms.max_depth = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using max depth = %d\n", rf_parms.max_depth) ;
  } else switch (toupper(*option)) {
    case 'V':
      verbose = !verbose ;
      break ;
    case 'N':
      if (sscanf(argv[2], "%d", &nclusters) != 1)
        ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                  Progname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      wsize = atoi(argv[2]) ;
      nargs = 1 ;
      printf("using window size = %d\n", wsize) ;
      break ;
    case 'F':
      if (sscanf(argv[2], "0x%x", &features) != 1)
        ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                  Progname, argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "using features 0x%x\n", features) ;
      break ;
    case 'P':
      strcpy(priors_fname, argv[2]) ;
      nargs = 1 ;
      if (verbose)
        fprintf(stderr, "using priors file %s\n", priors_fname) ;
      break ;
    case 'X':
      if (sscanf(argv[2], "%d", &extract) != 1)
        ErrorExit(ERROR_BADPARM, "%s: could not scan option from '%s'",
                  Progname, argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      printf("usage: %s <training file> <output file>\n", Progname) ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

#define NOT_WMSA     1
#define WMSA         2
#define FUTURE_WMSA  3

#if 0
#define WM       2
#define CAUDATE  3
static int aseg_labels[] = 
{
  Left_WM_hypointensities,  
  Right_WM_hypointensities,  
  Left_Cerebral_White_Matter,
  Right_Cerebral_White_Matter,
  Left_Caudate,
  Right_Caudate
} ;

static int training_labels[] = 
{
  WMSA,
  WMSA,
  WM,
  WM,
  CAUDATE,
  CAUDATE
} ;
static int naseg = sizeof(aseg_labels) / sizeof(aseg_labels[0]) ;
#endif

MRI *
label_voxels_with_wmsa(MRI **mri_long_seg, int nlong, int not_wmsa_label, int future_wmsa_label, int wmsa_label)
{
  MRI   *mri_seg ;
  int   x, y, z, l, label, nhypo ;

  mri_seg = MRIalloc(mri_long_seg[0]->width, mri_long_seg[0]->height, mri_long_seg[0]->depth,
		     mri_long_seg[0]->type) ;
  MRIcopyHeader(mri_long_seg[0], mri_seg) ;

  for (x = 0 ; x < mri_seg->width ; x++)
    for (y = 0 ; y < mri_seg->height ; y++)
      for (z = 0 ; z < mri_seg->depth ; z++)
      {
	if (x == Gx && y == Gy && z == Gz)
	  DiagBreak() ;
	label = MRIgetVoxVal(mri_long_seg[nlong-1], x, y, z, 0) ;

	// hmmm, what about labels that are wm at end but hypo in between?
	if (IS_HYPO(label) == 0) // only consider voxels that are WSMA at the end
	  continue ;
	for (nhypo = l = 0 ; l < nlong-1 ; l++) // make sure they weren't wmsa before the end
	{
	  label = MRIgetVoxVal(mri_long_seg[l], x, y, z, 0) ;
	  if (IS_HYPO(label))
	    nhypo++ ;
	  else if (!IS_WHITE_MATTER(label))  // must be either hypo or wm at every time point
	    break ;  
	}

	if (l < nlong-1)  // had some non wm and non-hypo label
	  continue ;
	if (nhypo > 0)
	  MRIsetVoxVal(mri_seg, x, y, z, 0, wmsa_label) ;
	else
	  MRIsetVoxVal(mri_seg, x, y, z, 0, future_wmsa_label) ;
      }

  return(mri_seg) ;
}

static MRI *
label_wmsa_and_future_wmsa(MRI *mri_mask_src, MRI *mri_mask_dst, MRI **mri_long_seg, int nlong,
			   int not_wmsa_label, int future_wmsa_label, int wmsa_label)
{
  int   x, y, z, l, label ;

  mri_mask_dst = MRIcopy(mri_mask_src, mri_mask_dst) ;

  for (x = 0 ; x < mri_mask_src->width ; x++)
    for (y = 0 ; y < mri_mask_src->height ; y++)
      for (z = 0 ; z < mri_mask_src->depth ; z++)
      {
	if (MRIgetVoxVal(mri_mask_src, x, y, z, 0) == 0) 
	  continue ;                    // not in the training mask
	for (l = 0 ; l < nlong ; l++) // see if it was wmsa before the end time
	{
	  label = MRIgetVoxVal(mri_long_seg[l], x, y, z, 0) ;
	  if (IS_HYPO(label))
	    break ;
	}
	label = (nlong - l) ;
	MRIsetVoxVal(mri_mask_dst, x, y, z, 0, label+1) ;
      }

  return(mri_mask_dst) ;
}
#if 0
static MRI *
remove_wmsa_voxels(MRI *mri_mask_src, MRI *mri_mask_dst, MRI **mri_long, int nlong)
{
  int   x, y, z, l, label ;

  mri_mask_dst = MRIcopy(mri_mask_src, mri_mask_dst) ;

  for (x = 0 ; x < mri_mask_src->width ; x++)
    for (y = 0 ; y < mri_mask_src->height ; y++)
      for (z = 0 ; z < mri_mask_src->depth ; z++)
      {
	for (l = 0 ; l < nlong-1 ; l++) // make sure they weren't wmsa before the end
	{
	  label = MRIgetVoxVal(mri_long[l], x, y, z, 0) ;
	  if (IS_HYPO(label))
	    break ;
	}
	if (l < nlong-1)  // it was wmsa before last scan, don't train on it
	  MRIsetVoxVal(mri_mask_dst, x, y, z, 0, 0) ;
      }

  return(mri_mask_dst) ;
}
#endif

static MRI *
make_mask_from_segmentation(const char *sdir, const char *subject, char **long_names, int nlong, 
					int min_dist, int max_dist) 
{
  MRI *mri_mask, *mri_min_dist, *mri_max_dist, *mri_long_seg[MAX_LONG] ;
  int i, l ;
  char fname[STRLEN] ;

  for (l = 0 ; l < nlong ; l++)
  {
    sprintf(fname, "%s/%s%s.%s_base/mri/%s", sdir, subject, long_names[l], subject, seg_name) ;
    mri_long_seg[l] = MRIread(fname) ;
    if (mri_long_seg[l] == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load segmentation from %s",
		Progname, fname) ;
  }
  mri_min_dist = label_voxels_with_wmsa(mri_long_seg, nlong, NOT_WMSA, FUTURE_WMSA, WMSA) ;
  mri_mask = MRIbinarize(mri_min_dist, NULL, 1, 0, WMSA) ;  // wmsa now or in the future
  MRIbinarize(mri_min_dist, mri_min_dist, 1, 0, NOT_WMSA) ; // will dilate into ring of NOT_WMSA

  for (i = 0 ; i < min_dist ; i++)
    MRIdilate(mri_min_dist, mri_min_dist) ;
  mri_max_dist = MRIcopy(mri_min_dist, NULL) ;
  for (i = min_dist ; i < max_dist ; i++)
    MRIdilate(mri_max_dist, mri_max_dist) ;

  MRIsubtract(mri_max_dist, mri_min_dist, mri_max_dist) ;
  MRIadd(mri_mask, mri_max_dist, mri_mask) ;

  label_wmsa_and_future_wmsa(mri_mask, mri_mask, mri_long_seg, nlong, NOT_WMSA, FUTURE_WMSA, WMSA) ;
  for (l = 0 ; l < nlong ; l++)
    MRIfree(&mri_long_seg[l]) ;
  MRIfree(&mri_min_dist) ; MRIfree(&mri_max_dist) ; 
  return(mri_mask) ;
}

#define MIN_WMSA_DIST   1
#define MAX_WMSA_DIST   2

static RANDOM_FOREST *
train_rforest(const char *training_file_name, RF_PARMS *parms, int wsize, int nlong, 
	      char **long_names, int nvols, const char **vol_names, const char *sdir, 
	      const char *seg_name, double  *paccuracy)
{
  RANDOM_FOREST *rf ;
  FILE          *fp ;
  double        **training_data ;
  int           *training_classes, x, y, z, x0, y0, z0, xk, yk, zk ;
  int           nsubjects, n, ntraining, l, v, label, whalf, t, tno ;
  char          line[MAX_LINE_LEN], *cp, fname[MAX_LINE_LEN], *subject ;
  MRI           *mri_intensity[MAX_LONG][MAX_VOLS],*mri_int,*mri_training_mask;
  MATRIX        *mX, *mXpinv ;
  VECTOR        *vP, *vY ;

  mX = MatrixAlloc(nlong-1, 2, MATRIX_REAL) ;
  vY = VectorAlloc(nlong-1, MATRIX_REAL) ;
  vP = VectorAlloc(2, 1) ;  // slope and offset
  mXpinv = NULL ;
  for (l = 0 ; l < nlong-1 ; l++)
    VECTOR_ELT(vY, l+1) = long_times[l] ;

  // use n-1 long intensity volumes to predict wmsa of the nth
  whalf = (wsize-1)/2 ;

  // one feature for each long vol, plus one for the slope, minus the last time point
  if (nlong > 2)
    parms->nfeatures = wsize*wsize*wsize*(nlong)*(nvols) ;
  else
    parms->nfeatures = wsize*wsize*wsize*(nlong-1)*(nvols) ;
  rf = RFalloc(parms->ntrees, parms->nfeatures, nlong+1, parms->max_depth, const_cast<char**>(wmsa_class_names), parms->nsteps) ;

  fp = fopen(training_file_name, "r") ;
  if (fp == NULL)
    ErrorReturn(NULL, (ERROR_NOFILE, "train_rforest(%s): could not open file",
		       training_file_name)) ;

  nsubjects = 0 ;
  cp = fgetl(line, MAX_LINE_LEN, fp) ;
  while (cp != NULL)
  {
    cp = fgetl(line, MAX_LINE_LEN, fp) ;
    nsubjects++ ;
  }
  rewind(fp) ;
  printf("%d subjects read from input file %s\n", nsubjects, training_file_name);

  for (ntraining = n = 0 ; n < nsubjects ; n++)
  {
    subject = fgetl(line, MAX_LINE_LEN, fp) ;
    printf("analyzing subject %s: %d of %d\n", subject, n+1, nsubjects) ;
    l = nlong-1 ;   // only use last one - voxels that will become WMSAs or not
    mri_training_mask = 
      make_mask_from_segmentation(sdir,subject,long_names,nlong,MIN_WMSA_DIST,MAX_WMSA_DIST);
    ntraining += MRItotalVoxelsOn(mri_training_mask, 1) ;
    MRIfree(&mri_training_mask) ; 
  }
  printf("allocating %d %2.1fM training vectors = %2.1lfM total\n", parms->nfeatures, ntraining/(1024.0*1024), (double)parms->nfeatures*ntraining / (1024*1024.0)) ;
  training_classes = (int *)calloc(ntraining, sizeof(int)) ;
  training_data = (double **)calloc(ntraining, sizeof(double *)) ;
  rewind(fp) ;

  for (tno = n = 0 ; n < nsubjects ; n++)
  {
    subject = fgetl(line, MAX_LINE_LEN, fp) ;
    printf("processing subject %s: %d of %d\n", subject, n+1, nsubjects) ;
    mri_training_mask = 
      make_mask_from_segmentation(sdir,subject,long_names,nlong,MIN_WMSA_DIST,MAX_WMSA_DIST);
    for (l = 0 ; l < nlong-1 ; l++)
    {
      for (v = 0 ; v < nvols ; v++)
      {
	sprintf(fname, "%s/%s%s.%s_base/mri/%s", sdir, 
		subject,long_names[l],subject,vol_names[v]);
	mri_intensity[l][v] = MRIread(fname) ;
	if (mri_intensity[l][v] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not intensity volume from %s",
		  Progname, fname) ;
      }
    }
    for (x0 = 0 ; x0 < mri_training_mask->width; x0++)
      for (y0 = 0 ; y0 < mri_training_mask->height; y0++)
	for (z0 = 0 ; z0 < mri_training_mask->depth; z0++)
	{
	  int l ;
	  
	  label = MRIgetVoxVal(mri_training_mask, x0, y0, z0, 0) ;
	  if (label == 0)
	    continue ;

	  if (tno >= ntraining)
	    break ;
	  training_data[tno] = (double *)calloc(parms->nfeatures, sizeof(double)) ;
	  if (training_data[tno] == NULL)
	    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %dth %d-long training vector",
		n, parms->nfeatures*sizeof(double)) ;
	  label-- ;   // classes to RF must be 0-based, but we want nonzero for visualization
	  training_classes[tno] = label ;
	  if (x0 == Gx && y0 == Gy && z0 == Gz)
	    printf("training voxel (%d, %d, %d): maps to tno %d\n", x0, y0, z0, tno) ;

	  for (t = 0, xk = -whalf ; xk <= whalf ; xk++)
	    for (yk = -whalf ; yk <= whalf ; yk++)
	      for (zk = -whalf ; zk <= whalf ; zk++)
		for (v = 0 ; v < nvols ; v++)
		{
		  float val ;
		  {
		    for (l = 0 ; l < nlong-1 ; l++)
		    {
		      mri_int = mri_intensity[l][v] ;
		      x = mri_int->xi[x0+xk] ;
		      y = mri_int->yi[y0+yk] ;
		      z = mri_int->zi[z0+zk] ;
		      val = MRIgetVoxVal(mri_int,x,y,z,0) ;
		      *MATRIX_RELT(mX, l+1, 1) = val ;
		      *MATRIX_RELT(mX, l+1, 2) = 1 ;
		      training_data[tno][t++] = val ;
		    }
		    if (nlong > 2) // otherwise not enough to estimate slope from nlong-1 tps
		    {
		      mXpinv = MatrixPseudoInverse(mX, mXpinv) ;
		      vP = MatrixMultiply(mXpinv, vY, vP) ;
		      training_data[tno][t++] = VECTOR_ELT(vP, 1) ;
		    }
		  }
		}
	  if (x0 == Gx && y0 == Gy && z0 == Gz)
	  {
	    int f ;
	    printf("training voxel (%d, %d, %d): maps to tno %d: \n", x0, y0, z0, tno) ;
	    for (f = 0 ; f < rf->nfeatures ; f++)
	      printf("%2.1f ", training_data[tno][f]) ;
	    printf("\n") ;
	    Gdiag_no = tno ;
	  }
	  tno++ ;
	}
    for (v = 0 ; v < nvols ; v++)
      for (l = 0 ; l < nlong-1 ; l++)
	MRIfree(&mri_intensity[l][v]) ;
    MRIfree(&mri_training_mask) ;
  }
  RFtrain(rf, parms->feature_fraction, parms->training_fraction, training_classes, training_data,ntraining);

  if (parms->training_fraction < 1)
  {
    int correct ;
    correct = RFcomputeOutOfBagCorrect(rf, training_classes, training_data,ntraining);
    printf("out of bag: %d of %d correct, %2.2f%%\n", correct, ntraining,
	   100.0*correct/ntraining) ;
    *paccuracy = (double)correct/ntraining ;
  }

  for (tno = 0 ; tno < ntraining ; tno++)
    free(training_data[tno]) ;
  free(training_data) ; free(training_classes) ;
  fclose(fp) ;
  if (mXpinv)
    MatrixFree(&mXpinv) ; 
  VectorFree(&vY) ; VectorFree(&vP) ; MatrixFree(&mX) ;
  return(rf) ;
}
static int
classify_subjects(RANDOM_FOREST *rf, const char *subject_list_file, int wsize, int nlong, 
		  char **long_names, int nvols, const char **vol_names, const char *sdir, 
		  const char *seg_name, const char *classify_name)
{
  FILE          *fp ;
  int           x, y, z, x0, y0, z0, xk, yk, zk, label ;
  int           nsubjects, l, v, whalf, t, classnum ;
  char          line[MAX_LINE_LEN], fname[MAX_LINE_LEN], *subject ;
  MRI           *mri_intensity[MAX_LONG][MAX_VOLS],*mri_int, *mri_aseg[MAX_LONG], *mri_labeled;
  MATRIX        *mX, *mXpinv ;
  VECTOR        *vP, *vY ;
  double        *feature_vec ;

  printf("classifying subjects...\n") ;
  mX = MatrixAlloc(nlong-1, 2, MATRIX_REAL) ;
  vY = VectorAlloc(nlong-1, MATRIX_REAL) ;
  vP = VectorAlloc(2, 1) ;  // slope and offset
  mXpinv = NULL ;
  for (l = 0 ; l < nlong-1 ; l++)
    VECTOR_ELT(vY, l+1) = long_times[l] ;

  // use n-1 long intensity volumes to predict wmsa of the nth
  whalf = (wsize-1)/2 ;

  fp = fopen(subject_list_file, "r") ;
  if (fp == NULL)
    ErrorReturn(ERROR_NOFILE, (ERROR_NOFILE, "train_rforest(%s): could not open file",
		       subject_list_file)) ;

  mri_labeled = NULL ;
  nsubjects = 0 ;
  feature_vec = (double *)calloc(rf->nfeatures, sizeof(double)) ;
  do
  {
    subject = fgetl(line, MAX_LINE_LEN, fp) ;
    if (subject == NULL)
      break ;
    printf("processing subject %s: %d\n", subject, nsubjects++) ;

    // read in all the volumes for each of the time points
    for (l = 0 ; l < nlong-1 ; l++)
    {
      sprintf(fname, "%s/%s%s.%s_base/mri/%s", sdir, 
	      subject,long_names[l],subject, aseg_name);
      mri_aseg[l] = MRIread(fname) ;
      if (mri_aseg[l] == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not aseg volume from %s",
		  Progname, fname) ;
      if (mri_labeled == NULL)
	mri_labeled = MRIclone(mri_aseg[l], NULL) ;
      for (v = 0 ; v < nvols ; v++)
      {
	sprintf(fname, "%s/%s%s.%s_base/mri/%s", sdir, 
		subject,long_names[l],subject,vol_names[v]);
	mri_intensity[l][v] = MRIread(fname) ;
	if (mri_intensity[l][v] == NULL)
	  ErrorExit(ERROR_NOFILE, "%s: could not intensity volume from %s",
		  Progname, fname) ;
      }
    }
    for (x0 = 0 ; x0 < mri_labeled->width; x0++)
      for (y0 = 0 ; y0 < mri_labeled->height; y0++)
	for (z0 = 0 ; z0 < mri_labeled->depth; z0++)
	{
	  int l, nzero ;
	  double pval ;

	  if (x0 == Gx && y0 == Gy && z0 == Gz)
	    DiagBreak() ;
	  for (label = nzero = l = 0 ; l < nlong-1 ; l++)
	  {
	    label = MRIgetVoxVal(mri_aseg[l], x0, y0, z0, 0) ;
	    if (!IS_WHITE_MATTER(label) && !IS_HYPO(label))
	      break ;
	    if (label == 0)
	      nzero++ ;
	  }
	  if (label == 0)
	    continue ;

#if 0
	  if (l < nlong-1)   // wasn't labeled wm or hypo
	    continue ;
#endif

	  for (t = 0, xk = -whalf ; xk <= whalf ; xk++)
	    for (yk = -whalf ; yk <= whalf ; yk++)
	      for (zk = -whalf ; zk <= whalf ; zk++)
          for (v = 0 ; v < nvols ; v++)
          {
            float val ;
            {
              for (l = 0 ; l < nlong-1 ; l++)
              {
                mri_int = mri_intensity[l][v] ;
                x = mri_int->xi[x0+xk] ;
                y = mri_int->yi[y0+yk] ;
                z = mri_int->zi[z0+zk] ;
                val = MRIgetVoxVal(mri_int,x,y,z,0) ;
                *MATRIX_RELT(mX, l+1, 1) = val ;
                *MATRIX_RELT(mX, l+1, 2) = 1 ;
                feature_vec[t++] = val ;
              }
              if (nlong>2)
              {
                mXpinv = MatrixPseudoInverse(mX, mXpinv) ;
                vP = MatrixMultiply(mXpinv, vY, vP) ;
                feature_vec[t++] = VECTOR_ELT(vP, 1) ;
              }
            }
          }
	  if (x0 == Gx && y0 == Gy && z0 == Gz)
	    DiagBreak() ;
	  classnum = RFclassify(rf, feature_vec, &pval, -1)+1 ;
#if 0
	  if (classnum == 0)
	    classnum = WM_hypointensities ;
	  else
	    classnum = MRIgetVoxVal(mri_aseg[nlong-2], x, y, z, 0) ;
#endif
	  MRIsetVoxVal(mri_labeled, x0, y0, z0, 0, classnum) ;
	}
    for (l = 0 ; l < nlong-1 ; l++)
    {
      MRIfree(&mri_aseg[l]) ;
      for (v = 0 ; v < nvols ; v++)
        MRIfree(&mri_intensity[l][v]) ;
    }
    sprintf(fname, "%s/%s%s.%s_base/mri/%s", sdir, 
	    subject,long_names[0],subject,classify_name);
    printf("writing classification to %s...\n", fname) ;
    MRIwrite(mri_labeled, fname) ;
    MRIfree(&mri_labeled) ;
  } while (subject != NULL);
  free(feature_vec) ;
  fclose(fp) ;
  if (mXpinv)
    MatrixFree(&mXpinv) ; 
  MatrixFree(&mX) ;
  VectorFree(&vY) ; VectorFree(&vP) ;
  return(NO_ERROR) ;
}
