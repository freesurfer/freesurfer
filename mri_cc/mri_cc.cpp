/**
 * @brief program for segmenting the corpus callosum
 *
 * segments the callosum into 5 parts divided along the primary
 * eigendirection (mostly anterior/posterior)
 * and writes the results into aseg_with_cc.mgz
 */
/*
 * Original Authors: Bruce Fischl and Peng Yu
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "mrimorph.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "transform.h"
#include "talairachex.h"
#include "matrix.h"
#include "voxlist.h"
#include "mriTransform.h"
#include "mrisegment.h"
#include "tritri.h"

// these were originally defined in cma.h but they created a big problem as
// they became out of sync with FreeSurferColorLUT.txt
#define MIN_CC_EXTRA 230 
#define MAX_CC_EXTRA 249

#define MAX_CC_ROT      RADIANS(7)
#define CC_ANGLE_DELTA  RADIANS(.25)
#define LABEL_IN_CC     128
#define LABEL_BORDER_CC 64
#define MAX_CENTRAL_SLICES      500

static int norm_thresh = 40 ;
static double max_cc_rot = MAX_CC_ROT ;
static int use_aseg = 1 ;
static int write_lta = 0 ;
static int force = 0 ;
static MRI *find_voxels_close_to_both_hemis(MRI *mri_aseg, int lh_label, int rh_label, int wsize) ;
int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
static void     print_usage();
static const char     *wmvolume = "mri/wm" ;
const char            *Progname ;
static int      dxi=2;  // thickness on either side of midline
static int      x_edge=0, y_edge=0;
static int      write_cc = 0 ;
static double cc_tal_x = 0.0 ;
static double cc_tal_y = 0.0 ;
static double cc_tal_z = 27.0 ;
static int fornix = 0 ;
static int lh_only = 0 ;
static int rh_only = 0 ;
static int skip = 0 ;
static LTA *lta = 0;
static char output_fname[STRLEN] = "aseg_with_cc.mgz";
static char norm_fname[STRLEN] = "norm.mgz" ;
static char aseg_fname[STRLEN] = "aseg.mgz" ;
static char lta_fname[STRLEN] = "" ;
static MRI *remove_fornix_new(MRI *mri_slice, MRI *mri_slice_edited) ;
static int cc_cutting_plane_correct(MRI *mri_aseg,
                                    double x0, double y0, double z0,
                                    double dx, double dy, double dz,
                                    VOXEL_LIST *vl_left_wm,
                                    VOXEL_LIST *vl_right_wm, int debug);
static MRI *find_cc_with_aseg(MRI *mri_aseg, MRI *mri_cc, LTA **plta,
                              double *pxc, double *pyc, double *pzc, int thick,
                              MRI *mri_norm, MRI *mri_fornix) ;
static int find_cc_slice(MRI *mri,
                         double *pccx, double *pccy, double *pccz,
                         const LTA *lta, MRI *mri_tal_cc) ;
static int find_corpus_callosum(MRI *mri,
                                double *ccx, double *ccy, double *ccz,
                                const LTA *lta, MRI *mri_tal_cc) ;
static MRI *remove_fornix(MRI *mri_filled, int xv, int yv, int zv);
static int edge_detection(MRI *mri_temp, int edge_count,int signal);
static int labels[] =
{
  THICKEN_FILL,
  NBHD_FILL,
  VENTRICLE_FILL,
  DIAGONAL_FILL,
  DEGENERATE_FILL
};
#define NLABELS  sizeof(labels) / (sizeof(labels[0]))
#define MAX_SLICES        21  /* 41*/
#define HALF_SLICES       ((MAX_SLICES-1)/2)
#define CUT_WIDTH         1
#define HALF_CUT          ((CUT_WIDTH-1)/2)
#define SEARCH_STEP       3
#define MAX_OFFSET        50
#define CC_VAL            100

/* aspect ratios are dy/dx */
#define MIN_CC_AREA       350  /* smallest I've seen is 389 */
#define MAX_CC_AREA      1400  /* biggest I've seen is 1154 */
#define MIN_CC_ASPECT     0.1
#define MAX_CC_ASPECT     0.75

static char sdir[STRLEN] = "" ;

double findMinSize(MRI *mri)
{
  double xsize, ysize, zsize, minsize;
  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;
  // there are 3! = 6 ways of ordering
  //             xy  yz  zx
  // x > y > z    z min
  // x > z > y    y min
  // z > x > y    y min
  //////////////////////////
  // y > x > z    z min
  // y > z > x    x min
  // z > y > x    x min
  if (xsize > ysize)
  {
    minsize = (ysize > zsize) ? zsize : ysize;
  }
  else
  {
    minsize = (zsize > xsize) ? xsize : zsize;
  }

  return minsize;
}


#define MAX_CC_DIVISIONS     50
#define DEFAULT_CC_DIVISIONS 5
static int cc_divisions = DEFAULT_CC_DIVISIONS ;

int
main(int argc, char *argv[])
{
  char        ifname[STRLEN], ofname[STRLEN],  data_dir[STRLEN], *cp ;
  int         nargs, msec, y, z, xi, temp, i, j, k ;
  double      xc,yc,zc, xv, yv, zv;
  MATRIX      *mrot = NULL, *mtrans = NULL;
  Timer then ;
  MRI         *mri_tal=NULL, *mri_talheader=NULL, *mri_header=NULL, *mri_cc;
  FILE        *fp=NULL;
  LTA         *lta2 = 0;
  MRI         *mri_wm = NULL,
               *mri_cc_tal = NULL, *mri_fornix = NULL, *mri_aseg = NULL ;
  double      means[3] ;
  float       volume[MAX_CC_DIVISIONS], evalues[3],
              ez_x, ez_y, ez_z, zf, zf_low=256, zf_high=0;
  MATRIX      *m_evectors ;

  nargs = handleVersionOption(argc, argv, "mri_cc");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  fflush(stdout);
  fflush(stderr);

  if ((argc < 2) || (argc > 2))
  {
    print_usage();
    ErrorExit(ERROR_BADPARM,
              "", Progname);
  }
  then.reset() ;

  if (!strlen(sdir))
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  strcpy(data_dir, sdir) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    int req = snprintf(ifname,STRLEN,"%s/%s/mri/cc_volume_%d.txt",data_dir,argv[1], dxi) ; 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if ((fp = fopen(ifname, "a")) == NULL)
    {
      ErrorReturn
      (ERROR_BADFILE,
       (ERROR_BADFILE,
        "could not open cc volume measurement file %s",
        ifname));
    }
    printf("writing results to %s\n",ifname);
  }

  if (use_aseg)
  {
    MRI *mri_norm ;

    int req = snprintf(ifname,STRLEN, "%s/%s/mri/%s",data_dir,argv[1], aseg_fname) ;  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    printf("reading aseg from %s\n", ifname);
    mri_aseg = MRIread(ifname) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s",
                Progname, ifname) ;
    if (lh_only || rh_only)
    {
      int req = snprintf(ofname,STRLEN,"%s/%s/mri/%s",data_dir,argv[1], output_fname) ; 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "copying aseg WITHOUT callosum to %s...\n", ofname) ;
      MRIwrite(mri_aseg, ofname) ;
      exit(0) ;
    }
    if (MRIvoxelsInLabel(mri_aseg, CC_Central) > 0)
    {
      if (force == 0)
      {
        fflush(stdout);
        fflush(stderr);
        ErrorExit
        (77,
         "%s: volume %s already has CC in it.\n"
         "Run with -force to reprocess (not recommended)\n",
         Progname, ifname) ;
      }
      // need to replace the cc labels with either lh or rh wm here...
    }

    req = snprintf(ifname,STRLEN,"%s/%s/mri/%s",data_dir,argv[1],norm_fname) ;  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("reading norm from %s\n", ifname);
    mri_norm = MRIread(ifname) ;
    if (mri_norm == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s",
                Progname, ifname) ;
    mri_fornix = MRIclone(mri_aseg, NULL) ;
    mri_cc = find_cc_with_aseg(mri_aseg, NULL, &lta,
                               &xc,&yc,&zc,dxi,mri_norm, mri_fornix) ;
    MRIfree(&mri_norm) ;

    // set cc center in xc, yc, zc
  }
  else
  {
    int req = snprintf(ifname, STRLEN, "%s/%s/%s",data_dir,argv[1],wmvolume) ;   
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("reading white matter volume from %s\n", ifname);
    mri_wm = MRIread(ifname) ;

    req = snprintf(ifname,STRLEN, "%s/%s/mri/transforms/talairach.xfm",data_dir,argv[1]) ;  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    lta = LTAreadEx(ifname);
    if (lta==0)
    {
      ErrorExit(ERROR_BADPARM,"ERROR: cound not load lta from %s.\n", ifname);
    }
    fprintf(stdout,
            "INFO: Using %s and its offset for Talairach volume ...\n",
            ifname);

    for (i = 0 ; i < NLABELS ; i++)
    {
      MRIreplaceValues(mri_wm, mri_wm, labels[i], 0) ;
    }

    mri_talheader =
      MRIallocHeader(mri_wm->width,
                     mri_wm->height,
                     mri_wm->depth,
                     mri_wm->type,1);
    MRIcopyHeader(mri_wm, mri_talheader); // not allocate memory, though

    ModifyTalairachCRAS(mri_talheader, lta);

    mri_tal =
      MRIalloc(mri_wm->width,
               mri_wm->height,
               mri_wm->depth,
               mri_wm->type);
    MRIcopyHeader(mri_talheader, mri_tal);
    // now fill the talairach volume values
    MRItoTalairachEx(mri_wm, mri_tal, lta);

    // binalize the talairach volume (mri_tal)
    MRIbinarize(mri_tal, mri_tal,
                DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      int req = snprintf(ofname,STRLEN,"%s/%s/mri/wm_tal.mgz",data_dir,argv[1]) ; 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout,
              "writing talairach transformed white matter volume to %s...\n",
              ofname) ;
      MRIwrite(mri_tal, ofname) ;
    }

    //find the transform matrix
    mtrans = MatrixAlloc(4, 4, MATRIX_REAL) ;
    mrot = MatrixAlloc(4, 4, MATRIX_REAL) ;

    //try method 2 to get the rotation matrix
    req = snprintf(ifname,STRLEN,"%s/%s/mri/transforms/talairach.xfm",data_dir,argv[1]) ; 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    lta2 = LTAreadEx(ifname);
    mtrans=lta2->xforms[0].m_L;
    Trns_ExtractRotationMatrix (mtrans,mrot);
    *MATRIX_RELT(mrot, 1, 4) = mtrans->rptr[1][4];
    *MATRIX_RELT(mrot, 2, 4) = mtrans->rptr[2][4];
    *MATRIX_RELT(mrot, 3, 4) = mtrans->rptr[3][4];
    lta2->xforms[0].m_L=mrot;

    //rotation wm volume to be upright, using cc volume temporarily
    mri_header = MRIallocHeader(mri_wm->width,
                                mri_wm->height,
                                mri_wm->depth,
                                mri_wm->type,1);
    MRIcopyHeader(mri_wm, mri_header);
    ModifyTalairachCRAS(mri_header, lta2);
    mri_cc = MRIcopy(mri_wm, NULL) ;
    MRIcopyHeader(mri_header, mri_cc);
    MRItoTalairachEx(mri_wm, mri_cc, lta2);
    // binalize the rotated wm  volume (mri_cc)
    MRIbinarize(mri_cc, mri_cc,
                DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      int req = snprintf(ofname,STRLEN,"%s/%s/mri/wm.mgz",data_dir,argv[1]) ; 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      fprintf(stdout,
              "writing rotated white matter volume to %s...\n", ofname) ;
      MRIwrite(mri_cc, ofname) ;
    }

    //now start cc segmentation in talairach space
    mri_cc_tal = MRIcopy(mri_tal, NULL) ;
    MRIcopyHeader(mri_talheader, mri_cc_tal);
    MRIvalueFill(mri_cc_tal, 0) ;

    //most of the work is done in find_corpus_callosum function
    find_corpus_callosum(mri_tal,
                         &cc_tal_x,&cc_tal_y,&cc_tal_z,
                         lta, mri_cc_tal);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      int req = snprintf(ofname,STRLEN,"%s/%s/mri/cc_tal.mgz",data_dir,argv[1]) ;   
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stdout, "writing output to %s...\n", ofname) ;
      MRIwrite(mri_cc_tal, ofname) ;
    }

    //starting the volume measurement

    //transform cc volume from talairach space to normal space
    MRIfromTalairachEx(mri_cc_tal, mri_wm, lta);
    // binalize the rotated cc volume (mri_wm)
    MRIbinarize(mri_wm, mri_wm, CC_VAL/2-1, 0, 100) ;
    req = snprintf(ofname,STRLEN,"%s/%s/mri/cc_org.mgz",data_dir,argv[1]) ;  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout,
            "writing corpus callosum in original space to %s...\n",
            ofname) ;
    MRIwrite(mri_wm, ofname) ;

    //trying to find the position of mid-sagital plane
    MRIcopyHeader(mri_talheader, mri_cc_tal);
    MRItalairachVoxelToVoxelEx(mri_cc_tal,
                               cc_tal_x, cc_tal_y, cc_tal_z,
                               &xv, &yv, &zv,
                               lta) ;

    //rotate the cc volume by the rotation matrix calculated above
    MRIcopyHeader(mri_header, mri_cc);
    MRItoTalairachEx(mri_wm, mri_cc, lta2);
    // binalize the rotated cc volume (mri_cc)
    MRIbinarize(mri_cc, mri_cc, CC_VAL/2-1, 0, 100) ;

    MRIvoxelToTalairachVoxelEx(mri_cc, xv, yv, zv, &xc, &yc, &zc, lta2) ;
  }

  //find the mid-sagital plane there
  xi=nint(xc);
  fprintf(stdout,"cc center is found at %d %d %d\n",xi, nint(yc),nint(zc));

  //count the number if equally segmented five parts
  m_evectors = MatrixAlloc(3, 3, MATRIX_REAL) ;
  MRIprincipleComponents(mri_cc, m_evectors, evalues, means, 1) ;
  printf("eigenvectors:\n") ;
  MatrixPrint(Gstdout, m_evectors) ;
  ez_x = *MATRIX_RELT(m_evectors, 1, 1) ;
  ez_y = *MATRIX_RELT(m_evectors, 2, 1) ;
  ez_z = *MATRIX_RELT(m_evectors, 3, 1) ;
  for (i = 2 ; i <= 3 ; i++)  // find eigenvector that is closest to z axis
  {
    if (fabs(*MATRIX_RELT(m_evectors, 3, i)) > fabs(ez_z))
    {
      ez_x = *MATRIX_RELT(m_evectors, 1, i) ;
      ez_y = *MATRIX_RELT(m_evectors, 2, i) ;
      ez_z = *MATRIX_RELT(m_evectors, 3, i) ;
    }
  }
  if (ez_z  < 0) // orient it anterior/posterior
  {
    ez_x *= -1 ;
    ez_y *= -1 ;
    ez_z *= -1 ;
  }
  MatrixFree(&m_evectors) ;

  //find the bounding box
  for (y = 0 ; y < mri_cc->height ; y++)
  {
    for (z = 0 ; z < mri_cc->depth ; z++)
    {
      if ( nint(MRIgetVoxVal(mri_cc, xi, y, z, 0)) )
      {
        zf = (xi-means[0])*ez_x + (y-means[1])*ez_y + (z-means[2])*ez_z ;
        if (zf < zf_low)
        {
          zf_low = zf ;
        }
        if (zf > zf_high)
        {
          zf_high = zf ;
        }
      }
    }
  }

  for (i = 0 ; i < cc_divisions ; i++)
  {
    volume[i] = 0.0 ;
  }

  for (i = xi-dxi ; i <= xi+dxi ; i++)
  {
    for ( j = 0 ; j < mri_cc->height ; j++)
    {
      for ( k = 0 ; k < mri_cc->depth ; k++)
      {
        if ( nint(MRIgetVoxVal(mri_cc, i, j, k,0))>0)
        {
          zf = (i-means[0])*ez_x + (j-means[1])*ez_y + (k-means[2])*ez_z ;
          if ( zf>=zf_low-10 )
          {
            temp = floor((zf-zf_low)/((zf_high-zf_low+1)/cc_divisions));
            if (temp < 0)
            {
              temp = 0;
            }
            if (temp >= cc_divisions)
            {
              temp = cc_divisions-1;
            }
            volume[temp] +=1 ;
            MRIsetVoxVal(mri_cc, i, j, k, 0, (temp+1)*20+10);
          }
        }
      }
    }
  }
  if (use_aseg)
  {
    int label ;
    MRI *mri_tmp ;

    mri_tmp = MRIlinearTransformInterp(mri_cc, NULL,
                                       lta->inv_xforms[0].m_L,
                                       SAMPLE_NEAREST) ;
    MRIcopy(mri_tmp, mri_cc) ;
    MRIfree(&mri_tmp) ;
    mri_tmp = MRIlinearTransformInterp(mri_fornix, NULL,
                                       lta->inv_xforms[0].m_L,
                                       SAMPLE_NEAREST) ;
    MRIcopy(mri_tmp, mri_fornix) ;
    MRIfree(&mri_tmp) ;

    // paste the cc and fornix labels into the aseg
    for (i = 0 ; i < mri_cc->width ; i++)
    {
      for ( j = 0 ; j < mri_cc->height ; j++)
      {
        for ( k = 0 ; k < mri_cc->depth ; k++)
        {
          if ( nint(MRIgetVoxVal(mri_cc, i, j, k, 0))>0)
          {
            temp = ((nint(MRIgetVoxVal(mri_cc, i, j, k,0))-10)/20) - 1 ;
            if (temp < 0)
            {
              temp = 0 ;
            }
            else if (temp >= cc_divisions)
            {
              temp = cc_divisions-1 ;
            }
	    if (cc_divisions != DEFAULT_CC_DIVISIONS)
	    {
	      label = MIN_CC_EXTRA + temp ;
	      if (label > 255)
		label = 255 ;
	    }
	    else
	    {
	      switch (temp)
	      {
	      default:
	      case 0:
		label = CC_Posterior ;
		break ;
	      case 1:
		label = CC_Mid_Posterior ;
		break ;
	      case 2:
		label = CC_Central ;
		break ;
	      case 3:
		label = CC_Mid_Anterior ;
		break ;
	      case 4:
		label = CC_Anterior ;
		break ;
	      }
	    }
            MRIsetVoxVal(mri_aseg, i, j, k, 0, label) ;
          }
          else if (nint(MRIgetVoxVal(mri_fornix, i, j, k, 0)) > 0 && fornix)
          {
            MRIsetVoxVal(mri_aseg, i, j, k, 0, Fornix) ;
          }
        }
      }
    }

    // fill small holes in cc due to interpolation
    for (i = 0 ; i < mri_cc->width ; i++)
    {
      int left_label, right_label, label ;

      for ( j = 0 ; j < mri_cc->height ; j++)
      {
        for ( k = 0 ; k < mri_cc->depth ; k++)
        {
          label = nint(MRIgetVoxVal(mri_aseg, i, j, k, 0)) ;
          if (label != Left_Cerebral_White_Matter &&
              label != Right_Cerebral_White_Matter)
          {
            continue ;
          }
          left_label = nint(MRIgetVoxVal(mri_aseg, i+1, j, k,0)) ;
          if (IS_CC(left_label) == 0)
          {
            continue ;
          }
          right_label = nint(MRIgetVoxVal(mri_aseg, i-1, j, k, 0)) ;
          if (IS_CC(right_label) == 0)
          {
            continue ;
          }
          MRIsetVoxVal(mri_aseg, i, j, k, 0, left_label) ;
        }
      }
    }

    // now make sure curved cc don't have wrong labels
    {
      MRI_SEGMENTATION *mseg ;
      MRI_SEGMENT      *seg ;
      int              i, j ;

      mseg = MRIsegment(mri_aseg, CC_Mid_Anterior, CC_Mid_Anterior) ;
      if (mseg->nsegments > 1)
      {
        j = MRIsegmentMax(mseg) ;
        for (i = 0 ; i < mseg->nsegments ; i++)
        {
          if (i == j)
          {
            continue ;
          }
          if (mseg->segments[i].y0 > mseg->segments[j].y1)
          {
            seg = &mseg->segments[i] ;
            printf("error in mid anterior detected - correcting...\n") ;
            int k ;
            for (k = 0 ; k < seg->nvoxels ; k++)
              MRIsetVoxVal(mri_aseg, seg->voxels[k].x, seg->voxels[k].y, seg->voxels[k].z, 0, CC_Anterior) ;
          }
        }
      }
      MRIsegmentFree(&mseg) ;
      mseg = MRIsegment(mri_aseg, CC_Mid_Posterior, CC_Mid_Posterior) ;
      if (mseg->nsegments > 1)
      {
        j = MRIsegmentMax(mseg) ;
        for (i = 0 ; i < mseg->nsegments ; i++)
        {
          if (mseg->segments[i].y0 > mseg->segments[j].y1)
          {
            seg = &mseg->segments[i] ;
            printf("error in mid anterior detected - correcting...\n") ;
            int k ;
            for (k = 0 ; k < seg->nvoxels ; k++)
              MRIsetVoxVal(mri_aseg, seg->voxels[k].x, seg->voxels[k].y, seg->voxels[k].z, 0, CC_Posterior) ;
          }
        }
      }
    }

    int req = snprintf(ofname,STRLEN,"%s/%s/mri/%s",data_dir,argv[1], output_fname) ;    
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing aseg with callosum to %s...\n", ofname) ;
    MRIwrite(mri_aseg, ofname) ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    fprintf(fp, "%s %d %d %d %d %d \n",
            argv[1],
            nint(volume[4]),
            nint(volume[3]),
            nint(volume[2]),
            nint(volume[1]),
            nint(volume[0]));
    fprintf(stdout, "%s %d %d %d %d %d \n",
            argv[1],
            nint(volume[4]),
            nint(volume[3]),
            nint(volume[2]),
            nint(volume[1]),
            nint(volume[0]));
  }

  if (write_cc)
  {
    int req = snprintf(ofname,STRLEN,"%s/%s/mri/cc_new.mgz",data_dir,argv[1]) ; 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing corpus callosum output to %s...\n", ofname) ;
    MRIwrite(mri_cc, ofname) ;
  }

  if (mri_fornix && fornix)
  {
    int req = snprintf(ofname,STRLEN,"%s/%s/mri/fornix.mgz",data_dir,argv[1]) ; 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stdout, "writing fornix to to %s...\n", ofname) ;
    MRIwrite(mri_fornix, ofname) ;
  }

  MRIfree(&mri_cc);
  if (mri_tal)
  {
    MRIfree(&mri_tal) ;
  }
  if (mri_wm)
  {
    MRIfree(&mri_wm);
  }
  if (mri_talheader)
  {
    MRIfree(&mri_talheader);
  }
  if (mri_header)
  {
    MRIfree(&mri_header);
  }
  if (mri_cc_tal)
  {
    MRIfree(&mri_cc_tal) ;
  }
  if (mtrans)
  {
    MatrixFree(&mtrans);
  }
  if (mrot)
  {
    MatrixFree(&mrot);
  }
  msec = then.milliseconds() ;
  printf("corpus callosum segmentation took %2.1f minutes\n",(float)msec/(1000.0f*60.0f));
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    fclose(fp);
  }
  printf("#VMPC# mri_cc VmPeak  %d\n",GetVmPeak());
  printf("mri_cc done\n");

  exit(0) ;
  return(0) ;
}

#define CC_SPREAD       10
#define MIN_THICKNESS   3
#define MAX_THICKNESS   20

static int
find_corpus_callosum(MRI *mri_tal,
                     double *pccx, double *pccy, double *pccz,
                     const LTA *lta, MRI *mri_cc_tal)
{
  int         xv, yv, zv, max_y, max_thick=0, thickness=0,
                                 y1, xcc, ycc, x, y,x0, extension=50 ;
  int         flag=0, counts=0;
  double      xr, yr, zr ;
  MRI_REGION  region ;
  int cc_spread, min_thickness, max_thickness, slice_size;
  double voxsize=findMinSize(mri_tal);

  cc_spread = ceil(CC_SPREAD/voxsize);
  min_thickness = ceil(MIN_THICKNESS/voxsize); // estimate bigger
  max_thickness = ceil(MAX_THICKNESS/voxsize);
  slice_size = mri_tal->width;

  MRIboundingBox(mri_tal, 1, &region) ;
  // bounding box center position in voxel coords
  x0 = region.x+region.dx/2 ;

  // this function is called with mri being talairached volume
  // get the talairach coords (0,0,0) in the voxel space
  if (
  mri_tal->linear_transform || 
  lta)
  {
    MRIworldToVoxel(mri_tal, 0.0, 0.0, 0.0,
                    &xr, &yr, &zr);   /* everything is now in tal coords */
    xv = nint(xr) ;
    yv = nint(yr) ;
    zv = nint(zr) ;
  }
  else
  {
    xv = x0;
    yv = region.y+region.dy/2;
    zv = region.z+region.dz/2;
  }

  x0 = xv ;  // BRF!!
  fprintf(stdout, "original seed found at x=%d, y=%d z=%d \n", xv, yv, zv );
  /* find the column with the lowest starting y value of any sign. thick. */
  xcc = ycc = max_y = 0 ;
  for (x = x0-cc_spread ; x <= x0+cc_spread ; x++)
  {
    /* search for first non-zero pixel */
    // in the talairach origin coronal slice from the top
    while (thickness==0 && yv-extension >= 0 && yv+extension <= 256)
    {
      for (y = yv-extension ; y < yv+extension ; y++)
      {
        if (nint(MRIgetVoxVal(mri_tal, x, y, zv, 0)) >= WM_MIN_VAL)
        {
          break ;
        }
      }
      // find y which is greater than WM_MIN_VAL
      /* check to make sure it as reasonably thick */
      if (y < yv+extension ) // within the region and bigger than so far
      {
        for (y1 = y, thickness = 0 ; y1 < slice_size ; y1++, thickness++)
          if (!nint(MRIgetVoxVal(mri_tal, x, y1, zv,0))) // if becomes zero, then break out
          {
            break ;
          }
        if ( thickness > min_thickness && thickness < max_thickness )
        {
          if ( y > max_y || (y == max_y && thickness > max_thick) )
          {
            // found the zero voxel at y1 -> thinckness
            xcc = x ;
            ycc = y+thickness/2 ;  /* in middle of cc */
            max_y = y ;            // mark starting y position
            max_thick = thickness;
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
              fprintf(stdout,
                      "potential cc found at (%d, %d), thickness = %d\n",
                      xcc, ycc, thickness) ;
          }
          else if ( y==max_y &&
                    thickness==max_thick &&
                    (x==xcc+1 ||
                     flag ==1 ))
          {
            flag = 1;
            counts ++;
          }
          else if (flag == 1)
          {
            xcc = xcc+ nint(counts/2);
            flag = 0;
            counts = 0;
          }
        }
      }
      extension += 10;
    }
    thickness = 0;
    extension = 50;
  }

  if (!max_y)
  {
    return(ERROR_BADPARM) ;
  }

  /* now convert the in-plane coords to Talairach coods */
  MRIvoxelToWorld(mri_tal, xcc, ycc, zv, pccx, pccy, pccz) ;
  fprintf(stdout, "%d, %d, %d\n", xcc, ycc, zv);

  find_cc_slice(mri_tal, pccx, pccy, pccz, lta, mri_cc_tal) ;

  return(NO_ERROR) ;
}


static int
find_cc_slice(MRI *mri_tal,
              double *pccx, double *pccy, double *pccz,
              const LTA *lta, MRI *mri_cc_tal)
{
  // here we can handle only up to .5 mm voxel size
  int         area[MAX_SLICES*2], flag[MAX_SLICES*2],
              min_area, min_slice, slice, offset,xv,yv,zv,
              i, total_area=0, left=0, right=0;
  MRI         *mri_slice, *mri_filled ;
  double      x_tal, y_tal, z_tal, x, y, z, xvv, yvv, zvv;
  MRI_REGION  region ;
  char        fname[STRLEN] ;
  int         half_slices, ii, jj;
  double      voxsize = findMinSize(mri_tal);
  int         slice_size = mri_tal->width;
  int         max_slices = ceil(MAX_SLICES/voxsize);
  int         max_cc_area = ceil(MAX_CC_AREA/(voxsize*voxsize));
  int         min_cc_area = floor(MIN_CC_AREA/(voxsize*voxsize));

  half_slices = floor(HALF_SLICES/voxsize);
  if ( half_slices <= 0)
  {
    half_slices = 1;
  }

  x_tal = *pccx ;
  y_tal = *pccy ;
  z_tal = *pccz ;
  offset = 0 ;
  xv = yv = zv = 0 ;
  for (slice = 0 ; slice < max_slices ; slice++)
  {
    offset = slice - half_slices ;

    i=0;
    area[slice]=0;
    while (area[slice]<100 && i<=5)
    {
      x = x_tal + offset ;
      y = y_tal ;
      z = z_tal ;
      MRIworldToVoxel(mri_tal, x, y,  z, &x, &y, &z) ;
      xv = nint(x) ;
      yv = nint(y)-i ;
      zv = nint(z) ;
      mri_slice = MRIextractPlane(mri_tal, NULL, MRI_SAGITTAL, xv);
      mri_filled =  MRIfillFG(mri_slice, NULL,
                              zv, yv,0,WM_MIN_VAL,CC_VAL,&area[slice]);
      MRIboundingBox(mri_filled, 1, &region) ;
      if (i++)
      {
        fprintf(stdout,"moved %d in slice %d \n", i-1, slice);
      }
    }
    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= slice_size -1 ||
        region.y+region.dy >= slice_size-1)
    {
      area[slice] = 0 ;
    }

    if ( !(area[slice]>1100 &&
           (nint(y)-region.y>11)) ||
         region.dy>=3.5*(y-region.y) )
    {
      mri_filled = remove_fornix(mri_filled,xv,yv,zv);

      area[slice] = 0;
      flag[slice] = 0;

      for (ii = 0 ; ii < mri_filled->width ; ii++)
      {
        for (jj = 0 ; jj < mri_filled->height ; jj++)
        {
          if ( nint(MRIgetVoxVal(mri_filled, ii, jj, 0,0))>0 )
          {
            area[slice]++;
          }
        }
      }
    }
    else
    {
      flag[slice] = 1;
      //if (offset>=-5&&offset<0) left++;
      //else if (offset<=5&&offset>0) right++;
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "slice[%d] @ (%d, %d, %d): area = %d\n",
              slice, xv, yv, zv, area[slice]) ;

    if ((Gdiag & DIAG_WRITE) && !(slice % 1) && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "cc_slice%d.mgz", slice);
      MRIwrite(mri_slice, fname) ;
      sprintf(fname, "cc_filled%d.mgz", slice);
      MRIwrite(mri_filled, fname) ;
    }
    MRIfillPlane(mri_filled, mri_cc_tal, MRI_SAGITTAL, xv, CC_VAL);

    MRIfree(&mri_filled) ;
    MRIfree(&mri_slice) ;
  }

#if 0
  min_area = 10000 ;
  min_slice = -1 ;
  for (slice = 1 ; slice < max_slices-1 ; slice++)
  {
    if (area[slice] < min_area &&
        (area[slice] >= min_cc_area && area[slice] <= max_cc_area))
    {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }
#else
  min_area = 10000*5 ;
  min_slice = -1 ;
  for (slice = 6 ; slice <= 14 ; slice++)
  {
    for (i=-2, total_area =0; i <=2; i++)
    {
      total_area += area[slice+i];
    }

    if (total_area < min_area &&
        (total_area >= min_cc_area*5 && total_area <= max_cc_area*5))
    {
      min_area = total_area ;
      min_slice = slice ;
    }
  }
#endif

  /* couldn't find a good slice - don't update estimate */
  if (min_slice < 0)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "%s: could not find valid seed for the cc",
                 Progname));

  offset = floor((min_slice - half_slices)/2) ;
  fprintf(stdout, "find offset as %d using area\n", offset);

  //another way to move the central slice
  for (slice = 0 ; slice < max_slices ; slice++)
  {
    if ( (slice-(half_slices+offset)>=-6) &&
         (slice<half_slices+offset) &&
         flag[slice]==1 )
    {
      left++;
    }
    else if ( (slice-(half_slices+offset)<=6) &&
              (slice>half_slices+offset) &&
              flag[slice]==1 )
    {
      right++;
    }
  }
  offset = offset+left-right;
  fprintf(stdout, "find offset as %d using shifting\n", left-right);
  if (abs(offset)>=5)
  {
    offset = 5*(offset/abs(offset));
  }
  *pccx = x = x_tal+floor(offset) ;
  *pccy = y = y_tal ;
  *pccz = z = z_tal ;

  // just for debugging
  MRIworldToVoxel(mri_tal, x, y,  z, &xvv, &yvv, &zvv) ;
  *pccx = xvv ;
  *pccy = yvv ;
  *pccz = zvv ;

  fprintf(stdout,
          "updating initial cc seed to Tal vol "
          "(%.2f, %.2f, %.2f) TAL (%.2f, %.2f, %.2f)\n",
          xvv, yvv, zvv, x, y, z);

  return(NO_ERROR) ;
}


static MRI *
remove_fornix(MRI *mri_filled, int xv, int yv, int zv)
{
  int    x, y, xi_low=255, xi_high=0, yi_low=255, yi_high=0,
               edge_count = 0, length=0;
  int    temp=0, temp2=0, old_temp =0;
  int    x1=0, y1=0, x2=0, y2=0,
         height=mri_filled->height, width=mri_filled->width, flag =0;
  int    i, section1[150];
  MRI *mri_temp1, *mri_temp2;
  char        fname[STRLEN] ;

  mri_temp1=MRIcopy(mri_filled,NULL);
  mri_temp2=MRIcopy(mri_filled,NULL);
  MRIvalueFill(mri_temp1, 0) ;
  MRIvalueFill(mri_temp2, 0) ;

  // mri_temp1 is dilated version of mri_filled
  for (x = 1 ; x < width-1 ; x++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      if ( nint(MRIgetVoxVal(mri_filled, x-1, y, 0, 0)) ||
           nint(MRIgetVoxVal(mri_filled, x, y-1, 0, 0)) ||
           nint(MRIgetVoxVal(mri_filled, x, y+1, 0, 0)) ||
           nint(MRIgetVoxVal(mri_filled, x+1, y, 0, 0)) ||
           nint(MRIgetVoxVal(mri_filled, x, y, 0, 0)) )
      {
        MRIsetVoxVal(mri_temp1,x,y,0,0,100);
      }
    }
  }

  // mri_temp2 is eroded version of mri_temp1
  for (x = 1 ; x < width-1 ; x++)
  {
    for (y = 1 ; y < height-1 ; y++)
    {
      if ( nint(MRIgetVoxVal(mri_temp1, x-1, y, 0, 0)) &&
           nint(MRIgetVoxVal(mri_temp1, x, y-1, 0, 0)) &&
           nint(MRIgetVoxVal(mri_temp1, x, y+1, 0, 0)) &&
           nint(MRIgetVoxVal(mri_temp1, x+1, y, 0, 0)) &&
           nint(MRIgetVoxVal(mri_temp1, x, y, 0, 0)) )
      {
        MRIsetVoxVal(mri_temp2,x,y,0, 0, 100);
      }
    }
  }

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      if ( nint(MRIgetVoxVal(mri_temp2, x, y, 0,0 )) )
      {
        if (x < xi_low)
        {
          xi_low = x ;
        }
        if (y < yi_low)
        {
          yi_low = y ;
        }
        if (x > xi_high )
        {
          xi_high = x ;
        }
        if (y > yi_high)
        {
          yi_high = y ;
        }
      }
    }
  }

  if (yi_high>yi_low+50)
  {
    yi_high=yi_low+40;
    for (x = 0 ; x < width ; x++)
    {
      for (y = yi_high ; y < height ; y++)
      {
        MRIsetVoxVal(mri_temp2, x, y, 0, 0, 0);
        MRIsetVoxVal(mri_temp1, x, y, 0, 0, 0) ;
      }
    }
  }

  sprintf(fname,
          "/space/neo/2/recon/buckner/001015_vc5442/mri/cc_dilation.mgz");
  //MRIwrite(mri_temp1, fname) ;
  mri_filled=MRIcopy(mri_temp2,NULL);
  sprintf(fname,
          "/space/neo/2/recon/buckner/001015_vc5442/mri/cc_filled.mgz");
  //MRIwrite(mri_temp2, fname) ;

  /*find the first edge of the spike */

  x_edge =0;
  y_edge =0;

  /*the first  edge of the spike */
  flag=0;
  for (x_edge=xi_high-nint((xi_high-xi_low)/4);
       x_edge>=xi_low+nint((xi_high-xi_low)/7);
       x_edge--)
  {
    edge_count=0;
    y_edge=0;
    while ( edge_count<3 && y_edge<height-2 )
    {
      length = edge_detection(mri_temp2,edge_count,0);
      if (length>1)
      {
        edge_count++ ;
      }
      if (edge_count==1&&length>1)
      {
        if (length>=20&&x_edge>xi_low+15&&x_edge<zv+20 )
        {
          temp=old_temp;
          if (x1<x_edge)
          {
            y1=old_temp;
            x1=x_edge;
            for (y=temp; y<256; y++)
            {
              MRIsetVoxVal(mri_filled, x_edge, y, 0, 0, 0);
            }
          }
          flag = 2;
        }
        else
        {
          temp=y_edge;
        }

        if ( length<=1 || y_edge<yv-3 )
        {
          //for (y=y_edge+1; y>0; y--)
          // MRIvox(mri_filled, x_edge, y, 0) =0;
          edge_count-=1;
        }
        else if (length>2&&x_edge>xi_low+15&&x_edge<zv+17)
        {
          for (y=y_edge+1; y<256; y++)
          {
            MRIsetVoxVal(mri_filled, x_edge, y, 0, 0, 0);
          }
        }
      }
      else if (length>1&&edge_count>1 && y_edge<yi_high+1 && y_edge>yv+15 )
      {
        edge_count -=1;
      }
    }
    if (edge_count>=2&&flag==0)
    {
      flag=1;
    }
    else if (edge_count<=1&&flag==1&&x_edge>xi_low+13)
    {
      flag = 0;
      y1=old_temp;
      x1=x_edge;
      //   if (x1<zv+13) break;
    }
    old_temp = temp;
  }
  //fprintf(stdout, "first point found at %d %d \n", x1, y1);
  x_edge =0;
  y_edge =0;

  /*the second edge of the spike */

  flag=0;
  //for (y_edge=yi_high-nint((yi_high-yi_low)/3); y_edge>=yi_low+4; y_edge--)
  for (y_edge=yv+20; y_edge>=yv; y_edge--)
  {
    edge_count=0;
    x_edge=0;
    i=yv+20-y_edge;
    section1[i]=0;
    while (x_edge<width-1)
    {
      length = edge_detection(mri_temp2,edge_count,1);
      if (length >=2)
      {
        edge_count++ ;
      }
      if (edge_count==1)
      {
        temp=x_edge;
        if (!section1[i])
        {
          section1[i]=length;
        }
      }
      if (edge_count==2)
      {
        temp2=x_edge-length;
      }
    }

    if (edge_count>=3&&flag==0)
    {
      flag=1;
    }
    else if (edge_count<=2&&flag==1)
    {
      flag = 0;
      x2=old_temp;
      y2=y_edge;
      if (y2<=yi_low+20)
      {
        break;
      }
    }
    else if ( x2==0 &&
              i>=0 &&
              (section1[i]>=19 ||
               (4*section1[i-1]<3*section1[i]))  )
    {
      x2=old_temp;
      y2=y_edge;
      if (y2<=yi_low+22)
      {
        break;
      }
    }
    if ( (edge_count>=4) && temp2<x1)
    {
      old_temp =  temp2-1;
    }
    else
    {
      old_temp = temp;
    }
  }
  //fprintf(stdout, "second point found at %d %d \n", x2, y2);

  if ( x2>0 && x1>xi_low+8 && x1>x2 && x1<zv+5)
  {
    if (x2<x1)
    {
      temp=x1;
      x1=x2;
      x2=temp;
      temp=y1;
      y1=y2;
      y2=temp;
    }
    for (x=x1; x<=x2; x++)
    {
      for (y=nint(y1+(y2-y1)*(x-x1)/(x2-x1))+1;
           y<height;
           y++)
      {
        MRIsetVoxVal(mri_filled, x, y, 0, 0, 0) ;
      }
    }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "/space/neo/2/recon/buckner/010601_vc6977/mri/cc_cut.mgz");
    MRIwrite(mri_filled, fname) ;
  }
  MRIfree(&mri_temp1);
  MRIfree(&mri_temp2);
  return(mri_filled) ;
}


static int
edge_detection(MRI *mri_temp, int edge_count, int signal)
{
  int length = 0, gap = 0;

  if (signal==1)
  {
    while (gap<=2)
    {
      gap=0;
      while ( x_edge < 256)
      {
        if (nint(MRIgetVoxVal(mri_temp, x_edge, y_edge, 0, 0)))
        {
          break ;
        }
        x_edge++;
      }

      while ( x_edge < 256 )
      {
        if (!nint(MRIgetVoxVal(mri_temp, x_edge, y_edge, 0, 0)))
        {
          break ;
        }
        else
        {
          length++ ;
        }
        x_edge++;
      }

      while ( x_edge < 256)
      {
        if (MRIgetVoxVal(mri_temp, x_edge, y_edge, 0, 0))
        {
          break ;
        }
        else
        {
          gap++;
        }
        x_edge++;
      }

      if (gap<=2&&x_edge<256)
      {
        length += gap;
      }
      else
      {
        x_edge -= gap;
        break;
      }
    }
  }
  else
  {
    while (y_edge < 256)
    {
      if (nint(MRIgetVoxVal(mri_temp, x_edge, y_edge, 0, 0)))
      {
        break ;
      }
      y_edge++;
    }

    while ( y_edge < 256 )
    {
      if (!nint(MRIgetVoxVal(mri_temp, x_edge, y_edge, 0, 0)))
      {
        break ;
      }
      else
      {
        length++ ;
      }
      y_edge++;
    }
  }
  return(length);
}

/*----------------------------------------------------------------------

 Parameters:

 Description:
 ----------------------------------------------------------------------*/
#include "mri_cc.help.xml.h"
static void
print_usage()
{
  outputHelpXml(mri_cc_help_xml,
                mri_cc_help_xml_len);
}
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "SDIR"))
  {
    strcpy(sdir, argv[2]) ;
    printf("using %s as SUBJECTS_DIR...\n", sdir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "-HELP")||!stricmp(option, "-USAGE"))
  {
    print_usage();
    exit(1) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else if (!stricmp(option, "ASEG"))
  {
    strcpy(aseg_fname, argv[2]) ;
    printf("will read input aseg from %s\n", aseg_fname);
    use_aseg = 1 ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "NORM"))
  {
    strcpy(norm_fname, argv[2]) ;
    printf("will read norm from %s\n", norm_fname);
    nargs = 1 ;
  }
  else if (!stricmp(option, "LTA"))
  {
    strcpy(lta_fname, argv[2]) ;
    printf("will write lta as %s\n", lta_fname);
    write_lta = 1 ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "force"))
  {
    force = 1 ;
    printf("processing regardless of existence of cc in input volume\n") ;
  }
  else if (!stricmp(option, "lh"))
  {
    lh_only = 1 ;
    printf("assuming only left hemisphere image\n") ;
  }
  else if (!stricmp(option, "rh"))
  {
    rh_only = 1 ;
    printf("assuming only right hemisphere image\n") ;
  }
  else switch (toupper(*option))
    {
    case 'M':
      max_cc_rot = RADIANS(atof(argv[2])) ;
      nargs = 1 ;
      printf("setting maximum rotation search to %2.1f deg\n",
             DEGREES(max_cc_rot)) ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage();
      exit(1) ;
      break ;
    case 'F':
      fornix = 1 ;
      printf("including fornix in segmentation\n") ;
      break ;
    case 'D':
      cc_divisions = atoi(argv[2]) ;
      if (cc_divisions > MAX_CC_DIVISIONS)
	ErrorExit(ERROR_UNSUPPORTED, "%s: too many CC divisions specified (max = %d)\n", cc_divisions, MAX_CC_DIVISIONS) ;
      nargs = 1 ;
      printf("subdividing the cc into %d compartments\n", cc_divisions) ;
      break ;
    case 'S':
      skip = atoi(argv[2]) ;
      nargs = 1 ;
      printf("skipping %d voxels in rotational alignment\n", skip) ;
      break ;
    case 'A':
      use_aseg = atoi(argv[2]) ;
      printf("%susing aseg as input instead of wm volume\n",
             use_aseg ? "" : "not ") ;
      nargs = 1 ;
      break ;
    case 'O':
      strcpy(output_fname, argv[2]) ;
      printf("writing aseg with cc labels to %s\n", output_fname) ;
      nargs = 1 ;
      break ;
    case 'T':
      dxi = atoi(argv[2]);
      fprintf(stdout,"setting callosum thickness to %d mm\n", 2*dxi+1);
      nargs = 1;
      break;

    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


static MRI *
find_cc_with_aseg(MRI *mri_aseg_orig, MRI *mri_cc, LTA **plta,
                  double *pxc, double *pyc, double *pzc,
                  int thick, MRI *mri_norm, MRI *mri_fornix)
{
  LTA         *lta ;
  VOXEL_LIST  *vl_left, *vl_right ;
  double      dx, dy, dz, yrot, zrot, yrot_best, zrot_best ;
  double      x0, x0_best, means[3] ;
  float       evalues[3] ;
  int         y0, z0, y0_best, z0_best, label,
              i, correct, max_correct, xmin, xmax, scale,
              slice_voxels[MAX_CENTRAL_SLICES], best_slice, changed,
              xleft_min, xright_max ;
  MATRIX      *m, *m_yrot, *m_zrot, *m_trans, *m_trans_inv, *m_tmp, *m_evectors ;
  VECTOR      *v1, *v2 ;
  MRI_SEGMENTATION *mseg ;
  MRI         *mri_slice, *mri_slice_edited, *mri_tmp = NULL, *mri_aseg, *mri_midline ;
  MRI_REGION  box  ;

  mri_aseg = MRIcopy(mri_aseg_orig, NULL) ;  // we'll modify it
  if (DIAG_VERBOSE_ON)
  {
    Gdiag |= (DIAG_VERBOSE | DIAG_WRITE) ;
  }

  *plta = lta = LTAalloc(1, NULL) ;
  lta->xforms[0].type = LINEAR_VOX_TO_VOX ;

  // find leftwards extent of right wm, and rightwards extent of left wm
  xleft_min = mri_aseg->width ;
  xright_max = 0 ;
  for (x0 = 0 ; x0 < mri_aseg->width  ; x0++)
  {
    for (y0 = 0 ; y0 < mri_aseg->height ; y0++)
      for (z0 = 0 ; z0 < mri_aseg->depth ; z0++)
      {
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x0, y0, z0, 0) ;
        if (label == Left_Cerebral_White_Matter ||
            label == Left_Cerebral_Cortex )
        {
          if (x0 < xleft_min)
          {
            if (x0 <= 116)
            {
              DiagBreak() ;
            }
            xleft_min = x0 ;
          }
        }
        else if (label == Right_Cerebral_White_Matter ||
                 label == Right_Cerebral_Cortex)
        {
          if (x0 > xright_max)
          {
            if (x0 >= 121)
            {
              DiagBreak() ;
            }
            xright_max = x0 ;
          }
        }
      }
  }
  xmin = MIN(xleft_min, xright_max)-1 ;
  xmax = MAX(xleft_min, xright_max)+1 ;
  box.y = box.z = 0 ;
  box.dy = mri_aseg->height ;
  box.dz = mri_aseg->depth ;
  box.x = xmin-1 ;
  box.dx = (xmax-xmin)+2 ;

  vl_left = VLSTcreateInRegion(mri_aseg, Left_Cerebral_White_Matter,
                               Left_Cerebral_Cortex, NULL, skip, 0, &box) ;
  vl_right = VLSTcreateInRegion(mri_aseg, Right_Cerebral_White_Matter,
                                Right_Cerebral_Cortex, NULL, skip, 0, &box) ;

  printf("%d voxels in left wm, %d in right wm, xrange [%d, %d]\n",
         vl_left->nvox, vl_right->nvox, xmin, xmax) ;

  // find principal components of regions close to both lh and rh wm.
  mri_midline = find_voxels_close_to_both_hemis(mri_aseg,
                Left_Cerebral_Cortex, Right_Cerebral_Cortex, 5) ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_midline, "m.mgz") ;
  }
  m_evectors = MatrixAlloc(3, 3, MATRIX_REAL) ;
  MRIprincipleComponents(mri_midline, m_evectors, evalues, means, 0) ;
  MRIfree(&mri_midline) ;

  dx = 1 ;
  dy = 0 ;
  dz = 0 ;
  yrot_best = 0 ;
  zrot_best = 0 ;
  x0 = x0_best = ((xmin+xmax)/2) ;
  y0 = y0_best = 128 ;
  z0 = z0_best = 128 ;
  max_correct = -1 ;

  dx = *MATRIX_RELT(m_evectors, 1, 3) ;
  dy = *MATRIX_RELT(m_evectors, 2, 3) ;
  dz = *MATRIX_RELT(m_evectors, 3, 3) ;
  if (dx < 0)
  {
    dx *= -1 ;
    dy *= -1 ;
    dz *= -1 ;
  }
  x0 = means[0] ;
  y0 = means[1] ;
  z0 = means[2] ;
  correct = cc_cutting_plane_correct(mri_aseg, x0, y0, z0, dx, dy, dz,
                                     vl_left, vl_right, 0) ;
  max_correct = correct ;
  x0_best = x0 ;
  y0_best = y0 ;
  z0_best = z0 ;
  zrot_best = -asin(dy) ;
  yrot_best = asin(dz/cos(zrot_best)) ;

  m = MatrixAlloc(4, 4, MATRIX_REAL) ;
  m_yrot = MatrixIdentity(4, NULL) ;
  m_zrot = MatrixIdentity(4, NULL) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ;
  v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ;
  VECTOR_ELT(v2, 4) = 1.0 ;
  for (scale = 1 ; scale >= 1 ; scale /= 2)
  {
    printf("searching rotation angles z=[%2.0f %2.0f], y=[%2.0f %2.0f]\n",
           DEGREES(zrot_best-scale*max_cc_rot),
           DEGREES(zrot_best+scale*max_cc_rot),
           DEGREES(yrot_best-scale*max_cc_rot),
           DEGREES(yrot_best+scale*max_cc_rot));
    for (zrot = zrot_best-scale*max_cc_rot ;
         zrot <= zrot_best+scale*max_cc_rot ;
         zrot += scale*CC_ANGLE_DELTA)
    {
      printf("\rsearching scale %d Z rot %2.1f  ", scale, DEGREES(zrot)) ;
      fflush(stdout) ;
      MatrixReallocRotation(4, zrot, Z_ROTATION, m_zrot) ;
      for (yrot = yrot_best-scale*max_cc_rot ;
           yrot <= yrot_best+scale*max_cc_rot ;
           yrot += scale*CC_ANGLE_DELTA)
      {
        // rotate vector away from 1,0,0
        MatrixReallocRotation(4, yrot, Y_ROTATION, m_yrot) ;
        MatrixMultiply(m_yrot, m_zrot, m) ;
        V3_X(v1) = 1 ;
        V3_Y(v1) = 0 ;
        V3_Z(v1) = 0 ;
        MatrixMultiply(m, v1, v2) ;
        dx = V3_X(v2) ;
        dy = V3_Y(v2) ;
        dz = V3_Z(v2) ;

        for (x0 = xmin ; x0 <= xmax ; x0 += 1.0)
        {
          correct = cc_cutting_plane_correct(mri_aseg, x0, y0, z0, dx, dy, dz,
                                             vl_left, vl_right, 0) ;
          if (correct > max_correct)
          {
            x0_best = x0 ;
            max_correct = correct ;
            yrot_best = yrot ;
            zrot_best = zrot ;
          }
        }
      }
    }

    MatrixReallocRotation(4, yrot_best, Y_ROTATION, m_yrot) ;
    MatrixReallocRotation(4, zrot_best, Z_ROTATION, m_zrot) ;
    MatrixMultiply(m_yrot, m_zrot, m) ;

    V3_X(v1) = 1 ;
    V3_Y(v1) = 0 ;
    V3_Z(v1) = 0 ;
    MatrixMultiply(m, v1, v2) ;
    dx = V3_X(v2) ;
    dy = V3_Y(v2) ;
    dz = V3_Z(v2) ;
    max_correct = -1 ;
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    correct = cc_cutting_plane_correct(mri_aseg, x0_best, y0_best, z0_best,
                                       dx, dy, dz, vl_left,vl_right,2);
  printf("global minimum found at slice %2.1f, rotations (%2.2f, %2.2f)\n",
         x0_best, DEGREES(yrot_best), DEGREES(zrot_best)) ;

  m_trans = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_trans, 1, 4) = x0_best ;
  *MATRIX_RELT(m_trans, 2, 4) = y0_best ;
  *MATRIX_RELT(m_trans, 3, 4) = z0_best ;
  m_trans_inv = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_trans_inv, 1, 4) = -128 ;
  *MATRIX_RELT(m_trans_inv, 2, 4) = -128 ;
  *MATRIX_RELT(m_trans_inv, 3, 4) = -128 ;

  MatrixMultiply(m_zrot, m_trans_inv, m) ;
  m_tmp = MatrixMultiply(m_yrot, m, NULL) ;
  MatrixMultiply(m_trans, m_tmp, m) ;
  MatrixInverse(m, lta->xforms[0].m_L) ;
  printf("final transformation (x=%2.1f, yr=%2.3f, zr=%2.3f):\n",x0_best,
         DEGREES(yrot_best), DEGREES(zrot_best)) ;
  MatrixPrint(Gstdout, lta->xforms[0].m_L) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    LTAwrite(lta, "test.lta") ;
  }
  MatrixFree(&m_trans) ;
  MatrixFree(&m_trans_inv) ;
  MatrixFree(&m_tmp) ;

  // remove aseg voxels with norm < norm_thresh (won't really be cc)
  if (mri_norm)
  {
    for (x0 = 0 ; x0 < mri_norm->width ; x0++)
      for (y0 = 0 ; y0 < mri_norm->height ; y0++)
        for (z0 = 0 ; z0 < mri_norm->depth ; z0++)
        {
          if (x0 == Gx && y0 == Gy && z0 == Gz)
          {
            DiagBreak() ;
          }
          if (MRIgetVoxVal(mri_norm, x0, y0, z0, 0) < norm_thresh)
          {
            MRIsetVoxVal(mri_aseg, x0, y0, z0, 0, 0) ;
          }
        }
  }

  mri_tmp = MRIlinearTransformInterp(mri_aseg, NULL,
                                     lta->xforms[0].m_L, SAMPLE_NEAREST) ;
  MRIcopy(mri_tmp, mri_aseg) ;
  MRIfree(&mri_tmp) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_aseg, "asegx.mgz") ;
  }

  if (mri_cc == NULL)
  {
    mri_cc = MRIclone(mri_aseg, NULL) ;
  }

  // update xmin and xmax to be in transformed coords
  xleft_min = mri_aseg->width ;
  xright_max = 0 ;
  for (x0 = 0 ; x0 < mri_aseg->width  ; x0++)
  {
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x0, y0, z0, 0) ;
        if (label == Left_Cerebral_White_Matter)
        {
          if (x0 < xleft_min)
          {
            xleft_min = x0 ;
          }
        }
        else if (label == Right_Cerebral_White_Matter)
        {
          if (x0 > xright_max)
          {
            xright_max = x0 ;
          }
        }
      }
  }
  if (xleft_min > xright_max)
    ErrorExit(ERROR_UNSUPPORTED, "%s: no WM voxels found with norm > %d -- check skull stripping\n",
	      Progname, norm_thresh) ;
  xmin = MIN(xleft_min, xright_max)-1 ;
  xmax = MAX(xleft_min, xright_max)+1 ;
  printf("updating x range to be [%d, %d] in xformed coordinates\n",xmin,xmax);


  for (x0 = xmin ; x0 <= xmax  ; x0++)
  {
    i = x0-xmin ; // slice_voxels index
    slice_voxels[i] = 0 ;
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x0, y0, z0, 0) ;
        if (label == Left_Cerebral_White_Matter)
        {
          label = Right_Cerebral_White_Matter ;
        }
        else if (label == Right_Cerebral_White_Matter)
        {
          label = Left_Cerebral_White_Matter ;
        }
        else
        {
          continue ;  // not wm
        }
        if (MRIlabelsInNbhd(mri_aseg, x0, y0, z0, 1, label) > 0)
        {
          MRIsetVoxVal(mri_cc, x0, y0, z0, 0, LABEL_IN_CC) ;
          slice_voxels[i]++ ;
        }
      }
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_cc, "s1.mgz") ;
  }

  best_slice = 0 ;
  for (x0 = xmin+1 ; x0 <= xmax ; x0++)
  {
    i = x0-xmin ;
    if (slice_voxels[i] > slice_voxels[best_slice])
    {
      best_slice = i ;
    }
  }
  best_slice += xmin ;

  printf("best xformed slice %d\n", best_slice) ;


  // erase voxels in other slices (for now)
  for (x0 = xmin ; x0 <= xmax  ; x0++)
  {
    if (x0 == best_slice)
    {
      continue ;
    }
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        MRIsetVoxVal(mri_cc, x0, y0, z0, 0, 0) ;
      }
  }

  // expand central slice to include all connected wm
  do
  {
    changed = 0 ;
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (best_slice == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        if (MRIgetVoxVal(mri_cc, best_slice, y0, z0, 0) > 0)  // already on
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_aseg, best_slice, y0, z0, 0) ;
        if (label != Left_Cerebral_White_Matter &&
            label != Right_Cerebral_White_Matter)
        {
          continue ;
        }
        if (MRIlabelsInPlanarNbhd(mri_cc, best_slice,
                                  y0, z0, 1, LABEL_IN_CC, MRI_SAGITTAL) > 1)
        {
          MRIsetVoxVal(mri_cc, best_slice, y0, z0, 0, LABEL_IN_CC) ;
          changed = 1 ;
        }
      }
  }
  while (changed > 0) ;

  // expand the cc one slice at a time
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_cc, "s2.mgz") ;
  }
  mri_tmp = MRIcopy(mri_cc, NULL) ;
  for (i = 1 ; i <= thick ; i++)
  {
    if (i == Gdiag_no)
    {
      DiagBreak() ;
    }
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        x0 = best_slice+i ;
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x0, y0, z0, 0) ;
        if (label == Left_Cerebral_White_Matter ||
            label == Right_Cerebral_White_Matter)
        {
          if (MRIlabelsInNbhd(mri_cc, x0, y0, z0, 1, LABEL_IN_CC) > 0)
          {
            MRIsetVoxVal(mri_tmp, x0, y0, z0, 0, LABEL_IN_CC) ;
          }
          else
          {
            MRIsetVoxVal(mri_tmp, x0, y0, z0, 0, 0) ;
          }
        }

        x0 = best_slice-i ;
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_aseg, x0, y0, z0, 0) ;
        if (label == Left_Cerebral_White_Matter ||
            label == Right_Cerebral_White_Matter)
        {
          if (MRIlabelsInNbhd(mri_cc, x0, y0, z0, 1, LABEL_IN_CC) > 0)
          {
            MRIsetVoxVal(mri_tmp, x0, y0, z0, 0, LABEL_IN_CC) ;
          }
          else
          {
            MRIsetVoxVal(mri_tmp, x0, y0, z0, 0, 0) ;
          }
        }
      }
    MRIcopy(mri_tmp, mri_cc) ;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_cc, "s3.mgz") ;
  }

  MRIdilate(mri_cc, mri_tmp) ;
  MRIdilate(mri_tmp, mri_tmp) ;
  MRIdilate(mri_tmp, mri_tmp) ;
  MRIdilate(mri_tmp, mri_tmp) ;
  mseg = MRIsegment(mri_tmp, 1, 255) ;
  i = MRIsegmentMax(mseg) ;
  MRIfree(&mri_tmp) ;
  mri_tmp = MRIsegmentToImage(mri_cc, NULL, mseg, i) ;
  MRIcopy(mri_tmp, mri_cc) ;
  MRIfree(&mri_tmp) ;
  MRIsegmentFree(&mseg) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_cc, "s4.mgz") ;
  }

  VLSTfree(&vl_left) ;
  VLSTfree(&vl_right) ;
  MatrixFree(&m) ;
  MatrixFree(&m_yrot) ;
  MatrixFree(&m_zrot) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  *plta = lta ;
  *pxc = (double)best_slice ;
  *pyc = (double)y0_best ;
  *pzc = (double)z0_best ;

  mri_tmp = MRIclone(mri_cc, NULL) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_cc, "cc_with_fornix.mgz") ;
  }

  MRIcopy(mri_cc, mri_fornix) ;
  for (x0 = best_slice-thick ; x0 <= best_slice+thick ; x0++)  // build fornix model slice-by-slice
  {
    if (x0 == Gdiag_no)
    {
      DiagBreak() ;
    }
    mri_slice = MRIextractPlane(mri_cc, NULL, MRI_SAGITTAL, x0);
    mri_slice_edited = remove_fornix_new(mri_slice, NULL) ;
    MRIfillPlane(mri_slice_edited, mri_tmp, MRI_SAGITTAL, x0, CC_VAL);
    MRIfree(&mri_slice) ;
    MRIfree(&mri_slice_edited) ;
  }
  MRIsubtract(mri_fornix, mri_tmp, mri_fornix) ;
  MRIbinarize(mri_fornix, mri_fornix, LABEL_IN_CC-1, 0, Fornix) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_fornix, "fornix.mgz") ;
  }
  MRIclose(mri_tmp, mri_cc) ;
  MRIfree(&mri_tmp) ;

  // add back in wm voxels that were removed because they were dark
  for (x0 = best_slice-thick ; x0 <= best_slice+thick ; x0++)
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (x0 == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        if (MRIgetVoxVal(mri_fornix, x0, y0, z0, 0) > 0)
        {
          continue ;
        }
        label = MRIgetVoxVal(mri_aseg_orig, x0, y0, z0, 0) ;
        if ((label == Left_Cerebral_White_Matter ||
             label == Right_Cerebral_White_Matter) &&
            (MRIgetVoxVal(mri_cc, x0, y0, z0, 0) == 0) &&
            (MRIlabelsInPlanarNbhd(mri_cc, x0, y0, z0, 1,
                                   LABEL_IN_CC, MRI_SAGITTAL) > 1))
        {
          MRIsetVoxVal(mri_cc, x0, y0, z0, 0, LABEL_IN_CC) ;
        }

      }

  changed = 0 ;
  i = 3*3*3 ;  // max nbhd size
  do
  {
    for (y0 = 0 ; y0 < mri_cc->height ; y0++)
      for (z0 = 0 ; z0 < mri_cc->depth ; z0++)
      {
        if (best_slice == Gx && y0 == Gy && z0 == Gz)
        {
          DiagBreak() ;
        }
        label = MRIgetVoxVal(mri_cc, best_slice, y0, z0, 0) ;
        if (label && MRIlabelsInNbhd(mri_cc, best_slice,
                                     y0, z0, 1, LABEL_IN_CC) >= i)
        {
          changed = 1 ;
          *pyc = y0 ;
          y0_best = y0;
          *pzc = z0 ;
          z0_best = z0;
          break ;
        }
      }
    i-- ;
  }
  while (changed == 0 && i > 0) ;

  if (write_lta)
  {
    LTA* lta2 = LTAalloc(1,mri_norm);
    lta2->xforms[0].m_L  = MatrixCopy(lta->xforms[0].m_L,lta2->xforms[0].m_L);
    lta2->xforms[0].type = LINEAR_VOX_TO_VOX ;
    getVolGeom(mri_norm, &lta2->xforms[0].src);
    getVolGeom(mri_norm, &lta2->xforms[0].dst);

    // adjust translation info so that best slice it at 128 y z
    *MATRIX_RELT(lta2->xforms[0].m_L, 1, 4) += 128 - best_slice;
    // *MATRIX_RELT(lta2->xforms[0].m_L, 2, 4) += 128 - y0_best ;
    // *MATRIX_RELT(lta2->xforms[0].m_L, 3, 4) += 128 - z0_best ;

    LTAwrite(lta2, lta_fname);
    LTAfree(&lta2);
  }


  LTAfillInverse(lta) ;

  return(mri_cc) ;
}

static int
cc_cutting_plane_correct(MRI *mri_aseg, double x0, double y0, double z0,
                         double dx, double dy, double dz,
                         VOXEL_LIST *vl_left, VOXEL_LIST *vl_right,
                         int debug)
{
  int     i, xv_left, yv_left, zv_left, xv_right, yv_right, zv_right,
          width, height, label_left, label_right, correct ;
  double  xf, yf, zf, e1[3], e2[3], n[3], len, x, y ;
  double  dot, dxi, dyi, dzi ;
  MRI     *mri_debug = NULL ;

  if (debug)
  {
    mri_debug = MRIclone(mri_aseg, NULL) ;
  }

  n[0] = dx ;
  n[1] = dy ;
  n[2] = dz ;
  e2[0] = dy ;
  e2[1] = dx ;
  e2[2] = dx ;
  CROSS(e1, n, e2) ;   // e1 is first basis vector
  len = VLEN(e1) ;
  SCALAR_MUL(e1, 1/len, e1) ;
  CROSS(e2, n, e1) ;   // e2 is second basis vector
  len = VLEN(e2) ;
  SCALAR_MUL(e2, 1/len, e2) ;

  width = mri_aseg->depth ;
  height = mri_aseg->height ;  // of extracted plane
  correct = 0 ;
  if (debug)
  {
    for (x = -width/2 ; x <= width/2 ; x += 0.5)
    {
      for (y = -height/2 ; y <= height/2 ; y += 0.5)
      {
        xf = x0 + e1[0]*x + e2[0]*y ;
        yf = y0 + e1[1]*x + e2[1]*y ;
        zf = z0 + e1[2]*x + e2[2]*y ;

#define SDIST 0.25
        if (ISINT(xf))
        {
          xv_left = (int)floor(xf + SDIST*n[0]) ;
          yv_left = (int)floor(yf + SDIST*n[1]) ;
          zv_left = (int)floor(zf + SDIST*n[2]) ;
          xv_right = (int)floor(xf - SDIST*n[0]) ;
          yv_right = (int)floor(yf - SDIST*n[1]) ;
          zv_right = (int)floor(zf - SDIST*n[2]) ;
        }
        else
        {
          xv_left = (int)floor(xf + SDIST*n[0]) ;
          yv_left = (int)floor(yf + SDIST*n[1]) ;
          zv_left = (int)floor(zf + SDIST*n[2]) ;
          xv_right = (int)ceil(xf - SDIST*n[0]) ;
          yv_right = (int)ceil(yf - SDIST*n[1]) ;
          zv_right = (int)ceil(zf - SDIST*n[2]) ;
        }
        if (xv_left == xv_right)
        {
          DiagBreak() ;
        }
        if (xv_left < 0 || xv_left >= mri_aseg->width ||
            yv_left < 0 || yv_left >= mri_aseg->height ||
            zv_left < 0 || zv_left >= mri_aseg->depth)
        {
          continue ;
        }
        label_left = MRIgetVoxVal(mri_aseg, xv_left, yv_left, zv_left, 0) ;

        if (xv_right < 0 || xv_right >= mri_aseg->width ||
            yv_right < 0 || yv_right >= mri_aseg->height ||
            zv_right < 0 || zv_right >= mri_aseg->depth)
        {
          continue ;
        }
        if (mri_debug)
        {
          MRIsetVoxVal(mri_debug, xv_left, yv_left, zv_left, 0, 160) ;
          MRIsetVoxVal(mri_debug, xv_right, yv_right, zv_right, 0, 100) ;
        }
        label_right = MRIgetVoxVal(mri_aseg, xv_right, yv_right, zv_right, 0) ;
        if (label_left == Left_Cerebral_White_Matter &&
            label_right == Right_Cerebral_White_Matter)
        {
          correct += 10 ;
        }
      }
    }
  }

  if (debug)
  {
    char fname[STRLEN] ;
    sprintf(fname, "plane%d.mgz", debug) ;
    MRIwrite(mri_debug, fname) ;
    MRIfree(&mri_debug) ;
  }
  for (i = 0 ; i < vl_left->nvox ; i++)
  {
    if ((vl_left->xi[i] == Gx) &&
        (vl_left->yi[i] == Gy) &&
        (vl_left->zi[i] == Gz))
    {
      DiagBreak() ;
    }
    // left should be on positive side of plane
    dxi = vl_left->xi[i] - x0 ;
    dyi = vl_left->yi[i] - y0 ;
    dzi = vl_left->zi[i] - z0 ;
    dot = dxi*dx + dyi*dy + dzi*dz ;
    if (dot > 0.5)
    {
      correct += 1 ;
    }
  }
  for (i = 0 ; i < vl_right->nvox ; i++)
  {
    if ((vl_right->xi[i] == Gx) &&
        (vl_right->yi[i] == Gy) &&
        (vl_right->zi[i] == Gz))
    {
      DiagBreak() ;
    }
    // right should be on negative side of plane
    dxi = vl_right->xi[i] - x0 ;
    dyi = vl_right->yi[i] - y0 ;
    dzi = vl_right->zi[i] - z0 ;
    dot = dxi*dx + dyi*dy + dzi*dz ;
    if (dot < -0.5)
    {
      correct += 1 ;
    }
  }
  return(correct) ;
}

#define FORNIX_VAL   32
#define LABEL_ERASE  96

static MRI *
remove_fornix_new(MRI *mri_slice, MRI *mri_slice_edited)
{
  int xmin, xmax, x, y, edges_found, val, last_val, i, x1, changed, found,ymax;
  MRI_SEGMENTATION *mseg ;
  MRI              *mri_tmp ;

  if (mri_slice_edited == NULL)
  {
    mri_slice_edited = MRIclone(mri_slice, NULL) ;
  }

  // This was just xmin = mri_slice->width, but can cause an error
  xmin = mri_slice->width - 1; 

  // find posterior/anterior extent of the cc
  ymax = xmax =  0 ;
  for (x = 0 ; x < mri_slice->width; x++)
  {
    for (y = 0 ; y < mri_slice->height; y++)
    {
      if (MRIgetVoxVal(mri_slice, x, y, 0, 0) > 0)
      {
        if (x < xmin)
        {
          xmin = x ;
        }
        if (x > xmax)
        {
          xmax = x ;
        }
        if (y > ymax)
        {
          ymax = y ;
        }
      }
    }
  }

  MRIcopy(mri_slice, mri_slice_edited) ;

  // start far enough in to avoid the curve at the head of the callosum
  //  for (x = xmin + (xmax-xmin)/4 ; x <= xmin + 2*(xmax-xmin)/3 ; x++)
  for (x = xmin ; x <= xmax ; x++)
  {
    // scan downwards from top.
    // First set of connected pixels will be callosum, 2nd will be fornix
    edges_found = 0 ;
    last_val = 0 ;
    for (y = 0 ; y < mri_slice->height ; y++)
    {
      val = (int)MRIgetVoxVal(mri_slice, x, y, 0, 0) ;
      if (val > 0 && last_val == 0)  // transition from off to on
      {
        edges_found++ ;
      }
      last_val = val ;
      if (val && edges_found > 1)  // vertical transition back to on
      {
        MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, FORNIX_VAL);  // mark fornix
      }
    }
  }

  // now use a 45deg line and look for on-off-on transitions


  // set seeds at either end of the cc and grow them inwards
  mri_tmp = MRIclone(mri_slice, NULL) ;
  for (x = xmax-1 ; x <= xmax ; x++)
    for (y = 0 ; y < mri_slice->height ; y++)
    {
      val = (int)MRIgetVoxVal(mri_slice, x, y, 0, 0) ;
      if (val > 0)  // change fornix vals to cc
      {
        MRIsetVoxVal(mri_tmp, x, y, 0, 0, LABEL_IN_CC) ;
      }
    }

  // now expand it posteriorly for 1/3 the length of
  // the cc to recover anterior end
  for (x = xmax-2 ; x >= xmax-2*(xmax-xmin)/3; x--)
  {
    do
    {
      changed = 0 ;
      for (y = 0 ; y < mri_slice->height ; y++)
      {
        val = (int)MRIgetVoxVal(mri_slice_edited, x, y, 0, 0) ;
        if (MRIlabelsInNbhd(mri_tmp, x, y, 0, 1, LABEL_IN_CC) > 0)
          if (val > 0)  // change fornix vals to cc
          {
            val = (int)MRIgetVoxVal(mri_tmp, x, y, 0, 0) ;
            MRIsetVoxVal(mri_tmp, x, y, 0, 0, LABEL_IN_CC) ;
            if (val != LABEL_IN_CC)
            {
              changed++ ;
            }
          }
      }
    }
    while (changed > 0) ;
  }

  // copy rest of definite cc into mri_tmp
  MRIcopyLabel(mri_slice_edited,mri_tmp,LABEL_IN_CC);
  // copy new non-ambiguous labels in
  MRIcopyLabel(mri_tmp, mri_slice_edited, LABEL_IN_CC) ;

  // search anterior posterior around fornix labels.
  // If no non-fornix within a few mm either way
  // mark it as to be erased
  for (x = xmin ; x <= xmax ; x++)
  {
    for (y = 0 ; y < mri_slice->height ; y++)
    {
      if (x == Gx && y == Gy)
      {
        DiagBreak() ;
      }
      val = (int)MRIgetVoxVal(mri_slice_edited, x, y, 0, 0) ;
      if (val != FORNIX_VAL)
      {
        continue ;
      }
      found = 0 ;
      for (i = -3 ; i <= 3 ; i++)
      {
        x1 = x + i ;
        if (x < 0 || x >= mri_slice_edited->width)
        {
          continue ;
        }
        val = (int)MRIgetVoxVal(mri_tmp, x1, y, 0, 0) ;
        if (val == LABEL_IN_CC)
        {
          //          MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, LABEL_IN_CC) ;
          found = 1 ;
          break ;
        }
      }
      if (!found)
      {
        MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, LABEL_ERASE) ;
      }
    }
  }

  /* now check all voxels that were labeled fornix and not yet erased.
     If they border something
     that was erased (i.e. definitely fornix) then they are fornix as well */
  for (x = xmin ; x <= xmax ; x++)
  {
    for (y = 0 ; y < mri_slice->height ; y++)
    {
      if (x == Gx && y == Gy)
      {
        DiagBreak() ;
      }
      val = (int)MRIgetVoxVal(mri_slice_edited, x, y, 0, 0) ;
      if (val != FORNIX_VAL)
      {
        continue ;
      }
      found = 0 ;
      if (MRIlabelsInNbhd(mri_slice_edited, x, y, 0, 1, LABEL_ERASE) > 0)
      {
        MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, LABEL_ERASE) ;
      }
      else
      {
        MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, LABEL_IN_CC) ;
      }
    }
  }
  /* compute the mean sup/inf thickness of the cc in the slices
     where fornix was removed, and use it to shave off remaining
     fornix posterior to what was found.
  */
  {
    double mean_thickness = 0.0, std_thickness = 0.0 ;
    int    num = 0, thickness = 256, first_on, last_on, min_x_fornix = mri_slice->width-1,
           max_thickness, first_off, ystart ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_slice, "s.mgz") ;
      MRIwrite(mri_slice_edited, "ed.mgz") ;
    }

    // add back stuff that is directly anterior to real CC
    for (x = xmin ; x <= xmax ; x++)
      for (y = ymax-10 ; y <= ymax ; y++)
      {
        if (nint(MRIgetVoxVal(mri_slice_edited, x, y, 0, 0)) != LABEL_ERASE)
        {
          continue ;
        }
        if ((nint(MRIgetVoxVal(mri_slice_edited, x-1, y, 0, 0)) == LABEL_IN_CC) &&
            (nint(MRIgetVoxVal(mri_slice_edited, x-2, y, 0, 0)) == LABEL_IN_CC) &&
            (nint(MRIgetVoxVal(mri_slice_edited, x-3, y, 0, 0)) == LABEL_IN_CC))
        {
          MRIsetVoxVal(mri_slice_edited, x, y, 0, 0, LABEL_IN_CC) ;
        }
      }


    int min_x_fornix_set = 0;
    for (x = xmin ; x <= xmax ; x++)
    {
      int val2, y1 ;

      for (y = 0 ; y < mri_slice->height ; y++)
      {
        if (x == Gx && y == Gy)
        {
          DiagBreak() ;
        }
        val = (int)MRIgetVoxVal(mri_slice, x, y, 0, 0) ;
        val2 = (int)MRIgetVoxVal(mri_slice_edited, x, y, 0, 0) ;
        if (val == LABEL_IN_CC && val2 == LABEL_ERASE)  // in fornix
        {
          thickness = 0 ;
          min_x_fornix = MIN(min_x_fornix, x) ;
	  min_x_fornix_set = 1;
          first_on = last_on = -1 ;
          for (y1 = y ; y1 >= 0 ; y1--)
          {
            val = (int)MRIgetVoxVal(mri_slice_edited, x, y1, 0, 0) ;
            if (val == LABEL_IN_CC)
            {
              if (first_on < 0)
              {
                first_on = y1 ;
              }
              else
              {
                last_on = y1 ;
              }
            }
          }
          thickness = first_on - last_on + 1 ;
          mean_thickness += thickness ;
          std_thickness += (thickness*thickness) ;
          num++ ;
          break ;   // don't bother with this x column any more
        }
      }
    }

    if(!min_x_fornix_set) printf("WARNING: min_x_fornix not set\n");
    else printf("min_x_fornix = %d\n",min_x_fornix);

    mean_thickness /= num ;
    std_thickness = sqrt(std_thickness/num - mean_thickness*mean_thickness);
    max_thickness = (int)ceil(mean_thickness+3*std_thickness); ;
    for (x = min_x_fornix - 3 ; x <= min_x_fornix ; x++)
    {
      first_on = last_on = first_off = -1 ;
      for (y = mri_slice->height-1 ; y >= 0 ; y--)
      {
        val = (int)MRIgetVoxVal(mri_slice_edited, x, y, 0, 0) ;
        if (val == LABEL_IN_CC)
        {
          if (first_on < 0)
          {
            first_on = y ;
          }
          else
          {
            last_on = y ;
          }
        }
        else if (first_on >= 0 && first_off < 0)
        {
          first_off = y ;
        }
      }
      thickness = first_on - last_on + 1 ;
      if (thickness >= max_thickness)  // erase inferior stuff
      {
        int y1 ;

        DiagBreak() ;
        if (first_off >
            last_on)  // two vertical runs of white - leftover fornix
        {
          ystart = first_off ;
        }
        else
        {
          ystart = last_on+max_thickness-1 ;
        }

        // boundary check: there was one oddball subject (BM)
        // where ystart was -2.  this handles that case:
        if (ystart < 0)
        {
          ystart=0;
        }

        for (y1 = ystart ; y1 < mri_slice->height ; y1++)
        {
          val = MRIgetVoxVal(mri_slice_edited, x, y1, 0, 0) ;
          if (val == LABEL_IN_CC)
          {
            MRIsetVoxVal(mri_slice_edited, x, y1, 0, 0, LABEL_ERASE) ;
          }
        }
      }
    }
  }

  MRIreplaceValues(mri_slice_edited, mri_slice_edited, LABEL_ERASE, 0) ;

  MRIfree(&mri_tmp) ;

  // remove unconnected stuff
  mseg = MRIsegment(mri_slice_edited, 1, 255) ;
  i = MRIsegmentMax(mseg) ;
  MRIclear(mri_slice_edited) ;
  MRIsegmentToImage(mri_slice, mri_slice_edited, mseg, i) ;
  MRIsegmentFree(&mseg) ;

  return(mri_slice_edited) ;
}

static MRI *
find_voxels_close_to_both_hemis(MRI *mri_aseg, int lh_label, int rh_label, int wsize)
{
  int  x, y, z, xk, yk, zk, xi, yi, zi, num_lh, num_rh, whalf, label ;
  MRI  *mri_dst ;

  mri_dst = MRIclone(mri_aseg, NULL) ;
  whalf = (wsize-1)/2 ;
  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        num_lh = num_rh= 0 ;
        for (xk = -whalf ; xk <= whalf ; xk++)
        {
          xi = mri_aseg->xi[x+xk] ;
          for (yk = -whalf ; yk <= whalf ; yk++)
          {
            yi = mri_aseg->yi[y+yk] ;
            for (zk = -whalf ; zk <= whalf ; zk++)
            {
              zi = mri_aseg->zi[z+zk] ;
              label = MRIgetVoxVal(mri_aseg, xi, yi, zi, 0) ;
              if (label == lh_label)
              {
                num_lh++ ;
              }
              else if (label == rh_label)
              {
                num_rh++ ;
              }
              if (num_lh > 0 && num_rh > 0)
              {
                break ;
              }
            }
            if (num_lh > 0 && num_rh > 0)
            {
              break ;
            }
          }
          if (num_lh > 0 && num_rh > 0)
          {
            break ;
          }
        }
        if (num_lh > 0 && num_rh > 0)
        {
          MRIsetVoxVal(mri_dst, x, y, z, 0, 1) ;
        }
      }

  return(mri_dst) ;
}

