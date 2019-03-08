/**
 * @file  mri_cc_medial_axis.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


///////////////////////////////////////////
// mri_cc_medial_axis.c
//
// written by Peng Yu
// date: 01/27/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:14 $
// Revision       : $Revision: 1.4 $
////////////////////////////////////////////

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
//#include "volume_io/geom_structs.h"
#include "transform.h"
#include "talairachex.h"
#include "matrix.h"
#include "mriTransform.h"

//static char vcid[] = "$Id: mri_cc_medial_axis.c,v 1.4 2011/03/02 00:04:14 nicks Exp $";


int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
static char     *wmvolume = "mri/wm" ;
static char     *talvolume = "fatalmodthr.img" ;
static int      WM = 0;
static int      TAL = 1;
static int      ROT = 1;
static int      NO_CORRECTION = 0;
static int      alert=0;
const char            *Progname ;
MRI             *mri_wm, *mri_cc_tal ;
MRI             *mri_intensity ;
int             dxi=0;
int             x_edge=0, y_edge=0;


typedef struct medial_axis_type_ {
  float x,y;            /* curr position */
  float nx,ny;         /* curr normal */
  float dx, dy ;     /* current change in position */
  float odx, ody ;  /* last change of position (for momentum) */
  float radius ;
}
atom_type, MEDATOM ;


static double cc_tal_x = 64 ;
static double cc_tal_y = 55 ;
static double cc_tal_z = 64 ;
static LTA *lta = 0;
static int rotatevolume();
static int find_mid_plane_wm(char *subject) ;
static int find_mid_plane_tal(char *subject, char *volume) ;

static MEDATOM *find_medial_axis(MEDATOM *medial_axis, int length, int n_sample, char *ofname) ;
static MRI     *sample_medial_axis(MEDATOM *medial_axis, MRI *cc_medial_axis, int length) ;
static MEDATOM *medial_axis_initialization(MEDATOM *medial_axis, int length) ;
static MEDATOM *medial_axis_integration(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length) ;
static MRI     *smooth_cc(MRI *mri_filled) ;
static MRI     *find_cc_boundary(MRI *cc_smoothed, MRI *cc_boundary) ;
static MEDATOM *compute_gradient_term(MEDATOM *medial_axis,MRI *cc_boundary,MRI *cc_smoothed,int length, float p);
static MEDATOM *compute_smooth_term(MEDATOM *medial_axis,int length, float k);
static MEDATOM *integrate_momentum(MEDATOM *medial_axis, MRI *cc_boundary,MRI *cc_smoothed, int length);
static double  compute_energy(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length);
static MEDATOM *compute_normal(MEDATOM *medial_axis, int length) ;
static MEDATOM *medial_axis_refine_leaf(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length);
static float    compute_error(MEDATOM *medial_axis,MRI *cc_boundary, MRI *cc_smoothed, int length, int new_length, int flag);
static float    sample_distance_map(MRI *cc_boundary, MRI *cc_smoothed, float x, float y) ;
static int      find_cc_slice(MRI *mri, double *pccx, double *pccy, double *pccz, const LTA *lta) ;
static int      find_corpus_callosum(MRI *mri, double *ccx, double *ccy, double *ccz, const LTA *lta) ;
static MRI     *cc_slice_correction(MRI *mri_filled, int xv, int yv, int zv);
static int      edge_detection(MRI *mri_temp, int edge_count,int signal);

static int labels[] = {
                        THICKEN_FILL, NBHD_FILL, VENTRICLE_FILL, DIAGONAL_FILL, DEGENERATE_FILL
                      };
#define NLABELS  sizeof(labels) / (sizeof(labels[0]))
#define MAX_SLICES        41
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


double findMinSize(MRI *mri) {
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
    minsize = (ysize > zsize) ? zsize : ysize;
  else
    minsize = (zsize > xsize) ? xsize : zsize;

  return minsize;
}



int
main(int argc, char *argv[]) {
  char        *cp, fname[STRLEN];
  int         nargs, msec, n_sample = 100, length = 100 ;
  Timer then ;
  MRI         *cc_medial_axis;
  MEDATOM     *medial_axis;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <subject name> <input volume>", Progname);

  cp = getenv("SUBJECTS_DIR");
  if (cp==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }

  then.reset() ;

  /* starting the searching of mid-saggital plane */
  if ( WM == 1 ) {
    find_mid_plane_wm(argv[1]) ;
  } else if ( TAL==1 ) {
    talvolume = argv[2];
    find_mid_plane_tal(argv[1], argv[2]) ;
  }

  /* find the medial axis */
  medial_axis = (MEDATOM *)calloc(length, sizeof(MEDATOM));
  sprintf(fname,"%s/%s/%s_ma.txt",cp,argv[1],argv[2]) ;
  find_medial_axis(medial_axis, length, n_sample, fname) ;
  cc_medial_axis = MRIcopy(mri_cc_tal, NULL) ;
  MRIvalueFill(cc_medial_axis, 0) ;
  sample_medial_axis(medial_axis, cc_medial_axis, length) ;
  sprintf(fname,"%s/%s/%s_med_axis.mgh",cp,argv[1],argv[2]) ;
  fprintf(stdout, "writing medial axis to %s...\n", fname) ;
  MRIwrite(cc_medial_axis, fname) ;

  MRIfree(&mri_cc_tal) ;
  MRIfree(&mri_intensity);
  MRIfree(&cc_medial_axis);
  free(medial_axis);
  msec = then.milliseconds() ;
  fprintf(stdout, "corpus callosum matter segmentation and processing took %2.1f minutes\n", (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}

static int find_mid_plane_wm(char *subject) {
  char        ifname[200], ofname[200],  data_dir[400], *cp ;
  int         y, z, xi, yi_low=256, yi_high=0, zi_low=256, zi_high=0, temp;
  double      xc,yc,zc;
  MATRIX      *mrot, *mtrans;
  int         i, j, k;
  MRI         *mri_talheader, *mri_header, *mri_cc, *mri_tal;
  double      xv, yv, zv;
  FILE        *fp;
  LTA         *lta2 = 0;
  float       volume[5];

  cp = getenv("SUBJECTS_DIR");
  if (cp==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }
  strcpy(data_dir, cp) ;

  sprintf(ifname,"%s/cc_volume_%d.txt",data_dir,dxi) ;
  if ((fp = fopen(ifname, "w")) == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "cc volume measurement: file %s does not exist!", ifname));
  }
  print("writing results to %s\n",ifname);

  sprintf(ifname,"%s/%s/%s",data_dir,subject,wmvolume) ;
  print("reading white matter volume from %s\n", ifname);
  mri_wm = MRIread(ifname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, ifname) ;

  sprintf(ifname,"%s/%s/mri/transforms/talairach.xfm",data_dir,subject) ;
  lta = LTAreadEx(ifname);
  if (lta==0)
    ErrorExit(ERROR_BADPARM,"ERROR: cound not load lta from %s.\n", ifname);
  fprintf(stdout, "INFO: Using %s and its offset for Talairach volume ...\n", ifname);

  for (i = 0 ; i < NLABELS ; i++) {
    MRIreplaceValues(mri_wm, mri_wm, labels[i], 0) ;
  }

  mri_talheader = MRIallocHeader(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_wm, mri_talheader); // not allocate memory, though

  ModifyTalairachCRAS(mri_talheader, lta);

  mri_tal = MRIalloc(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_talheader, mri_tal);
  // now fill the talairach volume values
  MRItoTalairachEx(mri_wm, mri_tal, lta);

  // binalize the talairach volume (mri_tal)
  MRIbinarize(mri_tal, mri_tal, DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;

  sprintf(ofname,"%s/%s/mri/wm_tal.mgh",cp,subject) ;
  fprintf(stdout, "writing talairach volume to %s...\n", ofname) ;
  MRIwrite(mri_tal, ofname) ;

  //find the transform matrix
  mtrans = MatrixAlloc(4, 4, MATRIX_REAL) ;
  mrot = MatrixAlloc(4, 4, MATRIX_REAL) ;

  //try method 2 to get the rotation matrix
  sprintf(ifname,"%s/%s/mri/transforms/talairach.xfm",data_dir,subject) ;
  lta2 = LTAreadEx(ifname);
  mtrans=lta2->xforms[0].m_L;
  Trns_ExtractRotationMatrix (mtrans,mrot);
  *MATRIX_RELT(mrot, 1, 4) = mtrans->rptr[1][4];
  *MATRIX_RELT(mrot, 2, 4) = mtrans->rptr[2][4];
  *MATRIX_RELT(mrot, 3, 4) = mtrans->rptr[3][4];
  lta2->xforms[0].m_L=mrot;

  //rotation wm volume to be upright, using cc volume temporarily
  mri_header = MRIallocHeader(mri_wm->width, mri_wm->height, mri_wm->depth, mri_wm->type);
  MRIcopyHeader(mri_wm, mri_header);
  ModifyTalairachCRAS(mri_header, lta2);
  mri_cc = MRIcopy(mri_wm, NULL) ;
  MRIcopyHeader(mri_header, mri_cc);
  MRItoTalairachEx(mri_wm, mri_cc, lta2);
  // binalize the rotated wm  volume (mri_cc)
  MRIbinarize(mri_cc, mri_cc, DEFAULT_DESIRED_WHITE_MATTER_VALUE/2-1, 0, 110) ;
  sprintf(ofname,"%s/%s/mri/wm.mgh",cp,subject) ;
  fprintf(stdout, "writing rotated white matter volume to %s...\n", ofname) ;
  MRIwrite(mri_cc, ofname) ;


  //now start cc segmentation in talairach space
  mri_cc_tal = MRIcopy(mri_tal, NULL) ;
  MRIcopyHeader(mri_talheader, mri_cc_tal);
  MRIvalueFill(mri_cc_tal, 0) ;

  //most of the work is done in find_corpus_callosum function
  find_corpus_callosum(mri_tal,&cc_tal_x,&cc_tal_y,&cc_tal_z, lta);
  sprintf(ofname,"%s/%s/cc_tal.mgh",cp,subject) ;
  fprintf(stdout, "writing output to %s...\n", ofname) ;
  MRIwrite(mri_cc_tal, ofname) ;

  //starting the volume measurement
  volume[0] = 0.0;
  volume[1] = 0.0;
  volume[2] = 0.0;
  volume[3] = 0.0;
  volume[4] = 0.0;

  //transform cc volume from talairach space to original space
  MRIfromTalairachEx(mri_cc_tal, mri_wm, lta);
  // binalize the rotated cc volume (mri_wm)
  MRIbinarize(mri_wm, mri_wm, CC_VAL/2-1, 0, 100) ;
  sprintf(ofname,"%s/%s/mri/cc_org.mgh",cp,subject) ;
  fprintf(stdout, "writing output to %s...\n", ofname) ;
  MRIwrite(mri_wm, ofname) ;

  //trying to find the position of mid-sagital plane
  MRIcopyHeader(mri_talheader, mri_cc_tal);
  MRItalairachVoxelToVoxelEx(mri_cc_tal, cc_tal_x, cc_tal_y, cc_tal_z, &xv, &yv, &zv, lta) ;

  //rotate the cc volume by the rotation matrix calculated above
  MRIcopyHeader(mri_header, mri_cc);
  MRItoTalairachEx(mri_wm, mri_cc, lta2);
  // binalize the rotated cc volume (mri_cc)
  MRIbinarize(mri_cc, mri_cc, CC_VAL/2-1, 0, 100) ;

  MRIvoxelToTalairachVoxelEx(mri_cc, xv, yv, zv, &xc, &yc, &zc, lta2) ;

  //find the mid-sagital plane there
  xi=nint(xc);
  fprintf(stdout,"cc center is found at %d %d %d\n",xi, nint(yc),nint(zc));

  //find the bounding box
  for (y = 0 ; y < mri_cc->height ; y++) {
    for (z = 0 ; z < mri_cc->depth ; z++) {
      if ( MRIvox(mri_cc, xi, y, z) ) {
        if (y < yi_low)
          yi_low = y ;
        if (z < zi_low)
          zi_low = z ;
        if (y > yi_high )
          yi_high = y ;
        if (z > zi_high)
          zi_high = z ;
      }
    }
  }

  //count the number if equally segmented five parts
  for (i = xi-dxi ; i <= xi+dxi ; i++) {
    for ( j = 0 ; j <= 255 ; j++) {
      for ( k = 0 ; k <= 255 ; k++) {
        if ( MRIvox(mri_cc, i, j, k)>0) {
          if ( k>=zi_low-10 ) {
            temp = floor((k-zi_low)/((zi_high-zi_low+1)/5));
            if (temp < 0) temp = 0;
            if (temp >= 5) temp = 4;
            volume[temp] +=1 ;
            MRIvox(mri_cc, i, j, k)=(temp+1)*20+10;
          }
        }
      }
    }
  }


  fprintf(fp, "%s %d %d %d %d %d %d %d %d %d \n", subject, nint(volume[4]), nint(volume[3]),nint(volume[2]), nint(volume[1]),nint(volume[0]), yi_low, yi_high, zi_low, zi_high);
  fprintf(stdout, "%s %d %d %d %d %d %d %d %d %d \n", subject, nint(volume[4]), nint(volume[3]),nint(volume[2]), nint(volume[1]),nint(volume[0]), yi_low, yi_high, zi_low, zi_high);

  sprintf(ofname,"%s/%s/mri/cc.mgh",cp,subject) ;
  fprintf(stdout, "writing output to %s...\n", ofname) ;
  MRIwrite(mri_cc, ofname) ;
  MRIfree(&mri_wm);
  MRIfree(&mri_tal);
  MRIfree(&mri_cc);
  MRIfree(&mri_talheader);
  MRIfree(&mri_header);
  MatrixFree(&mtrans);
  MatrixFree(&mrot);
  fclose(fp);

  return(NO_ERROR) ;
}



static int find_mid_plane_tal(char *subject, char *volume) {
  char        ifname[200], ofname[200],  data_dir[400], *cp ;
  MRI *mri_temp, *mri_filled, *mri_slice;
  int ii, jj, area;

  /* finish reading in the original volume registered with talaraich space*/
  cp = getenv("SUBJECTS_DIR");
  if (cp==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }
  strcpy(data_dir, cp) ;

  sprintf(ifname,"%s/%s/%s",data_dir,subject,talvolume) ;
  print("reading talairach transformed volume from %s\n", ifname);
  mri_intensity = MRIread(ifname) ;
  if (!mri_intensity)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, ifname) ;

  /*Rotate the original volume to find the slice with minimal white matter*/
  rotatevolume();

  if (NO_CORRECTION) {
    sprintf(ifname,"%s/%s/fatal.mgh",data_dir,subject) ;
    print("reading talairach transformed fa volume from %s\n", ifname);
    mri_intensity = MRIread(ifname) ;
  } else {
    sprintf(ofname,"%s/%s/%s_rotated.mgh",data_dir,subject,talvolume) ;
    fprintf(stdout, "writing rotated tal volume to %s...\n", ofname) ;
    MRIwrite(mri_intensity, ofname) ;

    /*now segment out corpus callosum in talairach space*/
    mri_slice = MRIextractPlane(mri_cc_tal, NULL, MRI_SAGITTAL, cc_tal_x);
    mri_filled =  MRIfillFG(mri_slice, NULL, cc_tal_z, cc_tal_y,0,WM_MIN_VAL,CC_VAL,&area);

    /*if needed, rotate the mri_filled 90 degree counterclockwise*/
    if (ROT) {
      mri_temp=MRIcopy(mri_filled,NULL);
      MRIvalueFill(mri_temp, 0) ;
      for (ii = 0 ; ii < mri_temp->width ; ii++) {
        for (jj = 0 ; jj < mri_temp->height ; jj++) {
          MRIvox(mri_temp,ii,jj,0) = MRIvox(mri_filled,mri_filled->width-jj-1,ii,0);
        }
      }
      //sprintf(fname,"/autofs/space/windu_006/users/salat/FR_1004/LEIK76_reconMS/callosum/cc_tal.mgh") ;
      //fprintf(stdout, "writing output to %s...\n", fname) ;
      //MRIwrite(mri_temp, fname) ;
      //used the rotated version for correction
      mri_temp = cc_slice_correction(mri_temp,cc_tal_x,cc_tal_z,cc_tal_y);
      //rotate it back and re-measure the area
      for (ii = 0 ; ii < mri_filled->width ; ii++) {
        for (jj = 0 ; jj < mri_filled->height ; jj++) {
          MRIvox(mri_filled,ii,jj,0) = MRIvox(mri_temp,jj,mri_temp->width-ii-1,0);
        }
      }
      MRIfree(&mri_temp) ;
    } else {
      //sprintf(fname,"/autofs/space/windu_006/users/salat/FR_1004/LEIK76_reconMS/callosum/cc_tal.mgh") ;
      //fprintf(stdout, "writing output to %s...\n", fname) ;
      //MRIwrite(mri_filled, fname) ;
      mri_filled = cc_slice_correction(mri_filled,cc_tal_x,cc_tal_y,cc_tal_z);
    }
    MRIvalueFill(mri_cc_tal, 0) ;
    MRIfillPlane(mri_filled, mri_cc_tal, MRI_SAGITTAL, cc_tal_x, CC_VAL);
    MRIfree(&mri_filled) ;
    MRIfree(&mri_slice) ;
  }

  sprintf(ofname,"%s/%s/%s_cc_tal.mgh",cp,subject,volume) ;
  fprintf(stdout, "writing segmented cc volume to %s...\n", ofname) ;
  MRIwrite(mri_cc_tal, ofname) ;

  return(NO_ERROR) ;
}

static int
rotatevolume() {
  int   width, height, depth,counter, number;
  int   y_rotate, z_rotate, x, y_alpha, z_alpha, xi, yi, zi, yk, zk ;
  int  whalf, p_x, p_y, p_z;
  int   x_center, y_center, z_center, xshift, i, j, k ;
  double      val=0;
  float       ey_x, ey_y, ey_z, ez_x, ez_y, ez_z, x1_x, x1_y, x1_z, xbase, ybase, zbase;
  /*  BUFTYPE     *pdst, *pptr ; */
  MRI   *mri_temp, *mri_tal;

  mri_tal=MRIcopy(mri_intensity,NULL);

  for (k = 0 ; k < mri_tal->depth ; k++) {
    for (j = 0 ; j < mri_tal->height ; j++) {
      for (i = 0 ; i < mri_tal->width ; i++) {
        MRIsampleVolumeFrameType(mri_tal, i, j, k, 0, SAMPLE_NEAREST, &val) ;
        if (val <= 0)
          val = 0 ;
        else {
          if (NO_CORRECTION) cc_tal_x = i;
          val = 110 ;
        }
        MRIsetVoxVal(mri_tal, i, j, k, 0, val) ;
      }
    }
  }

  /*Done with the construction of a binary volume mri_tal*/
  mri_tal = MRIchangeType(mri_tal, 0, 0, 0, 1) ;
  mri_cc_tal = MRIcopy(mri_tal, NULL) ;

  if (!NO_CORRECTION) {
    MRIvalueFill(mri_cc_tal, 0) ;

    width = mri_tal->width ;
    height = mri_tal->height ;
    depth = mri_tal->depth ;
    xshift = 3;
    y_rotate = 5;
    z_rotate = 5;
    whalf = mri_tal->width/2;
    p_x = 0;
    p_y = 0;
    p_z = 0;
    number = whalf*whalf;
    y_alpha = 0 ;
    z_alpha = 0 ;

    y_center = cc_tal_y ;
    z_center = cc_tal_z ;

    for ( x = -xshift;   x <= xshift  ; x++)  /* shift centoid point in x direction */
      for ( y_alpha = -1*y_rotate; y_alpha <= y_rotate; y_alpha++)
        for ( z_alpha = -1*z_rotate; z_alpha <= z_rotate; z_alpha++) {
          x_center = nint(cc_tal_x) + x;
          counter = 0;
          x1_x = cos( (double) y_alpha*M_PI/180) * cos( (double) z_alpha*M_PI/180 );  /* new x axis */
          x1_y = sin( (double) z_alpha*M_PI/180 ) ;
          x1_z = cos( (double) y_alpha*M_PI/180 ) * sin( (double) z_alpha*M_PI/180 ) ;

          ez_x = sin( (double) y_alpha*M_PI/180) ;  /* basis vectors for plane */
          ez_y = 0 ;
          ez_z = cos( (double) y_alpha*M_PI/180 );

          ey_x = ez_y*x1_z - ez_z*x1_y ;
          ey_y = ez_z*x1_x - ez_x*x1_z ;
          ey_z = ez_x*x1_y - ez_y*x1_x ;

          /* now find all the voxel in this plane */

          for (yk = -2*whalf ; yk <= 2*whalf ; yk++) {
            xbase = (float)x_center + (float)yk * ey_x ;
            ybase = (float)y_center + (float)yk * ey_y ;
            zbase = (float)z_center + (float)yk * ey_z ;

            for (zk = -2*whalf ; zk <= 2*whalf ; zk++) {
              /* in-plane vect. is linear combination of scaled basis vects */
              xi = nint(xbase + zk*ez_x) ;
              yi = nint(ybase + zk*ez_y) ;
              zi = nint(zbase + zk*ez_z) ;

              if ( (xi>=0) && (xi<width) && (yi>=0) && (yi<height) && (zi>=0) && (zi<depth))
                if ( MRIvox(mri_tal, xi, yi, zi)>0 ) counter++;
            }
          }
          if ( counter < number ) {
            number = counter;
            p_x = x_center;
            p_y = y_alpha;
            p_z = z_alpha;
          }
        }

    fprintf(stdout, "counts: %d rotation x center: %d y_alpha: %d z_alpha: %d \n", number, p_x, p_y, p_z) ;

    x_center = p_x ;
    x1_x = cos( (double) p_y*M_PI/180) * cos( (double) p_z*M_PI/180 );  /* new x axis */
    x1_y = sin( (double) p_z*M_PI/180) ;
    x1_z = cos( (double) p_y*M_PI/180) * sin( (double) p_z*M_PI/180 ) ;

    ez_x = sin( (double) p_y*M_PI/180) ;  /* basis vectors for plane */
    ez_y = 0 ;
    ez_z = cos( (double) p_y*M_PI/180 );

    ey_x = ez_y*x1_z - ez_z*x1_y ;
    ey_y = ez_z*x1_x - ez_x*x1_z ;
    ey_z = ez_x*x1_y - ez_y*x1_x ;

    mri_temp=MRIcopy(mri_intensity,NULL);
    MRIvalueFill(mri_temp, 0) ;

    /* now find all the voxel in this plane */
    for ( x=0; x < width; x++)
      for (yk = -whalf ; yk < whalf ; yk++) {
        xbase = (float)x + (float)yk * ey_x ;
        ybase = (float)y_center + (float)yk * ey_y ;
        zbase = (float)z_center + (float)yk * ey_z ;
        for (zk = -whalf ; zk < whalf ; zk++) {
          /* in-plane vect. is linear combination of scaled basis vects */
          xi = nint(xbase + zk*ez_x) ;
          yi = nint(ybase + zk*ez_y) ;
          zi = nint(zbase + zk*ez_z) ;
          if  ( (xi>=0) && (xi<width) && (yi>=0) && (yi<height) && (zi>=0) && (zi<depth)) {
            if (x==x_center)
              MRIvox(mri_cc_tal, x, yk+y_center, zk+z_center) = MRIvox(mri_tal, xi, yi, zi);
            MRIvox(mri_temp, x, yk+y_center, zk+z_center) = MRIvox(mri_intensity, xi, yi, zi);
          }
        }
      }
    cc_tal_x = x_center;
    mri_intensity = MRIcopy(mri_temp,NULL);
    MRIfree(&mri_temp);
  }

  MRIfree(&mri_tal);

  return(NO_ERROR) ;
}

#define CC_SPREAD       10
#define MIN_THICKNESS   3
#define MAX_THICKNESS   20

static int
find_corpus_callosum(MRI *mri_tal, double *pccx, double *pccy, double *pccz, const LTA *lta) {
  int         xv, yv, zv, max_y, max_thick=0, thickness=0, y1, xcc, ycc, x, y,x0, extension=50 ;
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
  if (mri_tal->linear_transform || lta) {
    MRIworldToVoxel(mri_tal, 0.0, 0.0, 0.0, &xr, &yr, &zr);   /* everything is now in tal coords */
    xv = nint(xr) ;
    yv = nint(yr) ;
    zv = nint(zr) ;
  } else {
    xv = x0;
    yv = region.y+region.dy/2;
    zv = region.z+region.dz/2;
  }

  fprintf(stdout, "original seed found at x=%d, y=%d z=%d \n", x0, yv, zv );
  /* find the column with the lowest starting y value of any sign. thick. */
  xcc = ycc = max_y = 0 ;
  for (x = x0-cc_spread ; x <= x0+cc_spread ; x++) {
    /* search for first non-zero pixel */
    // in the talairach origin coronal slice from the top
    while (thickness==0 && yv-extension >= 0 && yv+extension <= 256) {
      for (y = yv-extension ; y < yv+extension ; y++) {
        if (MRIvox(mri_tal, x, y, zv) >= WM_MIN_VAL)
          break ;
      }
      // find y which is greater than WM_MIN_VAL
      /* check to make sure it as reasonably thick */
      if (y < yv+extension ) // within the region and bigger than so far
      {
        for (y1 = y, thickness = 0 ; y1 < slice_size ; y1++, thickness++)
          if (!MRIvox(mri_tal, x, y1, zv)) // if becomes zero, then break out
            break ;
        if ( thickness > min_thickness && thickness < max_thickness ) {
          if ( y > max_y || (y == max_y && thickness > max_thick) ) {
            // found the zero voxel at y1 -> thinckness
            xcc = x ;
            ycc = y+thickness/2 ;  /* in middle of cc */
            max_y = y ;            // mark starting y position
            max_thick = thickness;
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
              fprintf(stdout, "potential cc found at (%d, %d), thickness = %d\n",
                      xcc, ycc, thickness) ;
          } else if ( y==max_y && thickness== max_thick && (x==xcc+1 || flag ==1 )) {
            flag = 1;
            counts ++;
          } else if (flag == 1) {
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
    return(ERROR_BADPARM) ;

  /* now convert the in-plane coords to Talairach coods */
  MRIvoxelToWorld(mri_tal, xcc, ycc, zv, pccx, pccy, pccz) ;
  fprintf(stdout, "%d, %d, %d\n", xcc, ycc, zv);

  find_cc_slice(mri_tal, pccx, pccy, pccz, lta) ;

  return(NO_ERROR) ;
}


static int
find_cc_slice(MRI *mri_tal, double *pccx, double *pccy, double *pccz, const LTA *lta) {
  // here we can handle only up to .5 mm voxel size
  int         area[MAX_SLICES*2], flag[MAX_SLICES*2], min_area, min_slice, slice, offset,xv,yv,zv,
  xo, yo ,i, total_area=0, left=0, right=0;
  MRI         *mri_slice, *mri_filled ;
  double      aspect, x_tal, y_tal, z_tal, x, y, z, xvv, yvv, zvv;
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
    half_slices = 1;

  x_tal = *pccx ;
  y_tal = *pccy ;
  z_tal = *pccz ;
  offset = 0 ;
  xo = yo = (slice_size-1)/2 ;  /* center point of the slice */
  for (slice = 0 ; slice < max_slices ; slice++) {
    offset = slice - half_slices ;

    i=0;
    area[slice]=0;
    while (area[slice]<100 && i<=5) {
      x = x_tal + offset ;
      y = y_tal ;
      z = z_tal ;
      MRIworldToVoxel(mri_tal, x, y,  z, &x, &y, &z) ;
      xv = nint(x) ;
      yv = nint(y)-i ;
      zv = nint(z) ;
      mri_slice = MRIextractPlane(mri_tal, NULL, MRI_SAGITTAL, xv);
      mri_filled =  MRIfillFG(mri_slice, NULL, zv, yv,0,WM_MIN_VAL,CC_VAL,&area[slice]);
      MRIboundingBox(mri_filled, 1, &region) ;
      aspect = (double)region.dy / (double)region.dx ;
      if (i++) fprintf(stdout,"moved %d in slice %d \n", i-1, slice);
    }
    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= slice_size -1 ||
        region.y+region.dy >= slice_size-1)
      area[slice] = 0 ;

    if ( !(area[slice]>1100&&(nint(y)-region.y>11)) || region.dy>=3.5*(y-region.y) ) {
      mri_filled = cc_slice_correction(mri_filled,xv,yv,zv);

      area[slice] = 0;
      flag[slice] = 0;

      for (ii = 0 ; ii < mri_filled->width ; ii++) {
        for (jj = 0 ; jj < mri_filled->height ; jj++) {
          if ( MRIvox(mri_filled, ii, jj, 0)>0 ) area[slice]++;
        }
      }
    } else {
      flag[slice] = 1;
      //if (offset>=-5&&offset<0) left++;
      //else if (offset<=5&&offset>0) right++;
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout, "slice[%d] @ (%d, %d, %d): area = %d\n",
              slice, xv, yv, zv, area[slice]) ;

    if ((Gdiag & DIAG_WRITE) && !(slice % 1) && DIAG_VERBOSE_ON) {
      sprintf(fname, "cc_slice%d.mgh", slice);
      MRIwrite(mri_slice, fname) ;
      sprintf(fname, "cc_filled%d.mgh", slice);
      MRIwrite(mri_filled, fname) ;
    }
    MRIfillPlane(mri_filled, mri_cc_tal, MRI_SAGITTAL, xv, CC_VAL);

    MRIfree(&mri_filled) ;
    MRIfree(&mri_slice) ;
  }

#if 0
  min_area = 10000 ;
  min_slice = -1 ;
  for (slice = 1 ; slice < max_slices-1 ; slice++) {
    if (area[slice] < min_area &&
        (area[slice] >= min_cc_area && area[slice] <= max_cc_area)) {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }
#else
  min_area = 10000*5 ;
  min_slice = -1 ;
  for (slice = 6 ; slice <= 14 ; slice++) {
    for (i=-2, total_area =0; i <=2; i++)
      total_area += area[slice+i];

    if (total_area < min_area &&
        (total_area >= min_cc_area*5 && total_area <= max_cc_area*5)) {
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
  for (slice = 0 ; slice < max_slices ; slice++) {
    if ( (slice-(half_slices+offset)>=-6) && (slice<half_slices+offset) && flag[slice]==1 ) left++;
    else if ( (slice-(half_slices+offset)<=6) && (slice>half_slices+offset) && flag[slice]==1 ) right++;
  }
  offset = offset+left-right;
  fprintf(stdout, "find offset as %d using shifting\n", left-right);
  if (abs(offset)>=5) offset = 5*(offset/abs(offset));
  *pccx = x = x_tal+floor(offset) ;
  *pccy = y = y_tal ;
  *pccz = z = z_tal ;
  // just for debugging
  MRIworldToVoxel(mri_tal, x, y,  z, &xvv, &yvv, &zvv) ;
  *pccx = xvv ;
  *pccy = yvv ;
  *pccz = zvv ;

  fprintf(stdout, "updating initial cc seed to Tal vol (%.2f, %.2f, %.2f) TAL (%.2f, %.2f, %.2f)\n",
          xvv, yvv, zvv, x, y, z);

  return(NO_ERROR) ;
}


static MRI *
cc_slice_correction(MRI *mri_filled, int xv, int yv, int zv) {
  int    x, y, xi_low=255, xi_high=0, yi_low=255, yi_high=0, edge_count = 0, length=0, ratio = 1;
  int    temp=0, temp2=0, old_temp =0;
  int    x1=0, y1=0, x2=0, y2=0, height=mri_filled->height, width=mri_filled->width, flag =0;
  int    i, section1[150];
  MRI *mri_temp1, *mri_temp2;

  mri_temp1=MRIcopy(mri_filled,NULL);
  mri_temp2=MRIcopy(mri_filled,NULL);
  MRIvalueFill(mri_temp1, 0) ;
  MRIvalueFill(mri_temp2, 0) ;

  ratio *= height/128 ;

  for (x = 1 ; x < width-1 ; x++) {
    for (y = 1 ; y < height-1 ; y++) {
      if ( MRIvox(mri_filled, x-1, y, 0) || MRIvox(mri_filled, x, y-1, 0) || MRIvox(mri_filled, x, y+1, 0) || MRIvox(mri_filled, x+1, y, 0) || MRIvox(mri_filled, x, y, 0) )
        MRIvox(mri_temp1,x,y,0)=100;
    }
  }

  for (x = 1 ; x < width-1 ; x++) {
    for (y = 1 ; y < height-1 ; y++) {
      if ( MRIvox(mri_temp1, x-1, y, 0) && MRIvox(mri_temp1, x, y-1, 0) && MRIvox(mri_temp1, x, y+1, 0) && MRIvox(mri_temp1, x+1, y, 0) && MRIvox(mri_temp1, x, y, 0) )
        MRIvox(mri_temp2,x,y,0)=100;
    }
  }


  for (x = 0 ; x < width ; x++) {
    for (y = 0 ; y < height ; y++) {
      if ( MRIvox(mri_temp2, x, y, 0) ) {
        if (x < xi_low)
          xi_low = x ;
        if (y < yi_low)
          yi_low = y ;
        if (x > xi_high )
          xi_high = x ;
        if (y > yi_high)
          yi_high = y ;
      }
    }
  }

  if ( (yi_high>yi_low+25*ratio) ) {
    if ((5*(yv-yi_low)>3*(yi_high-yi_low))) {
      yi_low=yi_high-ratio*20;
      for (x = 0 ; x < xi_high-nint((xi_high-xi_low)/4); x++) {
        for (y = yi_low ; y >0 ; y--) {
          MRIvox(mri_temp2, x, y, 0) = 0;
          MRIvox(mri_temp1, x, y, 0) = 0;
          MRIvox(mri_filled, x, y, 0) = 0;
        }
      }
    } else {
      yi_high=yi_low+15*ratio;
      for (x = 0 ; x < xi_high-nint((xi_high-xi_low)/4); x++) {
        for (y = yi_high ; y < height ; y++) {
          MRIvox(mri_temp2, x, y, 0) = 0;
          MRIvox(mri_temp1, x, y, 0) = 0;
          MRIvox(mri_filled, x, y, 0) = 0;
        }
      }
    }
  }

  //sprintf(fname,"/autofs/space/windu_006/users/salat/FR_1004/LEIK76_reconMS/callosum/cc_dilation.mgh") ;
  //MRIwrite(mri_temp1, fname) ;
  //mri_filled=MRIcopy(mri_temp2,NULL);
  //sprintf(fname,"/autofs/space/windu_006/users/salat/FR_1004/LEIK76_reconMS/callosum/cc_shrinkage.mgh") ;
  //MRIwrite(mri_temp2, fname) ;

  //mri_temp2=MRIcopy(mri_filled,NULL);
  /*find the first edge of the spike */

  x_edge =0;
  y_edge =0;

  /*the first  edge of the spike */
  flag=0;
  for (x_edge=xi_high-nint((xi_high-xi_low)/3.5); x_edge>=xi_low+nint((xi_high-xi_low)/7); x_edge--) {
    edge_count=0;
    y_edge=0;
    while ( edge_count<3 && y_edge<height-2 ) {
      length = edge_detection(mri_temp2,edge_count,0);
      if (length>1) edge_count++ ;
      if (edge_count==1&&length>1) {
        if (length>=10*ratio&&x_edge>xi_low+7*ratio&&x_edge<xi_high-10*ratio ) {
          temp=old_temp;
          if (x1<x_edge) {
            y1=old_temp;
            x1=x_edge;
            for (y=temp; y<mri_filled->width; y++)
              MRIvox(mri_filled, x_edge, y, 0) =0;
          }
          flag = 2;
        } else temp=y_edge;

        if ( length<=1 || y_edge<yi_low+3 ) {
          for (y=y_edge+1; y>0; y--)
            MRIvox(mri_filled, x_edge, y, 0) =0;
          edge_count-=1;
        } else if (length>2&&x_edge>xi_low+8*ratio&&x_edge<xi_high-10*ratio) {
          for (y=y_edge+1; y<mri_filled->width; y++)
            MRIvox(mri_filled, x_edge, y, 0) =0;
        }
      } else if (length>1&&edge_count>1 && y_edge<yi_high+1 && y_edge>yv+7.5*ratio )
        edge_count -=1;
    }
    if (edge_count>=2&&flag==0) flag=1;
    else if (edge_count<=1&&flag==1&&x_edge>xi_low+6) {
      flag = 0;
      y1=old_temp;
      x1=x_edge;
      //   if (x1<zv+6.5*ratio) break;
    }
    old_temp = temp;
  }
  fprintf(stdout, "first point found at %d %d \n", x1, y1);
  x_edge =0;
  y_edge =0;

  /*the second edge of the spike */

  flag=0;
  //for (y_edge=yi_high-nint((yi_high-yi_low)/3); y_edge>=yi_low+2; y_edge--)
  for (y_edge=yv+2; y_edge>=yi_low+2; y_edge--) {
    edge_count=0;
    x_edge=0;
    i=yv+2-y_edge;
    section1[i]=0;
    while (x_edge<width-1) {
      length = edge_detection(mri_temp2,edge_count,1);
      if (length >=2)     edge_count++ ;
      if (edge_count==1) {
        temp=x_edge;
        if (!section1[i]) section1[i]=length;
      }
      if (edge_count==2) temp2=x_edge-length;
    }

    if (edge_count>=3&&flag==0)     flag=1;
    else if (edge_count<=2&&flag==1) {
      flag = 0;
      x2=old_temp;
      y2=y_edge+1;
      if (y2<=yi_low+8*ratio) break;
    } else if ( x2==0&&i>=0&&(section1[i]>=12||(2*section1[i-1]<section1[i]))  ) {
      x2=old_temp;
      y2=y_edge;
      if (y2<=yi_low+8) break;
    }
    if ( (edge_count>=4) && temp2<x1) old_temp =  temp2-1;
    else  old_temp = temp;
  }
  fprintf(stdout, "second point found at %d %d \n", x2, y2);

  if ( x2>0 && x1>xi_low+4*ratio && x1>x2 ) {
    if (x2<x1) {
      temp=x1;
      x1=x2;
      x2=temp;
      temp=y1;
      y1=y2;
      y2=temp;
    }
    for (x=x1; x<=x2; x++) {
      for (y=nint(y1+(y2-y1)*(x-x1)/(x2-x1));y<height;y++) MRIvox(mri_filled, x, y, 0)=0;
    }
  }

  //sprintf(fname,"/space/solo/4/users/recon/ALBJ82_recon/callosum/cc_cut.mgh") ;
  //MRIwrite(mri_filled, fname) ;
  MRIfree(&mri_temp1);
  MRIfree(&mri_temp2);
  return(mri_filled) ;
}


static int
edge_detection(MRI *mri_temp, int edge_count, int signal) {
  int length = 0, gap = 0, width = mri_temp->width;

  if (signal==1) {
    while (gap<=2) {
      gap=0;
      while ( x_edge < width) {
        if (MRIvox(mri_temp, x_edge, y_edge, 0))
          break ;
        x_edge++;
      }

      while ( x_edge < width ) {
        if (!MRIvox(mri_temp, x_edge, y_edge, 0))
          break ;
        else length++ ;
        x_edge++;
      }

      while ( x_edge < width) {
        if (MRIvox(mri_temp, x_edge, y_edge, 0))
          break ;
        else gap++;
        x_edge++;
      }

      if (gap<=2&&x_edge<width) length += gap;
      else {
        x_edge -= gap;
        break;
      }
    }
  } else {
    while (y_edge < width) {
      if (MRIvox(mri_temp, x_edge, y_edge, 0))
        break ;
      y_edge++;
    }

    while ( y_edge < width ) {
      if (!MRIvox(mri_temp, x_edge, y_edge, 0))
        break ;
      else length++ ;
      y_edge++;
    }
  }
  return(length);
}

static MRI *smooth_cc(MRI *mri_filled) {
  int    x, y ;
  int    height=mri_filled->height, width=mri_filled->width;
  MRI *mri_temp1, *mri_temp2;

  mri_temp1=MRIcopy(mri_filled,NULL);
  mri_temp2=MRIcopy(mri_filled,NULL);
  MRIvalueFill(mri_temp1, 0) ;
  MRIvalueFill(mri_temp2, 0) ;

  for (x = 1 ; x < width-1 ; x++) {
    for (y = 1 ; y < height-1 ; y++) {
      if ( MRIvox(mri_filled, x-1, y, 0) || MRIvox(mri_filled, x, y-1, 0) || MRIvox(mri_filled, x, y+1, 0) || MRIvox(mri_filled, x+1, y, 0) || MRIvox(mri_filled, x, y, 0) )
        MRIvox(mri_temp1,x,y,0)=100;
    }
  }

  for (x = 1 ; x < width-1 ; x++) {
    for (y = 1 ; y < height-1 ; y++) {
      if ( MRIvox(mri_temp1, x-1, y, 0) && MRIvox(mri_temp1, x, y-1, 0) && MRIvox(mri_temp1, x, y+1, 0) && MRIvox(mri_temp1, x+1, y, 0) && MRIvox(mri_temp1, x, y, 0) )
        MRIvox(mri_temp2,x,y,0)=100;
    }
  }

  //sprintf(fname, "/space/yoda/5/users/salat/callosum/EVAC79_recon_dti_cc/cc_dilation.mgh");
  //MRIwrite(mri_temp1, fname) ;
  mri_filled=MRIcopy(mri_temp2,NULL);
  MRIwrite(mri_temp2, "/space/solo/4/users/recon/DONF55b_recon/callosum/cc_dilation.mgh") ;
  MRIfree(&mri_temp1);
  MRIfree(&mri_temp2);
  return(mri_filled);
}

static MRI *find_cc_boundary(MRI *cc_smoothed, MRI *cc_boundary) {
  int   width = cc_smoothed->width;
  int   height = cc_smoothed->height;
  int   depth = cc_smoothed->depth;
  int   k, j, i;

  cc_boundary = MRIcopy(cc_smoothed,NULL) ;
  MRIvalueFill(cc_boundary, 0) ;

  for (k=0;k<depth;k++)
    for (j=0;j<height;j++)
      for (i=0;i<width;i++)
        if (MRIvox(cc_smoothed,i,j,k)>0)
          if (MRIvox(cc_smoothed,i,cc_smoothed->yi[j-1],k)==0||MRIvox(cc_smoothed,i,cc_smoothed->yi[j+1],k)==0||MRIvox(cc_smoothed,cc_smoothed->xi[i-1],j,k)==0||MRIvox(cc_smoothed,cc_smoothed->xi[i+1],j,k)==0)
            //||MRIvox(cc_smoothed,cc_smoothed->xi[i-1],cc_smoothed->yi[j-1],k)==0||MRIvox(cc_smoothed,cc_smoothed->xi[i-1],cc_smoothed->yi[j+1],k)==0||MRIvox(cc_smoothed,cc_smoothed->xi[i+1],cc_smoothed->yi[j-1],k)==0||MRIvox(cc_smoothed,cc_smoothed->xi[i+1],cc_smoothed->yi[j+1],k)==0
          {
            MRIvox(cc_boundary,i,j,k)= 100;
          }
  //MRIwrite(cc_boundary,"/autofs/space/windu_006/users/salat/FR_1004/LEIK76_reconMS/callosum/cc_boundary.mgh" ) ;
  return(cc_boundary);
}

static MEDATOM *compute_normal(MEDATOM *medial_axis, int length) {
  int i, j, n;
  float dx, dy, sx, sy, dist;

#if 1

  /*compute the normal based on 2 neighbors*/
  for (i=1; i < length-1; i++) {
    sx=0;
    sy=0;
    dist = 0;
    j=0;
    n=0;
    do {
      j++;
      dx = medial_axis[i+j].x - medial_axis[i].x ;
      dy = medial_axis[i+j].y - medial_axis[i].y ;
      dist = sqrt(dx*dx+dy*dy);
      if (dist) {
        medial_axis[i+j].dx = dx/dist ;
        medial_axis[i+j].dy = dy/dist ;
      } else
        continue;

      dx = medial_axis[i-j].x - medial_axis[i].x ;
      dy = medial_axis[i-j].y - medial_axis[i].y ;
      dist = sqrt(dx*dx+dy*dy);
      if (dist) {
        medial_axis[i-j].dx = dx/dist ;
        medial_axis[i-j].dy = dy/dist ;
      } else
        continue;

      medial_axis[i].dx = medial_axis[i+j].dx - medial_axis[i-j].dx ;
      medial_axis[i].dy = medial_axis[i+j].dy - medial_axis[i-j].dy ;
      dist = sqrt(medial_axis[i].dx*medial_axis[i].dx+medial_axis[i].dy*medial_axis[i].dy);
      if (dist) {
        dx = medial_axis[i].dy/dist;
        dy = -1*medial_axis[i].dx/dist;
        n++;
        sx+=dx;
        sy+=dy;
      }
    } while ((n<2||dist==0)&&(i+j<length-1)&&(i-j>0));

    dist = sqrt(sx*sx+sy*sy);
    if (dist) {
      sx /= dist;
      sy /= dist;
    }
    medial_axis[i].nx = sx;
    medial_axis[i].ny = sy;
  }
#else
  /* compute normal using neighbors 1 pixel away*/

  for (i=1; i < length-1; i++) {
    for (j=1; i+j<(length-1)&&(medial_axis[i+j].y-medial_axis[i].y)<= 0.5 ; j++) {}
    dx = medial_axis[i+j].x - medial_axis[i].x ;
    dy = medial_axis[i+j].y - medial_axis[i].y ;
    dist = sqrt(dx*dx+dy*dy);
    if (dist) {
      medial_axis[i].dx = dx/dist ;
      medial_axis[i].dy = dy/dist ;
    }

    for (j=-1; i+j>=1 && (medial_axis[i].y-medial_axis[i+j].y)<= 0.5 ; j--) {}
    dx = medial_axis[i-j].x - medial_axis[i].x ;
    dy = medial_axis[i-j].y - medial_axis[i].y ;
    dist = sqrt(dx*dx+dy*dy);
    if (dist) {
      medial_axis[i].dx -= dx/dist ;
      medial_axis[i].dy -= dy/dist ;
    }

    dist = sqrt(medial_axis[i].dx*medial_axis[i].dx+medial_axis[i].dy*medial_axis[i].dy);
    if (dist) {
      dx = medial_axis[i].dy/dist;
      dy = -1*medial_axis[i].dx/dist;
    }
    medial_axis[i].nx = dx;
    medial_axis[i].ny = dy;
  }
#endif

  for (i=1; i < length-1; i++) {
    medial_axis[i].dx = 0;
    medial_axis[i].dy = 0;
  }
  return(medial_axis);
}

static MEDATOM *find_medial_axis(MEDATOM *medial_axis, int length, int n_sample, char *ofname) {
  int   i, j, xi_high=0, yi_high=0, n, count;
  MRI   *cc_boundary, *cc_smoothed ;
  int   width, height, yi_low, xi_low, xx, yy;
  float x1, x2, y1, y2, x, y, dx, dy ;
  FILE  *fp;
  double  dist, space, fraction;
  double val = 0, mean_val=0;

  cc_smoothed = MRIextractPlane(mri_cc_tal, NULL, MRI_SAGITTAL, cc_tal_x);
  cc_smoothed = smooth_cc(cc_smoothed);
  //MRIwrite(cc_smoothed, fname) ;
  cc_boundary = find_cc_boundary(cc_smoothed, cc_boundary) ;
  width=cc_boundary->width;
  height=cc_boundary->height;
  yi_low=height;
  xi_low=width;

  /*initialize the leaf nodes first*/
  if (ROT) {
    for (i = 0 ; i < width ; i++) {
      for (j = 0 ; j < height ; j++) {
        if ( MRIvox(cc_boundary, i, j, 0) ) {
          if (i < xi_low)
            xi_low = i ;
          if (j < yi_low) {
            yi_low = j ;
            medial_axis[0].x = i;
            medial_axis[0].y = j;
          }
          if (i > xi_high )
            xi_high = i ;
          if (j > yi_high) {
            yi_high = j ;
            medial_axis[length-1].x = i;
            medial_axis[length-1].y = j;
          }
        }
      }
    }
  } else {
    for (i = 0 ; i < width ; i++) {
      for (j = 0 ; j < height ; j++) {
        if ( MRIvox(cc_boundary, i, j, 0) ) {
          if (i < xi_low) {
            xi_low = i ;
            medial_axis[0].x = i;
            medial_axis[0].y = j;
          }
          if (j < yi_low)
            yi_low = j ;

          if (i > xi_high ) {
            xi_high = i ;
            medial_axis[length-1].x = i;
            medial_axis[length-1].y = j;
          }
          if (j > yi_high)
            yi_high = j ;
        }
      }
    }
  }
#if 1
  medial_axis = medial_axis_initialization(medial_axis, length) ;
  medial_axis = medial_axis_integration(medial_axis, cc_boundary, cc_smoothed, length) ;
  medial_axis = medial_axis_refine_leaf(medial_axis, cc_boundary, cc_smoothed, length) ;
#else
  { float val1, val2;
    length=102;
    fp = fopen(ofname, "r");
    for (j=0; j<length; j++) {
      fscanf(fp, "%d %f %f %f %f %f", &i, &medial_axis[j].x, &medial_axis[j].y, &medial_axis[j].radius, &medial_axis[j].dx, &medial_axis[j].dy);
    }
    fclose(fp);
  }
#endif

  if (0) {
    double xw,yw,zw;
    for (j=0; j<length; j++) {
      MRIvoxelToWorld(mri_cc_tal, cc_tal_x, medial_axis[j].y, medial_axis[j].x, &xw, &yw, &zw) ;
      fprintf(stdout, "0 %.2f %.2f %.2f 0\n", xw,yw,zw);
    }
    medial_axis = compute_normal(medial_axis, length);

    for (j=1; j<length-1; j++) {
      float radius;
      x = medial_axis[j].x;
      y = medial_axis[j].y;
      radius = medial_axis[j].radius;
      dx = medial_axis[j].nx;
      dy = medial_axis[j].ny;
      dx = dx/sqrt(dx*dx+dy*dy) ;
      dy = dy/sqrt(dx*dx+dy*dy) ;
      MRIvoxelToWorld(mri_cc_tal, cc_tal_x, y+radius*dy/2, x+radius*dx/2, &xw, &yw, &zw) ;
      fprintf(stdout, "0 %.2f %.2f %.2f 1\n", xw,yw,zw);
    }

    for (j=1; j<length-1; j++) {
      float radius;
      x = medial_axis[j].x;
      y = medial_axis[j].y;
      radius = medial_axis[j].radius;
      dx = medial_axis[j].nx;
      dy = medial_axis[j].ny;
      dx = dx/sqrt(dx*dx+dy*dy) ;
      dy = dy/sqrt(dx*dx+dy*dy) ;
      MRIvoxelToWorld(mri_cc_tal, cc_tal_x, y-radius*dy/2, x-radius*dx/2, &xw, &yw, &zw) ;
      fprintf(stdout, "0 %.2f %.2f %.2f 2\n", xw,yw,zw);
    }

  }


  medial_axis = compute_normal(medial_axis, length);
  for (i=0; i<length; i++)
    printf("%d  %f  %f  %f  %f\n", i+1, medial_axis[i].x, medial_axis[i].y, medial_axis[i].nx, medial_axis[i].ny);

  /* starting writing out medial axis to the text file */

  if ((fp = fopen(ofname, "w")) == NULL) {
    fprintf(stdout, "cc medial axis measurment: file %s does not exist!", ofname);
  }
  print("writing medial axis results to %s\n",ofname);
  // fprintf(fp, "This file contains the position and radius of 100 equally sampled medial atoms along the medial axis\n");
  MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, medial_axis[0].y, medial_axis[0].x, 0, SAMPLE_TRILINEAR, &val) ;
  fprintf(fp, "0   %f   %f   0   %f   %f   %f   %f\n", medial_axis[0].x, medial_axis[0].y, val, val, val, val);

  dist = 0;
  medial_axis[0].odx = 0;
  for (j=1; j<length; j++) {
    dist+= sqrt((medial_axis[j].x-medial_axis[j-1].x)*(medial_axis[j].x-medial_axis[j-1].x)+(medial_axis[j].y-medial_axis[j-1].y)*(medial_axis[j].y-medial_axis[j-1].y));
    medial_axis[j].odx = dist;
  }

  space = dist/(n_sample+1);
  j=0;
  for (n=1; n<=n_sample; n++) {
    double upmid, downmid;
    float radius;
    dist = n*space;
    for (;medial_axis[j].odx<dist&&j<length;j++) {
      ;
    }

    if (j==length) printf("sample error at %dth atoms\n", n);
    else {
      x1 = medial_axis[j-1].x;
      y1 = medial_axis[j-1].y;
      x2 = medial_axis[j].x;
      y2 = medial_axis[j].y;
      fraction = (dist-medial_axis[j-1].odx)/(medial_axis[j].odx-medial_axis[j-1].odx);
      x = fraction*x2+(1-fraction)*x1;
      y = fraction*y2+(1-fraction)*y1;

      dx = medial_axis[j-1].nx + medial_axis[j].nx;
      dy = medial_axis[j-1].ny + medial_axis[j].ny;
      dx = dx/sqrt(dx*dx+dy*dy) ;
      dy = dy/sqrt(dx*dx+dy*dy) ;
      xx = x;
      yy = y;
      count = 0;
      MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, yy, xx, 0, SAMPLE_NEAREST, &mean_val) ;
      for (i=1; i<=10; i++) {
        xx += dx;
        yy += dy;
        MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, yy, xx, 0, SAMPLE_NEAREST, &val) ;
        if (val) {
          count++;
          mean_val+=val;
        } else break;
      }
      xx = x;
      yy = y;
      for (i=1; i<=10; i++) {
        xx -= dx;
        yy -= dy;
        MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, yy, xx, 0, SAMPLE_NEAREST, &val) ;
        if (val) {
          count++;
          mean_val+=val;
        } else break;
      }
      mean_val /= count;

      MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, y, x, 0, SAMPLE_TRILINEAR, &val) ;
      radius = sample_distance_map(cc_boundary, cc_smoothed, x, y);
      MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, y+radius*dy/2, x+radius*dx/2, 0, SAMPLE_TRILINEAR, &upmid) ;
      MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, y-radius*dy/2, x-radius*dx/2, 0, SAMPLE_TRILINEAR, &downmid) ;
      fprintf(fp, "%d   %f   %f   %f   %f   %f   %f   %f\n", n, x, y, radius, val, mean_val, upmid, downmid);
    }
  }
  MRIsampleVolumeFrameType(mri_intensity, cc_tal_x, medial_axis[length-1].y, medial_axis[length-1].x , 0, SAMPLE_TRILINEAR, &val) ;
  fprintf(fp, "%d   %f   %f   0   %f   %f   %f   %f\n", n, medial_axis[length-1].x, medial_axis[length-1].y, val, val, val, val);

  fclose(fp);
  return (medial_axis) ;
}

static MEDATOM *medial_axis_initialization(MEDATOM *medial_axis, int length) {
  float dx,dy;
  int   i;

  /*draw a line between the two leaf nodes*/
  dx = (medial_axis[length-1].x - medial_axis[0].x)/(length-1);
  dy = (medial_axis[length-1].y - medial_axis[0].y)/(length-1);
  for (i=1; i < length; i++) {
    medial_axis[i].x = medial_axis[i-1].x + dx;
    medial_axis[i].y = medial_axis[i-1].y + dy;
    medial_axis[i].dx = 0;
    medial_axis[i].dy = 0;
    medial_axis[i].odx = 0;
    medial_axis[i].ody = 0;
    medial_axis[i].nx = 0 ;
    medial_axis[i].ny = 0 ;
  }
  /*compute normal*/
  medial_axis = compute_normal(medial_axis, length);
#if 0
  /*print out the medial axis*/
  for (i=0; i<length; i++) {
    fprintf(stdout,"%d th atom: x=%f  y=%f  nx=%f  ny=%f\n", i, medial_axis[i].x, medial_axis[i].y, medial_axis[i].nx, medial_axis[i].ny);
  }
#endif
  return (medial_axis) ;
}

static MEDATOM *medial_axis_integration(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length) {
  int      n=0, i=0, number=0, count=50;
  double   sse=0, old_sse=0, min_sse;
  float    p=1, k=1, iteration_number=500, tol=0.01, ratio=0, move=0;

  sse = compute_energy(medial_axis, cc_boundary, cc_smoothed, length);
  min_sse = sse;
  if (length<50) {
    iteration_number=200;
    count=20;
  }

  do {
    old_sse=sse;
    p = 1/(1 + 0.2*floor(n/100));
    k = 1/(1 + 0.2*floor(n/100));
    tol = 0.005/(1 + 0.1*floor(n/100));

    medial_axis = compute_gradient_term(medial_axis,cc_boundary,cc_smoothed,length,p);
    medial_axis = compute_smooth_term(medial_axis,length,k);
    move = 0;
    for (i=1; i<length-1; i++)
      move += (medial_axis[i].dx *medial_axis[i].dx) + (medial_axis[i].dy *medial_axis[i].dy) ;
    medial_axis = integrate_momentum(medial_axis,cc_boundary,cc_smoothed, length);
    sse = compute_energy(medial_axis, cc_boundary, cc_smoothed, length);
    if ( alert) iteration_number=1000;
    else if (length>50) iteration_number=500;
    else iteration_number=200;
    n++;
    ratio = (old_sse-sse)/sse;
    if (sse<min_sse) {
      min_sse=sse;
      number = 0;
    }
    if (ratio<0) ratio = 1;
    if ( (sse-min_sse)/min_sse < tol ) number++ ;
    else  number = 0 ;
    //fprintf(stdout, "%d iteration: min_sse= %f   sse=%f  old_sse=%f   ratio=%f tol=%f  p=%f move=%f number=%d\n", n, min_sse, sse, old_sse, ratio, tol,  p, move, number);
  } while (n<iteration_number&&move>0.01&&(number<count||ratio>tol||alert));

#if 0  //for debuging
  FILE     *fp;
  int      index;
  fp = fopen("/space/solo/4/users/recon/ALBJ82_recon/callosum/ma.txt", "r");
  for (i=0; i<length; i++) {
    fscanf(fp, "%d %f %f %f %f", &index, &medial_axis[i].x, &medial_axis[i].y, &medial_axis[i].nx, &medial_axis[i].ny);
  }
  fclose(fp);
#endif

#if 1
  /*print out the medial axis*/
  if ( (n%1)==0) {
    fprintf(stdout, "%d iteration: min_sse= %f   sse=%f  old_sse=%f   ratio=%f tol=%f  p=%f move=%f number=%d\n", n, min_sse, sse, old_sse, ratio, tol,  p, move, number);
    for (i=0; i<length; i++) {
      fprintf(stdout,"%f  %f  %f  %f\n", medial_axis[i].x, medial_axis[i].y, medial_axis[i].nx, medial_axis[i].ny);
    }
  }
#endif

  return (medial_axis) ;
}

static MEDATOM *compute_gradient_term(MEDATOM *medial_axis,MRI *cc_boundary,MRI *cc_smoothed,int length, float p) {
  int   i=0;
  float x, y, dx, dy, nx, ny, val_in, val_out;

  for (i=1; i<length-1; i++) {
    x = medial_axis[i].x;
    y = medial_axis[i].y;
    nx = medial_axis[i].nx;
    ny = medial_axis[i].ny;
    val_out = sample_distance_map(cc_boundary, cc_smoothed, x+nx, y+ny);
    val_in = sample_distance_map(cc_boundary, cc_smoothed, x-nx, y-ny);
    dx = (val_out-val_in)*p*nx/2;
    dy = (val_out-val_in)*p*ny/2;
    medial_axis[i].dx += dx;
    medial_axis[i].dy += dy;
  }
  return(medial_axis);
}

static float sample_distance_map(MRI *cc_boundary, MRI *cc_smoothed, float x, float y) {
  int width = cc_boundary->width,height = cc_boundary->height, i, j;
  float  dist, closest=1000;

  if (x<0) x=0;
  if (y<0) y=0;
  if (x>width) x=width;
  if (y>height) y=height;

  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      if (MRIvox(cc_boundary,i,j,0)) {
        dist = sqrt((x-i)*(x-i)+(y-j)*(y-j));
        if (dist<closest)
          closest = dist;
      }

  if (MRIvox(cc_smoothed,nint(x),nint(y),0)==0) closest *= -1;
  return(closest);
}

static MEDATOM *compute_smooth_term(MEDATOM *medial_axis,int length, float k) {
  int  i;
  float dx, dy, nx, ny, val;

  for (i=1; i<length-1; i++) {
    val = 0;
    nx = medial_axis[i].nx;
    ny = medial_axis[i].ny;
    dx =  medial_axis[i+1].x - medial_axis[i].x ;
    dy =  medial_axis[i+1].y - medial_axis[i].y ;
    val = dx*nx + dy*ny;
    dx =  medial_axis[i-1].x - medial_axis[i].x ;
    dy =  medial_axis[i-1].y - medial_axis[i].y ;
    val += dx*nx + dy*ny;
    medial_axis[i].dx += k*val*nx;
    medial_axis[i].dy += k*val*ny;
  }
  return(medial_axis);
}

static MEDATOM *integrate_momentum(MEDATOM *medial_axis, MRI *cc_boundary,MRI *cc_smoothed, int length) {
  int  i, j, k, count;
  float f, scale = 0, slope, x1, x2, y1, y2;
  double sse, min_sse, dist, space;
  MEDATOM *temp;

  temp = (MEDATOM *)calloc(length, sizeof(MEDATOM));
  for (i=0; i<length; i++) {
    temp[i].dx=medial_axis[i].dx;
    temp[i].dy=medial_axis[i].dy;
    temp[i].nx=medial_axis[i].nx;
    temp[i].ny=medial_axis[i].ny;
    temp[i].x=medial_axis[i].x;
    temp[i].y=medial_axis[i].y;
  }

  min_sse = compute_energy(medial_axis, cc_boundary, cc_smoothed, length);
  sse = compute_energy(temp, cc_boundary, cc_smoothed, length);
  for (f=-2; f<=3; f+=0.01) {
    for (i=1; i<length-1; i++) {
      temp[i].x = medial_axis[i].x + f*medial_axis[i].dx;
      temp[i].y = medial_axis[i].y + f*medial_axis[i].dy;
    }
    sse = compute_energy(temp, cc_boundary, cc_smoothed, length);
    if (sse<=min_sse) {
      scale = f;
      min_sse = sse;
    }
  }
  if (fabs(scale)<0.05) scale = 0.05;
  /* update the atoms */
  for (i=1; i<length-1; i++) {
    medial_axis[i].x += scale*medial_axis[i].dx;
    medial_axis[i].y += scale*medial_axis[i].dy;
    medial_axis[i].odx = medial_axis[i].dx;
    medial_axis[i].ody = medial_axis[i].dy;
    medial_axis[i].dx = 0;
    medial_axis[i].dy = 0;
  }

  /* check for intersection -simple version*/
  for (i=1; i<length-1; i++) {
    if ( (fabs(medial_axis[i].x-medial_axis[i-1].x)<0.001)&&(fabs(medial_axis[i].y-medial_axis[i-1].y)<0.001)) {
      medial_axis[i].x = (medial_axis[i+1].x - medial_axis[i-1].x)/2 + medial_axis[i-1].x;
      medial_axis[i].y = (medial_axis[i+1].y - medial_axis[i-1].y)/2 + medial_axis[i-1].y;
    }
  }


#if 0
  /* check for intersection - complicated version */
  for (i=1; i<length-1; i++) {
    if (medial_axis[i].y<medial_axis[i-1].y) {
      x1 = medial_axis[i-1].x;
      y1 = medial_axis[i-1].y;
      for (j=i+1; j<length && medial_axis[j].y<medial_axis[i-1].y; j++) {
        ;
      }
      if (j>=length) j=length-1;
      x2 = medial_axis[j].x;
      y2 = medial_axis[j].y;

      for (k=i ; k<j; k++) {
        slope = (medial_axis[j].y-medial_axis[i-1].y)/(j-i+1) ;
        medial_axis[k].y = (k-i+1)*slope+y1;
        slope = (medial_axis[j].x-medial_axis[i-1].x)/(j-i+1) ;
        medial_axis[k].x = (k-i+1)*slope+x1;
      }
    }
  }

  if (medial_axis[length-1].y<medial_axis[length-2].y) {
    for (j=length-2; j>0&&medial_axis[j].y+0.5<medial_axis[length-1].y; j--) {
      ;
    }
    for (k=j+1; k<length-1; k++) {
      slope = (medial_axis[length-1].y-medial_axis[j].y)/(length-1-j);
      medial_axis[k].y = (k-j)*slope+medial_axis[j].y ;
      slope = (medial_axis[length-1].x-medial_axis[j].x)/(length-1-j);
      medial_axis[k].x = (k-j)*slope+medial_axis[j].x ;
    }
  }
#endif

  /*check for small spacing */
  dist = sqrt((medial_axis[length-1].x-medial_axis[0].x)*(medial_axis[length-1].x-medial_axis[0].x)+(medial_axis[length-1].y-medial_axis[0].y)*(medial_axis[length-1].y-medial_axis[0].y));
  space = dist/length;
  for (i=0; i<length-1; i++) {
    dist = sqrt((medial_axis[i+1].x-medial_axis[i].x)*(medial_axis[i+1].x-medial_axis[i].x)+(medial_axis[i+1].y-medial_axis[i].y)*(medial_axis[i+1].y-medial_axis[i].y));
    if (dist<space/sqrt(2) && i<length-2) {
      count = 1;
      j=i+1;
      do {
        count++;
        j++;
        dist += sqrt((medial_axis[j].x-medial_axis[j-1].x)*(medial_axis[j].x-medial_axis[j-1].x)+(medial_axis[j].y-medial_axis[j-1].y)*(medial_axis[j].y-medial_axis[j-1].y));
      } while (dist<count*space&&j<length-1);
      if (j>i+1) {
        x1 = medial_axis[i].x;
        y1 = medial_axis[i].y;
        x2 = medial_axis[j].x;
        y2 = medial_axis[j].y;
        for (k=i+1;k<j;k++) {
          slope = (y2-y1)/(j-i) ;
          medial_axis[k].y = (k-i)*slope+y1;
          slope = (x2-x1)/(j-i) ;
          medial_axis[k].x = (k-i)*slope+x1;
        }
      }
      i=j;
    }
  }

  medial_axis = compute_normal(medial_axis, length);
  free(temp);
  return(medial_axis);
}

static double compute_energy(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length) {
  double sse=0;
  float  dx, dy, nx, ny;
  int    i, x, y, neg_count=0;

  for (i=1; i<length-1; i++) {
    x = medial_axis[i].x ;
    y = medial_axis[i].y ;
    nx = medial_axis[i].nx ;
    ny = medial_axis[i].ny ;

    if (sample_distance_map(cc_boundary, cc_smoothed, x, y)<0) neg_count++;
    sse += 10 - sample_distance_map(cc_boundary, cc_smoothed, x, y);
    dx = medial_axis[i+1].x + medial_axis[i-1].x - 2*medial_axis[i].x;
    dy = medial_axis[i+1].y + medial_axis[i-1].y - 2*medial_axis[i].y;
    sse += sqrt (dx*dx + dy*dy);
  }
  if (neg_count > length/50) {
    alert = 1;
    //fprintf(stdout, "found %d negative points\n", neg_count);
  } else alert = 0;
  return(sse);
}

static MEDATOM *medial_axis_refine_leaf(MEDATOM *medial_axis, MRI *cc_boundary, MRI *cc_smoothed, int length) {
  int     i, x,y, new_length=0, mark=0, ratio=1, height=cc_boundary->height;
  float   error=0, min_error=0, x1,y1,x2,y2;
  MEDATOM *temp;

  ratio = nint(height/128);
  if (ROT) {
    for (i=0; i<length && medial_axis[i].y-medial_axis[0].y-10*ratio<=0 ; i++ ) {
      ;
    }
    new_length = nint(length/5)>i+1 ? nint(length/5):i+1;

    temp = (MEDATOM *)calloc(new_length, sizeof(MEDATOM));

    /* fix the end leaf, adjust the head leaf*/
    error = compute_error(medial_axis,cc_boundary,cc_smoothed,length, new_length,0);
    min_error = 100;
    temp[new_length-1].x=medial_axis[new_length-1].x;
    temp[new_length-1].y=medial_axis[new_length-1].y;
    x1 = medial_axis[0].x ;
    y1 = medial_axis[0].y ;

    mark = nint(medial_axis[0].y)+1;
    for ( y=mark; y<=4*ratio+mark; y++) {
      for (x=0;x<temp[new_length-1].x;x++) {
        if (MRIvox(cc_boundary,x,y,0)) {
          temp[0].x = x ;
          temp[0].y = y ;
          temp = medial_axis_initialization(temp, new_length) ;
          temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
          for (i=0;i<new_length;i++) {
            medial_axis[i].x = temp[i].x;
            medial_axis[i].y = temp[i].y;
          }
          error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,0);
          printf("leaf node x=%d, y=%d, error=%f\n", x, y, error) ;
          if (error<min_error) {
            min_error = error;
            x1 = x;
            y1 =y;
          }
        }
      }
    }

    temp[0].x = x1 ;
    temp[0].y = y1 ;
    temp = medial_axis_initialization(temp, new_length) ;
    temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
    for (i=0;i<new_length;i++) {
      medial_axis[i].x = temp[i].x;
      medial_axis[i].y = temp[i].y;
    }
    free(temp);



    for (i=length-1; i>0 && medial_axis[length-1].y-medial_axis[i].y-ratio*10<=0 ; i-- ) {
      ;
    }
    new_length = nint(length/5)>(length-i+1) ? nint(length/5):(length-i+1);

    temp = (MEDATOM *)calloc(new_length, sizeof(MEDATOM));
    /* fix the head leaf, adjust the end leaf*/
    error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,1);
    //min_error = error;
    min_error=100;
    temp[0].x=medial_axis[length-new_length].x;
    temp[0].y=medial_axis[length-new_length].y;
    x2 = medial_axis[length-1].x;
    y2 = medial_axis[length-1].y;

    mark = nint(medial_axis[length-1].y)-1; //changed on Jan.6
    for (y=mark;y>=mark-12*ratio;y--) {
      for (x=0;x<temp[0].x;x++) {
        if (MRIvox(cc_boundary,x,y,0)) {
          temp[new_length-1].x = x ;
          temp[new_length-1].y = y ;
          temp = medial_axis_initialization(temp, new_length) ;
          temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
          for (i=0;i<new_length;i++) {
            medial_axis[length-new_length+i].x = temp[i].x;
            medial_axis[length-new_length+i].y = temp[i].y;
          }
          error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,1);
          printf("leaf node x=%d, y=%d, error=%f\n", x, y, error) ;
          if (error<min_error) {
            min_error = error;
            x2 = x;
            y2 =y;
          }
        }
      }
    }

    temp[new_length-1].x = x2 ;
    temp[new_length-1].y = y2 ;
    temp = medial_axis_initialization(temp, new_length) ;
    temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
    for (i=0;i<new_length;i++) {
      medial_axis[length-new_length+i].x = temp[i].x;
      medial_axis[length-new_length+i].y = temp[i].y;
    }
    /* end of re-estimation*/
  } else {
    for (i=0; i<length && medial_axis[i].x-medial_axis[0].x-10<=0 ; i++ ) {
      ;
    }
    new_length = nint(length/5)>i+1 ? nint(length/5):i+1;

    temp = (MEDATOM *)calloc(new_length, sizeof(MEDATOM));

    /* fix the end leaf, adjust the head leaf*/
    error = compute_error(medial_axis,cc_boundary,cc_smoothed,length, new_length,0);
    min_error = 100;
    temp[new_length-1].x=medial_axis[new_length-1].x;
    temp[new_length-1].y=medial_axis[new_length-1].y;
    x1 = medial_axis[0].x ;
    y1 = medial_axis[0].y ;

    mark = nint(medial_axis[0].x);
    for ( x=mark+2; x<=mark+20; x++) {
      for (y=temp[new_length-1].y;y<255;y++) {
        if (MRIvox(cc_boundary,x,y,0)) {
          temp[0].x = x ;
          temp[0].y = y ;
          temp = medial_axis_initialization(temp, new_length) ;
          temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
          for (i=0;i<new_length;i++) {
            medial_axis[i].x = temp[i].x;
            medial_axis[i].y = temp[i].y;
          }
          error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,0);
          printf("leaf node x=%d, y=%d, error=%f\n", x, y, error) ;
          if (error<min_error) {
            min_error = error;
            x1 = x;
            y1 =y;
          }
        }
      }
    }

    temp[0].x = x1 ;
    temp[0].y = y1 ;
    temp = medial_axis_initialization(temp, new_length) ;
    temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
    for (i=0;i<new_length;i++) {
      medial_axis[i].x = temp[i].x;
      medial_axis[i].y = temp[i].y;
    }
    free(temp);

    for (i=length-1; i>0 && medial_axis[length-1].x-medial_axis[i].x-10<=0 ; i-- ) {
      ;
    }
    new_length = nint(length/5)>(length-i+1) ? nint(length/5):(length-i+1);

    temp = (MEDATOM *)calloc(new_length, sizeof(MEDATOM));
    /* fix the head leaf, adjust the end leaf*/
    error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,1);
    min_error = 100;
    temp[0].x=medial_axis[length-new_length].x;
    temp[0].y=medial_axis[length-new_length].y;
    x2 = medial_axis[length-1].x ;
    y2 = medial_axis[length-1].y ;

    mark = nint(medial_axis[length-1].x);
    for (x=mark-2;x>=mark-25;x--) {
      for (y=temp[0].y;y<255;y++) {
        if (MRIvox(cc_boundary,x,y,0)) {
          temp[new_length-1].x = x ;
          temp[new_length-1].y = y ;
          temp = medial_axis_initialization(temp, new_length) ;
          temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
          for (i=0;i<new_length;i++) {
            medial_axis[length-new_length+i].x = temp[i].x;
            medial_axis[length-new_length+i].y = temp[i].y;
          }
          error = compute_error(medial_axis,cc_boundary,cc_smoothed,length,new_length,1);
          printf("leaf node x=%d, y=%d, error=%f\n", x, y, error) ;
          if (error<min_error) {
            min_error = error;
            x2 = x;
            y2 =y;
          }
        }
      }
    }

    temp[new_length-1].x = x2 ;
    temp[new_length-1].y = y2 ;
    temp = medial_axis_initialization(temp, new_length) ;
    temp = medial_axis_integration(temp, cc_boundary, cc_smoothed, new_length) ;
    for (i=0;i<new_length;i++) {
      medial_axis[length-new_length+i].x = temp[i].x;
      medial_axis[length-new_length+i].y = temp[i].y;
    }
    /* end of re-estimation*/
  }

  free(temp);
  return(medial_axis);
}

#if 0
static float compute_error(MEDATOM *medial_axis,MRI *cc_boundary, MRI *cc_smoothed, int length, int new_length, int flag) {
  float error=0, dist=0, closest=0, x, y;
  int   n, node, count=0, width = cc_boundary->width,height = cc_boundary->height, i, j;


  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      if (MRIvox(cc_boundary, i, j, 0)) {
        closest=1000;
        node=1000;
        for (n=0;n<length-1;n++) {
          dist = sqrt((medial_axis[n].x-i)*(medial_axis[n].x-i)+(medial_axis[n].y-j)*(medial_axis[n].y-j) );
          if ( dist<closest ) {
            closest = dist;
            node = n;
          } else if ( dist == closest ) {
            fprintf(stdout, "point %d %d has more than one cloesest point!\n", i,j);
          }
        }
        if ( (!flag&&node<new_length)||(flag&&node<length&&node>length-new_length-1) ) {
          dist = sample_distance_map(cc_boundary, cc_smoothed, medial_axis[node].x, medial_axis[node].y) ;
          error += (dist-closest)*(dist-closest);
          count++;
        }
      }
  // if(count>10)
  error = error/count;
  fprintf(stdout, "count number is %d\n", count);
  //else
  // error = 1000;
  MRIfree(&cc_new) ;
  return(error);
}
#else
static float compute_error(MEDATOM *medial_axis,MRI *cc_boundary, MRI *cc_smoothed, int length, int new_length, int flag) {
  float error=0, dist=0, closest=0, x, y, radius;
  int   node, count=0, width = cc_boundary->width,height = cc_boundary->height;
  int   i, j, m, n;
  MRI   *cc_new, *cc_newb;

  /* First build the new boundary */
  cc_new = MRIcopy(cc_boundary, NULL) ;
  MRIvalueFill(cc_new, 0) ;

  for (n=0;n<length;n++) {
    x = medial_axis[n].x ;
    y = medial_axis[n].y ;
    radius = sample_distance_map(cc_boundary, cc_smoothed, x, y) ;
    medial_axis[n].radius = radius ;
    for (i=0; i<width; i++)
      for (j=0; j<height; j++) {
        dist = sqrt((x-i)*(x-i)+(y-j)*(y-j)) ;
        if (dist<=radius+0.2)
          MRIvox(cc_new,i,j,0) += 1;
      }
  }
  cc_newb = find_cc_boundary(cc_new, cc_newb) ;
  MRIwrite(cc_new, "/space/solo/4/users/recon/DONF55b_recon/callosum/cc_new.mgh") ;

  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      if (MRIvox(cc_newb, i, j, 0)) {
        closest=1000;
        for (n=0;n<length;n++) {
          x = medial_axis[n].x ;
          y = medial_axis[n].y ;
          dist = sqrt((x-i)*(x-i)+(y-j)*(y-j)) ;
          dist = fabs (dist-medial_axis[n].radius);
          if ( dist < closest) {
            MRIvox(cc_newb,i,j,0)=n+1;
            closest = dist;
          }
        }
      }

  MRIwrite(cc_newb, "/space/solo/4/users/recon/DONF55b_recon/callosum/cc_newb.mgh") ;

  for (i=0; i<width; i++)
    for (j=0; j<height; j++)
      if (MRIvox(cc_boundary, i, j, 0)) {
        closest=1000;
        for (m=0; m<width; m++)
          for (n=0; n<height; n++)
            if (MRIvox(cc_newb, m, n, 0)) {
              dist = sqrt((m-i)*(m-i)+(n-j)*(n-j));
              if (dist<closest) {
                closest=dist;
                node=MRIvox(cc_newb, m, n, 0)-1;
              }
            }
        if ( (!flag&&node<new_length)||(flag&&node<length&&node>length-new_length-1) ) {
          error+=closest;
          count++;
        }
      }

  error = error/count;
  fprintf(stdout, "count number is %d\n", count);
  MRIfree(&cc_newb) ;
  MRIfree(&cc_new) ;
  return(error);
}

#endif

static MRI *sample_medial_axis(MEDATOM *medial_axis, MRI *cc_medial_axis, int length) {
  MRI   *cc_slice;
  int   j, n, yi_low, yi_high;
  float slope=0, x1, y1, x2, y2, x, y;

  cc_slice = MRIextractPlane(cc_medial_axis, NULL, MRI_SAGITTAL, cc_tal_x);
  yi_low = nint(medial_axis[0].y);
  yi_high = nint(medial_axis[length-1].y);
  MRIvox(cc_slice,nint(medial_axis[0].x),yi_low,0) = CC_VAL;
  MRIvox(cc_slice,nint(medial_axis[length-1].x),yi_high,0) = CC_VAL;

#if 0
  n = 0;
  for ( j=yi_low+1; j<yi_high; j++) {
    y = j ;

    for ( ; medial_axis[n].y<y && n<length ; n++) {}
    if (n==length) {
      fprintf(stdout, "ERROR in interpolation: position %d can not the sampled", j);
      n = 0;
    }

    y2 = medial_axis[n].y;
    x2 = medial_axis[n].x;

    if ( fabs(y-y2)<0.01)
      MRIvox(cc_slice,nint(x2),j,0) = CC_VAL;
    else {
      y1 = medial_axis[n-1].y;
      x1 = medial_axis[n-1].x;
      slope = (x2-x1)/(y2-y1);
      x = slope*(y-y1)+x1;
      MRIvox(cc_slice,nint(x),j,0) = CC_VAL;
    }
    if (n) n--;
  }
#else

  for (j=0;j<length-1;j++) {
    x1 = medial_axis[j].x;
    y1 = medial_axis[j].y;
    MRIvox(cc_slice,nint(x1),nint(y1),0) = CC_VAL;
    printf("voxel (%d %d) marked \n", nint(x1),nint(y1) );
    do {
      j++;
      x2 = medial_axis[j].x;
      y2 = medial_axis[j].y;
    } while ( j<length-1 && fabs(nint(x2)-nint(x1))<1 && fabs(nint(y2)-nint(y1) )<1 );

    if ( fabs(nint(y2)-nint(y1))>1 ) {
      if (nint(y2)<nint(y1)) {
        x = x1;
        y = y1;
        x1 = x2;
        y1 = y2;
        x2 = x;
        y2 = y;
      }
      slope = (x2-x1)/(y2-y1);
      n = nint(y1)+1;
      for (; n<nint(y2); n++) {
        x = slope*(n-y1)+x1;
        MRIvox(cc_slice,nint(x),n,0) = CC_VAL;
        printf("voxel (%d %d) y intepolated \n", nint(x),n );
      }
    } else if ( fabs(nint(x2)-nint(x1))>1 ) {
      if (nint(x2)<nint(x1)) {
        x = x1;
        y = y1;
        x1 = x2;
        y1 = y2;
        x2 = x;
        y2 = y;
      }
      slope = (y2-y1)/(x2-x1);
      n = nint(x1)+1;
      for (; n<nint(x2); n++) {
        y = slope*(n-x1)+y1;
        MRIvox(cc_slice,n,nint(y),0) = CC_VAL;
        printf("voxel (%d %d) x intepolated \n", n, nint(y) );
      }
    }


    if (j==length-1) break;
    else j--;
  }

#endif
  //MRIwrite(cc_slice, "/space/yoda/5/users/salat/callosum/EVAC79_recon_dti_cc/cc_medial_axis.mgh");

  MRIfillPlane(cc_slice, cc_medial_axis, MRI_SAGITTAL, cc_tal_x, CC_VAL);
  MRIfree(&cc_slice);

  return (cc_medial_axis) ;
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

  if (!stricmp(option, "wm")) {
    WM = 1 ;
    fprintf(stdout,"read in white matter volume %s\n", wmvolume);
  } else if (!stricmp(option, "tal")) {
    TAL = 1 ;
    fprintf(stdout,"read in talairach transformed volume\n");
  } else if (!stricmp(option, "seed")) {
    cc_tal_x = atof(argv[2]) ;
    cc_tal_y = atof(argv[3]) ;
    cc_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    ROT = 0;
    fprintf(stderr, "T1 volume: cc seed at (%f %f %f)\n",
            cc_tal_x, cc_tal_y, cc_tal_z) ;
  } else if (!stricmp(option, "NC")) {
    NO_CORRECTION = 1;
    fprintf(stdout, "No CC slice correction");
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      fprintf(stdout,
              "usage: %s <input volumes> <output volume>\n",
              Progname) ;
      exit(1) ;
      break ;
    case 'T':
      dxi = atoi(argv[2]);
      fprintf(stdout,"change thickness to %d mm\n", 2*dxi+1);
      nargs = 1;
      break;
    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}






