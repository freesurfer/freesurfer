#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "mri.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "mrimorph.h"
#include "timer.h"

static char vcid[] = "$Id: mri_fill.c,v 1.37 1999/11/05 00:17:39 fischl Exp $";

/*-------------------------------------------------------------------
                                CONSTANTS
-------------------------------------------------------------------*/

#define DEBUG_POINT(x,y,z)  (((x==140) && (y==74)) &&((z)==174))

#define PONS_LOG_FILE                   "pons.log"
#define CC_LOG_FILE                     "cc.log"

/* min # of neighbors which must be on to retain a point */
#define DEFAULT_NEIGHBOR_THRESHOLD      8


/* distance to search in each direction for a valid wm seed point */
#define SEED_SEARCH_SIZE                9

/* size of various orthogonal slices */
#define SLICE_SIZE                      128

/* anything below this is not white matter */
#define WM_MIN_VAL                       2 

/* the min # of neighbors on for a point to be allowed to be a seed point */
#define MIN_NEIGHBORS                    4

/* Talairach seed points - only used if heuristics fail */
#define CORPUS_CALLOSUM_TAL_X            0.0
#define CORPUS_CALLOSUM_TAL_Y            0.0
#define CORPUS_CALLOSUM_TAL_Z            27.0

#define PONS_TAL_X                      -2.0
#define PONS_TAL_Y                      -15.0  /* was -22.0 */
#define PONS_TAL_Z                      -17.0

#define WHITE_MATTER_RH_TAL_X            29.0
#define WHITE_MATTER_RH_TAL_Y            -12.0 ;
#define WHITE_MATTER_RH_TAL_Z            28.0 ;


#define WHITE_MATTER_LH_TAL_X            -29.0
#define WHITE_MATTER_LH_TAL_Y            -12.0 ;
#define WHITE_MATTER_LH_TAL_Z            28.0 ;


/* # if iterations for filling */
#define MAX_ITERATIONS  1000  /* was 10 */

/* min # of holes filled before quitting */
#define MIN_FILLED      0     /* was 100 */

#define FILLED_PRIORS_FNAME  "talairach_hemi_filled.mgh"

/*-------------------------------------------------------------------
                                GLOBAL DATA
-------------------------------------------------------------------*/

static FILE *log_fp = NULL ;
static int lh_fill_val = MRI_LEFT_HEMISPHERE ;
static int rh_fill_val = MRI_RIGHT_HEMISPHERE ;

static int ylim0,ylim1,xlim0,xlim1;
static int fill_holes_flag = TRUE;

/* corpus callosum seed point in Talairach coords */

static Real cc_tal_x = 0.0 ;
static Real cc_tal_y = 0.0 ;
static Real cc_tal_z = 27.0 ;

static Real pons_tal_x = -2.0 ;
static Real pons_tal_y = -15.0 /* -22.0 */ ;
static Real pons_tal_z = -17.0 ;


#if 0
/* test coords - not very close */
static Real cc_tal_x = -4.0 ;
static Real cc_tal_y = -32.0 ;
static Real cc_tal_z = 27.0 ;

static Real pons_tal_x = -10.0 ;
static Real pons_tal_y = -36.0 ;
static Real pons_tal_z = -20.0 ;

#endif

static int cc_seed_set = 0 ;
static int pons_seed_set = 0 ;

char *Progname ;

static int min_filled = 0 ;

static int neighbor_threshold = DEFAULT_NEIGHBOR_THRESHOLD ;

static MRI *mri_fill, *mri_im ;

static int logging = 0 ;

static int fill_val = 0 ;   /* only non-zero for generating images of planes */
static int cc_mask = 1 ;

static char *sname ;

static char *atlas_name = NULL ;
#if 0
static float blur_sigma = 0.25f ;
#endif

/*-------------------------------------------------------------------
                             STATIC PROTOTYPES
-------------------------------------------------------------------*/

static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;
int main(int argc, char *argv[]) ;

static int fill_holes(MRI *mri_fill) ;
static int fill_brain(MRI *mri_fill, MRI *mri_im, int threshold) ;
static MRI *find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,
                              int orientation, int *pxv, int *pvy, int *pzv,
                               int seed_set) ;
static int find_slice_center(MRI *mri,  int *pxo, int *pyo) ;
static int find_corpus_callosum(MRI *mri, Real *ccx, Real *ccy, Real *ccz) ;
static int find_pons(MRI *mri, Real *p_ponsx, Real *p_ponsy, Real *p_ponsz,
                     int x_cc, int y_cc, int z_cc, int which) ;
static int find_cc_slice(MRI *mri, Real *pccx, Real *pccy, Real *pccz) ;
static int neighbors_on(MRI *mri, int x0, int y0, int z0) ;
static int MRIfillVolume(MRI *mri_fill, MRI *mri_im, int x_seed, int y_seed, 
                         int z_seed, int fill_val) ;
static int   mriFindBoundingBox(MRI *mri_im) ;
MRI *MRIcombineHemispheres(MRI *mri_lh_fill, MRI *mri_rh_fill, MRI *mri_dst,
                      Real cc_tal_x) ;
static MRI *mriReadBinaryProbabilities(char *atlas_name, char *suffix, 
                                       M3D *m3d, char *subject_name, 
                                       MRI *mri_dst) ;
static MRI *mriReadConditionalProbabilities(MRI *mri_T1, char *atlas_name, 
                                            char *suffix, int offset, 
                                            M3D *m3d, MRI *mri_dst) ;
static MRI *MRIfillVentricle(MRI *mri_src, MRI *mri_prob, MRI *mri_T1,
                             MRI *mri_dst, float thresh, int out_label,
                             int xmin, int xmax) ;

/*-------------------------------------------------------------------
                                FUNCTIONS
-------------------------------------------------------------------*/

int
main(int argc, char *argv[])
{
  int     x, y, z, xd, yd, zd, xnew, ynew, znew, msec ;
  int     nargs, wm_rh_x, wm_rh_y, wm_rh_z, wm_lh_x, wm_lh_y, wm_lh_z ;
  char    input_fname[STRLEN],out_fname[STRLEN], fname[STRLEN] ;
  Real    xr, yr, zr, dist, min_dist, lh_tal_x, rh_tal_x ;
  MRI     *mri_cc, *mri_pons, *mri_lh_fill, *mri_rh_fill, *mri_lh_im, 
          *mri_rh_im /*, *mri_blur*/ ;
  int     x_pons, y_pons, z_pons, x_cc, y_cc, z_cc, xi, yi, zi ;
  MORPH_3D  *m3d ;
  struct timeb  then ;
  MRI    *mri_tmp ;

  TimerStart(&then) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  Progname = argv[0] ;

  for ( ; argc > 1 && (*argv[1] == '-') ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  
  if (argc < 3)
    print_help() ;   /* will exit */

  strcpy(input_fname, argv[1]) ;
  strcpy(out_fname, argv[2]) ;

  if (!Gdiag)
    fprintf(stderr, "reading input volume...") ;

  mri_im = MRIread(input_fname) ;
  if (!mri_im)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s", Progname, input_fname) ;

  if (atlas_name && 0)
  {
    MRI *mri_p_wm ;
    char subjects_dir[STRLEN], mri_dir[STRLEN], *cp ;

    cp = getenv("MRI_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: could not read MRI_DIR from environment",
                Progname) ;
    strcpy(mri_dir, cp) ;

    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not read SUBJECTS_DIR from environment",Progname) ;
    strcpy(subjects_dir, cp) ;


    sprintf(fname, "%s/%s/mri/transforms/%s.m3d", subjects_dir,sname,
            atlas_name) ;
    m3d = MRI3DreadSmall(fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not open transform file %s\n",
                Progname, fname) ;
    fprintf(stderr, "removing extraneous white matter using atlas %s...\n",
            atlas_name) ;
    mri_p_wm = mriReadBinaryProbabilities(atlas_name, "wm", m3d, sname, NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_p_wm, "p_wm.mgh") ;

    MRIprobabilityThresholdNeighborhoodOff(mri_im, mri_p_wm, mri_im, 0.0, 3);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_im, "wm_thresh.mgh") ;

#if 0
    fprintf(stderr, "filling high probability white matter...\n") ;
    mri_blur = MRIgaussian1d(blur_sigma, 0) ;
    if (!mri_blur)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not allocate blurring kernel with sigma=%2.3f",
                Progname, blur_sigma) ;
    mri_tmp = MRIconvolveGaussian(mri_p_wm, NULL, mri_blur) ;
    MRIfree(&mri_p_wm) ; mri_p_wm = mri_tmp ;
    MRIprobabilityThreshold(mri_im, mri_p_wm, mri_im, 98.0, 255) ;
#endif

    MRIfree(&mri_p_wm) ;
    MRI3DmorphFree(&m3d) ;
  }

  if (!Gdiag)
    fprintf(stderr, "done.\nsearching for cutting planes...") ;
  if (!cc_seed_set)
  {
    if (find_corpus_callosum(mri_im,&cc_tal_x,&cc_tal_y,&cc_tal_z) 
        != NO_ERROR)
      mri_cc = NULL ;
    else
      mri_cc =
        find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z,
                           MRI_SAGITTAL, &x_cc, &y_cc, &z_cc, cc_seed_set) ;
  }
  else
    mri_cc =
      find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z,
                         MRI_SAGITTAL, &x_cc, &y_cc, &z_cc, cc_seed_set) ;
  if (!mri_cc)  /* heuristic failed - use Talairach coordinates */
  {
    cc_tal_x = CORPUS_CALLOSUM_TAL_X ; 
    cc_tal_y = CORPUS_CALLOSUM_TAL_Y; 
    cc_tal_z = CORPUS_CALLOSUM_TAL_Z;
    mri_cc = 
      find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z, MRI_SAGITTAL,
                         &x_cc, &y_cc, &z_cc, cc_seed_set) ;
    if (!mri_cc)
      ErrorExit(ERROR_BADPARM, "%s: could not find corpus callosum", Progname);
  }
  
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON)
    MRIwrite(mri_cc, "cc.mgh") ;
  if (log_fp)
    fprintf(log_fp, "CC:   %d, %d, %d (TAL: %2.1f, %2.1f, %2.1f)\n",
            x_cc, y_cc, z_cc, cc_tal_x, cc_tal_y, cc_tal_z) ;
  if (!pons_seed_set)
  {
    MRI  *mri_mask ;

    if (cc_mask)
    {
      fprintf(stderr, 
              "masking possible pons locations using cc cutting plane\n") ;
      mri_tmp = MRIcopy(mri_im, NULL) ;
      mri_mask = MRIcopy(mri_cc, NULL) ;
      MRIdilate(mri_mask, mri_mask) ; MRIdilate(mri_mask, mri_mask) ;
      MRIreplaceValues(mri_mask, mri_mask, 0, 255) ;
      MRIreplaceValues(mri_mask, mri_mask, 1, 0) ;
      for (x = x_cc-64 ; x < x_cc+64 ; x++)
        MRIeraseTalairachPlaneNew(mri_tmp, mri_mask, MRI_SAGITTAL, x, y_cc, 
                                  z_cc, 2*SLICE_SIZE, 0);
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRIwrite(mri_tmp, "erased.mgh") ;
        MRIwrite(mri_mask, "mask.mgh") ;
      }
    }
    else
      mri_tmp = mri_im ;
    find_pons(mri_tmp, &pons_tal_x, &pons_tal_y, &pons_tal_z,x_cc,y_cc,z_cc,0);
    if (mri_tmp != mri_im)
    {
      MRIfree(&mri_tmp) ;
      MRIfree(&mri_mask) ;
    }
  }
  mri_pons = 
    find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                       MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons,
                       pons_seed_set);

  if (!pons_seed_set && !mri_pons)  /* first attempt failed - try different */
  {
    MRI  *mri_mask ;

    fprintf(stderr, 
            "initial attempt at finding brainstem failed - initiating backup "
            "plan A....\n") ;
    if (cc_mask)
    {
      mri_tmp = MRIcopy(mri_im, NULL) ;
      mri_mask = MRIcopy(mri_cc, NULL) ;
      MRIdilate(mri_mask, mri_mask) ; MRIdilate(mri_mask, mri_mask) ;
      MRIreplaceValues(mri_mask, mri_mask, 0, 255) ;
      MRIreplaceValues(mri_mask, mri_mask, 1, 0) ;
      for (x = x_cc-64 ; x < x_cc+64 ; x++)
        MRIeraseTalairachPlaneNew(mri_tmp, mri_mask, MRI_SAGITTAL, x, y_cc, 
                                  z_cc, 2*SLICE_SIZE, 0);
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        MRIwrite(mri_tmp, "erased.mgh") ;
        MRIwrite(mri_mask, "mask.mgh") ;
      }
    }
    else
      mri_tmp = mri_im ;
    find_pons(mri_tmp, &pons_tal_x, &pons_tal_y, &pons_tal_z,x_cc,y_cc,z_cc,1);
    if (mri_tmp != mri_im)
    {
      MRIfree(&mri_tmp) ;
      MRIfree(&mri_mask) ;
    }
    mri_pons = 
      find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                         MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons,
                         pons_seed_set);
  }


  if (!mri_pons) /* heuristic failed to find the pons - try Talairach coords */
  {
    pons_tal_x = PONS_TAL_X ; pons_tal_y = PONS_TAL_Y; pons_tal_z = PONS_TAL_Z;
    mri_pons = 
      find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                         MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons,
                         pons_seed_set);
    if (!mri_pons)
      ErrorExit(ERROR_BADPARM, "%s: could not find pons", Progname);
  }

  if (log_fp)
  {
    fprintf(log_fp, "PONS: %d, %d, %d (TAL: %2.1f, %2.1f, %2.1f)\n",
            x_pons, y_pons, z_pons, pons_tal_x, pons_tal_y, pons_tal_z) ;
    fclose(log_fp) ;
  }

  if (atlas_name && 0)
  {
    MRI *mri_p_filled ;
    char subjects_dir[STRLEN], mri_dir[STRLEN], *cp ;

    cp = getenv("MRI_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: could not read MRI_DIR from environment",
                Progname) ;
    strcpy(mri_dir, cp) ;

    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not read SUBJECTS_DIR from environment",Progname) ;
    strcpy(subjects_dir, cp) ;


    sprintf(fname, "%s/%s/mri/transforms/%s.m3d", subjects_dir,sname,
            atlas_name) ;
    m3d = MRI3DreadSmall(fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not open transform file %s\n",
                Progname, fname) ;
    fprintf(stderr, 
            "removing extraneous filled white matter using atlas %s...\n",
            atlas_name) ;
    mri_p_filled = 
      mriReadBinaryProbabilities(atlas_name, "filled", m3d, sname, NULL);
    MRIwrite(mri_p_filled, "p_filled.mgh") ;

    MRIprobabilityThresholdNeighborhoodOff(mri_im, mri_p_filled, mri_im,0.0,3);
#if 1
    fprintf(stderr, "filling high-probability regions...\n") ;
    MRIprobabilityThresholdNeighborhoodOn(mri_im, mri_p_filled, mri_im, 
                                          99.9, 5, 255);
#endif
    MRIwrite(mri_im, "filled_thresh.mgh") ;

    MRIfree(&mri_p_filled) ;
    MRI3DmorphFree(&m3d) ;
  }

  MRIeraseTalairachPlane(mri_im, mri_cc, MRI_SAGITTAL, x_cc, y_cc, z_cc, 
                         SLICE_SIZE, fill_val);
  MRIeraseTalairachPlane(mri_im, mri_pons, MRI_HORIZONTAL, 
                         x_pons, y_pons, z_pons, SLICE_SIZE, fill_val) ;
  if (fill_val)
  {
    fprintf(stderr,"writing out image with cutting planes to 'planes.mgh'.\n");
    MRIwrite(mri_im, "planes.mgh") ;
    fprintf(stderr, "done.\n") ;
    exit(0) ;
  }

  MRIfree(&mri_cc) ;
  MRIfree(&mri_pons) ;
  if (!Gdiag)
    fprintf(stderr, "done.\n") ;

  /* find white matter seed point for the left hemisphere */
  MRItalairachToVoxel(mri_im, cc_tal_x+SEED_SEARCH_SIZE,
                      cc_tal_y,cc_tal_z,&xr,&yr,&zr);

  wm_rh_x = nint(xr) ; wm_rh_y = nint(yr) ; wm_rh_z = nint(zr) ;
  if ((MRIvox(mri_im, wm_rh_x, wm_rh_y, wm_rh_z) <= WM_MIN_VAL) ||
      (neighbors_on(mri_im, wm_rh_x, wm_rh_y, wm_rh_z) < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for lh wm seed...") ;
    for (z = wm_rh_z-SEED_SEARCH_SIZE ; z <= wm_rh_z+SEED_SEARCH_SIZE ; z++)
    {
      zi = mri_im->zi[z] ;
      for (y = wm_rh_y-SEED_SEARCH_SIZE ; y <= wm_rh_y+SEED_SEARCH_SIZE ; y++)
      {
        yi = mri_im->yi[y] ;
        for (x = wm_rh_x-SEED_SEARCH_SIZE ;x <= wm_rh_x+SEED_SEARCH_SIZE ; x++)
        {
          xi = mri_im->xi[x] ;
          if ((MRIvox(mri_im, xi, yi, zi) >= WM_MIN_VAL) &&
              neighbors_on(mri_im, xi, yi, zi) >= MIN_NEIGHBORS)
          {
            xd = (xi - wm_rh_x) ; yd = (yi - wm_rh_y) ; zd = (zi - wm_rh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ; ynew = yi ; znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    wm_rh_x = xnew ; wm_rh_y = ynew ; wm_rh_z = znew ; 
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "found at (%d, %d, %d)\n", xnew, ynew, znew) ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "lh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_rh_x, wm_rh_y, wm_rh_z, 
            neighbors_on(mri_im, wm_rh_x, wm_rh_y, wm_rh_z)) ;

  /* find white matter seed point for the right hemisphere */
#if 0
  MRItalairachToVoxel(mri_im, wm_lh_tal_x,wm_lh_tal_y,wm_lh_tal_z,&xr,&yr,&zr);
#else
  MRItalairachToVoxel(mri_im, cc_tal_x-SEED_SEARCH_SIZE, 
                      cc_tal_y, cc_tal_z, &xr, &yr, &zr);
#endif
  wm_lh_x = nint(xr) ; wm_lh_y = nint(yr) ; wm_lh_z = nint(zr) ;
  if ((MRIvox(mri_im, wm_lh_x, wm_lh_y, wm_lh_z) <= WM_MIN_VAL) ||
      (neighbors_on(mri_im, wm_lh_x, wm_lh_y, wm_lh_z) < MIN_NEIGHBORS))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for rh wm seed...") ;
    for (z = wm_lh_z-SEED_SEARCH_SIZE ; z <= wm_lh_z+SEED_SEARCH_SIZE ; z++)
    {
      zi = mri_im->zi[z] ;
      for (y = wm_lh_y-SEED_SEARCH_SIZE ; y <= wm_lh_y+SEED_SEARCH_SIZE ; y++)
      {
        yi = mri_im->yi[y] ;
        for (x = wm_lh_x-SEED_SEARCH_SIZE ;x <= wm_lh_x+SEED_SEARCH_SIZE ; x++)
        {
          xi = mri_im->xi[x] ;
          if ((MRIvox(mri_im, xi, yi, zi) >= WM_MIN_VAL) &&
              (neighbors_on(mri_im, xi, yi, zi) >= MIN_NEIGHBORS))
          {
            xd = (xi - wm_lh_x) ; yd = (yi - wm_lh_y) ; zd = (zi - wm_lh_z) ;
            dist = xd*xd + yd*yd + zd*zd ;
            if (dist < min_dist)
            {
              xnew = xi ; ynew = yi ; znew = zi ;
              min_dist = dist ;
            }
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "found at (%d, %d, %d)\n", xnew, ynew, znew) ;
    wm_lh_x = xnew ; wm_lh_y = ynew ; wm_lh_z = znew ; 

  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "rh seed point at (%d, %d, %d): %d neighbors on.\n",
            wm_lh_x, wm_lh_y, wm_lh_z, 
            neighbors_on(mri_im, wm_lh_x, wm_lh_y, wm_lh_z)) ;

#if 0
  /* initialize the fill with the detected seed points */
  MRIvox(mri_fill, wm_lh_x, wm_lh_y, wm_lh_z) = lh_fill_val ;
  MRIvox(mri_fill, wm_rh_x, wm_rh_y, wm_rh_z) = rh_fill_val ;
#endif

  mri_lh_fill = MRIclone(mri_im, NULL) ;
  mri_rh_fill = MRIclone(mri_im, NULL) ;
  mri_lh_im = MRIcopy(mri_im, NULL) ;
  mri_rh_im = MRIcopy(mri_im, NULL) ;
  fprintf(stderr, "filling left hemisphere...\n") ;
  MRIfillVolume(mri_lh_fill, mri_lh_im, wm_lh_x, wm_lh_y, wm_lh_z,lh_fill_val);
  fprintf(stderr, "filling right hemisphere...\n") ;
  MRIfillVolume(mri_rh_fill, mri_rh_im, wm_rh_x, wm_rh_y, wm_rh_z,rh_fill_val);
  MRIfree(&mri_lh_im) ; MRIfree(&mri_rh_im) ; MRIfree(&mri_im) ;
  fprintf(stderr, "combining hemispheres...\n") ;
  MRIvoxelToTalairach(mri_lh_fill, (Real)wm_lh_x, (Real)wm_lh_y, (Real)wm_lh_z,
                      &xr, &yr, &zr) ;
  lh_tal_x = xr ;
  MRIvoxelToTalairach(mri_rh_fill, (Real)wm_rh_x, (Real)wm_rh_y, (Real)wm_rh_z,
                      &xr, &yr, &zr) ;
  rh_tal_x = xr ;
  mri_fill = MRIcombineHemispheres(mri_lh_fill, mri_rh_fill, NULL, cc_tal_x) ;
  MRIfree(&mri_lh_fill) ; MRIfree(&mri_rh_fill) ;

  if (atlas_name)
  {
    MRI *mri_p_ventricle, *mri_fg, *mri_bg, *mri_T1, *mri_ventricle ;
    char subjects_dir[STRLEN], mri_dir[STRLEN], *cp ;

    cp = getenv("MRI_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: could not read MRI_DIR from environment",
                Progname) ;
    strcpy(mri_dir, cp) ;

    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not read SUBJECTS_DIR from environment",Progname) ;
    strcpy(subjects_dir, cp) ;


    sprintf(fname, "%s/%s/mri/transforms/%s.m3d", subjects_dir,sname,
            atlas_name) ;
    m3d = MRI3DreadSmall(fname) ;
    if (!m3d)
      ErrorExit(ERROR_NOFILE, 
                "%s: could not open transform file %s\n",
                Progname, fname) ;
    fprintf(stderr, "filling ventricles using atlas %s...\n", atlas_name) ;
    mri_p_ventricle = 
      mriReadBinaryProbabilities(atlas_name, "left_ventricle",m3d,sname,NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_p_ventricle, "p_left_ventricle.mgh") ;

    sprintf(fname, "%s/%s/mri/T1", subjects_dir,sname) ;
    mri_T1 = MRIread(fname) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume %s\n", 
                Progname, fname) ;
#if 0
#define TMP_FILL_VAL   5
    MRIprobabilityThreshold(mri_fill, mri_p_ventricle, mri_fill, 50.0, 
                            TMP_FILL_VAL) ;
    MRIdilateLabel(mri_fill, mri_fill, TMP_FILL_VAL, 2) ;
    MRIreplaceValues(mri_fill, mri_fill, TMP_FILL_VAL, lh_fill_val) ;
#else
    MRItalairachToVoxel(mri_fill,cc_tal_x-10, cc_tal_y, cc_tal_z, &xr,&yr,&zr);
    mri_ventricle = 
      MRIfillVentricle(mri_fill, mri_p_ventricle, mri_T1, NULL, 90.0f, 
                     lh_fill_val, nint(xr), mri_fill->width-1) ;
    MRIwrite(mri_ventricle, "left_ventricle.mgh") ;
    MRIdilate(mri_ventricle, mri_ventricle) ;
    MRIdilate(mri_ventricle, mri_ventricle) ;
    MRIunion(mri_fill, mri_ventricle, mri_fill) ;
    MRIfree(&mri_ventricle) ;
#endif
                            
    MRIfree(&mri_p_ventricle) ;
    
    mri_p_ventricle = 
      mriReadBinaryProbabilities(atlas_name,"right_ventricle",m3d,sname,NULL);
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_p_ventricle, "p_right_ventricle.mgh") ;
#if 0
    MRIprobabilityThreshold(mri_fill, mri_p_ventricle, mri_fill, 50.0, 
                            TMP_FILL_VAL);
    MRIdilateLabel(mri_fill, mri_fill, TMP_FILL_VAL, 2) ;
    MRIreplaceValues(mri_fill, mri_fill, TMP_FILL_VAL, rh_fill_val) ;
#else
    MRItalairachToVoxel(mri_fill,cc_tal_x+10,cc_tal_y, cc_tal_z, &xr,&yr,&zr);
    mri_ventricle = 
      MRIfillVentricle(mri_fill, mri_p_ventricle, mri_T1, NULL, 90.0f, 
                       rh_fill_val, 0, nint(xr)) ;
    MRIwrite(mri_ventricle, "right_ventricle.mgh") ;
    MRIdilate(mri_ventricle, mri_ventricle) ;
    MRIdilate(mri_ventricle, mri_ventricle) ;
    MRIunion(mri_fill, mri_ventricle, mri_fill) ;
#endif

    MRIfree(&mri_T1) ; MRIfree(&mri_ventricle) ; MRIfree(&mri_p_ventricle) ;
    MRI3DmorphFree(&m3d) ;

    mri_fg = MRIclone(mri_fill, NULL) ;
    MRIvox(mri_fg, wm_rh_x, wm_rh_y, wm_rh_z) = rh_fill_val ;
    MRIvox(mri_fg, wm_lh_x, wm_lh_y, wm_lh_z) = lh_fill_val ;
    MRIgrowLabel(mri_fill, mri_fg, lh_fill_val, lh_fill_val) ;
    MRIgrowLabel(mri_fill, mri_fg, rh_fill_val, rh_fill_val) ;
    mri_bg = MRIclone(mri_fill, NULL) ;
    MRIvox(mri_bg, 0, 0, 0) = 1 ;
    MRIgrowLabel(mri_fill, mri_bg, 0, 1) ;
    MRIturnOnFG(mri_fill, mri_fg, mri_bg) ;
    MRIturnOffBG(mri_fill, mri_bg) ;
  }
  
  fprintf(stderr, "writing output to %s...\n", out_fname) ;
  MRIwrite(mri_fill, out_fname) ;
  msec = TimerStop(&then) ;
  fprintf(stderr,"filling took %2.1f minutes\n", (float)msec/(60*1000.0f));
  exit(0) ;
  return(0) ;
}
static int 
fill_brain(MRI *mri_fill, MRI *mri_im, int threshold)
{
  int dir = -1, nfilled = 10000, ntotal = 0,iter = 0;
  int im0,im1,j0,j1,i0,i1,imnr,i,j;
  int v1,v2,v3,vmax ;

  mriFindBoundingBox(mri_im) ;
  while (nfilled>min_filled && iter<MAX_ITERATIONS)
  {
    iter++;
    nfilled = 0;
    dir = -dir;
    if (dir==1)   /* filling foreground */
    {
      im0=1;
      im1=mri_fill->depth-1;
      i0=j0=1;
      i1=mri_fill->height ; j1=mri_fill->width ;
    } 
    else          /* filling background */
    {
      im0=mri_fill->depth-2;
      im1= -1;
      i0 = mri_fill->height - 2 ; j0 = mri_fill->width - 2 ;
      i1=j1= -1;
    }
    for (imnr=im0;imnr!=im1;imnr+=dir)
    {
      for (i=i0;i!=i1;i+=dir)
      {
        for (j=j0;j!=j1;j+=dir)
        {
          if ((j == 143) && (i == 106) && (imnr == 128))
              DiagBreak() ;
          if (MRIvox(mri_fill, j, i, imnr) ==0)   /* not filled yet */
          {
            if ((threshold<0 &&   /* previous filled off */
                 MRIvox(mri_im, j, i, imnr)<-threshold)  ||  
                (threshold>=0 &&
                 MRIvox(mri_im, j, i, imnr) >threshold))   /* wm is on */
            {
              /* three inside 6-connected nbrs */
              v1=MRIvox(mri_fill, j, i, imnr-dir);
              v2=MRIvox(mri_fill, j, i-dir, imnr);
              v3=MRIvox(mri_fill, j-dir, i, imnr) ;
              if (v1>0||v2>0||v3>0)       /* if any are on */
              {

                /* set vmax to biggest of three interior neighbors */
                vmax = (v1>=v2&&v1>=v3)?v1:((v2>=v1&&v2>=v3)?v2:v3);

                MRIvox(mri_fill, j, i, imnr) = vmax;
                nfilled++;
                ntotal++;
              }
            }
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%d voxels filled\n",nfilled);
  } 
  fprintf(stderr, "total of %d voxels filled...",ntotal);
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "\n") ;
  return(NO_ERROR) ;
}

static int
fill_holes(MRI *mri_fill)
{
  int  nfilled, ntotal = 0, cnt, cntmax= neighbor_threshold ;
  int im0,x0,i0,z,i,x;
  int v,vmax;

  do
  {
    nfilled = 0;
    for (z=1;z!=mri_fill->depth-1;z++)
    for (i=1;i!=mri_fill->height-1;i++)
    for (x=1;x!=mri_fill->width-1;x++)
    if (MRIvox(mri_fill, x, i, z)==0 && 
        i>ylim0-10 && i<ylim1+10 && x>xlim0-10 && x<xlim1+10)
    {
      cnt = 0;
      vmax = 0;
      for (im0= -1;im0<=1;im0++)
      for (i0= -1;i0<=1;i0++)
      for (x0= -1;x0<=1;x0++)
      {
        v = MRIvox(mri_fill, x+x0, i+i0, z+im0) ;
        if (v>vmax) vmax = v;
        if (v == 0) cnt++;              /* count # of nbrs which are off */
        if (cnt>cntmax) im0=i0=x0=1;    /* break out of all 3 loops */
      }
      if (x == 79 && i == 86 && z == 119)
        DiagBreak() ;
      if (cnt<=cntmax)   /* toggle pixel (off to on, or on to off) */
      {
        MRIvox(mri_fill, x, i, z) = vmax; 
        nfilled++;
        ntotal++;
      }
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%d holes filled\n",nfilled);
  } while (nfilled > 0) ;
  
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "total of %d holes filled\n",ntotal);
  return(NO_ERROR) ;
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
  if (!strcmp(option, "-help"))
    print_help() ;
  else if (!strcmp(option, "-version"))
    print_version() ;
  else if (!strcmp(option, "rval"))   /* sorry */
  {
    rh_fill_val = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d as fill val for right hemisphere.\n",
            rh_fill_val);
  }
  else if (!strcmp(option, "lval"))   /* sorry */
  {
    lh_fill_val = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using %d as fill val for left hemisphere.\n",lh_fill_val);
  }
  else if (!strcmp(option, "ccmask"))
  {
    cc_mask = 1 /*!cc_mask*/ ;
    fprintf(stderr,"%susing corpus callosum to mask possible location of "
            "pons.\n", cc_mask ? "" : "not ");
  }
  else if (!strcmp(option, "atlas"))
  {
    atlas_name = argv[2] ; sname = argv[3] ;
    nargs = 2 ;
    fprintf(stderr, "using atlas %s for auto filling...\n", atlas_name) ;
  }
  else switch (toupper(*option))
  {
  case 'P':
    pons_tal_x = atof(argv[2]) ;
    pons_tal_y = atof(argv[3]) ;
    pons_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "using pons seed point (%2.0f, %2.0f, %2.0f)\n",
            pons_tal_x, pons_tal_y, pons_tal_z) ;
    pons_seed_set = 1 ;
    break ;
  case 'C':
    cc_tal_x = atof(argv[2]) ;
    cc_tal_y = atof(argv[3]) ;
    cc_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "using corpus callosum seed point (%2.0f, %2.0f, %2.0f)\n",
            cc_tal_x, cc_tal_y, cc_tal_z) ;
    cc_seed_set = 1 ;
    break ;
  case 'T':
    if (sscanf(argv[2], "%d", &neighbor_threshold) < 1)
    {
      fprintf(stderr, "fill: could not scan threshold from '%s'", argv[2]) ;
      exit(1) ;
    }
    nargs = 1 ;
    break ;
  case 'L':
    log_fp = fopen(argv[2], "w") ;
    if (!log_fp)
      ErrorExit(ERROR_BADPARM, "%s: could not open cutting plane log file %s",
                Progname, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "logging cutting plane coordinates to %s...\n", argv[2]) ;
    break ;
  case 'F':
    fill_val = atoi(argv[2]) ;
    /*    min_filled = atoi(argv[2]) ;*/
    nargs = 1 ;
    break ;
  case 'D':
    logging = 1 ;
    break ;
  case '?':
  case 'U':
    print_help() ;
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
print_version(void)
{
  fprintf(stderr, "fill version %s\n", vcid) ;
  exit(0) ;
}
static void
print_help(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input MR directory> <output MR directory>\n",
          Progname) ;
  fprintf(stderr, "where options are:\n") ;
  fprintf(stderr, 
          "\t-T <threshold> - specify fill_holes threshold (default=%d)\n",
          DEFAULT_NEIGHBOR_THRESHOLD) ;
  fprintf(stderr, 
          "\t-L <x> <y> <z> - the Talairach coords of the white matter seed\n"
          "\t                 for the left hemisphere\n") ;
  fprintf(stderr, 
          "\t-R <x> <y> <z> - the Talairach coords of the white matter seed\n"
          "\t                 for the right hemisphere\n") ;
  fprintf(stderr, 
          "\t-P <x> <y> <z> - the Talairach coords of the seed for the "
          "pons\n");
  fprintf(stderr, 
          "\t-C <x> <y> <z> - the Talairach coords of the seed for the\n"
          "\t                 corpus callosum\n") ;
  exit(0) ;
}

/* build a set of slices in Talairach space, and find the one in which 
   the central connected component has the smallest cross-sectional
   area. Then use this as a mask for erasing ('cutting') stuff out
   of the input image.
*/

#define MAX_SLICES        15  /* 41*/
#define HALF_SLICES       ((MAX_SLICES-1)/2)
#define CUT_WIDTH         1
#define HALF_CUT          ((CUT_WIDTH-1)/2)
#define SEARCH_STEP       3
#define MAX_OFFSET        50

/* aspect ratios are dy/dx */
#define MIN_CC_AREA       350  /* smallest I've seen is 389 */
#define MAX_CC_AREA      1400  /* biggest I've seen is 1154 */
#define MIN_CC_ASPECT     0.1
#define MAX_CC_ASPECT     0.65

#define MIN_PONS_AREA     350  /* smallest I've seen is 385 */
#define MAX_PONS_AREA     700  /* biggest I've seen is 620 */  
#define MIN_PONS_ASPECT   0.6  /* was .8 */
#define MAX_PONS_ASPECT   1.2



static MRI *
find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,int orientation,
                   int *pxv, int *pyv, int *pzv, int seed_set)
{
  MRI        *mri_slices[MAX_SLICES], *mri_filled[MAX_SLICES], *mri_cut=NULL ;
  Real       dx, dy, dz, x, y, z, aspect,MIN_ASPECT,MAX_ASPECT,
             aspects[MAX_SLICES] ;
  int        slice, offset, area[MAX_SLICES], min_area, min_slice,xo,yo,
             xv, yv, zv, x0, y0, z0, xi, yi, zi, MIN_AREA, MAX_AREA, done ;
  FILE       *fp = NULL ;   /* for logging pons and cc statistics */
  char       fname[100], *name ;
  MRI_REGION region ;

  switch (orientation)
  {
  default:
  case MRI_SAGITTAL:   
    dx = 1.0 ; dy = dz = 0.0 ; 
    name = "corpus callosum" ;
    MIN_AREA = MIN_CC_AREA ; MAX_AREA = MAX_CC_AREA ;
    MIN_ASPECT = MIN_CC_ASPECT ; MAX_ASPECT = MAX_CC_ASPECT ;
    if (logging)
      fp = fopen(CC_LOG_FILE, "a") ;
    break ;
  case MRI_HORIZONTAL: 
    dz = 1.0 ; dx = dy = 0.0 ; 
    name = "pons" ;
    MIN_AREA = MIN_PONS_AREA ; MAX_AREA = MAX_PONS_AREA ;
    MIN_ASPECT = MIN_PONS_ASPECT ; MAX_ASPECT = MAX_PONS_ASPECT ;
    if (logging)
      fp = fopen(PONS_LOG_FILE, "a") ;
    break ;
  case MRI_CORONAL:    
    dy = 1.0 ; dx = dz = 0.0 ; 
    MIN_AREA = MIN_CC_AREA ; MAX_AREA = MAX_CC_AREA ;
    MIN_ASPECT = MIN_CC_ASPECT ; MAX_ASPECT = MAX_CC_ASPECT ;
    name = "coronal" ;
    break ;
  }

  xo = yo = (SLICE_SIZE-1)/2 ;  /* center point of the slice */

  /* search for valid seed point */
  MRItalairachToVoxel(mri, x_tal, y_tal,  z_tal, &x, &y, &z) ;
  xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
  mri_slices[0] = 
    MRIextractTalairachPlane(mri, NULL, orientation, xv, yv, zv, SLICE_SIZE);
  mri_filled[0] =MRIfillFG(mri_slices[0],NULL,xo,yo,0,WM_MIN_VAL,127,&area[0]);
  MRIboundingBox(mri_filled[0], 1, &region) ;

#if 0
#undef DIAG_VERBOSE_ON
#define DIAG_VERBOSE_ON  1
#endif
  
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s_seed.mgh", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_slices[0], fname) ;
    sprintf(fname, "%s_seed_fill.mgh", orientation==MRI_SAGITTAL?"cc":"pons");
    MRIwrite(mri_filled[0], fname) ;
  }
  
  /* now check to see if it could be a valid seed point based on:
     1) bound on the area
     2) the connected component is completely contained in the slice.
     */
  aspect = (Real)region.dy / (Real)region.dx ;
  done =  
    ((area[0] >= MIN_AREA) &&
     (area[0] <= MAX_AREA) &&
     (aspect  >= MIN_ASPECT) &&
     (aspect  <= MAX_ASPECT) &&
     (region.y > 0) &&
     (region.x > 0) &&
     (region.x+region.dx < SLICE_SIZE-1) &&
     (region.y+region.dy < SLICE_SIZE-1)) ;

  if (seed_set)  /* just fill around the specified seed */
  {
    done = 1 ;
    MRIboundingBox(mri_filled[0], 1, &region) ;
    if (orientation == MRI_SAGITTAL)/* extend to top and bottom of slice */
      region.dy = SLICE_SIZE - region.y ;
    MRItalairachToVoxel(mri, x_tal, y_tal,  z_tal, &x, &y, &z) ;
    xv = *pxv = nint(x) ; yv = *pyv = nint(y) ; zv = *pzv = nint(z) ;
    fprintf(stderr, "using seed (%d, %d, %d), TAL = (%2.1f, %2.1f, %2.1f)\n",
            *pxv, *pyv, *pzv, x_tal, y_tal, z_tal) ;
    mri_cut = MRIcopy(mri_filled[0], NULL) ;
    for (yv = region.y ; yv < region.y+region.dy ; yv++)
    {
      for (xv = region.x ; xv < region.x+region.dx ; xv++)
      {
        MRIvox(mri_cut, xv, yv, 0) = 1 ;
      }
    }
    xv = *pxv ; yv = *pyv ; zv = *pzv ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "%s_filled.mgh", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_slices[0], fname) ;
      sprintf(fname, "%s_cut.mgh", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_cut, fname) ;
    }
  }

  if (!seed_set && done)  /* center the seed */
  {
    find_slice_center(mri_filled[0], &xi, &yi) ;
    switch (orientation)   
    {
    default:
    case MRI_HORIZONTAL: xv += xi - xo ; zv += yi - yo ; break ;
    case MRI_SAGITTAL:   zv += xi - xo ; yv += yi - yo ; break ;
    }
    x = (Real)xv ; y = (Real)yv ; z = (Real)zv ;
    MRIvoxelToTalairach(mri, x, y, z, &x_tal, &y_tal, &z_tal) ;
  }

  MRIfree(&mri_slices[0]) ;
  MRIfree(&mri_filled[0]) ;

  offset = 0 ;
  if (!seed_set) while (!done)
  {
    offset += SEARCH_STEP ;   /* search at a greater radius */
    if (offset >= MAX_OFFSET)
      return(NULL) ;
    for (z0 = zv-offset ; !done && z0 <= zv+offset ; z0 += offset)
    {
      zi = mri->zi[z0] ;
      for (y0 = yv-offset ; !done && y0 <= yv+offset ; y0 += offset)
      {
        yi = mri->yi[y0] ;
        for (x0 = xv-offset ; !done && x0 <= xv+offset ; x0 += offset)
        {
          xi = mri->xi[x0] ;
          mri_slices[0] = 
            MRIextractTalairachPlane(mri, NULL, orientation, xi, yi, zi, 
                                     SLICE_SIZE);
          mri_filled[0] = 
            MRIfillFG(mri_slices[0], NULL, xo, yo, 0,WM_MIN_VAL,127, &area[0]);
          MRIboundingBox(mri_filled[0], 1, &region) ;

          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
          {
            sprintf(fname, "%s_seed.mgh", 
                    orientation == MRI_SAGITTAL ? "cc":"pons");
            MRIwrite(mri_slices[0], fname) ;
            sprintf(fname, "%s_seed_fill.mgh", 
                    orientation == MRI_SAGITTAL ? "cc":"pons");
            MRIwrite(mri_filled[0], fname) ;
          }
    
          /* now check to see if it could be a valid seed point based on:
             1) bound on the area
             2) the connected component is completely contained in the slice.
             */

          aspect = (Real)region.dy / (Real)region.dx ;
          if ((area[0] >= MIN_AREA) &&
              (area[0] <= MAX_AREA) &&
              (aspect  >= MIN_ASPECT) &&
              (aspect  <= MAX_ASPECT) &&
              (region.y > 0) &&
              (region.x > 0) &&
              (region.x+region.dx < SLICE_SIZE-1) &&
              (region.y+region.dy < SLICE_SIZE-1))
          {
            /* center the seed */
            find_slice_center(mri_filled[0], &xv, &yv) ;
            switch (orientation)   
            {
            default:
            case MRI_HORIZONTAL: xi += xv - xo ; zi += yv - yo ; break ;
            case MRI_SAGITTAL:   zi += xv - xo ; yi += yv - yo ; break ;
            }

            x = (Real)xi ; y = (Real)yi ; z = (Real)zi ;
            MRIvoxelToTalairach(mri, x, y, z, &x_tal, &y_tal, &z_tal) ;
            done = 1 ;
            xv = xi ; yv = yi ; zv = zi ;
          }
          MRIfree(&mri_slices[0]) ;
          MRIfree(&mri_filled[0]) ;
        }
      }
    }
  }

  if (Gdiag & DIAG_SHOW)
  {
    MRIvoxelToTalairach(mri, xv, yv, zv, &x, &y, &z) ;
    fprintf(stderr, 
          "%s seed point found at (%d, %d, %d), TAL: (%2.1f, %2.1f, %2.1f)\n", 
          name, xv,yv,zv, x, y, z);
  }

  if (!seed_set)   /* find slice with smallest cross-section */
  {
    for (slice = 0 ; slice < MAX_SLICES ; slice++)
    {
      offset = slice - HALF_SLICES ;
      x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ;
      MRItalairachToVoxel(mri, x, y,  z, &x, &y,&z) ;
      xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
      mri_slices[slice] = 
        MRIextractTalairachPlane(mri,NULL,orientation, xv, yv,zv,SLICE_SIZE);
      mri_filled[slice] = 
        MRIfillFG(mri_slices[slice],NULL,xo,yo,0,WM_MIN_VAL,127,&area[slice]);
      MRIboundingBox(mri_filled[slice], 1, &region) ;
      aspects[slice] = (Real)region.dy / (Real)region.dx ;

#if 0      
      /* don't trust slices that extend to the border of the image */
      if (!region.x || !region.y || region.x+region.dx >= SLICE_SIZE-1 ||
          region.y+region.dy >= SLICE_SIZE-1)
        area[slice] = 0 ;
#endif
      
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stderr, "slice[%d] @ (%d, %d, %d): area = %d\n", 
                slice, xv, yv, zv, area[slice]) ;
      
      if (orientation == MRI_SAGITTAL)/* extend to top and bottom of slice */
        region.dy = SLICE_SIZE - region.y ;
      
      /*    for (yv = region.y ; yv < region.y+region.dy ; yv++)*/
      for (yv = region.y ; yv < region.y+region.dy ; yv++)
      {
        for (xv = region.x ; xv < region.x+region.dx ; xv++)
        {
          MRIvox(mri_filled[slice], xv, yv, 0) = 1 ;
        }
      }
      
      if ((Gdiag & DIAG_WRITE) && !(slice % 1) && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s_slice%d.mgh", 
                orientation == MRI_SAGITTAL ? "cc":"pons", slice);
        MRIwrite(mri_slices[slice], fname) ;
        sprintf(fname, "%s_filled%d.mgh", 
                orientation == MRI_SAGITTAL ? "cc":"pons", slice);
        MRIwrite(mri_filled[slice], fname) ;
      }
    }
    
    if (orientation == MRI_HORIZONTAL)
    {
      min_area = 10000 ; min_slice = -1 ;
      for (slice = 1 ; slice < MAX_SLICES-1 ; slice++)
      {
        if (area[slice] < min_area && 
            (area[slice] >= MIN_AREA && area[slice] <= MAX_AREA))
        {
          min_area = area[slice] ;
          min_slice = slice ;
        }
      }
    }
    else    /* search for middle of corpus callosum */
    {
      int valid[MAX_SLICES], num_on, max_num_on, max_on_slice_start ;
      
      for (slice = 1 ; slice < MAX_SLICES-1 ; slice++)
      {
        valid[slice] =
          ((area[slice]    >= MIN_AREA)   && (area[slice]    <= MAX_AREA) &&
           (aspects[slice] >= MIN_ASPECT) && (aspects[slice] <= MAX_ASPECT));
      }
    
      max_on_slice_start = max_num_on = num_on = 0 ;
      for (slice = 1 ; slice < MAX_SLICES-1 ; slice++)
      {
        if (valid[slice])
          num_on++ ;
        else
        {
          if (num_on > max_num_on)
          {
            max_num_on = num_on ;
            max_on_slice_start = slice - num_on ;
          }
        }
      }
      if (num_on > max_num_on)   /* last slice was on */
      {
        max_num_on = num_on ;
        max_on_slice_start = slice - num_on ;
      }
      
      min_slice = max_on_slice_start + (max_num_on-1)/2 ;
      min_area = area[min_slice] ;
    }
    
    if (min_slice < 0)
      return(NULL) ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "%s_slice.mgh", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_slices[min_slice], fname) ;
      sprintf(fname, "%s_filled.mgh", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_filled[min_slice], fname) ;
    }

    if (logging)
    {
      fprintf(fp, "area %d, aspect %2.2f\n", 
              area[min_slice],aspects[min_slice]);
      fclose(fp) ;
    }
    
    offset = min_slice - HALF_SLICES ;
    x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
    MRItalairachToVoxel(mri, x, y,  z, &x, &y, &z) ;
    *pxv = nint(x) ; *pyv = nint(y) ; *pzv = nint(z) ;
    mri_cut = MRIcopy(mri_filled[min_slice], NULL) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, 
              "%s: cutting at slice %d, area %d, (%d, %d, %d)\n", 
              name, min_slice, min_area, *pxv, *pyv, *pzv) ;

    for (slice = 0 ; slice < MAX_SLICES ; slice++)
    {
      MRIfree(&mri_slices[slice]) ;
      MRIfree(&mri_filled[slice]) ;
    }
  }
  else   /* seed specified by user */
  {
  }

  return(mri_cut) ;
}

static int
find_slice_center(MRI *mri,  int *pxo, int *pyo)
{
  int  x, y, dist, min_total_dist, total_dist, width, height, x1, y1, xo, yo,
       border, xi, yi ;

  xo = yo = 0 ;
  width = mri->width ; height = mri->height ;
  min_total_dist = width*height*width*height ;
  for (y = 0 ; y < height ; y++)
  {
    for (x = 0 ; x < width ; x++)
    {
      /* find total distance to all other points */
      if (MRIvox(mri, x, y, 0) >= WM_MIN_VAL)
      {
        total_dist = 0 ;
        for (y1 = 0 ; y1 < height ; y1++)
        {
          for (x1 = 0 ; x1 < width ; x1++)
          {
            if (MRIvox(mri, x1, y1, 0) >= WM_MIN_VAL) 
            {
              dist = (x1-x)*(x1-x) + (y1-y)*(y1-y) ;
              total_dist += dist ;
            }
          }
        }
        if (total_dist < min_total_dist)
        {
          min_total_dist = total_dist ;
          xo = x ; yo = y ;
        }
      }
    }
  }

  /* now find a point which is not on the border */
  border = 1 ;
  for (y = yo-1 ; border && y <= yo+1 ; y++)
  {
    for (x = xo-1 ; border && x <= xo+1 ; x++)
    {
      if (MRIvox(mri, x, y, 0) >= WM_MIN_VAL) /* see if it is a border pixel */
      {
        border = 0 ;
        for (y1 = y-1 ; y1 <= y+1 ; y1++)
        {
          yi = mri->yi[y1] ;
          for (x1 = x-1 ; x1 <= x+1 ; x1++)
          {
            xi = mri->xi[x1] ;
            if (!MRIvox(mri, xi, yi, 0))
              border = 1 ;
          }
        }
      }
      if (!border)  /* (x,y) is not a border pixel */
      {
        xo = x ; yo = y ;
      }
    }
  }

  *pxo = xo ; *pyo = yo ;
  return(NO_ERROR) ;
}

#define CC_SPREAD       17
#define MIN_THICKNESS   3
#define MAX_THICKNESS   20

static int
find_corpus_callosum(MRI *mri, Real *pccx, Real *pccy, Real *pccz)
{
  MRI         *mri_slice ;
  int         xv, yv, zv, xo, yo, max_y, thickness, y1, xcc, ycc, x, y,x0 ;
  Real        xr, yr, zr ;
  MRI_REGION  region ;


  MRItalairachToVoxel(mri, 0.0, 0.0, 0.0, &xr, &yr, &zr);
  xv = nint(xr) ; yv = nint(yr) ; zv = nint(zr) ;
  mri_slice = 
    MRIextractTalairachPlane(mri, NULL, MRI_CORONAL,xv,yv,zv,SLICE_SIZE) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_slice, "cor.mgh") ;
  }
  MRIboundingBox(mri_slice, 1, &region) ;
  x0 = region.x+region.dx/2 ;
  xo = mri_slice->width/2 ;  yo = mri_slice->height/2 ;

  /* find the column with the lowest starting y value of any sign. thick. */
  xcc = ycc = max_y = 0 ; 
  for (x = x0-CC_SPREAD ; x <= x0+CC_SPREAD ; x++)
  {
    /* search for first non-zero pixel */
    for (y = 0 ; y < SLICE_SIZE ; y++)
    {
      if (MRIvox(mri_slice, x, y, 0) >= WM_MIN_VAL)
        break ;
    }
    
    /* check to make sure it as reasonably thick */
    if ((y < SLICE_SIZE) && (y > max_y))
    {
      for (y1 = y, thickness = 0 ; y1 < SLICE_SIZE ; y1++, thickness++)
        if (!MRIvox(mri_slice, x, y1, 0))
          break ;
      if ((thickness > MIN_THICKNESS) && (thickness < MAX_THICKNESS))
      {
        xcc = x ; ycc = y+thickness/2 ;  /* in middle of cc */
        max_y = y ;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stderr, "potential cc found at (%d, %d), thickness = %d\n",
                  xcc, ycc, thickness) ;
      }
    }
  }

  MRIfree(&mri_slice) ;
  if (!max_y)
    return(ERROR_BADPARM) ;

  /* now convert the in-plane coords to Talairach coods */
  xv += xcc - xo ; yv += ycc - yo ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "cc found at (%d, %d, %d)\n", xv, yv, zv) ;
  xr = (Real)xv ; yr = (Real)yv ; zr = (Real)zv ;
  MRIvoxelToTalairach(mri, xr, yr, zr, pccx, pccy, pccz) ;

  find_cc_slice(mri, pccx, pccy, pccz) ;
  return(NO_ERROR) ;
}

#define MIN_BRAINSTEM_THICKNESS         15
#define MAX_BRAINSTEM_THICKNESS         35
#define MIN_DELTA_THICKNESS             6
#define MIN_BRAINSTEM_HEIGHT            25

/* if can't find a really thick piece of brainstem (15 voxels wide),
   then try to find some contiguous slices that are fairly wide.
*/
#define MIN_CONT_BRAINSTEM_THICKNESS    8
#define MIN_CONT_BRAINSTEM_HEIGHT       3

static int
find_pons(MRI *mri, Real *p_ponsx, Real *p_ponsy, Real *p_ponsz, 
          int x_cc, int y_cc, int z_cc, int method)
{
  MRI   *mri_slice, *mri_filled, *mri_closed ;
  Real  xr, yr, zr ;
  int   xv, yv, zv, x, y, width, height, thickness, xstart, xo, yo, area,xmax,
        bheight ;
  MRI_REGION region ;

  thickness = xstart = -1 ;
  MRIboundingBox(mri, 10, &region) ;
  xr = (Real)(region.x+region.dx/2) ;
  yr = (Real)(region.y+region.dy/2) ;
  zr = (Real)(region.z+region.dz/2) ;
  xv = (int)xr ; yv = (int)yr ; zv = (int)zr ;
  xv = x_cc ; yv = y_cc ; zv = z_cc ;  /* use corpus callosum!! */
  
  MRIvoxelToTalairach(mri, xr, yr, zr, &xr, &yr, &zr);


  mri_slice = 
    MRIextractTalairachPlane(mri, NULL, MRI_SAGITTAL,xv,yv,zv,SLICE_SIZE) ;


  /* (xv,yv,zv) are the coordinates of the center of the slice */

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_slice, "brainstem_slice.mgh") ;

  switch (method)
  {
  case 1:
    mri_closed = MRIclose(mri_slice, NULL) ;  /* remove small holes */
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_closed, "brainstem_closed.mgh") ;

/*
  search from the front of the head at the bottom of a sagittal slice
  for the first set of contiguous slices that are significantly thick, 
  this will (hopefully) be the brainstem. 
*/
    
    /* first find the brainstem */
    width = mri_slice->width ; height = mri_slice->height ;
    bheight = 0 ;
    for (y = height-1 ; y >= 0 ; y--)
    {
      for (x = width-1 ; x >= 0 ; x--)
      {
        thickness = 0 ;
        while (MRIvox(mri_closed, x, y, 0) && x >= 0)
        {
          if (!thickness)
            xstart = x ;
          thickness++ ; x-- ;
        }
        if (thickness >= MIN_CONT_BRAINSTEM_THICKNESS && 
            thickness <= MAX_BRAINSTEM_THICKNESS)
          break ;
        else
          thickness = 0 ;
      }
      if (thickness > 0)   /* found the slice */
        bheight++ ;
      else
        bheight = 0 ;
      if (bheight >= MIN_CONT_BRAINSTEM_HEIGHT)
        break ;  /* found it */
    }
    MRIfree(&mri_closed) ;
    break ;
  default:
  case 0:
/*
  search from the front of the head at the bottom of a sagittal slice
  for the first significantly thick slice of white matter, which will
  be the brainstem. Then, follow the contour of the brainstem upwards
  until the first 'backwards' indentation which is about where we
  want to make the cut.
*/
    
    /* first find the brainstem */
    width = mri_slice->width ; height = mri_slice->height ;
    for (y = height-1 ; y >= 0 ; y--)
    {
      for (x = width-1 ; x >= 0 ; x--)
      {
        thickness = 0 ;
        while (MRIvox(mri_slice, x, y, 0) && x >= 0)
        {
          if (!thickness)
            xstart = x ;
          thickness++ ; x-- ;
        }
        if (thickness >= MIN_BRAINSTEM_THICKNESS && 
            thickness <= MAX_BRAINSTEM_THICKNESS)
          break ;
        else
          thickness = 0 ;
      }
      if (thickness > 0)   /* found the slice */
        break ;
    }
    break ;
  }

  /*
    Then, follow the contour of the brainstem upwards until the first 
    'backwards' indentation which is about where we want to make the cut.
  */

  mri_filled = 
    MRIfillFG(mri_slice, NULL,xstart-(thickness-1)/2,y,0,WM_MIN_VAL,127,&area);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "found brainstem at y=%d, x = %d --> %d, area=%d\n",
            y, xstart, xstart-thickness+1, area) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_filled, "brainstem_filled.mgh") ;

/* 
   now, starting from that slice, proceed upwards until we find a
   'cleft', which is about where we want to cut.
*/
  xmax = x = xstart ;  /* needed for initialization */
  bheight = 0 ;
  do
  {
    xstart = x ;
    y-- ; bheight++ ;
    for (x = width-1 ; !MRIvox(mri_filled,x,y,0) && (x >= 0) ; x--)
    {}

    if (x < 0)/* sometimes pons is disconnected at cleft: assume prior slice */
    {
      y++ ; bheight-- ;
      for (x = width-1 ; !MRIvox(mri_filled,x,y,0) && (x >= 0) ; x--)
      {}
      break ;
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice %d, xstart %d\n", y, x) ;
    if (x > xmax)
    {
      /* check to see if we've gone too far */
      if ((bheight > MIN_BRAINSTEM_HEIGHT) && (x-xmax >= MIN_DELTA_THICKNESS))
      {
        y-- ;   /* gone too far */
        break ;
      }
      xmax = x ;
    }
  } while ((xmax - x < MIN_DELTA_THICKNESS) && (y >0)) ;

  /* now search forward to find center of slice */
  for (thickness = 0, xstart = x ; 
       (x >= 0) && 
       (thickness < MAX_BRAINSTEM_THICKNESS) &&
       MRIvox(mri_slice, x, y, 0) ; 
       x--)
    thickness++ ;
  x = xstart - (thickness-1)/2 ;
  xo = mri_slice->width/2 ;  yo = mri_slice->height/2 ;
  zv += x - xo ; yv += y - yo ;  /* convert to voxel coordinates */
  xr = (Real)xv ; yr = (Real)yv ; zr = (Real)zv ;
  MRIvoxelToTalairach(mri, xr, yr, zr, p_ponsx, p_ponsy, p_ponsz) ;
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "putative pons found at (%d, %d) --> (%d, %d, %d)\n", 
            xstart-thickness/2, y, xv, yv, zv) ;

  MRIfree(&mri_slice) ;
  MRIfree(&mri_filled) ;
  return(NO_ERROR) ;
}

static int
find_cc_slice(MRI *mri, Real *pccx, Real *pccy, Real *pccz)
{
  int         area[MAX_SLICES], min_area, min_slice, slice, offset,xv,yv,zv,
              xo, yo ;
  MRI         *mri_slice, *mri_filled ;
  Real        aspect, x_tal, y_tal, z_tal, x, y, z ;
  MRI_REGION  region ;
  char        fname[100] ;

  x_tal = *pccx ; y_tal = *pccy ; z_tal = *pccz ;
  offset = 0 ;
  xo = yo = (SLICE_SIZE-1)/2 ;  /* center point of the slice */
  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    offset = slice - HALF_SLICES ;
    x = x_tal + offset ; y = y_tal ; z = z_tal ;
    MRItalairachToVoxel(mri, x, y,  z, &x, &y, &z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    mri_slice = 
      MRIextractTalairachPlane(mri, NULL, MRI_SAGITTAL, xv, yv, zv,SLICE_SIZE);
    mri_filled = 
      MRIfillFG(mri_slice, NULL, xo, yo,0,WM_MIN_VAL,127,&area[slice]);
    MRIboundingBox(mri_filled, 1, &region) ;
    aspect = (Real)region.dy / (Real)region.dx ;

    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= SLICE_SIZE-1 ||
        region.y+region.dy >= SLICE_SIZE-1)
      area[slice] = 0 ;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice[%d] @ (%d, %d, %d): area = %d\n", 
              slice, xv, yv, zv, area[slice]) ;

    if ((Gdiag & DIAG_WRITE) && !(slice % 1) && DIAG_VERBOSE_ON)
    {
      sprintf(fname, "cc_slice%d.mgh", slice);
      MRIwrite(mri_slice, fname) ;
      sprintf(fname, "cc_filled%d.mgh", slice);
      MRIwrite(mri_filled, fname) ;
    }
    MRIfree(&mri_filled) ; MRIfree(&mri_slice) ;
  }

  min_area = 10000 ; min_slice = -1 ;
  for (slice = 1 ; slice < MAX_SLICES-1 ; slice++)
  {
    if (area[slice] < min_area && 
        (area[slice] >= MIN_CC_AREA && area[slice] <= MAX_CC_AREA))
    {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }

  /* couldn't find a good slice - don't update estimate */
  if (min_slice < 0)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: could not find valid seed for the cc",
                 Progname));

  offset = min_slice - HALF_SLICES ;
  *pccx = x = x_tal + offset ; *pccy = y = y_tal ; *pccz = z = z_tal ;
  MRItalairachToVoxel(mri, x, y,  z, &x, &y, &z) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, 
            "updating initial cc seed to (%d, %d, %d)\n", 
            nint(x), nint(y), nint(z)) ;
  return(NO_ERROR) ;
}

static int
neighbors_on(MRI *mri, int x0, int y0, int z0)
{
  int   nbrs = 0 ;

  if (MRIvox(mri,x0-1,y0,z0) >= WM_MIN_VAL)
    nbrs++ ;
  if (MRIvox(mri,x0+1,y0,z0) >= WM_MIN_VAL)
    nbrs++ ;
  if (MRIvox(mri,x0,y0+1,z0) >= WM_MIN_VAL)
    nbrs++ ;
  if (MRIvox(mri,x0,y0-1,z0) >= WM_MIN_VAL)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0+1) >= WM_MIN_VAL)
    nbrs++ ;
  if (MRIvox(mri,x0,y0,z0-1) >= WM_MIN_VAL)
    nbrs++ ;
  return(nbrs) ;
}
#if 0
static int mriExtendMaskDownward(MRI *mri) ;
static int
mriExtendMaskDownward(MRI *mri)
{
  int     width, height, depth, x, y, z, y1 ;
  BUFTYPE *psrc, val ;

  width = mri->width ;
  height = mri->height ;
  depth = mri->depth ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      psrc = &MRIvox(mri, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if ((val = *psrc++) > 0)
        {
          for (y1 = y ; y1 < mri->height ; y1++)
            MRIvox(mri, x, y1, z) = val ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

#endif

/*
  take the volume in mri_im and find the 'on' connected component
  containing the point (x_seed, y_seed, z_seed), then fill interior holds
  and remove exterior islands to ensure that it contains one connected
  mass (this will be one hemisphere).
*/
static int
MRIfillVolume(MRI *mri_fill, MRI *mri_im, int x_seed, int y_seed, int z_seed,
              int fill_val)
{
  int   x, y, z ;

#if 0
  MRIvox(mri_fill, wm_lh_x, wm_lh_y, wm_lh_z) = lh_fill_val ;
  MRIvox(mri_fill, wm_rh_x, wm_rh_y, wm_rh_z) = rh_fill_val ;
#else
  MRIvox(mri_fill, x_seed, y_seed, z_seed) = fill_val ;
#endif
  if (!Gdiag)
    fprintf(stderr, "filling volume: pass 1 of 3...") ;
  /* fill from 2 seeds in wm outwards */
  fill_brain(mri_fill, mri_im, WM_MIN_VAL);

  /*MRIwrite(mri_fill, "fill1.mgh") ;*/
  
  /* Set im to initial outward fill */
  MRIcopy(mri_fill, mri_im) ;
  MRIclear(mri_fill) ;
  MRIvox(mri_fill,1,1,1) = 255;              /* seed for background */

  /* fill in connected component of background */  
  if (!Gdiag)
    fprintf(stderr, "\rfilling volume: pass 2 of 3...") ;
  fill_brain(mri_fill, mri_im, -WM_MIN_VAL);

  if (fill_holes_flag)              
    fill_holes(mri_fill);   /* fill in islands with fewer than 10 nbrs on */
  
  /* put complement into im (im==on means part of background) */
  MRIcopy(mri_fill, mri_im) ;
  MRIclear(mri_fill) ;

  MRIvox(mri_fill, x_seed, y_seed, z_seed) = fill_val ;

  /* fill in background of complement (also sometimes called the foreground) */
  if (!Gdiag)
    fprintf(stderr, "\rfilling volume: pass 3 of 3...") ;
  fill_brain(mri_fill, mri_im, -WM_MIN_VAL);   
  mriFindBoundingBox(mri_fill) ;
  
  if (fill_holes_flag)
    fill_holes(mri_fill);
  
  if (!Gdiag)
    fprintf(stderr, "done.\n") ;
  for (z = 0 ; z < mri_fill->depth ; z++)
    for (y = 0; y < mri_fill->height ; y++)
      for (x = 0; x < mri_fill->width ; x++)
      {
        if (z==0||z==mri_fill->depth-1) 
          MRIvox(mri_fill, x, y, z) = 0;
      }
  return(NO_ERROR) ;
}

/*
  NOTE: this procedure sets some global variable (xlim0, ylim0, xlim1, ylim1)
  which are used in a couple of places!!!! Sorry, but don't blame me -
  it was from Anders' and Marty's original code....
*/
static int
mriFindBoundingBox(MRI *mri_im)
{
  int  x, y, z ;

  /* find bounding box for valid data */
  ylim0 = mri_im->height ; xlim0 = mri_im->width ;
  ylim1 = xlim1 = 0;
  for (z = 0 ; z < mri_im->depth ; z++)
  {
    for (y = 0 ; y < mri_im->height; y++)
    {
      for (x = 0 ; x < mri_im->width ; x++)
      {
        if (MRIvox(mri_im, x, y, z) >= WM_MIN_VAL)
        {
          if (y<ylim0) ylim0=y;
          if (y>ylim1) ylim1=y;
          if (x<xlim0) xlim0=x;
          if (x>xlim1) xlim1=x;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

/*
  combine the two filled hemispheres into a single filled volume.
  Overlap in the volumes is resolved by filling from the most lateral
  part of each hemisphere, and making the fill variable so that wide 
  connections are filled before narrow ones (defined by the # of neighbors
  that are on). Note that this routine modifies mri_lh and mri_rh so
  that they do not overlap.
*/

MRI *
MRIcombineHemispheres(MRI *mri_lh, MRI *mri_rh, MRI *mri_dst, Real cc_tal_x)
{
  int     width, height, depth, x, y, z, v, ncorrected, ambiguous, lh, rh,
          nfilled, xi, yi, zi, xk, yk, zk, nleft, nright, nambiguous, i ;
  BUFTYPE *plh, *prh, *pdst, vlh, vrh ;
  Real    xtal, ytal, ztal, dist ;

  width = mri_lh->width ;
  height = mri_lh->height ;
  depth = mri_lh->depth ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_lh, NULL) ;

  ambiguous = rh = lh = 0 ;
  for (ncorrected = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      plh = &MRIvox(mri_lh, 0, y, z) ;
      prh = &MRIvox(mri_rh, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        vlh = *plh++ ;
        vrh = *prh++ ;
        if (vlh >= WM_MIN_VAL && vrh >= WM_MIN_VAL)
        {
          if (ambiguous == 0)
          {
            ambiguous = (vlh+vrh) / 2 ;
            lh = vlh ;
            rh = vrh ;
          }
          /* both on - use talairach coordinate to arbitrate */
          MRIvoxelToTalairach(mri_lh, (Real)x, (Real)y,(Real)z,
                              &xtal,&ytal,&ztal) ;
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
            fprintf(stderr, "arbitrating overlap at (%d,%d,%d) [tal "
                    "(%2.1f,%2.1f,%2.1f] - ",x,y,z,(float)xtal,(float)ytal,
                    (float)ztal) ;

          dist = (xtal-cc_tal_x) ;
          if (fabs(dist) < 30)
            v = ambiguous ;
          else if (dist < 0)
            v = vlh ;
          else
            v = vrh ;
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
            fprintf(stderr, "%s\n", v == vlh ? "lh" : "rh") ;
          ncorrected++ ;
        }
        else
          v = MAX(vlh, vrh) ;
        *pdst++ = v ;
      }
    }
  }
  if (ncorrected > 0)   /* correct ambiguous voxels using diffusion distance */
  {
    MRI *mri_tmp = NULL ;
    int max_on = 0, min_on = (3*3*3)-1 ;

    fprintf(stderr, 
            "using variable coefficient diffusion "
            "to correct hemispheric overlap...\n") ;

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      fprintf(stderr, "writing ambiguous voxel image to amb.mgh\n") ;
      MRIwrite(mri_dst, "amb.mgh") ;
    }
    mri_tmp = MRIcopy(mri_dst, mri_tmp) ;
    do
    {
      max_on = i = nambiguous = nfilled = 0 ;
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          pdst = &MRIvox(mri_dst, 0, y, z) ;
          for (x = 0 ; x < width ; x++)
          {
            v = *pdst++ ;
            if (v == ambiguous)
            {
              nleft = nright = 0 ;
              for (zk = -1 ; zk <= 1 ; zk++)
              {
                zi = mri_dst->zi[z+zk] ;
                for (yk = -1 ; yk <= 1 ; yk++)
                {
                  yi = mri_dst->yi[y+yk] ;
                  for (xk = -1 ; xk <= 1 ; xk++)
                  {
                    xi = mri_dst->xi[x+xk] ;
                    v = MRIvox(mri_dst, xi, yi, zi) ;
                    if (v == lh)
                      nleft++ ;
                    else if (v == rh)
                      nright++ ;
                  }
                }
              }
              if (MAX(nright,nleft) > max_on)
                max_on = MIN(min_on, MAX(nright,nleft)) ;
              if (nright >= min_on && nleft >= min_on && nright == nleft)
              {
                if (ODD(i))
                  nleft-- ;
                else
                  nright-- ;
                i++ ;
              }
              if (nright > nleft && nright >= min_on)
              {
                nfilled++ ;
                MRIvox(mri_tmp, x, y, z) = rh ;
                MRIvox(mri_lh, x, y, z) = 0 ;
              }
              else if (nleft > nright && nleft >= min_on)
              {
                nfilled++ ;
                MRIvox(mri_tmp, x, y, z) = lh ;
                MRIvox(mri_rh, x, y, z) = 0 ;
              }
              else
                nambiguous++ ;
            }
          }
        }
      }
      MRIcopy(mri_tmp, mri_dst) ;
      if (min_on != max_on || DIAG_VERBOSE_ON)
        fprintf(stderr, "%d:  %d ambiguous voxels remaining\n", 
                min_on, nambiguous) ;
      min_on = max_on ;
    } while ((nfilled > 0) || (min_on > 0)) ;
  }

  if (Gdiag & DIAG_SHOW && ncorrected > 0)
    fprintf(stderr, "%d overlapping voxels corrected\n", ncorrected) ;
  return(mri_dst) ;
}


static MRI *
mriReadConditionalProbabilities(MRI *mri_T1, char *atlas_name, char *suffix, 
                                int offset, M3D *m3d, MRI *mri_dst)
{
  MRI  *mri_mean, *mri_std, *mri_tmp ;
  char *mri_dir, fname[STRLEN] ;

  mri_dir = getenv("MRI_DIR") ;
  if (!mri_dir)
    ErrorExit(ERROR_BADPARM, "%s: could not read MRI_DIR from environment\n") ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_T1, NULL) ;

  /* read and transform the mean volume */
  sprintf(fname, "%s/average/%s_%s.mgh#%d",
          mri_dir, atlas_name, suffix, offset) ;
  fprintf(stderr, "reading atlas means from %s...\n", fname) ;
  mri_mean = MRIread(fname) ;
  if (!mri_mean)
    ErrorExit(ERROR_NOFILE, "%s: could not read mean file %s",
              Progname, fname) ;
  if (m3d)
  {
    mri_tmp = MRIapplyInverse3DMorph(mri_mean, m3d, NULL) ;
    MRIfree(&mri_mean) ; mri_mean = mri_tmp ;
  }

  /* read and transform the standard deviation volume */
  sprintf(fname, "%s/average/%s_%s.mgh#%d",
          mri_dir, atlas_name, suffix, offset+1) ;
  fprintf(stderr, "reading atlas sigmas from %s...\n", fname) ;
  mri_std = MRIread(fname) ;
  if (!mri_std)
    ErrorExit(ERROR_NOFILE, "%s: could not read std file %s",
              Progname, fname) ;
  if (m3d)
  {
    mri_tmp = MRIapplyInverse3DMorph(mri_std, m3d, NULL) ;
    MRIfree(&mri_std) ; mri_std = mri_tmp ;
  }
  mri_dst = 
    MRIcomputeConditionalProbabilities(mri_T1, mri_mean, mri_std, NULL) ;
  MRIfree(&mri_mean) ; MRIfree(&mri_std) ;
  return(mri_dst) ;
}

static MRI *
mriReadBinaryProbabilities(char *atlas_name, char *suffix, M3D *m3d, 
                           char *subject_name, MRI *mri_dst)
{
  MRI  *mri_priors, *mri_T1, *mri_on_conditional, *mri_off_conditional, *mri_tmp ;
  char *mri_dir, *subjects_dir, fname[STRLEN] ;

  mri_dir = getenv("MRI_DIR") ;
  if (!mri_dir)
    ErrorExit(ERROR_BADPARM, "%s: could not read MRI_DIR from environment\n") ;
  subjects_dir = getenv("SUBJECTS_DIR") ;
  if (!subjects_dir)
    ErrorExit(ERROR_BADPARM, 
              "%s: could not read SUBJECTS_DIR from environment\n") ;

  /* read in subject's T1 volume */
  sprintf(fname, "%s/%s/mri/T1", subjects_dir, subject_name) ;
  fprintf(stderr, "reading T1 volume from %s...\n", fname) ;
  mri_T1 = MRIread(fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 volume from %s",
              Progname, fname) ;

  /* read and transform the standard deviation volume */
  sprintf(fname, "%s/average/%s_%s.mgh#0", mri_dir, atlas_name, suffix) ;
  fprintf(stderr, "reading atlas priors from %s offset 0...\n", fname) ;
  mri_priors = MRIread(fname) ;
  if (!mri_priors)
    ErrorExit(ERROR_NOFILE, "%s: could not read priors file %s",
              Progname, fname) ;
  if (m3d)
  {
    mri_tmp = MRIapplyInverse3DMorph(mri_priors, m3d, NULL) ;
    MRIfree(&mri_priors) ; mri_priors = mri_tmp ;
  }

  mri_on_conditional = 
    mriReadConditionalProbabilities(mri_T1, atlas_name, suffix, 1, m3d, NULL) ;
  mri_off_conditional = 
    mriReadConditionalProbabilities(mri_T1, atlas_name, suffix, 3, m3d, NULL) ;
  mri_dst = 
    MRIapplyBayesLaw(mri_priors,mri_on_conditional,mri_off_conditional,NULL);
  MRIfree(&mri_T1) ; MRIfree(&mri_on_conditional) ;
  MRIfree(&mri_priors) ; MRIfree(&mri_off_conditional) ;
  return(mri_dst) ;
}

static MRI *
MRIfillVentricle(MRI *mri_src, MRI *mri_prob, MRI *mri_T1, MRI *mri_dst, 
                 float threshold, int out_label, int xabs_min, int xabs_max)
{
  BUFTYPE   *pprob, *pdst, *psrc, out_val, prob, in_val ;
  int       width, height, depth, x, y, z, ymin, ymax, zmin, zmax, 
            nfilled, xk, yk, zk, niter, xmin, xmax, expanded, expansions,
            xi, yi, zi, ventricle_voxels, total_expanded ;

  if (mri_prob->type != MRI_UCHAR)
    ErrorReturn(NULL, (ERROR_UNSUPPORTED, 
                       "MRI3Dthreshold: prob must be MRI_UCHAR")) ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ; 
  /* now apply the inverse morph to build an average wm representation
     of the input volume 
     */


  total_expanded = ventricle_voxels = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      pprob = &MRIvox(mri_prob, 0, y, z) ;
      psrc = &MRIvox(mri_src, 0, y, z) ;
      pdst = &MRIvox(mri_dst, 0, y, z) ;
      for (x = 0 ; x < width ; x++)
      {
        if (DEBUG_POINT(x,y,z))
          DiagBreak() ;
        out_val = 0 ;
        prob = *pprob++ ;   /* value from inverse morphed volume */
        in_val = *psrc++ ;
        if (prob >= threshold)        /* probably on */
          out_val = out_label ;
        else                         /* not sure, use original val */
          out_val = 0 ;

        if (out_val)
          ventricle_voxels++ ;
        *pdst++ = out_val ;
      }
    }
  }

  nfilled = niter = 0 ;
  xmin = width ; ymin = height ; zmin = depth ;
  xmax = ymax = zmax = 0 ;
  do
  {
    expansions = expanded = 0 ;

    /* compute bounding box */
    if (!nfilled)  /* couldn't fill anymore in this bbox - expand it */
    {
      xmin = width ; ymin = height ; zmin = depth ;
      xmax = ymax = zmax = 0 ;
      for (z = 0 ; z < depth ; z++)
      {
        for (y = 0 ; y < height ; y++)
        {
          for (x = 0 ; x < width ; x++)
          {
            if (x == 151 && y == 139 && z == 69)
              DiagBreak() ;
            if (MRIvox(mri_dst, x, y, z) == out_label)
            {
              if (x < xmin)
                xmin = x ;
              if (y < ymin)
                ymin = y ;
              if (z < zmin)
                zmin = z ;
              if (x > xmax)
                xmax = x ;
              if (y > ymax)
                ymax = y ;
              if (z > zmax)
                zmax = z ;
            }
          }
        }
      }
      expanded = 1 ; expansions++ ;
    }
    else
      expanded = 0 ;

    nfilled = 0 ;

    /* don't let bounding box grow outside limits imposed by caller.
       This will (hopefully) prevent the filling from expanding medially
       and escaping into the third ventricle.
    */
    if (xmin < xabs_min)
      xmin = xabs_min ;
    if (xmax > xabs_max)
      xmax = xabs_max ;

    if (!niter++)  /* first time through - shrink filled region */
    {
#define MAX_EXPANSIONS 5
#define BUFFER         10

      for (z = zmin ; z <= zmax ; z++)
      {
        for (y = ymin ; y <= ymax ; y++)
        {
          for (x = xmin ; x <= xmax ; x++)
          {
            /* if close to one of the borders */
            if ((((x-xmin) <= BUFFER) || ((xmax-x) <= BUFFER)) ||
                (((z-zmin) <= BUFFER) || ((zmax-z) <= BUFFER)) ||
                (((y-ymin) <= BUFFER) || ((ymax-y) <= BUFFER)))
            {
              if (x == 151 && y == 139 && z == 69)
                DiagBreak() ;
              if (MRIvox(mri_dst, x, y, z) == out_label)
                MRIvox(mri_dst, x, y, z) = 0 ;
            }
          }
        }
      }
      xmin += BUFFER ; xmax -= BUFFER ;
      ymin += BUFFER ; ymax -= BUFFER ;
      zmin += BUFFER ; zmax -= BUFFER ;
#if 1
      if (out_label == lh_fill_val)
        MRIwrite(mri_dst, "left_ventricle_p_fill.mgh") ;
      else
        MRIwrite(mri_dst, "right_ventricle_p_fill.mgh") ;
#endif
    }

    fprintf(stderr, "ventricle bounding box [%d, %d, %d] --> [%d, %d, %d]\n",
            xmin, ymin, zmin, xmax, ymax, zmax) ;
    for (z = zmin ; z <= zmax ; z++)
    {
      for (y = ymin ; y <= ymax ; y++)
      {
        for (x = xmin ; x <= xmax ; x++)
        {
          if (x == 122 && y == 109 && z == 117)
            DiagBreak() ;
          if (MRIvox(mri_dst, x, y, z) == out_label)
          {   /* a filled voxel - expand fill */
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = mri_dst->zi[z+zk] ;
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                yi = mri_dst->yi[y+yk] ;
                for (xk = -1 ; xk <= 1 ; xk++)
                {
                  xi = mri_dst->xi[x+xk] ;
                  if (x == 122 && y == 108 && z == 117)
                    DiagBreak() ;
                  if ((MRIvox(mri_T1, xi, yi, zi) < 50) &&
                      (MRIvox(mri_dst, xi, yi, zi) != out_label))
                  {
                    MRIvox(mri_dst, xi, yi, zi) = out_label ;
                    nfilled++ ;
                  }
                }
              }
            }
          }
        }
      }
    }
    fprintf(stderr, "%d ventricular voxels filled\n", nfilled) ;
    total_expanded += nfilled ;

    if ((ventricle_voxels + total_expanded > 20000) ||
        (++expansions >= MAX_EXPANSIONS))
      break ;

  } while (!expanded || (nfilled > 0)) ;
  
  ventricle_voxels += total_expanded ;
  fprintf(stderr, "%5d ventricular voxels, %5d from atlas %5d from fill.\n", 
          ventricle_voxels, ventricle_voxels-total_expanded, total_expanded);

  return(mri_dst) ;
}
