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

/*-------------------------------------------------------------------
                                CONSTANTS
-------------------------------------------------------------------*/

/* min # of neighbors which must be on to retain a point */
#define DEFAULT_NEIGHBOR_THRESHOLD      8

/* seed values */
#define LEFT_HEMISPHERE_WHITE_MATTER    127
#define RIGHT_HEMISPHERE_WHITE_MATTER   255

/* distance to search in each direction for a valid wm seed point */
#define SEED_SEARCH_SIZE                9

/* size of various orthogonal slices */
#define SLICE_SIZE                      128

/* anything below this is not white matter */
#define WM_MIN_VAL                       2 

/* Talairach seed points - only used if heuristics fail */
#define CORPUS_CALLOSUM_TAL_X            0.0
#define CORPUS_CALLOSUM_TAL_Y            0.0
#define CORPUS_CALLOSUM_TAL_Z            27.0

#define PONS_TAL_X                      -2.0
#define PONS_TAL_Y                      -15.0  /* was -22.0 */
#define PONS_TAL_Z                      -17.0

/* # if iterations for filling */
#define MAX_ITERATIONS  1000  /* was 10 */

/* min # of holes filled before quitting */
#define MIN_FILLED      0     /* was 100 */

/*-------------------------------------------------------------------
                                GLOBAL DATA
-------------------------------------------------------------------*/

static int ilim0,ilim1,jlim0,jlim1;
static int fill_holes_flag = TRUE;

/* Talairach seed points for white matter in left and right hemispheres */
static Real wm_lh_tal_x = 29.0 ;
static Real wm_lh_tal_y = -12.0 ;
static Real wm_lh_tal_z = 28.0 ;

static Real wm_rh_tal_x = -29.0 ;
static Real wm_rh_tal_y = -12.0 ;
static Real wm_rh_tal_z = 28.0 ;

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

/*-------------------------------------------------------------------
                             STATIC PROTOTYPES
-------------------------------------------------------------------*/

static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;
void main(int argc, char *argv[]) ;

static int fill_holes(void) ;
static int fill_brain(int threshold) ;
void main(int argc, char *argv[]) ;
static MRI *find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,
                              int orientation, int *pxv, int *pvy, int *pzv) ;
static int find_slice_center(MRI *mri,  int *pxo, int *pyo) ;
static int find_corpus_callosum(MRI *mri, Real *ccx, Real *ccy, Real *ccz) ;
static int find_pons(MRI *mri, Real *p_ponsx, Real *p_ponsy, Real *p_ponsz) ;

/*-------------------------------------------------------------------
                                FUNCTIONS
-------------------------------------------------------------------*/

void
main(int argc, char *argv[])
{
  int     i,j,k, xd, yd, zd, xnew, ynew, znew ;
  int     nargs, wm_lh_x, wm_lh_y, wm_lh_z, wm_rh_x, wm_rh_y, wm_rh_z ;
  char    input_fname[200],out_fname[200],*data_dir;
  Real    x, y, z, dist, min_dist ;
  MRI     *mri_cc, *mri_pons ;
  int     x_pons, y_pons, z_pons, x_cc, y_cc, z_cc, xi, yi, zi ;
#if 0
  int     imnr_seed[MAXSEED],i_seed[MAXSEED],j_seed[MAXSEED],val_seed[MAXSEED];
  int     snum, nseed ;
#endif

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

#if 0
  if (argc > 2)
    snum = atoi(argv[2]);
  else
    snum = 0 ;
#endif
  
  data_dir = getenv("SUBJECTS_DIR");
  if (data_dir==NULL)
#if 1
    data_dir = "." ;
#else
  {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(0);
  }
#endif
#if 0
  sprintf(input_fname,"%s/%s/mri/wm/COR-",data_dir,argv[1]);
  sprintf(out_fname,"%s/%s/mri/filled/COR-",data_dir,argv[1]);
#else
  strcpy(input_fname, argv[1]) ;
  strcpy(out_fname, argv[2]) ;
#endif

  mri_im = MRIread(input_fname) ;
  if (!mri_im)
    ErrorExit(ERROR_NOFILE, "%s: could not read %s", Progname, input_fname) ;
  mri_fill = MRIclone(mri_im, NULL) ;

  find_pons(mri_im, &pons_tal_x, &pons_tal_y, &pons_tal_z) ;
  mri_pons = 
    find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                       MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons);
  if (!mri_pons) /* heuristic failed to find the pons - try Talairach coords */
  {
    pons_tal_x = PONS_TAL_X ; pons_tal_y = PONS_TAL_Y; pons_tal_z = PONS_TAL_Z;
    mri_pons = 
      find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                         MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons);
    if (!mri_pons)
      ErrorExit(ERROR_BADPARM, "%s: could not find pons", Progname);
  }

  if (find_corpus_callosum(mri_im,&cc_tal_x,&cc_tal_y,&cc_tal_z) != NO_ERROR)
    mri_cc = NULL ;
  else
    mri_cc =
      find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z, MRI_SAGITTAL,
                         &x_cc, &y_cc, &z_cc) ;
  if (!mri_cc)  /* heuristic failed - use Talairach coordinates */
  {
    cc_tal_x = CORPUS_CALLOSUM_TAL_X ; 
    cc_tal_y = CORPUS_CALLOSUM_TAL_Y; 
    cc_tal_z = CORPUS_CALLOSUM_TAL_Z;
    mri_cc = 
      find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z, MRI_SAGITTAL,
                         &x_cc, &y_cc, &z_cc) ;
    if (!mri_cc)
      ErrorExit(ERROR_BADPARM, "%s: could not find corpus callosum", Progname);
  }

  MRIeraseTalairachPlane(mri_im, mri_cc, MRI_SAGITTAL, x_cc, y_cc, z_cc, 
                         SLICE_SIZE);
  MRIeraseTalairachPlane(mri_im, mri_pons, MRI_HORIZONTAL, 
                         x_pons, y_pons, z_pons, SLICE_SIZE) ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_pons, "pons.mnc") ;
    MRIwrite(mri_cc, "cc.mnc") ;
  }
  MRIfree(&mri_cc) ;
  MRIfree(&mri_pons) ;

  /* find bounding box for valid data */
  ilim0 = mri_im->height ; jlim0 = mri_im->width ;
  ilim1 = jlim1 = 0;
  for (k = 0 ; k < mri_im->depth ; k++)
  {
    for (i = 0 ; i < mri_im->height; i++)
    {
      for (j = 0 ; j < mri_im->width ; j++)
      {
        if (MRIvox(mri_im, j, i, k) >0)
        {
          if (i<ilim0) ilim0=i;
          if (i>ilim1) ilim1=i;
          if (j<jlim0) jlim0=j;
          if (j>jlim1) jlim1=j;
        }
      }
    }
  }

  /* find white matter seed point for the left hemisphere */
  MRItalairachToVoxel(mri_im, wm_lh_tal_x, wm_lh_tal_y, wm_lh_tal_z,&x, &y,&z);
  wm_lh_x = nint(x) ; wm_lh_y = nint(y) ; wm_lh_z = nint(z) ;
  if (!MRIvox(mri_im, wm_lh_x, wm_lh_y, wm_lh_z))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for lh wm seed") ;
    for (k = wm_lh_z-SEED_SEARCH_SIZE ; k <= wm_lh_z+SEED_SEARCH_SIZE ; k++)
    {
      zi = mri_im->zi[k] ;
      for (j = wm_lh_y-SEED_SEARCH_SIZE ; j <= wm_lh_y+SEED_SEARCH_SIZE ; y++)
      {
        yi = mri_im->yi[j] ;
        for (i = wm_lh_x-SEED_SEARCH_SIZE ;i <= wm_lh_x+SEED_SEARCH_SIZE ; i++)
        {
          xi = mri_im->xi[i] ;
          if (MRIvox(mri_im, xi, yi, zi))
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
    wm_lh_x = xnew ; wm_lh_y = ynew ; wm_lh_z = znew ; 
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "found at (%d, %d, %d)\n", xnew, ynew, znew) ;
  }

  /* find white matter seed point for the right hemisphere */
  MRItalairachToVoxel(mri_im, wm_rh_tal_x, wm_rh_tal_y, wm_rh_tal_z,&x, &y,&z);
  wm_rh_x = nint(x) ; wm_rh_y = nint(y) ; wm_rh_z = nint(z) ;
  if (!MRIvox(mri_im, wm_rh_x, wm_rh_y, wm_rh_z))
  {
    xnew = ynew = znew = 0 ;
    min_dist = 10000.0f ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "searching for rh wm seed") ;
    for (k = wm_rh_z-SEED_SEARCH_SIZE ; k <= wm_rh_z+SEED_SEARCH_SIZE ; k++)
    {
      zi = mri_im->zi[k] ;
      for (j = wm_rh_y-SEED_SEARCH_SIZE ; j <= wm_rh_y+SEED_SEARCH_SIZE ; j++)
      {
        yi = mri_im->yi[j] ;
        for (i = wm_rh_x-SEED_SEARCH_SIZE ;i <= wm_rh_x+SEED_SEARCH_SIZE ; i++)
        {
          xi = mri_im->xi[i] ;
          if (MRIvox(mri_im, xi, yi, zi))
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
    wm_rh_x = xnew ; wm_rh_y = ynew ; wm_rh_z = znew ; 
  }

  /* initialize the fill with the detected seed points */
  MRIvox(mri_fill, wm_rh_x, wm_rh_y, wm_rh_z) = RIGHT_HEMISPHERE_WHITE_MATTER ;
  MRIvox(mri_fill, wm_lh_x, wm_lh_y, wm_lh_z) = LEFT_HEMISPHERE_WHITE_MATTER ;

  fill_brain(WM_MIN_VAL);    /* fill from 2 seeds in wm outwards */

  /*MRIwrite(mri_fill, "fill1.mnc") ;*/
  
  /* set im to initial outward fill */
  MRIcopy(mri_fill, mri_im) ;
  MRIclear(mri_fill) ;
  MRIvox(mri_fill,1,1,1) = 255;              /* seed for background */
  
  fill_brain(-WM_MIN_VAL);/* fill in connected component of background */
  /*MRIwrite(mri_fill, "fill2.mnc") ;*/
  if (fill_holes_flag)              
    fill_holes();         /* fill in islands with fewer than 10 nbrs on */
  
  /*MRIwrite(mri_fill, "fill3.mnc") ;*/

  /* put complement into im (im==on means part of background) */
  MRIcopy(mri_fill, mri_im) ;
  MRIclear(mri_fill) ;

#if 0 
  for (i=0;i<nseed;i++)
    if (snum==0 || snum==i+1)
    {
      if ((j_seed[i] >= 0 && j_seed[i] < mri_im->width) &&
          (i_seed[i] >= 0 && i_seed[i] < mri_im->height) &&
          (imnr_seed[i] >= 0 && imnr_seed[i] < mri_im->depth))
        MRIvox(mri_fill, j_seed[i], i_seed[i], imnr_seed[i]) = val_seed[i];
      else
        fprintf(stderr, "seed[%d] = (%d, %d, %d) is out of bounds\n",
                i, j_seed[i], i_seed[i], imnr_seed[i]) ;
    }
#else
  MRIvox(mri_fill, wm_rh_x, wm_rh_y, wm_rh_z) = RIGHT_HEMISPHERE_WHITE_MATTER ;
  MRIvox(mri_fill, wm_lh_x, wm_lh_y, wm_lh_z) = LEFT_HEMISPHERE_WHITE_MATTER ;
#endif


  /* fill in background of complement (also sometimes called the foreground) */
  fill_brain(-WM_MIN_VAL);   
  /*MRIwrite(mri_fill, "fill4.mnc") ;*/
  ilim0 = mri_im->height ; jlim0 = mri_im->width ;
  ilim1 = jlim1 = 0;
  for (k = 0 ; k < mri_fill->depth ; k++)
  {
    for (i = 0 ; i < mri_fill->height ; i++)
    {
      for (j = 0 ; j < mri_fill->width;j++)
      {
        if (MRIvox(mri_fill, j, i, k) >0)
        {
          if (i<ilim0) ilim0=i;
          if (i>ilim1) ilim1=i;
          if (j<jlim0) jlim0=j;
          if (j>jlim1) jlim1=j;
        }
      }
    }
  }
  
  if (fill_holes_flag)
    fill_holes();
  
  /*MRIwrite(mri_fill, "fill5.mnc") ;*/

  for (k = 0 ; k < mri_fill->depth ; k++)
    for (i = 0; i < mri_fill->height ; i++)
      for (j = 0; j < mri_fill->width ; j++)
      {
        if (k==0||k==mri_fill->depth-1) 
          MRIvox(mri_fill, j, i, k) = 0;
      }

  MRIwrite(mri_fill, out_fname) ;
}
static int 
fill_brain(int threshold)
{
  int dir = -1, nfilled = 10000, ntotal = 0,iter = 0;
  int im0,im1,j0,j1,i0,i1,imnr,i,j;
  int v1,v2,v3,vmax;

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
  printf("total of %d voxels filled\n",ntotal);
  return(NO_ERROR) ;
}

static int
fill_holes(void)
{
  int  nfilled, ntotal = 0, cnt, cntmax= neighbor_threshold ;
  int im0,j0,i0,imnr,i,j;
  int v,vmax;

  do
  {
    nfilled = 0;
    for (imnr=1;imnr!=mri_fill->depth-1;imnr++)
    for (i=1;i!=mri_fill->height-1;i++)
    for (j=1;j!=mri_fill->width-1;j++)
    if (MRIvox(mri_fill, j, i, imnr)==0 && 
        i>ilim0-10 && i<ilim1+10 && j>jlim0-10 && j<jlim1+10)
    {
      cnt = 0;
      vmax = 0;
      for (im0= -1;im0<=1;im0++)
      for (i0= -1;i0<=1;i0++)
      for (j0= -1;j0<=1;j0++)
      {
        v = MRIvox(mri_fill, j+j0, i+i0, imnr+im0) ;
        if (v>vmax) vmax = v;
        if (v == 0) cnt++;              /* count # of nbrs which are off */
        if (cnt>cntmax) im0=i0=j0=1;    /* break out of all 3 loops */
      }
      if (cnt<=cntmax)   /* toggle pixel (off to on, or on to off) */
      {
        if (imnr == 149 && (i == 96) && (j == 147))
          printf("filled to %d", vmax) ;
        MRIvox(mri_fill, j, i, imnr) = vmax; 
        nfilled++;
        ntotal++;
      }
    }
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "%d holes filled\n",nfilled);
  } while (nfilled > 0) ;
  
  printf("total of %d holes filled\n",ntotal);
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
  else switch (toupper(*option))
  {
  case 'L':
    wm_lh_tal_x = atof(argv[2]) ;
    wm_lh_tal_y = atof(argv[3]) ;
    wm_lh_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "using lh wm seed point (%2.0f, %2.0f, %2.0f)\n",
            wm_lh_tal_x, wm_lh_tal_y, wm_lh_tal_z) ;
    break ;
  case 'R':
    wm_rh_tal_x = atof(argv[2]) ;
    wm_rh_tal_y = atof(argv[3]) ;
    wm_rh_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "using rh wm seed point (%2.0f, %2.0f, %2.0f)\n",
            wm_rh_tal_x, wm_rh_tal_y, wm_rh_tal_z) ;
    break ;
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
  case 'F':
    min_filled = atoi(argv[2]) ;
    nargs = 1 ;
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
  fprintf(stderr, "fill version 1\n") ;
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
          "\t-L <x,y,z>     - the Talairach coords of the wm seed for the\n"
          "\t                 left hemisphere\n") ;
  fprintf(stderr, 
          "\t-R <x,y,z>     - the Talairach coords of the wm seed for the\n"
          "\t                 right hemisphere\n") ;
  fprintf(stderr, 
          "\t-P <x,y,z>     - the Talairach coords of the seed for the pons") ;
  fprintf(stderr, 
          "\t-C <x,y,z>     - the Talairach coords of the seed for the\n"
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
#define MIN_CC_AREA       350
#define MAX_CC_AREA       850
#define MIN_CC_ASPECT     0.1
#define MAX_CC_ASPECT     0.575

#define MIN_PONS_AREA     400
#define MAX_PONS_AREA     800
#define MIN_PONS_ASPECT   0.8
#define MAX_PONS_ASPECT   1.2



static MRI *
find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,int orientation,
                   int *pxv, int *pyv, int *pzv)
{
  MRI        *mri_slices[MAX_SLICES], *mri_filled[MAX_SLICES], *mri_cut ;
  Real       dx, dy, dz, x, y, z, aspect,MIN_ASPECT,MAX_ASPECT;
  int        slice, offset, area[MAX_SLICES], min_area, min_slice,xo,yo,
             xv, yv, zv, x0, y0, z0, xi, yi, zi, MIN_AREA, MAX_AREA, done ;
             
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
    break ;
  case MRI_HORIZONTAL: 
    dz = 1.0 ; dx = dy = 0.0 ; 
    name = "pons" ;
    MIN_AREA = MIN_PONS_AREA ; MAX_AREA = MAX_PONS_AREA ;
    MIN_ASPECT = MIN_PONS_ASPECT ; MAX_ASPECT = MAX_PONS_ASPECT ;
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
    sprintf(fname, "%s_seed.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_slices[0], fname) ;
    sprintf(fname, "%s_seed_fill.mnc", orientation==MRI_SAGITTAL?"cc":"pons");
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

  if (!done &&
      ((pons_seed_set && orientation == MRI_HORIZONTAL) ||
       (cc_seed_set && orientation == MRI_SAGITTAL)))
      done = 1 ;

  if (done)  /* center the seed */
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
  while (!done)
  {
    offset += SEARCH_STEP ;   /* search at a greater radius */
    if (offset >= MAX_OFFSET)
      ErrorReturn(NULL,
                  (ERROR_BADPARM, "%s: could not find valid seed for the %s",
                   Progname, name)) ;
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
            sprintf(fname, "%s_seed.mnc", 
                    orientation == MRI_SAGITTAL ? "cc":"pons");
            MRIwrite(mri_slices[0], fname) ;
            sprintf(fname, "%s_seed_fill.mnc", 
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
    fprintf(stderr, "%s seed point found at (%d, %d, %d)\n", name, xv, yv, zv);

  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    offset = slice - HALF_SLICES ;
    x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
    MRItalairachToVoxel(mri, x, y,  z, &x, &y,&z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    mri_slices[slice] = 
      MRIextractTalairachPlane(mri, NULL, orientation, xv, yv, zv, SLICE_SIZE);
    mri_filled[slice] = 
      MRIfillFG(mri_slices[slice], NULL, xo, yo,0,WM_MIN_VAL,127,&area[slice]);
    MRIboundingBox(mri_filled[slice], 1, &region) ;

    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= SLICE_SIZE-1 ||
        region.y+region.dy >= SLICE_SIZE-1)
      area[slice] = 0 ;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice[%d] @ (%d, %d, %d): area = %d\n", 
              slice, xv, yv, zv, area[slice]) ;

    if (orientation == MRI_SAGITTAL)  /* extend to top and bottom of slice */
      region.dy = SLICE_SIZE - region.y ;

    /*    for (yv = region.y ; yv < region.y+region.dy ; yv++)*/
    for (yv = region.y ; yv < region.y+region.dy ; yv++)
    {
      for (xv = region.x ; xv < region.x+region.dx ; xv++)
      {
        MRIvox(mri_filled[slice], xv, yv, 0) = 1 ;
      }
    }

    if ((Gdiag & DIAG_WRITE) && !(slice % 1))
    {
      sprintf(fname, "%s_slice%d.mnc", 
              orientation == MRI_SAGITTAL ? "cc":"pons", slice);
      MRIwrite(mri_slices[slice], fname) ;
      sprintf(fname, "%s_filled%d.mnc", 
              orientation == MRI_SAGITTAL ? "cc":"pons", slice);
      MRIwrite(mri_filled[slice], fname) ;
    }
  }

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

  if (min_slice < 0)
    ErrorReturn(NULL, 
                (ERROR_BADPARM, "%s: could not find valid seed for the %s",
                 Progname, name));
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s_slice.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_slices[min_slice], fname) ;
    sprintf(fname, "%s_filled.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_filled[min_slice], fname) ;
  }


  offset = min_slice - HALF_SLICES ;
  x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
  MRItalairachToVoxel(mri, x, y,  z, &x, &y, &z) ;
  *pxv = nint(x) ; *pyv = nint(y) ; *pzv = nint(z) ;
  mri_cut = MRIcopy(mri_filled[min_slice], NULL) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, 
            "%s: min area %d at slice %d, (%d, %d, %d)\n", 
            name, min_area, min_slice, *pxv, *pyv, *pzv) ;



  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    MRIfree(&mri_slices[slice]) ;
    MRIfree(&mri_filled[slice]) ;
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
      if (MRIvox(mri, x, y, 0))   /* find total distance to all other points */
      {
        total_dist = 0 ;
        for (y1 = 0 ; y1 < height ; y1++)
        {
          for (x1 = 0 ; x1 < width ; x1++)
          {
            if (MRIvox(mri, x1, y1, 0)) 
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
      if (MRIvox(mri, x, y, 0))  /* see if it is a border pixel */
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
  MRI  *mri_slice ;
  int  xv, yv, zv, xo, yo, max_y, thickness, y1, xcc, ycc, x, y ;
  Real xr, yr, zr ;

  MRItalairachToVoxel(mri, 0.0, 0.0, 0.0, &xr, &yr, &zr);
  xv = nint(xr) ; yv = nint(yr) ; zv = nint(zr) ;
  mri_slice = 
    MRIextractTalairachPlane(mri, NULL, MRI_CORONAL,xv,yv,zv,SLICE_SIZE) ;
  if (Gdiag & DIAG_WRITE)
  {
    MRIwrite(mri_slice, "cor.mnc") ;
  }

  xo = mri_slice->width/2 ;  yo = mri_slice->height/2 ;

  /* find the column with the lowest starting y value of any sign. thick. */
  xcc = ycc = max_y = 0 ; 
  for (x = xo-CC_SPREAD ; x <= xo+CC_SPREAD ; x++)
  {
    /* search for first non-zero pixel */
    for (y = 0 ; y < SLICE_SIZE ; y++)
    {
      if (MRIvox(mri_slice, x, y, 0))
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
  return(NO_ERROR) ;
}

#define MIN_BRAINSTEM_THICKNESS    15
#define MAX_BRAINSTEM_THICKNESS    35
#define MIN_DELTA_THICKNESS         3

static int
find_pons(MRI *mri, Real *p_ponsx, Real *p_ponsy, Real *p_ponsz)
{
  MRI   *mri_slice, *mri_filled ;
  Real  xr, yr, zr ;
  int   xv, yv, zv, x, y, width, height, thickness, xstart, xo, yo, area ;

  thickness = xstart = 0 ;  /* for compiler warning */
  MRItalairachToVoxel(mri, 0.0, 0.0, 0.0, &xr, &yr, &zr);
  xv = nint(xr) ; yv = nint(yr) ; zv = nint(zr) ;
  mri_slice = 
    MRIextractTalairachPlane(mri, NULL, MRI_SAGITTAL,xv,yv,zv,SLICE_SIZE) ;

  if (Gdiag & DIAG_WRITE)
    MRIwrite(mri_slice, "pons.mnc") ;

/*
   search from the front of the head at the bottom of a sagittal slice
   for the first significantly thick slice of white matter, which will
   be the brainstem. Then, follow the contour of the brainstem upwards
   until the first 'backwards' indentation which is about where we
   want to make the cut.
*/

  /* first find the brainstem */
  width = mri->width ; height = mri->height ;
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

  mri_filled = 
    MRIfillFG(mri_slice, NULL,xstart-(thickness-1)/2,y,0,WM_MIN_VAL,127,&area);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stderr, "found brainstem at y=%d, x = %d --> %d, area=%d\n",
            y, xstart, xstart-thickness+1, area) ;

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_filled, "pons_filled.mnc") ;

/* 
   now, starting from that slice, proceed upwards until we find a
   'cleft', which is about where we want to cut.
*/
  x = xstart ;  /* needed for initialization */
  do
  {
    xstart = x ;
    y-- ;
    for (x = width-1 ; !MRIvox(mri_filled,x,y,0) && (x >= 0) ; x--)
    {}

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice %d, xstart %d\n", y, x) ;
  } while ((xstart - x < MIN_DELTA_THICKNESS) && (y >0)) ;

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
    fprintf(stderr, "pons found at (%d, %d) --> (%d, %d, %d)\n", 
            xstart-thickness/2, y, xv, yv, zv) ;

  MRIfree(&mri_slice) ;
  MRIfree(&mri_filled) ;
  return(NO_ERROR) ;
}

