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

#define LEFT_HEMISPHERE_WHITE_MATTER    127
#define RIGHT_HEMISPHERE_WHITE_MATTER   255

#define MAXSEED 20
#define IMAGE_DIR "/usr3/people/dale/MRI/IMAGES"

#define SLICE_SIZE 128

#define WM_MIN_VAL    2.0 /* anything below this is not white matter */

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
#if 1
static Real cc_tal_x = -1.0 ;
static Real cc_tal_y = -1.0 ;
static Real cc_tal_z = 27.0 ;

static Real pons_tal_x = -2.0 ;
static Real pons_tal_y = -22.0 ;
static Real pons_tal_z = -17.0 ;
#else

static Real cc_tal_x = -4.0 ;
static Real cc_tal_y = -32.0 ;
static Real cc_tal_z = 27.0 ;

static Real pons_tal_x = -10.0 ;
static Real pons_tal_y = -36.0 ;
static Real pons_tal_z = -20.0 ;

#endif

char *Progname ;

#define MAX_ITERATIONS  1000  /* was 10 */
#define MIN_FILLED      100     /* was 100 */


static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;
void main(int argc, char *argv[]) ;

/* used in fill_holes */
#define DEFAULT_NEIGHBOR_THRESHOLD    8
static int neighbor_threshold = DEFAULT_NEIGHBOR_THRESHOLD ;

static MRI *mri_fill, *mri_im ;

static int fill_holes(void) ;
static int fill_brain(float threshold) ;
void main(int argc, char *argv[]) ;
static MRI *find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,
                              int orientation, int *pxv, int *pvy, int *pzv) ;

void
main(int argc, char *argv[])
{
  int     i,j,k,snum ;
  int     imnr_seed[MAXSEED],i_seed[MAXSEED],j_seed[MAXSEED],val_seed[MAXSEED];
  int     nseed, nargs, wm_lh_x, wm_lh_y, wm_lh_z, wm_rh_x, wm_rh_y, wm_rh_z ;
  char    input_fname[200],out_fname[200],*data_dir;
  Real    x, y, z ;
  MRI     *mri_cc, *mri_pons ;
  int     x_pons, y_pons, z_pons, x_cc, y_cc, z_cc ;
#if 0
  int     timin ;
  char    fname[200];
  FILE    *fp;
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
  
  if (argc < 2)
  {
    printf("Usage: %s <name> <seed_number>\n",argv[0]);
    exit(0);
  }
#if 0
  if (argc > 2)
    snum = atoi(argv[2]);
  else
#endif
    snum = 0 ;
  
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
    ErrorExit(ERROR_NOFILE, "%s: could not read %s", input_fname) ;
  mri_fill = MRIclone(mri_im, NULL) ;

  /* 
     note: must do CC first to avoid disconnecting brainstem. Otherwise,
     if the seed point winds up in the brainstem it may be mistaken
     for the CC (it will be an in-plane-disconnected structure in
     about the right place, with about the right area.
     */
  mri_cc = 
    find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z, MRI_SAGITTAL,
                       &x_cc, &y_cc, &z_cc) ;
  mri_pons = 
    find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,
                       MRI_HORIZONTAL, &x_pons, &y_pons, &z_pons);

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

#if 0
  sprintf(fname,"fill.dat");
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("file %s not found.\n",fname);
    exit(0);
  }
  
  fscanf(fp,"%*s %d",&nseed);
  for (i=0;i<nseed;i++)
  {
    fscanf(fp,"%*s %d %d %d %d",
           &imnr_seed[i],&i_seed[i],&j_seed[i],&val_seed[i]);
    imnr_seed[i] += mri_im->zstart ;
    j_seed[i] += mri_im->xstart ;
    i_seed[i] = mri_im->yend-i_seed[i]+1;
  }

  fscanf(fp,"%*s %d %d",&immin1,&immax1);
  fscanf(fp,"%*s %d %d",&imin1,&imax1);
  fscanf(fp,"%*s %d %d",&jmin1,&jmax1);
  fscanf(fp,"%*s %d %d",&immin2,&immax2);
  fscanf(fp,"%*s %d %d",&imin2,&imax2);
  fscanf(fp,"%*s %d %d",&jmin2,&jmax2);
  timin = imin1;
  imin1 = 255-imax1+mri_im->ystart;
  imax1 = 255-timin+mri_im->ystart;
  timin = imin2;
  imin2 = 255-imax2+mri_im->ystart ;
  imax2 = 255-timin+mri_im->ystart;
  jmin1 += mri_im->xstart ; jmin2 += mri_im->xstart ;
  jmax1 += mri_im->xstart ; jmax2 += mri_im->xstart ;
  immin1 += mri_im->zstart ; immin2 += mri_im->zstart ;
  immax1 += mri_im->zstart ; immax2 += mri_im->zstart ;
#endif

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

  MRItalairachToVoxel(mri_im, wm_lh_tal_x, wm_lh_tal_y, wm_lh_tal_z,&x, &y,&z);
  wm_lh_x = nint(x) ; wm_lh_y = nint(y) ; wm_lh_z = nint(z) ;
  MRItalairachToVoxel(mri_im, wm_rh_tal_x, wm_rh_tal_y, wm_rh_tal_z,&x, &y,&z);
  wm_rh_x = nint(x) ; wm_rh_y = nint(y) ; wm_rh_z = nint(z) ;
  MRIvox(mri_fill, wm_rh_x, wm_rh_y, wm_rh_z) = RIGHT_HEMISPHERE_WHITE_MATTER ;
  MRIvox(mri_fill, wm_lh_x, wm_lh_y, wm_lh_z) = LEFT_HEMISPHERE_WHITE_MATTER ;

  i_seed[0] = wm_rh_y ; j_seed[0] = wm_rh_x ; imnr_seed[0] = wm_rh_z ;
  val_seed[0] = RIGHT_HEMISPHERE_WHITE_MATTER ;
  i_seed[1] = wm_lh_y ; j_seed[1] = wm_lh_x ; imnr_seed[1] = wm_lh_z ;
  val_seed[1] = LEFT_HEMISPHERE_WHITE_MATTER ;
  nseed = 2 ;

#if 0
  for (i=0;i<nseed;i++)
    if (snum==0 || snum==i+1)
    {
      MRIvox(mri_fill, j_seed[i], i_seed[i], imnr_seed[i]) = val_seed[i];
    }
#endif
  
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
fill_brain(float threshold)
{
  int dir = -1, nfilled = 10000, ntotal = 0,iter = 0;
  int im0,im1,j0,j1,i0,i1,imnr,i,j;
  int v1,v2,v3,vmax;

  while (nfilled>MIN_FILLED && iter<MAX_ITERATIONS)
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
    printf("%d voxels filled\n",nfilled);
  } 
  printf("total of %d voxels filled\n",ntotal);
  return(NO_ERROR) ;
}

static int
fill_holes(void)
{
  int nfilled = 10000, ntotal = 0, cnt, cntmax= neighbor_threshold ;
  int im0,j0,i0,imnr,i,j;
  int v,vmax;

  while (nfilled>MIN_FILLED)
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
    printf("%d holes filled\n",nfilled);
  } 
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
  case 'T':
    if (sscanf(argv[2], "%d", &neighbor_threshold) < 1)
    {
      fprintf(stderr, "fill: could not scan threshold from '%s'", argv[2]) ;
      exit(1) ;
    }
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
          "\t-T <threshold>  - specify fill_holes threshold (default=%d)\n",
          DEFAULT_NEIGHBOR_THRESHOLD) ;
  exit(0) ;
}

/* build a set of slices in Talairach space, and find the one in which 
   the central connected component has the smallest cross-sectional
   area. Then use this as a mask for erasing ('cutting') stuff out
   of the input image.
*/

#define MAX_SLICES     15
#define HALF_SLICES    ((MAX_SLICES-1)/2)
#define CUT_WIDTH      1
#define HALF_CUT      ((CUT_WIDTH-1)/2)
#define MIN_CC_AREA    1500
#define MAX_CC_AREA    2500
#define MIN_PONS_AREA  500
#define MAX_PONS_AREA  1500
#define SEARCH_STEP    5
#define MAX_OFFSET     50

static MRI *
find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,int orientation,
                   int *pxv, int *pyv, int *pzv)
{
  MRI        *mri_slices[MAX_SLICES], *mri_filled[MAX_SLICES], *mri_cut ;
  Real       dx, dy, dz, x, y, z, dist, min_dist ;
  int        slice, offset, area[MAX_SLICES], min_area, min_slice,xo,yo,
             xv, yv, zv, x0, y0, z0, xi, yi, zi, MIN_AREA, MAX_AREA, done,
             x1, y1 ;
  char       fname[100], *name ;
  MRI_REGION region ;

  switch (orientation)
  {
  default:
  case MRI_SAGITTAL:   
    dx = 1.0 ; dy = dz = 0.0 ; 
    name = "corpus callosum" ;
    MIN_AREA = MIN_CC_AREA ; MAX_AREA = MAX_CC_AREA ;
    break ;
  case MRI_HORIZONTAL: 
    dz = 1.0 ; dx = dy = 0.0 ; 
    name = "pons" ;
    MIN_AREA = MIN_PONS_AREA ; MAX_AREA = MAX_PONS_AREA ;
    break ;
  case MRI_CORONAL:    
    dy = 1.0 ; dx = dz = 0.0 ; 
    MIN_AREA = MIN_CC_AREA ; MAX_AREA = MAX_CC_AREA ;
    name = "coronal" ;
    break ;
  }

  xo = yo = (SLICE_SIZE-1)/2 ;  /* center point of the slice */

  /* search for valid seed point */
  MRItalairachToVoxel(mri, x_tal, y_tal,  z_tal, &x, &y, &z) ;
  xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
  mri_slices[0] = 
    MRIextractTalairachPlane(mri, NULL, orientation, xv, yv, zv, SLICE_SIZE);
  mri_filled[0] = MRIfillFG(mri_slices[0],NULL,xo,yo, 0, (int)WM_MIN_VAL,127);
  MRIboundingBox(mri_filled[0], 1, &region) ;
  area[0] = region.dx * region.dy ;

#if 0
#undef DIAG_VERBOSE_ON
#define DIAG_VERBOSE_ON  1
#endif
  
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    sprintf(fname, "%s_slice.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_slices[0], fname) ;
    sprintf(fname, "%s_filled.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_filled[0], fname) ;
  }
  
  /* now check to see if it could be a valid seed point based on:
     1) bound on the area
     2) the connected component is completely contained in the slice.
     */
  done =  
    ((area[0] >= MIN_AREA) &&
     (area[0] <= MAX_AREA) &&
     (region.y > 0) &&
     (region.x > 0) &&
     (region.x+region.dx < SLICE_SIZE-1) &&
     (region.y+region.dy < SLICE_SIZE-1)) ;

  if (done)  /* center the seed */
  {
    x0 = region.x + region.dx / 2 ;
    y0 = region.y + region.dy / 2 ;
    min_dist = 10000 ; xi = yi = 0 ;
    for (y1 = 0 ; y1 < SLICE_SIZE ; y1++)
    {
      for (x1 = 0 ; x1 < SLICE_SIZE ; x1++)
      {
        if (MRIvox(mri_filled[0], x1, y1, 0))
        {
          dist = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) ;
          if (dist < min_dist)
          {
            min_dist = dist ;
            xi = x1 ; yi = y1 ;
          }
        }
      }
    }
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
      ErrorExit(ERROR_BADPARM, "%s: could not find valid seed for the %s",
                Progname, name) ;
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
            MRIfillFG(mri_slices[0], NULL, xo, yo, 0, (int)WM_MIN_VAL,127);
          MRIboundingBox(mri_filled[0], 1, &region) ;
          area[0] = region.dx * region.dy ;

          if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
          {
            sprintf(fname, "%s_slice.mnc", 
                    orientation == MRI_SAGITTAL ? "cc":"pons");
            MRIwrite(mri_slices[0], fname) ;
            sprintf(fname, "%s_filled.mnc", 
                    orientation == MRI_SAGITTAL ? "cc":"pons");
            MRIwrite(mri_filled[0], fname) ;
          }
    
          /* now check to see if it could be a valid seed point based on:
             1) bound on the area
             2) the connected component is completely contained in the slice.
             */
          if ((area[0] >= MIN_AREA) &&
              (area[0] <= MAX_AREA) &&
              (region.y > 0) &&
              (region.x > 0) &&
              (region.x+region.dx < SLICE_SIZE-1) &&
              (region.y+region.dy < SLICE_SIZE-1))
          {
            /* center the seed */
            x0 = region.x + region.dx / 2 ;
            y0 = region.y + region.dy / 2 ;
            min_dist = 10000 ;
            for (y1 = 0 ; y1 < SLICE_SIZE ; y1++)
            {
              for (x1 = 0 ; x1 < SLICE_SIZE ; x1++)
              {
                if (MRIvox(mri_filled[0], x1, y1, 0))
                {
                  dist = (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) ;
                  if (dist < min_dist)
                  {
                    min_dist = dist ;
                    xv = x1 ; yv = y1 ;
                  }
                }
              }
            }
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
      MRIfillFG(mri_slices[slice], NULL, xo, yo, 0, (int)WM_MIN_VAL, 127);
    MRIboundingBox(mri_filled[slice], 1, &region) ;

    /* don't trust slices that extend to the border of the image */
    if (!region.x || !region.y || region.x+region.dx >= SLICE_SIZE ||
        region.y+region.dy >= SLICE_SIZE)
      area[slice] = 0 ;
    else
      area[slice] = region.dx * region.dy ;

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

    if ((Gdiag & DIAG_WRITE) && !(slice % 5))
    {
      sprintf(fname, "%s_slice.mnc", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_slices[slice], fname) ;
      sprintf(fname, "%s_filled.mnc", 
              orientation == MRI_SAGITTAL ? "cc":"pons");
      MRIwrite(mri_filled[slice], fname) ;
    }
  }

  min_area = 10000 ; min_slice = -1 ;
  for (slice = 1 ; slice < MAX_SLICES-1 ; slice++)
  {
    if (area[slice] < min_area && 
        area[slice] >= MIN_AREA && area[slice] <= MAX_AREA)
    {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }

  if (min_slice < 0)
    ErrorExit(ERROR_BADPARM, "%s: could not find valid seed for the %s",
              Progname, name);
  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "%s_slice.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_slices[min_slice], fname) ;
    sprintf(fname, "%s_filled.mnc", orientation == MRI_SAGITTAL ? "cc":"pons");
    MRIwrite(mri_filled[min_slice], fname) ;
  }


  offset = min_slice - HALF_SLICES ;
  x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
  MRItalairachToVoxel(mri, x_tal, y,  z, &x, &y, &z) ;
  *pxv = nint(x) ; *pyv = nint(y) ; *pzv = nint(z) ;
  mri_cut = MRIcopy(mri_filled[min_slice], NULL) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, 
            "%s: min area %d at slice %d, coord = (%d, %d, %d)\n", 
            name, min_area, min_slice, *pxv, *pyv, *pzv) ;



  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    MRIfree(&mri_slices[slice]) ;
    MRIfree(&mri_filled[slice]) ;
  }

  return(mri_cut) ;
}

