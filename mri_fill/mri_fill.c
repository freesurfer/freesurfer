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

#if 0
static int immin1,immax1,imin1,imax1,jmin1,jmax1;
static int immin2,immax2,imin2,imax2,jmin2,jmax2;
#endif
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
static Real cc_tal_x = -1.0 ;
static Real cc_tal_y = -1.0 ;
static Real cc_tal_z = 27.0 ;

static Real pons_tal_x = -2.0 ;
static Real pons_tal_y = -22.0 ;
static Real pons_tal_z = -17.0 ;

char *Progname ;

#define MAX_ITERATIONS  1000  /* was 10 */
#define MIN_FILLED      0     /* was 100 */


static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;
void main(int argc, char *argv[]) ;

/* used in fill_holes */
#define DEFAULT_NEIGHBOR_THRESHOLD    8
static int neighbor_threshold = DEFAULT_NEIGHBOR_THRESHOLD ;

static MRI *mri_fill, *mri_im ;

static int fill_holes(void) ;
static int fill_brain(float threshold,int clip_flag) ;
void main(int argc, char *argv[]) ;
static int find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,
                              int orientation) ;

void
main(int argc, char *argv[])
{
  int     i,j,k,snum ;
  int     imnr_seed[MAXSEED],i_seed[MAXSEED],j_seed[MAXSEED],val_seed[MAXSEED];
  int     nseed, nargs, wm_lh_x, wm_lh_y, wm_lh_z, wm_rh_x, wm_rh_y, wm_rh_z ;
  char    input_fname[200],out_fname[200],*data_dir;
  Real    x, y, z ;
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

  find_cutting_plane(mri_im, pons_tal_x,pons_tal_y, pons_tal_z,MRI_HORIZONTAL);
  find_cutting_plane(mri_im, cc_tal_x, cc_tal_y, cc_tal_z, MRI_SAGITTAL) ;


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
  
  fill_brain(WM_MIN_VAL,TRUE);    /* fill from 2 seeds in wm outwards */

  /*MRIwrite(mri_fill, "fill1.mnc") ;*/
  
  /* set im to initial outward fill */
  MRIcopy(mri_fill, mri_im) ;
  MRIclear(mri_fill) ;
  MRIvox(mri_fill,1,1,1) = 255;              /* seed for background */
  
  fill_brain(-WM_MIN_VAL,FALSE);/* fill in connected component of background */
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
  fill_brain(-WM_MIN_VAL,TRUE);   
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
fill_brain(float threshold,int clip_flag)
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
#if 0
                if (clip_flag&&((imnr>immin1&&imnr<immax1&&
                                 i>imin1&&i<imax1&&
                                 j>jmin1&&j<jmax1) ||
                                (imnr>immin2&&imnr<immax2&&
                                 i>imin2&&i<imax2&&
                                 j>jmin2&&j<jmax2)));
                else
#endif
                {
                  MRIvox(mri_fill, j, i, imnr) = vmax;
                  nfilled++;
                  ntotal++;
                }
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
  fprintf(stderr, "fill version -1\n") ;
  exit(0) ;
}
static void
print_help(void)
{
  fprintf(stderr, "usage: fill [options] <name> <seed>\n") ;
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

#define MAX_SLICES   15
#define MIN_AREA     100
#define SEARCH_SIZE  7
#define SHALF        ((SEARCH_SIZE-1)/2)
#define CUT_WIDTH    1
#define HWIDTH       ((CUT_WIDTH-1)/2)

static int
find_cutting_plane(MRI *mri, Real x_tal, Real y_tal,Real z_tal,int orientation)
{
  MRI        *mri_slices[MAX_SLICES], *mri_filled[MAX_SLICES] ;
  Real       dx, dy, dz, x, y, z, dist, min_dist ;
  int        slice, offset, shalf, area[MAX_SLICES], min_area, min_slice,
             xv, yv, zv, x0, y0, z0, xi, yi, zi ;
  char       fname[100] ;
  MRI_REGION region ;

  /* check to make sure seed point is in white matter */
  MRItalairachToVoxel(mri, x_tal, y_tal,  z_tal, &x, &y,&z) ;
  xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
  min_dist = SEARCH_SIZE*SEARCH_SIZE*2 ;
  if (MRIvox(mri, xv, yv, zv) < WM_MIN_VAL)   /* find a valid seed point */
  {
    for (z0 = -SHALF ; z0 <= SHALF ; z0++)
    {
      zi = mri->zi[zv+z0] ;
      for (y0 = -SHALF ; y0 <= SHALF ; y0++)
      {
        yi = mri->yi[yv+y0] ;
        for (x0 = -SHALF ; x0 <= SHALF ; x0++)
        {
          xi = mri->xi[xv+x0] ;
          dist = x0*x0 + y0*y0 + z0*z0 ;
          if ((MRIvox(mri, xi, yi, zi) >= WM_MIN_VAL) && (dist < min_dist))
          {
            /* found a new seed point */
            x = (Real)xi ; y = (Real)yi ; z = (Real)zi ;
            min_dist = dist ;
            MRIvoxelToTalairach(mri, x, y, z, &x_tal, &y_tal, &z_tal) ;
          }
        }
      }
    }
    if (min_dist > SEARCH_SIZE*SEARCH_SIZE)
      ErrorExit(ERROR_BADPARM, 
                "%s: could not find valid seed around (%d, %d, %d)",
                Progname, xv, yv, zv) ;
  }
  

  switch (orientation)
  {
  default:
  case MRI_SAGITTAL:   dx = 1.0 ; dy = dz = 0.0 ; break ;
  case MRI_HORIZONTAL: dz = 1.0 ; dx = dy = 0.0 ; break ;
  case MRI_CORONAL:    dy = 1.0 ; dx = dz = 0.0 ; break ;
  }

  shalf = (MAX_SLICES-1)/2 ;
  x0 = y0 = (SLICE_SIZE-1)/2 ;
  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    offset = slice - shalf ;
    x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
    MRItalairachToVoxel(mri, x, y,  z, &x, &y,&z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;
    mri_slices[slice] = 
      MRIextractTalairachPlane(mri, NULL, orientation, xv, yv, zv, SLICE_SIZE);
    mri_filled[slice] = 
      MRIfillFG(mri_slices[slice], NULL, x0, y0, 0, (int)WM_MIN_VAL, 127);
    MRIboundingBox(mri_filled[slice], 1, &region) ;
    area[slice] = region.dx * region.dy ;

    if (orientation == MRI_SAGITTAL)  /* extend to top and bottom of slice */
    {
#if 0
      region.dy = SLICE_SIZE ;region.y = 0 ; 
#else
      region.dy = SLICE_SIZE - region.y ;
#endif
    }

    /*    for (yv = region.y ; yv < region.y+region.dy ; yv++)*/
    for (yv = region.y ; yv < region.y+region.dy ; yv++)
    {
      for (xv = region.x ; xv < region.x+region.dx ; xv++)
      {
        MRIvox(mri_filled[slice], xv, yv, 0) = 1 ;
      }
    }

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stderr, "slice[%d]: area = %d\n", slice, area[slice]) ;

    if ((Gdiag & DIAG_WRITE) && !(slice % 5))
    {
      sprintf(fname, "slice%d.mnc", slice) ;
      MRIwrite(mri_slices[slice], fname) ;
      sprintf(fname, "filled%d.mnc", slice) ;
      MRIwrite(mri_filled[slice], fname) ;
    }
  }

  min_area = area[1] ; min_slice = 0 ;
  for (slice = 2 ; slice < MAX_SLICES-1 ; slice++)
  {
    if (area[slice] < min_area && area[slice] > MIN_AREA)
    {
      min_area = area[slice] ;
      min_slice = slice ;
    }
  }

  if (Gdiag & DIAG_WRITE)
  {
    sprintf(fname, "slice.mnc") ;
    MRIwrite(mri_slices[min_slice], fname) ;
    sprintf(fname, "filled.mnc") ;
    MRIwrite(mri_filled[min_slice], fname) ;
  }

  /* now make the cut(s) */
  for (slice = min_slice-HWIDTH ; slice <= min_slice+HWIDTH ; slice++)
  {
    offset = slice - shalf ;
    x = x_tal + dx*offset ; y = y_tal + dy*offset ; z = z_tal + dz*offset ; 
    MRItalairachToVoxel(mri, x, y,  z, &x, &y,&z) ;
    xv = nint(x) ; yv = nint(y) ; zv = nint(z) ;

    if (Gdiag & DIAG_SHOW && slice == min_slice)
      fprintf(stderr, 
              "min area %d at slice %d, coord = (%d, %d, %d)\n", 
              min_area, min_slice, xv, yv, zv) ;

    MRIeraseTalairachPlane(mri, mri_filled[slice], 
                           orientation, xv, yv, zv, SLICE_SIZE);
  }

  for (slice = 0 ; slice < MAX_SLICES ; slice++)
  {
    MRIfree(&mri_slices[slice]) ;
    MRIfree(&mri_filled[slice]) ;
  }

  return(NO_ERROR) ;
}

