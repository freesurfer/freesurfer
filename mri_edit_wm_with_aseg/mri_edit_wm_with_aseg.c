//
// mri_auto_fill.c
//
// written by Bruce Fischl
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/01/19 14:59:47 $
// Revision       : $Revision: 1.1 $
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "cma.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "const.h"
#include "transform.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label) ;
static int neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label);
static int edit_segmentation(MRI *mri_im, MRI *mri_seg) ;


char *Progname ;

static int fillven = 0 ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[])
{
	MRI    *mri_wm, *mri_aseg ;
  struct timeb  then ;
	int    msec, nargs ;

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
    usage_exit(1) ;   /* will exit */


	printf("reading wm segmentation from %s...\n", argv[1]) ;
	mri_wm = MRIread(argv[1]) ;
	if (!mri_wm)
		ErrorExit(ERROR_NOFILE, "%s: could not read wm volume from %s",
							Progname, argv[1]) ;
	mri_aseg = MRIread(argv[2]) ;
	if (!mri_aseg)
		ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s",
							Progname, argv[2]) ;


	edit_segmentation(mri_wm, mri_aseg) ;
	printf("writing edited volume to %s....\n", argv[3]) ;
	MRIwrite(mri_wm, argv[3]) ;

  msec = TimerStop(&then) ;
  fprintf(stderr, "auto filling took %2.2f minutes\n",
          (float)msec/(1000.0f*60.0f));
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
  if (!stricmp(option, "dilate"))
  {
  }
  else if (!strcmp(option, "fillven"))
  {
		fillven = atoi(argv[2]) ;
		printf("%sfilling ventricles\n", fillven ? "not " : "") ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
  {
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
static void
usage_exit(int code)
{
  printf("usage: %s <input wm volume> <aseg volume> <output wm volume> "
         "<output volume>\n", 
         Progname) ;
  exit(code) ;
}
static int
edit_segmentation(MRI *mri_wm, MRI *mri_seg)
{
  int   width, height, depth, x, y, z, label, non, noff, xi, yi, zi,  xk, yk, zk, nchanged, wsize ;
  MRI   *mri_filled ;

  mri_filled =  MRIclone(mri_wm,  NULL);

  width = mri_wm->width ; height = mri_wm->height ; depth = mri_wm->depth ;

  non = noff = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
        label = MRIvox(mri_seg, x, y, z) ;
        
        switch (label)
        {
				case Unknown:
					wsize=5 ;
					if (MRIlabelsInNbhd(mri_seg, x, y, z, (wsize-1)/2, Unknown) < (wsize*wsize*wsize-1))
						break ;
	  
					/* !!! no break - erase unknown if it is surrounded by only  unknowns */
	  
					/* erase these  labels */
        case Left_Cerebellum_White_Matter:
        case Left_Cerebellum_Exterior:
        case Left_Cerebellum_Cortex:
        case Right_Cerebellum_White_Matter:
        case Right_Cerebellum_Exterior:
        case Right_Cerebellum_Cortex:
        case Right_Cerebral_Cortex:
#if 0
					/* otherwise will never be able to find pons */
        case Brain_Stem:
        case Left_VentralDC:
        case Right_VentralDC:
        case Left_Substancia_Nigra:
        case Right_Substancia_Nigra:
#endif
          if ((neighborLabel(mri_seg, x,y,z,1,Left_Cerebral_Cortex) == 0) &&
              (neighborLabel(mri_seg, x,y,z,1,Right_Cerebral_Cortex) == 0))
          {
            if (MRIvox(mri_wm, x, y, z) >= WM_MIN_VAL)
            {
              MRIvox(mri_wm, x, y, z) = 0 ;
              noff++ ;
            }
          }
          break ;
	  
					/* fill these */
				case Left_Lesion:
				case Right_Lesion:
				case WM_hypointensities:
				case Left_WM_hypointensities:
				case Right_WM_hypointensities:
				case non_WM_hypointensities:
				case Left_non_WM_hypointensities:
				case Right_non_WM_hypointensities:
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_Cortex) >= 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_Cortex) >= 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            MRIvox(mri_wm, x, y, z) = 255 ;
            MRIvox(mri_filled, x, y, z) = 255 ;
            non++ ;  
          }
					break ;
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
					if (fillven == 0)
						break ;
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_Cortex) >= 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_Cortex) >= 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
            MRIvox(mri_wm, x, y, z) = 255 ;
            MRIvox(mri_filled, x, y, z) = 255 ;
            non++ ;  
          }
          yi = mri_wm->yi[y+1] ;
          label = MRIvox(mri_seg, x,yi, z) ;
          if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
							 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter) ||
							 (label == Unknown))
              && (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
          {
            MRIvox(mri_wm, x, yi, z) = 255 ;
            MRIvox(mri_filled, x, yi, z) = 255 ;
            non++ ;
          }
					if (label == Left_Inf_Lat_Vent || label ==  Right_Inf_Lat_Vent)  /* fill inferior wm */
					{
						int xi,  olabel ;
	    
						xi = label ==  Left_Inf_Lat_Vent ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
						olabel = MRIvox(mri_seg, xi, y, z) ;
						if (olabel != label)  /* voxel lateral to this one is not hippocampus   */
						{
							MRIvox(mri_wm, xi, y, z) = 255 ;
							MRIvox(mri_filled, xi, y, z) = 255 ;
							non++ ;
						}
	    
						yi = mri_wm->yi[y+1] ;
						label = MRIvox(mri_seg, x,yi, z) ;
						if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
								 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
								&& (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
						{
							MRIvox(mri_wm, x, yi, z) = 255 ;
							MRIvox(mri_filled, x, yi, z) = 255 ;
							yi = mri_wm->yi[y+2] ;
							MRIvox(mri_wm, x, yi, z) = 255 ;
							MRIvox(mri_filled, x, yi, z) = 255 ;
							non += 2 ;
						}
					}
          break ;
        case Left_Hippocampus:
        case Right_Hippocampus:
					{
						int xi,  olabel ;
	    
						xi = label == Right_Hippocampus ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
						olabel = MRIvox(mri_seg, xi, y, z) ;
						if (olabel != label)  /* voxel lateral to this one is not hippocampus   */
						{
							MRIvox(mri_wm, xi, y, z) = 255 ;
							MRIvox(mri_filled, xi, y, z) = 255 ;
							non++ ;
						}
	    
						yi = mri_wm->yi[y+1] ;
						label = MRIvox(mri_seg, x,yi, z) ;
						if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
								 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
								&& (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
						{
							MRIvox(mri_wm, x, yi, z) = 255 ;
							MRIvox(mri_filled, x, yi, z) = 255 ;
							yi = mri_wm->yi[y+2] ;
							MRIvox(mri_wm, x, yi, z) = 255 ;
							MRIvox(mri_filled, x, yi, z) = 255 ;
							non += 2 ;
						}
						break ;
					}
        case Left_Accumbens_area:
        case Right_Accumbens_area:
        case Left_Caudate:
        case Right_Caudate:
        case Left_Putamen:
        case Right_Putamen:
        case Left_Pallidum:
        case Right_Pallidum:
          if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
          {
            MRIvox(mri_wm, x, y, z) = 255 ;
            MRIvox(mri_filled, x, y, z) = 255 ;
            non++ ; 
          }
          break ;
        default:
          break ;
        }
      }
    }
  }

  /*
    fill in voxels that were labeled wm by the aseg, but not by  wmfilter, and are
    neighbors  of  voxels that  have been already been filled .
  */
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
				for (x = 0 ; x < width ; x++)
				{
					if (x == Gx && y == Gy && z == Gz)  
						DiagBreak() ;
					if  (MRIvox(mri_filled, x,  y, z) == 0)
						continue  ;
					for (xk = -1 ; xk <= 1 ; xk++)
					{
						xi = mri_filled->xi[x+xk] ;
						for (yk = -1 ; yk <= 1 ; yk++)
						{
							yi = mri_filled->yi[y+yk] ;
							for (zk = -1 ; zk <= 1 ; zk++)
							{
								zi = mri_filled->zi[z+zk] ;
								if (xi == Gx && yi == Gy && zi == Gz)  
									DiagBreak() ;
								label = MRIvox(mri_seg, xi, yi, zi) ;
								if (IS_WM(label) &&  (MRIvox(mri_wm, xi, yi, zi) < WM_MIN_VAL))
								{
									nchanged++ ;
									MRIvox(mri_wm, xi, yi, zi) = 255 ;
#if 0
									MRIvox(mri_filled, xi, yi, zi) = 255 ;
#endif
									non++ ;  
								}
							}
						}
					}
				}
      }
    }
    printf("%d additional wm voxels added\n", nchanged)  ;
  } while (nchanged >  0) ;
  
  printf("SEG EDIT: %d voxels turned on, %d voxels turned off.\n", non, noff) ;
	MRIfree(&mri_filled) ;
  return(NO_ERROR) ;
}

#if 0
static int
neighbors(MRI *mri, int x, int y,int z,int whalf,int label)
{
  int xi, yi, zi, xk, yk, zk, nbrs ;
  
  for (nbrs = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        if (abs(xk)+abs(yk)+abs(zk) > 1) /* only 6-connected neighbors */
          continue ;
        xi = mri->xi[x+xk] ;
        if (MRIvox(mri, xi, yi, zi) == label)
          nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
#endif

static int
neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label)
{
  int xi, yi, zi, xk, yk, zk ;
  
  for (zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        if (abs(xk)+abs(yk)+abs(zk) > 1) /* only 6-connected neighbors */
          continue ;
        xi = mri->xi[x+xk] ;
        if (MRIvox(mri, xi, yi, zi) == label)
          return(1) ;
      }
    }
  }
  return(0) ;
}

static int
MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label)
{
  int xi, yi, zi, xk, yk, zk, count ;
  
  for (count = 0, zk = -whalf ; zk <= whalf ; zk++)
  {
    zi = mri->zi[z+zk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (xk = -whalf ; xk <= whalf ; xk++)
      {
        xi = mri->xi[x+xk] ;
        if (MRIvox(mri, xi, yi, zi) == label)
	  count++;
      }
    }
  }
  return(count) ;
}
