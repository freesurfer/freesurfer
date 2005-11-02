//
// mri_edith_wm_with_aseg.c
//
// written by Bruce Fischl
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/11/02 21:38:34 $
// Revision       : $Revision: 1.10 $
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
#include "tags.h"
#include "mrimorph.h"
#include "const.h"
#include "transform.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"

static void
DiagBreak2(){}

static int MRInonzeroInNbhd(MRI *mri, int x, int y, int z, int whalf) ;
static int MRInonfilledInNbhd(MRI *mri, int x, int y, int z, int whalf) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static int MRIlabelsInNbhd(MRI *mri, int x, int y, int z, int whalf, int label) ;
static int neighborLabel(MRI *mri, int x, int y, int z, int whalf, int label);
static int edit_segmentation(MRI *mri_im, MRI *mri_T1, MRI *mri_seg) ;
static int distance_to_label(MRI *mri_labeled, int label, int x,  
														 int y, int z, int dx, int dy, 
                             int dz, int max_dist) ;

static int distance_to_nonzero(MRI *mri_wm, int x,  
															 int y, int z, int dx, int dy, 
															 int dz, int max_dist) ;

static int distance_to_zero(MRI *mri_wm, int x,  
														int y, int z, int dx, int dy, 
														int dz, int max_dist) ;


char *Progname ;

static int fillven = 1 ;
static void usage_exit(int code) ;


int
main(int argc, char *argv[])
{
	MRI    *mri_wm, *mri_aseg, *mri_T1 ;
  struct timeb  then ;
	int    msec, nargs ;
	char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string (argc, argv, "$Id: mri_edit_wm_with_aseg.c,v 1.10 2005/11/02 21:38:34 fischl Exp $", "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_edit_wm_with_aseg.c,v 1.10 2005/11/02 21:38:34 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);

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
  
  if (argc < 5)
    usage_exit(1) ;   /* will exit */


	printf("reading wm segmentation from %s...\n", argv[1]) ;
	mri_wm = MRIread(argv[1]) ;
	if (!mri_wm)
		ErrorExit(ERROR_NOFILE, "%s: could not read wm volume from %s",
							Progname, argv[1]) ;


	MRIaddCommandLine(mri_wm, cmdline) ;
	mri_T1 = MRIread(argv[2]) ;
	if (!mri_T1)
		ErrorExit(ERROR_NOFILE, "%s: could not read T1/brain volume from %s",
							Progname, argv[2]) ;


	mri_aseg = MRIread(argv[3]) ;
	if (!mri_aseg)
		ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s",
							Progname, argv[3]) ;


	edit_segmentation(mri_wm, mri_T1, mri_aseg) ;
	printf("writing edited volume to %s....\n", argv[4]) ;
	MRIwrite(mri_wm, argv[4]) ;

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
		printf("%sfilling ventricles\n", fillven == 0 ? "not " : "") ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "debug_voxel"))
  {
		Gx = atoi(argv[2]) ;
		Gy = atoi(argv[3]) ;
		Gz = atoi(argv[4]) ;
		printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
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
  printf("usage: %s <input wm volume> <input T1/brain volume> <aseg volume> <output wm volume>\n", 
         Progname) ;
  exit(code) ;
}
static int
edit_segmentation(MRI *mri_wm, MRI *mri_T1, MRI *mri_seg)
{
  int   i, width, height, depth, x, y, z, label, non, noff, xi, yi, zi,  xk, yk, zk, nchanged, 
		    wsize, alabel, hlabel, slabel, olabel, left;
  MRI   *mri_filled ;

  mri_filled =  MRIclone(mri_wm,  NULL);

  width = mri_wm->width ; height = mri_wm->height ; depth = mri_wm->depth ;

  non = noff = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = height-2 ; y > 0 ; y--)
    {
      for (x = 1 ; x < width-1 ; x++)
      {
        if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
        label = MRIvox(mri_seg, x, y, z) ;

        left = 0 ;
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
						if (x == Gx && y == Gy && z == Gz)  
							DiagBreak2() ;
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;  
          }
					break ;
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_White_Matter) > 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_White_Matter) > 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
					{
						if (x == Gx && y == Gy && z == Gz)  
						{
							printf("filling ventricle adjacent to wm at (%d, %d, %d)\n", x,y,z) ;
							DiagBreak2() ;
						}
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;
						break ;
					}
					if (fillven == 0)
						break ;
					if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
					{
						MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
						MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
						non++ ;
					}
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
					xi = (label ==  Left_Inf_Lat_Vent) ?  mri_wm->xi[x+1] :  mri_wm->xi[x-1] ; // lateral
					olabel = MRIvox(mri_seg, xi, y, z) ;

					/* don't allow cortex to be directly lateral to inf-lat-vent - should be some wm there
					   also don't allow it to be diagonally connected */
					if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, y, z) < WM_MIN_VAL))
					{
						if (xi == Gx && y == Gy && z == Gz)
							printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", xi, y, z);
            MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
            non++ ;  
					}

					yi = mri_wm->yi[y+1] ; // inferior
					olabel = MRIvox(mri_seg, xi, yi, z) ;
					if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, yi, z) < WM_MIN_VAL))
					{
						if (xi == Gx && yi == Gy && z == Gz)
							printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", xi, yi, z);
            MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, xi, yi, z) = AUTO_FILL ;
            non++ ;  
					}

					// for spackling, don't do it if we are too close to superior wm
					if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
						continue ;

					// check diagonally anterior/posterior and inferior
					zi = mri_wm->zi[z-1] ; // posterior
					olabel = MRIvox(mri_seg, x, yi, zi) ;
					if (IS_CORTEX(olabel) && (MRIvox(mri_wm, x, yi, zi) < WM_MIN_VAL))
					{
						if (x == Gx && yi == Gy && zi == Gz)
							printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", x, yi, zi);
            MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, zi) = AUTO_FILL ;
            non++ ;  
					}
					zi = mri_wm->zi[z+1] ; // anterior
					olabel = MRIvox(mri_seg, x, yi, zi) ;
					if (IS_CORTEX(olabel) && (MRIvox(mri_wm, x, yi, zi) < WM_MIN_VAL))
					{
						if (x == Gx && yi == Gy && zi == Gz)
							printf("changing label (%d, %d, %d) to wm (gm lateral to inf-lat-vent)\n", x, yi, zi);
            MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, zi) = AUTO_FILL ;
            non++ ;  
					}

					hlabel = ((label == Left_Lateral_Ventricle) || (label == Left_Inf_Lat_Vent)) ? Left_Hippocampus : Right_Hippocampus ;
					if (distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 10) < 10)
						continue ;  // don't fill ventricular voxels superior to hippo
          if ((neighborLabel(mri_seg, x, y, z,1,Left_Cerebral_Cortex) > 0) &&
              (neighborLabel(mri_seg, x, y, z,1,Right_Cerebral_Cortex) > 0) &&
              (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL))
          {
						if (x == Gx && y == Gy && z == Gz)  
							DiagBreak2() ;
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ;  
          }
          yi = mri_wm->yi[y+1] ;
          label = MRIvox(mri_seg, x,yi, z) ;
          if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
							 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter) ||
							 (label == Unknown))
              && (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
          {
						if (x == Gx && yi == Gy && z == Gz)  
							DiagBreak2() ;
            MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
            non++ ;
          }
					if (label == Left_Inf_Lat_Vent || label ==  Right_Inf_Lat_Vent)  /* fill inferior wm */
					{
						int xi ;
	    
						xi = label ==  Left_Inf_Lat_Vent ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
						olabel = MRIvox(mri_seg, xi, y, z) ;
						/* voxel lateral to this one is not hippocampus   */
						if ((olabel != label) && (MRIvox(mri_wm, xi, y, z) < WM_MIN_VAL) && !IS_AMYGDALA(olabel) &&
								((distance_to_label(mri_seg, label ==  Left_Inf_Lat_Vent ? Left_Cerebral_White_Matter :
																		Right_Cerebral_White_Matter, xi, y, z, 0, 1, 0, 5) < 3) ||
								 (distance_to_label(mri_seg, label ==  Left_Inf_Lat_Vent ? Left_Cerebral_Cortex :
																		Right_Cerebral_Cortex, xi, y, z, 0, 1, 0, 5) < 3)) &&
								!IS_LAT_VENT(olabel))
								
						{
							if (xi == Gx && y == Gy && z == Gz)  
								DiagBreak2() ;
							MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
							MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
							non++ ;
						}
	    
						yi = mri_wm->yi[y+1] ;
						label = MRIvox(mri_seg, x,yi, z) ;
						if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
								 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
								&& (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
						{
							if (x == Gx && yi == Gy && z == Gz)  
								DiagBreak2() ;
							MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
							MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
							non++ ;
							yi = mri_wm->yi[y+2] ;
							if (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL)
							{
								if (x == Gx && yi == Gy && z == Gz)  
									DiagBreak2() ;
								MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
								MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
								non++ ;
							}
						}
					}
          break ;
        case Left_Hippocampus:
					left = 1 ;
        case Right_Hippocampus:
					{
						int xi ;

						// don't mess around with voxels near the medial or superior edge of hippocampus
						if (left)
						{
							if (distance_to_label(mri_seg, Unknown, x, y, z, -1,0,0, 10) <8)
								continue ;
						}
						else
							if (distance_to_label(mri_seg, Unknown, x, y, z, +1,0,0, 10) <8)
								continue ;
						if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
							continue ;

						// make sure we're not close to superior edge of hippo 
						olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;
						olabel = left ? Left_Thalamus_Proper : Right_Thalamus_Proper ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;
						olabel = left ? Left_VentralDC : Right_VentralDC ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;


						xi = label == Right_Hippocampus ?  mri_wm->xi[x-1] :  mri_wm->xi[x+1] ;
						yi = mri_wm->yi[y+1] ;
						olabel = MRIvox(mri_seg, xi, y, z) ;
						/* voxel lateral to this one is not hippocampus, and not
							 superior to hippocampus, and not far from wm or gm */
						if (olabel != label && (MRIvox(mri_wm, xi, y, z) < MIN_WM_VAL) &&
								(distance_to_label(mri_seg, label, xi, y, z, 0, 1, 0, 10) >= 10) &&
								((distance_to_label(mri_seg, label == Left_Hippocampus ? Left_Cerebral_Cortex : Right_Cerebral_Cortex, xi, y, z, 0, 1, 0, 5) < 3) ||
								 (distance_to_label(mri_seg, label == Left_Hippocampus ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter, xi, y, z, 0, 1, 0, 5) < 3) ||
								 (distance_to_label(mri_seg, Unknown, xi, y, z, 0, 1, 0, 5) < 3)))
								
																																									
						{
							if (xi == Gx && y == Gy && z == Gz)  
								DiagBreak2() ;
							MRIvox(mri_wm, xi, y, z) = AUTO_FILL ;
							MRIvox(mri_filled, xi, y, z) = AUTO_FILL ;
							non++ ;
						}
	    
						yi = mri_wm->yi[y+1] ;
						olabel = MRIvox(mri_seg, xi, yi, z) ;  // diagonally lateral and inferior 
						if (IS_CORTEX(olabel) && (MRIvox(mri_wm, xi, yi, z) < WM_MIN_VAL))
						{
							if (xi == Gx && yi == Gy && z == Gz)  
								DiagBreak2() ;
							MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
							MRIvox(mri_filled, xi, yi, z) = AUTO_FILL ;
							non++ ;
						}


						label = MRIvox(mri_seg, x,yi, z) ;

						if (((label == Left_Cerebral_Cortex || label == Right_Cerebral_Cortex) ||
								 (label == Left_Cerebral_White_Matter || label == Right_Cerebral_White_Matter))
								&& (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL))
						{
							if (x == Gx && yi == Gy && z == Gz)  
								DiagBreak2() ;
							MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
							MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;
							yi = mri_wm->yi[y+2] ;
							non++ ;
							if (MRIvox(mri_wm, x, yi, z) < WM_MIN_VAL)
							{
								if (x == Gx && yi == Gy && z == Gz)  
									DiagBreak2() ;
								MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
								MRIvox(mri_filled, x, yi, z) = AUTO_FILL ;

								non++ ;
							}
						}
						break ;
					}
        case Left_Accumbens_area:
        case Right_Accumbens_area:
        case Left_Caudate:
				case Left_vessel:
				case Right_vessel:
        case Right_Caudate:
        case Left_Putamen:
        case Right_Putamen:
        case Left_Pallidum:
        case Right_Pallidum:
				case Right_Thalamus_Proper:
				case Left_Thalamus_Proper:
				case Right_Thalamus:
				case Left_Thalamus:
				case Left_VentralDC:
				case Right_VentralDC:
          if (MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
          {
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ; 
          }
          break ;
				case Left_Cerebral_White_Matter:
				case Right_Cerebral_White_Matter:
					yi = mri_wm->yi[y-1] ;
					slabel = MRIvox(mri_seg, x, yi, z) ;
					if (IS_INF_LAT_VENT(slabel) && MRIvox(mri_wm, x, y, z) < WM_MIN_VAL)
					{
						if (x == Gx && y == Gy && z == Gz)
							printf("changing voxel (%d, %d, %d) to WM, due to superior inf-lat-vent\n",
										 x, y, z) ;
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ; 
					}						
					break ;
        default:
          break ;
        }
      }
    }
  }


	/* fill in the borders of the ventricle - 2mm thick. This shouldn't affect the folds
		 but will prevent small holes from ventricle into wm
	*/
	for (z = 0 ; z < depth ; z++)
	{
		for (y = 0 ; y < height ; y++)
		{
			for (x = 0 ; x < width ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Unknown:
					if (((neighborLabel(mri_seg, x, y, z, 1, Left_Lateral_Ventricle) > 0) &&
							 (neighborLabel(mri_seg, x, y, z, 1, Left_Cerebral_White_Matter) > 0)) ||
							((neighborLabel(mri_seg, x, y, z, 1, Right_Lateral_Ventricle) > 0) &&
							 (neighborLabel(mri_seg, x, y, z, 1, Right_Cerebral_White_Matter) > 0)))
					{
						if (x == Gx && y == Gy && z == Gz)
							DiagBreak2() ;
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
            non++ ; 
					}
					break ;

				case Left_Lateral_Ventricle:
					left = 1 ;
				case Right_Lateral_Ventricle:
					olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					if (neighborLabel(mri_seg, x, y, z, 2, olabel) > 0)
					{
						if (x == Gx && y == Gy && z == Gz)
							DiagBreak2() ;
            MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
            MRIvox(mri_filled, x, y, z) = AUTO_FILL ;
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
								label = MRIvox(mri_seg, xi, yi, zi) ;
								if (IS_WM(label) &&  (MRIvox(mri_wm, xi, yi, zi) < WM_MIN_VAL))
								{
									if (xi == Gx && yi == Gy && zi == Gz)  
										DiagBreak2() ;
									nchanged++ ;
									MRIvox(mri_wm, xi, yi, zi) = AUTO_FILL ;
#if 0
									MRIvox(mri_filled, xi, yi, zi) = AUTO_FILL ;
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


	// now look for voxels in which there is no wm inferior to hippo and spackle them.
	for (z = 1 ; z < depth-1 ; z++)
	{
		for (y = height-2 ; y > 0 ; y--)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
					continue ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Cerebral_Cortex:
					left = 1 ;
				case Right_Cerebral_Cortex:
				case Unknown:
					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
					xi = left ? x-1 : x+1 ;  // lateral
					olabel = left ? Left_Amygdala : Right_Amygdala ;
					if (MRIvox(mri_seg, xi, y, z) != olabel)
						continue ;
					if (distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 8) < 8)
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				case Left_Amygdala:
				case Left_Inf_Lat_Vent:
				case Left_Hippocampus:
					left = 1 ;
				case Right_Amygdala:
				case Right_Inf_Lat_Vent:
				case Right_Hippocampus:
					if (MRIvox(mri_seg, x, y+1, z) == label)  // not at inferior border
						continue ;
					if (IS_INF_LAT_VENT(label))  // check to make sure it's not all hippo inferior
					{
						olabel = MRIvox(mri_seg, x, y+1, z) ;
						if (IS_HIPPO(olabel) || IS_AMYGDALA(olabel))
							continue ;
					}
					if (IS_HIPPO(label))  // check to make sure it's not all hippo inferior
					{
						olabel = MRIvox(mri_seg, x, y+1, z) ;
						if (IS_INF_LAT_VENT(olabel))
							continue ;
					}
					if (IS_AMYGDALA(label))  // check to make sure it's not all hippo inferior
					{
						olabel = MRIvox(mri_seg, x, y+1, z) ;
						if (IS_INF_LAT_VENT(olabel) || IS_HIPPO(olabel))
							continue ;
					}

					// first make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;

					// see if there is any wm inferior to this label
					olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					if (distance_to_nonzero(mri_wm, x, y, z, 0, 1, 0, 4) <= 4)
						continue ;  // found wm inferior

					// change either this voxel or the one inferior to it to wm
					if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z))
						yi = y ;
					else 
						yi = y+1 ;
					nchanged++ ;
					MRIvox(mri_wm, x, yi, z) = AUTO_FILL ;
#if 0
					MRIvox(mri_filled, xi, yi, zi) = AUTO_FILL ;
#endif
					non++ ;  
					break ;
				}
			}
		}
	}
  
	/* more spackling. Look for hippocampal or ventricular voxels that have wm 
		 inferior, but are diagonally connected to non-wm inferior. Fill these.
	*/
	for (i = 0 ; i < 1 ; i++)
	{
		for (z = 0 ; z < depth ; z++)
		{
			for (y = height-1 ; y > 0 ; y--)
			{
				for (x = 2 ; x < width-2 ; x++)
				{
					if (x == Gx && y == Gy && z == Gz)  
						DiagBreak() ;
					if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
						continue ;
					label = MRIvox(mri_seg, x, y, z) ;
					left = 0 ;
					switch (label)
					{
					case Left_Inf_Lat_Vent:
					case Left_Hippocampus:
					case Left_Amygdala:
						left = 1 ;
					case Right_Inf_Lat_Vent:
					case Right_Hippocampus:
					case Right_Amygdala:
#if 0  // don't need because of next line
						if (MRIvox(mri_seg, x, y+1, z) == label)  // not at inferior border
							continue ;
#endif
						if (MRIvox(mri_wm, x, y+1, z) < MIN_WM_VAL)
							continue ;   // no white matter inferior

						// if on top of a thick piece of wm, this isn't needed
						if (distance_to_zero(mri_wm, x, y, z, 0, 1, 0, 4) > 2)
							continue ;

						// only if we are close to cortex or unknowns below
						olabel = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
						if ((distance_to_label(mri_seg, Unknown, x, y, z, 0, 1, 0, 10) > 5) &&
								(distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 10) > 5))
							continue ;

						// but not if there is wm close above (can cause defects!)
						if (distance_to_nonzero(mri_wm, x, y, z, 0, -1, 0, 5) <= 3)
							continue ;
						olabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;
						olabel = left ? Left_Thalamus_Proper : Right_Thalamus_Proper ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;
						olabel = left ? Left_VentralDC : Right_VentralDC ;
						if (distance_to_label(mri_seg, olabel, x, y, z, 0, -1, 0, 5) < 3)
							continue ;

#if 0  // don't need to do this check, because wm is inferior
						if (IS_INF_LAT_VENT(label))  // check to make sure it's not all hippo inferior
						{
							olabel = MRIvox(mri_seg, x, y+1, z) ;
							if (IS_HIPPO(olabel) || IS_AMYGDALA(olabel))
								continue ;
						}
						if (IS_HIPPO(label))  // check to make sure it's not all ventricle inferior
						{
							olabel = MRIvox(mri_seg, x, y+1, z) ;
							if (IS_INF_LAT_VENT(olabel))
								continue ;
						}
						if (IS_AMYGDALA(label))  // check to make sure it's not all hippo inferior
						{
							olabel = MRIvox(mri_seg, x, y+1, z) ;
							if (IS_INF_LAT_VENT(olabel) || IS_HIPPO(olabel))
								continue ;
						}
#endif

						// make sure we aren't on the medial edge of hippo
						if (left)
						{
							if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
								continue ;
						}
						else
							if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
								continue ;

						if (MRIvox(mri_wm, x, y+1, z-1) < MIN_WM_VAL) // fill one of these
						{
							if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z-1))
							{
								yi = y ;
								zi = z ;
							}
							else 
							{
								yi = y+1 ;
								zi = z-1 ;
							}
							if (x == Gx && yi == Gy && zi == Gz)
								DiagBreak2() ;
							MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
							nchanged++ ;
							non++ ;  
						}

						if (MRIvox(mri_wm, x, y+1, z+1) < MIN_WM_VAL) // fill one of these
						{
							if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x, y+1, z+1))
							{
								yi = y ;
								zi = z ;
							}
							else 
							{
								yi = y+1 ;
								zi = z+1 ;
							}
							if (x == Gx && yi == Gy && zi == Gz)
								DiagBreak2() ;
							MRIvox(mri_wm, x, yi, zi) = AUTO_FILL ;
							nchanged++ ;
							non++ ;  
						}

						if (MRIvox(mri_wm, x-1, y+1, z) < MIN_WM_VAL) // fill one of these
						{
							if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x-1, y+1, z))
							{
								xi = x ;
								yi = y ;
							}
							else 
							{
								xi = x-1 ;
								yi = y+1 ;
							}
							if (xi == Gx && yi == Gy && z == Gz)
								DiagBreak2() ;
							MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
							nchanged++ ;
							non++ ;  
						}

						if (MRIvox(mri_wm, x+1, y+1, z) < MIN_WM_VAL) // fill one of these
						{
							if (MRIvox(mri_T1, x, y, z) > MRIvox(mri_T1, x+1, y+1, z))
							{
								xi = x ;
								yi = y ;
							}
							else 
							{
								xi = x-1 ;
								yi = y+1 ;
							}
							if (xi == Gx && yi == Gy && z == Gz)
								DiagBreak2() ;
							MRIvox(mri_wm, xi, yi, z) = AUTO_FILL ;
							nchanged++ ;
							non++ ;  
						}

						break ;
					}
				}
			}
		}
	}

	// spackle the amygdala
	for (z = 0 ; z < depth ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Cerebral_Cortex:
					left = 1 ;
				case Right_Cerebral_Cortex:
					if (MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL)
						continue ;
					// make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;

					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
					xi = left ? x-1 : x+1 ;  // lateral
					olabel = left ? Left_Amygdala : Right_Amygdala ;
					if (MRIvox(mri_seg, xi, y, z) != olabel)
						continue ;
					if (distance_to_label(mri_seg, olabel, x, y, z, 0, 1, 0, 8) < 8)
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				}
			}
		}
	}

	// spackle diagonal connectivity topology flaws
	for (z = 1 ; z < depth-1 ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Unknown:
					if (MRIvox(mri_wm, x, y-1,z) >= MIN_WM_VAL)
						continue ;
					if (MRIvox(mri_wm, x, y,z) >= MIN_WM_VAL)
						continue ;
					if ((MRIvox(mri_wm, x, y-1, z-1) >= MIN_WM_VAL) && (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL))
						continue ;
					if (!IS_HIPPO(MRIvox(mri_seg, x, y-1, z-1)) && !IS_AMYGDALA(MRIvox(mri_seg, x, y-1, z-1)) &&
							!IS_INF_LAT_VENT(MRIvox(mri_seg, x, y-1, z-1)) &&
							!IS_HIPPO(MRIvox(mri_seg, x, y-1, z+1)) && !IS_AMYGDALA(MRIvox(mri_seg, x, y-1, z+1)) &&
							!IS_INF_LAT_VENT(MRIvox(mri_seg, x, y-1, z+1)))
						continue ;

					if (distance_to_label(mri_seg, Left_Inf_Lat_Vent, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (distance_to_label(mri_seg, Right_Inf_Lat_Vent, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (distance_to_label(mri_seg, Right_Hippocampus, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (distance_to_label(mri_seg, Right_Amygdala, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (distance_to_label(mri_seg, Left_Hippocampus, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (distance_to_label(mri_seg, Left_Amygdala, x, y, z, 0, 1, 0, 6) < 5)
						continue ;
					if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 10)  // only if there is enough wm in vicinity
						continue ;
					if (MRIvox(mri_wm, x, y-1, z-1) < MIN_WM_VAL)
					{
						if (x == Gx && y-1 == Gy && z-1 == Gz)
							DiagBreak() ;
						non++ ;
						nchanged++ ;
						MRIvox(mri_wm, x, y-1, z-1) = AUTO_FILL ;
					}
					if (MRIvox(mri_wm, x, y-1, z+1) < MIN_WM_VAL)
					{
						if (x == Gx && y-1 == Gy && z+1 == Gz)
							DiagBreak() ;
						non++ ;
						nchanged++ ;
						MRIvox(mri_wm, x, y-1, z+1) = AUTO_FILL ;
					}
					break ;
				case Left_Cerebral_Cortex:
					left = 1 ;
				case Right_Cerebral_Cortex:
					// make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;
					slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
#endif
					xi = left ? x-1 : x+1 ;  // lateral
					if (MRIvox(mri_wm, xi, y-1, z) >= MIN_WM_VAL)
						continue ;   // only if diagonal voxel isn't on
					alabel = left ? Left_Amygdala : Right_Amygdala ;
					hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
					if ((MRIvox(mri_seg, xi, y-1, z) != alabel) && (MRIvox(mri_seg, xi, y-1, z) != hlabel))
						continue ;
					// check to make sure we are at inferior border of hippo or amygdala
					if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
							(distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
						continue ;
					if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				}
			}
		}
	}
	for (z = 1 ; z < depth-1 ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Cerebral_Cortex:
					left = 1 ;
				case Right_Cerebral_Cortex:
					// make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;
					slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
#endif
					if (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL)
						continue ;   // only if diagonal voxel isn't on
					alabel = left ? Left_Amygdala : Right_Amygdala ;
					hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
					if ((MRIvox(mri_seg, x, y-1, z+1) != alabel) && (MRIvox(mri_seg, x, y-1, z+1) != hlabel))
						continue ;
					// check to make sure we are at inferior border of hippo or amygdala
					if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
							(distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
						continue ;
					if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				}
			}
		}
	}
	for (z = 1 ; z < depth-1 ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Cerebral_Cortex:
					left = 1 ;
				case Right_Cerebral_Cortex:
					// make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;
					slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
#endif
					if (MRIvox(mri_wm, x, y-1, z-1) >= MIN_WM_VAL)
						continue ;   // only if diagonal voxel isn't on
					alabel = left ? Left_Amygdala : Right_Amygdala ;
					hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
					if ((MRIvox(mri_seg, x, y-1, z-1) != alabel) && (MRIvox(mri_seg, x, y-1, z-1) != hlabel))
						continue ;
					// check to make sure we are at inferior border of hippo or amygdala
					if ((distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8) ||
							(distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8))
						continue ;
					if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				}
			}
		}
	}


	// look for unknown voxels with hippo diagonal connectivity.
	for (z = 1 ; z < depth-1 ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Unknown:
					if (MRIvox(mri_seg, x-1, y-1, z) == Left_Hippocampus)
					{
						xi = x-1 ;
						left = 1 ;
						hlabel = Left_Hippocampus ;
					}
					else if (MRIvox(mri_seg, x+1, y-1, z) != Right_Hippocampus)
						continue ;
					else
					{
						xi = x+1 ;
						hlabel = Right_Hippocampus ;
					}
					if (left)
					{
						i =  distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10)  ;
						if (i < 10 && MRIvox(mri_wm, x-i, y, z) < WM_MIN_VAL)
							continue ;
					}
					else
					{
						i = distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) ;
						if (i < 10 && MRIvox(mri_wm, x+i, y, z) < WM_MIN_VAL)
							continue ;
					}

					if ((distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))  // no hippo inferior
						continue ;
					if (distance_to_nonzero(mri_wm, x, y, z, 0, 1, 0, 8) < 7)
						continue ;
					non++ ;
					nchanged++ ;
					if (xi == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, xi, y-1, z) = AUTO_FILL ;
					break ;
					
				case Right_Cerebral_Cortex:
					// make sure we aren't on the medial edge of hippo
					if (left)
					{
						if (distance_to_label(mri_seg, Unknown, x, y, z, -1, 0, 0, 10) < 10)
							continue ;
					}
					else
						if (distance_to_label(mri_seg, Unknown, x, y, z, 1, 0, 0, 10) < 10)
							continue ;
					slabel = left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
					// look for voxels that are lateral to amygdala, and inf to wm. Should be filled
#if 0
					if (MRIvox(mri_wm, x, y-1, z) < MIN_WM_VAL)
						continue ;
#endif
					if (MRIvox(mri_wm, x, y-1, z+1) >= MIN_WM_VAL)
						continue ;   // only if diagonal voxel isn't on
					alabel = left ? Left_Amygdala : Right_Amygdala ;
					hlabel = left ? Left_Hippocampus : Right_Hippocampus ;
					if ((MRIvox(mri_seg, x, y-1, z+1) != alabel) && (MRIvox(mri_seg, x, y-1, z+1) != hlabel))
						continue ;
					// check to make sure we are at inferior border of hippo or amygdala
					if ((distance_to_label(mri_seg, alabel, x, y, z, 0, 1, 0, 8) < 8) ||
							(distance_to_label(mri_seg, hlabel, x, y, z, 0, 1, 0, 8) < 8))
						continue ;
					if (MRInonfilledInNbhd(mri_wm, x, y, z, 1) < 2)  // only if there is enough wm in vicinity
						continue ;
					non++ ;
					if (x == Gx && y == Gy && z == Gz)
						DiagBreak2() ;
					MRIvox(mri_wm, x, y, z) = AUTO_FILL ;
					nchanged++ ;
					break ;
				}
			}
		}
	}

  // remove uneeded filled voxels by checking for ones that are hippo/amy/inf lat and have lots of wm inferior
	for (z = 1 ; z < depth ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				if (MRIvox(mri_wm, x, y, z) != AUTO_FILL)
					continue ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Amygdala:
				case Left_Inf_Lat_Vent:
				case Left_Hippocampus:
					left = 1 ;
				case Right_Inf_Lat_Vent:
				case Right_Hippocampus:
				case Right_Amygdala:
					{
						int erase = 1 ;
						
						// if either of the 2 inferior planes have all 9 voxels on, erase this one
						for (yi = y+1 ; yi<=y+2 ; yi++)
						{
							erase = 1 ;
							for (xi = x-1 ; erase && xi <= x+1 ; xi++)
							{
								for (zi = z-1 ; erase && zi <= z+1 ; zi++)
								{
									if (MRIvox(mri_wm, xi, yi, zi) < MIN_WM_VAL)
									{
										erase = 0 ;
										break ;
									}
								}
							}
							if (erase)  // just need one of the two planes to be full
								break ;
						}
						if (erase)
						{
							if (x == Gx && y == Gy && z == Gz)
								DiagBreak2() ;
              MRIvox(mri_wm, x, y, z) = 0 ;
              noff++ ;
						}
					}
					break ;
				}
			}
		}
	}

  // add voxels that are wm in the aseg, but not in the wm vol
	for (z = 1 ; z < depth ; z++)
	{
		for (y = 1 ; y < height-1 ; y++)
		{
			for (x = 2 ; x < width-2 ; x++)
			{
				int done ;

				if (x == Gx && y == Gy && z == Gz)  
					DiagBreak() ;
				if ((MRIvox(mri_wm, x, y, z) >= MIN_WM_VAL) ||
						(!IS_WHITE_CLASS(MRIvox(mri_seg, x, y, z))))
					continue ;
				label = MRIvox(mri_seg, x, y, z) ;
				left = 0 ;
				switch (label)
				{
				case Left_Cerebral_White_Matter:
					left = 1 ;
				case Right_Cerebral_White_Matter:
					hlabel = left ?  Left_Hippocampus : Right_Hippocampus ;
					// if there is any hippo superior to this label, turn it on

					yi = y-1 ;
					done = 0 ;
					for (xi = x-1 ; !done && xi <= x+1 ; xi++)
					{
						for (zi = z-1 ; !done && zi <= z+1 ; zi++)
						{
							if (MRIvox(mri_seg, xi, yi, zi) == hlabel)
							{
								done = 1 ;
								if (x == Gx && y == Gy && z == Gz)
									DiagBreak2() ;
								MRIvox(mri_wm, x,y,z) = AUTO_FILL ;
								non++ ; 
								nchanged++ ;
								break ;
							}
						}
					}
					break ;
				}
			}
		}
	}
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
#if 0
        if (abs(xk)+abs(yk)+abs(zk) > 1) /* only 6-connected neighbors */
          continue ;
#endif
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
static int
MRInonzeroInNbhd(MRI *mri, int x, int y, int z, int whalf)
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
        if (MRIvox(mri, xi, yi, zi) >= MIN_WM_VAL)
					count++;
      }
    }
  }
  return(count) ;
}
static int
distance_to_label(MRI *mri_labeled, int label, int x, int y, int z, int dx, 
                  int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ; yi = y + d * dy ; zi = z + d * dz ;
    xi = mri_labeled->xi[xi] ; 
    yi = mri_labeled->yi[yi] ; 
    zi = mri_labeled->zi[zi];
    if (MRIvox(mri_labeled, xi, yi, zi) == label)
      break ;
  }

  return(d) ;
}

static int
distance_to_nonzero(MRI *mri_wm, int x, int y, int z, int dx, 
										int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ; yi = y + d * dy ; zi = z + d * dz ;
    xi = mri_wm->xi[xi] ; 
    yi = mri_wm->yi[yi] ; 
    zi = mri_wm->zi[zi];
    if (MRIvox(mri_wm, xi, yi, zi) >= MIN_WM_VAL)
      break ;
  }

  return(d) ;
}
static int
distance_to_zero(MRI *mri_wm, int x, int y, int z, int dx, 
										int dy, int dz, int max_dist)
{
  int   xi, yi, zi, d ;

  for (d = 1 ; d <= max_dist ; d++)
  {
    xi = x + d * dx ; yi = y + d * dy ; zi = z + d * dz ;
    xi = mri_wm->xi[xi] ; 
    yi = mri_wm->yi[yi] ; 
    zi = mri_wm->zi[zi];
    if (MRIvox(mri_wm, xi, yi, zi) < MIN_WM_VAL)
      break ;
  }

  return(d) ;
}

static int
MRInonfilledInNbhd(MRI *mri, int x, int y, int z, int whalf)
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
        if ((MRIvox(mri, xi, yi, zi) >= MIN_WM_VAL) && (MRIvox(mri, xi, yi, zi) != AUTO_FILL))
					count++;
      }
    }
  }
  return(count) ;
}
