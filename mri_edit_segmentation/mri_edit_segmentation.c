#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "timer.h"
#include "proto.h"
#include "mrinorm.h"
#include "mri_conform.h"
#include "cma.h"

static int change_label(MRI *mri_T1, MRI *mri_out_labeled, int x, int y, 
                        int z, int wsize, int left) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void  usage_exit(void) ;

static int distance_to_label(MRI *mri_labeled, int label, int x, 
                             int y, int z, int dx, int dy, 
                             int dz, int max_dist) ;
static MRI *edit_hippocampus(MRI *mri_in_labeled, MRI *mri_T1, 
                             MRI *mri_out_labeled) ;

static MRI *edit_amygdala(MRI *mri_in_labeled, MRI *mri_T1, 
                          MRI *mri_out_labeled) ;


char *Progname ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri_in_labeled, *mri_T1, *mri_out_labeled ;
  char   *in_fname, *T1_fname, *out_fname ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  in_fname = argv[1] ;
  T1_fname = argv[2] ;
  out_fname = argv[3] ;

  printf("reading from %s...\n", in_fname) ;
  mri_in_labeled = MRIread(in_fname) ;
  if (!mri_in_labeled)
    ErrorExit(ERROR_NO_FILE, "%s: could not open labeled file %s", 
              Progname, in_fname) ;

  mri_T1 = MRIread(T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_NO_FILE, "%s: could not open T1 file %s", 
              Progname, in_fname) ;

  TimerStart(&start) ;

  mri_out_labeled = edit_hippocampus(mri_in_labeled, mri_T1, NULL);
  edit_amygdala(mri_out_labeled, mri_T1, mri_out_labeled);

  MRIwrite(mri_out_labeled, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("segmentation adjustment took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
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
  if (!stricmp(option, "no1d"))
  {
    printf("disabling 1d normalization...\n") ;
  }
  else switch (toupper(*option))
  {
  case '?':
  case 'U':
    printf("usage: %s [input directory] [output directory]\n", argv[0]) ;
    exit(1) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
static void
usage_exit(void)
{
  printf("usage: %s [options] "
          "<input segmentation> <T1 volume> <output segmentation>\n", 
          Progname) ;
  exit(0) ;
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
static MRI *
edit_hippocampus(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled)
{
  int   width, height, depth, x, y, z, nchanged, dleft, label,  
        dright, dpos, dant, dup, ddown, i, left ;
  MRI   *mri_tmp ;

  nchanged = 0 ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ; height = mri_T1->height ; depth = mri_T1->depth ;

  /* change gray to hippocampus based on wm */
  for (i = 0 ; i < 3 ; i++)
  {
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 94 && y == 125 && z == 114)  
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;

          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            left = 1 ;
          case Right_Cerebral_Cortex:
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 1)
            {
              label = left ? 
                Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }

            /*
              if the current label is gray, and there is white matter below
              and hippocampus above, change to hippocampus.
            */
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2)
            {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }
            
            /*
              if the current label is gray, and there is white matter above
              and hippocampus below, change to hippocampus.
            */
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Cerebral_White_Matter :
                                    Right_Cerebral_White_Matter,
                                    x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2)
            {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,0,1,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,
                                    x,y,z,0,-1,0,2);
            if (dup <= 2 && ddown <= 2)
            {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }
            
            dleft = distance_to_label(mri_out_labeled, left ? 
                                      Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,-1,0,0,3);
            dright = distance_to_label(mri_out_labeled, left ? 
                                      Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                       x,y,z,1,0,0,3);
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Hippocampus :
                                    Right_Hippocampus,x,y,z,0,-1,0,2);
            if (dleft <= 2 && dright <= 2 && dup <= 1)
            {
              label = left ? 
                Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
            dright = distance_to_label(mri_out_labeled, left ? 
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,1,0,0,3);
            dleft = distance_to_label(mri_out_labeled, left ? 
                                      Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,-1,0,0,3);
            if (dleft <= 1 && dright <= 1)
            {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            dleft = distance_to_label(mri_out_labeled, left ? 
                                       Left_Cerebral_White_Matter :
                                       Right_Cerebral_White_Matter,
                                       x,y,z,-1,0,0,3);
            dright = distance_to_label(mri_out_labeled, left ? 
                                      Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,1,0,0,3);
            if (dleft <= 1 && dright <= 1)
            {
              change_label(mri_T1, mri_tmp, x, y, z, 9, left) ;
              nchanged++ ;
              continue ;
            }

            dright = distance_to_label(mri_out_labeled, left ? 
                                       Left_Hippocampus :
                                       Right_Hippocampus,
                                       x,y,z,1,0,0,3);
            dleft = distance_to_label(mri_out_labeled, left ? 
                                      Left_Hippocampus :
                                      Right_Hippocampus,
                                      x,y,z,-1,0,0,3);
            if (dleft <= 1 && dright <= 1)
            {
              label = left ? 
                Left_Hippocampus : Right_Hippocampus;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }

            dpos = distance_to_label(mri_out_labeled,
                                   Right_Cerebral_White_Matter,x,y,z,0,0,-1,3);
            dant = distance_to_label(mri_out_labeled,
                                   Right_Cerebral_White_Matter,x,y,z,0,0,1,3);


            break ;
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
  }

  MRIfree(&mri_tmp) ;
  printf("%d hippocampal voxels changed.\n", nchanged) ;
  return(mri_out_labeled) ;
}
static MRI *
edit_amygdala(MRI *mri_in_labeled, MRI *mri_T1, MRI *mri_out_labeled)
{
  int   width, height, depth, x, y, z, nchanged, label, total_changed, 
        dup, ddown, left ;
  MRI   *mri_tmp ;

  mri_out_labeled = MRIcopy(mri_in_labeled, mri_out_labeled) ;
  mri_tmp = MRIcopy(mri_out_labeled, NULL) ;

  width = mri_T1->width ; height = mri_T1->height ; depth = mri_T1->depth ;

  total_changed = 0 ;
  do
  {
    nchanged = 0 ;

    /* change gray to wm if near amygdala */
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 95 && y == 127 && z == 119)  
            DiagBreak() ;
          label = MRIvox(mri_tmp, x, y, z) ;
          
          left = 0 ;
          switch (label)
          {
          case Left_Cerebral_Cortex:
            left = 1 ;
          case Right_Cerebral_Cortex:
            dup = distance_to_label(mri_out_labeled,
                                    left ?  Left_Amygdala :
                                    Right_Amygdala,x,y,z,0,-1,0,2);
            ddown = distance_to_label(mri_out_labeled,
                                      left ? Left_Cerebral_White_Matter :
                                      Right_Cerebral_White_Matter,
                                      x,y,z,0,1,0,3);
            if (dup <= 1 && ddown <= 1)
            {
              label = left ? 
                Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
              MRIvox(mri_tmp, x, y, z) = label ;
              nchanged++ ;
              continue ;
            }
          default:
            break ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_out_labeled) ;
    total_changed += nchanged ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d amygdala voxels changed.\n", total_changed) ;
  return(mri_out_labeled) ;
}
float
label_mean(MRI *mri_T1, MRI *mri_labeled, int x, int y, int z, int wsize, 
           int label)
{
  int   xi, yi, zi, xk, yk, zk, whalf, nvox ;
  float mean ;


  whalf = (wsize-1)/2 ;

  for (mean = 0.0, nvox = 0, xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri_T1->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri_T1->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri_T1->zi[z+zk] ;
        if (MRIvox(mri_labeled, xi, yi, zi) != label)
          continue ;
        nvox++ ;
        mean += (float)MRIvox(mri_T1, xi, yi, zi) ;
      }
    }
  }

  if (nvox > 0)
    mean /= (float)nvox ;
  else
    mean = 0.0f ;
  return(mean) ;
}
static int 
change_label(MRI *mri_T1,MRI *mri_labeled,int x,int y,int z,int wsize,int left)
{
  float  wm_mean, hippo_mean, val ;
  int    label ;

  if (x == 95 && y == 127 && z == 119) /* dark wm (68) */
    DiagBreak() ;
  if (x == 94 && y == 126 && z == 119)  /* bright hippo (104) */
    DiagBreak() ;
  val = (float)MRIvox(mri_T1, x, y, z) ;
  wm_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize,
                       left ? Left_Cerebral_White_Matter :
                       Right_Cerebral_White_Matter) ;
  hippo_mean = label_mean(mri_T1, mri_labeled, x, y, z, wsize,
                          left ? Left_Hippocampus :
                          Right_Hippocampus) ;
  if (fabs(wm_mean-val) < fabs(hippo_mean-val))
    label = left ? 
      Left_Cerebral_White_Matter : Right_Cerebral_White_Matter;
  else
    label = left ? 
      Left_Hippocampus : Right_Hippocampus;
  
  MRIvox(mri_labeled, x, y, z) = label ;
  return(NO_ERROR) ;
}

