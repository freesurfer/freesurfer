/**
 * @file  mri_simulate_atrophy.c
 * @brief fuse a set of segmentations (asegs)
 *
 * program to fuse a group of cross sectional segmentations into 
 * an initial estimate of a longitudinal one
 * See Sabuncu et al., MICCA 2009 (SNIP paper).
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.5 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "transform.h"


int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;
MRI *MRIsimulateAtrophy(MRI *mri_norm, MRI *mri_aseg,  int target_label, 
                        int *border_labels, int nlabels, float atrophy_pct, 
                        MRI *mri_norm_atrophy)  ;

static float noise_sigma = 4 ;
static float atrophy_pct = 0.05 ;
static int nlabels = 1 ;
#define MAX_LABELS 100
static int target_labels[MAX_LABELS] = { Left_Hippocampus } ;
static int border_labels[MAX_LABELS][MAX_LABELS] ;

int
main(int argc, char *argv[])
{
  char        **av, *out_fname ;
  int          ac, nargs, msec, minutes, seconds ;
  struct timeb start ;
  MRI         *mri_aseg, *mri_norm, *mri_norm_atrophy, *mri_noise ;

  /* rkt: check for and handle version tag */
  nargs = 
    handle_version_option
    (argc, argv,
     "$Id: mri_simulate_atrophy.c,v 1.5 2011/03/02 00:04:24 nicks Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  memset(border_labels, -1, sizeof(border_labels)) ;
  border_labels[0][0] = Left_Lateral_Ventricle ;
  border_labels[0][1] = Left_Inf_Lat_Vent ;
  border_labels[0][2] = Unknown ;
  nlabels = 3 ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  mri_aseg = MRIread(argv[1]) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume %s",
              Progname, argv[1]) ;
  mri_norm = MRIread(argv[2]) ;
  if (mri_norm == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read norm volume %s",
              Progname, argv[2]) ;
  out_fname = argv[argc-1] ;

  mri_noise = MRIrandn(mri_norm->width, mri_norm->height, mri_norm->depth,1,
                       0, noise_sigma, NULL) ;
  MRImaskZero(mri_noise, mri_norm, mri_noise) ;
  mri_norm_atrophy = MRIsimulateAtrophy(mri_norm, mri_aseg,  target_labels[0], border_labels[0], nlabels, atrophy_pct, NULL) ;

  MRIadd(mri_norm_atrophy, mri_noise, mri_norm_atrophy) ;
  printf("writing simulated atrophy image to %s\n", out_fname) ;
  MRIwrite(mri_norm_atrophy, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr,
          "simulation took %d minutes and %d seconds.\n", minutes, seconds) ;

  exit(0) ;
  return(0) ;
}



static int get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else switch (toupper(*option))
  {
  case 'N':
    noise_sigma = atof(argv[2]) ;
    printf("using noise level %f\n", noise_sigma) ;
    nargs = 1 ;
    break ;
  case 'A':
    atrophy_pct = atof(argv[2]) ;
    printf("simulating %f atrophy level\n", atrophy_pct) ;
    nargs = 1 ;
    if (atrophy_pct > 1)
    {
      if (atrophy_pct > 100)
        ErrorExit(ERROR_BADPARM, "%s: must specify atrophy in [0 1]",
                  Progname) ;
      atrophy_pct /= 100 ; // assume it was specified in pct
    }
    break ;
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

static void usage_exit(int code)
{
  printf("usage: %s [options] <input aseg> <input norm> <output norm>\n",
         Progname) ;
  printf("\noptions:\n");
  printf("  -a <atrophy %%> - %% atrophy to simulate of structure\n");
  printf("  -n <sigma>     - gaussian noise level to add\n") ;
  exit(code) ;
}
MRI *
MRIsimulateAtrophy(MRI *mri_norm, MRI *mri_aseg,  int target_label, 
                   int *border_labels, int nlabels, float atrophy_pct, 
                   MRI *mri_norm_atrophy) 
{
  float  total_volume, pv, border_volume, vox_vol, pct_border_reduction ;
  MRI    *mri_mixing = MRIcloneDifferentType(mri_aseg, MRI_FLOAT),
         *mri_nbr_labels = MRIclone(mri_aseg, NULL);
  int    x, y, z ;
  int    nbr_label, vox_label, i ;
  float  mean_label, mean_nbr, val ;

  vox_vol = mri_norm->xsize * mri_norm->ysize * mri_norm->zsize ;
  mri_norm_atrophy = MRIcopy(mri_norm, NULL) ;

  MRIsetValues(mri_nbr_labels, Left_undetermined) ;
  total_volume = 
    MRIvoxelsInLabelWithPartialVolumeEffects(mri_aseg, mri_norm, target_label,
                                             mri_mixing, mri_nbr_labels) ;
  border_volume = 0 ;
  for (x = 0 ; x < mri_norm->width ; x++)
  {
    for (y = 0 ; y < mri_norm->height ; y++)
    {
      for (z = 0 ; z < mri_norm->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        vox_label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        nbr_label = MRIgetVoxVal(mri_nbr_labels, x, y, z, 0) ;
        for (i = 0 ; i < nlabels ; i++)
        {
          if (nbr_label == border_labels[i])
            break ;
        }
        if (i >= nlabels)
          continue ;   // not an allowable border label

        pv = MRIgetVoxVal(mri_mixing,x, y, z, 0) ;
        if (pv < 0)
          continue ;
        if (pv > 1)
        {
          pv = 1 ;
          DiagBreak() ;
        }
        border_volume += pv ;
      }
    }
  }
  pct_border_reduction = total_volume *  atrophy_pct / border_volume ;
  printf("total volume = %2.0fmm^3, border volume = %2.0f mm^3,reduction=%2.1f%%\n",
         total_volume, border_volume, 100*pct_border_reduction) ;
  pct_border_reduction = 1-pct_border_reduction ;
  for (x = 0 ; x < mri_norm->width ; x++)
  {
    for (y = 0 ; y < mri_norm->height ; y++)
    {
      for (z = 0 ; z < mri_norm->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        vox_label = MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        nbr_label = MRIgetVoxVal(mri_nbr_labels, x, y, z, 0) ;
        for (i = 0 ; i < nlabels ; i++)
        {
          if (nbr_label == border_labels[i])
            break ;
        }
        if (i >= nlabels)
          continue ;   // not an allowable border label

        // compute mean of this label
        pv = MRIgetVoxVal(mri_mixing, x, y, z, 0) ;
        mean_label = MRImeanInLabelInRegion(mri_norm, mri_aseg, target_label,
                                            x, y, z, 7) ;
        mean_nbr = MRImeanInLabelInRegion(mri_norm, mri_aseg, nbr_label,
                                          x, y, z, 7) ;
        val = MRIgetVoxVal(mri_norm, x, y, z, 0) ;
        pv = (val - mean_nbr) / (mean_label - mean_nbr) ;
        if (pv <= 0)
          continue ;
        if (pv > 1)
          pv = 1 ;
        pv *= pct_border_reduction ;
        val = mean_label * pv + mean_nbr * (1-pv) ;
        MRIsetVoxVal(mri_norm_atrophy, x, y, z, 0, val) ;
      }
    }
  }
  MRIfree(&mri_mixing) ; MRIfree(&mri_nbr_labels) ;
  return(mri_norm_atrophy) ;
}
