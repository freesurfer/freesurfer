/**
 * @brief extract cortical ribbon
 *
 * Usage: mri_ribbon \
 *              inner_surface_fname outer_surface_fname
 *              input_volume_pref output_volume_pref
 */
/*
 * Original Author: Andre van der Kouwe
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
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "mri.h"
#include "diag.h"
#include "macros.h"
#include "mrisurf.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "error.h"
#include "mrishash.h"

const char *Progname;

MRI *MRIcropVolumeToLabel(MRI *mri_src, 
                          MRI *mri_dst, 
                          LABEL *area, 
                          MRI_SURFACE *mris_white, 
                          MRI_SURFACE *mris_pial);

int main(int argc, char *argv[]) {
  char *inner_mris_fname,*outer_mris_fname,*input_mri_pref,*output_mri_pref;
  MRI *mri=0,*mri_src=0;
  MRI_SURFACE *inner_mris=0,*outer_mris=0;
  int nargs;
  LABEL *area = NULL ;

  nargs = handleVersionOption(argc, argv, "mri_ribbon");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname=argv[0];
  if (argc > 1)
    while (*argv[1] == '-')
    {
      int nargs = 0 ;
      switch (toupper(argv[1][1]))
      {
      case 'L':
        printf("cropping ribbon to label file %s\n", argv[2]) ;
        nargs = 1 ;
        area = LabelRead(NULL, argv[2]) ;
        if (area == NULL)
          ErrorExit(ERROR_NOFILE, 
                    "%s: could not read label file %s\n", 
                    argv[2]) ;
        break ;
      default:
        break ;
      }
      argc -= (nargs+1) ;
      if (argc == 1) break;
      argv += (nargs+1) ;
    }

  /* Set command-line parameters */
  if (argc!=5) {
    printf("Usage: mri_ribbon [-l fname.label] \\ \n"
           "       inner_surface_fname \\ \n"
           "       outer_surface_fname \\ \n"
           "       input_volume_pref \\ \n"
           "       output_volume_pref\n");
    exit(1);
  }

  inner_mris_fname=argv[1];
  outer_mris_fname=argv[2];
  input_mri_pref=argv[3];
  output_mri_pref=argv[4];

  /* Read surface information from inner surface file */
  printf("Reading surface file %s.\n",inner_mris_fname);
  inner_mris=MRISread(inner_mris_fname);
  if (!inner_mris) {
    printf("Could not read surface file %s.\n",inner_mris_fname);
    exit(1);
  }

  /* Read surface information from outer surface file */
  printf("Reading surface file %s.\n",outer_mris_fname);
  outer_mris=MRISread(outer_mris_fname);
  if (!outer_mris) {
    printf("Could not read surface file %s.\n",outer_mris_fname);
    exit(1);
  }

  /* Read example volume from file */
  printf("Reading MRI volume directory %s.\n",input_mri_pref);
  mri_src=MRIread(input_mri_pref);
  if (!mri_src) {
    printf("Could not read MRI volume directory %s.\n",input_mri_pref);
    exit(1);
  }

  /* Extract ribbon */
  printf("Extracting ribbon.\n");
  mri=MRISribbon(inner_mris,outer_mris,mri_src,NULL);

  if (area)
    MRIcropVolumeToLabel(mri, mri, area, inner_mris, outer_mris) ;
  /* Save MRI volume to directory */
  printf("Writing volume file %s.\n",output_mri_pref);
  MRIwrite(mri,output_mri_pref);

  MRIfree(&mri);
  MRIfree(&mri_src);
  MRISfree(&inner_mris);
  MRISfree(&outer_mris);

  return 0;
}

MRI *
MRIcropVolumeToLabel(MRI *mri_src, 
                     MRI *mri_dst, 
                     LABEL *area, 
                     MRI_SURFACE *mris_white, 
                     MRI_SURFACE *mris_pial)
{
  MHT    *mht_white, *mht_pial ;
  int    x, y, z ;
  VERTEX *v_white, *v_pial ;
  double xs, ys, zs ;

  mht_white = MHTcreateVertexTable_Resolution(mris_white, CURRENT_VERTICES, 5.0) ;
  mht_pial = MHTcreateVertexTable_Resolution (mris_pial,  CURRENT_VERTICES, 5.0) ;
  if (mri_dst == NULL)
    mri_dst = MRIclone(mri_src, NULL) ;

  MRISclearMarks(mris_white) ; MRISclearMarks(mris_pial) ;
  LabelMarkSurface(area, mris_white) ; LabelMarkSurface(area, mris_pial) ;
  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
              
        if (MRIgetVoxVal(mri_src, x, y, z, 0) > 0)
        {
          MRISsurfaceRASFromVoxel(mris_white, mri_dst, 
                                  x, y, z, 
                                  &xs, &ys, &zs) ;
          v_white = MHTfindClosestVertexInTable(mht_white, mris_white, 
                                                xs, ys, zs, 1) ;
          v_pial = MHTfindClosestVertexInTable(mht_pial, mris_pial, 
                                               xs, ys, zs, 1) ;
          if ((v_white && v_white->marked == 1) ||
              (v_pial && v_pial->marked == 1))
          {
            MRIsetVoxVal(mri_dst, x, y, z, 0, 255) ;
          }
          else
            MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
        }
      }
    }
  }

  MHTfree(&mht_white) ; MHTfree(&mht_pial) ;
  return(mri_dst) ;
}
