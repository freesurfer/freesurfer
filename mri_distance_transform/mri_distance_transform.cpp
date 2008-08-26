/**
 * @file  mri_distance_transform.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Florent Segonne 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2008/08/26 02:19:06 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

extern "C" {
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "timer.h"
#include "cma.h"
}
#include "fastmarching.h"

char *Progname ;

static int get_option(int argc, char *argv[]) ;
static int MRIexpandDistancesIntoGrayMatter(MRI *mri_distance, MRI *mri_aseg, MRI *mri_white, int target_label) ;

static MRI *mri_white = NULL, *mri_aseg = NULL ;
static char *surf_name = NULL ;

static int ndilations = 0 ;

int main(int argc, char *argv[]) {
  MRI *mri,*mri_distance;
  int label,mode;
  float max_distance;
  int   ac, nargs ;
  char  **av ;

  max_distance=10;
  mode=1;

  Progname=argv[0];

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  fprintf(stderr,"mri_distance_transform <input volume> <label> <max_distance> <mode[=1]> <output volume>\n");
  fprintf(stderr,"mode : 1 = outside , mode : 2 = inside , mode : 3 = both, mode : 4 = both unsigned \n");

  if (argc < 5)
    exit(0) ;

  mri=MRIread(argv[1]);
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s", Progname, argv[1]) ;
  label=atoi(argv[2]);
  max_distance=atof(argv[3]);
  mode=atoi(argv[4]);

  if (ndilations > 0)
      MRIdilateLabel(mri, mri, label, ndilations) ;
  if (mri->type != MRI_FLOAT)
    {
      MRI *mri_tmp;
      mri_tmp = MRIchangeType(mri, MRI_FLOAT, 0, 1, 1) ;
      MRIfree(&mri) ;
      mri = mri_tmp ;
    }
  fprintf(stderr,"label=%d distance=%.0f mode=%d\n",label,max_distance,mode);

  mri_distance=MRIalloc(mri->width,mri->height,mri->depth,MRI_FLOAT);
  MRIcopyHeader(mri, mri_distance) ;

  if (surf_name)
    {
      MRI_SURFACE *mris ;
      int         x, y, z ;

      mris = MRISread(surf_name) ;
      if (mris == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, surf_name) ;
      mri_white = MRIclone(mri, NULL) ;
      MRISfillInterior(mris, mri_white->xsize, mri_white) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_white, "w.mgz") ;
      mri_aseg = MRIdilate(mri_white, NULL) ;
      for (x = 0 ; x < mri_aseg->width ; x++)
        for (y = 0 ; y < mri_aseg->height ; y++)
          for (z = 0 ; z < mri_aseg->depth ; z++)
          {
            if (MRIgetVoxVal(mri_white, x, y, z, 0) != 0)
              MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_Cerebral_White_Matter) ;
            else if ((MRIgetVoxVal(mri_aseg, x, y, z, 0) != 0)) 
              MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_Cerebral_Cortex) ;
          }
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_aseg, "a.mgz") ;
    }

  mri_distance=MRIextractDistanceMap(mri,mri_distance,label, max_distance, mode, mri_white);

  if (mri_aseg)
    {
      MRIexpandDistancesIntoGrayMatter(mri_distance, mri_aseg, mri_white, label) ;
    }

  MRIwrite(mri_distance,argv[5]);

  MRIfree(&mri);
  MRIfree(&mri_distance);

  return 0;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int wm_labels[] =
  {
    Left_Cerebral_White_Matter,
    Right_Cerebral_White_Matter,
    Left_Cerebellum_White_Matter,
    Right_Cerebellum_White_Matter,
    Left_Thalamus_Proper,
    Right_Thalamus_Proper,
    Brain_Stem,
    Left_VentralDC,
    Right_VentralDC
  } ;
#define NWM_LABELS (int)((sizeof(wm_labels) / sizeof(wm_labels[0])))

static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  /*  StrUpper(option) ;*/
  if (!stricmp(option, "wm"))
    {
      int i ;
      nargs = 1 ;
      mri_aseg = MRIread(argv[2]) ;
      if (mri_aseg == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s", Progname, argv[2]) ;
      mri_white = MRIclone(mri_aseg, NULL);
      for (i = 0 ; i < NWM_LABELS ; i++)
        MRIcopyLabel(mri_aseg, mri_white, wm_labels[i]) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_white, "white.mgz");
    }
  else if (!stricmp(option, "wsurf"))
    {
      nargs = 1 ;
      surf_name = argv[2] ;
    }
  else if (!stricmp(option, "dilate"))
    {
      ndilations = atoi(argv[2]) ;
      nargs = 1;
      printf("performing %d dilations on labeled volume before computing distance transform\n", 
             ndilations) ;
    }
  return(nargs) ;
}

static int
MRIexpandDistancesIntoGrayMatter(MRI *mri_distance, MRI *mri_aseg, MRI *mri_white, int target_label)
{
  int   x, y, z, xi, yi, zi, xk, yk, zk, label ;
  float dist, min_dist, vdist ;

  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (label == target_label)
          MRIsetVoxVal(mri_distance, x, y, z, 0, 0) ;
      }

  // sample from wm into adjacent gm to prevent jumping the banks of a sulcus
  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (IS_GM(label) == 0)
          continue ;
        min_dist = MRIgetVoxVal(mri_distance, x, y, z, 0) ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = mri_aseg->xi[x+xk] ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = mri_aseg->yi[y+yk] ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = mri_aseg->zi[z+zk] ;
              /* only allow distances to propagate through the wm.
                 This will prevent the distances from jumping
                 across the banks of a sulcus (which is the whole point)
              */
              if (MRIgetVoxVal(mri_white, xi, yi, zi, 0) == 0)
                continue ;  
              vdist = sqrt(xk*xk + yk*yk + zk*zk);
              zi = mri_aseg->zi[z+zk] ;
              dist = MRIgetVoxVal(mri_distance, xi, yi, zi, 0) + vdist ;
              if (dist < min_dist)
                min_dist = dist ;
            }
          }
        }
        MRIsetVoxVal(mri_distance, x, y, z, 0, min_dist);
      }

  return(NO_ERROR) ;
}

