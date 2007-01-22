/**
 * @file  mris_wm_volume.c
 * @brief computes the interior volume of the ?h.white surfaces (excluding subcortical labels)
 *
 * computes the interior volume of the ?h.white surfaces (excluding subcortical labels),
 * uses the aseg.mgz and the ?h.white surfaces.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/01/22 14:26:38 $
 *    $Revision: 1.1 $
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


//

// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2007/01/22 14:26:38 $
// Revision       : $Revision: 1.1 $
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#include "cma.h"


double MRIScomputeWhiteVolume(MRI_SURFACE *mris, MRI *mri_aseg, double resolution) ;
int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;

static void usage_exit(int code) ;
static char sdir[STRLEN] = "" ;
static char *white_name = "white" ;
static char *aseg_name = "aseg.mgz" ;

static double resolution = 1.0/4.0 ;

int
main(int argc, char *argv[]) 
{
  char          **av, *hemi, fname[STRLEN], *cp, *subject;
  int           ac, nargs;
  MRI_SURFACE   *mris;
  int           msec, minutes, seconds;
  struct timeb start ;
  double        total_volume ;
  MRI           *mri_aseg ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_wm_volume.c,v 1.1 2007/01/22 14:26:38 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(sdir)) 
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }
  subject = argv[1] ;
  hemi = argv[2] ;
  if (argc < 3)
    usage_exit(1) ;

  //  printf("command line parsing finished\n");

  /*** Read in the input surfaces ***/
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject, hemi, white_name) ;
  mris = MRISread(fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, fname) ;

  /*** Read in the aseg volume ***/
  sprintf(fname, "%s/%s/mri/%s", sdir, subject, aseg_name) ;
  mri_aseg = MRIread(fname) ;
  if (mri_aseg == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read segmentation volume %s", Progname, fname) ;

  total_volume = MRIScomputeWhiteVolume(mris, mri_aseg, resolution) ;
  printf("%g\n", total_volume);
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;


  MRISfree(&mris);
  MRIfree(&mri_aseg) ;

  exit(0);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "white"))
  {
    white_name = argv[2] ;
    fprintf(stderr, "using %s as white surface name", white_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "aseg"))
  {
    aseg_name = argv[2] ;
    fprintf(stderr, "using %s as aseg volume name", aseg_name) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    fprintf(stderr, "using %s as SUBJECTS_DIR", sdir) ;
    nargs = 1 ;
  }
  else switch (*option) {
  default:
    printf("unknown option %s\n", argv[1]);
    exit(1);
    break;
  }

  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s <subject> <hemi>\n", Progname) ;
  printf("\t This program computes the volume of the enclosed ?h.white surface,\n\tignoring non-wm voxels in the aseg.\n");
  printf("\t use -v option to output more messages\n");
  printf("  -SDIR SUBJECTS_DIR \n");
  printf("  -white whitesurfname \n");
  printf("  -aseg  asegname \n");
  exit(code) ;
}

double
MRIScomputeWhiteVolume(MRI_SURFACE *mris, MRI *mri_aseg, double resolution)
{
  MRI    *mri_filled ;
  MATRIX *m_vox2vox ;
  double total_volume, vox_volume ;
  int    x, y, z, label, xa, ya, za ;
  VECTOR *v1, *v2 ;
  Real   val ;

  mri_filled = MRISfillInterior(mris, resolution, NULL) ;
  m_vox2vox = MRIgetVoxelToVoxelXform(mri_filled, mri_aseg) ;
  v1 = VectorAlloc(4, MATRIX_REAL) ; v2 = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v1, 4) = 1.0 ; VECTOR_ELT(v2, 4) = 1.0 ;
  vox_volume = mri_filled->xsize * mri_filled->ysize * mri_filled->zsize ;

  for (x = 0 ; x < mri_filled->width ; x++)
  {
    V3_X(v1) = x ;
    for (y = 0 ; y < mri_filled->height ; y++)
    {
      V3_Y(v1) = y ;
      for (z = 0 ; z < mri_filled->depth ; z++)
      {
        val = MRIgetVoxVal(mri_filled, x, y, z, 0) ;
        if (FZERO(val))
          continue ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        V3_Z(v1) = z ;
        MatrixMultiply(m_vox2vox, v1, v2) ;
        xa = nint(V3_X(v2)) ; ya = nint(V3_Y(v2)) ; za = nint(V3_Z(v2)) ;
        if (xa < 0 || xa >= mri_aseg->width ||
            ya < 0 || ya >= mri_aseg->height ||
            za < 0 || za >= mri_aseg->depth)
          continue ;
        label = (int)MRIgetVoxVal(mri_aseg, xa, ya, za, 0) ;
        if (xa == Gx && ya == Gy && za == Gz)
          DiagBreak() ;
        switch (label)
        {
        case Left_Cerebral_Cortex:
        case Right_Cerebral_Cortex:
        case Left_Cerebral_White_Matter:
        case Right_Cerebral_White_Matter:
        case Left_WM_hypointensities:
        case Right_WM_hypointensities:
          total_volume += vox_volume ;
          break ;
        }
      }
    }
  }

  MatrixFree(&m_vox2vox) ; MatrixFree(&v1) ; MatrixFree(&v2) ;
  MRIfree(&mri_filled) ;
  return(total_volume) ;
}

