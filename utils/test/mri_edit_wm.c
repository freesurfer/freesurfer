/*
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "cma.h"
#include "version.h"

void usage(int exit_val);
MRI *fliplr(MRI *src);

const char *Progname;

/* Modified from mri_fill.c */
static int edit_segmentation(MRI *mri_im, MRI *mri_seg) ;

int main(int argc, char *argv[])
{

  MRI *mri_wm, *mri_seg;
  int nargs;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_edit_wm");
  argc -= nargs ;
  if (1 == argc)
    exit (0);

  if (argc != 4)
    usage(1);

  mri_wm = MRIread(argv[1]) ;
  if (!mri_wm)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;
  mri_seg = MRIread(argv[2]) ;
  if (!mri_seg)
    ErrorExit(ERROR_BADPARM, "%s: could not read mask volume %s",
              Progname, argv[1]) ;

  if (mri_wm->width != mri_seg->width ||
      mri_wm->height != mri_seg->height ||
      mri_wm->depth != mri_seg->depth)
  {
    ErrorExit(ERROR_BADPARM, "%s: the two input volumes have unequal size. Exit.",
              Progname) ;
  }

  edit_segmentation(mri_wm, mri_seg);

  printf("writing edited volume to %s...\n", argv[3]) ;
  MRIwrite(mri_wm, argv[3]);

  MRIfree(&mri_wm);
  MRIfree(&mri_seg);

  exit(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <in vol> <aseg vol> <out vol>\n", Progname);
  fprintf(fout, "this program fills in subcortical regions using aseg \n") ;

  exit(exit_val);

}  /*  end usage()  */
/*  EOF  */

static int
edit_segmentation(MRI *mri_wm, MRI *mri_seg)
{
  int   width, height, depth, x, y, z, label, non;

  width = mri_wm->width ;
  height = mri_wm->height ;
  depth = mri_wm->depth ;

  non = 0 ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        label = (int) MRIgetVoxVal(mri_seg, x, y,z,0);

        switch (label)
        {

          /* fill these */
        case Left_Thalamus_Proper:
        case Right_Thalamus_Proper:
        case Left_Caudate:
        case Right_Caudate:
        case Left_Putamen:
        case Right_Putamen:
        case Left_Lateral_Ventricle:
        case Right_Lateral_Ventricle:
        case Left_Inf_Lat_Vent:
        case Right_Inf_Lat_Vent:
          if (MRIgetVoxVal(mri_wm, x, y, z, 0) < 110)
          {
            MRIsetVoxVal(mri_wm, x, y, z, 0, 110);
            non++ ;
          }
          break ;
        default:
          break ;
        }
      }
    }
  }

  printf("SEG EDIT: %d voxels turned on.\n", non) ;

  return(NO_ERROR) ;
}
