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


/* Compute the mean and std of a labeled region */

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
#include "version.h"

void usage(int exit_val);
MRI *fliplr(MRI *src);

const char *Progname;

static int label1 = 0;
static int label2 = 0;

/* Compute CNR inside a window only */

static int debug_flag = 0;
static int window_flag = 0;
static int window_size = 30;

static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{

  char **av;
  MRI *mri_src, *mri_label, *mri_mask;
  int ac, nargs;
  int width, height, depth, x, y, z;
  int cnrflag = 0;
  int count;

  float v1, v2;
  double meanV1, stdV1, meanV2, stdV2;
  int total1, total2;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_compute_stats");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if ( argc != 2 && argc != 4 && argc != 5)
    usage(1);

  mri_src = MRIread(argv[1]) ;
  if (!mri_src)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;

  if (argc == 2)
  {
    count = 0;
    for (z=0; z < mri_src->depth; z++)
      for (y=0; y< mri_src->height; y++)
        for (x=0; x < mri_src->width; x++)
        {
          if (MRIgetVoxVal(mri_src,x, y, z, 0) > 1e-30)
            count++;
        }
    printf("# foreground voxels = %d\n", count);
    MRIfree(&mri_src);
    exit(0);
  }

  mri_label = MRIread(argv[2]) ;
  if (!mri_label)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume %s",
              Progname, argv[2]) ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  mri_mask = MRIclone(mri_label, NULL);

#if 0
  printf("This program is borrowed to create things");
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        if (MRIvox(mri_label, x, y, z) == 17 || MRIvox(mri_label, x, y, z) == 53)
          MRIvox(mri_mask, x, y,z) = 255;
        else
          MRIvox(mri_mask, x, y,z) = 0;
      }

  MRIwrite(mri_mask, "hippo.mgz");
#endif

  if ((mri_src->width != mri_label->width) ||
      (mri_src->height != mri_label->height) ||
      (mri_src->depth != mri_label->depth))
    ErrorExit(ERROR_BADPARM, "%s: source (type %d) and label (type %d) volumes don't match (%d x %d x %d) vs (%d x %d x %d)\n",
              Progname, mri_src->type, mri_label->type, mri_src->width,
              mri_src->height, mri_src->depth,
              mri_label->width, mri_label->height, mri_label->depth) ;


  label1 = atoi(argv[3]);

  printf("Region of interest has label = %d\n", label1);

  if (argc == 5)
  {
    label2 = atoi(argv[4]);
    cnrflag = 1;
    printf("Second region has label = %d\n", label2);
  }
  else
    cnrflag = 0;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;

  mri_mask = MRIalloc(mri_src->width, mri_src->height, mri_src->depth, MRI_UCHAR);
  MRIcopyHeader(mri_src, mri_mask);

  /* Simply set mask to be 1 everywhere */
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++)
      {
        MRIvox(mri_mask, x, y,z) = 1;
      }

  if (debug_flag && window_flag)
  {
    /* Limit LDA to a local window */
    printf("Local window size = %d\n", window_size);
    window_size /= 2;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++)
        {

          if (z < (Gz - window_size) || z >(Gz + window_size)
              || y <(Gy - window_size) || y > (Gy + window_size)
              || x < (Gx - window_size) || x > (Gx + window_size))
            MRIvox(mri_mask, x, y,z) = 0;
        }

  }


  total1 = 0;
  meanV1 = 0.0;
  total2 = 0;
  meanV2 = 0.0;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_mask, x, y,z) == 0) continue;

        v2 = MRIgetVoxVal(mri_label,x,y,z,0);
        v1 = MRIgetVoxVal(mri_src,x,y,z,0);
        if ( v2 < label1 + 0.000001 && v2 > label1 -0.000001)
        {
          total1++;
          meanV1 += v1;
        }

        if (cnrflag)
        {
          if ( v2 < label2 + 0.000001 && v2 > label2 -0.000001)
          {
            total2++;
            meanV2 += v1;
          }
        }

      }
    }
  }

  if (total1 == 0)
  {
    printf("The region with label %d is empty! \n", label1);

    MRIfree(&mri_src);
    MRIfree(&mri_label);
    exit(0);
  }

  meanV1 /= (float)total1;
  if (cnrflag)
    meanV2 /= (float)(total2 + 0.000000000001);

  stdV1 = 0.0;
  stdV2 = 0.0;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_mask, x, y,z) == 0) continue;

        v2 = MRIgetVoxVal(mri_label,x,y,z,0);
        v1 = MRIgetVoxVal(mri_src,x,y,z,0);
        if ( v2 < label1 + 0.000001 && v2 > label1 -0.000001)
        {
          stdV1 += (v1-meanV1)*(v1-meanV1);
        }

        if (cnrflag)
        {
          if ( v2 < label2 + 0.000001 && v2 > label2 -0.000001)
          {
            stdV2 += (v1-meanV2)*(v1-meanV2);
          }
        }
      }
    }
  }

  stdV1  = sqrt(stdV1/total1);
  if (cnrflag)
    stdV2  = sqrt(stdV2/(total2 + 0.00000001));
  printf("Region with label %d has size = %d, mean = %g, std = %g\n", label1, total1, meanV1, stdV1);
  if (cnrflag)
  {
    printf("Region with label %d has size = %d, mean = %g, std = %g\n", label2, total2, meanV2, stdV2);

    stdV1 = (meanV1 - meanV2)/sqrt((stdV1*stdV1 + stdV2*stdV2)*0.5 + 0.0000000001);
    if (stdV1 < 0) stdV1 = -stdV1;

    printf("CNR of regions %d and %d is equal to %g\n", label1, label2, stdV1);
  }

  MRIfree(&mri_src);
  MRIfree(&mri_label);

  exit(0);

}  /*  end main()  */

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    debug_flag = 1;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "window"))
  {
    window_flag = 1 ;
    window_size = atoi(argv[2]) ;
    nargs = 1;
    printf("interpolating volume to be isotropic 1mm^3\n") ;
  }

  return(nargs) ;
}

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <in vol> <label vol> label [label2] \n", Progname);
  fprintf(fout, "this program computes the mean and std of a labeled region. \n") ;
  fprintf(fout, "If label2 is also specified, the program computes the CNR between the two labelled regions. \n") ;
  fprintf(fout, "If only <in vol> is given, it will output the volume of all foreground voxels in it. \n") ;

  exit(exit_val);

}  /*  end usage()  */
/*  EOF  */
