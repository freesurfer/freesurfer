/**
 * @brief extract the largest connected component from a volume and write it into another volume
 *
 * extract the largest connected component from a volume and write it into another volume
 */
/*
 * Original Author: Florent Segonne
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


//Extract the largest connected component of a binary segmentation

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
#include "connectcomp.h"
#include "mrisurf.h"

static int lh_label = LH_LABEL ;
static int rh_label = RH_LABEL ;

void usage(int exit_val);
static int get_option(int argc, char *argv[]) ;

static float threshold = 90;

static char hemi[80] = "lh";
static int inverse = 0 ;

const char *Progname;
static MRI *mri_orig = NULL ;
static int target_value = 255;

int main(int argc, char *argv[]) {
  char   **av;
  MRI *mri_seg;
  int ac, nargs;

  int x, y, z;


  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_extract_largest_CC");
  argc -= nargs ;
  if (1 >= argc)
    usage (0);

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 3)
    usage(1);

  if (target_value == 255)
  {
    if (!stricmp(hemi, "lh")) {
      target_value = lh_label;
    } else {
      target_value = rh_label;
    }
  }

  mri_seg = MRIread(argv[1]) ;
  if (!mri_seg)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume %s",
              Progname, argv[1]) ;
  if (mri_orig == NULL)
    mri_orig = MRIcopy(mri_seg, NULL) ;

  for (z = 0; z < mri_seg->depth; z++)
    for (y = 0; y < mri_seg->height; y++)
      for (x = 0; x < mri_seg->width; x++) {
	if (target_value != rh_label && target_value != lh_label)
	{
	  if (MRIgetVoxVal(mri_seg, x, y, z, 0) == target_value)
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 255);
	  else
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 0);
	} 
	else if (inverse == 0)
	{
	  if (MRIgetVoxVal(mri_seg, x, y, z, 0) < threshold)
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 0);
	  else
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 255);
	}
	else  // extract background
	{
	  if (MRIgetVoxVal(mri_seg, x, y, z, 0) > threshold)
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 0);
	  else
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 255);
	}
      }

  GetLargestCC6(mri_seg);

  for (z = 0; z < mri_seg->depth; z++)
    for (y = 0; y < mri_seg->height; y++)
      for (x = 0; x < mri_seg->width; x++) {
	if (inverse == 0)
	{
	  if (MRIgetVoxVal(mri_seg, x, y, z, 0) > 0)
	    MRIsetVoxVal(mri_seg, x, y, z, 0, target_value);
	}
	else
	{
	  if (MRIgetVoxVal(mri_seg, x, y, z, 0) == 0)
	  {
	    MRIsetVoxVal(mri_seg, x, y, z, 0, MRIgetVoxVal(mri_orig, x, y, z, 0));
	  }
	  else
	    MRIsetVoxVal(mri_seg, x, y, z, 0, 0) ;
	}
      }

  printf("writing output to %s\n", argv[2]) ;
  MRIwrite(mri_seg, argv[2]);

  exit(0);

}  /*  end main()  */

void usage(int exit_val) {

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <input vol> <output vol>\n", Progname);
  fprintf(fout, "this program extracts the largest connected component of the input volume \n") ;
  fprintf(fout, "\t Options are: \n") ;
  fprintf(fout, "\t\t -T #: threshold for object \n") ;
  fprintf(fout, "\t\t -hemi lh/rh: set the target value corresponding to lh (255) or rh (127) \n") ;
  fprintf(fout, "\t\t -I : find the largest CC in the background\n") ;
  fprintf(fout, "\t\t -O <orig volume> : clone values from <orig volume> into output (used with -I)\n") ;
  fprintf(fout, "\t\t -L <label val> : perform connected components on voxels with value <label val>\n") ;
  exit(exit_val);

}  /*  end usage()  */


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "T")) {
    threshold = atof(argv[2]);
    printf("threshold = %g\n", threshold);
    nargs = 1;
  } else if (!stricmp(option, "hemi")) {
    strcpy(hemi,argv[2]);
    printf("hemisphere = %s\n", hemi);
    nargs = 1;
  } else switch (toupper(*option)) {
  case 'O':
    mri_orig = MRIread(argv[2]) ;
    if (mri_orig == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not open orig volume %s", Progname, argv[2]) ;
    nargs = 1 ;
    break ;
    case 'L':
      target_value = atoi(argv[2]) ;
      printf("using target value %d\n", target_value) ;
      nargs = 1 ;
      break ;
  case 'I':
    inverse = 1 ;
    printf("extracting CC from background instead of foreground\n") ;
    break ;
    case '?':
    case 'U':
      usage(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

/*  EOF  */
