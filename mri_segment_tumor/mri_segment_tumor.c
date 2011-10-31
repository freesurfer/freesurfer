/**
 * @file  mri_segment_tumor.c
 * @brief program for segmenting tumors from multispectral data.
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2011/10/31 18:30:41 $
 *    $Revision: 1.3 $
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



/*!
\file mri_segment_tumor.c
\brief Example c file that can be used as a template.
\author Douglas Greve

*/


// $Id: mri_segment_tumor.c,v 1.3 2011/10/31 18:30:41 fischl Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "gca.h"
#include "transform.h"

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_segment_tumor.c,v 1.3 2011/10/31 18:30:41 fischl Exp $";
char *Progname = NULL;
//char *cmdline, cwd[2000];
//struct utsname uts;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;
static float mthresh = 2 ;

MRI *
GCAsegmentTumor(GCA* gca, TRANSFORM  *transform, MRI *mri_inputs, MRI *mri_dst, float mthresh)
{
  int             x, y, z, xn, yn, zn, label, n ;
  GCA_NODE        *gcan ;
  GCA_PRIOR       *gcap ;
  GC1D            *gc ;
  float           vals[MAX_GCA_INPUTS], mdist, min_mdist ;
  MRI             *mri_mask ;

  mri_mask = MRIalloc(mri_inputs->width, mri_inputs->height, mri_inputs->depth, mri_inputs->type) ;
  MRIcopyHeader(mri_inputs, mri_mask) ;
  MRIsetVoxelsWithValue(mri_inputs, mri_mask, 0, 1) ;
  MRIdilate(mri_mask, mri_mask) ;
  

  if (mri_dst == NULL)
  {
    mri_dst = MRIalloc(mri_inputs->width, mri_inputs->height, mri_inputs->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  for (x = 0 ; x < mri_inputs->width ; x++)
  {
    for (y = 0 ; y < mri_inputs->height ; y++)
    {
      for (z = 0 ; z < mri_inputs->depth ; z++)
      {
	if (MRIgetVoxVal(mri_mask, x, y, z, 0))
	  continue ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        if (!GCAsourceVoxelToNode(gca, mri_inputs,
                                  transform, x, y, z, &xn, &yn, &zn))
        {
          load_vals(mri_inputs, x, y, z, vals, gca->ninputs);
          gcan = &gca->nodes[xn][yn][zn] ;
          gcap = getGCAP(gca, mri_inputs, transform, x, y, z) ;
          MRIsetVoxVal(mri_dst, x, y, z, 0, 0) ;
          if (gcap==NULL)
            continue;
	  for (n = 0 ; n < gca->ninputs ; n++)
	    if (!FZERO(vals[n]))
	      break ;
	  if (n == gca->ninputs) // all inputs are 0 - set output to 0
	    continue ;
	  if (gcan->nlabels == 0)
	    continue ;

          min_mdist = 10000 ;
          for (label = 1, n = 0 ; n < gcan->nlabels ; n++)
          {
            gc = &gcan->gcs[n] ;
            mdist = sqrt(GCAmahDist(gc, vals, gca->ninputs)) ;
            if (mdist < mthresh)  // there is some label that is not unlikely
              label = 0 ;
            if (mdist < min_mdist)
              min_mdist = mdist ;
          }

          MRIsetVoxVal(mri_dst, x, y, z, 0, min_mdist) ;
        }
      }
    }
  }
  MRIfree(&mri_mask) ;
  return(mri_dst) ;
}
/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  char      **av ;
  int       nargs, ac, i, ninputs ;
  MRI       *mri_in = NULL, *mri_tmp, *mri_labeled ;
  GCA       *gca ;
  TRANSFORM *xform ;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");

  if (nargs && argc - nargs == 1) 
    exit (0);
  argc -= nargs;
  //  cmdline = argv2cmdline(argc,argv);
  //  uname(&uts);
  //  getcwd(cwd,2000);

  Progname = argv[0] ;
#if 0
  argc --;
  argv++;
#endif
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    usage_exit(1) ;

  ninputs = argc - 4 ;

  for (i = 0 ; i < ninputs ; i++)
  {
    printf("reading input %d: %s\n", i, argv[i+1]) ;
    mri_tmp = MRIread(argv[i+1]) ;
    if (mri_tmp == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s", Progname, argv[i+1]) ;
    if (i == 0)
    {
      mri_in = MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type, ninputs) ;
      MRIcopyHeader(mri_tmp, mri_in) ;
    }
    MRIcopyFrame(mri_tmp, mri_in, 0, i) ;
    MRIfree(&mri_tmp) ;
  }
  printf("reading gca from %s\n",argv[ninputs+1]) ;fflush(stdout) ;
  gca = GCAread(argv[ninputs+1]) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read gca from %s", Progname, argv[ninputs+1]) ;
  xform = TransformRead(argv[2+ninputs]) ;
  if (xform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, argv[2+ninputs]) ;

  mri_labeled = GCAsegmentTumor(gca, xform, mri_in, NULL, mthresh) ;
  printf("writing labeled tumor to %s\n", argv[3+ninputs]) ;
  MRIwrite(mri_labeled, argv[3+ninputs]) ;
  return 0;
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } 
  else switch (toupper(*option)) {
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  case 'T':
    mthresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("using Mahalanobis distance %2.1f\n", mthresh) ;
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
usage_exit(int code) {
  printf("usage: %s [options] <input vol 1> <input vol 2> ... <gca> <output labeling>",
         Progname) ;
  exit(code) ;
}
