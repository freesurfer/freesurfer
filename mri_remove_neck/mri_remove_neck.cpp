/**
 * @brief removes neck fromm t1 (stuff below brain)
 *
 */
/*
 * Original Author: Bruce Fischl
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
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"

const char         *Progname ;

static int fill_val = 0 ;
static int radius = 25 ;

static double TR = 0.0 ;
static double TE = 0.0 ;
static double alpha = 0.0 ;

static MRI * fill_brain_volume(MRI *mri, GCA *gca, TRANSFORM *transform, int radius) ;
MRI *MRIremoveNonBrain(MRI *mri_src, MRI *mri_dst, TRANSFORM *transform, GCA *gca, int radius, int fill_val) ;
static int get_option(int argc, char *argv[]) ;

static void usage_exit(int code);

/*
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/


int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, **av, *transform_fname ;
  MRI          *mri_in, *mri_out ;
  GCA          *gca ;
  int          ac, nargs ;
  int          msec, minutes, seconds,err ;
  Timer start ;
  TRANSFORM    *transform ;

  Progname = argv[0] ;
  setRandomSeed(-1L) ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  nargs = handleVersionOption(argc, argv, "mri_remove_neck");
  argc -= nargs ;
  if (1 == argc)
  {
    exit (0);
  }

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
  {
    usage_exit(0);
  }

  in_fname = argv[1] ;
  transform_fname = argv[2] ;
  gca_fname = argv[3] ;
  out_fname = argv[4] ;

  start.reset() ;
  printf("reading atlas '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n", Progname, gca_fname) ;
  }

  printf("reading input volume '%s'...\n", in_fname) ;
  fflush(stdout) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n", Progname, in_fname) ;
  }

  printf("reading transform '%s'...\n", transform_fname) ;
  fflush(stdout) ;
  transform = TransformRead(transform_fname) ;
  if (!transform)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not open transform %s.\n", Progname, transform_fname) ;
  }
  TransformInvert(transform, mri_in) ;

  printf("removing structures at least %d mm from brain...\n", radius) ;
  fflush(stdout) ;
  mri_out = MRIremoveNonBrain(mri_in, NULL, transform, gca, radius, fill_val) ;
  printf("writing output to %s...\n", out_fname) ;
  fflush(stdout) ;

  GCAfree(&gca) ;
  MRIfree(&mri_in) ;

  err=MRIwrite(mri_out, out_fname);
  if (err)
  {
    exit(1);
  }

  MRIfree(&mri_out)  ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("nonbrain removal took %d minutes and %d seconds.\n", minutes, seconds) ;
  fflush(stdout) ;
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
  StrUpper(option) ;
  if (!stricmp(option, "FILL"))
  {
    fill_val = atoi(argv[2]) ;
    nargs = 1 ;
    printf("filling defaced regions with %d\n", fill_val) ;
  }
  else if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage_exit(0);
  }
  else if (!stricmp(option, "RADIUS"))
  {
    radius = atoi(argv[2]) ;
    nargs = 1 ;
    printf("erasing everything more than %d mm from possible brain\n", radius) ;
  }
  else if (!stricmp(option, "DEBUG_LABEL"))
  {
    Ggca_label = atoi(argv[2]) ;
    printf("debugging label %s (%d)\n", cma_label_to_name(Ggca_label), Ggca_label) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!stricmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!stricmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else switch (*option)
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      usage_exit(0);
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

MRI *
MRIremoveNonBrain(MRI *mri_src, MRI *mri_dst, TRANSFORM *transform, GCA *gca, int radius, int fill_val)
{
  int       x, y, z, frame, nerased = 0 ;
  MRI       *mri_brain ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_brain = fill_brain_volume(mri_dst, gca, transform, radius) ;

  for (x = 0 ; x < mri_src->width ; x++)
  {
    for (y = 0 ; y < mri_src->height ; y++)
    {
      for (z = 0 ; z < mri_src->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        if (MRIvox(mri_brain, x, y, z)  > 0)
        {
          continue ;
        }

        for (frame = 0 ; frame < mri_dst->nframes ; frame++)
        {
          MRIseq_vox(mri_dst, x, y, z, frame) = fill_val ;
        }
        nerased++ ;
      }
    }
  }

  printf("%d nonbrain voxels erased\n", nerased) ;
  MRIfree(&mri_brain) ;
  return(mri_dst)  ;
}

static MRI *
fill_brain_volume(MRI *mri, GCA *gca, TRANSFORM *transform, int radius)
{
  int       i, x, y, z, n ;
  MRI       *mri_brain ;
  GCA_PRIOR *gcap ;

  mri_brain = MRIclone(mri, NULL) ;

  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
        {
          DiagBreak() ;
        }
        MRIvox(mri_brain, x, y, z) = 0 ;
        gcap = getGCAP(gca, mri, transform, x, y, z) ;
        if (gcap == NULL)
        {
          continue ;
        }
        for (n = 0 ; n < gcap->nlabels ; n++)
        {
          if (IS_BRAIN(gcap->labels[n]) && (gcap->labels[n] != Brain_Stem))
          {
            MRIvox(mri_brain, x, y, z) = 128 ;
            break ;
          }
        }
      }
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_brain, "brain_mask.mgz") ;
  }
  for (i = 0 ; i < radius ; i++)
  {
    MRIdilate(mri_brain, mri_brain) ;
  }
  return(mri_brain) ;
}

#include "mri_remove_neck.help.xml.h"
static void
usage_exit(int code)
{
  outputHelpXml(mri_remove_neck_help_xml,mri_remove_neck_help_xml_len);
  exit(code) ;
}
