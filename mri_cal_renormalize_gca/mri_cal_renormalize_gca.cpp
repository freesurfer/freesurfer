/**
 * @brief Normalize a set of longituindal volumes making use of subcortical atlas data
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
#include "utils.h"
#include "gca.h"
#include "tags.h"
#include "cma.h"
#include "mrinorm.h"
#include "version.h"


const char *Progname ;


static char *mask_fname = NULL ;
static int novar = 0 ;

static void usage_exit(int code) ;
static int get_option(int argc, char *argv[]) ;

static int longinput = 0;


/*
  command line consists of three inputs:

  argv[1]  - input volume
  argv[2]  - atlas (gca)
  argv[3]  - transform (lta/xfm/m3d)
  argv[4]  - output atlas
*/


#define MAX_TIMEPOINTS 2000
static char *subjects[MAX_TIMEPOINTS] ;
int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, **av, *xform_fname, fname[STRLEN] ;
  MRI          *mri_in, *mri_tmp ;
  GCA          *gca ;
  int          ac, nargs, msec, minutes, seconds;
  int          input, ninputs ;
  Timer start ;
  TRANSFORM    *transform ;
  char         line[STRLEN], *cp, sdir[STRLEN], base_name[STRLEN] ;
  FILE         *fp ;

  std::string cmdline = getAllInfo(argc, argv, "mri_cal_renormalize_gca");

  nargs = handleVersionOption(argc, argv, "mri_cal_renormalize_gca");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  setRandomSeed(-1L) ;
  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    ErrorExit
      (ERROR_BADPARM,
       "usage: %s [<options>] <longitudinal time point file> <in vol> <input atlas> <transform file> <output atlas> \n",
       Progname) ;
  in_fname = argv[2] ;
  gca_fname = argv[3] ;
  xform_fname = argv[4] ;
  out_fname = argv[5] ;

  printf("reading atlas from '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",Progname, gca_fname) ;
  GCAregularizeConditionalDensities(gca, .5) ;

  transform = TransformRead(xform_fname) ;
  if (transform == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s", Progname, xform_fname) ;

  start.reset() ;

  FileNamePath(argv[1], sdir) ;
  cp = strrchr(sdir, '/') ; 
  if (cp)
  {
    strcpy(base_name, cp+1) ;
    *cp = 0 ;  // remove last component of path, which is base subject name
  }
  ninputs = 0 ;
  fp = fopen(argv[1], "r") ;
  if (fp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read time point file %s", Progname, argv[1]) ;

  do
  {
    cp = fgetl(line, STRLEN-1, fp) ;
    if (cp != NULL && strlen(cp) > 0)
    {
      subjects[ninputs] = (char *)calloc(strlen(cp)+1, sizeof(char)) ;
      strcpy(subjects[ninputs], cp) ;
      ninputs++ ;
    }
  } while (cp != NULL && strlen(cp) > 0) ;
  fclose(fp) ;
  printf("processing %d timepoints in SUBJECTS_DIR %s...\n", ninputs, sdir) ;
  for (input = 0 ; input < ninputs ; input++)
  {
    printf("reading subject %d of %d: %s \n", input+1, ninputs, subjects[input]) ;
    if (longinput)
      sprintf(fname, "%s/%s.long.%s/mri/%s", sdir, subjects[input], base_name, in_fname) ;
    else
      sprintf(fname, "%s/%s/longtp/%s/%s", sdir, base_name, subjects[input], in_fname) ;

    mri_tmp = MRIread(fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read input MR volume from %s",
                Progname, fname) ;
    MRImakePositive(mri_tmp, mri_tmp) ;
    if (input == 0)
    {
      mri_in =
        MRIallocSequence(mri_tmp->width, mri_tmp->height, mri_tmp->depth,
                         mri_tmp->type, ninputs) ;
      if (!mri_in)
        ErrorExit(ERROR_NOMEMORY,
                  "%s: could not allocate input volume %dx%dx%dx%d",
                  mri_tmp->width,mri_tmp->height,mri_tmp->depth,ninputs) ;
      MRIcopyHeader(mri_tmp, mri_in) ;
    }

    if (mask_fname)
    {
      int i ;
      MRI *mri_mask ;

      mri_mask = MRIread(mask_fname) ;
      if (!mri_mask)
        ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                  Progname, mask_fname) ;

      for (i = 1 ; i < WM_MIN_VAL ; i++)
        MRIreplaceValues(mri_mask, mri_mask, i, 0) ;
      MRImask(mri_tmp, mri_mask, mri_tmp, 0, 0) ;
      MRIfree(&mri_mask) ;
    }
    MRIcopyFrame(mri_tmp, mri_in, 0, input) ;
    MRIfree(&mri_tmp) ;
  }
  MRIaddCommandLine(mri_in, cmdline) ;



  FileNameOnly(out_fname, line) ;
  FileNameRemoveExtension(line, line) ;
  strcpy(base_name, line) ;
  TransformInvert(transform, mri_in) ;
  GCAmapRenormalizeWithAlignmentLongitudinal
      (gca,
       mri_in,
       transform,
       NULL,
       base_name,
       NULL,
       0) ;

  printf("writing updated atlas to %s\n", out_fname) ;
  GCAwrite(gca, out_fname) ;
  MRIfree(&mri_in) ;

  printf("freeing GCA...") ;
  if (gca)
    GCAfree(&gca) ;
  printf("done.\n") ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("normalization took %d minutes and %d seconds.\n",
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
  StrUpper(option) ;
  if (!strcmp(option, "-HELP")||!strcmp(option, "-USAGE"))
    usage_exit(0) ;
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else if (!strcmp(option, "DEBUG_NODE"))
  {
    Ggca_x = atoi(argv[2]) ;
    Ggca_y = atoi(argv[3]) ;
    Ggca_z = atoi(argv[4]) ;
    printf("debugging node (%d, %d, %d)\n", Ggca_x, Ggca_y, Ggca_z) ;
    nargs = 3 ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else if (!strcmp(option, "LONGINPUT"))
  {
    longinput = 1;
    printf("reading inputs from longitudinal directories\n") ; 
  }
  else switch (*option)
  {
  case 'W':
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'H':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(int code)
{
  //  outputHelp(Progname);
#ifdef GREGT
  printf("usage: %s [<options>] <in volume> <atlas> <transform> "
         "<normalized atlas>\n\n",
         Progname) ;
  printf("\t<atlas>                      path/to/some.gca file "
         "(or 'noatlas') \n");
  printf("\t<transform>                  ex. transforms/talairach.lta "
         "(or 'noxform') \n");
  printf("\noptions:\n");
  printf("\t-seg <filename>              aseg file, to help normalization\n");
  printf("\t-sigma <bias sigma>          smoothing sigma for bias field if control points specified (def=4)\n");
  printf("\t-fsamples <filename>         write control points to filename\n");
  printf("\t-nsamples <filename>         write transformed "
         "normalization control points to filename\n");
  printf("\t-mask <mri_vol>              use mri_vol to mask input\n");
  printf("\t-f <filename>                define control points "
         "from filename\n");
  printf("\t-fonly <filename>            only use control points "
         "from filename\n");
  printf("\t-longinput                   load inputs from <tp>.long.<base> dirs, "
         "instead of <base>/longtp subdirs.\n");
  printf("\t-diag <filename>             write to log file\n");
  printf("\t-debug_voxel <x> <y> <z>     debug voxel\n");
  printf("\t-debug_node <x> <y> <z>      debug node\n");
  printf("\t-tr <float n>                set TR in msec\n");
  printf("\t-te <float n>                set TE in msec\n");
  printf("\t-alpha <float n>             set alpha in radians\n");
  printf("\t-example <mri_vol> <segmentation> use T1 (mri_vol) "
         "and segmentation as example\n");
  printf("\t-novar                       do not use variance estimates\n");
  printf("\t-renorm <mri_vol>            renormalize using "
         "predicted intensity values in mri_vol\n");
  printf("\t-flash                       use FLASH forward model "
         "to predict intensity values\n");
  printf("\t-prior <float t>             use prior threshold t "
         "(default=.6)\n");
  printf("\t-w                           write normalized volume "
         "each nregion iteration to norm(n).mgh(see -n)\n");
  printf("\t-n <int n>                   use n regions/struct "
         "for normalization\n");
  printf("\t-v <int n>                   does nothing as far "
         "as i can tell, but an option\n");
  printf("\t-p <float p>                 use top p percent(default=.25) "
         "white matter points as control points\n");
#endif
  exit(code) ;
}
