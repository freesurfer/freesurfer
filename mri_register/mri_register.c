

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
#include "transform.h"

char         *Progname ;
static MORPH_PARMS  parms ;

static int get_option(int argc, char *argv[]) ;
static LTA  *register_mri(MRI *mri_in, MRI *mri_ref) ;
static unsigned char thresh_low = 90 ;
static unsigned char thresh_hi = 120 ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;

static int linear = 0 ;
static char *transform_fname = NULL ;

/* 
   command line consists of three inputs:

   argv[1]  - input brain (to be registered)
   argv[2]  - 'canonical' brain (target of registration)
   argv[3]  - transform fle (output)
*/

int
main(int argc, char *argv[])
{
  char         ref_fname[100], *in_fname, *out_fname, fname[100], **av ;
  MRI          *mri_ref, *mri_in ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  MORPH_3D     *m3d ;
  MRI          *mri_in_reduced, *mri_ref_reduced/*, *mri_reg*/ ;

  parms.niterations = 100 ;
  parms.l_intensity = 0.1 ;
  parms.levels = -1 ;   /* use default */
  parms.l_dist = 1.0 ; parms.l_area = 1.0 ;
  parms.l_narea = 100 ;
  parms.dt = .5 ;
  parms.tol = 1.5e-2 ;
  parms.sigma = 0.0f ;
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

  if (argc < 4)
    ErrorExit(ERROR_BADPARM, 
              "usage: %s <in brain> <ref brain> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  strcpy(ref_fname, argv[2]) ;
  if (strchr(ref_fname, '#') == NULL)
    strcat(ref_fname, "#-2") ;   /* only read in 2 frames - I know, a hack */
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  /*
     check to make sure output directory exists, and, if not, create it.
     */
  fprintf(stderr, "reading '%s'...", ref_fname) ;
  fflush(stderr) ;
  mri_ref = MRIread(ref_fname) ;
  if (!mri_ref)
    ErrorExit(ERROR_NOFILE, "%s: could not open reference volume %s.\n",
              Progname, ref_fname) ;
#if 0
  if (mri_ref->nframes > 2)
    MRIfreeFrames(mri_ref, 2) ;
#endif
  fprintf(stderr, "done.\n") ;
  fflush(stderr) ;

  fprintf(stderr, "reading '%s'...", in_fname) ;
  fflush(stderr) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;
  fprintf(stderr, "done.\n") ;
  fflush(stderr) ;

  if (linear || !parms.lta)    /* find optimal linear transformation */
  {
    fprintf(stderr, "computing optimal linear transformation...\n") ;
    parms.lta = register_mri(mri_in, mri_ref) ;
    
    fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
    LTAwrite(parms.lta, out_fname) ;
  }

  if (mri_in->width > 128)
  {
    mri_in_reduced = MRIreduceByte(mri_in, NULL) ;
    MRIfree(&mri_in) ;
  }
  else
    mri_in_reduced = mri_in ;
  if (mri_ref->width > 128)
  {
    mri_ref_reduced = MRIreduceMeanAndStdByte(mri_ref, NULL) ;
    MRIfree(&mri_ref) ;
  }
  else
    mri_ref_reduced = mri_ref ;

  m3d = MRI3Dmorph(mri_in_reduced, mri_ref_reduced, &parms) ;
  MRIfree(&mri_ref_reduced) ;
  MRIfree(&mri_in_reduced) ;
#if 1
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing 3d morph transform to %s...\n", out_fname) ;
  MRI3Dwrite(m3d, out_fname) ;
#else
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "applying 3d morph...\n") ;
  mri_reg = MRIapply3DMorph(mri_in, mri_ref, m3d, NULL) ;
  
  fprintf(stderr, "writing transformed volume to %s...\n", out_fname) ;
  MRI3DmorphFree(&m3d) ;
  MRIwrite(mri_reg, out_fname) ;

  if (mri_reg)
    MRIfree(&mri_reg) ;
  if (mri_ref)
    MRIfree(&mri_ref) ;
  if (mri_in)
    MRIfree(&mri_in) ;
#endif
  LTAfree(&parms.lta) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "morphing took %d minutes and %d seconds.\n", 
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
  if (!stricmp(option, "TOL"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE"))
  {
    parms.l_dist = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_dist = %2.2f\n", parms.l_dist) ;
  }
  else if (!strcmp(option, "DT"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "dt = %2.2e\n", parms.dt) ;
  }
  else if (!strcmp(option, "AREA"))
  {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.2f\n", parms.l_area) ;
  }
  else if (!strcmp(option, "NAREA"))
  {
    parms.l_narea = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_narea = %2.2f\n", parms.l_narea) ;
  }
  else if (!strcmp(option, "LEVELS"))
  {
    parms.levels = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "levels = %d\n", parms.levels) ;
  }
  else if (!strcmp(option, "INTENSITY") || !strcmp(option, "CORR"))
  {
    parms.l_intensity = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_intensity = %2.2f\n", parms.l_intensity) ;
  }
  else switch (*option)
  {
  case 'S':
    parms.sigma = atof(argv[2]) ;
    fprintf(stderr, "using sigma=%2.3f as upper bound on blurring.\n", 
            parms.sigma) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [image file name]\n", argv[0]) ;
    exit(1) ;
    break ;
#if 0
  case 'T':
    thresh_low = atoi(argv[2]) ;
    thresh_hi = atoi(argv[3]) ;
    fprintf(stderr, "thresholds set to %d --> %d\n", thresh_low, thresh_hi) ;
    nargs = 2 ;
    break ;
#endif
  case 'T':
    {
      transform_fname = argv[2] ;
      fprintf(stderr, "reading optimal linear alignment from %s...\n",
              transform_fname) ;
      nargs = 1 ;
      parms.lta = LTAread(transform_fname) ;
      if (!parms.lta)
        ErrorExit(ERROR_NOFILE, "%s: could not open transform file %s\n",
                  Progname, transform_fname) ;
    }
    break ;
  case 'L':
    linear = 1 ;
    fprintf(stderr, "finding optimal linear alignment...\n") ;
    break ;
  case 'N':
    parms.niterations = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "niterations = %d\n", parms.niterations) ;
    break ;
  case 'W':
    parms.write_iterations = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "write iterations = %d\n", parms.write_iterations) ;
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "momentum = %2.2f\n", parms.momentum) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static LTA *
register_mri(MRI *mri_in, MRI *mri_ref)
{
  MORPH_PARMS parms ;

  memset(&parms, 0, sizeof(parms)) ;
  parms.mri_ref = mri_ref ; parms.mri_in = mri_in ;
  parms.dt = 5e-6 ; parms.tol = 1e-2 ; parms.niterations = 10 ;
  parms.l_intensity = 1.0f ;
  parms.lta = LTAalloc(1, mri_in) ;
  if (MRIremoveNeck(mri_ref, mri_ref, thresh_low, thresh_hi, &parms, 1) == 
      NULL)
    ErrorExit(ERROR_BADPARM, "%s: registration failed", Progname) ;
  if (MRIremoveNeck(mri_in, mri_in, thresh_low, thresh_hi, &parms, -1)  == 
      NULL)
    ErrorExit(ERROR_BADPARM, "%s: registration failed", Progname) ;

  MRIlinearAlign(mri_in, mri_ref, &parms) ;

  return(parms.lta) ;
}
char *
FileNameRemoveExtension(char *in_fname, char *out_fname)
{
  char *dot ;

  if (out_fname != in_fname)
    strcpy(out_fname, in_fname) ;
  dot = strrchr(out_fname, '.') ;
  if (dot)
    *dot = 0 ;
  return(out_fname) ;
}

