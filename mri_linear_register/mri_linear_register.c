

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

char         *Progname ;
static MORPH_PARMS  parms ;

static int get_option(int argc, char *argv[]) ;
static int register_mri(MRI *mri_in, MRI *mri_ref, MP *parms) ;
static unsigned char thresh_low = 90 ;
static unsigned char thresh_hi = 120 ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;

static int num_xforms = 1 ;
static int transform_loaded = 0 ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

int
main(int argc, char *argv[])
{
  char         ref_fname[100], *in_fname, *out_fname, fname[100], **av ;
  MRI          *mri_ref, *mri_in ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  parms.l_intensity = 1.0f ;
  parms.niterations = 100 ;
  parms.levels = -1 ;   /* use default */
  parms.dt = 1e-6 ;  /* was 5e-6 */
  parms.tol = INTEGRATION_TOL*5 ;

  parms.dt = 5e-6 ;  /* was 5e-6 */
  parms.tol = 1e-3 ;
  parms.momentum = 0.8 ;
  parms.niterations = 15 ;
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
              "usage: %s <in brain> <template> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  strcpy(ref_fname, argv[2]) ;
  if (strchr(ref_fname, '#') == NULL)
    strcat(ref_fname, "#-2") ;
  out_fname = argv[3] ;
  FileNameOnly(out_fname, fname) ;
  FileNameRemoveExtension(fname, fname) ;
  strcpy(parms.base_name, fname) ;
  Gdiag |= DIAG_WRITE ;
  fprintf(stderr, "logging results to %s.log\n", parms.base_name) ;

  TimerStart(&start) ;
  fprintf(stderr, "reading '%s'...", ref_fname) ;
  fflush(stderr) ;
  mri_ref = MRIread(ref_fname) ;
  if (!mri_ref)
    ErrorExit(ERROR_NOFILE, "%s: could not open reference volume %s.\n",
              Progname, ref_fname) ;
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
  if (!transform_loaded)   /* wasn't preloaded */
    parms.lta = LTAalloc(1, mri_in) ;

  register_mri(mri_in, mri_ref, &parms) ;
  
  fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
  LTAwrite(parms.lta, out_fname) ;
  if (mri_ref)
    MRIfree(&mri_ref) ;
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "registration took %d minutes and %d seconds.\n", 
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
  if (!strcmp(option, "DIST") || !strcmp(option, "DISTANCE"))
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
  else if (!strcmp(option, "TOL"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "tol = %2.2e\n", parms.tol) ;
  }
  else if (!strcmp(option, "NUM"))
  {
    num_xforms = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "finding a total of %d linear transforms\n", num_xforms) ;
  }
  else if (!strcmp(option, "AREA"))
  {
    parms.l_area = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_area = %2.2f\n", parms.l_area) ;
  }
  else if (!strcmp(option, "NLAREA"))
  {
    parms.l_nlarea = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "l_nlarea = %2.2f\n", parms.l_nlarea) ;
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
  case 'T':
    parms.lta = LTAread(argv[2]) ;
    if (!parms.lta)
      ErrorExit(ERROR_BADFILE, "%s: could not read transform file %s",
                Progname, argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using previously computed transform %s\n", argv[2]) ;
    transform_loaded = 1 ;
    break ;
  case 'S':
    parms.sigma = atof(argv[2]) ;
    fprintf(stderr, "using sigma=%2.3f as upper bound on blurring.\n", 
            parms.sigma) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s <in volume> <template volume> <output transform>\n", 
           argv[0]) ;
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

static int
register_mri(MRI *mri_in, MRI *mri_ref, MORPH_PARMS *parms)
{
  MRI  *mri_in_red, *mri_ref_red ;

  parms->mri_ref = mri_ref ; parms->mri_in = mri_in ;  /* for diagnostics */
  mri_in_red = MRIreduceByte(mri_in, NULL) ;
  mri_ref_red = MRIreduceMeanAndStdByte(mri_ref,NULL);

  /*  parms->write_iterations = 0 ; */
  if (!parms->niterations)
    parms->niterations = 1000 ;
  if (transform_loaded)  /* don't recompute rotation based on neck */
  {
    MRIremoveNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, NULL, 1) ; 
    MRIremoveNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, NULL, -1) ; 
  }
  else
  {
    MRIremoveNeck(mri_ref_red, mri_ref_red, thresh_low, thresh_hi, parms, 1) ; 
    MRIremoveNeck(mri_in_red, mri_in_red, thresh_low, thresh_hi, parms, -1) ; 
  }

  while (parms->lta->num_xforms < num_xforms)
    LTAdivide(parms->lta, mri_in_red) ;
  fprintf(stderr,"computing %d linear transformation%s...\n",
          parms->lta->num_xforms, parms->lta->num_xforms>1?"s":"");
  MRIlinearAlign(mri_in_red, mri_ref_red, parms) ;

  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;
#if 0
  mri_reg = MRIlinearTransform(mri_in_red, mri_reg, parms->m_L) ;

  return(mri_reg) ;
#else
  return(NO_ERROR) ;
#endif
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

