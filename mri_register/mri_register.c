

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
static MRI  *register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg, MP *parms) ;
static unsigned char thresh_low = 90 ;
static unsigned char thresh_hi = 120 ;
char *FileNameRemoveExtension(char *in_fname, char *out_fname) ;

static int linear = 0 ;
static char *transform_fname = NULL ;

/* 
   command line consists of three inputs:

   argv[1]  - directory containing 'canonical' brain
   argv[2]  - directory containing brain to be registered
   argv[3]  - directory in which to write out registered brain.
*/

int
main(int argc, char *argv[])
{
  char         *ref_fname, *in_fname, *out_fname, fname[100], **av ;
  MRI          *mri_ref, *mri_in, *mri_reg ;
  int          ac, nargs ;
  FILE         *fp ;
  int          msec, minutes, seconds ;
  struct timeb start ;

  parms.niterations = 100 ;
  parms.l_intensity = 0.1 ;
  parms.levels = -1 ;   /* use default */
  parms.l_dist = 0.0f ;
  parms.l_area = 0.1 ;
  parms.l_narea = 1000 ;
  parms.dt = .5 ;
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
              "usage: %s <ref brain> <in brain> <output file name>\n",
              Progname) ;

  ref_fname = argv[1] ;
  in_fname = argv[2] ;
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

  if (linear || !parms.m_L)    /* find optimal linear transformation */
  {
    parms.mri_ref = mri_ref ; parms.mri_in = mri_in ;
    parms.dt = 5e-6 ;
    fprintf(stderr, "computing optimal linear transformation...\n") ;
    mri_reg = register_mri(mri_in, mri_ref, NULL, &parms) ;
    
    fprintf(stderr, "writing output transformation to %s...\n", out_fname) ;
    fp = fopen(out_fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open output file %s\n",
                Progname, out_fname) ;
    fprintf(fp, "#%s %s %s\n", Progname, ref_fname, in_fname) ;
    MatrixAsciiWriteInto(fp, parms.m_L) ;
    fclose(fp) ;
    parms.niterations = 15 ;
    parms.l_intensity = parms.l_dist = 1.0f ;
    parms.dt = .005 ;
  }

  if (!linear)   /* do the 3d morph */
  {
    MORPH_3D *m3d ;
    MRI      *mri_in_reduced, *mri_ref_reduced ;

    if (mri_in->width > 128)
      mri_in_reduced = MRIreduceByte(mri_in, NULL) ;
    else
      mri_in_reduced = mri_in ;
    if (mri_ref->width > 128)
      mri_ref_reduced = MRIreduceByte(mri_ref, NULL) ;
    else
      mri_ref_reduced = mri_ref ;
    
    m3d = MRI3Dmorph(mri_in_reduced, mri_ref_reduced, &parms) ;
    if (mri_ref_reduced != mri_ref)
      MRIfree(&mri_ref_reduced) ;
    if (mri_in_reduced != mri_in)
      MRIfree(&mri_in_reduced) ;
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr, "applying 3d morph...\n") ;
    mri_reg = MRIapply3DMorph(mri_in, mri_ref, m3d, NULL) ;

    fprintf(stderr, "writing transformed volume to %s...\n", out_fname) ;
    MRI3DmorphFree(&m3d) ;
    MRIwrite(mri_reg, out_fname) ;
  }
  if (mri_reg)
    MRIfree(&mri_reg) ;
  if (mri_ref)
    MRIfree(&mri_ref) ;
  if (mri_in)
    MRIfree(&mri_in) ;
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
      parms.m_L = MatrixAsciiRead(transform_fname, NULL) ;
      if (!parms.m_L)
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

static MRI *
register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg, MORPH_PARMS *parms)
{
  /*  parms->write_iterations = 0 ; */
  if (!parms->niterations)
    parms->niterations = 1000 ;
  parms->l_intensity = parms->l_area = parms->l_dist = 1.0f ;
  parms->m_L = MatrixIdentity(4, NULL) ;
  MRIremoveNeck(mri_ref, mri_ref, thresh_low, thresh_hi, parms, 1) ; 
  MRIremoveNeck(mri_in, mri_in, thresh_low, thresh_hi, parms, -1) ; 
  MRIlinearAlign(mri_in, mri_ref, parms) ;
  mri_reg = MRIlinearTransform(mri_in, mri_reg, parms->m_L) ;

  return(mri_reg) ;
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

