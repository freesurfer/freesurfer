
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"

static char vcid[] = "$Id: mris_average_curvature.c,v 1.3 1998/07/02 22:56:37 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int normalize_flag = 0 ;
static int condition_no = 0 ;
static int stat_flag = 0 ;
static char *output_surf_name = NULL ;
static float sigma = 0.0f ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, *out_fname, *surf_name, fname[200], *sdir, 
               *hemi ;
  int          ac, nargs, i ;
  MRI_SURFACE  *mris ;
  MRI_SP       *mrisp, *mrisp_total ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir)
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n", Progname);
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  in_fname = argv[1] ;
  hemi = argv[2] ;
  surf_name = argv[3] ;
  out_fname = argv[argc-1] ;

  mrisp_total = MRISPalloc(1, 3) ;
  for (i = 4 ; i < argc-1 ; i++)
  {
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, argv[i], hemi, surf_name) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
    if (MRISreadCurvatureFile(mris, in_fname) != NO_ERROR)
      ErrorExit(ERROR_BADFILE,"%s: could not read curvature file %s.\n",
                Progname, in_fname);
#if 0
    if (normalize_flag)
      MRISnormalizeCurvature(mris) ;
#endif
    mrisp = MRIStoParameterization(mris, NULL, 1, 0) ;
    if (!FZERO(sigma))
    {
      MRI_SP  *mrisp_blur ;

#if 0
      mrisp_blur = MRISPblur(mrisp, NULL, sigma, 0) ;
#else
      mrisp_blur = MRISPconvolveGaussian(mrisp, NULL, sigma,mris->radius,0);
#endif
      MRISPfree(&mrisp) ;
      mrisp = mrisp_blur ;
    }
    MRISPcombine(mrisp, mrisp_total, 0) ;
    MRISPfree(&mrisp) ;
    if (i < argc-2)
      MRISfree(&mris) ;
  }

  if (output_surf_name)
  {
    sprintf(fname, "%s/%s/surf/%s.%s", sdir,output_surf_name,hemi,surf_name);
    fprintf(stderr, "reading output surface %s...\n", fname) ;
    MRISfree(&mris) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  }
  MRISfromParameterization(mrisp_total, mris, 0) ;
  if (stat_flag)    /* write out summary statistics files */
  {
    int    vno ;
    VERTEX *v ;
    float  dof ;
    FILE   *fp ;

    sprintf(fname, "%s/sigavg%d-%s.w", out_fname, condition_no, hemi);
    fprintf(stderr, "writing output means to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;

    MRISfromParameterization(mrisp_total, mris, 1) ;
    
    /* change variances to squared standard errors */
    dof = *IMAGEFseq_pix(mrisp_total->Ip, 0, 0, 2) ;
    if (!FZERO(dof)) for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
          continue ;
        v->curv /= dof ;   /* turn it into a standard error */
      }

    sprintf(fname, "%s/sigvar%d-%s.w", out_fname, condition_no, hemi);
    fprintf(stderr, "writing output variances to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;

    /* write out dof file */
    sprintf(fname, "%s/sigavg%d.dof", out_fname, condition_no) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open dof file %s\n",
                Progname,fname);
    fprintf(stderr, "writing dof file %s\n", fname) ;
    fprintf(fp, "%d\n", (int)dof) ;
    fclose(fp) ;
  }
  else
  {
    if (Gdiag & DIAG_SHOW)
      fprintf(stderr,"writing blurred pattern to surface to %s\n",out_fname);
    MRISwriteCurvature(mris, out_fname) ;
  }

  MRISfree(&mris) ;
  MRISPfree(&mrisp_total) ;
  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
  case 'A':
    sigma = atof(argv[2]) ;
    fprintf(stderr, "blurring thickness measures with sigma=%2.3f\n",sigma);
    nargs = 1 ;
    break ;
  case 'O':
    output_surf_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "painting output onto subject %s.\n", output_surf_name);
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'S':   /* write out stats */
    stat_flag = 1 ;
    condition_no = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "writing out summary statistics as condition %d\n",
            condition_no) ;
    break ;
  case 'N':
    normalize_flag = 1 ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input curv. file> <hemi> <surface> <subject> ... "
          " <output curv file >\n", Progname) ;
  fprintf(stderr, "the output curvature file will be painted onto the last "
          "subject name specified\non the command line.\n"
          "if the -s flag is specified then the last parameter specifies\n"
          "the directory in which to write the statistical maps.\n"
          "if the -o flag is specified then it overrides the last subject\n"
          "as the output surface.\n") ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-s <cond #>     generate summary statistics and write\n"
                  "                them into sigavg<cond #>-<hemi>.w and\n"
                  "                sigvar<cond #>-<hemi>.w.\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

