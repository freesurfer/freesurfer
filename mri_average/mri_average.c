#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static MRI *align_with_average(MRI *mri_src, MRI *mri_avg) ;

char *Progname ;
static int align = 0 ;
static int conform = 0 ;
static MORPH_PARMS  parms ;

static void usage_exit(int code) ;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs, i ;
  MRI    *mri_src, *mri_avg = NULL, *mri_tmp ;
  char   *in_fname, *out_fname ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  parms.dt = 1e-7 ;
  parms.tol = 1e-4 ;
  parms.momentum = 0.8 ;
  parms.niterations = 50 ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  out_fname = argv[argc-1] ;

  for (i = 1 ; i < argc-1 ; i++)
  {
    in_fname = argv[i] ;
    fprintf(stderr, "reading %s...\n", in_fname) ;
    mri_src = MRIread(in_fname) ;
    if (!mri_src)
      ErrorExit(Gerror, "%s: could not read MR volume %s",Progname,in_fname);
		if (conform)
		{
			MRI *mri_tmp ;

			mri_tmp = MRIconform(mri_src) ;
			MRIfree(&mri_src) ;
			mri_src = mri_tmp ;
		}
		else
		{
			mri_src->xsize = mri_src->ysize=mri_src->zsize = mri_src->thick = 1.0f ;
			mri_src->imnr0 = 1 ; mri_src->imnr1 = mri_src->depth ;
		}
    if (align && mri_avg)  /* don't align the first time */
    {
      mri_tmp = align_with_average(mri_src, mri_avg) ;
      MRIfree(&mri_src) ;
      mri_src = mri_tmp ;
    }

    mri_avg = MRIaverage(mri_src, i-1, mri_avg) ;
    MRIfree(&mri_src) ;
  }
  fprintf(stderr, "writing to %s..\n", out_fname) ;
  MRIwrite(mri_avg, out_fname) ;
  MRIfree(&mri_avg) ;
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
  if (!stricmp(option, "dt"))
  {
    parms.dt = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using dt = %e\n", parms.dt) ;
  }
  else if (!stricmp(option, "tol"))
  {
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using tol = %e\n", parms.tol) ;
  }
  else if (!stricmp(option, "conform"))
  {
		conform = 1 ;
    fprintf(stderr, "interpolating and padding input volumes to 256^3 "
						"isotopic voxels\n") ;
  }
  else switch (toupper(*option))
  {
  case 'M':
    parms.momentum = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using momentum = %2.3f\n", parms.momentum) ;
    break ;
  case 'A':
    align = 1 ;
    fprintf(stderr, "aligning volumes before averaging...\n") ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
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
usage_exit(int code)
{
  printf("usage: %s [options] <volume> ... <output volume>\n", Progname) ;
  printf("\t-a    linearly align input volumes before averaging\n") ;
  exit(code) ;
}

static MRI *
align_with_average(MRI *mri_src, MRI *mri_avg)
{
  MRI          *mri_aligned, *mri_in_red, *mri_ref_red ;

  parms.lta = LTAalloc(1, NULL) ;
  parms.in_np.neck_x0 = parms.in_np.neck_y0 = parms.in_np.neck_z0 = 1000 ;
  parms.in_np.neck_dx = parms.in_np.neck_dy = parms.in_np.neck_dz = 1.0 ;
#if 1
  mri_in_red = MRIreduceByte(mri_src, NULL) ;
  mri_ref_red = MRIreduceByte(mri_avg,NULL);
  parms.mri_ref = mri_ref_red ;
  parms.mri_in = mri_in_red ;  /* for diagnostics */
  fprintf(stderr,"computing %d linear transformation%s...\n",parms.lta->num_xforms, parms.lta->num_xforms>1?"s":"");
  MRIlinearAlign(mri_in_red, mri_ref_red, &parms) ;
d263 13

  mri_aligned = align_pca(mri_src, mri_avg) ;
  if (Gdiag & DIAG_WRITE)
    MRIwriteImageViews(mri_aligned, "after_pca", 400) ;

  fprintf(stderr, "aligning volume with average...\n") ;
#define WINDOW_IMAGES 0
#if WINDOW_IMAGES
  mri_in_windowed = MRIwindow(mri_aligned, NULL, WINDOW_HANNING,127,127,127,100.0f);
  mri_ref_windowed = MRIwindow(mri_avg,NULL,WINDOW_HANNING,127,127,127,100.0f);
#else
  mri_in_windowed = MRIcopy(mri_aligned, NULL) ;
  mri_ref_windowed = MRIcopy(mri_avg, NULL) ;
#endif
  mri_aligned = MRIlinearTransform(mri_src, NULL, parms.lta->xforms[0].m_L) ;

  MRIfree(&mri_in_red) ; MRIfree(&mri_ref_red) ;
#else
  parms.mri_ref = mri_src ;
  parms.mriNULL) ;
  mri_re /* for diagnostics */
  fprintf(stderr,"computing %d linear transformation%s...\n",parms.lta->num_xforms, parms.lta->num_xforms>1?"s":"");

  MRIlinearAlign(mri_src, mri_avg, &parms) ;

  mri_aligned = MRIlinearTransform(mri_src, NULL, parms.lta->xforms[0].m_L) ;

#endif
  return(mri_aligned) ;
}

