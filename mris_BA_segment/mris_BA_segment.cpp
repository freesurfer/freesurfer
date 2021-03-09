/**
 * @brief program for segmenting a Brodmann area from MRI.
 *
 * This program uses laminar intensity profiles to segment Brodmann areas from
 * imaging data, using MT as the initial target.
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
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int make_profiles_zero_baseline(MRI *mri) ;
static int compute_spherical_distances_to_vertex(MRI_SURFACE *mris, int vno0) ;
static double compute_MT_log_likelihood(MRI_SURFACE *mris, MRI *mri_profiles, int vno0, double r,
                                        double exterior_mm) ;

const char *Progname ;

static int nbhd_size = 2 ;
static int navgs = 0 ;

static double MT_radius_mean = 7.79 ; // mm from Zilles data set
static double MT_radius_std = 1.13 ; // mm from Zilles data set
static LABEL *segment_MT(MRI_SURFACE *mris, MRI *mri, LABEL *lprior, 
                         double MT_radius_mean, double MT_radius_std) ;

static const char *white_name = "white" ;
int
main(int argc, char *argv[]) 
{
  char          **av, *out_name, *surf_name, *profile_name, *prior_name ;
  int           ac, nargs ;
  MRI_SURFACE   *mris ;
  LABEL         *lprior, *lout ;
  MRI           *mri ;

  nargs = handleVersionOption(argc, argv, "mris_BA_segment");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  surf_name = argv[1] ;
  profile_name = argv[2] ;
  prior_name = argv[3] ;
  out_name = argv[4] ;
  if (argc < 5)
    usage_exit() ;

  mris = MRISread(surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_name) ;
  if (MRISreadOriginalProperties(mris, white_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, white_name) ;
  MRISscaleBrain(mris, mris, sqrt(mris->orig_area / mris->total_area)) ;
  mris->radius = MRISaverageRadius(mris) ;
  MRIScenter(mris, mris) ;   // put center of sphere at (0,0,0)

  mri = MRIread(profile_name) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read profiles from %s",
              Progname, profile_name) ;
  make_profiles_zero_baseline(mri) ;

  lprior = LabelRead(NULL, prior_name) ;
  if (lprior == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read prior label from %s",
              Progname, prior_name) ;


  LabelCopyStatsToSurface(lprior, mris, VERTEX_VALS) ;
  if (navgs > 0)
  {
    printf("smoothing priors %d times (sigma=%2.2f)\n", navgs, sqrt(navgs*2/M_PI));
    MRISaverageVals(mris, navgs) ;
  }
  lout = segment_MT(mris, mri, lprior, MT_radius_mean, MT_radius_std) ;

  LabelWrite(lout, out_name) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_usage() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
  case 'A':
    navgs = atoi(argv[2]) ;
    printf("smoothing prior distribution with %d averages\n", navgs) ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    nargs = 1 ;
    break ;
  case 'N':
    nbhd_size = atoi(argv[2]) ;
    fprintf(stderr, "using neighborhood size=%d\n", nbhd_size) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
    "usage: %s [options] <surface> <profiles> <prior label> <output label>\n", Progname) ;
  print_help() ;
}

static void
print_help(void) {
  fprintf(stderr,
    "\nThis program segments a Brodmann area (MT currently) from a laminar"
    " intensity profile overlay\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#define EXTERIOR_MM 2   // mm for outside stuff
static LABEL *
segment_MT(MRI_SURFACE *mris, MRI *mri_profiles, LABEL *lprior, 
           double MT_radius_mean, double MT_radius_std)
{
  LABEL  *area ;
  int    *vertices, nvertices, vno, i, vno_best ;
  VERTEX *v ;
  double  ll, max_ll, rmin, rmax, r, rstep = 0.25, rbest, pradius ;

  for (vno = nvertices = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].val > 0.01)
      nvertices++ ;
  vertices = (int *)calloc(nvertices, sizeof(int)) ;
  if (vertices == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d vertex array", Progname, nvertices) ;
  for (vno = nvertices = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].val > 0.01)
      vertices[nvertices++] = vno ;

  printf("searching over %d vertices for MT center\n", nvertices) ;

  ll = max_ll = -1e8; ;
  rmin = MT_radius_mean-2*MT_radius_std ;
  rmax = MT_radius_mean+2*MT_radius_std;
  vno_best = -1 ; rbest = -1 ;
  for (i = 0 ; i < nvertices ; i++)
  {
    vno = vertices[i] ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    compute_spherical_distances_to_vertex(mris, vno) ;

    for (r = rmin ; r <= rmax ; r += rstep)
    {
      pradius = exp(-SQR(r-MT_radius_mean) / 2*SQR(MT_radius_std)) ;
      pradius *= (1.0 / (MT_radius_std * sqrt(2*M_PI))) ;
      ll = compute_MT_log_likelihood(mris, mri_profiles, vno, r, EXTERIOR_MM) ;
      ll += log(pradius) ;
      if (ll > max_ll)
      {
        rbest = r ;
        max_ll = ll ;
        vno_best = vno ;
        printf("new max %2.2f found at vno %d, r = %2.2f\n", ll, vno, r) ;
      }
    }
  }

  // build best label
  compute_spherical_distances_to_vertex(mris, vno_best) ;
  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->d < rbest)
      v->marked = 1 ;
  }
  area = LabelFromMarkedSurface(mris) ;
  {
    FILE *fp ;
    int  num, i ;
    double profile[100] ;

    memset(profile, 0, sizeof(profile)) ;
    for (num = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked == 0)
        continue ;
      for (i = 0 ; i < mri_profiles->nframes ; i++)
        profile[i] += MRIgetVoxVal(mri_profiles, vno, 0, 0, i) ;
      num++ ;
    }
    fp = fopen("profile.dat", "w") ;
    for (i = 0 ; i < mri_profiles->nframes ; i++)
    {
      profile[i] /= num ;
      fprintf(fp, "%2.3f\n", profile[i]) ;
    }
    fclose(fp) ;
  }
    
  return(area) ;
}


static double
compute_MT_log_likelihood(MRI_SURFACE *mris, MRI *mri_profiles, int vno0, double radius,
                          double exterior_mm)
{
  int  vno, i, nsamples = mri_profiles->nframes, num ;
  VERTEX *v, *v0 ;
  double val, ll, ll_vno, ll_exterior ;
  static double *u_interior = NULL, *v_interior, *v_exterior ;
 
  if (vno0 == Gdiag_no)
    DiagBreak() ;
  if (u_interior == NULL)
  {
    u_interior = (double *)calloc(nsamples, sizeof(double)) ;
    v_exterior = (double *)calloc(nsamples, sizeof(double)) ;
    v_interior = (double *)calloc(nsamples, sizeof(double)) ;
    if (u_interior == NULL || v_interior == NULL || v_exterior == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d elt vector", Progname,mri_profiles->nframes) ;
  }
  else
  {
    memset(u_interior, 0, nsamples*sizeof(double)) ;
    memset(v_interior, 0, nsamples*sizeof(double)) ;
    memset(v_exterior, 0, nsamples*sizeof(double)) ;
  }
  
  v0 = &mris->vertices[vno0] ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->d > radius)
      continue ;   // not in interior of putative MT
    for (i = 0 ; i < nsamples ; i++)
    {
      val = MRIgetVoxVal(mri_profiles, vno, 0,0, i) ;
      u_interior[i] += val ;
      v_interior[i] += val*val ;
    }
    num++ ;
  }
  for (i = 0 ; i < nsamples ; i++)
  {
    u_interior[i] /= num ;
    v_interior[i] = v_interior[i]/num - u_interior[i]*u_interior[i] ;
  }

  // compute log likelihood of interior 
  for (ll = 0.0,  vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->d > radius)
      continue ;   // not in interior of putative MT
    for (ll_vno = 0.0, i = 1 ; i < nsamples ; i++)
    {
      val = MRIgetVoxVal(mri_profiles, vno, 0,0, i) ;
      ll_vno += -SQR(u_interior[i]-val) / (2*v_interior[i]) ;
    }
    ll += ll_vno ;
    if (!std::isfinite(ll))
      DiagBreak() ;
  }
  if (num == 0)
    return(-1e10) ;
  ll /= num ;

  for (ll_exterior = 0.0, num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->d <= radius || v->d > radius+exterior_mm)
      continue ;   // not in immediate exterior of putative MT
    for (ll_vno = 0.0 ; i < nsamples ; i++)
    {
      val = MRIgetVoxVal(mri_profiles, vno, 0,0, i) ;
      v_exterior[i] += SQR(val-u_interior[i]) ;
      ll_vno += SQR(val-u_interior[i]) ;
    }
    ll_vno /= (nsamples-1) ;
#define EXT_K 10 // the bigger this is, the gentler the slope of the pdf is
    ll_exterior += log(.5 / (.5 + EXT_K * exp(-ll_vno))) ;
    num++ ;
  }

  ll_exterior /= num ;
  for (i = 1 ; i < nsamples ; i++)
    v_exterior[i] /= num ;
  return(ll + ll_exterior) ;
}
static int
compute_spherical_distances_to_vertex(MRI_SURFACE *mris, int vno0)
{
  int    vno ;
  VERTEX *v0, *v ;
  double x0, y0, z0, dot, x, y, z, theta, d ;

  v0 = &mris->vertices[vno0] ;
  x0 = v0->x ;  y0 = v0->y ;  z0 = v0->z ; 
  d = sqrt(x0*x0+y0*y0+z0*z0) ;
  x0 /= d ; y0 /= d ; z0 /= d ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    x = v->x ;  y = v->y ;  z = v->z ; 
    d = sqrt(x*x+y*y+z*z) ;
    x /= d ; y /= d ; z /= d ;
    dot = x*x0 + y*y0 + z*z0 ;
    theta = acos(dot) ;
    v->d = mris->radius * theta ;
  }
  return(NO_ERROR) ;
}

static int
make_profiles_zero_baseline(MRI *mri)
{
  int   vno, i ;
  double val0 ;

  for (vno = 0 ; vno < mri->width ; vno++)
  {
    val0 = MRIgetVoxVal(mri, vno, 0, 0, 0) ;
    for (i = 0 ; i < mri->nframes ; i++)
      MRIsetVoxVal(mri, vno, 0, 0, i, MRIgetVoxVal(mri, vno, 0, 0, i) - val0) ;
  }
  return(NO_ERROR) ;
}

