/**
 * @file  mris_smooth.c
 * @brief iterative averaging of vertex positions to smooth a surface.
 *
 * See (Fischl et al, NeuroImage, 1999)
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:34 $
 *    $Revision: 1.28 $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"

static char vcid[] =
  "$Id: mris_smooth.c,v 1.28 2011/03/02 00:04:34 nicks Exp $";

int main(int argc, char *argv[]) ;

static int count_big_curvatures(MRI_SURFACE *mris, double thresh) ;
static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int which_norm = NORM_MEAN ;
static int normalize_flag = 0 ;
static char curvature_fname[STRLEN] = "curv" ;
static char area_fname[STRLEN] = "area" ;
static int nbrs = 2 ;
static int normalize_area = 0 ;
static int navgs = 10 ;
static int niterations = 10 ;
static int rescale = 0 ;
static int write_iterations = 0 ;
static double l_spring = 1.0 ;
static float momentum = 0.0 ;
static int no_write = 0 ;

// -g 20 8 works well for hippo
static double gaussian_norm = 0 ;
static int gaussian_avgs = 0 ;
int MRIShistoThresholdGaussianCurvatureToMarked(MRI_SURFACE *mris, double pct) ;
int MRISthresholdPrincipalCurvatures(MRI_SURFACE *mris, double thresh) ;
int MRISthresholdGaussianCurvatureToMarked(MRI_SURFACE *mris, double low_thresh, double hi_thresh);

int
main(int argc, char *argv[])
{
  char               **av, *in_fname, *out_fname, fname[STRLEN], path[STRLEN] ;
  int                ac, nargs, start_t, pass ;
  MRI_SURFACE        *mris ;

  char cmdline[CMD_LINE_LEN] ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_smooth.c,v 1.28 2011/03/02 00:04:34 nicks Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_smooth.c,v 1.28 2011/03/02 00:04:34 nicks Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    print_help() ;
  }

  in_fname = argv[1] ;
  out_fname = argv[2] ;
  FileNamePath(out_fname, path) ;

  mris = MRISfastRead(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;

  MRISaddCommandLine(mris, cmdline) ;
  MRISremoveTriangleLinks(mris) ;
  fprintf(stderr, "smoothing surface tessellation for %d iterations...\n",
          niterations);

  MRIScomputeMetricProperties(mris) ;
  MRISstoreMetricProperties(mris) ;
  MRISsetNeighborhoodSize(mris, nbrs) ;
#define DT 0.5
  if (gaussian_norm > 0)
  {
    int i, done, start_avgs = gaussian_avgs, j ;

    done = 0;
    start_t = 0 ;
    pass = 0 ;
    do
    {
      for (i = start_t ; i < niterations+start_t ; i++)
      {
        MRIScomputeMetricProperties(mris) ;
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        for (j = 0 ; j < 5 ; j++)
        {
          MRISaverageVertexPositions(mris, 2) ; // turn flat spikes into tubular ones
          MRIScomputeMetricProperties(mris) ;
          MRIScomputeSecondFundamentalForm(mris) ;
          MRIShistoThresholdGaussianCurvatureToMarked(mris, (float)(mris->nvertices-20)/mris->nvertices) ;
        }
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
        MRIScomputeMetricProperties(mris) ;
        MRISsmoothSurfaceNormals(mris, gaussian_avgs) ;
        MRISclearMarks(mris) ;
        MRISthresholdGaussianCurvatureToMarked(mris, 10, 50);
        MRIScomputeSecondFundamentalForm(mris) ;
        MRIShistoThresholdGaussianCurvatureToMarked(mris, (float)(mris->nvertices-20)/mris->nvertices) ;
        MRISthresholdGaussianCurvatureToMarked(mris, 10, 50);
        if ((write_iterations > 0) && ((i % write_iterations) == 0))
        {
          char fname[STRLEN] ;

          sprintf(fname, "%s%04d", out_fname, i) ;
          printf("writing snapshot to %s...\n", fname) ;
          MRISwrite(mris, fname) ;
          if (Gdiag & DIAG_WRITE)
          {
            MRISuseGaussianCurvature(mris) ;
            sprintf(fname, "%s_K%04d", out_fname, i) ;
            printf("writing curvature to %s...\n", fname) ;
            MRISwriteCurvature(mris, fname) ;
            sprintf(fname, "%s_marked%04d", out_fname, i) ;
            printf("writing marks to %s...\n", fname) ;
            MRISwriteMarked(mris, fname) ;
          }
        }
        for (j = 0 ; j <= 5*nint(1/DT) ; j++)
        {
          MRISmarkedSpringTerm(mris, l_spring) ;
          MRISaverageGradients(mris, gaussian_avgs) ;
          MRISmomentumTimeStep(mris, momentum, DT, 1, gaussian_avgs) ;
          MRISclearGradient(mris) ;
          MRIScomputeMetricProperties(mris) ;
          MRISsmoothSurfaceNormals(mris, gaussian_avgs) ;
          {
            int vno ;
            VERTEX *v ;

            for (vno = 0 ; vno < mris->nvertices ; vno++)
            {
              v = &mris->vertices[vno] ;
              if (v->marked > 0)
              {
                v->K = 1.0/(v->marked) ;
              }
              else
              {
                v->K = 0 ;
              }
            }
          }
        }
      }
      MRISclearGradient(mris) ;
      if (gaussian_avgs == 2)
      {
        if (pass++ > 4)
        {
          done = 1 ;
        }
        else
        {
          int num = count_big_curvatures(mris, 2) ;
          printf("------------------------------------------------------\n") ;
          printf("------------------------------------------------------\n") ;
          printf("------------------ pass %d (num=%d) ------------------\n",
                 pass, num) ;
          printf("------------------------------------------------------\n") ;
          printf("------------------------------------------------------\n") ;
          gaussian_avgs = start_avgs ;
        }
      }
      else
      {
        gaussian_avgs /= 2 ;
        if (done ==0)
        {
          printf("----------------- setting avgs to %d -----------------\n", gaussian_avgs) ;
        }
      }
      start_t = i ;
    }
    while (!done) ;

#if 0
    // more smoothing with principal curvatures
    gaussian_avgs = start_avgs ;
    printf("--------------------------------------------------------------------------\n") ;
    printf("--------------------------------------------------------------------------\n") ;
    printf("---------------------- starting threshold smoothing ----------------------\n") ;
    printf("--------------------------------------------------------------------------\n") ;
    printf("--------------------------------------------------------------------------\n") ;
    do
    {
      for (i = start_t ; i < niterations+start_t ; i++)
      {
        MRIScomputeMetricProperties(mris) ;
        MRIScomputeSecondFundamentalForm(mris) ;
        MRISsmoothSurfaceNormals(mris, 16) ;
#define KTHRESH 1.5  // everything with kmin less than this will not move
        MRISthresholdPrincipalCurvatures(mris, KTHRESH) ;
        MRISspringTermWithGaussianCurvature(mris, gaussian_norm, l_spring) ;
        MRISaverageGradients(mris, gaussian_avgs) ;
        MRISmomentumTimeStep(mris, 0, 0.1, 1, gaussian_avgs) ;
        MRISclearGradient(mris) ;
        if ((write_iterations > 0) && (((i+1) % write_iterations) == 0))
        {
          char fname[STRLEN] ;

          sprintf(fname, "%s%04d", out_fname, i+1) ;
          printf("writing snapshot to %s...\n", fname) ;
          MRISwrite(mris, fname) ;
          if (Gdiag & DIAG_WRITE/* && DIAG_VERBOSE_ON*/)
          {
            MRISuseGaussianCurvature(mris) ;
            sprintf(fname, "%s_K%04d", out_fname, i+1) ;
            printf("writing curvature to %s...\n", fname) ;
            MRISwriteCurvature(mris, fname) ;
          }
        }
      }
      MRISclearGradient(mris) ;
      done = (gaussian_avgs == 2) ;
      gaussian_avgs /= 2 ;
      if (done ==0)
      {
        printf("---------------------- setting avgs to %d ----------------------\n", gaussian_avgs) ;
      }
      start_t = i ;
    }
    while (!done) ;
#endif
  }
  else
  {
    MRISaverageVertexPositions(mris, niterations) ;
  }

  fprintf(stderr, "smoothing complete - recomputing first and second "
          "fundamental forms...\n") ;
  MRIScomputeMetricProperties(mris) ;

  if (rescale)
  {
    MRISscaleBrainArea(mris) ;
  }
  MRIScomputeSecondFundamentalForm(mris) ;
  MRISuseMeanCurvature(mris) ;
  MRISaverageCurvatures(mris, navgs) ;
  if (normalize_flag)
  {
    MRISnormalizeCurvature(mris, which_norm) ;
  }
  sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
          curvature_fname);
  if (no_write == 0)
  {
    fprintf(stderr, "writing smoothed curvature to %s/%s\n", path,fname) ;
    MRISwriteCurvature(mris, fname) ;
    sprintf(fname, "%s.%s", mris->hemisphere == LEFT_HEMISPHERE?"lh":"rh",
            area_fname);
    fprintf(stderr, "writing smoothed area to %s/%s\n", path, fname) ;
    MRISwriteArea(mris, fname) ;
  }

  if (Gdiag & DIAG_SHOW)
  {
    fprintf(stderr, "writing smoothed surface to %s\n", out_fname) ;
  }
  MRISwrite(mris, out_fname) ;
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
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size = %d\n", nbrs) ;
  }
  else if (!stricmp(option, "normalize"))
  {
    normalize_flag = 1 ;
  }
  else if (!stricmp(option, "nw"))
  {
    no_write = 1 ;
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "area"))
  {
    normalize_area = 1 ;
    printf("normalizing area after smoothing\n") ;
  }
  else switch (toupper(*option))
    {
    case 'M':
      momentum = atof(argv[2]) ;
      printf("using momentum = %2.2f\n", momentum) ;
      nargs = 1 ;
      break ;
    case 'G':
      gaussian_norm = atof(argv[2]) ;
      gaussian_avgs = atoi(argv[3]) ;
      printf("using Gaussian curvature smoothing with norm %2.2f with %d smooth steps\n",
             gaussian_norm, gaussian_avgs) ;
      nargs = 2 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      nargs = 1 ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      write_iterations = atoi(argv[2]) ;
      printf("writing out snapshots every %d iterations\n", write_iterations) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'R':
      rescale = 1 ;
      fprintf(stderr, "rescaling brain area after smoothing...\n") ;
      break ;
    case 'B':
      strcpy(area_fname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'C':
      strcpy(curvature_fname, argv[2]) ;
      nargs = 1 ;
      break ;
    case 'N':
      niterations = atoi(argv[2]) ;
      fprintf(stderr, "smoothing for %d iterations\n", niterations) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging curvature for %d iterations\n", navgs) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

#include "mris_smooth.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_smooth_help_xml,mris_smooth_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

#define NBINS 256
int
MRIShistoThresholdGaussianCurvatureToMarked(MRI_SURFACE *mris, double pct)
{
  HISTOGRAM  *h ;
  double     min_curv, max_curv, K, bin_size, total, mode = 0.0, mode_peak, mean, std, dmean,
                                                     dmin, dmax, dsigma ;
  int        vno, num, b, bin_no, bin_thresh = 0, vno_min, vno_max, skipped, nvertices ;
  VERTEX     *v ;

  dmean = MRIScomputeVertexSpacingStats(mris, &dsigma,&dmin, &dmax, &vno_min, &vno_max,
                                        CURRENT_VERTICES);
  min_curv = 1e8 ;
  max_curv = -min_curv ;

  // compute mean nbr spacing for each vertex so we can ignore tangles
  for (dmean = 0.0, num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    int    n ;

    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    for (v->d = 0.0, n = 0 ; n < v->vtotal ; n++)
    {
      v->d += v->dist[n] ;
    }
    v->d /= v->vtotal ;
    dmean += v->d ;
    num++ ;
  }

  dmean /= num ;
  mean = std = 0.0 ;
  for (skipped = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (v->d < 0.01*dmean)
    {
      skipped++ ;
      continue ;
    }
    K = MIN(fabs(v->k1), fabs(v->k2)) ;
    K = fabs(v->K) ;
    if (K < min_curv)
    {
      min_curv = K ;
    }
    if (K > max_curv)
    {
      max_curv = K ;
    }
    mean += K ;
    std += K*K ;
  }

  nvertices = mris->nvertices - skipped ;
  max_curv = MIN(max_curv, 1000) ;
  mean /= (float)nvertices ;
  std = sqrt(std/(float)nvertices - mean*mean) ;
  bin_size = (max_curv - min_curv + 1) / NBINS ;
  h = HISTOalloc(NBINS) ;
  h->bin_size = bin_size ;
  for (b = 0 ; b < NBINS ; b++)
  {
    h->bins[b] = (float)b*bin_size + min_curv ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->d < 0.01*dmean)
    {
      continue ;
    }
    K = MIN(fabs(v->k1), fabs(v->k2)) ;
    K = fabs(v->K) ;
    bin_no = (int)((float)(K - min_curv) / (float)bin_size) ;
    if (bin_no > NBINS-1 || bin_no < 0)
    {
      bin_no = NBINS-1 ;
    }
    h->counts[bin_no]++ ;
  }

  mode_peak = 0 ;
  for (total = 0, b = 0 ; b < NBINS-1 ; b++)
  {
    if (h->counts[b] > mode_peak)
    {
      mode_peak = h->counts[b];
      mode = h->bins[b] ;
    }
    total += h->counts[b] ;
    if (total/nvertices >= pct)
    {
      bin_thresh = b+1 ;
      break ;
    }
#if 0
    if (h->bins[b] > mean+10*std)
    {
      printf("setting threshold based on mean %2.2f + 4 * std %2.2f\n",
             mean, std) ;
      bin_thresh = b ;
      break ;
    }
#endif
  }

#if 0
  {
    int b1, b2, prev_count, next_count, bthresh ;

#define BWIN 3
    for (b = BWIN ; b < NBINS ; b++)
    {
      prev_count = next_count = 0 ;
      // if sum of previous BWIN bins is <= sum of next BWIN bins, set thresh here
      for (b1 = MAX(b-BWIN,0) ; b1 < b ; b1++)
      {
        prev_count += h->counts[b1] ;
      }
      for (b2 = b+1 ; b2 < MIN(NBINS, b+BWIN) ; b2++)
      {
        next_count += h->counts[b2] ;
      }
      if (prev_count-3 <= next_count)
      {
        bthresh = b ;
        break ;
      }
    }
    if (bthresh > bin_thresh)
    {
      bin_thresh = bthresh ;  // use more conservature of the two
    }
  }
#endif

  if (bin_thresh == 0)
  {
    return(0) ;
  }
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("mode at %2.2f, max %2.2f\n", mode, max_curv) ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (v->ripflag)
    {
      continue ;
    }
    K = MIN(fabs(v->k1), fabs(v->k2)) ;
    K = fabs(v->K) ;
    bin_no = (int)((float)(K - min_curv) / (float)bin_size) ;
    if (bin_no >= bin_thresh)
    {
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      v->marked = 1 ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->marked == 1)
    {
      int    n ;
      VERTEX *vn ;

      for (n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (v->v[n] == Gdiag_no)
        {
          DiagBreak() ;
        }
        if (vn->marked == 0)
        {
          vn->marked = 2 ;
        }
      }
    }
  }

  HISTOfree(&h) ;
  return(NO_ERROR) ;
}

int
MRISthresholdPrincipalCurvatures(MRI_SURFACE *mris, double thresh)
{
  int     vno, num ;
  VERTEX  *v ;
  double  K ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    K = fabs(v->K) ;
    if (fabs(v->k1) < thresh || fabs(v->k2) < thresh || (v->k1 > 0 && v->k2 > 0))
    {
      v->K = 0 ;
    }
    else
    {
      num++ ;
    }

  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->K >= thresh)
    {
      int n ;

      v->marked = 1 ;
      for (n = 0 ; n < v->vtotal ; n++)
      {
        mris->vertices[v->v[n]].marked = 2 ;
        mris->vertices[v->v[n]].K = 0.5 ;
      }
    }
  }

  printf("%d vertices over threshold\n", num) ;
  return(NO_ERROR) ;
}

static int
count_big_curvatures(MRI_SURFACE *mris, double thresh)
{
  int    num, vno ;
  VERTEX *v ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (fabs(v->k1) > thresh && fabs(v->k2) > thresh)
    {
      num++ ;
    }
  }
  return(num) ;
}



int
MRISthresholdGaussianCurvatureToMarked(MRI_SURFACE *mris, double low_thresh, double hi_thresh)
{
  int    vno ;
  VERTEX *v, *vn ;
  double K ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }
    if (v->ripflag || v->marked == 2) // ignore if it was the nbr of one
    {
      continue ;
    }
    K = fabs(v->k1 * v->k2) ;
    if (K < low_thresh)
    {
      int n ;
      v->marked = 0 ;

      for (n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked == 2)
        {
          vn->marked = 0 ;
        }
      }
    }
    else if (K > hi_thresh)
    {
      v->marked = 1 ;
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (v->marked == 1)
    {
      int n ;

      for (n = 0 ; n < v->vtotal ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (v->v[n] == Gdiag_no)
        {
          DiagBreak() ;
        }
        if (vn->marked == 0)
        {
          vn->marked = 2 ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

