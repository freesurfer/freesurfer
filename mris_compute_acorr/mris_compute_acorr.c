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
#include "macros.h"
#include "fio.h"
#include "mrishash.h"

static char vcid[] = "$Id: mris_compute_acorr.c,v 1.1 2000/04/06 17:45:59 fischl Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

#define BIN_SIZE  1
#define MAX_DIST  20

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static double *MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris, 
                                                  float bin_size, 
                                                  float max_dist, int *pn) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;


/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  MRI_SURFACE  *mris ;
  char         **av, *curv_name, *surf_name, *hemi, fname[STRLEN],
               *cp, *subject_name, subjects_dir[STRLEN] ;
  int          ac, nargs, n ;
  double       *acorr ;

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

  /* subject_name hemi surface curvature */
  if (argc < 4)
    usage_exit() ;
  
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
              Progname) ;
  strcpy(subjects_dir, cp) ;

  subject_name = argv[1] ;
  hemi = argv[2] ;
  surf_name = argv[3] ;
  curv_name = argv[4] ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjects_dir,subject_name,hemi,surf_name);
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  if (strchr(curv_name, '/') != NULL)
    strcpy(fname, curv_name) ;  /* full path specified */
  else
    sprintf(fname,"%s/%s/surf/%s.%s",subjects_dir,subject_name,hemi,curv_name);
  if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
    ErrorExit(Gerror, "%s: could no read curvature file %s",Progname,fname) ;

  acorr = MRIScomputeCurvatureAutocorrelation(mris, BIN_SIZE, MAX_DIST, &n) ;
  free(acorr) ;
  MRISfree(&mris) ;

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
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, "usage: %s [options] <subject> <hemi> <surf> <curv>\n",
          Progname) ;
  fprintf(stderr, "where surf must be a spherical surface suitable for "
          "computing geodesics\n") ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will compute the autocorrelation function of"
          " a curvature file\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

static double *
MRIScomputeCurvatureAutocorrelation(MRI_SURFACE *mris, float bin_size, 
                                    float max_dist, int *pn)
{
  MHT     *mht ;
  int     vno, n, i, index, *counts, nbins ;
  VERTEX  *v, *vn ;
  MHBT    *bucket ;
  MHB     *bin ;
  float   x, y, z, radius, dist ;
  double  angle, circumference, *acorr ;
  VECTOR  *v1, *v2 ;

  fprintf(stderr, "building spatial LUT...\n") ;
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;
  radius = MRISaverageRadius(mris) ;
  circumference = M_PI * 2.0 * radius ;
  mht = MHTfillVertexTableRes(mris, NULL, CURRENT_VERTICES, 2*max_dist) ;

  nbins = max_dist/bin_size+1 ;
  acorr = (double *)calloc(nbins, sizeof(double)) ;
  counts = (int *)calloc(nbins, sizeof(int)) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (!(vno % 10000))
      fprintf(stderr, "%d of %d vertices processed\n", vno, mris->nvertices) ;
    x = v->x ; y = v->y ; z = v->z ;
    bucket = MHTgetBucket(mht, x, y, z) ;
    VECTOR_LOAD(v1, v->x, v->y, v->z) ;  /* radius vector */
    for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++)
    {
      n = bin->fno ; vn = &mris->vertices[n] ;
      VECTOR_LOAD(v2, vn->x, vn->y, vn->z) ;  /* radius vector */
      angle = fabs(Vector3Angle(v1, v2)) ;
#if 0
      xd = v->x - vn->x ; yd = v->y - vn->y ; zd = v->z - vn->z ;
      dist = sqrt(xd*xd + yd*yd + zd*zd) ;
#else
      dist = circumference * angle / (2.0 * M_PI) ;
#endif
      if (dist < max_dist)
      {
        index = (int)((float)dist/bin_size) ;
        counts[index]++ ;
        acorr[index] += (double)(v->curv * vn->curv) ;
      }
    }

  }
  
  MHTfree(&mht) ;

  for (i = 0 ; i < nbins ; i++)
  {
    if (counts[i])
    {
      acorr[i] /= (float)counts[i] ;
      printf("%2.4f  %2.4f  %d\n",
             (float)i*bin_size, acorr[i], counts[i]) ;
    }
#if 0
    else
      printf("0  0  0\n") ;
#endif
  }
  *pn = nbins ;
  free(counts) ;
  VectorFree(&v1) ; VectorFree(&v2) ;
  return(acorr) ;
}

