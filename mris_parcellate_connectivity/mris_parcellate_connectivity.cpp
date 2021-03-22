/**
 * @brief tool for parcellating a cortical model into relatively uniform connectivity regions.
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
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"
#include "mrishash.h"
#include "icosahedron.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

static int nbrs = 1 ;
static int navgs = 0 ;
static int nsub = 2 ;

static int ico_no = 2 ;

int
main(int argc, char *argv[]) {
  char         **av, *surf_fname, *parcellation_fname, *cmatrix_fname ;
  int          ac, nargs;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI_SURFACE  *mris, *mris_ico ;
  MRI          *mri_cmatrix ;
  char         fname[STRLEN], *cp ;
  MHT          *mht ;
  ICP          *icp ;

  setRandomSeed(1L) ;

  nargs = handleVersionOption(argc, argv, "mris_parcellate_connectivity");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit(1) ;

  icp = ICPread(ico_no, 6) ;
  exit(0) ;
  surf_fname = argv[1] ;
  cmatrix_fname = argv[2] ;
  parcellation_fname = argv[3] ;

  mris = MRISread(surf_fname) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s\n", Progname, surf_fname) ;
  MRISresetNeighborhoodSize(mris, nbrs) ;
  mri_cmatrix = MRIread(cmatrix_fname) ;
  if (mri_cmatrix == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s\n", Progname, cmatrix_fname) ;

  if (mri_cmatrix->width != mris->nvertices || mri_cmatrix->nframes != mris->nvertices)
    ErrorExit(ERROR_BADFILE, "%s: cmatrix must be %d x 1 x 1 x %d",Progname,mris->nvertices,mris->nvertices);

  cp = getenv("FREESURFER_HOME") ;
  if (cp == NULL)
    ErrorExit(ERROR_NOFILE, "%s: FREESURFER_HOME must be defined in the env",
              Progname) ;
  sprintf(fname, "%s/lib/bem/ic%d.tri", cp, ico_no) ;
  mris_ico = ICOreadOverAlloc(fname, 100, 1.0f) ;
                //
                // This is overallocating by 100x 
                // The parameter used to be called pct_ but was multiplied by without dividing by 100
                // so I don't know if the 100 is deliberate or not.  Left unchanged to not change the behaviour.
                
  printf("ico %d read, %d faces, %d vertices (%d, %d max)\n",
         ico_no, mris_ico->nvertices, mris_ico->nfaces,
         mris_ico->max_vertices, mris_ico->max_faces) ;

  MRIScomputeMetricProperties(mris_ico) ;
  MRIScomputeMetricProperties(mris) ;
  printf("scaling ico by %2.3f to radius %2.3f\n", mris->radius/mris_ico->radius, mris->radius) ;
  MRISscaleBrain(mris_ico, mris_ico, mris->radius/mris_ico->radius) ;

  if (navgs > 0)
    MRISsmoothMRI(mris, mri_cmatrix, navgs, NULL, mri_cmatrix) ;

  MRISprojectOntoSphere(mris_ico, mris_ico, DEFAULT_RADIUS) ;
  MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;

  mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, 1.0) ;

  MRISwrite(mris_ico, "lh.ico") ; printf("dividing icosahedral edges...\n") ;

  // figure out which vertices are on boundary between regions
  {
    int     vno, fno, ico_down, n, bad_triangle, add_vertex, vno2 ;
    MRIS    *mris_ico_down, *mris_ico_base ;
    FACE    *f ;

    sprintf(fname, "%s/lib/bem/ic%d.tri", cp, ico_no-(nsub-1)) ;
    mris_ico_base = ICOread(fname) ;
    MRISclearMarks(mris_ico);
    for (vno = 0 ; vno < mris_ico_base->nvertices ; vno++)
    {
      mris_ico->vertices[vno].border = mris_ico->vertices[vno].marked = 1 ;
    }
    // cannot be part of a face that is *only* at this level
    for (ico_down = ico_no - (nsub-2) ; ico_down <= ico_no ; ico_down++) 
    {
      sprintf(fname, "%s/lib/bem/ic%d.tri", cp, ico_down) ;
      mris_ico_down = ICOread(fname) ;
      for (fno = 0 ; fno < mris_ico->nfaces ; fno++)
      {
        f = &mris_ico->faces[fno] ;
        if (fno == Gdiag_no)
          DiagBreak() ;
        for (bad_triangle = 1, n = 0 ; n < VERTICES_PER_FACE ; n++)
        {
          VERTEX const * const v = &mris_ico->vertices[f->v[n]] ;
          if (f->v[n] == Gdiag_no)
            DiagBreak() ;
          if (v->border == 1)
          {
            bad_triangle = 0 ;
            break ;
          }
        }

        if (bad_triangle)
          for (n = 0 ; n < VERTICES_PER_FACE ; n++)
          {
            if (f->v[n] == Gdiag_no)
              DiagBreak() ;
            VERTEX * const v = &mris_ico->vertices[f->v[n]] ;
            v->marked = 2 ; // part of a bad triangle
          }
      }
      for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
      {
        VERTEX * const v = &mris_ico->vertices[vno] ;
        if (vno == Gdiag_no)
          DiagBreak() ;
        if (v->marked == 0)
        {
          v->marked = 1 ;
          v->border = 1 ;
        }
        else if (v->marked == 2) // part of an higher-order triangle
          v->marked = 0 ;
      }
      MRISfree(&mris_ico_base) ; mris_ico_base = mris_ico_down ;
    }


    sprintf(fname, "%s/lib/bem/ic%d.tri", cp, ico_no-(nsub)) ;
    mris_ico_base = ICOread(fname) ;

    for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris_ico->vertices_topology[vno];
      VERTEX                * const v  = &mris_ico->vertices         [vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      if (v->border == 1)  // already added
        continue ;
      for (add_vertex = 0, n = 0 ; n < vt->vnum ; n++)
      {
        vno2 = vt->v[n] ;
        if (vno2 == Gdiag_no)
          DiagBreak() ;
        if (vno2 < mris_ico_base->nvertices)
          add_vertex = 1;
      }

      if (add_vertex)
        v->border = v->marked = 1 ;
    }
    MRISfree(&mris_ico_base) ; mris_ico_base = mris_ico_down ;

    MRISfree(&mris_ico_down) ;
    MRISwriteMarked(mris_ico, "lh.border.mgz") ;
  }
#if 0
  MRISdivideEdges(mris_ico, nsub) ;
#endif
  MRISprojectOntoSphere(mris_ico, mris_ico, DEFAULT_RADIUS) ;
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.sub%d.ico", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", nsub) ;
    MRISwrite(mris_ico, fname) ;
  }

  MHTfree(&mht) ;
#if 0
  {
    MRI *mri_laplacian = MRISlaplacian(mris, mri_cmatrix, 1.0, 3.0) ;
    printf("writing output to %s\n", parcellation_fname) ;
    MRIwrite(mri_laplacian, parcellation_fname) ;
    MRIfree(&mri_laplacian) ;
  }
#endif


  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "parcellation took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  switch (toupper(*option)) {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    printf("debugging vertex %d\n", Gdiag_no) ;
    break ;
  case 'S':
    nsub = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting # of subdivisions = %d (default = %d)\n", nsub, 3) ;
    break ;
  case 'N':
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting nbhd size to %d\n", nbrs) ;
    break ;
  case 'A':
    navgs = atoi(argv[2]) ;
    nargs = 1 ;
    printf("smoothing correlation matrix %d times\n", navgs) ;
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
usage_exit(int code) {
  printf("usage: %s [options] <input surface> <input correlations> <output parcellation>\n",
         Progname) ;
  printf("\tn <avgs> - iteratively smooth correlation matrix\n") ;
  exit(code) ;
}

double
MRIScomputeUniformityEnergy(MRI_SURFACE *mris)
{
  return(0.0) ;
}
