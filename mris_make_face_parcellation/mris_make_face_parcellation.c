/**
 * @file  mris_make_face_parcellation.c
 * @brief create a parcellation of equal areas
 *
 * make a parcellation where each unit is assigned based on the face 
 * it maps to in a specified spherical surface (usually an icosahedral one)
 * 
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2009/10/09 20:04:27 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2009,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "timer.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "mrishash.h"


static char vcid[] =
  "$Id: mris_make_face_parcellation.c,v 1.6 2009/10/09 20:04:27 fischl Exp $";

typedef struct
{
  int    write_iterations ;
  char   base_name[STRLEN] ;
  int    max_iterations ;
  double tol ;
  double l_markov ;
} PARMS ;

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int write_snapshot(MRI_SURFACE *mris, PARMS *parms, int n) ;
static double markov_energy(MRI_SURFACE *mris) ;
static int adjust_parcellation_boundaries(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, MRI *mri_cmatrix,
                                          PARMS *parms) ;

static int
update_parcellation_statistics(MRI_SURFACE *mris, int vno, int old_parcel, int new_parcel,
                               MRI *mri_cmatrix, MRI *mri_means, MRI *mri_vars, 
                               MRI *mri_stats);
char *Progname ;

static char *cmatrix_fname = NULL ;
static MRI  *mri_cmatrix ;
static PARMS  parms ;
int
main(int argc, char *argv[]) {
  char               **av, *in_fname, *ico_fname, *out_fname, path[STRLEN], ico_name[STRLEN] ;
  int                ac, nargs ;
  float              scale ;
  MRI_SURFACE        *mris, *mris_ico ;
  float              radius ;
  int                fno, r, g, b, vno, annot ;
  double             fdist ;
  char               cmdline[CMD_LINE_LEN] ;
  FACE               *face ;
  MHT                *mht ;
  VERTEX             *v ;
  int                msec, minutes, seconds ;
  struct timeb       start ;

  parms.max_iterations = 100 ;
  parms.tol = 1e-4 ;
  parms.l_markov = 50 ;
  TimerStart(&start) ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_make_face_parcellation.c,v 1.6 2009/10/09 20:04:27 fischl Exp $",
   "$Name:  $", cmdline);

  setRandomSeed(1L) ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_make_face_parcellation.c,v 1.6 2009/10/09 20:04:27 fischl Exp $",
     "$Name:  $");
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

  if (argc < 4)
    print_help() ;

  in_fname = argv[1] ;
  ico_fname = argv[2] ;
  out_fname = argv[3] ;
  FileNamePath(out_fname, path) ;
  FileNameOnly(ico_fname, ico_name) ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  mris_ico = MRISread(ico_fname) ;
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not read ico file %s",
              Progname, ico_fname) ;
  if (mris_ico->nfaces < 256)
    scale = 256 / mris_ico->nfaces ;
  else
    scale = 1 ;

  if (cmatrix_fname)  // read in a correlation matrix for refining the parcellation
  {
    mri_cmatrix = MRIread(cmatrix_fname) ;
    if (mri_cmatrix == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface %s\n", Progname, cmatrix_fname) ;
    
    if (mri_cmatrix->width != mris->nvertices || mri_cmatrix->nframes != mris->nvertices)
      ErrorExit(ERROR_BADFILE, "%s: cmatrix must be %d x 1 x 1 x %d",Progname,mris->nvertices,mris->nvertices);
  }

  MRISaddCommandLine(mris, cmdline) ;

  mris->ct = CTABalloc(mris_ico->nfaces) ;
  strcpy (mris->ct->fname, ico_fname);

  printf("parcellating hemisphere into %d units\n", mris_ico->nfaces) ;
  for (fno = 0 ; fno < mris_ico->nfaces ; fno++)
  {
    int f = nint(scale * fno), m ;

    // don't let the color be too close to 0 by scaling them
    r = (nint(scale*f) / (256*256)) ;
    g = (nint(scale*f) / 256) ;
    b = (nint(scale*f) % 256) ;
    m = fno % 10 ;

    // try to avoid having adjacent triangles with similar colors
    switch (m)
    {
    default:
    case 0: r += 128 ; break ;
    case 1: g += 128 ; break ;
    case 2: b += 128 ; break ;
    case 3: r = 255-r ; break ;
    case 4: g = 255-g ; break ;
    case 5: b = 255-b ; break ;
    case 6: b += 128 ; r = 255-r ; break ;
    case 7: g += 128 ; b = 255-b ;break ;
    case 8: g += 128 ; b = 255-b ;break ; r += 128 ;
    case 9: g += 128 ; b = 255-b ;break ; r = 255-r ;
    }

    if (r < 0)
      r = 0 ;
    if (g < 0)
      g = 0 ;
    if (b < 0)
      b = 0 ;
    r = r % 256 ; g = g % 256 ; b = b %256 ;
    if (r > 255 || g > 255 || b > 255)
      DiagBreak() ;
    sprintf (mris->ct->entries[fno]->name, "%s face %d", ico_name, fno);
    mris->ct->entries[fno]->ri = r ;
    mris->ct->entries[fno]->gi = g ;
    mris->ct->entries[fno]->bi = b ;
    mris->ct->entries[fno]->ai = 255;
    
    /* Now calculate the float versions. */
    mris->ct->entries[fno]->rf =
      (float)mris->ct->entries[fno]->ri / 255.0;
    mris->ct->entries[fno]->gf =
      (float)mris->ct->entries[fno]->gi / 255.0;
    mris->ct->entries[fno]->bf =
      (float)mris->ct->entries[fno]->bi / 255.0;
    mris->ct->entries[fno]->af =
      (float)mris->ct->entries[fno]->ai / 255.0;
  }

  radius = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris_ico, mris_ico, radius / mris_ico->radius) ;
  MRIScomputeMetricProperties(mris_ico) ;
  mht = MHTfillTableAtResolution(mris_ico, NULL, CURRENT_VERTICES, 1.0);
  //  mht = MHTfillTableAtResolution(mris_ico, NULL, CURRENT_VERTICES, mris_ico->avg_vertex_dist/10);

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if ((((vno % (mris->nvertices/10))) == 0) && DIAG_VERBOSE_ON)
      printf("%2.1f%% done\n", 100.0*(float)vno / mris->nvertices) ;
    v = &mris->vertices[vno] ;
    MHTfindClosestFaceGeneric(mht, mris, 
                              v->x, v->y, v->z, 
                              mris_ico->avg_vertex_dist/10, 
                              mris_ico->avg_vertex_dist/10, 
                              &face, &fno, &fdist) ;
    if (fno < 0)
    {
      fno = mhtBruteForceClosestFace(mris_ico, v->x, v->y, v->z, 
                                      CURRENT_VERTICES, NULL);    
      if (fno  < 0)
      {
        printf("warning: v %d not found in MHT\n", vno) ;
        continue ;
      }
    }
    CTABannotationAtIndex(mris->ct, fno, &annot);
    v->annotation = annot ;
    v->marked = fno ;
  }

  if (mri_cmatrix)
  {
    char fname[STRLEN], *cp ;
    FileNameOnly(out_fname, fname) ;
    strcpy(parms.base_name, fname) ;
    cp = strrchr(parms.base_name, '.') ;
    if (cp)
      *cp = 0 ; // take out extension

    adjust_parcellation_boundaries(mris, mris_ico, mri_cmatrix, &parms) ;
  }
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing annotation to %s\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "parcellation took %d minutes and %d seconds.\n", minutes, seconds) ;
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
    print_help() ;
  else if (!stricmp(option, "markov")){
    parms.l_markov = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting markov coefficient to %2.3f\n", parms.l_markov) ;
  }
  else if (!stricmp(option, "tol")){
    parms.tol = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting tol to %2.3f\n", parms.tol) ;
  }
  else if (!stricmp(option, "-version")){
    print_version() ;
  } else switch (toupper(*option)) {
  case 'C':
    cmatrix_fname = argv[2] ;
    nargs = 1 ;
    printf("reading correlation matrix from %s\n", cmatrix_fname) ;
    break ;
  case 'M':
    parms.max_iterations = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting max iterations to %d\n", parms.max_iterations) ;
    break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      nargs = 1 ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      parms.write_iterations = atoi(argv[2]) ;
      nargs = 1 ;
      printf("setting write iterations to %d\n", parms.write_iterations) ;
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
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <ico file> <output annot>\n"
          "example: %s lh.inflated $FREESURFER_HOME/lib/bem/ic3.tri "
          "./lh.ic3.annot",
          Progname, Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
    "\nThis generates a parcellation based on which icosahedral face each vertex maps to.\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
static double compute_parcellation_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, double *penergy, MRI *mri_means, MRI *mri_vars) ;
static int build_parcellation_border_permutation(MRI_SURFACE *mris, int *vertex_permutation, 
                                                 int *pnborder) ;
static int compute_parcellation_statistics(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_stats, MRI *mri_means, MRI *mri_vars);

static int
adjust_parcellation_boundaries(MRI_SURFACE *mris, MRI_SURFACE *mris_ico, MRI *mri_cmatrix,
                               PARMS *parms)
{
  MRI     *mri_stats, *mri_means, *mri_vars ;
  int     *vertex_permutation, parcel, nborder, done, nchanged, iter, index, *nbrs, nparcels,  
    min_parcel,vno, n, nframes;
  double  energy, min_energy = 1e10, last_energy, parc_energy ;
  VERTEX  *v, *vn ;

  nframes = mri_cmatrix->nframes ;
  nparcels = mris_ico->nfaces ;
  mri_stats = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, 3) ;
  mri_means = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, nframes) ;
  mri_vars = MRIallocSequence(nparcels, 1, 1, MRI_FLOAT, nframes) ;

  nbrs = (int *)calloc(nparcels, sizeof(int)) ;
  vertex_permutation = (int *)calloc(mris->nvertices, sizeof(int)) ;
  if (vertex_permutation == NULL)
    ErrorExit(ERROR_NOMEMORY, "adjust_parcellation_boundaries: could not allocated permutation");
  iter = done = 0 ; 
  do
  {
    if (parms->write_iterations > 0 && ((iter % parms->write_iterations) == 0))
      write_snapshot(mris, parms, iter) ;
    nchanged = 0 ;
    build_parcellation_border_permutation(mris, vertex_permutation, &nborder) ;

    compute_parcellation_statistics(mris, mri_cmatrix, mri_stats, mri_means, mri_vars) ;
    last_energy = energy = compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, mri_means, mri_vars) ;
    if (iter == 0)
      printf("iter %d: nchanged = %d, nborder = %d, energy = %2.3f (%2.3f)\n", iter, nchanged,nborder,energy,parc_energy);
    for (index = 0 ; index < nborder ; index++)
    {
      vno = vertex_permutation[index] ;
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      parcel = v->marked ;

      // build list of all neihboring parcels at this point
      memset(nbrs, 0, nparcels*sizeof(int)) ;
      min_parcel = parcel ; min_energy = energy ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (v->v[n] == Gdiag_no)
          DiagBreak() ;
        if (vn->marked != parcel && nbrs[vn->marked] != 1) // not already done
        {
          // try changing the parcel and see if it decreases the energy
          nbrs[vn->marked] = 1 ;  // only do it once
          update_parcellation_statistics(mris, vno, parcel, vn->marked,
                                         mri_cmatrix, mri_means, mri_vars, mri_stats);
          v->marked = vn->marked ;
          energy = compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, mri_means, mri_vars) ;
          update_parcellation_statistics(mris, vno, vn->marked, parcel,
                                         mri_cmatrix, mri_means, mri_vars, mri_stats);
          if (energy < min_energy)
          {
            min_energy = energy ;
            min_parcel = vn->marked ;
          }
          v->marked = parcel ;
        }
      }
      if (min_parcel != parcel)
      {
        int annot ;
        v->marked = min_parcel ;
        energy = min_energy ;
        nchanged++ ;
        CTABannotationAtIndex(mris->ct, v->marked, &annot);
        v->annotation = annot ;
        update_parcellation_statistics(mris, vno, parcel, min_parcel,
                                       mri_cmatrix, mri_means, mri_vars, mri_stats);
      }
      energy = compute_parcellation_energy(mris, mri_stats, parms, &parc_energy, mri_means, mri_vars) ;
      if (!FZERO(energy- min_energy))
        DiagBreak() ;
    }
    done = (nchanged == 0) || (iter++ > parms->max_iterations) ||
      ((last_energy-energy)/last_energy < parms->tol) ;
        
    printf("iter %d: nchanged = %d, nborder = %d, energy = %2.3f (%2.3f)\n", iter, nchanged,nborder,min_energy,parc_energy);
  }  while (!done) ;
  if (parms->write_iterations > 0 && ((iter % parms->write_iterations) == 0))
    write_snapshot(mris, parms, iter) ;
  MRIfree(&mri_stats) ; MRIfree(&mri_means) ; MRIfree(&mri_vars) ;
  free(nbrs) ; free(vertex_permutation) ;
  return(NO_ERROR) ;
}

#define PARCELLATION_EPSILON 0.1
static double
compute_parcellation_energy(MRI_SURFACE *mris, MRI *mri_stats, PARMS *parms, double *penergy, MRI *mri_means, MRI *mri_vars)
{
  double  energy, var_within, var_between, mean, mean_nbr, parc_energy ;
  int     **nbrs, parcel, vno, n, nparcels, n_nbrs ;
  VERTEX  *v, *vn ;

  nparcels = mri_stats->width ;
  nbrs = (int **)calloc(nparcels, sizeof(int *)) ;
  for (parcel = 0 ; parcel < nparcels ; parcel++)
    nbrs[parcel] = (int *)calloc(nparcels, sizeof(int)) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ; 
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue;
    parcel = v->marked ;
    for (n = 0 ; n < v->v[n] ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      nbrs[parcel][vn->marked] = 1 ;
    }
  }
  for (parcel = 0, energy = 0.0 ; parcel < nparcels ; parcel++)
  {
    var_within = MRIgetVoxVal(mri_stats, parcel, 0, 0, 1) ;
    var_between = 0.0 ;
    mean = MRIgetVoxVal(mri_stats, parcel, 0, 0, 0) ;
    for (n_nbrs = n = 0 ; n < nparcels ; n++)
    {
      if (n != parcel && nbrs[parcel][n] > 0)
      {
        n_nbrs++ ;
        mean_nbr = MRIgetVoxVal(mri_stats, n, 0, 0, 0) ;
        var_between += SQR(mean-mean_nbr) ;
      }
    }
    if (n_nbrs > 0)
      var_between /= n_nbrs ;
    energy += var_between / (var_within+PARCELLATION_EPSILON) ;
  }
  energy /= nparcels ;
  for (parcel = 0 ; parcel < nparcels ; parcel++)
    free(nbrs[parcel]) ;
  free(nbrs) ; 

  parc_energy = 10000*energy ;
  energy = parc_energy + parms->l_markov * markov_energy(mris) ;
  *penergy = parc_energy ;
  return(energy) ;
}

static int
build_parcellation_border_permutation(MRI_SURFACE *mris, int *vertex_permutation, int *pnborder)
{
  int    vno, nborder, n, tmp, index ;
  VERTEX *v, *vn ;

  for (nborder = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (v->v[n] == Gdiag_no)
        DiagBreak() ;
      if (vn->marked != v->marked)
        break ;
    }
    if (n < v->vnum)   // is a border vertex
      vertex_permutation[nborder++] = vno ;
  }

  for (n = 0 ; n < nborder ; n++)
  {
    index = randomNumber(0.0, (double)(nborder-0.0001)) ;
    tmp = vertex_permutation[n] ;
    vertex_permutation[n] = vertex_permutation[index] ;
    vertex_permutation[index] = tmp ;
  }

  *pnborder = nborder ;
  return(NO_ERROR) ;
}

/*
  frame 0 of mri_stats will be the mean for the ith parcellation unit
  frame 1 of mri_stats will be the variance for the ith parcellation unit
  frame 2 of mri_stats will be the # of vertices in the ith parcellation unit
*/
static int
compute_parcellation_statistics(MRI_SURFACE *mris, MRI *mri_cmatrix, MRI *mri_stats,
                                MRI *mri_means, MRI *mri_vars)
{
  
  int    vno, parcel, nparcels, frame, nframes ;
  VERTEX *v ;
  double  var, mean,val, dof, mean_total, var_total ;
  
  if (Gdiag_no >= 0) 
    printf("computing parcellation statistics\n") ;
  nframes = mri_cmatrix->nframes ;
  nparcels = mri_stats->width ;
  MRIclear(mri_stats) ; MRIclear(mri_means) ; MRIclear(mri_vars) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    dof = MRIgetVoxVal(mri_stats, v->marked, 0, 0, 2) ;
    MRIsetVoxVal(mri_stats, v->marked, 0, 0, 2, dof+1) ;
    if (v->marked == Gdiag_no)
      printf("adding %d to parcel %d, dofs = %d, ", vno, v->marked, (int)dof+1) ;

    for (frame = 0 ; frame < nframes ; frame++)
    {
      val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
      mean = MRIgetVoxVal(mri_means, v->marked, 0, 0, frame) ;
      var = MRIgetVoxVal(mri_vars, v->marked, 0, 0, frame) ;
      mean += val ;
      var += val*val ;
      MRIsetVoxVal(mri_means, v->marked, 0, 0, frame, mean) ;
      MRIsetVoxVal(mri_vars, v->marked, 0, 0, frame, var) ;
      if (v->marked == Gdiag_no && frame == 2)
        printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
    }
  }

  for (parcel = 0 ; parcel < nparcels ; parcel++)
  {
    dof = MRIgetVoxVal(mri_stats, parcel, 0, 0, 2) ;
    if (dof <= 0)
      continue ;
    for (mean_total = var_total = 0.0, frame = 0 ; frame < nframes ; frame++)
    {
      mean = MRIgetVoxVal(mri_means, parcel, 0, 0, frame) ;
      var = MRIgetVoxVal(mri_vars, parcel, 0, 0, frame) ;
      mean /= dof ;
      var = var / dof - mean*mean ;
      if (var < 0)
      {
        FILE *fp = fopen("val.log", "w") ;

        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          v = &mris->vertices[vno] ;
          if (v->marked != parcel)
            continue ;
          val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;
          fprintf(fp, "%f\n", val) ;
        }
        fclose(fp) ;
        DiagBreak() ;
      }
      //      MRIsetVoxVal(mri_means, parcel, 0, 0, frame, mean) ;
      //      MRIsetVoxVal(mri_vars, parcel, 0, 0, frame, var) ;
      mean_total += mean ; var_total += var ;
    }
    mean_total /= nframes ;  var_total /= nframes ;
    MRIsetVoxVal(mri_stats, parcel, 0,0, 0, mean_total) ;
    MRIsetVoxVal(mri_stats, parcel, 0,0, 1, var_total) ;
  }
  return(NO_ERROR) ;
}
static int
write_snapshot(MRI_SURFACE *mris, PARMS *parms, int n)
{
  char fname[STRLEN] ;

  sprintf(fname, "%s.%3.3d.annot", parms->base_name, n) ;
  printf("writing snapshot to %s\n", fname) ;
  MRISwriteAnnotation(mris, fname) ;
  return(NO_ERROR) ;
}

static double
markov_energy(MRI_SURFACE *mris)
{
  int    vno, nborders, n ;
  VERTEX *v, *vn ;
  double energy ;

  for (energy = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    for (nborders = n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked != v->marked)
        nborders++ ;
    }
    energy += exp((double)nborders/(double)v->vnum)-1 ;
  }
  energy /= mris->nvertices ;
  return(energy) ;
}

static int
update_parcellation_statistics(MRI_SURFACE *mris, int vno, int old_parcel, int new_parcel,
                               MRI *mri_cmatrix, MRI *mri_means, MRI *mri_vars, 
                               MRI *mri_stats)
{
  int    frame, dofs, nframes ;
  double mean, var, val, total_mean, total_var ;

  nframes = mri_means->nframes ;
  dofs = MRIgetVoxVal(mri_stats, old_parcel, 0, 0, 2) - 1 ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 2, dofs) ;
  if (old_parcel == Gdiag_no)
    printf("removing %d from parcel %d, dofs = %d, ", vno, old_parcel, dofs) ;

  for (total_mean = total_var = 0.0, frame = 0 ; frame < nframes ; frame++)
  {
    val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;

    mean = MRIgetVoxVal(mri_means, old_parcel, 0, 0, frame) ;
    mean -= val ;
    var = MRIgetVoxVal(mri_vars, old_parcel, 0, 0, frame) ;
    var -= (val*val) ;
    MRIsetVoxVal(mri_means, old_parcel, 0, 0, frame, mean) ;
    MRIsetVoxVal(mri_vars, old_parcel, 0, 0, frame, var) ;
    mean /= dofs ;
    var = var /  dofs - mean*mean ;
    if (var < 0)
      DiagBreak() ;
    total_mean += mean ;
    total_var += var ;
    if (old_parcel == Gdiag_no && frame == 2)
      printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
  }
  total_mean /= nframes ; total_var /= nframes ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 0, total_mean) ;
  MRIsetVoxVal(mri_stats, old_parcel, 0, 0, 1, total_var) ;
  
  dofs = MRIgetVoxVal(mri_stats, new_parcel, 0, 0, 2)+1 ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 2, dofs) ;
  if (new_parcel == Gdiag_no)
    printf("adding %d to parcel %d, dofs = %d, ", vno, new_parcel, dofs) ;
  for (total_mean = total_var = 0.0, frame = 0 ; frame < nframes ; frame++)
  {
    val = MRIgetVoxVal(mri_cmatrix, vno, 0, 0, frame) ;

    mean = MRIgetVoxVal(mri_means, new_parcel, 0, 0, frame) ;
    mean += val ;
    var = MRIgetVoxVal(mri_vars, new_parcel, 0, 0, frame) ;
    var += (val*val) ;
    MRIsetVoxVal(mri_means, new_parcel, 0, 0, frame, mean) ;
    MRIsetVoxVal(mri_vars, new_parcel, 0, 0, frame, var) ;
    mean /= dofs ;
    var = var /  dofs - mean*mean ;
    if (var < 0)
      DiagBreak() ;
    total_mean += mean ;
    total_var += var ;
    if (new_parcel == Gdiag_no && frame == 2)
      printf("\tval = %2.3f, mean = %2.3f, var = %2.5f\n", val, mean, var) ;
  }
  total_mean /= nframes ; total_var /= nframes ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 0, total_mean) ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 1, total_var) ;
  MRIsetVoxVal(mri_stats, new_parcel, 0, 0, 2, dofs) ;

  return(NO_ERROR) ;
}
