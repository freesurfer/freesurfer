#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief program for computing various curvature metrics of a surface
 *
 * program for computing various curvature metrics of a surface.
 * This includes gaussian, mean, 1st and 2nd principle ones, and others.
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
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static char output_type[STRLEN] = "" ;
static const char *suffix = "" ;
static int write_flag = 0 ;
static int nbrs = 2 ;
static double cthresh = -1.0;
static int navgs = 0 ;
static char *param_file = NULL ;
static int normalize = 0 ;
static int diff_flag = 0 ;
static int max_flag = 0 ;
static int min_flag = 0 ;
static int stretch_flag = 0 ;
static int patch_flag = 0 ;
static int neg_flag = 0 ;
static int param_no = 0 ;
static int normalize_param = 0 ;
static int ratio_flag = 0 ;
static int contrast_flag = 0 ;

#define MAX_NBHD_SIZE 5000
static int nbhd_size = 0 ;
static int nbrs_per_distance = 0 ;
static float max_mm = 0 ;

static int which_norm = NORM_MEAN;

int MRIScomputeNeighbors(MRI_SURFACE *mris, float max_mm)  ;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname,fname[STRLEN],hemi[10], path[STRLEN],
               name[STRLEN],*cp ;
  int          ac, nargs, nhandles,err ;
  MRI_SURFACE  *mris ;
  double       ici, fi, var ;

  nargs = handleVersionOption(argc, argv, "mris_curvature");
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

  if (argc < 2)
  {
    usage_exit() ;
  }

  in_fname = argv[1] ;

  FileNamePath(in_fname, path) ;
  FileNameOnly(in_fname, name) ;
  cp = strchr(name, '.') ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: could not scan hemisphere from '%s'",
              Progname, fname) ;
  strncpy(hemi, cp-2, 2) ;
  hemi[2] = 0 ;

  if (patch_flag)  /* read the orig surface, then the patch file */
  {
    int req = snprintf(fname, STRLEN, "%s/%s.orig", path, hemi) ;
    if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    mris = MRISfastRead(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
    if (Gdiag & DIAG_SHOW)
    {
      fprintf(stderr, "reading patch file %s...\n", in_fname) ;
    }
    if (MRISreadPatch(mris, in_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;

  }
  else     /* just read the surface normally */
  {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }

  MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;

  if (nbhd_size > 0)
  {
    MRISsampleAtEachDistance(mris, nbhd_size, nbrs_per_distance) ;
  }
  if (max_mm > 0)
  {
    float ratio ;

    MRISstoreMetricProperties(mris) ;
    if (MRISreadCanonicalCoordinates(mris, "sphere") != NO_ERROR)
    {
      ErrorExit(ERROR_NOFILE,
                "%s: could not read canonical coordinates from ?h.sphere",
                Progname);
    }

    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    ratio = mris->orig_area / M_PI * mris->radius * mris->radius * 4.0 ;
    ratio = mris->orig_area / mris->total_area ;
    MRISscaleBrain(mris, mris, sqrt(ratio)) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeNeighbors(mris, max_mm) ;
  }

  if (param_file)
  {
    MRI_SP *mrisp ;
    mrisp = MRISPread(param_file) ;
    if (normalize_param)
    {
      MRISnormalizeFromParameterization(mrisp, mris, param_no) ;
    }
    else
    {
      MRISfromParameterization(mrisp, mris, param_no) ;
    }
    MRISPfree(&mrisp) ;
    if (normalize)
    {
      MRISnormalizeCurvature(mris,which_norm) ;
    }
    int req = snprintf(fname, STRLEN, "%s/%s%s.param", path,name,suffix) ;
    if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "writing parameterized curvature to %s...\n", fname) ;
    MRISwriteCurvature(mris, fname) ;
    fprintf(stderr, "done.\n") ;
  }
  else
  {
    MRIScomputeSecondFundamentalFormThresholded(mris, cthresh) ;

    nhandles = nint(1.0 - mris->Ktotal / (4.0*M_PI)) ;
    fprintf(stderr, "total integrated curvature = %2.3f*4pi (%2.3f) --> "
            "%d handles\n", (float)(mris->Ktotal/(4.0f*M_PI)),
            (float)mris->Ktotal, nhandles) ;

#if 0
    fprintf(stderr, "0: k1 = %2.3f, k2 = %2.3f, H = %2.3f, K = %2.3f\n",
            mris->vertices[0].k1, mris->vertices[0].k2,
            mris->vertices[0].H, mris->vertices[0].K) ;
    fprintf(stderr, "0: vnum = %d, v2num = %d, total=%d, area=%2.3f\n",
            mris->vertices[0].vnum, mris->vertices[0].v2num,
            mris->vertices[0].vtotal,mris->vertices[0].area) ;
#endif
    MRIScomputeCurvatureIndices(mris, &ici, &fi);
    var = MRIStotalVariation(mris) ;
    fprintf(stderr,"ICI = %2.1f, FI = %2.1f, variation=%2.3f\n", ici, fi, var);

    if (diff_flag)
    {
      MRISuseCurvatureDifference(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      int req = snprintf(fname, STRLEN, "%s/%s%s.diff", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing curvature difference to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (ratio_flag)
    {
      MRISuseCurvatureRatio(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int req = snprintf(fname, STRLEN, "%s/%s%s.ratio", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "writing curvature ratio to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (contrast_flag)
    {
      MRISuseCurvatureContrast(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int req = snprintf(fname, STRLEN, "%s/%s%s.contrast", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing curvature contrast to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }
    if (neg_flag)
    {
      int neg ;
      if (mris->patch)
      {
        mris->status = MRIS_PLANE ;
      }
      MRIScomputeMetricProperties(mris) ;
      neg = MRIScountNegativeTriangles(mris) ;
      MRISuseNegCurvature(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      int req = snprintf(fname, STRLEN, "%s/%s%s.neg", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing negative vertex curvature to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "%d negative triangles\n", neg) ;
      fprintf(stderr, "done.\n") ;
      {
        int    vno, fno ;
        FACE   *f ;
        for (vno = 0 ; vno < mris->nvertices ; vno++)
        {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
          VERTEX          const * const v  = &mris->vertices         [vno];
          if (v->ripflag)
          {
            continue ;
          }
          neg = 0 ;
          for (fno = 0 ; fno < vt->num ; fno++)
          {
            f = &mris->faces[vt->f[fno]] ;
            if (f->area < 0.0f)
            {
              neg = 1 ;
            }
          }
          if (neg)
          {
            fprintf(stdout, "%d\n", vno) ;
          }
        }
      }
    }

    if (max_flag)
    {
      MRISuseCurvatureMax(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int req = snprintf(fname, STRLEN, "%s/%s%s.max", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing curvature maxima to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (min_flag)
    {
      MRISuseCurvatureMin(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int req = snprintf(fname, STRLEN, "%s/%s%s.min", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing curvature minima to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (stretch_flag)
    {
      MRISreadOriginalProperties(mris, NULL) ;
      MRISuseCurvatureStretch(mris) ;
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int req = snprintf(fname, STRLEN, "%s/%s%s.stretch", path,name,suffix) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
      fprintf(stderr, "writing curvature stretch to %s...\n", fname) ;
      MRISwriteCurvature(mris, fname) ;
      fprintf(stderr, "done.\n") ;
    }

    if (write_flag)
    {
      MRISuseGaussianCurvature(mris) ;
      if (cthresh > 0)
      {
        MRIShistoThresholdCurvature(mris, cthresh) ;
      }
      MRISaverageCurvatures(mris, navgs) ;
      int req = snprintf(fname, STRLEN, "%s/%s%s.K%s", path,name, suffix, output_type) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "writing Gaussian curvature to %s...\n", fname) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      int len_curv_name = strlen(name)+strlen(suffix)+strlen(output_type)+2;
      char curv_name[len_curv_name+1];
      snprintf(curv_name, len_curv_name+1, "%s%s.K%s", name, suffix,output_type) ;
      MRISwriteCurvature(mris, fname,curv_name) ;
      MRISwriteCurvature(mris, fname) ;
      MRISuseMeanCurvature(mris) ;
      if (cthresh > 0)
      {
        MRIShistoThresholdCurvature(mris, cthresh) ;
      }
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize)
      {
        MRISnormalizeCurvature(mris,which_norm) ;
      }
      req = snprintf(fname, STRLEN, "%s/%s%s.H%s", path,name, suffix,output_type) ;
      if( req >= STRLEN ) {
	    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "done.\nwriting mean curvature to %s...\n", fname) ;
      snprintf(curv_name, len_curv_name+1, "%s%s.H%s", name, suffix,output_type) ;
      printf("%s\n",curv_name);
      err = MRISwriteCurvature(mris, fname,curv_name) ;
      if(err) exit(1);

      // write k1 and k2
      //MRI *mritmp = MRIcopyMRIS(NULL, mris, 0, "k1");
      //MRIwrite(mritmp,"k1.mgz");
      //mritmp = MRIcopyMRIS(mritmp, mris, 0, "k2");
      //MRIwrite(mritmp,"k2.mgz");

    }
  }
  printf("mris_curvature done.\n") ;
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
  else if (!stricmp(option, "distances") || !stricmp(option, "vnum"))
  {
    nbhd_size = atoi(argv[2]) ;
    nbrs_per_distance = atoi(argv[3]) ;
    fprintf(stderr, "sampling %d neighbors out to a distance of %d mm\n",
            nbrs_per_distance, nbhd_size) ;
    nargs = 2 ;
  }
  else if (!stricmp(option, "distance"))
  {
    max_mm = atof(argv[2]) ;
    fprintf(stderr, "sampling neighbors out to %f mm\n", max_mm) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "diff"))
  {
    diff_flag = 1 ;
  }
  else if (!stricmp(option, "mgh") || !stricmp(option, "-mgh") ||
	   !stricmp(option, "mgz") || !stricmp(option, "-mgz"))
  {
    strcpy(output_type, ".mgz");
    printf("appending %s to output files to change type\n", output_type) ;
  }
  else if (!stricmp(option, "ratio") || !stricmp(option, "defect"))
  {
    ratio_flag = 1 ;
  }
  else if (!stricmp(option, "contrast"))
  {
    contrast_flag = 1 ;
  }
  else if (!stricmp(option, "suffix"))
  {
    suffix = argv[2] ;
    nargs = 1 ;
    printf("appending suffix %s to output names\n", suffix) ;
  }
  else if (!stricmp(option, "neg"))
  {
    neg_flag = 1 ;
  }
  else if (!stricmp(option, "max"))
  {
    max_flag = 1 ;
  }
  else if (!stricmp(option, "min"))
  {
    min_flag = 1 ;
  }
  else if (!stricmp(option, "stretch"))
  {
    stretch_flag = 1 ;
  }
  else if (!stricmp(option, "median"))
  {
    which_norm = NORM_MEDIAN ;
    printf("using median normalization for curvature\n") ;
  }
  else if (!stricmp(option, "thresh"))
  {
    cthresh = atof(argv[2]) ;
    nargs = 1 ;
    printf("thresholding curvature at %2.2f%% level\n",
           cthresh*100) ;
  }
  else if (!stricmp(option, "param"))
  {
    param_file = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using parameterization file %s\n", param_file) ;
  }
  else if(!stricmp(option,"curvs")||!stricmp(option,"H")||!stricmp(option,"K")||!stricmp(option,"k1")||!stricmp(option,"k2")||!stricmp(option,"k1k2"))
  {
    MRIS *surf = MRISread(argv[2]);
    if(!surf) exit(1);
    int err, nbrhdsize;
    sscanf(argv[3],"%d",&nbrhdsize);
    printf("nbrhdsize %d\n",nbrhdsize);
    MRISsetNeighborhoodSizeAndDist(surf, nbrhdsize) ;
    MRIScomputeSecondFundamentalFormDiscrete(surf,0); // consistent with mris_curvature_stats
    char tmpstr[2000];
    MRI *mri=NULL;
    if(!stricmp(option,"curvs")||!stricmp(option,"H")){
      sprintf(tmpstr,"%s.H.mgz",argv[4]);
      mri = MRIcopyMRIS(mri, surf, 0, "H");
      err = MRIwrite(mri,tmpstr);
      if(err)exit(err);
    }
    if(!stricmp(option,"curvs")||!stricmp(option,"K")){
      sprintf(tmpstr,"%s.K.mgz",argv[4]);
      mri = MRIcopyMRIS(mri, surf, 0, "K");
      err = MRIwrite(mri,tmpstr);
      if(err)exit(err);
    }
    if(!stricmp(option,"curvs")||!stricmp(option,"k1")||!stricmp(option,"k1k2")){
      sprintf(tmpstr,"%s.k1.mgz",argv[4]);
      mri = MRIcopyMRIS(mri, surf, 0, "k1");
      err = MRIwrite(mri,tmpstr);
      if(err)exit(err);
    }
    if(!stricmp(option,"curvs")||!stricmp(option,"k2")||!stricmp(option,"k1k2")){
      sprintf(tmpstr,"%s.k2.mgz",argv[4]);
      mri = MRIcopyMRIS(mri, surf, 0, "k2");
      err = MRIwrite(mri,tmpstr);
      if(err)exit(err);
    }
    exit(0);
  }
  else if (!stricmp(option, "nparam"))
  {
    char *cp ;
    cp = strchr(argv[2], '#') ;
    if (cp)   /* # explicitly given */
    {
      param_no = atoi(cp+1) ;
      *cp = 0 ;
    }
    else
    {
      param_no = 0 ;
    }
    normalize_param = 1 ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "nbrs"))
  {
    nbrs = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
  }
  else if (!stricmp(option, "seed"))
  {
    setRandomSeed(atol(argv[2])) ;
    fprintf(stderr,"setting seed for random number generator to %d\n",
            atoi(argv[2])) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'N':
      fprintf(stderr, "normalizing curvature values.\n") ;
      normalize = 1 ;
      break ;
    case 'P':
      patch_flag = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "averaging curvature patterns %d times.\n", navgs) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'W':
      write_flag = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_help() ;
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
  print_help() ;
  exit(1) ;
}

#include "mris_curvature.help.xml.h"
static void
print_help(void)
{
  outputHelpXml(mris_curvature_help_xml,
                mris_curvature_help_xml_len);
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


int
MRIScomputeNeighbors(MRI_SURFACE *mris, float max_mm)
{
  int    vno, n, vlist[MAX_NBHD_SIZE], done, found, m, nbhd, first =1 ;
  size_t nbrs=0;
  float  dist, dx, dy, dz ;


  MRISresetNeighborhoodSize(mris, -1) ;  /* back to max */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno] ;
    VERTEX          * const v  = &mris->vertices         [vno] ;
    if (v->ripflag)
    {
      continue ;
    }
    if (vno == Gdiag_no)
    {
      DiagBreak() ;
    }

    for (n = 0 ; n < vt->vtotal ; n++)
    {
      mris->vertices[vt->v[n]].marked = 1 ;
      vlist[n] = vt->v[n] ;
    }

    cheapAssert(vt->nsizeCur == vt->nsizeMax);

    nbhd = vt->nsizeCur+1 ;
    nbrs = vt->vtotal ;
    do
    {
      found = 0 ;
      for (n = 0 ; n < nbrs ; n++)
      {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[n]] ;
        for (m = 0 ; m < vnt->vnum ; m++)
        {
          VERTEX * const vn2 = &mris->vertices[vnt->v[m]] ;
          if (vn2->marked)  // already in the nbhd list
          {
            continue ;
          }
          dx = vn2->cx - v->cx ;
          dy = vn2->cy - v->cy ;
          dz = vn2->cz - v->cz ;
          dist = sqrt(dx*dx + dy*dy + dz*dz) ;
          if (dist < max_mm)
          {
            vlist[nbrs] = vnt->v[m] ;
            nbrs++ ;
            found++ ;
            vn2->marked = 1 ;
            if (nbrs >= MAX_NBHD_SIZE)
            {
              if (first)
              {
                printf("max nbrs %zu exceeded at vertex %d\n", nbrs, vno) ;
              }
              first = 0 ;
              break ;
            }
          }
        }
      }
      done = ((found == 0) || (nbrs >= MAX_NBHD_SIZE)) ;
    }
    while (!done) ;
    free(vt->v) ;
    vt->v = (int *)calloc(nbrs, sizeof(int)) ;
    if (vt->v == NULL)
      ErrorExit(ERROR_NOMEMORY,
                "%s: vno %d could not allocate %zu vertex array",
                Progname, vno, nbrs) ;
    memmove(vt->v, vlist, sizeof(vlist[0])*nbrs) ;
    vt->vtotal = nbrs ;
    for (n = 0 ; n < nbrs ; n++)
    {
      mris->vertices[vlist[n]].marked = 0 ;
    }
  }

  mrisCheckVertexFaceTopology(mris);
  
  return(NO_ERROR) ;
}
