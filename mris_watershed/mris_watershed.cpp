/*
From Wikipedia:


S. Beucher and F. Meyer introduced an algorithmic inter-pixel implementation of the watershed method,[5] given the following procedure:

1. Label each minimum with a distinct label. Initialize a set S with the labeled nodes.

2. Extract from S a node x of minimal altitude F, that is to say F(x) = min{F(y)|y ∈ S}. Attribute the label of x to each non-labeled node y adjacent to x, and insert y in S.

3. Repeat Step 2 until S is empty.


[5] Serge Beucher and Fernand Meyer. The morphological approach to segmentation: the watershed transformation. In Mathematical Morphology in Image Processing (Ed. E. R. Dougherty), pages 433481 (1993).

*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "romp_support.h"

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri.h"
#include "version.h"
#include "colortab.h"
#include "fsinit.h"


const char *Progname ;
int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISfindMaxLabel(MRI_SURFACE *mris)  ;
int  MRISinitWatershed(MRI_SURFACE *mris, MRI *mri)  ;
int  MRISwatershed(MRI_SURFACE *mris, MRI *mri, int max_clusters, int merge_type) ;
static int MRISfindMostSimilarBasins(MRI_SURFACE *mris, MRI *mri, int *pb2) ;
int MRISmergeBasins(MRI_SURFACE *mris, int b1, int b2) ;

#define MERGE_MOST_SIMILAR 1
#define MERGE_SMALLEST     2

static int merge_type = MERGE_SMALLEST ; //MERGE_MOST_SIMILAR ;
static int nbrs = 3 ;
static int max_clusters = 60 ;

static LABEL *mask_label= NULL ;

int
main(int argc, char *argv[])
{
  MRI_SURFACE *mris ;
  MRI         *mri ;
  int         nargs, nlabels ;
  char        *out_fname ;

  nargs = handleVersionOption(argc, argv, "mris_watershed");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  setRandomSeed(-1L) ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  mris = MRISread(argv[1]) ; 
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, argv[1]) ;
  mri = MRIread(argv[2]) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read intensity volume %s", Progname, argv[2]) ;
  out_fname = argv[3] ;

  MRISresetNeighborhoodSize(mris, nbrs) ;
  if (mask_label)
  {
    int vno ;
    VERTEX *v ;
    
    LabelMark(mask_label, mris) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->marked == 0)
	MRIsetVoxVal(mri, vno, 0, 0, 0, 0);
    }
    LabelUnmark(mask_label, mris) ;
  }

  MRISinitWatershed(mris, mri) ;
  MRISwatershed(mris, mri, max_clusters, merge_type) ;
  nlabels = MRISfindMaxLabel(mris)+1 ;

  printf("writing to %s with %d labels\n", out_fname, nlabels) ;

  {
    int vno, annot ;
    VERTEX *v ;
    
    mris->ct = CTABalloc(nlabels) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
	continue ;
      CTABannotationAtIndex(mris->ct, v->annotation, &annot);
      v->annotation = annot ;
    }
  }

  MRISwriteAnnotation(mris, out_fname) ;
  return(0) ;
}
  
static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <input gradient field> <output annotation>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program computes the watershed transform on the surface of an intensity gradient\n"
          "and writes the resulting measurement into a .annot file\n"
          "<output file>.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "\t-M <max clusters>: set the number of clusters\n") ;
  fprintf(stderr, "\t-mask_label <label>: read in and mask the input volume that is not in the specified label\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "openmp")) {
    char str[STRLEN] ;
    sprintf(str, "OMP_NUM_THREADS=%d", atoi(argv[2]));
    putenv(str) ;
#ifdef HAVE_OPENMP
    omp_set_num_threads(atoi(argv[2]));
#else
    fprintf(stderr, "Warning - built without openmp support\n");
#endif
    nargs = 1 ;
    fprintf(stderr, "Setting %s\n", str) ;
  } else if (!stricmp(option, "mask_label"))
  {
    nargs = 1 ;
    printf("reading masking label from %s\n", argv[2]) ;
    mask_label = LabelRead(NULL,argv[2]) ;
    if (mask_label == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load label %s for masking\n", Progname, argv[2]) ;
  }
  else if (!stricmp(option, "dilate") || !stricmp(option, "label_dilate") ||
	   !stricmp(option, "dilate_label"))
  {
    nargs = 1 ;
//    label_dilate = atoi(argv[2]) ;
//    printf("dilating label %d times\n", label_dilate) ;
  }
  else if (!stricmp(option, "erode") || !stricmp(option, "label_erode") ||
	   !stricmp(option, "erode_label"))
  {
    nargs = 1 ;
//    label_erode = atoi(argv[2]) ;
//    printf("eroding label %d times\n", label_erode) ;
  }
  else switch (toupper(*option)) {
    case 'N':
      nbrs = atoi(argv[2]) ;
      fprintf(stderr, "using neighborhood size=%d\n", nbrs) ;
      nargs = 1 ;
      break ;
    case 'M':
      max_clusters = atoi(argv[2]) ;
      printf("using max clusters=%d\n", max_clusters) ;
      nargs = 1 ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
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

int 
MRISinitWatershed(MRI_SURFACE *mris, MRI *mri) 
{
  int     vno, local_min, n, label_no = 1 ;
  float   val0, val ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    val0 = MRIgetVoxVal(mri, vno, 0, 0, 0) ;
    for (local_min = 1, n = 0 ; n < vt->vtotal ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->ripflag)
	continue ;
      val = MRIgetVoxVal(mri, vt->v[n], 0, 0, 0) ;
      if (val <= val0)
      {
	local_min = 0 ;
	break ;
      }
    }
    if (local_min)
      v->annotation = label_no++ ;
    else
      v->annotation = 0 ;
  }
  return(label_no) ;
}

typedef struct
{
  int    vno ;
  int    label; 
  float  val ;
} VERTEX_VALUE ;

static int
compare_vertex_vals(const void *v1, const void *v2)
{
  VERTEX_VALUE *vv1, *vv2 ;

  vv1 = (VERTEX_VALUE *)v1 ;
  vv2 = (VERTEX_VALUE *)v2 ;
  if (vv1->val < vv2->val) 
    return(-1) ;
  if (vv1->val > vv2->val) 
    return(1) ;
  return(0) ;  // equal
}
static VERTEX_VALUE *
MRISgetSortedVertexValues(MRI_SURFACE *mris, MRI *mri, VERTEX_VALUE *vv, int *pnvert) 
{
  int    label, vno, n ;

  if (vv == NULL)
  {
    vv = (VERTEX_VALUE *)calloc(mris->nvertices, sizeof(VERTEX_VALUE)) ;
    if (vv == NULL)
      ErrorExit(ERROR_NOMEMORY, "MRISwatershed: could not allocated %d-elt val array",
		mris->nvertices) ;
  }

  label = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (v->annotation > 0)
    {
      for (n = 0 ; n < vt->vnum ; n++)
      {
	VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
	if (vn->annotation <= 0)  // only if at least one neighbor is unlabeled
	{
	  vv[label].vno = vno ;
	  vv[label].val = MRIgetVoxVal(mri, vno, 0, 0, 0) ;
	  vv[label].label = v->annotation ;
	  label++ ;
	  break ;
	}
      }
    }
  }
  qsort(vv, label, sizeof(VERTEX_VALUE),compare_vertex_vals);
  *pnvert = label ;
  return(vv) ;
}

static int
MRISfindMostSimilarBasin(MRI_SURFACE *mris, MRI *mri, int min_basin) 
{
  int    best_basin, vno, n, basin, *nbr_vertices, nbr_basin, max_basin ;
  double *avg_grad, min_grad ;

  if( mris->nvertices >= 0 ) {
    nbr_vertices = (int *)calloc(mris->nvertices, sizeof(*nbr_vertices)) ;
    avg_grad = (double *)calloc(mris->nvertices, sizeof(*avg_grad)) ;
  } else {
    ErrorExit(ERROR_BADPARM, "%s: mris->nvertices<0\n", Progname) ;
    abort();
  }

  max_basin = 0 ;
  for (vno = 0 ;  vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag || v->annotation != min_basin)
      continue ;
    for (n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      nbr_basin = vn->annotation ;
      if (vn->ripflag || nbr_basin == min_basin || nbr_basin == 0)
	continue ;
      nbr_vertices[nbr_basin]++ ;
      avg_grad[nbr_basin] += MRIgetVoxVal(mri, vt->v[n], 0, 0, 0) ;
      if (nbr_basin > max_basin)
	max_basin = nbr_basin ;
    }
  }

  min_grad = 1e10 ;
  best_basin = -1 ;
  for (basin = 1 ; basin <= max_basin ; basin++)
  {
    if (nbr_vertices[basin] <= 0)
      continue ;
    avg_grad[basin] /= nbr_vertices[basin] ;
    if (avg_grad[basin] < min_grad)
    {
      best_basin = basin ;
      min_grad = avg_grad[basin] ;
    }
  }

  free(nbr_vertices) ;
  free(avg_grad) ;
  return(best_basin) ;
}
static int
MRISfindMostSimilarBasins(MRI_SURFACE *mris, MRI *mri, int *pb2) 
{
  int    best_basin, vno, basin, *nbr_vertices, max_basin, nbasins, best_b1, best_b2, basin1 ;
  double *avg_grad, min_grad, best_grad ;

  nbasins = MRISfindMaxLabel(mris) + 1 ;
  nbr_vertices = (int *)calloc(nbasins, sizeof(*nbr_vertices)) ;
  avg_grad = (double *)calloc(nbasins, sizeof(*avg_grad)) ;
  best_grad = 1e20;
  best_b1 = best_b2 = best_basin = -1;
  for (basin1 = 1 ; basin1 < nbasins ; basin1++)
  {
    memset(nbr_vertices, 0, nbasins*sizeof(*nbr_vertices)) ;
    memset(avg_grad, 0, nbasins*sizeof(*avg_grad)) ;
    max_basin = 0 ;
// reductions for min and max aren't available in earlier openmp
    ROMP_PF_begin
#if defined(HAVE_OPENMP) && GCC_VERSION > 40408
    #pragma omp parallel for if_ROMP(experimental) reduction(max:max_basin)
#endif
    for (vno = 0 ;  vno < mris->nvertices ; vno++)
    {
      ROMP_PFLB_begin
      int    n, nbr_basin ;

      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno] ;
      VERTEX          const * const v  = &mris->vertices         [vno] ;
      if (v->ripflag || v->annotation != basin1)
	continue ;
      for (n = 0 ; n < vt->vnum ; n++)
      {
	VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
	nbr_basin = vn->annotation ;
	if (vn->ripflag || nbr_basin == basin1 || nbr_basin == 0)
	  continue ;
	nbr_vertices[nbr_basin]++ ;
	avg_grad[nbr_basin] += MRIgetVoxVal(mri, vt->v[n], 0, 0, 0) ;
	if (nbr_basin > max_basin)
	  max_basin = nbr_basin ;
      }
      ROMP_PFLB_end
    }
    ROMP_PF_end
    
    min_grad = 1e10 ;
    for (basin = 1 ; basin <= max_basin ; basin++)
    {
      if (nbr_vertices[basin] <= 0)
	continue ;
      avg_grad[basin] /= nbr_vertices[basin] ;
      if (avg_grad[basin] < min_grad)
      {
	best_basin = basin ;
	min_grad = avg_grad[basin] ;
      }
    }
    if (min_grad < best_grad)
    {
      best_grad = min_grad ;
      best_b1 = basin1 ;
      best_b2 = best_basin ;
    }
  }
  free(nbr_vertices) ;
  free(avg_grad) ;
  *pb2 = best_b2 ;
  return(best_b1) ;
}
int
MRISfindMaxLabel(MRI_SURFACE *mris) 
{
  int    vno, max_label ;
  VERTEX *v ;
  
  max_label = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->annotation > max_label)
      max_label = v->annotation ;
  }
  return(max_label) ;
}
static int
MRISmergeSmallestBasin(MRI_SURFACE *mris, MRI *mri) 
{
  int    nbasins, vno, basin, min_basin, neighbor_basin ;
  double *basin_area, min_basin_area ;
  VERTEX *v ;

  nbasins = MRISfindMaxLabel(mris)+1 ;
  basin_area = (double *)calloc(nbasins, sizeof(*basin_area)) ;
  for (vno = 0 ;  vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation == 0)
      continue ;
    basin_area[v->annotation] += v->area ;
  }
  min_basin = 1 ; min_basin_area = basin_area[min_basin] ;
  for (basin = 1 ; basin < nbasins ; basin++)
    if (basin_area[basin] < min_basin_area)
    {
      min_basin_area = basin_area[basin] ;
      min_basin = basin ;
    }

  neighbor_basin = MRISfindMostSimilarBasin(mris, mri, min_basin) ;
  MRISmergeBasins(mris, min_basin, neighbor_basin) ;
  free(basin_area) ;
  return(NO_ERROR) ;
}

static int
MRISmergeMostSimilarBasin(MRI_SURFACE *mris, MRI *mri) 
{
  int    min_basin, neighbor_basin ;

  neighbor_basin = MRISfindMostSimilarBasins(mris, mri, &min_basin) ;

  MRISmergeBasins(mris, min_basin, neighbor_basin) ;
  return(NO_ERROR) ;
}

int
MRISmergeBasins(MRI_SURFACE *mris, int b1, int b2)
{
  int    vno, new_basin, other_basin ;
  VERTEX *v ;

  new_basin = MIN(b1, b2) ;    // basin to merge into
  other_basin = MAX(b1, b2) ;  // basin to be eliminated
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->annotation == 2349)
      DiagBreak() ;
    if (v->ripflag || v->annotation < other_basin)
      continue ;
    if (v->annotation == other_basin)
      v->annotation = new_basin ;
    else
      v->annotation = v->annotation - 1 ;   // we removed one basin - compact the rest
  }
  return(NO_ERROR) ;
}
int
MRISwatershed(MRI_SURFACE *mris, MRI *mri, int max_clusters, int merge_type) 
{
  int           vno, nvert, niter ;
  VERTEX_VALUE *vv = NULL ;
  double       max_val ;

  printf("flooding surface...\n") ;
  niter = 0 ;
  do
  {
    vv = MRISgetSortedVertexValues(mris, mri, vv, &nvert) ;
    max_val = 0 ;
//ifdef HAVE_OPENMP
//pragma omp parallel for if_ROMP(experimental) reduction(max:max_val)
//endif
    for (vno = 0 ; vno < nvert ; vno++)
    {
      double       min_nbr_val ;
      int          n ;

      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vv[vno].vno] ;
      VERTEX          const * const v  = &mris->vertices         [vv[vno].vno] ;
      if (v->ripflag)
	continue ;

      if (vv[vno].val > max_val)
	max_val = vv[vno].val ;   // new flooding height
      min_nbr_val = max_val+1 ;
      for (n = 0 ; n < vt->vnum ; n++)
      {
	VERTEX * const vn = &mris->vertices[vt->v[n]] ;
	if (vn->annotation <= 0)  
	{
	  float val ;
	  vn->annotation = v->annotation ;
	  val = MRIgetVoxVal(mri, vt->v[n], 0, 0, 0) ;
	  if (val < min_nbr_val)
	    min_nbr_val = val ;
	}
      }
      if (min_nbr_val < max_val)  // have to redo sorting
	vno = nvert ;
//	break ;
    }
    if (++niter % 1000 == 0)
      printf("iter %3.3d: nvertices = %d\n", niter,nvert) ;
  } while (nvert > 0) ;
      
  free(vv) ;

  // now do merging
  printf("merging basins...\n") ;
  niter = 0 ;
  while (MRISfindMaxLabel(mris) > max_clusters)
  {
    if (++niter % 200 == 0)
      printf("iter %3.3d: nclusters = %d\n", niter, MRISfindMaxLabel(mris)) ;
    switch (merge_type)
    {
    case MERGE_SMALLEST:
      MRISmergeSmallestBasin(mris, mri) ;
      break ;
    case MERGE_MOST_SIMILAR:
      MRISmergeMostSimilarBasin(mris, mri) ;
      break ;
    default:
      break ;
    }
  }

  return(NO_ERROR) ;
}

