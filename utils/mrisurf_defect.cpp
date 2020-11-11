#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
#define COMPILING_MRISURF_METRIC_PROPERTIES_FRIEND

/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mrisurf_defect.h"

#include "mrisurf_metricProperties.h"
#include "mrisurf_base.h"

//==================================================================
// Utilities for editing a surface
// Some of these need to be refactored to move portions of them into mrisurf_topology.c
//
void mrisDivideEdge(MRIS * const mris, int const vno1, int const vno2)
{
  int const vnew_no = mrisDivideEdgeTopologically(mris, vno1, vno2);

  VERTEX const * const v1   = &mris->vertices[vno1];
  VERTEX const * const v2   = &mris->vertices[vno2];
  VERTEX       * const vnew = &mris->vertices[vnew_no];

  vnew->x = (v1->x + v2->x) / 2;
  vnew->y = (v1->y + v2->y) / 2;
  vnew->z = (v1->z + v2->z) / 2;
  vnew->tx = (v1->tx + v2->tx) / 2;
  vnew->ty = (v1->ty + v2->ty) / 2;
  vnew->tz = (v1->tz + v2->tz) / 2;

  vnew->infx = (v1->infx + v2->infx) / 2;
  vnew->infy = (v1->infy + v2->infy) / 2;
  vnew->infz = (v1->infz + v2->infz) / 2;

  vnew->pialx = (v1->pialx + v2->pialx) / 2;
  vnew->pialy = (v1->pialy + v2->pialy) / 2;
  vnew->pialz = (v1->pialz + v2->pialz) / 2;

  vnew->cx = (v1->cx + v2->cx) / 2;
  vnew->cy = (v1->cy + v2->cy) / 2;
  vnew->cz = (v1->cz + v2->cz) / 2;
  vnew->x = (v1->x + v2->x) / 2;
  vnew->y = (v1->y + v2->y) / 2;
  vnew->z = (v1->z + v2->z) / 2;
  vnew->odx = (v1->odx + v2->odx) / 2;
  vnew->ody = (v1->ody + v2->ody) / 2;
  vnew->odz = (v1->odz + v2->odz) / 2;
  vnew->val = (v1->val + v2->val) / 2;
  
  MRISsetOriginalXYZ(mris, vnew_no, 
    (v1->origx + v2->origx) / 2,
    (v1->origy + v2->origy) / 2,
    (v1->origz + v2->origz) / 2);
}


/*
  nsubs = 1 --> divide edge in half
        = 2 --> divide edge in half twice, add 3 vertices
        = 3 --> divide edge in half three times, add 7 vertices
*/
#define MAX_SURFACE_FACES 500000
int MRISdivideEdges(MRIS *mris, int nsubdivisions)
{
  int nadded, sub, nfaces, fno, nvertices, faces[MAX_SURFACE_FACES], index;
  FACE *f;

  for (nadded = sub = 0; sub < nsubdivisions; sub++) {
    nfaces = mris->nfaces;
    nvertices = mris->nvertices;  // before adding any
    for (fno = 0; fno < nfaces; fno++) {
      faces[fno] = fno;
    }
    for (fno = 0; fno < nfaces; fno++) {
      int tmp;

      index = (int)randomNumber(0.0, (double)(nfaces - 0.0001));
      tmp = faces[fno];
      if (faces[fno] == Gdiag_no || faces[index] == Gdiag_no) {
        DiagBreak();
      }
      faces[fno] = faces[index];
      faces[index] = tmp;
    }

    for (index = 0; index < nfaces; index++) {
      fno = faces[index];
      f = &mris->faces[fno];
      if (fno == Gdiag_no) {
        DiagBreak();
      }

      if (f->v[0] < nvertices && f->v[1] < nvertices) { mrisDivideEdge(mris, f->v[0], f->v[1]); nadded++; }
      if (f->v[0] < nvertices && f->v[2] < nvertices) { mrisDivideEdge(mris, f->v[0], f->v[2]); nadded++; }
      if (f->v[1] < nvertices && f->v[2] < nvertices) { mrisDivideEdge(mris, f->v[1], f->v[2]); nadded++; }
    }
  }

  if (Gdiag & DIAG_SHOW && nadded > 0) {
    fprintf(stdout,
            "MRISdivideEdges(%d): %d vertices added: # of vertices=%d, # of faces=%d.\n",
            nsubdivisions,
            nadded,
            mris->nvertices,
            mris->nfaces);
#if 0
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
#endif
  }
  return (nadded);
}


int MRISdivideLongEdges(MRIS *mris, double thresh)
{
  int nadded = 0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    float const x = v->x;
    float const y = v->y;
    float const z = v->z;

    /*
      only add vertices if average neighbor vector is in normal direction, 
      that is, if the region is concave or sulcal.
    */
    int n;
    for (n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      double dist = sqrt(SQR(vn->x - x) + SQR(vn->y - y) + SQR(vn->z - z));
      if (dist > thresh) {
        mrisDivideEdge(mris, vno, vt->v[n]);
        nadded++;
      }
    }
  }

  if (Gdiag & DIAG_SHOW && nadded > 0) {
    fprintf(stdout,
            "%2.2f mm: %d vertices added: # of vertices=%d, # of faces=%d.\n",
            thresh,
            nadded,
            mris->nvertices,
            mris->nfaces);
#if 0
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
#endif
  }
  
  return (nadded);
}


//==================================================================
// Saving and restoring portions of the surface
// This holds a strong hint to what fields of the surface the rest of this code should be allowed to modify!
//
static void mrisFreeDefectVertexState(DEFECT_VERTEX_STATE *dvs)
{
  int i;
  for (i = 0; i < dvs->nvertices; i++) {
    freeAndNULL(dvs->vs[i].v);
    freeAndNULL(dvs->vs[i].f);
    freeAndNULL(dvs->vs[i].n);
  }
  freeAndNULL(dvs->vs);
  freeAndNULL(dvs);
}


static DEFECT_VERTEX_STATE* mrisRecordVertexState(MRIS const * const mris, DEFECT * const defect, int * const vertex_trans)
{
  DEFECT_VERTEX_STATE* const dvs = (DVS *)calloc(1, sizeof(DVS));
  if (!dvs) ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate dvs");

  dvs->defect       = defect;
  dvs->vertex_trans = vertex_trans;

/* in theory, the convex hull vertices should not be changed */
/* however, numerical errors might generate some unexpected 'bugs' */
/* including the convex hull vertices in DVS  limits these errors */

  // following code assumes these are zeroed
  //
#if MRIS_FIX_TOPOLOGY_ERROR_MODE
  dvs->nvertices = defect->nvertices + defect->nchull;
  dvs->vs = (VS *)calloc(dvs->nvertices, sizeof(VS));
#else
  dvs->nvertices = defect->nvertices + defect->nborder;
  dvs->vs = (VS *)calloc(dvs->nvertices, sizeof(VS));
#endif
  if (!dvs->vs) ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate %d dvs->vs", dvs->nvertices);

  /* keep the # of faces before retessellation */
  dvs->nfaces = mris->nfaces;

  /* keep the vertex state */
  int n;
  for (n = 0; n < defect->nvertices; n++) {
    dvs->vs[n].vno = vertex_trans[defect->vertices[n]];
  }
#if MRIS_FIX_TOPOLOGY_ERROR_MODE
  for (n = 0; n < defect->nchull; n++) {
    dvs->vs[defect->nvertices + n].vno = vertex_trans[defect->chull[n]];
  }
#else
  for (n = 0; n < defect->nborder; n++) {
    dvs->vs[defect->nvertices + n].vno = vertex_trans[defect->border[n]];
  }
#endif

  int i;
  for (i = 0; i < dvs->nvertices; i++) {
    VERTEX_STATE * const vs = &dvs->vs[i];
    int const vno = vs->vno;

    if (vno < 0) continue;

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (v->ripflag) {
      vs->vtotal = -1;
    } else {

      // save the topology
      //
      vs->nsizeCur = vt->nsizeCur;
      vs->nsizeMax = vt->nsizeMax;
      vs->vtotal   = vt->vtotal;
      vs->vnum     = vt->vnum;
      vs->v2num    = vt->v2num;
      vs->v3num    = vt->v3num;

      int const vsize = mrisVertexVSize(mris, vno);   // vtotal is based on nsizeCur, but may be bigger
      if (vsize) {                                    // and we need to restore up to nsizeMax
        int* vv = (int *)malloc(vs->vtotal*sizeof(int));
        if (!vv) ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate %dth array of %d elts", i, vsize);
        memcpy(vv,vt->v,vsize*sizeof(int));
        vs->v = vv;
      }

      if (vt->num > 0) {
        int num = vt->num;
        vs->f = (int           *)malloc(num*sizeof(int));
        vs->n = (unsigned char *)malloc(num*sizeof(unsigned char));
        if (!vs->f || !vs->n) ErrorExit(ERROR_NOMEMORY, "mrisRecordVertexState: could not allocate %dth array of %d elts", i, vsize);
        memcpy(vs->f,vt->f,num*sizeof(int));
        memcpy(vs->n,vt->n,num*sizeof(unsigned char));
        vs->num = num;
      }

      // Save the original position
      //
      vs->origx = v->origx; 
      vs->origy = v->origy;
      vs->origz = v->origz;

      vs->nx = v->nx;
      vs->ny = v->ny;
      vs->nz = v->nz;
    }
        
    // Save the checksum to test all put back correctly.
    //
    // For performance reasons, can't afford to check the hash every time, so vs->hash.hash is usually left 0
    //
    static size_t count, limit = 1; // no locking needed - doesn't matter if it gets it wrong...
    if (count++ > limit) {
      limit *= 2;
      if (limit > 1024*1024) limit = 1;
      mrisVertexHash(&vs->hash, mris, vno);
    }
  }

  return (dvs);
}


static void mrisRestoreOneVertexFaceState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs, int const i) {
  VERTEX_STATE const * const vs = &dvs->vs[i];
  int const vno = vs->vno;
  if (vno < 0) return;

  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];

  if (vs->vtotal == -1) return;

  int const num = vs->num;
  vt->num = num;
    
  if (num > 0) {
    vt->n = (unsigned char *)realloc(vt->n, num  *sizeof(unsigned char));
    vt->f = (int *)          realloc(vt->f, num  *sizeof(int));
      // use realloc for efficiency
        
    if (!vt->f || !vt->n)
      ErrorExit(ERROR_NOMEMORY, "mrisRestoreVertexState: could not reallocate v.f v.n arrays for num:%d", vt->num);
  }
  
  memcpy(vt->n,vs->n,num  *sizeof(unsigned char));
  memcpy(vt->f,vs->f,num  *sizeof(int));
}


static void mrisRestoreOneVertexState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs, int const i) {
  VERTEX_STATE const * const vs = &dvs->vs[i];
  int vsize;
  int const vno = vs->vno; 
  if (vno < 0) return;

  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
  VERTEX          * const v  = &mris->vertices         [vno];

  if (vs->vtotal == -1) goto Done;     // was ripped

  // Do the face
  //
  mrisRestoreOneVertexFaceState(mris, dvs, i);
    
  // Restore original position
  //
  v->nx = vs->nx;
  v->ny = vs->ny;
  v->nz = vs->nz;

  MRISsetOriginalXYZ(mris, vno, vs->origx, vs->origy, vs->origz);

  // Restore the topology
  //
  vt->nsizeMax = vs->nsizeMax;
  modVnum(mris, vno, vs->vnum, true);
  vt->v2num    = vs->v2num;
  vt->v3num    = vs->v3num;
  MRIS_setNsizeCur(mris, vno, vs->nsizeCur);
  cheapAssert(vt->vtotal == vs->vtotal);
  
  vsize = mrisVertexVSize(mris, vno);

  if (vsize > 0) {  // keep the zero-size ones for faster future realloc's, realloc the non-zero ones now
  
    vt->v = (int*)realloc(vt->v, vsize*sizeof(int));
      // use realloc for efficiency
        
    if (!vt->v)
      ErrorExit(ERROR_NOMEMORY, "mrisRestoreVertexState: could not reallocate topology arrays for num:%d vsize:%d", vt->num, vsize);
    
    memcpy(vt->v,vs->v,vsize*sizeof(int));
  }
  
Done:
  if (vs->hash.hash) {
    MRIS_HASH hash;
    mrisVertexHash(&hash, mris, vno);
    cheapAssert(hash.hash == vs->hash.hash);
  }
}

static void mrisRestoreFaceVertexState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs)
{
  MRIStruncateNFaces(mris, dvs->nfaces);    // remove the added faces
  int i;
  for (i = 0; i < dvs->nvertices; i++) {
    mrisRestoreOneVertexFaceState(mris, dvs, i);
  }
}

static void mrisRestoreVertexState(MRI_SURFACE *mris, DEFECT_VERTEX_STATE *dvs)
{
  MRIStruncateNFaces(mris, dvs->nfaces);  // remove the added faces
  int i;
  for (i = 0; i < dvs->nvertices; i++) {
    mrisRestoreOneVertexState(mris, dvs, i);
  }
  mrisCheckVertexFaceTopology(mris);
}


//
//


//==================================================================
//
#define VOLUME_SCALE 2.0f

#define MAX_SEGMENTS 10
#define MAX_SEGMENT_EDGES 100

#define USE_SOME_VERTICES 0  // FLO
#define USE_ALL_VERTICES 1

#define USE_NO_OVERLAP 0
#define USE_OVERLAP 1
#define USE_OVERLAP_IN_MUTATION 2

#define MAX_PATCHES 1000
#define SELECTION_PCT 0.50
#define ELITISM_PCT 0.10
#define MUTATION_PCT 0.1
static double MUTATION_PCT_INIT = (MUTATION_PCT * 1.0);
#define REPLACEMENT_PCT         \
  0.1 /* replace this many with \
mutated versions of best */
#define MAX_UNCHANGED 3
static int max_unchanged = MAX_UNCHANGED;

#define AREA_THRESHOLD 35.0f


// #define mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM

#define MARK_AMBIGUOUS_RETAIN   1
#define MARK_AMBIGUOUS_DISCARD  2
#define MARK_SEGMENTED          3

#define MARK_RETAIN     MARK_AMBIGUOUS_RETAIN
#define MARK_AMBIGUOUS  MARK_AMBIGUOUS_RETAIN
#define MARK_DISCARD    MARK_AMBIGUOUS_DISCARD



typedef struct DefectFacesCache {
    int  size;
    int  capacity;
    int* fnos;
} DefectFacesCache;

static void initDefectFacesCache(DefectFacesCache* p) { bzero(p, sizeof(*p)); }

static void finiDefectFacesCache(DefectFacesCache* p) { freeAndNULL(p->fnos); }

static void insertIntoDefectFacesCache(DefectFacesCache* p, int fno) {
    if (p->size == p->capacity) { 
        if (!(p->capacity *= 2)) p->capacity = 1000;
        p->fnos = (int*)realloc(p->fnos, p->capacity * sizeof(int));
    }
    p->fnos[p->size++] = fno;
}



static int mrisDumpDefectiveEdge(MRI_SURFACE *mris, int vno1, int vno2);
static SEGMENTATION *SEGMENTATIONalloc(int max_segments, int max_edges);
static void SEGMENTATIONfree(SEGMENTATION **segmentation);
static int findSegment(SEGMENTATION *segmentation);
static int compatibility(SEGMENTATION *segmentation, int segment_n, int edge_n);
static int addEdgeToSegment(SEGMENTATION *segmentation, int segment_n, int edge_n);
static void segmentEdge(MRIS *mris, SEGMENTATION *segmentation, int edge_n);
static int mergeSegmentToSegment(SEGMENTATION *segmentation, int segment_1, int segment_2);
static SEGMENTATION *segmentIntersectingEdges(MRIS *mris, DEFECT *defect, int *vertex_trans, ES *es, int *nedges);
static void saveSegmentation(
    MRIS *mris, MRIS *mris_corrected, DEFECT *defect, int *vertex_trans, ES *es, int nes, char *fname);
// generate an ordering based on the segmented overlapping edges
static void generateOrdering(DP *dp, SEGMENTATION *segmentation, int i);
static void savePatch(MRI *mri, MRIS *mris, MRIS *mris_corrected, DVS *dvs, DP *dp, char *fname, TOPOLOGY_PARMS *parms);
// compute statistics of the surface
static void mrisComputeSurfaceStatistics(
    MRIS *mris, MRI *mri, HISTOGRAM *h_k1, HISTOGRAM *h_k2, MRI *mri_k1_k2, MRI *mri_gray_white, HISTOGRAM *h_dot);
static void computeDefectStatistics(MRI *mri,
                                    MRIS *mris,
                                    DEFECT *defect,
                                    HISTOGRAM *h_white,
                                    HISTOGRAM *h_gray,
                                    MRI *mri_gw,
                                    HISTOGRAM *h_k1,
                                    HISTOGRAM *h_k2,
                                    MRI *mri_k1_k2,
                                    int verbose);
// initialization, free for the tessellated patch structure
static void TPinit(TP *tp);
static void TPfree(TP *tp);
static void TPprint(TP *tp);
static int mrisAddAllDefectFaces(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vertex_trans);
static int mrisFindDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect);
static int mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vtrans);
static int mrisMarkDefect(MRI_SURFACE *mris, DEFECT *defect, int mark);
static int mrisMarkDefectBorder(MRI_SURFACE *mris, DEFECT *defect, int mark);
static int mrisMarkDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect, int mark);

FACE_DEFECT_LIST *MRISmarkAmbiguousVertices(MRI_SURFACE *mris, int mark);
DEFECT_LIST *MRISsegmentDefects(MRI_SURFACE *mris, int mark_ambiguous, int mark_segmented);
static int mrisSegmentDefect(MRI_SURFACE *mris, int vno, DEFECT *defect, int mark_ambiguous, int mark_segmented);
//#if FIND_ENCLOSING_LOOP
static int mrisSimplyConnectedDefect(MRI_SURFACE *mris, DEFECT *defect, int mark_ambiguous, int mark_segmented);
static int mrisSegmentConnectedComponents(MRIS *mris);
//#endif
static int mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris,
                                        DEFECT *defect,
                                        FACE_DEFECT_LIST *fdl,
                                        float area_threshold,
                                        int mark_retain,
                                        int mark_discard,
                                        MHT *mht,
                                        int mode);
static int mrisTessellateDefect(MRI_SURFACE *mris,
                                MRI_SURFACE *mris_corrected,
                                DEFECT *defect,
                                int *vertex_trans,
                                MRI *mri,
                                HISTOGRAM *h_k1,
                                HISTOGRAM *h_k2,
                                MRI *mri_k1_k2,
                                HISTOGRAM *h_white,
                                HISTOGRAM *h_gray,
                                HISTOGRAM *h_border,
                                HISTOGRAM *h_grad,
                                MRI *mri_gray_white,
                                HISTOGRAM *h_dot,
                                TOPOLOGY_PARMS *parms);
static int mrisDefectRemoveDegenerateVertices(MRI_SURFACE *mris, float min_sphere_dist, DEFECT *defect);
static int mrisDefectRemoveProximalVertices(MRI_SURFACE *mris, float min_orig_dist, DEFECT *defect);
static int mrisDefectRemoveNegativeVertices(MRI_SURFACE *mris, DEFECT *defect);

static MRI *mriDefectVolume(MRIS *mris, EDGE_TABLE *etable, TOPOLOGY_PARMS *parms);
static void defectVolumeLikelihood(MRI *mri,
                                   MRI *mri_defect,
                                   MRI *mri_white,
                                   MRI *mri_gray,
                                   HISTOGRAM *h_white,
                                   HISTOGRAM *h_gray,
                                   float white_mean,
                                   float gray_mean,
                                   int type,
                                   int contrast);
static int vertexNeighbor(MRI_SURFACE *mris, int vno1, int vno2);
static int isVertexInsideFace(MRIS *mris, int vno, int fno);
static int isDiscarded(MRIS *mris, int vno1, int vno2);
static void removeVertex(MRIS *mris, int vno);
static int updateVertexTriangle(MRIS *mris, int vno, int fno);
// static void updateTriangle(MRIS *mris,int fno);
static void possiblyAddNewFaces(MRIS *mris, int vno1, int vno2);
static int possiblyAddEdgesAndFaces(MRIS *mris, int vno1, int vno2, int mode);
static int retessellateDefect(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, DVS *dvs, DP *dp);
static int mrisRetessellateDefect(MRI_SURFACE *mris,
                                  MRI_SURFACE *mris_corrected,
                                  DEFECT *defect,
                                  int *vertex_trans,
                                  EDGE *et,
                                  int nedges,
                                  int *ordering,
                                  EDGE_TABLE *etable);
static void defectMatch(MRI *mri, MRI_SURFACE *mris, DP *dp, int smooth, int match);
static void defectSmooth(MRI_SURFACE *mris, DP *dp, int niter, double alpha, int type);
static void defectMaximizeLikelihood(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha);
static void detectDefectFaces(MRIS *mris, DEFECT_PATCH *dp);
static int computePatchEulerNumber(MRIS *mris, DP *dp);
static void orientDefectFaces(MRIS *mris, DP *dp);

// segmentation of the edges into clusters
static void computeDefectFaceNormals(MRIS *mris, DP *dp, DefectFacesCache* dfc);
static void computeDefectVertexNormals(MRIS *mris, DP *dp);
static void printDefectStatistics(DP *dp);
static void computeDisplacement(MRI_SURFACE *mris, DP *dp);
static void updateVertexStatistics(
    MRIS *mris, MRIS *mris_corrected, DVS *dvs, RP *rp, DP *dp, int *vertex_trans, float fitness);
static int deleteWorstVertices(MRIS *mris, RP *rp, DEFECT *defect, int *vertex_trans, float fraction, int count);
static double mrisDefectPatchFitness(
    ComputeDefectContext* computeDefectContext,
    MRI_SURFACE *mris,
    MRI_SURFACE *mris_corrected,
    MRI *mri,
    DEFECT_PATCH *dp,
    int *vertex_trans,
    DEFECT_VERTEX_STATE *dvs,
    RP *rp,
    HISTOGRAM *h_k1,
    HISTOGRAM *h_k2,
    MRI *mri_k1_k2,
    HISTOGRAM *h_white,
    HISTOGRAM *h_gray,
    HISTOGRAM *h_border,
    HISTOGRAM *h_grad,
    MRI *mri_gray_white,
    HISTOGRAM *h_dot,
    TOPOLOGY_PARMS *parms);
static void vertexPseudoNormal(MRIS *mris1, int vn1, MRIS *mris2, int vn2, float norm[3]);

static int mrisMutateDefectPatch(DEFECT_PATCH *dp, EDGE_TABLE *etable, double pmutation);
static int mrisCrossoverDefectPatches(DEFECT_PATCH *dp1, DEFECT_PATCH *dp2, DEFECT_PATCH *dp_dst, EDGE_TABLE *etable);
static int defectPatchRank(DEFECT_PATCH *dps, int index, int npatches);
static int mrisCopyDefectPatch(DEFECT_PATCH *dp_src, DEFECT_PATCH *dp_dst);
static int mrisComputeOptimalRetessellation(MRI_SURFACE *mris,
                                            MRI_SURFACE *mris_corrected,
                                            MRI *mri,
                                            DEFECT *defect,
                                            int *vertex_trans,
                                            EDGE *et,
                                            int nedges,
                                            ES *es,
                                            int nes,
                                            HISTOGRAM *h_k1,
                                            HISTOGRAM *h_k2,
                                            MRI *mri_k1_k2,
                                            HISTOGRAM *h_white,
                                            HISTOGRAM *h_gray,
                                            HISTOGRAM *h_border,
                                            HISTOGRAM *h_grad,
                                            MRI *mri_gray_white,
                                            HISTOGRAM *h_dot,
                                            TOPOLOGY_PARMS *parms);
static int mrisComputeRandomRetessellation(MRI_SURFACE *mris,
                                           MRI_SURFACE *mris_corrected,
                                           MRI *mri,
                                           DEFECT *defect,
                                           int *vertex_trans,
                                           EDGE *et,
                                           int nedges,
                                           ES *es,
                                           int nes,
                                           HISTOGRAM *h_k1,
                                           HISTOGRAM *h_k2,
                                           MRI *mri_k1_k2,
                                           HISTOGRAM *h_white,
                                           HISTOGRAM *h_gray,
                                           HISTOGRAM *h_border,
                                           HISTOGRAM *h_grad,
                                           MRI *mri_gray_white,
                                           HISTOGRAM *h_dot,
                                           TOPOLOGY_PARMS *parms);
static OPTIMAL_DEFECT_MAPPING *mrisFindOptimalDefectMapping(MRIS *mris, DEFECT *defect);

static int mrisComputeGrayWhiteBorderDistributions(MRI_SURFACE *mris,
                                                   MRI *mri,
                                                   DEFECT *defect,
                                                   HISTOGRAM *h_white,
                                                   HISTOGRAM *h_gray,
                                                   HISTOGRAM *h_border,
                                                   HISTOGRAM *h_grad);


static int isVertexInsideFace(MRIS *mris, int vno, int fno)
{
  int i;
  VERTEX *V0, *V1, *V2, *V, *Vn;
  double sign, test, normal[3], orgn[3], tangent[3];

  /* the vertices constituting the triangle fno */
  V0 = &mris->vertices[mris->faces[fno].v[0]];
  V1 = &mris->vertices[mris->faces[fno].v[1]];
  V2 = &mris->vertices[mris->faces[fno].v[2]];

  /* if one of the vertices is outside of the defect, then vno is outside of fno */
  if (V0->marked == OUTSIDE_VERTEX || V1->marked == OUTSIDE_VERTEX || V2->marked == OUTSIDE_VERTEX) {
    return 0;
  }

  /* the vertex to be tested */
  Vn = &mris->vertices[vno];

  for (i = 0; i < 3; i++) {
    if (i) {
      // circular rotation
      V = V0;
      V0 = V1;
      V1 = V2;
      V2 = V;
    }

    // compute normal for edge V1<-->V2
    orgn[0] = (V1->cx + V2->cx) / 2.0;
    orgn[1] = (V1->cy + V2->cy) / 2.0;
    orgn[2] = (V1->cz + V2->cz) / 2.0;

    tangent[0] = V2->cx - V1->cx;
    tangent[1] = V2->cy - V1->cy;
    tangent[2] = V2->cz - V1->cz;

    // normal to edge in the planar basis
    F_CROSS(orgn, tangent, normal);

    tangent[0] = V0->cx - orgn[0];
    tangent[1] = V0->cy - orgn[1];
    tangent[2] = V0->cz - orgn[2];

    sign = F_DOT(tangent, normal);

    tangent[0] = Vn->cx - orgn[0];
    tangent[1] = Vn->cy - orgn[1];
    tangent[2] = Vn->cz - orgn[2];

    test = F_DOT(tangent, normal);

    if (sign * test < 0) {
      return 0; /* outside */
    }
  }
  return 1;
}

static void vertexPseudoNormal(MRIS *mris1, int vn1, MRIS *mris2, int vn2, float norm[3])
{
  int n, n0, n1, n2;
  float v1[3], v2[3], alpha;

  // fprintf(stderr,"-vpn: %d and %d - ",vn1,vn2);

  norm[0] = norm[1] = norm[2] = 0;

  {
    VERTEX_TOPOLOGY const * const vt = &mris1->vertices_topology[vn1];
    VERTEX          const * const v  = &mris1->vertices         [vn1];

    for (n = 0; n < vt->num; n++) {
      int const fno = vt->f[n];
      FACE * const face = &mris1->faces[fno];
      if (face->marked) {
        continue;  // defect
      }

      n0 = vt->n[n];
      n1 = (n0 == 2) ? 0 : n0 + 1;
      n2 = (n0 == 0) ? 2 : n0 - 1;

  #if 1
      if ((face->v[n0] != vn1) || (face->v[n1] == vn1) || (face->v[n2] == vn1) || (face->v[n2] == face->v[n1])) {
        if (1) {
          // verbose>=VERBOSE_MODE_MEDIUM){
          if (face->v[n0] != vn1) {
            fprintf(WHICH_OUTPUT, "error for vno in face %d", vt->f[n]);
          }
          if (face->v[n1] == vn1) {
            fprintf(WHICH_OUTPUT, "error for vn1 in face %d", vt->f[n]);
          }
          if (face->v[n2] == vn1) {
            fprintf(WHICH_OUTPUT, "error for vn2 in face %d", vt->f[n]);
          }
          if (face->v[n2] == face->v[n1]) {
            fprintf(WHICH_OUTPUT, "error for vn in face %d", vt->f[n]);
          }

          fprintf(WHICH_OUTPUT, "face %d (%d,%d,%d) != (%d)\n", vt->f[n], face->v[n0], face->v[n1], face->v[n2], vn1);

          if (1)  // verbose==VERBOSE_MODE_MEDIUM)
          {
            fprintf(stderr, "vertexPseudoNormal: SHOULD NOT HAPPEN\n");
          }

          //                      MRISwrite(mris,"rh.testdebug1");
          // MRISrestoreVertexPositions(mris,CANONICAL_VERTICES);
          // MRISwrite(mris,"rh.testdebug2");
          if (1)  // verbose==VERBOSE_MODE_HIGH)
          {
            ErrorExit(ERROR_BADPARM, "vertexPseudoNormal: SHOULD NOT HAPPEN\n");
          }
        }
      }
  #endif
      v1[0] = mris1->vertices[face->v[n1]].origx - v->origx;
      v1[1] = mris1->vertices[face->v[n1]].origy - v->origy;
      v1[2] = mris1->vertices[face->v[n1]].origz - v->origz;

      v2[0] = mris1->vertices[face->v[n2]].origx - v->origx;
      v2[1] = mris1->vertices[face->v[n2]].origy - v->origy;
      v2[2] = mris1->vertices[face->v[n2]].origz - v->origz;

      alpha = MAX(0.0, MIN(1.0, F_DOT(v1, v2) / NORM3(v1) / NORM3(v2)));
      alpha = acos(alpha);

      FaceNormCacheEntry const * fNorm = getFaceNorm(mris1, fno);

      norm[0] += alpha * fNorm->nx;
      norm[1] += alpha * fNorm->ny;
      norm[2] += alpha * fNorm->nz;
    }
  }
  
  {
    VERTEX_TOPOLOGY const * const vt = &mris2->vertices_topology[vn2];
    VERTEX          const * const v  = &mris2->vertices         [vn2];
    for (n = 0; n < vt->num; n++) {
      int const fno = vt->f[n];
      FACE* const face = &mris2->faces[fno];
      n0 = vt->n[n];
      n1 = (n0 == 2) ? 0 : n0 + 1;
      n2 = (n0 == 0) ? 2 : n0 - 1;
#if 1
      if ((face->v[n0] != vn2) || (face->v[n1] == vn2) || (face->v[n2] == vn2) || (face->v[n2] == face->v[n1])) {
        if (1) {
          // verbose>=VERBOSE_MODE_MEDIUM){
          if (face->v[n0] != vn2) {
            fprintf(WHICH_OUTPUT, "error for vno in face %d", vt->f[n]);
          }
          if (face->v[n1] == vn2) {
            fprintf(WHICH_OUTPUT, "error for vn1 in face %d", vt->f[n]);
          }
          if (face->v[n2] == vn2) {
            fprintf(WHICH_OUTPUT, "error for vn2 in face %d", vt->f[n]);
          }
          if (face->v[n2] == face->v[n1]) {
            fprintf(WHICH_OUTPUT, "error for vn in face %d", vt->f[n]);
          }

          fprintf(WHICH_OUTPUT, "face %d (%d,%d,%d) != (%d)\n", vt->f[n], face->v[n0], face->v[n1], face->v[n2], vn2);

          if (1)  // verbose==VERBOSE_MODE_MEDIUM)
          {
            fprintf(stderr, "vertexPseudoNormal: SHOULD NOT HAPPEN\n");
          }

          //                      MRISwrite(mris,"rh.testdebug1");
          // MRISrestoreVertexPositions(mris,CANONICAL_VERTICES);
          // MRISwrite(mris,"rh.testdebug2");
          if (1)  // verbose==VERBOSE_MODE_HIGH)
          {
            ErrorExit(ERROR_BADPARM, "vertexPseudoNormal: SHOULD NOT HAPPEN\n");
          }
        }
      }
  #endif
      v1[0] = mris2->vertices[face->v[n1]].origx - v->origx;
      v1[1] = mris2->vertices[face->v[n1]].origy - v->origy;
      v1[2] = mris2->vertices[face->v[n1]].origz - v->origz;

      v2[0] = mris2->vertices[face->v[n2]].origx - v->origx;
      v2[1] = mris2->vertices[face->v[n2]].origy - v->origy;
      v2[2] = mris2->vertices[face->v[n2]].origz - v->origz;

      alpha = MAX(0.0, MIN(1.0, F_DOT(v1, v2) / NORM3(v1) / NORM3(v2)));
      alpha = acos(alpha);

      FaceNormCacheEntry const * fNorm = getFaceNorm(mris2, fno);
      norm[0] += alpha * fNorm->nx;
      norm[1] += alpha * fNorm->ny;
      norm[2] += alpha * fNorm->nz;
    }
  }
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static void computeDefectTangentPlaneAtVertex(MRIS *mris, int vno)
{
  float a[3], b[3], c[3], len;

  VERTEX *v;

  v = &mris->vertices[vno];

  a[0] = v->nx;
  a[1] = v->ny;
  a[2] = v->nz;

  /* pick one random vector */
  b[0] = 1.0f;
  b[1] = 0.0f;
  b[2] = 0.0f;

  F_CROSS(a, b, c);
  len = SQR(c[0]) + SQR(c[1]) + SQR(c[2]);

  if (FZERO(len)) {
    /* the vector b was parallel to a */
    b[0] = 0.0f;
    b[1] = 1.0f;
    b[2] = 0.0f;
    F_CROSS(a, b, c);
    len = SQR(c[0]) + SQR(c[1]) + SQR(c[2]);
  }
  /* normalize */
  len = sqrt(len);

  if (FZERO(len)) {
    fprintf(WHICH_OUTPUT, "first tangent vector of length zero at vertex %d\n", vno);
    len = 1;
  }
  v->e1x = c[0] / len;
  v->e1y = c[1] / len;
  v->e1z = c[2] / len;

  F_CROSS(a, c, b);

  /* normalize */
  len = sqrt(SQR(b[0]) + SQR(b[1]) + SQR(b[2]));
  if (FZERO(len)) {
    fprintf(WHICH_OUTPUT, "second tangent vector of length zero at vertex %d\n", vno);
    len = 1;
  }
  v->e2x = b[0] / len;
  v->e2y = b[1] / len;
  v->e2z = b[2] / len;
}


static void computeDefectSecondFundamentalForm(MRIS *mris, TP *tp)
{
  int p, vno, i, n, nbad = 0;
  MATRIX *m_U, *m_Ut, *m_tmp1, *m_tmp2, *m_inverse, *m_eigen, *m_Q;
  VECTOR *v_c, *v_z, *v_n, *v_e1, *v_e2, *v_yi;
  float k1, k2, evalues[3], a11, a12, a21, a22, cond_no, kmax, kmin, rsq, k;
  double ui, vi;

  /* allocation of diverse vectors matrices*/
  v_c = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL);
  v_e1 = VectorAlloc(3, MATRIX_REAL);
  v_e2 = VectorAlloc(3, MATRIX_REAL);
  v_yi = VectorAlloc(3, MATRIX_REAL);
  m_Q = MatrixAlloc(2, 2, MATRIX_REAL); /* the quadratic form */
  m_eigen = MatrixAlloc(2, 2, MATRIX_REAL);

  /* compute faces only for modified vertices */
  for (p = 0; p < tp->nvertices; p++) {
    vno = tp->vertices[p];

    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];

    /* compute tangent plane */
    computeDefectTangentPlaneAtVertex(mris, vno);

// FLO TO BE CHECKED !!!
#if 0
    mrisFindSecondNeighborhood(mris, vno, nbrs, &num_nbrs) ;
    if (num_nbrs < 3)
    {
      continue ;
    }
    nvertices++;

    /* need ex,ey,ez */
    MRIScomputeSecondFundamentalFormAtVertex(mris, vno, nbrs, num_nbrs) ;
#endif

    VECTOR_LOAD(v_n,  vertex->nx,  vertex->ny,  vertex->nz);
    VECTOR_LOAD(v_e1, vertex->e1x, vertex->e1y, vertex->e1z);
    VECTOR_LOAD(v_e2, vertex->e2x, vertex->e2y, vertex->e2z);

    if (vertext->vtotal <= 0) {
      continue;
    }

    m_U = MatrixAlloc(vertext->vtotal, 3, MATRIX_REAL);
    v_z = VectorAlloc(vertext->vtotal, MATRIX_REAL);

    /* fit a quadratic form to the surface at this vertex */
    kmin = 10000.0f;
    kmax = -kmin;
    for (n = i = 0; i < vertext->vtotal; i++) {
      VERTEX const * const vnb = &mris->vertices[vertext->v[i]];

      /* calculate the projection of this vertex onto the local tangent plane */
      VECTOR_LOAD(v_yi, vnb->origx - vertex->origx, vnb->origy - vertex->origy, vnb->origz - vertex->origz);
      ui = V3_DOT(v_yi, v_e1);
      vi = V3_DOT(v_yi, v_e2);
      *MATRIX_RELT(m_U, n + 1, 1) = ui * ui;
      *MATRIX_RELT(m_U, n + 1, 2) = 2 * ui * vi;
      *MATRIX_RELT(m_U, n + 1, 3) = vi * vi;
      VECTOR_ELT(v_z, n + 1) = V3_DOT(v_n, v_yi); /* height above TpS */
      rsq = ui * ui + vi * vi;

      if (!FZERO(rsq)) {
        k = VECTOR_ELT(v_z, n + 1) / rsq;
        if (k > kmax) {
          kmax = k;
        }
        if (k < kmin) {
          kmin = k;
        }
      }
      n++;
    }
    
    MatrixBuffer m_Ut_buffer, m_tmp2_buffer;
    m_Ut   = MatrixAlloc2(m_U->cols, m_U->rows, MATRIX_REAL, &m_Ut_buffer); 
    m_Ut   = MatrixTranspose(m_U, m_Ut);        /* Ut */
    
    m_tmp2 = MatrixAlloc2(m_Ut->rows, m_U->cols, MATRIX_REAL, &m_tmp2_buffer);
    m_tmp2 = MatrixMultiply(m_Ut, m_U, m_tmp2); /* Ut U */
    cond_no = MatrixConditionNumber(m_tmp2);

    MatrixBuffer m_inverse_buffer;
    m_inverse = MatrixAlloc2(m_U->cols, m_U->cols, MATRIX_REAL, &m_inverse_buffer);
    m_inverse = MatrixSVDInverse(m_tmp2, m_inverse); /* (Ut U)^-1 */

    if (!m_inverse) {
      /* singular matrix - must be planar?? */
      nbad++;
      evalues[0] = evalues[1] = 0.0;
    }
    else {
      MatrixBuffer m_tmp1_buffer;
      m_tmp1 = MatrixAlloc2(m_Ut->rows, v_z->cols, MATRIX_REAL, &m_tmp1_buffer); 
      m_tmp1 = MatrixMultiply(m_Ut, v_z, m_tmp1); /* Ut z */
      MatrixMultiply(m_inverse, m_tmp1, v_c);   /* (Ut U)^-1 Ut z */

      /* now build Hessian matrix */
      *MATRIX_RELT(m_Q, 1, 1) = 2 * VECTOR_ELT(v_c, 1);
      *MATRIX_RELT(m_Q, 1, 2) = *MATRIX_RELT(m_Q, 2, 1) = 2 * VECTOR_ELT(v_c, 2);
      *MATRIX_RELT(m_Q, 2, 2) = 2 * VECTOR_ELT(v_c, 3);

      if (cond_no >= ILL_CONDITIONED) {
        vertex->k1 = k1 = kmax;
        vertex->k2 = k2 = kmin;

        vertex->K = k1 * k2;
        vertex->H = (k1 + k2) / 2;
        MatrixFree(&m_Ut);
        MatrixFree(&m_tmp2);
        MatrixFree(&m_U);
        VectorFree(&v_z);
        MatrixFree(&m_tmp1);
        MatrixFree(&m_inverse);
        continue;
      }

      /* the columns of m_eigen will be the eigenvectors of m_Q */
      if (MatrixEigenSystem(m_Q, evalues, m_eigen) == NULL) {
        nbad++;
        MatrixSVDEigenValues(m_Q, evalues);
        vertex->k1 = k1 = evalues[0];
        vertex->k2 = k2 = evalues[1];
        vertex->K = k1 * k2;
        vertex->H = (k1 + k2) / 2;
        MatrixFree(&m_Ut);
        MatrixFree(&m_tmp2);
        MatrixFree(&m_U);
        VectorFree(&v_z);
        MatrixFree(&m_tmp1);
        MatrixFree(&m_inverse);
        continue;
      }

      MatrixFree(&m_tmp1);
      MatrixFree(&m_inverse);
    }
    k1 = evalues[0];
    k2 = evalues[1];
    vertex->k1 = k1;
    vertex->k2 = k2;
    vertex->K = k1 * k2;
    vertex->H = (k1 + k2) / 2;

    /* now update the basis vectors to be the principal directions */
    a11 = *MATRIX_RELT(m_eigen, 1, 1);
    a12 = *MATRIX_RELT(m_eigen, 1, 2);
    a21 = *MATRIX_RELT(m_eigen, 2, 1);
    a22 = *MATRIX_RELT(m_eigen, 2, 2);
    vertex->e1x = V3_X(v_e1) * a11 + V3_X(v_e2) * a21;
    vertex->e1y = V3_Y(v_e1) * a11 + V3_Y(v_e2) * a21;
    vertex->e1z = V3_Z(v_e1) * a11 + V3_Z(v_e2) * a21;
    vertex->e2x = V3_X(v_e1) * a12 + V3_X(v_e2) * a22;
    vertex->e2y = V3_Y(v_e1) * a12 + V3_Y(v_e2) * a22;
    vertex->e2z = V3_Z(v_e1) * a12 + V3_Z(v_e2) * a22;

    MatrixFree(&m_Ut);
    MatrixFree(&m_tmp2);
    MatrixFree(&m_U);
    VectorFree(&v_z);
  }

  MatrixFree(&m_eigen);
  VectorFree(&v_e1);
  VectorFree(&v_e2);
  VectorFree(&v_c);
  VectorFree(&v_n);
  VectorFree(&v_yi);
  MatrixFree(&m_Q);
}


double computeSurfaceLikelihood(MRIS *mris, DP *dp, int verbose)
{
  int n;
  double Km, Hm;
  VERTEX *v;

  computeDefectSecondFundamentalForm(mris, &dp->tp);
  Km = Hm = 0.0;
  for (n = 0; n < dp->tp.ninside; n++) {
    v = &mris->vertices[dp->tp.vertices[n]];
    Km += v->K * v->K;
    Hm += v->H * v->H;
  }
  Km /= dp->tp.ninside;
  Hm /= dp->tp.ninside;

  if (verbose) {
    fprintf(WHICH_OUTPUT, "CURV: K=%f H=%f\n", Km, Hm);
  }

  return -Hm;
}

double computeVolumeLikelihood(DP *dp, int contrast, int verbose)
{
  double int_g, int_w, val, val2, dval, delta, delta_ll;
  int i, j, k;
  double white_mean, gray_mean, white_vol, gray_vol;
  double scale, valvox, delta_NRG;

  // double sigma,dwhite,dgray,Dw,Dg;

  MRI *mri_defect, *mri_sign, *mri_white, *mri_gray, *mri_init;
  TP *tp;

  tp = &dp->tp;
  mri_defect = dp->mri_defect;
  mri_sign = dp->mri_defect_sign;
  mri_white = dp->mri_defect_white;
  mri_gray = dp->mri_defect_gray;
  mri_init = dp->mri_defect_initial_sign;

  white_mean = dp->defect->white_mean;
  gray_mean = dp->defect->gray_mean;
  // sigma = fabs(white_mean-gray_mean)/2.0;

  scale = mri_defect->xsize;
  scale = 1 / (scale * scale * scale);

  white_vol = gray_vol = 0.0001;
  int_g = gray_vol * dp->defect->gray_mean;
  int_w = white_vol * dp->defect->white_mean;
  // dwhite = dgray = 0.0;
  // Dw     = Dg    = 0.0;
  delta_NRG = 0.0;
  delta = delta_ll = 0.0;

  /* compute differences with mri_defect_initial_sign */
  for (k = 1; k < mri_sign->depth - 1; k++)
    for (j = 1; j < mri_sign->height - 1; j++)
      for (i = 1; i < mri_sign->width - 1; i++) {
        val = MRIFvox(mri_sign, i, j, k);   // new signed volume
        val2 = MRIFvox(mri_init, i, j, k);  // old signed volume

        if ((val >= 1.0f && val2 >= 1.0f) || (val <= -1.0f && val2 <= -1.0f)) {
          continue;  // compute the difference with the orig segmentation only
        }

        dval = fabs(val2 - val) / 2.0f;  // change in volume
        delta_NRG += (val2 - val) / 2.0f * (MRIFvox(mri_white, i, j, k) - MRIFvox(mri_gray, i, j, k));
        valvox = MRIvox(mri_defect, i, j, k);  // intensity at location (i,j,k)

        if (val2 > val) {
          white_vol += dval;
          int_w += dval * valvox;
        }
        else {
          gray_vol += dval;
          int_g += dval * valvox;
        }

        /*
            if(valvox < white_mean) Dw     += (1-val2)/2.0f   * SQR((valvox - white_mean)/sigma);
            if(valvox > gray_mean)  Dg     += (1.0+val2)/2.0f * SQR((valvox - gray_mean)/sigma);
            if(valvox < white_mean) dwhite += (1.0-val)/2.0f  * SQR((valvox - white_mean)/sigma);
            if(valvox > gray_mean)  dgray  += (1.0+val)/2.0f  * SQR((valvox - gray_mean)/sigma);
            delta_ll += MRIFvox(mri_white,i,j,k)-MRIFvox(mri_gray,i,j,k);
        */
      }
  int_w /= white_vol;
  int_g /= gray_vol;

  /*
   dwhite *= scale;
   dgray *= scale;
   Dw *= scale;
   Dg *= scale;
   delta=(dwhite+dgray-Dw-Dg)/(gray_vol+white_vol);
  */

  delta_NRG *= scale;
  delta_ll = delta_NRG / (white_vol + gray_vol);
  if (verbose) {
    if (gray_vol > white_vol) {
      fprintf(stderr, "CUTTING HANDLE\n");
    }
    else {
      fprintf(stderr, "FILLING HANDLE\n");
    }

    fprintf(stderr,
            "Gray [ Vol=%f - Int = %f ]- White [ Vol = %f - Int = %f ] %f (%f)\n",
            gray_vol,
            int_g,
            white_vol,
            int_w,
            delta_NRG,
            delta_ll);
  }

  // dp->tp.unmri_ll=(white_ll+gray_ll);
  dp->tp.unmri_ll = delta_ll;
  // return (white_ll+gray_ll);
  return (-delta_NRG);
}

double l_mri = 1.0;
double l_unmri = 1.0;
double l_curv = 1.0;
double l_qcurv = 1.0;

double l_vol = 1.0;
double l_surf = 1.0;
double l_wm = 1.0;

double MRIScomputeFitness(MRIS *mris, TOPOFIX_PARMS *parms, int verbose)
{
  static int first_time = 1;
  double fitness, unmri_ll, curv_ll;
  DP *dp;
  TP *tp;

#if 0
  static int now = 0;
  static int def = -1;
  static int when = 0;
#endif

  fitness = unmri_ll = curv_ll = 0.0;

  MRIScomputeNormals(mris);
  MRIScomputeTriangleProperties(mris);
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);

  dp = (DP *)parms->dp;
  tp = &dp->tp;

  dp->tp.face_ll = 0.0f;
  dp->tp.vertex_ll = 0.0f;
  dp->tp.curv_ll = 0.0f;
  dp->tp.qcurv_ll = 0.0f;
  dp->tp.unmri_ll = 0.0f;

  if (first_time) {
    l_vol = parms->l_unmri;
    l_surf = parms->l_curv;
    l_wm = 0;  // parms->l_wm;
    first_time = 0;
  }

  // compute the likelihood of the topologically-corrected volume
  if (!FZERO(l_vol) && dp->mri_defect->width >= 5 && dp->mri_defect->height >= 5 && dp->mri_defect->depth >= 5) {
    /* compute the signed distance volume */
    MRIScomputeDistanceVolume(parms, 2.0);
#if 0
    if (def == parms->defect_number && when == 1)
    {
      now = 1;
    }
    if (def != parms->defect_number)
    {
      def = parms->defect_number;
      when = 1;
    }
    if (now)
    {
      MRISsaveLocal(parms->mris_defect,parms,"./defect");
      now=0;
      when=0;
    }
#endif
    unmri_ll = computeVolumeLikelihood(dp, parms->contrast, 0);
    fitness += l_vol * unmri_ll;
  }
  // compute the likelihood of the surface
  if (!FZERO(l_surf)) {
    curv_ll = computeSurfaceLikelihood(mris, dp, 0);
    fitness += l_surf * curv_ll;
  }
  unmri_ll = -unmri_ll;
  curv_ll = -curv_ll;
  dp->tp.unmri_ll = unmri_ll;
  dp->tp.curv_ll = curv_ll;
  dp->tp.qcurv_ll = 2 * curv_ll;
  dp->tp.cll = curv_ll;
  dp->tp.qcll = 2 * curv_ll;

  return (fitness);
}


static int mrisCheckDefectFaces(MRI_SURFACE *mris, DEFECT_PATCH *dp)
{
  int fno, n1, n2, vno1, vno2, fshared;

  for (n1 = 0; n1 < dp->tp.nvertices; n1++) {
    vno1 = dp->tp.vertices[n1];
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno1];
    for (n2 = 0; n2 < v->vnum; n2++) {
      vno2 = v->v[n2];
      for (fshared = fno = 0; fno < v->num; fno++)
        if (vertexInFace(mris, vno2, v->f[fno])) fshared++;
      if (fshared != 2) return (-1);
    }
  }
  return (1);
}

float mrisDefectFaceMRILogLikelihood(
    MRI_SURFACE *mris, MRI *mri, TP *tp, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white)
{
  int n, vno0, vno1, vno2;
  double x, y, z, nx, ny, nz, xv, yv, zv, white_val, gray_val, val;
  double int_w, int_g;
  double fll, t_area, tf_area;
  VERTEX *v0, *v1, *v2;

#if 1

  t_area = tf_area = 0.0;

  fll = 0.0;
  int_w = int_g = 0.0;
  for (n = 0; n < tp->nfaces; n++) {
    int const fno = tp->faces[n];
    FACE * const face = &mris->faces[fno];
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    
    vno0 = face->v[0];
    vno1 = face->v[1];
    vno2 = face->v[2];

    v0 = &mris->vertices[vno0];
    v1 = &mris->vertices[vno1];
    v2 = &mris->vertices[vno2];

    /* find face centroid */
    x = (v0->origx + v1->origx + v2->origx) / 3.0f;
    y = (v0->origy + v1->origy + v2->origy) / 3.0f;
    z = (v0->origz + v1->origz + v2->origz) / 3.0f;

    /* face normal */
    nx = fNorm->nx;
    ny = fNorm->ny;
    nz = fNorm->nz;

#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &white_val);

    int_w += white_val;

#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &gray_val);

    int_g += gray_val;

    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);
    fll += log(val);

    t_area += face->area;
    tf_area += log(val) * face->area;
  }

  if (tp->nfaces) {
    tp->face_ll = (float)fll / (float)tp->nfaces;
    tp->fll = tf_area / t_area;
    int_w /= (double)tp->nfaces;
    int_g /= (double)tp->nfaces;
    // fs_topo_fixer_test
    //  fprintf(stderr,"face : gray = %f and white = %f\n",int_g,int_w);
  }
  else {
    tp->fll = 0.0;
    tp->face_ll = 0.0;
  }

  return tp->face_ll;

#else
  double x, y z, xa, ya, za, xc, yc, zc, t0, t1, adx, ady, adz, dx, dy, dz, grad, cdx, cdy, cdz, alen, clen, delta_t0,
      delta_t1, len, nx, ny, nz, xv, yv, zv, white_val, gray_val, cnx, cny, cnz, dot, val;
  double ll = 0.0, jll = 0.0;

  adx = v1->x - v0->x;
  ady = v1->y - v0->y;
  adz = v1->z - v0->z;
  alen = sqrt(SQR(adx) + SQR(ady) + SQR(adz));
  cdx = v2->x - v0->x;
  cdy = v2->y - v0->y;
  cdz = v2->z - v0->z;
  clen = sqrt(SQR(cdx) + SQR(cdy) + SQR(cdz));

  /*
  sample along legs of the triangle making sure the maximum spacing
  between samples (along the longer leg) is SAMPLE_DIST.
  */

  /*
  move along v0->v1 and v3->v2 lines and draw in crossing line to fill face
  t0 parameterizes lines from v0->v1 and v0->v2
  */
  if (FZERO(alen) && FZERO(clen)) {
    delta_t0 = 0.99;
  }
  else {
    delta_t0 = (alen > clen) ? (SAMPLE_DIST / alen) : (SAMPLE_DIST / clen);
  }
  if (FZERO(delta_t0))
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "mrisDefectFaceMRILogLikelihood: face has infinite leg (%d, %d)\n", alen, clen));

  if (delta_t0 >= 1.0) {
    delta_t0 = 0.99;
  }

  /* delta_t0 is % of alen or clen (whichever is bigger) of SAMPLE_DIST */
  for (nsamples = 0, ll = 0.0f, t0 = 0; t0 <= 1.0f; t0 += delta_t0) {
    /* compute points (xa,ya,za) and (xc,yc,zc) on the a and c lines resp. */
    xa = v0->x + t0 * adx;
    ya = v0->y + t0 * ady;
    za = v0->z + t0 * adz;
    xc = v0->x + t0 * cdx;
    yc = v0->y + t0 * cdy;
    zc = v0->z + t0 * cdz;
    dx = xc - xa;
    dy = yc - ya;
    dz = zc - za;
    len = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
    if (FZERO(len)) {
      delta_t1 = 0.99;
    }
    else {
      delta_t1 = SAMPLE_DIST / len; /* sample at SAMPLE_DIST intervals */
      if (delta_t1 >= 1.0f) {
        delta_t1 = 0.99;
      }
    }

    /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
    for (t1 = 0; t1 <= 1.0f; t1 += delta_t1, nsamples++) {
      /* compute a point on the line connecting a and c */
      x = xa + t1 * dx;
      y = ya + t1 * dy;
      z = za + t1 * dz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif

      MRIsampleVolume(mri, xv, yv, zv, &gray_val);
      ll += log(h_gray->counts[nint(gray_val)]);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &white_val);
      ll += log(h_white->counts[nint(white_val)]);
      grad = white_val - gray_val;
      bin = nint((grad - h_grad->bins[0]) / h_grad->bin_size);
      if (bin < 0) {
        bin = 0;
      }
      else if (bin >= h_grad->nbins) {
        bin = h_grad->nbins - 1;
      }
      ll += log(h_grad->counts[bin]);
      MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);
      jll += log(val);
    }
    /* compute last point on line */
    t1 = 1.0f;
    x = xa + t1 * dx;
    y = ya + t1 * dy;
    z = za + t1 * dz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &gray_val);
    ll += log(h_gray->counts[nint(gray_val)]);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &white_val);
    ll += log(h_white->counts[nint(white_val)]);
    grad = white_val - gray_val;
    bin = nint((grad - h_grad->bins[0]) / h_grad->bin_size);
    if (bin < 0) {
      bin = 0;
    }
    else if (bin >= h_grad->nbins) {
      bin = h_grad->nbins - 1;
    }
    ll += log(h_grad->counts[bin]);
    nsamples++;
    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);
    jll += log(val);
  }

  /* compute last line on the a and c lines resp. */
  t0 = 1.0f;
  xa = v0->x + t0 * adx;
  ya = v0->y + t0 * ady;
  za = v0->z + t0 * adz;
  xc = v0->x + t0 * cdx;
  yc = v0->y + t0 * cdy;
  zc = v0->z + t0 * cdz;
  dx = xc - xa;
  dy = yc - ya;
  dz = zc - za;
  len = sqrt(SQR(dx) + SQR(dy) + SQR(dz));
  if (FZERO(len)) {
    delta_t1 = 0.99;
  }
  else {
    delta_t1 = SAMPLE_DIST / len; /* sample at SAMPLE_DIST intervals */
    if (delta_t1 >= 1.0f) {
      delta_t1 = 0.99;
    }
  }

  /* now draw a line from (xa,ya,za) to (xc, yc, zc) */
  for (t1 = 0; t1 <= 1.0f; t1 += delta_t1, nsamples++) {
    /* compute a point on the line connecting a and c */
    x = xa + t1 * dx;
    y = ya + t1 * dy;
    z = za + t1 * dz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &gray_val);
    ll += log(h_gray->counts[nint(gray_val)]);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &white_val);
    ll += log(h_white->counts[nint(white_val)]);
    grad = white_val - gray_val;
    bin = nint((grad - h_grad->bins[0]) / h_grad->bin_size);
    if (bin < 0) {
      bin = 0;
    }
    else if (bin >= h_grad->nbins) {
      bin = h_grad->nbins - 1;
    }
    ll += log(h_grad->counts[bin]);
    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);
    jll += log(val);
  }
  /* compute last point on line */
  t1 = 1.0f;
  x = xa + t1 * dx;
  y = ya + t1 * dy;
  z = za + t1 * dz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
  mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
  MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
  MRIsampleVolume(mri, xv, yv, zv, &gray_val);
  ll += log(h_gray->counts[nint(gray_val)]);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
  mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
  MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
  MRIsampleVolume(mri, xv, yv, zv, &white_val);
  ll += log(h_white->counts[nint(white_val)]);
  grad = white_val - gray_val;
  bin = nint((grad - h_grad->bins[0]) / h_grad->bin_size);
  if (bin < 0) {
    bin = 0;
  }
  else if (bin >= h_grad->nbins) {
    bin = h_grad->nbins - 1;
  }
  ll += log(h_grad->counts[bin]);
  MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);
  jll += log(val);
  nsamples++;

#if 0
  return(ll/nsamples) ;
#else
  return (jll / nsamples);
#endif
#endif
}

float mrisDefectVertexMRILogLikelihood(
    MRI_SURFACE *mris, MRI *mri, TP *tp, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white)
{
  int n;
  double x, y, z, nx, ny, nz, xv, yv, zv, white_val, gray_val, val;
  double int_w, int_g;
  double v_ll, total_ll, t_area, tv_area;
  VERTEX *v;

  total_ll = 0.0;
  t_area = tv_area = 0.0;
  int_w = int_g = 0.0;

  for (n = 0; n < tp->nvertices; n++) {
    v = &mris->vertices[tp->vertices[n]];

    x = v->origx;
    y = v->origy;
    z = v->origz;

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &white_val);

    int_w += white_val;

#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &gray_val);

    int_g += gray_val;

    MRIsampleVolume(mri_gray_white, white_val, gray_val, 0, &val);

    v_ll = log(val);
    total_ll += v_ll;
    tv_area += v_ll * v->area;
    t_area += v->area;
  }

  if (tp->nvertices) {
    tp->vertex_ll = total_ll / (double)tp->nvertices;
    tp->vll = tv_area / t_area;
    int_w /= (double)tp->nvertices;
    int_g /= (double)tp->nvertices;
    // fs_topo_fixer_test
    //  fprintf(stderr,"vertex : gray = %f and white = %f\n",int_g,int_w);
  }
  else {
    tp->vll = 0.0;
    tp->vertex_ll = 0.0;
  }

  return tp->vertex_ll;
}

static double mrisComputeDefectMRILogLikelihood(
    MRI_SURFACE *mris, MRI *mri, TP *tp, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white)
{
  double mll;

  mll = mrisDefectVertexMRILogLikelihood(mris, mri, tp, h_white, h_gray, h_grad, mri_gray_white);
  mll += mrisDefectFaceMRILogLikelihood(mris, mri, tp, h_white, h_gray, h_grad, mri_gray_white);
  return mll;
}


double mrisComputeDefectLogLikelihood(
    ComputeDefectContext* computeDefectContext, 
    MRI_SURFACE *mris,
    MRI *mri,
    DEFECT_PATCH *dp,
    HISTOGRAM *h_k1,
    HISTOGRAM *h_k2,
    MRI *mri_k1_k2,
    HISTOGRAM *h_white,
    HISTOGRAM *h_gray,
    HISTOGRAM *h_border,
    HISTOGRAM *h_grad,
    MRI *mri_gray_white,
    HISTOGRAM *h_dot,
    TOPOLOGY_PARMS *parms)
{
  static int first_time = 1;
  double ll = 0.0;

  dp->tp.face_ll = 0.0f;
  dp->tp.vertex_ll = 0.0f;
  dp->tp.curv_ll = 0.0f;
  dp->tp.qcurv_ll = 0.0f;
  dp->tp.unmri_ll = 0.0f;

  if (first_time) {
    l_mri = parms->l_mri;
    l_unmri = parms->l_unmri;
    l_curv = parms->l_curv;
    l_qcurv = parms->l_qcurv;

    first_time = 0;

    /*      if (!FZERO(l_mri))
            fprintf(WHICH_OUTPUT,"l_mri = %2.2f ", l_mri) ;
            if (!FZERO(l_unmri))
            fprintf(WHICH_OUTPUT,"l_unmri = %2.2f ", l_unmri) ;
            if (!FZERO(l_curv))
            fprintf(WHICH_OUTPUT,"l_curv = %2.2f ", l_curv) ;
            if (!FZERO(l_qcurv))
            fprintf(WHICH_OUTPUT,"l_qcurv = %2.2f ", l_qcurv) ;
            fprintf(WHICH_OUTPUT,"\n") ;*/
  }

  if (!FZERO(l_unmri) && (dp->mri_defect->width <= 5 || dp->mri_defect->height <= 5 || dp->mri_defect->depth <= 5)) {
    l_unmri = 0;
  }

  if (!FZERO(l_mri)) {
    ll += l_mri * mrisComputeDefectMRILogLikelihood(mris, mri, &dp->tp, h_white, h_gray, h_grad, mri_gray_white);
  }
  if (!FZERO(l_unmri)) {
    ll += l_unmri * mrisComputeDefectMRILogUnlikelihood(computeDefectContext, mris, dp, h_border);
  }
  if (!FZERO(l_qcurv)) {
    /*compute the second fundamental form */
    computeDefectSecondFundamentalForm(mris, &dp->tp);
    ll += l_qcurv * mrisComputeDefectCurvatureLogLikelihood(mris, &dp->tp, h_k1, h_k2, mri_k1_k2);
  }
  if (!FZERO(l_curv)) {
    ll += l_curv * mrisComputeDefectNormalDotLogLikelihood(mris, &dp->tp, h_dot);
  }

  l_unmri = parms->l_unmri;
  if (mrisCheckDefectFaces(mris, dp) < 0) ll -= 10000000;

  return (ll);
}


double mrisComputeDefectNormalDotLogLikelihood(MRI_SURFACE *mris, TP *tp, HISTOGRAM *h_dot)
{
  double v_ll, total_ll = 0.0, nx, ny, nz, x, y, z, dx, dy, dz, dot;
  int vno, n, i, bin;

  double t_area, tc_area;

  t_area = tc_area = 0.0;

  /* compute faces only for modified vertices */
  for (i = 0; i < tp->nvertices; i++) {
    vno = tp->vertices[i];
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    x = v->origx;
    y = v->origy;
    z = v->origz;

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;

    for (v_ll = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX * const vn = &mris->vertices[vt->v[n]];
      dx = vn->origx - x;
      dy = vn->origy - y;
      dz = vn->origz - z;
      dot = dx * nx + dy * ny + dz * nz;

      bin = nint((dot - h_dot->bins[0]) / h_dot->bin_size);
      if (bin < 0) {
        bin = 0;
      }
      else if (bin >= h_dot->nbins) {
        bin = h_dot->nbins - 1;
      }

      v_ll += log(h_dot->counts[bin]);
    }
    total_ll += v_ll / vt->vnum;
    t_area += v->area;
    tc_area += v->area * v_ll / vt->vnum;
  }

  if (tp->nvertices) {
    tp->cll = tc_area / t_area;
    tp->curv_ll = total_ll / (double)tp->nvertices;
  }
  else {
    tp->cll = 0.0;
    tp->curv_ll = 0.0;
  }

  return tp->curv_ll;
}


double mrisComputeDefectCurvatureLogLikelihood(
    MRI_SURFACE *mris, TP *tp, HISTOGRAM *h_k1, HISTOGRAM *h_k2, MRI *mri_k1_k2)
{
  double v_ll, total_ll = 0.0, new_total;
  int i, vno, bin, bink1, bink2;
  VERTEX *v;
  double t_area, tc_area;

  new_total = 0.0;

  t_area = tc_area = 0.0;

  /* compute faces only for modified vertices */
  for (i = 0; i < tp->nvertices; i++) {
    vno = tp->vertices[i];
    v = &mris->vertices[vno];

    bin = nint((v->k1 - h_k1->bins[0]) / h_k1->bin_size);
    if (bin < 0) {
      bin = 0;
    }
    else if (bin >= h_k1->nbins) {
      bin = h_k1->nbins - 1;
    }
    v_ll = log(h_k1->counts[bin]);

    total_ll += v_ll;

    bin = nint((v->k2 - h_k2->bins[0]) / h_k2->bin_size);
    if (bin < 0) {
      bin = 0;
    }
    else if (bin >= h_k2->nbins) {
      bin = h_k2->nbins - 1;
    }
    v_ll += log(h_k2->counts[bin]);

    total_ll += v_ll;

    bink1 = MIN(mri_k1_k2->width - 1, MAX(0, (int)((v->k1 - mri_k1_k2->xstart) / mri_k1_k2->xsize)));
    bink2 = MIN(mri_k1_k2->height - 1, MAX(0, (int)((v->k2 - mri_k1_k2->ystart) / mri_k1_k2->ysize)));

    new_total += log(MRIFvox(mri_k1_k2, bink1, bink2, 0));

    t_area += v->area;
    tc_area += v->area * log(MRIFvox(mri_k1_k2, bink1, bink2, 0));
  }
  if (tp->nvertices) {
    new_total /= (double)tp->nvertices;
  }

  if (tp->nvertices) {
    tp->qcurv_ll = total_ll / (double)tp->nvertices;
    tp->qcll = tc_area / t_area;
  }
  else {
    tp->qcll = 0.0;
    tp->qcurv_ll = 0.0;
  }

  tp->qcurv_ll = new_total; /* using mri_k1_k2 */

  return tp->qcurv_ll;
}

static void get_origxyz(VERTEX const * vertex, float* x, float* y, float* z) {
    *x = vertex->origx;
    *y = vertex->origy;
    *z = vertex->origz;
}

static double mrisComputeDefectMRILogUnlikelihood_wkr(
    ComputeDefectContext* computeDefectContext,
    MRI_SURFACE  * const mris_nonconst, 			        // various subcomponents of these structures get updated
    DEFECT_PATCH * const dp_nonconst, 
    HISTOGRAM    * const h_border_nonconst);
    
double mrisComputeDefectMRILogUnlikelihood(
    ComputeDefectContext* computeDefectContext,
    MRI_SURFACE  * const mris_nonconst, 			        // various subcomponents of these structures get updated
    DEFECT_PATCH * const dp_nonconst, 
    HISTOGRAM    * const h_border_nonconst) {

    static int once;
    static int suppress_usecomputeDefectContext = 0;
    if (!once++) {
        if (getenv("FREESURFER_SUPPRESS_using_computeDefectContext")) {
            fprintf(stderr, "Suppressing using computeDefectContext\n");
            suppress_usecomputeDefectContext = 1;
        }
    }
    if (suppress_usecomputeDefectContext) computeDefectContext = NULL;
    
    int saved_noteVnoMovedInActiveRealmTreesCount = noteVnoMovedInActiveRealmTreesCount++;
    
    //  TIMER_INTERVAL_BEGIN(A)

    double result = 
        mrisComputeDefectMRILogUnlikelihood_wkr(
            computeDefectContext,
            mris_nonconst,
            dp_nonconst, 
            h_border_nonconst);

    //  TIMER_INTERVAL_END(A)

    if (0) 
        printf("noteVnoMovedInActiveRealmTrees called:%d\n",
            noteVnoMovedInActiveRealmTreesCount-saved_noteVnoMovedInActiveRealmTreesCount);

    return result;
}


// PerThreadVertexPseudoNormalCache
// Since each vertex is in two or more faces, it makes sense to cache this
// Sizing it will be important - based on test_mris_fix_topology
//
#define PerThreadVertexPseudoNormalCacheSize (1<<12)	    	    	    	    	// below code assumes this is divisible by (sizeof(long) * 8) 
#define PerThreadVertexPseudoNormalCacheMask (PerThreadVertexPseudoNormalCacheSize-1)
typedef struct PerThreadVertexPseudoNormalCacheEntry {
  int   vno;
  float norm[3];
} PerThreadVertexPseudoNormalCacheEntry;

#define PerThreadVertexPseudoNormalCacheInitedSize (PerThreadVertexPseudoNormalCacheSize / (sizeof(long) * 8))

typedef struct PerThreadVertexPseudoNormalCache {
  long                                  inited [PerThreadVertexPseudoNormalCacheInitedSize];
  PerThreadVertexPseudoNormalCacheEntry entries[PerThreadVertexPseudoNormalCacheSize      ];
} PerThreadVertexPseudoNormalCache;

static PerThreadVertexPseudoNormalCache* makePerThreadVertexPseudoNormalCache() {
  PerThreadVertexPseudoNormalCache* p = (PerThreadVertexPseudoNormalCache*)malloc(sizeof(PerThreadVertexPseudoNormalCache));
  bzero(&p->inited, sizeof(p->inited));
  return p;
}

static void freePerThreadVertexPseudoNormalCache(PerThreadVertexPseudoNormalCache** pp) {
  PerThreadVertexPseudoNormalCache* p = *pp;
  *pp = NULL;
  if (!p) return;
  freeAndNULL(p);
}

static void cachedOrComputeVertexPseudoNormal(
  PerThreadVertexPseudoNormalCache** pp,
  MRI_SURFACE const  * const mris,
  int                  const vno,
  float*               const norm,	    	    // output!
  DEFECT_PATCH const * const dp) {
    
  static volatile bool once;
  static bool use_cache, test_cache;
  if (!once) 
#ifdef HAVE_OPENMP
  #pragma omp critical
#endif
  if (!once) {
    use_cache  =  !getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_ComputeVertexPseudoNormalCache_suppress");	
    test_cache = !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_ComputeVertexPseudoNormalCache_test");	
    use_cache |= test_cache;
    if (!use_cache)  
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood not using ComputeVertexPseudoNormalCache\n");
    if (test_cache) 
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood testing ComputeVertexPseudoNormalCache\n");
    once = true;
  }

  static long count, limit = 1, hits, invalidates;
  count++;
  
  PerThreadVertexPseudoNormalCacheEntry* entry = NULL;
  
  bool valid = false;
  if (use_cache) {
    PerThreadVertexPseudoNormalCache* p = *pp;
    if (!p) {
    	*pp = p = (PerThreadVertexPseudoNormalCache*)malloc(sizeof(PerThreadVertexPseudoNormalCache));
	bzero(&p->inited, sizeof(p->inited));
    }

    bool inited = false;

    int index = (vno*1234567) & PerThreadVertexPseudoNormalCacheMask;
    int tries;
    for (tries = 0; tries < 2; tries++) {
      int const initedIndex =       (index / ( sizeof(long) * 8)   );
      int const initedMask  = 1L << (index & ((sizeof(long) * 8)-1));
      
      if (0) {
      	cheapAssert(0 <= index);
      	cheapAssert(index       < PerThreadVertexPseudoNormalCacheSize);
      	cheapAssert((unsigned)initedIndex < PerThreadVertexPseudoNormalCacheInitedSize);
      }
      
      inited = p->inited[initedIndex] & initedMask;
      p->inited[initedIndex] |= initedMask;   	    	// will be inited soon
      
      
      entry = &p->entries[index];
      if (!inited) break;	    	    	    	// use this entry, since currently empty

      valid = inited && (entry->vno == vno);
      if (valid) {  // use this entry, since is the right one
	hits++; 
	break; 
      }
      
      index = (index + 1) & PerThreadVertexPseudoNormalCacheMask;
    }
    
    if (!valid) {	    	    	// either empty, or too many tries made
      if (inited) invalidates++;  	// count evictions
      entry->vno = vno;     	    	// use it
    }
  }
      
  // Fetch from, or update, or test cache
  //
  if (valid && !test_cache) {

    memmove(norm, &entry->norm, sizeof(entry->norm));

  } else {

    computeVertexPseudoNormal(mris, vno, norm, dp->verbose_mode);

    if (valid && memcmp(&entry->norm, norm, sizeof(entry->norm))) {
      static int count;
      if (count++ < 10) {
	fprintf(stdout,"%s:%d PerThreadVertexPseudoNormalCacheEntry wrong norm for vno:%d tid:%d<<<<<<<<<<<<<\n", __FILE__, __LINE__, 
	    vno, omp_get_thread_num());
	fprintf(stdout,"cache:%g %g %g\n", entry->norm[0], entry->norm[1], entry->norm[2]);
	fprintf(stdout,"norm :%g %g %g\n",        norm[0],        norm[1],        norm[2]);
	exit(1);
      }
    }
    
    if (entry) memmove(&entry->norm, norm, sizeof(entry->norm));
  }
  
  if (test_cache) {
    if (count >= limit) {
      if (limit < 100000) limit *= 2; else limit += 100000;
      fprintf(stdout, "cachedOrComputeVertexPseudoNormal count:%g hits:%g invalidates:%g\n", (float)count, (float)hits, (float)invalidates);
    }
  }
}


static void useComputeDefectContextRealmTree(
    ComputeDefectContext* computeDefectContext,
    MRIS  const * const   mris,
    GetXYZ_FunctionType   getXYZ) 
{
    if (computeDefectContext->realmTree == NULL) {
    
#ifdef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
        fprintf(stderr, "%s:%d useComputeDefectContextRealmTree making realmTree\n",__FILE__,__LINE__);
#endif
        computeDefectContext->realmTree = makeRealmTree(mris, getXYZ);
        mrisurf_orig_clock++;
        
        insertActiveRealmTree(mris, computeDefectContext->realmTree, getXYZ);

    } else {
        updateRealmTree(computeDefectContext->realmTree, mris, getXYZ);
#ifdef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
        if (checkRealmTree(computeDefectContext->realmTree, mris, getXYZ)) {
            fprintf(stderr, "%s:%d useComputeDefectContextRealmTree failed checking realmTree\n",__FILE__,__LINE__);
            exit(1);
        }
#endif
    }
}


static double mrisComputeDefectMRILogUnlikelihood_wkr(
    ComputeDefectContext* computeDefectContext,
    MRI_SURFACE  * const mris_nonconst, 			            	// various subcomponents of these structures get updated
    DEFECT_PATCH * const dp_nonconst, 
    HISTOGRAM    * const h_border_nonconst)
{
#if 1
  // The tests themselves are expensive, so eliminate them except when developing
  //
  static const bool 
    keep_sign_bug = false,  	    	    	    	    	    	    	    	// no test
    do_new_MRIDistance = true,     	    	    	    	    	    	    	// if both set, that tests the new
    do_old_MRIDistance = false,
    do_new_loop3 = true,   	    	    	    	    	    	    	    	// if both set, that tests the new
    do_old_loop3 = false, 
    	    	    	            	        use_fast_avoidable_prediction = true,	// no test
    test_avoidable_prediction = false, 	        use_avoidable_prediction = true,
    test_sharedVertexPseudoNormalCache = false, use_sharedVertexPseudoNormalCache = true;
#else
  // Allow comparing different options quickly, but the times are affected compared to unconditional
  //
  static bool 
    keep_sign_bug,  	    	    	    	    	    	    	    	// no test
    do_new_MRIDistance,     	    	    	    	    	    	    	// if both set, that tests the new
    do_old_MRIDistance,
    do_new_loop3,   	    	    	    	    	    	    	    	// if both set, that tests the new
    do_old_loop3, 
    	    	    	            	    use_fast_avoidable_prediction,	// no test
    test_avoidable_prediction,      	    use_avoidable_prediction,
    test_sharedVertexPseudoNormalCache,     use_sharedVertexPseudoNormalCache;
    
  static bool once;
  if (!once) { once = true;
    keep_sign_bug             	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_dont_fix_sign_bug"); 
    do_old_loop3              	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_loop3_old");
    do_new_loop3              	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_loop3_new") || !do_old_loop3;
    do_old_MRIDistance  	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_MRIDistance_old");	
    do_new_MRIDistance  	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_MRIDistance_new") || !do_old_MRIDistance;
    test_avoidable_prediction 	    	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_test_avoidable_prediction");
    use_avoidable_prediction  	    	=  !getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_dont_use_avoidable_prediction") 
    	    	    	    	    	    || test_avoidable_prediction;
    use_fast_avoidable_prediction   	=  !getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_dont_use_fast_avoidable_prediction") 
    	    	    	    	    	    && use_avoidable_prediction;
    
    test_sharedVertexPseudoNormalCache	= !!getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_test_sharedVertexPseudoNormalCache"),
    use_sharedVertexPseudoNormalCache	=  !getenv("FREESURFER_mrisComputeDefectMRILogUnlikelihood_dont_use_sharedVertexPseudoNormalCache")
    	    	    	    	    	    || test_sharedVertexPseudoNormalCache;

    if (do_old_loop3) {
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood using old loop3 algorithm with the %s\n",keep_sign_bug?"sign bug still there":"sign bug fixed");
      if (do_new_loop3) 
        fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood comparing old and new loop3 algorithm\n");
    }
    if (do_old_MRIDistance) {
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood using old mris_distance\n");
      if (do_new_MRIDistance) {
        fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood comparing old and new mris_distance\n");
      }
    }
    if (!use_avoidable_prediction) {
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood using old tests only\n");
    }
    if (!use_fast_avoidable_prediction) {
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood using slower prediction code only\n");
    }
    if (!use_sharedVertexPseudoNormalCache) {
      fprintf(stdout, "mrisComputeDefectMRILogUnlikelihood not using sharedVertexPseudoNormalCache%s\n",
        test_sharedVertexPseudoNormalCache?" only, but testing it":"");
    }
  }
#endif

  MRI_SURFACE  const * const mris     = mris_nonconst;			// find where the modifiers are
  DEFECT_PATCH const * const dp       = dp_nonconst;
  //HISTOGRAM  const * const h_border = h_border_nonconst;		// unused

  MRI          * const mri_distance_nonconst = dp->mri_defect_sign;
  MRI    const * const mri_distance = dp->mri_defect_sign;
  MRI    const * const mri_defect   = dp->mri_defect;
  DEFECT const * const dp_defect    = dp->defect;

  //  TIMER_INTERVAL_BEGIN(getRealmTree)
  if (computeDefectContext) {
    useComputeDefectContextRealmTree(computeDefectContext, mris, get_origxyz);
  }
  //  TIMER_INTERVAL_END(getRealmTree)

  int const maxThreads = omp_get_max_threads();
  if (computeDefectContext) {
    computeDefectContext->mris_deferred_norms = mris_nonconst;  // MODIFIER
    mrisurf_deferSetFaceNorms(mris_nonconst);
  } else {
    mrisurf_recomputeFaceNorms(mris_nonconst);                         // old code 
  }


  /* look at approximately +/- 2mm */
  int const delta = 2.0 * mri_defect->xsize;

/* optimize by listing the concerned faces */

  /* find the distance to each surface voxels */
  {int p;
  

//#define BEVIN_COUNT_EXITS
#ifdef BEVIN_COUNT_EXITS
  int noexits = 0;
  int exit_counters[6]; { int i; for (i = 0; i < 6; i++) exit_counters[i] = 0; }
#endif


  // calculate the tasks to do
  //
  typedef struct Entry {
#define ENTRY_MEMBERS \
	ELT(int, fno) \
	ELT(FACE const *, face) \
	ELT(int, imin_nobnd) \
	ELT(int, imax_nobnd) \
	ELT(int, jmin_nobnd) \
	ELT(int, jmax_nobnd) \
	ELT(int, kmin_nobnd) \
	ELT(int, kmax_nobnd) \
	ELT(float, x0) \
	ELT(float, x1) \
	ELT(float, x2) \
	ELT(float, y0) \
	ELT(float, y1) \
	ELT(float, y2) \
	ELT(float, z0) \
	ELT(float, z1) \
	ELT(float, z2) \
	// end of macro
#define ELT(T,X) T X; 
    ENTRY_MEMBERS
#undef ELT
  } Entry;
  Entry* buffer  = (Entry*)malloc(mris->nfaces * sizeof(Entry));
  int bufferSize = 0; 

  //  TIMER_INTERVAL_BEGIN(taskCreation)

  // the {i,j,k}VOL just inside the following loop implements a linear function of the .orig{x,y,z}
  // so it should be possible to avoid searching all the faces to find the ones that fit in the box
  //
  // jVOL does the following mapping
  //        #define zVOL(mri,z) (mri->zsize*(z-mri->zstart)) 
  //        #define kVOL(mri,z) ((int)(zVOL(mri,z)+0.5))
  // so
  //        kmax_nobnd = ((int)(mri->zsize*(z - mri->zstart))+0.5) + delta
  //
  // we want to ignore when  kmax_nobnd < 0, but for safety, make it -1
  // i.e.   ((int)(mri->zsize*(z - mri->zstart))+0.5) + delta <  -1
  //        ((int)(mri->zsize*(z - mri->zstart))+0.5)         <  -delta - 1                  subtract delta from both sides
  //               mri->zsize*(z - mri->zstart))+0.5          <  -delta - 2                  conversion to int rounds towards zero
  //               mri->zsize*(z - mri->zstart))              <  -delta - 2.5                subtract 0.5 from both sides
  //                           z - mri->zstart                < (-delta - 2.5) / mri->zsize
  //                           z                              < (-delta - 2.5) / mri->zsize + mri->zstart         
  //
  // 
  // we want to ignore when  kmin_nobnd > mri_defect->depth - 1, but for safety, make it -0
  // i.e.        kVOL(mri_defect, z) - delta                  >  mri_defect->depth
  //        ((int)(mri->zsize*(z - mri->zstart))+0.5)         >  mri_defect->depth + delta
  //              (mri->zsize*(z - mri->zstart))+0.5)         >  mri_defect->depth + delta + 1
  //               mri->zsize*(z - mri->zstart)               >  mri_defect->depth + delta + 1.5
  //                           z - mri->zstart                > (mri_defect->depth + delta + 1.5) / mri->zsize
  //                           z                              > (mri_defect->depth + delta + 1.5) / mri->zsize + mri->zstart
  //            
  //
  float const realm_xLo = (0                  - delta - 2.5) / mri_defect->xsize + mri_defect->xstart;
  float const realm_xHi = (mri_defect->width  + delta + 1.5) / mri_defect->xsize + mri_defect->xstart;

  float const realm_yLo = (0                  - delta - 2.5) / mri_defect->ysize + mri_defect->ystart;
  float const realm_yHi = (mri_defect->height + delta + 1.5) / mri_defect->ysize + mri_defect->ystart;

  float const realm_zLo = (0                  - delta - 2.5) / mri_defect->zsize + mri_defect->zstart;
  float const realm_zHi = (mri_defect->depth  + delta + 1.5) / mri_defect->zsize + mri_defect->zstart;

  // Get the list of interesting fno's in ascending order
  // to make this worthwhile, there is a precomputation shared across all the defects
  // This precomputation is the computeDefectContext->realmTree
  //
  Realm* realm =        // This is the realm the face must possibly intersect to be interesting below
    !computeDefectContext
    ? NULL
    : makeRealm(
        computeDefectContext->realmTree, 
        realm_xLo, realm_xHi, 
        realm_yLo, realm_yHi,
        realm_zLo, realm_zHi);
  
  int   const fnosCapacity = !realm ? 0    : realmNumberOfMightTouchFno(realm);
  int * const fnos         = !realm ? NULL : (int*)malloc(fnosCapacity*sizeof(int));
  int   const fnosSize     = !realm ? 0    : realmMightTouchFno(realm, fnos, fnosCapacity);
#if 0
  qsort(fnos, fnosSize, sizeof(int), int_compare);
#else
  sort_int(fnos, fnosSize, true);
#endif
  
  bool const   iterateOverFnos =
#ifndef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
    !!computeDefectContext;
#else
    false;
#endif
  
  int const fnosToIterateOverSize = 
    iterateOverFnos ? fnosSize : mris->nfaces;

#ifdef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
  static int tellCheckingFnosOnce = 0;
    // shared but doesn't need to be locked because worst result is an extra info message
#endif
 
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible) 
#endif
  for (p = 0; p < fnosToIterateOverSize; p++) {
    ROMP_PFLB_begin

    int const fno = iterateOverFnos ? fnos[p] : p;

    FACE const * const face = &mris->faces[fno];

    // calculate three vertices
    float const
      y0 = mris->vertices[face->v[0]].origy,
      y1 = mris->vertices[face->v[1]].origy,
      y2 = mris->vertices[face->v[2]].origy;
    float const
      min_y012 = MIN3(y0, y1, y2);
    int const jmin_nobnd = jVOL(mri_defect, min_y012) - delta;

#ifndef BEVIN_COUNT_EXITS
    if (jmin_nobnd > mri_defect->height - 1) continue;			// hottest
#endif      
    float const
      z0 = mris->vertices[face->v[0]].origz,
      z1 = mris->vertices[face->v[1]].origz,
      z2 = mris->vertices[face->v[2]].origz;
    float const
      min_z012 = MIN3(z0, z1, z2);
    int const kmin_nobnd = kVOL(mri_defect, min_z012) - delta;
#ifndef BEVIN_COUNT_EXITS
    if (kmin_nobnd > mri_defect->depth  - 1) continue;			// hot
#endif
    float const
      x0 = mris->vertices[face->v[0]].origx,
      x1 = mris->vertices[face->v[1]].origx,
      x2 = mris->vertices[face->v[2]].origx;
    float const
      max_x012 = MAX3(x0, x1, x2);
    int const imax_nobnd = iVOL(mri_defect, max_x012) + delta;
#ifndef BEVIN_COUNT_EXITS
    if (imax_nobnd < 0)                     continue;			// hot
#endif

    float const max_y012 = MAX3(y0, y1, y2);
    int const jmax_nobnd = jVOL(mri_defect, max_y012) + delta;
    
    float const max_z012 = MAX3(z0, z1, z2);
    int const kmax_nobnd = kVOL(mri_defect, max_z012) + delta;

    float const min_x012 = MIN3(x0, x1, x2);
    int const imin_nobnd = iVOL(mri_defect, min_x012) - delta;
      
    
    /* find the bounding box */


    // TO BE CHECKED
    /* we don't count faces that extend outside the volume -
       should not change the sign */
    // TO BE CHECKED it some defects are close from each other!
    
#ifdef BEVIN_COUNT_EXITS
    if (imin_nobnd > mri_defect->width  - 1) exit_counters[0]++;
    if (jmin_nobnd > mri_defect->height - 1) exit_counters[1]++;	// hottest
    if (kmin_nobnd > mri_defect->depth  - 1) exit_counters[2]++;	// hot
    if (imax_nobnd < 0) exit_counters[3]++;				// hot
    if (jmax_nobnd < 0) exit_counters[4]++;
    if (kmax_nobnd < 0) exit_counters[5]++;
#endif
    
    if (imin_nobnd > mri_defect->width  - 1 
#ifndef BEVIN_COUNT_EXITS
     || jmin_nobnd > mri_defect->height - 1 
     || kmin_nobnd > mri_defect->depth  - 1 
     || imax_nobnd < 0 
#endif
     || jmax_nobnd < 0 
     || kmax_nobnd < 0) {
      ROMP_PF_continue;
    }

#ifdef BEVIN_COUNT_EXITS
    noexits++;
#endif

#ifdef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
    // Make sure fno was in the list of interesting fnos
    //
    if (fnos) {
        if (!tellCheckingFnosOnce++) fprintf(stderr,"Checking all relevant fnos were found\n");
        
        if (!bsearch(&fno, fnos, fnosSize, sizeof(int), int_compare)) {
            
            fprintf(stdout,"%s:%d fno:%d not in the interesting list\n", __FILE__, __LINE__, fno);
            fprintf(stdout,"  face  x:%f..%f y:%f..%f z:%f..%f\n",min_x012,max_x012,min_y012,max_y012,min_z012,max_z012);
            fprintf(stdout,"  realm x:%f..%f y:%f..%f z:%f..%f\n",realm_xLo,realm_xHi,realm_yLo,realm_yHi,realm_zLo,realm_zHi);
            fprintf(stdout,"  fnosSize:%d\n", fnosSize);
            
            // info leading to the realm_zLo bound
            if (1) {
                fprintf(stdout,"0 - delta:%f + 1) / mri_defect->zsize:%f + mri_defect->zstart:%f     realm_zLo:%f\n",
                    (float)(delta), (float)(mri_defect->zsize), (float)(mri_defect->zstart), realm_zLo);

                fprintf(stdout,"kmax_nobnd:%d ok if >= 0\n", kmax_nobnd);
                fprintf(stdout,"   kVOL(mri_defect, max_z012:%f):%d + delta:%d\n", max_z012, kVOL(mri_defect, max_z012), delta);
                fprintf(stdout,"   ((int)(zVOL(mri_defect,max_z012:%f):%f + 0.5))\n", max_z012, zVOL(mri_defect,max_z012));
                fprintf(stdout,"   (mri_defect->zsize:%f * (z:%f - mri_defect->zstart:%f))\n", mri_defect->zsize, max_z012, mri_defect->zstart);
            }

            // info leading to the realm_zHi bound
            if (0) {
                fprintf(stdout,"mri_defect->depth:%f - 1 + delta:%f - 1) / mri_defect->zsize:%f + mri_defect->zstart:%f\n",
                    (float)(mri_defect->depth), (float)(delta), (float)(mri_defect->zsize), (float)(mri_defect->zstart));
            }

            if (1) {
                summarizeRealmTreeFno(computeDefectContext->realmTree, fno);
                summarizeRealmTree   (computeDefectContext->realmTree);
            }
            
            *(int*)(-1) = 0;
        }
    }
#endif
    
    // Add an entry for this fno
    //
    Entry* entry;
#ifdef HAVE_OPENMP
    #pragma omp critical
#endif
    entry = &buffer[bufferSize++];

#define ELT(T,X) entry->X = X; 
    ENTRY_MEMBERS
#undef ELT

    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  //  TIMER_INTERVAL_END(taskCreation)
  //  TIMER_INTERVAL_BEGIN(realmFree)
  
  if (!iterateOverFnos) { 
    ; // printf("Only searching %d fnos, instead of %d to create %d tasks\n", mris->nfaces, mris->nfaces, bufferSize);
  } else {
    ; // printf("Only searching %d fnos, instead of %d to create %d tasks\n", fnosSize,     mris->nfaces, bufferSize);

#ifdef mrisComputeDefectMRILogUnlikelihood_CHECK_USE_OF_REALM
    if (fnosSize > mris->nfaces/4) 
    #pragma omp critical
    {
      float rtXLo, rtXHi, rtYLo, rtYHi, rtZLo, rtZHi;
      getRealmTreeBnds(
        computeDefectContext->realmTree, &rtXLo, &rtXHi, &rtYLo, &rtYHi, &rtZLo, &rtZHi);
      printf("%s:%d Not many filtered when delta:%d\n", __FILE__, __LINE__, delta);
      printf("  rtree x:%8.2f..%8.2f y:%8.2f..%8.2f z:%8.2f..%8.2f\n",rtXLo, rtXHi, rtYLo, rtYHi, rtZLo, rtZHi);
      printf("  realm x:%8.2f..%8.2f y:%8.2f..%8.2f z:%8.2f..%8.2f\n",realm_xLo,realm_xHi,realm_yLo,realm_yHi,realm_zLo,realm_zHi);
      summarizeRealmTree(computeDefectContext->realmTree);
      printf("%s:%d exit(1) called\n", __FILE__, __LINE__);
      exit(1);
    }
#endif

    freeRealm(&realm);
  }

  
#ifdef BEVIN_COUNT_EXITS
  if (1) { 
    fprintf(stderr, "exit counters: ");
    int i; for (i= 0; i < 6; i++) fprintf(stderr,"%d:%d ", i, exit_counters[i]);
    fprintf(stderr," noexits:%d\n",noexits);
  }
#endif

  //  TIMER_INTERVAL_END(realmFree)

  //  TIMER_INTERVAL_BEGIN(taskExecution)

  // The following is a very complex loop
  // We need to understand its behavior better before starting to optimize it
  //
  bool traceTid0 = false;
  if (0) {
    static long count, limit=1, sumBufferSize;
    count++;
    sumBufferSize += bufferSize;
    if (count >= limit) {
      limit *= 2;
      fprintf(stdout, "%s:%d count:%ld avg bufferSize:%g\n", __FILE__, __LINE__, count, (float)sumBufferSize/(float)count);
      traceTid0 = (count == 1024);
    }
  }
  
  
  // do the tasks
  //
  PerThreadVertexPseudoNormalCache* perThreadVertexPseudoNormalCache[_MAX_FS_THREADS];
  PerThreadMRIDistance*             perThreadMRIDistances           [_MAX_FS_THREADS];
  { int tid; 
    for (tid = 0; tid < maxThreads; tid++) {
      perThreadMRIDistances           [tid] = NULL; 	// allocate in the loop to get in the right cache and to parallelize the init
      perThreadVertexPseudoNormalCache[tid] = NULL;
    }
  }
  
  if (do_old_MRIDistance) 
  { int i,j,k;
    for (i = 0; i < mri_distance->width; i++) {
      for (j = 0; j < mri_distance->height; j++) {
        for (k = 0; k < mri_distance->depth; k++) {
          MRIFvox(mri_distance_nonconst, i, j, k) = NPY;		// MODIFIER THAT WASN'T DETECTED
        }
      }
    }
  }

  typedef struct Float3 { float elts[3]; } Float3;
  Float3* sharedVertexPseudoNormalCache = NULL;
  if (use_sharedVertexPseudoNormalCache) {
  
    sharedVertexPseudoNormalCache = (Float3*)malloc(sizeof(Float3) * mris->nvertices);

    volatile char* done = (char*)calloc(mris->nvertices,sizeof(char)); 
    
    ROMP_PF_begin
    int bufferIndex;
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(shown_reproducible) schedule(guided)
#endif
    for (bufferIndex = 0; bufferIndex < bufferSize; bufferIndex++) {
      Entry const * const entry = &buffer[bufferIndex];
      FACE  const * const face  = entry->face;
      int vi;
      for (vi = 0; vi < 3; vi++) {
        int const vno = face->v[vi];
        char alreadyDone = done[vno];
        if (alreadyDone) continue;

#ifdef HAVE_OPENMP
#if GCC_VERSION >= 50400
        #pragma omp atomic capture
#else
        #pragma omp critical
#endif
#endif
        {  alreadyDone = done[vno]; 
           done[vno] = 1; 
        }
        if (alreadyDone) continue;

        computeVertexPseudoNormal(mris, vno, sharedVertexPseudoNormalCache[vno].elts, dp->verbose_mode);
      }
    }
    ROMP_PF_end
    
    free((char*)done);
  }

  ROMP_PF_begin

  int bufferIndex;
#ifdef HAVE_OPENMP
  int const chunkSize = MAX(1,bufferSize/MAX(omp_get_max_threads()*4,32));
  #pragma omp parallel for if_ROMP(shown_reproducible) schedule(dynamic, chunkSize )
#endif
  for (bufferIndex = 0; bufferIndex < bufferSize; bufferIndex++) {
    ROMP_PFLB_begin
    int const tid = omp_get_thread_num();
    
    bool const trace = false; // traceTid0 && (tid == 0);

    Entry const * entry = &buffer[bufferIndex];

#define ELT(T,X) T const X = entry->X;
    ENTRY_MEMBERS
#undef ELT

    long const imin = MAX(imin_nobnd, 0);
    long const imax = MIN(imax_nobnd, mri_defect->width - 1);

    long const jmin = MAX(jmin_nobnd, 0);
    long const jmax = MIN(jmax_nobnd, mri_defect->height - 1);

    long const kmin = MAX(kmin_nobnd, 0);
    long const kmax = MIN(kmax_nobnd, mri_defect->depth - 1);

    // Get the pseudonormals
    // My previous belief was they would not many would be needed by more than one threads, but maybe they are
    // I should measure this!  Or maybe we should use a snoopy cache!
    //
    float n_v0[3], n_v1[3], n_v2[3];
    
    if (!use_sharedVertexPseudoNormalCache || test_sharedVertexPseudoNormalCache) {
      cachedOrComputeVertexPseudoNormal(&perThreadVertexPseudoNormalCache[tid], mris, face->v[0], n_v0, dp);
      cachedOrComputeVertexPseudoNormal(&perThreadVertexPseudoNormalCache[tid], mris, face->v[1], n_v1, dp);
      cachedOrComputeVertexPseudoNormal(&perThreadVertexPseudoNormalCache[tid], mris, face->v[2], n_v2, dp);
    } 

    if (use_sharedVertexPseudoNormalCache) {
      if (test_sharedVertexPseudoNormalCache) {
    	if (memcmp(n_v0, sharedVertexPseudoNormalCache[face->v[0]].elts, sizeof(n_v0)) 
	||  memcmp(n_v1, sharedVertexPseudoNormalCache[face->v[1]].elts, sizeof(n_v1))
	||  memcmp(n_v2, sharedVertexPseudoNormalCache[face->v[2]].elts, sizeof(n_v2))) {
	  fprintf(stdout, "%s:%d sharedVertexPseudoNormalCache wrong\n", __FILE__, __LINE__);
	  float* e0 = sharedVertexPseudoNormalCache[face->v[0]].elts;
	  float* e1 = sharedVertexPseudoNormalCache[face->v[1]].elts;
	  float* e2 = sharedVertexPseudoNormalCache[face->v[2]].elts;
	  fprintf(stdout, " v0 %g %g %g  ?=  %g %g %g \n", n_v0[0], n_v0[1], n_v0[2], e0[0], e0[1], e0[2]);
	  fprintf(stdout, " v1 %g %g %g  ?=  %g %g %g \n", n_v1[0], n_v1[1], n_v1[2], e1[0], e1[1], e1[2]);
	  fprintf(stdout, " v2 %g %g %g  ?=  %g %g %g \n", n_v2[0], n_v2[1], n_v2[2], e2[0], e2[1], e2[2]);
	  exit(1);
	}
      } else {
    	memmove(n_v0, sharedVertexPseudoNormalCache[face->v[0]].elts, sizeof(n_v0)); 
	memmove(n_v1, sharedVertexPseudoNormalCache[face->v[1]].elts, sizeof(n_v1));
	memmove(n_v2, sharedVertexPseudoNormalCache[face->v[2]].elts, sizeof(n_v2));
      }
    } 
    
    // Generate the edges and their normals
    //
    float n_f [3];
    float   n0[3],   n1[3],   n2[3];
    float   e0[3],   e1[3],   e2[3];
    float n_e0[3], n_e1[3], n_e2[3];
    {
      // Do these first, so the optimizer has plenty of common subexpressions knowing the pointed-to items aren't changing
      //
      int fnoV01 = findOtherEdgeFace(mris, fno, face->v[0], face->v[1]);
      int fnoV12 = findOtherEdgeFace(mris, fno, face->v[1], face->v[2]);
      int fnoV20 = findOtherEdgeFace(mris, fno, face->v[2], face->v[0]);

      if (trace) {
      	fprintf(stdout, "  fno:%d abuts fnos:%d %d %d\n", fno, fnoV01, fnoV12, fnoV20);
      }
      
      FaceNormCacheEntry const * fNorm  = getFaceNorm(mris, fno);
      FaceNormCacheEntry const * fNorm0 = getFaceNorm(mris, fnoV01);
      FaceNormCacheEntry const * fNorm1 = getFaceNorm(mris, fnoV12);
      FaceNormCacheEntry const * fNorm2 = getFaceNorm(mris, fnoV20);

      n_f[0] = fNorm->nx;
      n_f[1] = fNorm->ny;
      n_f[2] = fNorm->nz;

      /* edge0: x0 <--> x1 */
      e0[0] = x1 - x0;
      e0[1] = y1 - y0;
      e0[2] = z1 - z0;
    
      F_CROSS(n_f, e0, n0);

      n_e0[0] = fNorm->nx + fNorm0->nx;
      n_e0[1] = fNorm->ny + fNorm0->ny;
      n_e0[2] = fNorm->nz + fNorm0->nz;

      /* edge1: x1 <--> x2 */
      e1[0] = x2 - x1;
      e1[1] = y2 - y1;
      e1[2] = z2 - z1;

      F_CROSS(n_f, e1, n1);

      n_e1[0] = fNorm->nx + fNorm1->nx;
      n_e1[1] = fNorm->ny + fNorm1->ny;
      n_e1[2] = fNorm->nz + fNorm1->nz;

      /* edge2: x2 <--> x0 */
      e2[0] = x0 - x2;
      e2[1] = y0 - y2;
      e2[2] = z0 - z2;

      F_CROSS(n_f, e2, n2);

      n_e2[0] = fNorm->nx + fNorm2->nx;
      n_e2[1] = fNorm->ny + fNorm2->ny;
      n_e2[2] = fNorm->nz + fNorm2->nz;

    }

    float const SQR3_e0 = SQR3(e0);
    float const SQR3_e1 = SQR3(e1);
    float const SQR3_e2 = SQR3(e2);


    if (trace) {
      fprintf(stdout, "  computeVertexPseudoNormal for vnos:%d %d %d\n", face->v[0], face->v[1], face->v[2]);
    }

    if (0) {
    	// Find out the sizes of these boxes
    	static int count, limit = 128;
	if (count++ > limit) 
    	#pragma omp critical
	{
	    limit *= 2;
	    float 
	    	xmin = xSURF(mri_defect, imin),
	    	xmax = xSURF(mri_defect, imax),
	    	ymin = xSURF(mri_defect, jmin),
	    	ymax = xSURF(mri_defect, jmax),
	    	zmin = xSURF(mri_defect, kmin),
	    	zmax = xSURF(mri_defect, kmax);
	    fprintf(stdout, "%s:%d box %g..%g,%g..%g,%g..%g %g*%g*%g\n", __FILE__, __LINE__, 
	    	xmin,xmax,ymin,ymax,zmin,zmax,xmax-xmin,ymax-ymin,zmax-zmin);
	    fprintf(stdout, "%ld..%ld,%ld..%ld,%ld..%ld %ld*%ld*%ld=%ld\n",
	    	imin,imax,jmin,jmax,kmin,kmax,imax+1-imin,  jmax+1-jmin,  kmax+1-kmin,
		                             (imax+1-imin)*(jmax+1-jmin)*(kmax+1-kmin));
	}
    }
    
    // Bevin has an ideas for making the following faster.
    // It is looping over a box, and for each point in the box deciding...
    // 
    // (a) which aspect of the face this point is closest to
    //	    - the face
    //	    - one of the six corners
    //	    - one of the six edges
    //
    // (b) what the distance is to that aspect
    //
    // (c) the minimum of the distance to that aspect of this face and to all the similar situation on all the other faces
    //
    // However there are two facts that can speed this up
    //
    // 1.   It is possible to get a lower bound on the distance to this face, and - when it is too large - don't bother refining it
    //	    This estimate can be based on the average of the vertices, and the distance from it to the furtherest vertex.
    //	    because we know
    //	    	let IJK be any point in the box, and FP be any point on the face
    //	    	    distance(IJK, Center) <= distance(IJK,FP) + distance(FP,Center) 	hence
    //	    	    distance(IJK, Center) - distance(FP,Center) <= distance(IJK,FP)
    //	    so when the LHS is larger than the current least distance for IJK, it is not worth further consideration of this face 
    //
    // 2.   If two points in the box are nearest to the same aspect, then all points between them are also nearest to it
    //      which avoids deciding which of the branches below to take.  This is not yet exploited
    //
    // In addition an initial estimate of the minimum distance can be obtained from the realm tree, which will ignore even more faces.
    //	    TBD
    //
    // Calculate the center of the face - note, the precise location is not important.
    //    
    float cx = 0, cy = 0, cz = 0;
    {
    	int vi;
	for (vi = 0; vi < 3; vi++) {
	    VERTEX const * v = &mris->vertices[face->v[vi]];
	    cx += v->origx;
	    cy += v->origy;
	    cz += v->origz;
	}
	cx *= 0.33f; cy *= 0.33f; cz *= 0.33f;
    }
    
    // Calculate the square of the furtherest distance
    //
    float furtherestVertexDistanceSquared = 0;
    {
    	int vi;
	for (vi = 0; vi < 3; vi++) {
	    VERTEX const * v = &mris->vertices[face->v[vi]];
	    float x = cx - v->origx;
	    float y = cy - v->origy;
	    float z = cz - v->origz;
	    furtherestVertexDistanceSquared = MAX(furtherestVertexDistanceSquared, x*x + y*y + z*z);
	}
    }
    float const furtherestVertexDistance = sqrtf(furtherestVertexDistanceSquared) * 1.01f;  // margin for error included
    
    // For each point in the box, find the distance, and update the leasts accordingly
    //
    // GCC does not find many of the loop invariants in the following, hence the explicit hoisting
    //
    float vec[3], vec0[3], vec1[3], vec2[3];

    PerThreadMRIDistance* ptd = perThreadMRIDistances[tid];
    if (do_new_MRIDistance && !ptd) ptd = perThreadMRIDistances[tid] = makePerThreadMRIDistance(mri_distance);

    long k,j,i;

    long   jToYMapSize = jmax - jmin + 1;
    float  jToYMapBuffer            [128];
    float* jToYMap = (jToYMapSize <= 128) ? jToYMapBuffer : (float*)malloc(jToYMapSize * sizeof(float));
    for (j = jmin; j <= jmax; j++) jToYMap[j - jmin] = ySURF(mri_defect, j);
    
    long   kToZMapSize = kmax - kmin + 1;
    float  kToZMapBuffer            [128];
    float* kToZMap = (kToZMapSize <= 128) ? kToZMapBuffer : (float*)malloc(kToZMapSize * sizeof(float));
    for (k = kmin; k <= kmax; k++) kToZMap[k - kmin] = zSURF(mri_defect, k);
    
    long   kToDoBuffer[128];
    long*  kToDo = (kToZMapSize <= 128) ? kToDoBuffer : (long*)malloc(kToZMapSize * sizeof(long));
    long   kToDoSize = kToZMapSize;

    for (k = kmin; k <= kmax; k++) kToDo[k - kmin] = k;
        //
        // By putting these in a buffer, I hope to avoid a lot of branch mispredicts that were resulting
        // in lots of cancelled floating point ops below
           
    for (i = imin; i <= imax; i++) {
      	  float const x = xSURF(mri_defect, i);
          vec0[0] = x - x0;
          vec1[0] = x - x1;
          vec2[0] = x - x2;
          vec [0] = (vec0[0] + vec1[0] + vec2[0]) * 0.3333333f;

      for (j = jmin; j <= jmax; j++) {
          float const y = jToYMap[j - jmin];
          vec0[1] = y - y0;
          vec1[1] = y - y1;
          vec2[1] = y - y2;
          vec [1] = (vec0[1] + vec1[1] + vec2[1]) * 0.3333333f;

        float const partialPythagorasSum = squaref(x-cx) + squaref(y-cy);
	
        float* const ijkDistances =
            do_new_MRIDistance ? perThreadMRIDistanceElt(ptd, i,j,0) : &MRIFvox(mri_distance_nonconst, i, j, 0);
                  
        if (use_fast_avoidable_prediction && !test_avoidable_prediction) {

	    // This path is so important I have hand-optimized it.
            //
	    // It skips over all the predictedIrrelevant k's
	    // This path does not yet have a specific correctness test
            //
            // By moving it out of the following loop, there may be more overlapping of the floating point operations and the memory traffic
            // and fewer branch mispredicts
            // and more opportunities for unrolling and for vector operations and for cse's
	    //
            kToDoSize = 0;
            for (k = kmin; k <= kmax; k++) {
                kToDo[kToDoSize] = k;
                float zt                    = kToZMap[k - kmin];
                float ijkDistance           = ijkDistances[k];
                float distanceToCxyzSquared = partialPythagorasSum + squaref(zt-cz);
                kToDoSize += ((distanceToCxyzSquared*0.98f <= (squaref(fabsf(ijkDistance) + furtherestVertexDistance))));
            }
        }
        
        long kToDoIndex;
        for (kToDoIndex = 0; kToDoIndex < kToDoSize; kToDoIndex++) {
          k = kToDo[kToDoIndex];

    	  float const z = kToZMap[k - kmin];
	  
    	  // Here is the minimum distance to update
	  //
    	  float* const ijkDistanceElt = &ijkDistances[k];
	  
    	  // The following code has four different ways of getting the answer
	  // and can check them against each other
	  //	      
	  float distanceToCxyzSquared = 0.0; bool predictedIrrelevant = false;
	    
      	  static long pointCount, pointLimit = 1000000, irrelevantPointCount;
#if 0
	  pointCount++;
#endif

          // Calculate the distance to center, and ignore when no chance of providing a new least distance
	  //
    	  if (use_avoidable_prediction && test_avoidable_prediction) {
	      distanceToCxyzSquared = partialPythagorasSum + squaref(z-cz);
	      predictedIrrelevant   =  distanceToCxyzSquared*0.98f > squaref(fabsf(*ijkDistanceElt) + furtherestVertexDistance);
	      if (predictedIrrelevant) irrelevantPointCount++;
	      if (!test_avoidable_prediction && predictedIrrelevant) 
	      	continue;
	  }

       	  float new_distanceSquared = -1.0f;
          float new_distanceSignArg = 0.0f;
          float new_distance        = 0.0f;
	  float old_distance        = 0.0f;

    	  // Complete the three vectors
	  //
          vec0[2] = z - z0;
          vec1[2] = z - z1;
          vec2[2] = z - z2;
	  
          if (true) {
	  
            vec[2] = (vec0[2] + vec1[2] + vec2[2]) * 0.3333333f;

            /* compute distance to face */
            /* where is the point */
            float const
	      val0 = F_DOT(vec0, n0),
              val1 = F_DOT(vec1, n1),
              val2 = F_DOT(vec2, n2);

    	    if (do_new_loop3) {

              if ((val0 >= 0) && (val1 >= 0) && (val2 >= 0)) {

		// the projection of the vertex is inside
        	//
		new_distance = F_DOT(n_f, vec); 	// n_f is already normalized

              } else {

		float  least_distance_squared = NPY*NPY;
		float* least_sign_lhs = NULL;
		float* least_sign_rhs = NULL;
#define THIS_CASE_DEF       	float* this_sign_lhs,*this_sign_rhs;
#define LEAST_CASE(LHS,RHS)     { least_sign_lhs = (LHS); least_sign_rhs = (RHS); }
#define THIS_CASE(LHS,RHS,CASE) { this_sign_lhs  = (LHS); this_sign_rhs = (RHS); }
#define MAKE_LEAST_THIS         { least_sign_lhs = this_sign_lhs; least_sign_rhs = this_sign_rhs; }

        	if (val0 <= 0) {
        	  /* compute distance to edge0 */
		  float distance_squared;
        	  float val = F_DOT(vec0, e0);
	          THIS_CASE_DEF
        	  if (val < 0) {
                    /* closer to x0 */
		    THIS_CASE(n_v0, vec0, 2);
                    distance_squared = SQR3(       vec0);
        	  }
        	  else if (val < SQR3_e0) {
                    /* closer to edge0 */
		    THIS_CASE(n_e0, vec0, 3);
                    distance_squared = MAX(0, SQR3(vec0) - SQR(val) / SQR3_e0);
        	  }
        	  else {
                    /* closer to x1 */
		    THIS_CASE(n_v1, vec1, 4);
                    distance_squared = SQR3 (vec1);
        	  }
		  least_distance_squared = distance_squared;
	    	  MAKE_LEAST_THIS;
        	}

        	if (val1 <= 0) {
		  float distance_squared;
        	  float val = F_DOT(vec1, e1);
	          THIS_CASE_DEF
        	  if (val < 0) {
                    /* closer to x1 */
		    THIS_CASE(n_v1, vec1, 5);
                    distance_squared = SQR3(       vec1);
        	  }
        	  else if (val < SQR3_e1) {
                    /* closer to edge1 */
		    THIS_CASE(n_e1, vec1, 6);
                    distance_squared = MAX(0, SQR3(vec1) - SQR(val) / SQR3_e1);
        	  }
        	  else {
                    /* closer to x2 */
  		    THIS_CASE(n_v2, vec2, 7);
                    distance_squared = SQR3(vec2);
        	  }
		  if (least_distance_squared > distance_squared) {
		    least_distance_squared = distance_squared;
		    MAKE_LEAST_THIS;
		  }
        	}

        	if (val2 <= 0) {
		  float distance_squared;
        	  float val = F_DOT(vec2, e2);
	          THIS_CASE_DEF
        	  if (val < 0) {
                    /* closer to x2 */
		    THIS_CASE(n_v2, vec2, 8);
                    distance_squared = SQR3(       vec2);
        	  }
        	  else if (val < SQR3(e2)) {
                    /* closer to edge2 */
		    THIS_CASE(n_e2, vec2, 9);
                    distance_squared = MAX(0, SQR3(vec2) - SQR(val) / SQR3_e2);
        	  }
        	  else {
                    /* closer to x0 */
		    THIS_CASE(n_v0, vec0, 10);
                    distance_squared = SQR3(vec0);
        	  }
		  if (least_distance_squared > distance_squared) {
		    least_distance_squared = distance_squared;
		    MAKE_LEAST_THIS;
		  }
        	}
#undef THIS_CASE_DEF
#undef LEAST_CASE
#undef THIS_CASE
#undef MAKE_LEAST_THIS

    		float least_sign = F_DOT(least_sign_lhs, least_sign_rhs);
		new_distanceSignArg = least_sign;
                new_distanceSquared = least_distance_squared;
              }
    	    }
    	    // end if (do_new_loop3)
	  
            if (do_old_loop3) {
	      if (keep_sign_bug) {
		float val, valu, sign, distance;

		// THIS CODE HAS A FUNDAMENTAL PROBLEM
		//
		// IT FINDS THE LEAST ABS DISTANCE, BUT THE SIGN OF THE DISTANCE IS THE SIGN OF THE LAST TESTED VALUE
		// WHICH CAN'T POSSIBLY BE RIGHT!  CONCLUSION - MAYBE SIGN DOESNT MATTER?  IT IS EXPENSIVE TO COMPUTE...

        	if ((val0 >= 0) && (val1 >= 0) && (val2 >= 0)) {
        	  /* the projection of the vertex is inside */
        	  val = F_DOT(n_f, vec);
        	  valu     = 1;
        	  sign     = val;
        	  distance = val; /* n_f is already normalized */
        	}
        	else {
        	  distance = NPY;
        	  sign = 0;
        	  valu = 0;

        	  if (val0 <= 0) {
        	    /* compute distance to edge0 */
        	    val = F_DOT(vec0, e0);
        	    if (val < 0) {
                      /* closer to x0 */
                      sign = F_DOT(n_v0, vec0);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
        	    }
        	    else if (val < SQR3(e0)) {
                      /* closer to edge0 */
                      sign = F_DOT(n_e0, vec0);
                      valu = 3;
                      distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec0) - SQR(val) / SQR3(e0))));
        	    }
        	    else {
                      /* closer to x1 */
                      sign = F_DOT(n_v1, vec1);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
        	    }
        	  };

        	  if (val1 <= 0) {
        	    val = F_DOT(vec1, e1);
        	    if (val < 0) {
                      /* closer to x1 */
                      sign = F_DOT(n_v1, vec1);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
        	    }
        	    else if (val < SQR3(e1)) {
                      /* closer to edge1 */
                      sign = F_DOT(n_e1, vec1);
                      valu = 3;
                      distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec1) - SQR(val) / SQR3(e1))));
        	    }
        	    else {
                      /* closer to x2 */
                      sign = F_DOT(n_v2, vec2);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
        	    }
        	  };

        	  if (val2 <= 0) {
        	    val = F_DOT(vec2, e2);
        	    if (val < 0) {
                      /* closer to x2 */
                      sign = F_DOT(n_v2, vec2);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
        	    }
        	    else if (val < SQR3(e2)) {
                      /* closer to edge2 */
                      sign = F_DOT(n_e2, vec2);
                      valu = 3;
                      distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec2) - SQR(val) / SQR3(e2))));
        	    }
        	    else {
                      /* closer to x0 */
                      sign = F_DOT(n_v0, vec0);
                      valu = 2;
                      distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
        	    }
        	  };
        	}

    		old_distance = distance;
              } else {
		float val, valu, sign, distance;

		// THIS CODE KEEPS THE SIGN FROM THE LEAST OF THE DISTANCES

        	if ((val0 >= 0) && (val1 >= 0) && (val2 >= 0)) {
        	  /* the projection of the vertex is inside */
        	  val = F_DOT(n_f, vec);
        	  valu     = 1;
        	  sign     = val;
        	  distance = val; /* n_f is already normalized */
        	}
        	else {
        	  distance = NPY;
        	  sign = 0;
        	  valu = 0;

    		  float trialDistance;

        	  if (val0 <= 0) {
        	    /* compute distance to edge0 */
        	    val = F_DOT(vec0, e0);
        	    if (val < 0) {
                      /* closer to x0 */
                      sign = F_DOT(n_v0, vec0);
                      valu = 2;
		      trialDistance = NORM3(vec0);
        	    }
        	    else if (val < SQR3(e0)) {
                      /* closer to edge0 */
                      sign = F_DOT(n_e0, vec0);
                      valu = 3;
		      trialDistance = sqrt(MAX(0, SQR3(vec0) - SQR(val) / SQR3(e0)));
        	    }
        	    else {
                      /* closer to x1 */
                      sign = F_DOT(n_v1, vec1);
                      valu = 2;
		      trialDistance = NORM3(vec1);
        	    }
                    if (fabs(distance) > trialDistance) distance = SIGN(sign) * trialDistance;
        	  };

        	  if (val1 <= 0) {
        	    val = F_DOT(vec1, e1);
        	    if (val < 0) {
                      /* closer to x1 */
                      sign = F_DOT(n_v1, vec1);
                      valu = 2;
                      trialDistance = NORM3(vec1);
        	    }
        	    else if (val < SQR3(e1)) {
                      /* closer to edge1 */
                      sign = F_DOT(n_e1, vec1);
                      valu = 3;
                      trialDistance = sqrt(MAX(0, SQR3(vec1) - SQR(val) / SQR3(e1)));
        	    }
        	    else {
                      /* closer to x2 */
                      sign = F_DOT(n_v2, vec2);
                      valu = 2;
                      trialDistance = NORM3(vec2);
        	    }
    	    	    if (fabs(distance) > trialDistance) distance = SIGN(sign) * trialDistance;
        	  };

        	  if (val2 <= 0) {
        	    val = F_DOT(vec2, e2);
        	    if (val < 0) {
                      /* closer to x2 */
                      sign = F_DOT(n_v2, vec2);
                      valu = 2;
                      trialDistance = NORM3(vec2);
        	    }
        	    else if (val < SQR3(e2)) {
                      /* closer to edge2 */
                      sign = F_DOT(n_e2, vec2);
                      valu = 3;
                      trialDistance = sqrt(MAX(0, SQR3(vec2) - SQR(val) / SQR3(e2)));
        	    }
        	    else {
                      /* closer to x0 */
                      sign = F_DOT(n_v0, vec0);
                      valu = 2;
                      trialDistance = NORM3(vec0);
        	    }
    	    	    if (fabs(distance) > trialDistance) distance = SIGN(sign) * trialDistance;
        	  };
        	}

    		old_distance = distance;
              }
	    }
            // end if (do_old_loop3)

    	  }
          // end of if (true)

    	  // Show the stats about how each point was calculated
	  //
	  if (false && (pointCount >= pointLimit)) {
	    if (pointLimit > 10000000) pointLimit += 10000000; else pointLimit *= 2;
	    fprintf(stdout, "pointCount:%g irrelevantPointCount:%g remainder:%g\n", 
	      (float) pointCount, (float)irrelevantPointCount,
	      (float)(pointCount -       irrelevantPointCount));
	  }

	  // Compare the various answers when testing
	  //
    	  if (test_avoidable_prediction && predictedIrrelevant) {
            if (new_distanceSquared >= 0.0f) new_distance = SIGN(new_distanceSignArg) * sqrtf(new_distanceSquared);
	    float distance = do_new_loop3 ? new_distance : old_distance;

     	    if (fabsf(*ijkDistanceElt) > fabsf(distance)) {
	      fprintf(stdout, "%s:%d prediction failed!\n", __FILE__, __LINE__);
	      fprintf(stdout, "do_new_MRIDistance:%d\n", do_new_MRIDistance);
	      fprintf(stdout, "c (%g, %g, %g)\n", cx,cy,cz);
	      fprintf(stdout, "furtherestVertexDistance %g\n", furtherestVertexDistance);
	      fprintf(stdout, "distanceToCxyz %g\n", sqrt(distanceToCxyzSquared));
	      fprintf(stdout, "distanceToCxyz - furtherestVertexDistance:%g\n", sqrt(distanceToCxyzSquared) - furtherestVertexDistance);
	      fprintf(stdout, "ijk (%g, %g, %g)\n", x,y,z);
	      fprintf(stdout, "fabsf(distance) %g\n", fabsf(distance));
	      fprintf(stdout, "fabsf(*ijkDistanceElt) %g\n\n", fabsf(*ijkDistanceElt));
	      exit(1);
	    }
	  }
	  
	  if (do_new_loop3 && do_old_loop3) {
            if (new_distanceSquared >= 0.0f) new_distance = SIGN(new_distanceSignArg) * sqrtf(new_distanceSquared);
	  
	    // Sadly the old code compares the distances rather than the distances-squared
	    // forcing it to take a sqrt before doing the comparison.  But sqrt can make unequal
	    // things equal, so the old code might select a different sign distance than the new.
	    //
	    if (!closeEnough(fabsf(new_distance), fabsf(old_distance))) {
	      fprintf(stdout, "%s:%d new_distance:%g not near old_distance:%g,  magnitude diff:%g\n", __FILE__, __LINE__, 
	      	new_distance, old_distance, fabsf(new_distance)-fabsf(old_distance));
	      exit(1);
	    }
	  }
	  
          // update distance map
          //
          if (trace) fprintf(stdout, "  update distance for i:%ld j:%ld j:%ld\n", i,j,k);

          if (new_distanceSquared >= 0) {
            bool const lockNeeded =
#ifdef HAVE_OPENMP
	      do_old_MRIDistance;
#else
    	      false;
#endif
     	    if (lockNeeded) updateDistanceEltFromSignArgAndSquareLockNeeded  (ijkDistanceElt, new_distanceSignArg, new_distanceSquared);
            else            updateDistanceEltFromSignArgAndSquareNoLockNeeded(
                                ijkDistanceElt, 
                                new_distanceSignArg, 
                                new_distanceSquared,
                                sqrtf(new_distanceSquared));    // start this as soon as possible - the unit is idle, and the result
                                                                // is discarded if not needed
	    
          } else {
    	    updateDistanceElt(ijkDistanceElt, do_new_loop3 ? new_distance : old_distance, 
#ifdef HAVE_OPENMP
	      do_old_MRIDistance
#else
    	      false
#endif
	      );
	  }
          
	  if (do_old_MRIDistance && do_new_MRIDistance) {
            if (new_distanceSquared >= 0.0f) new_distance = SIGN(new_distanceSignArg) * sqrtf(new_distanceSquared);
	    float distance = do_new_loop3 ? new_distance : old_distance;
	    updateDistanceElt(&MRIFvox(mri_distance_nonconst, i, j, 0), distance, false);
	  }
	  
        } // k
      } // j
    } // i

    if (kToDo != kToDoBuffer) freeAndNULL(kToDo);
    if (kToZMap != kToZMapBuffer) freeAndNULL(kToZMap);
    if (jToYMap != jToYMapBuffer) freeAndNULL(jToYMap);
    
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  free(buffer);

  if (do_new_MRIDistance) {
    ROMP_PF_begin
    long i;
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(shown_reproducible) 
#endif
    for (i = 0; i < mri_distance->width; i++) {
      ROMP_PFLB_begin
      long j,k;
      for (j = 0; j < mri_distance->height; j++) {
        for (k = 0; k < mri_distance->depth; k++) {
	
	  float distance = NPY;
          int tid;
          for (tid = 0; tid < maxThreads; tid++) {
	    PerThreadMRIDistance* ptd = perThreadMRIDistances[tid];
	    if (ptd) updateDistanceElt(&distance, *perThreadMRIDistanceElt(ptd, i,j,k), false);
          }
	  
	  if (do_old_MRIDistance) {
	    float old_distance = MRIFvox(mri_distance_nonconst, i, j, k);
	    if (old_distance != distance) {
	      fprintf(stdout, "%s:%d diff distances at i:%ld j:%ld k:%ld old:%g new:%g\n", __FILE__, __LINE__,
	      	i,j,k,old_distance,distance);
	      exit(1);
	    }
	  }
	  
	  MRIFvox(mri_distance_nonconst, i, j, k) = distance;		// MODIFIER THAT WASN'T DETECTED
        }
      }
      ROMP_PFLB_end
    }
    ROMP_SCOPE_end
    
  }
  
  int tid;
  for (tid = 0; tid < maxThreads; tid++) {
    freePerThreadVertexPseudoNormalCache(&perThreadVertexPseudoNormalCache[tid]);
    freePerThreadMRIDistance            (&perThreadMRIDistances           [tid]);
  }

  freeAndNULL(sharedVertexPseudoNormalCache);

  //  TIMER_INTERVAL_END(taskExecution)

  free(fnos);
  } // int p;

  //  TIMER_INTERVAL_BEGIN(volumeLikelihood)

  /* compute the volumeLikelihood */
  /* init log values */
  MRI const * const mri_white = dp->mri_defect_white;
  MRI const * const mri_gray  = dp->mri_defect_gray;

    
  float white_ll = dp_defect->white_mean_ll;
  float gray_ll  = dp_defect->gray_mean_ll;
  int   nwhite   = 1;
  int   ngray    = 1;
  
  float int_g = dp_defect->gray_mean;
  float int_w = dp_defect->white_mean;

  float max_distance = 0.0;

  int k;
  for (k = 3; k < mri_distance->depth - 3; k++) {
    int j;
    for (j = 3; j < mri_distance->height - 3; j++) {
      int i;
      for (i = 3; i < mri_distance->width - 3; i++) {
      	// Note: The above order is how the _old code does it
	//       even though it is inefficient.
      
        if (!MRIvox(mri_defect, i, j, k)) {
          continue;
        }
        float val = MRIFvox(mri_distance, i, j, k);
        if (val == NPY) {
          continue;
        }
        if (fabs(val) > max_distance) {
          max_distance = fabs(val);
        }

        if (val > 0.5) {
          /* gray matter */
          gray_ll  += MRIFvox(mri_gray,   i, j, k);
          int_g    += MRIvox (mri_defect, i, j, k);
          ngray++;
        }
        if (val < -0.5) {
          white_ll += MRIFvox(mri_white,  i, j, k);
          int_w    += MRIvox (mri_defect, i, j, k);
          nwhite++; 
        }
      }
    }
  }
  
  if (nwhite) {
    white_ll /= nwhite;
    int_w    /= nwhite;
  }
  if (ngray) {
    gray_ll /= ngray;
    int_g   /= ngray;
  }

  {
    max_distance = MIN(2.0f, MAX(1.0f, max_distance));

    int p;
    for (p = 0; p < dp_defect->nvertices; p++) {
      if (dp_defect->status[p] == DISCARD_VERTEX) {
        continue;
      }
      int const vno = dp_defect->vertex_trans[dp_defect->vertices[p]];
      VERTEX * vertex_nonconst = &mris_nonconst->vertices[vno];

      int const i = iVOL(mri_defect, vertex_nonconst->fx);
      int const j = jVOL(mri_defect, vertex_nonconst->fy);
      int const k = kVOL(mri_defect, vertex_nonconst->fz);

      if ((i < 0) || (i >= mri_defect->width) || (j < 0) || (j >= mri_defect->height) || (k < 0) ||
          (k >= mri_defect->depth)) {
        continue;
      }

      // fprintf(WHICH_OUTPUT,"*%d-%d-%d");
      float val = MRIFvox(mri_distance, i, j, k);

      // FLO : max_distance=MIN(2.0f, max_distance) ?
      if (val == NPY) {
        val = 1.0f;
      }
      else {
        val = fabs(val) / max_distance;
      }

      val = MIN(1.0f, val);

      vertex_nonconst->curvbak *= val;			// MODIFIER
    }
  }

  dp_nonconst->tp.unmri_ll = (white_ll + gray_ll);	// MODIFIER

  //  TIMER_INTERVAL_END(volumeLikelihood)

  return (white_ll + gray_ll);
}





/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISripDefectiveFaces(MRI_SURFACE *mris)
{
  FACE *f;
  int fno, flist[MAX_INT_FACES], retained_i, i, j, nfaces, nripped, n;
  double dist, min_face_dist, max_dist, r;
  float dx, dy, dz, cx[MAX_INT_FACES], cy[MAX_INT_FACES], cz[MAX_INT_FACES];
  MHT *mht;

  MRISclearCurvature(mris);
  r = MRISaverageRadius(mris);
  MRISscaleBrain(mris, mris, 100.0 / r);  // TO BE CHECKED
  mht = MHTcreateFaceTable(mris);

/*
  first remove all faces with negative areas as they will be
  the 'underside' of a defect.
*/
#if 0
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    if (mris->faces[177616].v[2] > mris->nfaces)
    {
      DiagBreak() ;
    }
    f = &mris->faces[fno] ;
    if (f->ripflag)
    {
      continue ;
    }
    if (f->area < 0)
    {
      f->ripflag = 1 ;  /* part of a defect! */
      if (Gdiag & DIAG_SHOW)
      {
        fprintf(stdout, "ripping face %d\n", fno) ;
      }
      nripped++ ;
      for (i = 0 ; i < VERTICES_PER_FACE ; i++)
      {
        mris->vertices[f->v[i]].curv = -2.0f ;
      }
    }
  }
#endif
  for (nripped = fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    nfaces = mrisFindAllOverlappingFaces(mris, mht, fno, flist);
    if (nfaces > 1) /* overlapping faces - rip all but one of them */
    {
      /* retain the face that is furthest from any negative face */
      max_dist = 0.0f;
      for (i = 0; i < nfaces; i++) {
        mrisCalculateFaceCentroid(mris, flist[i], cx + i, cy + i, cz + i);
      }
      for (retained_i = -1, i = 0; i < nfaces; i++) {
        if (mris->faces[flist[i]].area < 0) {
          continue; /* don't ever retain a negative face */
        }

        /* find the distance to the closest negative face */
        for (min_face_dist = 1000000.0, j = 0; j < nfaces; j++) {
          if (mris->faces[flist[j]].area > 0) {
            continue; /* only consider distances to negative faces */
          }
          dx = cx[j] - cx[i];
          dy = cy[j] - cy[i];
          dz = cz[j] - cz[i];
          dist = (dx * dx + dy * dy + dz * dz);
          if (dist < min_face_dist) {
            min_face_dist = dist;
          }
        }

        /*
          if this face is more distant than any other face to a negative
          face (so far), tentatively mark it as the one to keep.
        */
        if (min_face_dist > max_dist) {
          max_dist = min_face_dist;
          retained_i = i;
        }
      }

      if (retained_i >= 0)
        for (i = 0; i < nfaces; i++) {
          VERTEX *v;
          FACE *f;

          f = &mris->faces[flist[i]];
          if (i == retained_i) {
            for (n = 0; n < VERTICES_PER_FACE; n++) {
              v = &mris->vertices[f->v[n]];
              if (FZERO(v->curv)) {
                v->curv = 1.0; /* good part of a defect! */
              }
            }
          }

          if (i == retained_i || f->ripflag) {
            continue;
          }
          f->ripflag = 1; /* part of a defect! */
          if (Gdiag & DIAG_SHOW) {
            fprintf(stdout, "ripping face %d\n", flist[i]);
          }
          nripped++;
          for (n = 0; n < VERTICES_PER_FACE; n++) {
            v = &mris->vertices[f->v[n]];
            if (v->curv >= 0.0) {
              v->curv = f->area > 0 ? -1.0 : -2.0;
            }
          }
        }
    }
  }

  MRISscaleBrain(mris, mris, r / 100.0);
  fprintf(stdout, "removing %d ripped faces from tessellation\n", nripped);
#if 0
  mrisRipVertices(mris) ;
#else
  {
#if 0
    int vno ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->curv = 0 ;
      for (fno = 0 ; fno < v->num ; fno++)
      {
        f = &mris->faces[v->f[fno]] ;
        if (f->ripflag == 1)
        {
          v->curv = -1 ;
          break ;
        }
        else if (f->ripflag == 2)
        {
          v->curv = 1.0 ;
        }
      }
    }
#endif
    MRISwriteCurvature(mris, "defects");
  }
#endif
  /*  MRISremoveRipped(mris) ;*/
  return (NO_ERROR);
}


static void destructComputeDefectContext(ComputeDefectContext* computeDefectContext) {

    if (computeDefectContext->mris_deferred_norms) {
        mrisurf_undeferSetFaceNorms(computeDefectContext->mris_deferred_norms);
        computeDefectContext->mris_deferred_norms = NULL;
    }
    
    if (computeDefectContext->realmTree) {
        removeActiveRealmTree(computeDefectContext->realmTree);
        freeRealmTree(&computeDefectContext->realmTree);
    }
}


// IntersectDefectEdgesContext are used to speed up intersectDefectEdges and intersectDefectConvexHullEdges
// by sharing computations across multiple calls
//

//#define COPE_WITH_VERTEX_MOVEMENT
typedef struct IntersectDefectEdgesContext_Entry {
  int vno, n;
#ifdef COPE_WITH_VERTEX_MOVEMENT
  int next;
  float cx,cy,cz,cx2,cy2,cz2;   // the values when created so movement can be detected
#endif
} IntersectDefectEdgesContext_Entry;

typedef struct IntersectDefectEdgesContext {
    int                                 entriesCapacity;
    int                                 entriesSize;
    IntersectDefectEdgesContext_Entry*  entries;
    GreatArcSet*                        greatArcSet;
    bool                                obsoleted;
    int                                 nvertices_seen;
#ifdef COPE_WITH_VERTEX_MOVEMENT
    int*                                vnoToFirstNPlus1;   // 0 is end of list
#else
    int*                                vnosHighestNSeen;   // Detect added vertices
#endif
} IntersectDefectEdgesContext;

static int intersectDefectEdges          (MRI_SURFACE *mris, DEFECT *defect, EDGE *e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2);

static int intersectDefectConvexHullEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2);


static void initIntersectDefectEdgesContext(IntersectDefectEdgesContext* ctx, MRI_SURFACE* mris) {
    ctx->entriesCapacity  = 0;
    ctx->entriesSize      = 0;
    ctx->entries          = NULL;
    ctx->greatArcSet      = NULL;
    ctx->obsoleted        = false;
    ctx->nvertices_seen   = mris->nvertices;
#ifdef COPE_WITH_VERTEX_MOVEMENT
    ctx->vnoToFirstNPlus1 = (int*)calloc(ctx->nvertices_seen,sizeof(int));
#else
    ctx->vnosHighestNSeen = (int*)calloc(ctx->nvertices_seen,sizeof(int));
#endif
}

static void finiIntersectDefectEdgesContext(IntersectDefectEdgesContext* ctx) {
#ifdef COPE_WITH_VERTEX_MOVEMENT
    freeAndNULL(ctx->vnoToFirstNPlus1); ctx->nvertices_seen = 0;
#else
    freeAndNULL(ctx->vnosHighestNSeen);
#endif
    freeGreatArcSet(&ctx->greatArcSet);
    free(ctx->entries); ctx->entries = NULL;
    ctx->entriesCapacity = ctx->entriesSize = 0;
}

static void obsoleteIntersectDefectEdgesContext(IntersectDefectEdgesContext* ctx) {
    ctx->obsoleted = true;
}


DEFECT_LIST *mrisSegmentDefects(MRI_SURFACE *mris, int mark_ambiguous, int mark_segmented)
{
  DEFECT_LIST *dl;
  int vno, nadded;
  VERTEX *v;
  DEFECT *defect;

  dl = (DEFECT_LIST *)calloc(1, sizeof(DEFECT_LIST));
  if (!dl) ErrorExit(ERROR_NO_MEMORY, "MRISsegmentDefects: could allocate defect list");

  /* uses fixedval to mark border vertices */
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].fixedval = 0;
  }

  for (nadded = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked != mark_ambiguous) {
      continue;
    }

    defect = &dl->defects[dl->ndefects];
    defect->defect_number = dl->ndefects++;

    /* segment defect #defect->defect_number */
    nadded += mrisSegmentDefect(mris, vno, defect, mark_ambiguous, mark_segmented);

    /* update the defect so it becomes simply connected */
    mrisSimplyConnectedDefect(mris, defect, mark_ambiguous, mark_segmented);
  }
  if (nadded) fprintf(stderr, "   total of %d vertices have been added to the surface\n", nadded);

  return (dl);
}

static int mrisMarkAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag);

static DEFECT_LIST *mrisRemoveOverlappingDefects(MRIS *mris, DEFECT_LIST *dl)
{
  int i, n, found, removed;
  DEFECT *defect, *removed_defect;
  DEFECT_LIST *new_dl;

  fprintf(stderr, "analyzing defects for overlaps...\n");

  found = 1;
  while (found) {
    found = 0;

    /* initiliaze flags */
    for (n = 0; n < mris->nvertices; n++) {
      mris->vertices[n].ripflag = 0;
      mris->vertices[n].marked = 0;
      mris->vertices[n].fixedval = 0;
      mris->vertices[n].undefval = 0;
      mris->vertices[n].old_undefval = 0;
    }

    for (i = 0; i < dl->ndefects; i++) {
      if (found) {
        break;
      }
      defect = &dl->defects[i];
      for (n = 0; n < defect->nvertices; n++) {
        if (mris->vertices[defect->vertices[n]].marked) {
          /* defect i overlaps defect (marked-1) */
          /* by construction defect (marked-1) is inside defect i */
          /* remove defect (marked-1) */

          fprintf(WHICH_OUTPUT,
                  "  -removing defect %d (overlapping defect %d)\n",
                  mris->vertices[defect->vertices[n]].marked - 1,
                  i);
          // fprintf(WHICH_OUTPUT,"defect %d has %d vertices - ",i,defect->nvertices);
          // fprintf(WHICH_OUTPUT,"defect %d has %d vertices :
          // %d\n",mris->vertices[defect->vertices[n]].marked-1,dl->defects[mris->vertices[defect->vertices[n]].marked-1].nvertices,defect->vertices[n]);

          removed = mris->vertices[defect->vertices[n]].marked - 1;
          removed_defect = &dl->defects[removed];

          /* free inside vertices */
          free(removed_defect->vertices);
          free(removed_defect->status);
          removed_defect->vertices = NULL;
          removed_defect->status = NULL;
          removed_defect->nvertices = 0;
          /* free border */
          free(removed_defect->border);
          removed_defect->border = NULL;
          removed_defect->nborder = 0;

          // newflo
          // free edge list!!!

          /* clean defect from empty spaces */
          new_dl = (DEFECT_LIST *)calloc(1, sizeof(DEFECT_LIST));
          for (i = 0; i < dl->ndefects; i++) {
            if (removed) {
              continue;
            }
            defect = &dl->defects[i];
            if (defect->nvertices) {
              /* this defect is not empty */
              memmove(&new_dl->defects[new_dl->ndefects++], defect, sizeof(DEFECT));
            }
          }
          free(dl);
          dl = new_dl;

          found = 1;
          break;
        }
        mris->vertices[defect->vertices[n]].marked = i + 1;
      }
    }
  }

  return dl;
}

static int mrisRipAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int ripflag);
static int mrisRipDefect(MRI_SURFACE *mris, DEFECT *defect, int ripflag);

void MRISidentifyDefects(MRIS *mris, TOPOFIX_PARMS *parms)
{
  int fno;
  FACE_DEFECT_LIST *fdl;
  DEFECT_LIST *dl;

  /* using CANONICAL_VERTICES */
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);

  /* marking intersecting edges */
  fdl = MRISmarkAmbiguousVertices(mris, MARK_AMBIGUOUS);

  fprintf(WHICH_OUTPUT, "segmenting defects...\n");

  dl = mrisSegmentDefects(mris, MARK_AMBIGUOUS, MARK_SEGMENTED);

  dl = mrisRemoveOverlappingDefects(mris, dl);

  MRISclearMarks(mris);

  /* free structures */
  for (fno = 0; fno < mris->nfaces; fno++)
    if (fdl->nfaces[fno] > 0) {
      free(fdl->faces[fno]);
    }
  free(fdl->faces);
  free(fdl->nfaces);
  free(fdl);

  // save structure into parms
  parms->defect_list = (void *)dl;
}

static int mrisComputeJointGrayWhiteBorderDistributions(MRI_SURFACE *mris, MRI *mri, MRI *mri_gray_white, MRI *mri_wm)
{
  int vno, x, y;
  VERTEX *v;
  float norm;
  double nx, ny, nz, xv, yv, zv, xw, yw, zw, white_val, gray_val;

  MRIScomputeMetricProperties(mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked || v->ripflag) /* part of a defect - ignore it */
    {
      continue;
    }

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    xw = v->x;
    yw = v->y;
    zw = v->z;

// MRIworldToVoxel(mri, xw+.5*nx, yw+.5*ny, zw+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(xw + .5 * nx, yw + .5 * ny, zw + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, xw + .5 * nx, yw + .5 * ny, zw + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolumeType(mri, xv, yv, zv, &gray_val, SAMPLE_NEAREST);
// MRIworldToVoxel(mri, xw-.5*nx, yw-.5*ny, zw-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(xw - .5 * nx, yw - .5 * ny, zw - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, xw - .5 * nx, yw - .5 * ny, zw - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolumeType(mri, xv, yv, zv, &white_val, SAMPLE_NEAREST);

#if 0
    //to be checked
#if 0
    if (gray_val >= MIN_WM_VAL &&
        white_val >= MIN_WM_VAL)  /* white on both sides */
    {
      continue ;
    }
#endif

    MRIFvox(mri_gray_white, nint(white_val), nint(gray_val), 0) += 1.0f ;
#else
    // set the value
    if (nint(v->val2) < 0 || nint(v->val2bak < 0) || nint(v->val2) >= mri_gray_white->width ||
        nint(v->val2bak) >= mri_gray_white->height)
      ErrorExit(ERROR_UNSUPPORTED, "gray/white vals out of [0 255] range (%d, %d)\n", nint(v->val2), nint(v->val));
    MRIFvox(mri_gray_white, nint(v->val2), nint(v->val2bak), 0) += 1.0f;

//      if ((nint(nint(v->val2)) == 110) && (nint(v->val2bak) == 110)) DiagBreak() ;
#endif
  }

  for (x = 0; x < 256; x++)
    for (y = 0; y < 256; y++) {
      if (FZERO(MRIFvox(mri_gray_white, x, y, 0))) {
        MRIFvox(mri_gray_white, x, y, 0) = 0.1;
      }
    }
  for (norm = 0.0, x = 0; x < 256; x++)
    for (y = 0; y < 256; y++) {
      norm += MRIFvox(mri_gray_white, x, y, 0);
    }

  for (x = 0; x < 256; x++)
    for (y = 0; y < 256; y++) {
      MRIFvox(mri_gray_white, x, y, 0) = MRIFvox(mri_gray_white, x, y, 0) / norm;
    }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_gray_white, "gw.mgh");
  }
  return (NO_ERROR);
}

static int mrisFindGrayWhiteBorderMean(MRI_SURFACE *mris, MRI *mri)
{
  double x, y, z, xv, yv, zv, gray_val, white_val, nx, ny, nz;
  int vno;
  VERTEX *v;

  MRIScomputeNormals(mris);
  /*  MRISsmoothSurfaceNormals(mris, 10) ;*/
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->marked) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    x = v->x;
    y = v->y;
    z = v->z;

// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &gray_val);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
    MRIsampleVolume(mri, xv, yv, zv, &white_val);
    v->val2 = white_val;
    v->val2bak = gray_val;
    v->val = (white_val + gray_val) / 2;
  }
#if 0
  MRISaverageVals(mris, 10) ;
  MRISaverageVal2s(mris, 10) ;
  MRISaverageVal2baks(mris, 10) ;
#else
  MRISmedianFilterVals(mris, 2);
  MRISmedianFilterVal2s(mris, 2);
  MRISmedianFilterVal2baks(mris, 2);
#endif
  return (NO_ERROR);
}



void MRISinitTopoFixParameters(MRIS *mris, TOPOFIX_PARMS *parms)
{
#if MATRIX_ALLOCATION
  /* allocation of the transform matrix */
  VoxelFromSRASmatrix = GetSurfaceRASToVoxelMatrix(parms->mri);
  parms->transformation_matrix = VoxelFromSRASmatrix;
#endif

  //    MRISsmoothSurfaceNormals(mris,10);

  /* v->val = border, v->val2 = white, v->val2bak = gray */
  mrisRipAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 1);
  mrisFindGrayWhiteBorderMean(mris, parms->mri);

  // find if there is a contrast inversion
  if (parms->contrast == -2) {
    parms->contrast = detectContrast(mris);
  }
  mrisRipAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 0);

  // computing curvature statistics
  MRISsetNeighborhoodSizeAndDist(mris, 2);
  parms->h_k1 = HISTOalloc(100);
  parms->h_k2 = HISTOalloc(100);
  parms->mri_k1_k2 = MRIalloc(100, 100, 1, MRI_FLOAT);
  parms->h_dot = HISTOalloc(100);
  mrisComputePrincipalCurvatureDistributions(mris, parms->h_k1, parms->h_k2, parms->mri_k1_k2);
  mrisComputeNormalDotDistribution(mris, parms->h_dot);
  MRISsetNeighborhoodSizeAndDist(mris, 1);

  // computing mri statistics
  MRIScomputeMetricProperties(mris);
  MRISsmoothSurfaceNormals(mris, 10);

  parms->h_gray = HISTOalloc(256);
  parms->h_white = HISTOalloc(256);
  parms->h_border = HISTOalloc(256);
  parms->h_grad = HISTOalloc(256);
  parms->mri_gray_white = MRIalloc(256, 256, 1, MRI_FLOAT);

  mrisMarkAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 1);
  mrisComputeJointGrayWhiteBorderDistributions(mris, parms->mri, parms->mri_gray_white, parms->mri_wm);

  /* compute statistics on original */
  mrisComputeSurfaceStatistics(
      mris, parms->mri, parms->h_k1, parms->h_k2, parms->mri_k1_k2, parms->mri_gray_white, parms->h_dot);
  mrisMarkAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 0);
}

MRI *mriInitDefectVolume(MRIS *mris, TOPOFIX_PARMS *parms)
{
  MRI *mri;
  VERTEX *v;

  int k, l, p, q;
  int width, height, depth;
  float xmin, xmax, ymin, ymax, zmin, zmax, scale;

  if (parms->volume_resolution == -1) {
    scale = VOLUME_SCALE;
  }
  else {
    scale = parms->volume_resolution;
  }

  /* find dimension of the volume */
  xmin = ymin = zmin = 1000.0;
  xmax = ymax = zmax = -1000.0;
  for (k = 0; k < mris->nvertices; k++) {
    v = &mris->vertices[k];
    if (v->marked2 != parms->defect_number) {
      continue;
    }
    if (v->origx > xmax) {
      xmax = v->origx;
    }
    if (v->origy > ymax) {
      ymax = v->origy;
    }
    if (v->origz > zmax) {
      zmax = v->origz;
    }
    if (v->origx < xmin) {
      xmin = v->origx;
    }
    if (v->origy < ymin) {
      ymin = v->origy;
    }
    if (v->origz < zmin) {
      zmin = v->origz;
    }
  }

  xmin -= 1.0f;
  ymin -= 1.0f;
  zmin -= 1.0f;
  xmax += 1.0f;
  ymax += 1.0f;
  zmax += 1.0f;

  /* allocate the volume */
  width = ceil(scale * (xmax - xmin));
  height = ceil(scale * (ymax - ymin));
  depth = ceil(scale * (zmax - zmin));
  mri = MRIalloc(width, height, depth, MRI_UCHAR);

  if (parms->verbose == VERBOSE_MODE_HIGH)
    fprintf(WHICH_OUTPUT, "      defect volume : %d by %d by %d (scale = %d)\n", width, height, depth, (int)scale);

  mri->xstart = xmin;
  mri->xsize = scale;

  mri->ystart = ymin;
  mri->ysize = scale;

  mri->zstart = zmin;
  mri->zsize = scale;

  for (q = 0; q < depth; q++)
    for (p = 0; p < height; p++)
      for (l = 0; l < width; l++) {
        MRIvox(mri, l, p, q) = 1;
      }

  return mri;
}


// used for fs_topo_fixer
void MRISsaveLocal(MRIS *mris, TOPOFIX_PARMS *parms, char *name)
{
  int static n_br = 0;
  char fname[512];
  int n;
  double x, y, z, xv, yv, zv;
  MRI *mri;

  mri = ((DP *)parms->dp)->mri_defect;

  n_br = parms->defect_number;

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  for (n = 0; n < mris->nvertices; n++) {
    x = xVOL(mri, mris->vertices[n].x);
    y = yVOL(mri, mris->vertices[n].y);
    z = zVOL(mri, mris->vertices[n].z);
    MRIvoxelToWorld(mri, x, y, z, &xv, &yv, &zv);
    mris->vertices[n].x = xv;
    mris->vertices[n].y = yv;
    mris->vertices[n].z = zv;
    mris->vertices[n].curv = mris->vertices[n].H;
  }
  sprintf(fname, "./def_%d.asc", n_br);
  MRISwrite(mris, fname);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  sprintf(fname, "./def_o_%d.asc", n_br);
  MRISwrite(mris, fname);

  mris = parms->mrip->mris_source;
  sprintf(fname, "./source_o_%d.asc", n_br);
  MRISwrite(mris, fname);
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  for (n = 0; n < mris->nvertices; n++) {
    x = xVOL(mri, mris->vertices[n].x);
    y = yVOL(mri, mris->vertices[n].y);
    z = zVOL(mri, mris->vertices[n].z);
    MRIvoxelToWorld(mri, x, y, z, &xv, &yv, &zv);
    mris->vertices[n].x = xv;
    mris->vertices[n].y = yv;
    mris->vertices[n].z = zv;
    mris->vertices[n].curv = mris->vertices[n].H;
  }
  sprintf(fname, "./source_%d.asc", n_br);
  MRISwrite(mris, fname);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);

  sprintf(fname, "./mri_sign_%d.mgz", n_br);
  MRIwrite(((DP *)parms->dp)->mri_defect_sign, fname);
  sprintf(fname, "./mri_sign_init_%d.mgz", n_br);
  MRIwrite(((DP *)parms->dp)->mri_defect_initial_sign, fname);
  // sprintf(fname,"def.curv_%d",n_br++);
  // MRISwriteCurvature(mris,fname);
}

static void defectVolumeWM(MRI *mri, MRI *mri_defect, MRI *mri_wm);

void mrisInitDefectMRIParameters(MRIS *mris, TOPOFIX_PARMS *parms)
{
  MRI *mri_defect, *mri_defect_white, *mri_defect_gray, *mri_defect_sign, *mri_defect_initial_sign, *mri_defect_wm;
  DP *dp;
  DEFECT *defect;

  dp = (DP *)parms->dp;
  defect = dp->defect;

  mri_defect = mri_defect_white = mri_defect_gray = mri_defect_sign = mri_defect_initial_sign = mri_defect_wm = NULL;

  if (!FZERO(parms->l_unmri)) {
    // initialization
    mri_defect = mriInitDefectVolume(mris, parms);
    // allocate volumes
    mri_defect_white = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_gray = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_sign = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_initial_sign = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_wm = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    // compute likelihood
    defectVolumeLikelihood(parms->mri,
                           mri_defect,
                           mri_defect_white,
                           mri_defect_gray,
                           parms->h_white,
                           parms->h_gray,
                           defect->white_mean,
                           defect->gray_mean,
                           1,
                           parms->contrast);
    defectVolumeWM(parms->mri_wm, mri_defect, mri_defect_wm);
  }
  dp->mri_defect = mri_defect;
  dp->mri_defect_white = mri_defect_white;
  dp->mri_defect_gray = mri_defect_gray;
  dp->mri_defect_sign = mri_defect_sign;
  dp->mri_defect_initial_sign = mri_defect_initial_sign;
  dp->mri_defect_wm = mri_defect_wm;
  dp->mri = parms->mri;
}

void MRISinitDefectParameters(MRIS *mris, TOPOFIX_PARMS *parms)
{
  int n, nvertices, nchull, defect_number;
  DEFECT *defect;
  DP *dp;

#if MATRIX_ALLOCATION
//  VoxelFromSRASmatrix=parms->transformation_matrix;
#endif

  defect_number = parms->defect_number;
  defect = (DEFECT *)calloc(1, sizeof(DEFECT));
  defect->defect_number = defect_number;

  // reset marks to zero
  MRISclearMarks(mris);
  nvertices = 0;
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked2 == defect_number) {
      nvertices++;
      mris->vertices[n].marked = 1;
    };
  if (nvertices == 0) {
    // empty defect!
    free(defect);
    return;
  }
  defect->nvertices = 0;
  defect->vertices = (int *)malloc(nvertices * sizeof(int));
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked) {
      defect->vertices[defect->nvertices++] = n;
    }
  // expand marks
  for (n = 0; n < 3; n++) {
    MRISexpandMarked(mris);
  }
  mrisMarkDefect(mris, defect, 2);
  nchull = 0;
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked == 1) {
      nchull++;
    }
  defect->chull = (int *)malloc(nchull * sizeof(int));
  defect->nchull = 0;
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked == 1) {
      defect->chull[defect->nchull++] = n;
    }

  // computing statistics
  MRISclearMarks(mris);

  mrisMarkAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 1);
  mrisComputeGrayWhiteBorderDistributions(
      mris, parms->mri, defect, parms->h_white, parms->h_gray, parms->h_border, parms->h_grad);

  computeDefectStatistics(parms->mri,
                          mris,
                          defect,
                          parms->h_white,
                          parms->h_gray,
                          parms->mri_gray_white,
                          parms->h_k1,
                          parms->h_k2,
                          parms->mri_k1_k2,
                          parms->verbose);
  mrisMarkAllDefects(mris, (DEFECT_LIST *)parms->defect_list, 0);

  // HISTOplot(parms->h_white,"w.plt") ;
  // HISTOplot(parms->h_gray,"g.plt") ;
  // MRIwrite(parms->mri_k1_k2,"./k1k2.mgz");
  // MRIwrite(parms->mri_gray_white,"gw.mgz");

  // initialize the Defect Patch associated with the current defect
  dp = (DP *)calloc(1, sizeof(DP));
  free(defect->chull);
  defect->chull = NULL;
  defect->nchull = 0;
  free(defect->vertices);
  defect->vertices = NULL;
  defect->nvertices = 0;
  dp->defect = defect;
  parms->dp = (void *)dp;

  // now initialize the volume associated with the defect
  mrisInitDefectMRIParameters(mris, parms);
}

void TOPOFIXfreeDP(TOPOFIX_PARMS *parms)
{
  DP *dp;
  dp = (DP *)parms->dp;
  if (dp == NULL) {
    return;
  }
  if (dp->defect) {
    free(dp->defect);
  }
  TPfree(&dp->tp);
  if (dp->mri_defect) {
    MRIfree(&dp->mri_defect);
  }
  if (dp->mri_defect_white) {
    MRIfree(&dp->mri_defect_white);
  }
  if (dp->mri_defect_gray) {
    MRIfree(&dp->mri_defect_gray);
  }
  if (dp->mri_defect_sign) {
    MRIfree(&dp->mri_defect_sign);
  }
  if (dp->mri_defect_initial_sign) {
    MRIfree(&dp->mri_defect_initial_sign);
  }
  if (dp->mri_defect_wm) {
    MRIfree(&dp->mri_defect_wm);
  }

  free(dp);

  parms->dp = NULL;
}

void MRISmarkPatchVertices(MRIS *mris, TOPOFIX_PARMS *parms, int ninitvertices)
{
  int n;
  DP *dp;
  TP *tp;

  dp = (DP *)parms->dp;
  tp = &dp->tp;

  MRISclearMarks(mris);

  for (n = ninitvertices; n < mris->nvertices; n++) {
    mris->vertices[n].marked = 1;
  }
}

void MRISmarkBorderVertices(MRIS *mris, TOPOFIX_PARMS *parms, int mark)
{
  int n;
  DP *dp;
  TP *tp;

  dp = (DP *)parms->dp;
  tp = &dp->tp;

  for (n = tp->ninside; n < tp->nvertices; n++) {
    mris->vertices[tp->vertices[n]].marked = mark;
  }
}

void MRISinitDefectPatch(MRIS *mris, TOPOFIX_PARMS *parms)
{
  int n;
  DP *dp;
  TP *tp;

  dp = (DP *)parms->dp;
  tp = &dp->tp;

  TPfree(tp);
  TPinit(tp);

  tp->nvertices = mris->nvertices;
  tp->vertices = (int *)malloc(tp->nvertices * sizeof(int));
  // first find the inside vertices
  MRISclearMarks(mris);
  tp->ninside = 0;
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[n];
    VERTEX                * const v  = &mris->vertices         [n];
    if (vt->vnum == vt->num) {
      v->marked = 1;
      tp->vertices[tp->ninside++] = n;
    }
  }
  // then the remaining ones
  tp->nvertices = tp->ninside;
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX *v = &mris->vertices[n];
    if (v->marked == 1) {
      continue;
    }
    tp->vertices[tp->nvertices++] = n;
  }
  MRISclearMarks(mris);
  // finally the faces
  tp->nfaces = mris->nfaces;
  tp->faces = (int *)malloc(tp->nfaces * sizeof(int));
  for (n = 0; n < tp->nfaces; n++) {
    tp->faces[n] = n;
  }
  tp->nedges = 0;
}

void MRIScomputeInitialFitness(MRIS *mris, TOPOFIX_PARMS *parms)
{
  parms->mris_defect = mris;
  MRIScomputeFitness(mris, parms, 0);
  parms->initial_fitness = 0.0;
  // save initial sign into mri_defect_initial_sign
  MRIcopy(((DP *)parms->dp)->mri_defect_sign, ((DP *)parms->dp)->mri_defect_initial_sign);
  //  fs_topo_fixer_test
  //  MRIwrite(((DP*)parms->dp)->mri_defect_initial_sign,"./mri_init_sign.mgz"); //
  //  MRIwrite(((DP*)parms->dp)->mri_defect_wm,"./mri_defect_wm.mgz");
}

static void MRISdefectMaximizeLikelihood(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha, int mode);

void MRISdefectMatch(MRIS *mris, TOPOFIX_PARMS *parms)
{
  // fs_topo_fixer_test
  // defectMatch(parms->mri,mris,(DP*)parms->dp,2,1);//parms->smooth, parms->match);
  MRISdefectMaximizeLikelihood(parms->mri, mris, (DP *)parms->dp, 40, 0.5, 1);
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
}

void MRISprintInfo(TOPOFIX_PARMS *parms)
{
  TP *tp;
  tp = &((DP *)parms->dp)->tp;
  TPprint(tp);
}

// ------------ Definition of the static functions -------------------- //
#define DO_NOT_USE_AREA 0

static void TPprint(TP *tp)
{
  double mri, curv;
#if DO_NOT_USE_AREA
  mri = (tp->face_ll + tp->vertex_ll) / 4.0;
  curv = (tp->qcurv_ll + tp->curv_ll) / 3.0;
#else
  mri = (tp->fll + tp->vll) / 4.0;
  curv = (tp->qcll + tp->cll) / 3.0;
#endif
  fprintf(WHICH_OUTPUT, "         mri =%3.3f   curv = %3.3f unmri = %3.3f\n", mri, curv, tp->unmri_ll);
  fprintf(WHICH_OUTPUT,
          "         ( f=%2.2f , v=%2.2f , c=%2.2f , q= %2.2f  ) \n",
          tp->face_ll,
          tp->vertex_ll,
          tp->curv_ll,
          tp->qcurv_ll);
  fprintf(WHICH_OUTPUT, "         ( f=%2.2f , v=%2.2f , c=%2.2f , q= %2.2f ) \n", tp->fll, tp->vll, tp->cll, tp->qcll);
}

static void TPinit(TP *tp)
{
  tp->vertices = NULL;
  tp->faces = NULL;
  tp->edges = NULL;
}

static void TPfree(TP *tp)
{
  if (tp->vertices) {
    free(tp->vertices);
  }
  if (tp->faces) {
    free(tp->faces);
  }
  if (tp->edges) {
    free(tp->edges);
  }
  TPinit(tp);
}

static SEGMENTATION *SEGMENTATIONalloc(int max_segments, int max_edges)
{
  int n;
  SEGMENTATION *seg;
  SEGMENT *segment;

  seg = (SEGMENTATION *)calloc(1, sizeof(SEGMENTATION));
  seg->segments = (SEGMENT *)calloc(max_segments, sizeof(SEGMENT));
  seg->nsegments = 0;
  seg->max_segments = max_segments;

  for (n = 0; n < max_segments; n++) {
    segment = &seg->segments[n];
    segment->edges = (int *)malloc(max_edges * sizeof(int));
    segment->max_edges = max_edges;
    segment->nedges = 0;

    segment->xedges = (int *)malloc(max_edges * sizeof(int));
    segment->max_xedges = max_edges;
    segment->nxedges = 0;
  }

  return seg;
}

static void SEGMENTATIONfree(SEGMENTATION **segmentation)
{
  int n;
  SEGMENTATION *seg;

  seg = *segmentation;
  *segmentation = NULL;

  for (n = 0; n < seg->max_segments; n++) {
    if (seg->segments[n].edges) {
      free(seg->segments[n].edges);
    }
    if (seg->segments[n].xedges) {
      free(seg->segments[n].xedges);
    }
  }
  if (seg->segments) {
    free(seg->segments);
  }
  free(seg);
}

static int findSegment(SEGMENTATION *segmentation)
{
  int n;
  SEGMENT *tmp;

  if (segmentation->nsegments == segmentation->max_segments) {
    segmentation->max_segments += 10;
    tmp = segmentation->segments;
    segmentation->segments = (SEGMENT *)calloc(segmentation->max_segments, sizeof(SEGMENT));
    memmove(segmentation->segments, tmp, segmentation->nsegments * sizeof(SEGMENT));
    free(tmp);
  }
  for (n = 0; n < segmentation->nsegments; n++)
    if (segmentation->segments[n].nedges == 0) {
      break;
    }

  return n;
}

static int compatibility(SEGMENTATION *segmentation, int segment_n, int edge_n)
{
  SEGMENT *segment;
  ES *es;
  int n, m;

  if (segment_n < 0) {
    return 0;
  }

  segment = &segmentation->segments[segment_n];

  for (n = 0; n < segment->nxedges; n++)
    if (segment->xedges[n] == edge_n) {
      return 0;
    }

  es = &segmentation->edges[edge_n];

  for (m = 0; m < es->nxedges; m++)
    for (n = 0; n < segment->nedges; n++)
      if (es->xedges[m] == segment->edges[n]) {
        return 0;
      }

  return 1;
}

/* add edge n into segment n */
static int addEdgeToSegment(SEGMENTATION *segmentation, int segment_n, int edge_n)
{
  SEGMENT *segment;
  int n, m, inside, val, *tmp, nedges;
  ES *es;

  segment = &segmentation->segments[segment_n];
  /* first add edge */
  if (segment->nedges == segment->max_edges) {
    tmp = segment->edges;
    segment->max_edges += 10;
    segment->edges = (int *)malloc(segment->max_edges * sizeof(int));
    memmove(segment->edges, tmp, segment->nedges * sizeof(int));
    free(tmp);
  }
  segment->edges[segment->nedges++] = edge_n;

  /* modify edge */
  es = &segmentation->edges[edge_n];
  es->segment = segment_n;

  /* add excluded edges */
  nedges = segment->nxedges + es->nxedges;
  if (nedges >= segment->max_xedges) {
    tmp = segment->xedges;
    segment->max_xedges = nedges + 10;
    segment->xedges = (int *)malloc(segment->max_xedges * sizeof(int));
    memmove(segment->xedges, tmp, segment->nxedges * sizeof(int));
    free(tmp);
  }
  for (n = 0; n < es->nxedges; n++) {
    val = es->xedges[n];
    for (inside = m = 0; m < segment->nxedges; m++)
      if (val == segment->xedges[m]) {
        inside = 1;
        break;
      }
    if (inside) {
      continue;
    }
    /* add edge */
    segment->xedges[segment->nxedges++] = val;
  }

  return NO_ERROR;
}

/* merge segment 2 into segment 1 */
static int mergeSegmentToSegment(SEGMENTATION *segmentation, int segment_1, int segment_2)
{
  int n, m, val, inside, nedges, *tmp;
  SEGMENT *segment1, *segment2;

  segment1 = &segmentation->segments[segment_1];
  segment2 = &segmentation->segments[segment_2];

  /* check compatibility of the two segments */
  for (n = 0; n < segment2->nxedges; n++)
    for (m = 0; m < segment1->nedges; m++)
      if (segment1->edges[m] == segment2->xedges[n]) {
        return NO_ERROR;
      }

  for (n = 0; n < segment1->nxedges; n++)
    for (m = 0; m < segment2->nedges; m++)
      if (segment2->edges[m] == segment1->xedges[n]) {
        return NO_ERROR;
      }

  nedges = segment1->nedges + segment2->nedges;
  /* reallocating if necessary */
  if (nedges >= segment1->max_edges) {
    tmp = segment1->edges;
    segment1->max_edges = nedges + 10;
    segment1->edges = (int *)malloc(segment1->max_edges * sizeof(int));
    memmove(segment1->edges, tmp, segment1->nedges * sizeof(int));
    free(tmp);
  }

  /* copy edges into segment 1 and updating edges */
  for (n = 0; n < segment2->nedges; n++) {
    val = segment2->edges[n];
    segment1->edges[n + segment1->nedges] = val;
    segmentation->edges[val].segment = segment_1;
  }
  segment1->nedges = nedges;
  segment2->nedges = 0;

  /* overlapping edges */
  nedges = segment1->nxedges + segment2->nxedges;
  if (nedges >= segment1->max_xedges) {
    tmp = segment1->xedges;
    segment1->max_xedges = nedges + 10;
    segment1->xedges = (int *)malloc(segment1->max_xedges * sizeof(int));
    memmove(segment1->xedges, tmp, segment1->nxedges * sizeof(int));
    free(tmp);
  }
  /* copy edges into segment1 */
  for (n = 0; n < segment2->nxedges; n++) {
    val = segment2->xedges[n];
    for (inside = m = 0; m < segment1->nxedges; m++)
      if (segment1->xedges[m] == val) {
        /* already inside segment1 */
        inside = 1;
        break;
      }
    if (inside) {
      continue;
    }
    segment1->xedges[segment1->nxedges++] = val;
  }
  segment2->nxedges = 0;
  segmentation->nsegments--;

  return NO_ERROR;
}

#define MAX_CLUSTERS 40

static void segmentEdge(MRIS *mris, SEGMENTATION *segmentation, int edge_n)
{
  int n, m, vno, stop;
  int clusters[MAX_CLUSTERS], nclusters, segment_n;
  EP *ep;

  /* check if the neighbors are already segmented */
  vno = segmentation->edges[edge_n].vno1;
  ep = (EP *)mris->vertices[vno].vp;
  for (nclusters = 0, n = 0; n < ep->nedges; n++) {
    if (ep->edges[n] == edge_n) {
      continue;
    }
    segment_n = segmentation->edges[ep->edges[n]].segment;
    for (stop = m = 0; m < nclusters; m++)
      if (segment_n == clusters[m]) {
        stop = 1;
        break;
      }
    if (stop) {
      continue;
    }
    if (compatibility(segmentation, segment_n, edge_n)) {
      clusters[nclusters++] = segment_n;
    }
  }
  vno = segmentation->edges[edge_n].vno2;
  ep = (EP *)mris->vertices[vno].vp;
  for (n = 0; n < ep->nedges; n++) {
    if (ep->edges[n] == edge_n) {
      continue;
    }
    segment_n = segmentation->edges[ep->edges[n]].segment;
    for (stop = m = 0; m < nclusters; m++)
      if (segment_n == clusters[m]) {
        stop = 1;
        break;
      }
    if (stop) {
      continue;
    }
    if (compatibility(segmentation, segment_n, edge_n)) {
      clusters[nclusters++] = segment_n;
    }
  }

  /* cluster this edge */
  switch (nclusters) {
    case 0: /* create new cluster */
      // fprintf(stderr,",0");
      n = findSegment(segmentation);
      addEdgeToSegment(segmentation, n, edge_n);
      segmentation->nsegments++;
      break;
    case 1: /* add edge to cluster */
      // fprintf(stderr,",1");
      addEdgeToSegment(segmentation, clusters[0], edge_n);
      break;
    default: /* add edge and merge clusters if possible :-) */
      //    fprintf(stderr,",2");
      /* add edge to smallest cluster */
      for (n = 1; n < nclusters; n++)
        if (segmentation->segments[clusters[1]].nedges < segmentation->segments[clusters[0]].nedges) {
          /* swap */
          m = clusters[0];
          clusters[0] = clusters[1];
          clusters[1] = m;
        }

      addEdgeToSegment(segmentation, clusters[0], edge_n);
      for (n = 1; n < nclusters; n++) {
        mergeSegmentToSegment(segmentation, clusters[0], clusters[n]);
      }
      break;
  }

  //  fprintf(stderr," %d",segmentation->nsegments);
}

static SEGMENTATION *segmentIntersectingEdges(MRIS *mris, DEFECT *defect, int *vertex_trans, ES *es, int *nedges)
{
  int n, m, vno1, vno2, ndiscarded, nes, vtn;
  int nsegments, *segments, *sizes, changed, tmp;
  EDGE e1, e2;
  SEGMENTATION *segmentation, *new_segmentation;

  nes = *nedges;

  if (!nes) {
    return NULL;
  }

  /* mark the inside vertices */
  for (n = 0; n < defect->nvertices; n++) {
    vtn = vertex_trans[defect->vertices[n]];
    if (vtn < 0 || vtn >= mris->nvertices) {
      fprintf(stderr, "segmentIntersectingEdges: vertex < 0 : should not happen");
      continue;
    }
    mris->vertices[vtn].marked = 1;
  }
  for (n = 0; n < defect->nborder; n++) {
    mris->vertices[vertex_trans[defect->border[n]]].marked = 0;
  }

/* first consider the inside edges */
#if 0

  for ( ndiscarded = n = 0 ; n < nes ; n++ )
  {
    if (mris->vertices[es[n].vno1].marked==0)
    {
      if (n<nes-1)
      {
        memmove(&es[n], &es[n+1], (nes-n-1)*sizeof(ES)) ;
      }
      nes--;
      n--;
      ndiscarded++;
      continue;
    }
    if (mris->vertices[es[n].vno2].marked==0)
    {
      if (n<nes-1)
      {
        memmove(&es[n], &es[n+1], (nes-n-1)*sizeof(ES)) ;
      }
      nes--;
      n--;
      ndiscarded++;
      continue;
    }
  }
  fprintf(WHICH_OUTPUT,
          "%d out of %d edges were discarded \n",ndiscarded,nes+ndiscarded);
#else
  {
    ES *new_es;
    new_es = (ES *)malloc(nes * sizeof(ES));
    for (m = n = 0; n < nes; n++) {
      if (mris->vertices[es[n].vno1].marked == 0) {
        continue;
      }
      if (mris->vertices[es[n].vno2].marked == 0) {
        continue;
      }

      new_es[m].vno1 = es[n].vno1;
      new_es[m].vno2 = es[n].vno2;
      new_es[m].n = es[n].n;
      new_es[m].segment = es[n].segment;
      new_es[m].xedges = es[n].xedges;
      new_es[m++].nxedges = es[n].nxedges;
    }
    for (n = 0; n < nes; n++) {
      if (mris->vertices[es[n].vno1].marked && mris->vertices[es[n].vno2].marked) {
        continue;
      }

      new_es[m].vno1 = es[n].vno1;
      new_es[m].vno2 = es[n].vno2;
      new_es[m].n = es[n].n;
      new_es[m].segment = es[n].segment;
      new_es[m].xedges = es[n].xedges;
      new_es[m++].nxedges = es[n].nxedges;
    }
    memmove(es, new_es, nes * sizeof(ES));
    free(new_es);
  }
#endif

  /* unmark the inside vertices */
  for (n = 0; n < defect->nvertices; n++) {
    vtn = vertex_trans[defect->vertices[n]];
    if (vtn < 0 || vtn >= mris->nvertices) {
      fprintf(stderr, "segmentIntersectingEdges: vertex < 0: should not happen");
      continue;
    }
    mris->vertices[vtn].marked = 0;
  }
  for (n = 0; n < defect->nborder; n++) {
    mris->vertices[vertex_trans[defect->border[n]]].marked = 0;
  }

  //  fprintf(WHICH_OUTPUT,"update surface\n");

  /* update surface */
  for (n = 0; n < defect->nvertices; n++) {
    mris->vertices[defect->vertices[n]].vp = NULL;
  }
  for (n = 0; n < defect->nborder; n++) {
    mris->vertices[defect->border[n]].vp = NULL;
  }

  for (n = 0; n < nes; n++) {
    EP *ep;

    vno1 = es[n].vno1;

    if (mris->vertices[vno1].vp == NULL) {
      mris->vertices[vno1].vp = (void *)calloc(1, sizeof(EP));
    }
    ep = (EP *)mris->vertices[vno1].vp;
    if (ep->nedges == 0) {
      ep->edges = (int *)malloc((25 + mris->vertices_topology[vno1].vnum) * sizeof(int));
    }
    ep->edges[ep->nedges++] = n;

    vno2 = es[n].vno2;

    if (mris->vertices[vno2].vp == NULL) {
      mris->vertices[vno2].vp = (void *)calloc(1, sizeof(EP));
    }
    ep = (EP *)mris->vertices[vno2].vp;
    if (ep->nedges == 0) {
      ep->edges = (int *)malloc((25 + mris->vertices_topology[vno2].vnum) * sizeof(int));
    }
    ep->edges[ep->nedges++] = n;
  }

  // fprintf(WHICH_OUTPUT,"finding overlapping edges \n");

  /* finding overlapping edges */
  for (n = 0; n < nes; n++) {
    es[n].nxedges = 0;
    es[n].xedges = (int *)malloc((nes - 1) * sizeof(int));
    e1.vno1 = es[n].vno1;
    e1.vno2 = es[n].vno2;
    for (m = 0; m < nes; m++) {
      if (n == m) {
        continue;
      }
      e2.vno1 = es[m].vno1;
      e2.vno2 = es[m].vno2;
      if (edgesIntersect(mris, &e1, &e2)) {
        es[n].xedges[es[n].nxedges++] = m;
      }
    }
#if 0
    fprintf(stdout,"%d and %d:",n,es[n].nxedges);
    if (es[n].nxedges)
      for ( m = 0  ; m <es[n].nxedges ; m++)
      {
        fprintf(stdout," %d",es[n].xedges[m]);
      }
    fprintf(stdout,"\n");
#endif
  }

  /* allocating structure */
  segmentation = SEGMENTATIONalloc(MAX_SEGMENTS, MAX_SEGMENT_EDGES);
  segmentation->edges = es;
  segmentation->nedges = nes;

  *nedges = nes;

  /* ********************************************* */
  /*                 MAIN LOOP                     */

  /* segmenting edges */
  for (n = 0; n < nes; n++) {
    segmentEdge(mris, segmentation, n);
  }

  /* ********************************************* */

  /* sorting of the segments */
  nsegments = segmentation->nsegments;
  segments = (int *)malloc(nsegments * sizeof(int));
  sizes = (int *)malloc(nsegments * sizeof(int));
  for (nsegments = n = 0; n < segmentation->max_segments; n++)
    if (segmentation->segments[n].nedges > 0) {
      segments[nsegments] = n;
      sizes[nsegments++] = segmentation->segments[n].nedges;
    }

  changed = 1;
  while (changed) {
    changed = 0;
    for (n = 0; n < nsegments - 1; n++)
      if (sizes[n] < sizes[n + 1]) {
        tmp = sizes[n];
        sizes[n] = sizes[n + 1];
        sizes[n + 1] = tmp;
        tmp = segments[n];
        segments[n] = segments[n + 1];
        segments[n + 1] = tmp;
        changed = 1;
      }
  }

  free(sizes);

  new_segmentation = (SEGMENTATION *)calloc(1, sizeof(SEGMENTATION));
  new_segmentation->nsegments = nsegments;
  new_segmentation->max_segments = nsegments;
  new_segmentation->segments = (SEGMENT *)calloc(nsegments, sizeof(SEGMENT));
  new_segmentation->edges = segmentation->edges;
  new_segmentation->nedges = segmentation->nedges;
  new_segmentation->mris = segmentation->mris;

  for (n = 0; n < nsegments; n++) {
    SEGMENT *sd, *ss;
    sd = &new_segmentation->segments[n];
    ss = &segmentation->segments[segments[n]];

    memmove(sd, ss, sizeof(SEGMENT));

    ss->edges = NULL;
    ss->xedges = NULL;

    for (m = 0; m < sd->nedges; m++) {
      new_segmentation->edges[sd->edges[m]].segment = n;
    }
  }

  free(segments);

  SEGMENTATIONfree(&segmentation);
  segmentation = new_segmentation;

  /* remove segments with less than 5 edges -
     keep at least 3 segments , at max 10 */
  for (ndiscarded = 0, n = 3; n < segmentation->max_segments; n++) {
    if (segmentation->segments[n].nedges <= 0) {
      continue;
    }
    if ((n >= 10) || segmentation->segments[n].nedges < 5) {
      for (m = 0; m < segmentation->segments[n].nedges; m++) {
        segmentation->edges[segmentation->segments[n].edges[m]].segment = -1;
      }

      segmentation->segments[n].nedges = 0;
      segmentation->segments[n].nxedges = 0;

      segmentation->nsegments--;

      ndiscarded++;
    }
  }

  if (DIAG_VERBOSE_ON)
    fprintf(WHICH_OUTPUT,
            "Edge Clustering: %d segments were found (%d were discarded )\n",
            segmentation->nsegments,
            ndiscarded);

#if 0
  for ( n = 0 ; n < segmentation->max_segments ; n++)
    if (segmentation->segments[n].nedges>0)
      fprintf(WHICH_OUTPUT,
              "#%d (%d edges) -",n,segmentation->segments[n].nedges);
#endif
  //  fprintf(WHICH_OUTPUT,"\nclustering...\n");

  /* merging edges into the smallest neighboring component */
  changed = 1;
  while (changed) {
    int vno, nbh, nbh_size;
    EP *ep;
    changed = 0;
    for (n = 0; n < segmentation->nedges; n++) {
      if (segmentation->edges[n].segment >= 0) {
        continue;
      }
      /* find neighbors of this edge */
      nbh = -1;
      nbh_size = -1;
      vno = segmentation->edges[n].vno1;
      ep = (EP *)mris->vertices[vno].vp;
      for (m = 0; m < ep->nedges; m++)
        if (segmentation->edges[ep->edges[m]].segment >= 0) {
          if (nbh_size == -1) {
            nbh = segmentation->edges[ep->edges[m]].segment;
            nbh_size = segmentation->segments[nbh].nedges;
          }
          else if (segmentation->segments[segmentation->edges[ep->edges[m]].segment].nedges < nbh_size) {
            nbh = segmentation->edges[ep->edges[m]].segment;
            nbh_size = segmentation->segments[nbh].nedges;
          }
        }
      vno = segmentation->edges[n].vno2;
      ep = (EP *)mris->vertices[vno].vp;
      for (m = 0; m < ep->nedges; m++)
        if (segmentation->edges[ep->edges[m]].segment >= 0) {
          if (nbh_size == -1) {
            nbh = segmentation->edges[ep->edges[m]].segment;
            nbh_size = segmentation->segments[nbh].nedges;
          }
          else if (segmentation->segments[segmentation->edges[ep->edges[m]].segment].nedges < nbh_size) {
            nbh = segmentation->edges[ep->edges[m]].segment;
            nbh_size = segmentation->segments[nbh].nedges;
          }
        }
      if (nbh >= 0) {
        SEGMENT *seg;
        changed = 1;
        /* add this edge into the segment nbh */
        segmentation->edges[n].segment = nbh;
        seg = &segmentation->segments[nbh];
        if (seg->nedges == seg->max_edges) {
          int *tp;
          tp = seg->edges;
          seg->max_edges += 10;
          seg->edges = (int *)malloc(seg->max_edges * sizeof(int));
          memmove(seg->edges, tp, seg->nedges * sizeof(int));
          free(tp);
        }
        seg->edges[seg->nedges++] = n;
      }
    }
  }
  for (n = 0; n < segmentation->max_segments; n++)
    if (segmentation->segments[n].nedges > 0) {
      if (DIAG_VERBOSE_ON)
        fprintf(WHICH_OUTPUT, "                 cluster %d has %d edges\n", n, segmentation->segments[n].nedges);
    }

  for (n = 0; n < nes; n++)
    if (es[n].xedges) {
      free(es[n].xedges);
    }

  for (n = 0; n < nes; n++) {
    EP *ep;
    vno1 = es[n].vno1;

    if (mris->vertices[vno1].vp) {
      ep = (EP *)mris->vertices[vno1].vp;
      if (ep->edges) {
        free(ep->edges);
      }
      free(mris->vertices[vno1].vp);
      mris->vertices[vno1].vp = NULL;
    }

    vno2 = es[n].vno2;

    if (mris->vertices[vno2].vp) {
      ep = (EP *)mris->vertices[vno2].vp;
      if (ep->edges) {
        free(ep->edges);
      }
      free(mris->vertices[vno2].vp);
      mris->vertices[vno2].vp = NULL;
    }
  }

  return segmentation;
}

static void saveSegmentation(
    MRIS *mris, MRIS *mris_corrected, DEFECT *defect, int *vertex_trans, ES *es, int nes, char *fname)
{
  int n, val, v_rgb, vtn;
  int r, g, b, rgb[10];
  char name[500];

  MRISRGBToAnnot(0, 225, 225, rgb[0]);
  MRISRGBToAnnot(205, 62, 78, rgb[1]);
  MRISRGBToAnnot(120, 62, 78, rgb[2]);
  MRISRGBToAnnot(196, 58, 250, rgb[3]);
  MRISRGBToAnnot(0, 148, 0, rgb[4]);
  MRISRGBToAnnot(220, 248, 164, rgb[5]);
  MRISRGBToAnnot(230, 148, 34, rgb[6]);
  MRISRGBToAnnot(0, 118, 14, rgb[7]);
  MRISRGBToAnnot(12, 48, 255, rgb[8]);
  MRISRGBToAnnot(122, 186, 220, rgb[9]);

  for (n = 0; n < nes; n++) {
    val = es[n].segment;
    r = 100 + 20 * val;
    g = 150 - 20 * val;
    b = 200 - 20 * val;
    // v_rgb=(int)(((int)((r<<24)+(g<<16)+((b)<<8)))>>8);

    if (val >= 0) {
      // v_rgb=(int)((int)(rgb[val%5]<<8)>>8);
      v_rgb = rgb[val % 10];
      mris_corrected->vertices[es[n].vno1].annotation = v_rgb;
      mris_corrected->vertices[es[n].vno2].annotation = v_rgb;
    }
  }

  MRISclearAnnotations(mris);

  for (n = 0; n < defect->nvertices; n++) {
    vtn = vertex_trans[defect->vertices[n]];
    if (vtn < 0 || vtn >= mris_corrected->nvertices) {
      fprintf(stderr, "saveSegmentation: vertex < 0 : should not happen");
      continue;
    }
    v_rgb = mris_corrected->vertices[vtn].annotation;
    mris->vertices[defect->vertices[n]].annotation = (int)((int)(v_rgb << 8) >> 8);
  }
  for (n = 0; n < defect->nborder; n++) {
    v_rgb = mris_corrected->vertices[vertex_trans[defect->border[n]]].annotation;
    mris->vertices[defect->border[n]].annotation = (int)((int)(v_rgb << 8) >> 8);
  }
  sprintf(name, "%s/rh.annotation%d", fname, defect->defect_number);
  fprintf(WHICH_OUTPUT,
          "writting annotation file for edge clustering "
          "segmentation of defect %d\n",
          defect->defect_number);

  MRISwriteAnnotation(mris, name);
}

static void generateOrdering(DP *dp, SEGMENTATION *segmentation, int i)
{
  int n, m, val, r;
  int *ordering, *counter, nedges;
  int *seg_order, nseg;

  /* nothing to be done for the first segment */
  if (i == 0) {
    return;
  }

  if (segmentation == NULL) {
    mrisMutateDefectPatch(dp, dp->etable, MUTATION_PCT_INIT);
    return;
  }

  seg_order = (int *)malloc(sizeof(int) * segmentation->nsegments);
  nseg = segmentation->nsegments;
  for (n = 0; n < nseg; n++) {
    seg_order[n] = n;
  }

  r = i % segmentation->nsegments + 1;

  /* find the nth segment */
  for (n = 0; n < segmentation->max_segments; n++) {
    if (segmentation->segments[n].nedges > 0) {
      r--;
    }
    if (r == 0) {
      break;
    }
  }
  counter = (int *)calloc(dp->nedges, sizeof(int));
  ordering = (int *)calloc(dp->nedges, sizeof(int));

  /* use the nth segment */
  nedges = 0;
  for (m = 0; m < segmentation->segments[n].nedges; m++) {
    val = segmentation->edges[segmentation->segments[n].edges[m]].n;
    ordering[nedges++] = val;
    counter[val] = 1;
  }

  /* use all the other segments */
  /* first generate random ordering of the segments */
  for (m = 0; m < 11; m++)
    for (n = 0; n < nseg; n++) {
      fflush(stdout);  // nicknote: prevents segfault on Linux PowerPC
      // when -O2 optimization is used w/gcc 3.3.3

      r = nint(randomNumber(0.0, (double)nseg - 1));

      val = seg_order[n];
      seg_order[n] = seg_order[r];
      seg_order[r] = val;
    }

  for (n = 0; n < nseg; n++)
    for (m = 0; m < segmentation->segments[seg_order[n]].nedges; m++) {
      val = segmentation->edges[segmentation->segments[seg_order[n]].edges[m]].n;
      if (counter[val]) {
        continue;
      }
      ordering[nedges++] = val;
      counter[val] = 1;
    }

  for (m = 0; m < dp->nedges; m++) {
    val = dp->ordering[m];
    if (counter[val] == 0) {
      counter[val] = 1;
      ordering[nedges++] = val;
    }
  }

  free(dp->ordering);
  dp->ordering = ordering;
  free(seg_order);

  if (r != i + 1) {
    mrisMutateDefectPatch(dp, dp->etable, MUTATION_PCT_INIT);
  }
}

static MRIS *extractDefect(MRIS *mris, DEFECT *defect)
{
  int n, vno, nvertices, nfaces, *vtrans, *vertex_trans;
  VERTEX *vdst, *vsrc;
  FACE *fdst, *fsrc;
  MRIS *mris_small;

  vertex_trans = defect->vertex_trans;

  /* marking vertices */
  nvertices = 0;
  for (n = 0; n < defect->nvertices; n++) {
    if (defect->status[n] == DISCARD_VERTEX) {
      continue;
    }
    vno = vertex_trans[defect->vertices[n]];
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    mris->vertices[vno].marked = 110;
    nvertices++;
  }
  for (n = 0; n < defect->nchull; n++) {
    vno = vertex_trans[defect->chull[n]];
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    mris->vertices[vno].marked = 110;
    nvertices++;
  }

  nfaces = 0;
  for (n = 0; n < mris->nfaces; n++) {
    fsrc = &mris->faces[n];
    if (mris->vertices[fsrc->v[0]].marked == 110 && mris->vertices[fsrc->v[1]].marked == 110 &&
        mris->vertices[fsrc->v[2]].marked == 110) {
      nfaces++;
    }
  }

  mris_small = MRISalloc(nvertices, nfaces);
  mris_small->type = MRIS_TRIANGULAR_SURFACE;
  mris_small->useRealRAS = mris->useRealRAS;

  vtrans = (int *)malloc(mris->nvertices * sizeof(int));

  /* vertex positions */
  for (nvertices = n = 0; n < defect->nvertices + defect->nchull; n++) {
    if (n < defect->nvertices) {
      if (defect->status[n] == DISCARD_VERTEX) {
        continue;
      }
      vno = vertex_trans[defect->vertices[n]];
    }
    else {
      vno = vertex_trans[defect->chull[n - defect->nvertices]];
    }

    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }

    vdst = &mris_small->vertices[nvertices];
    vsrc = &mris->vertices[vno];
    vdst->x = vsrc->x;
    vdst->y = vsrc->y;
    vdst->z = vsrc->z;
    if (vsrc->old_undefval) {
      vdst->curv = 1;
    }
    vtrans[vno] = nvertices++;
  }

  /* now the faces */
  nfaces = 0;
  for (n = 0; n < mris->nfaces; n++) {
    fsrc = &mris->faces[n];
    if (mris->vertices[fsrc->v[0]].marked == 110 && mris->vertices[fsrc->v[1]].marked == 110 &&
        mris->vertices[fsrc->v[2]].marked == 110) {
      fdst = &mris_small->faces[nfaces];
      fsrc = &mris->faces[n];
      fdst->v[0] = vtrans[fsrc->v[0]];
      fdst->v[1] = vtrans[fsrc->v[1]];
      fdst->v[2] = vtrans[fsrc->v[2]];
      nfaces++;
    }
  }

  /* unmarking vertices */
  for (n = 0; n < defect->nvertices; n++) {
    if (defect->status[n] == DISCARD_VERTEX) {
      continue;
    }
    vno = vertex_trans[defect->vertices[n]];
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    mris->vertices[vno].marked = 0;
  }

  for (n = 0; n < defect->nchull; n++) {
    vno = vertex_trans[defect->chull[n]];
    if (vno < 0 || vno >= mris->nvertices) {
      continue;
    }
    mris->vertices[vno].marked = 0;
  }

  free(vtrans);

  mrisCheckVertexFaceTopology(mris_small);
  
  return mris_small;
}

#define EXTRACT_SMALL_SURFACE 0

static void savePatch(MRI *mri, MRIS *mris, MRIS *mris_corrected, DVS *dvs, DP *dp, char *fname, TOPOLOGY_PARMS *parms)
{
  int i;
  VERTEX *vsrc, *vdst;
  MRIS *mris_small;

  retessellateDefect(mris, mris_corrected, dvs, dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* detect the new set of faces */
  detectDefectFaces(mris_corrected, dp);

  /* orient the patch faces */
  orientDefectFaces(mris_corrected, dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  if (parms->verbose == VERBOSE_MODE_LOW)
    fprintf(WHICH_OUTPUT,
            "(%d , %d , %d ) - %d vertices were discarded \n",
            dp->tp.ninside,
            dp->tp.nedges,
            dp->tp.nfaces,
            dp->tp.ndiscarded);

  // before smoothing and after
  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);

  if (EXTRACT_SMALL_SURFACE) {
    /* extract 'small' surface */
    mris_small = extractDefect(mris_corrected, dp->defect);
    /* save surface */
    MRISwrite(mris_small, fname);
  }
  else
  /* save surface */
  {
    MRISwrite(mris_corrected, fname);
  }

  /* smooth original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, dp, parms->smooth, 0);

  strcat(fname, "s");

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);

  if (EXTRACT_SMALL_SURFACE) {
    int nvertices;
    /* transfer current coord to 'small' surface */
    for (nvertices = i = 0; i < dp->defect->nvertices + dp->defect->nchull; i++) {
      int vno;
      if (i < dp->defect->nvertices) {
        if (dp->defect->status[i] == DISCARD_VERTEX) {
          continue;
        }
        vno = dp->defect->vertex_trans[dp->defect->vertices[i]];
      }
      else
        vno = dp->defect->vertex_trans[dp->defect->chull[i - dp->defect->nvertices]];
      if (vno < 0 || vno >= mris_corrected->nvertices) {
        continue;
      }
      vdst = &mris_small->vertices[nvertices++];
      vsrc = &mris_corrected->vertices[vno];
      vdst->x = vsrc->x;
      vdst->y = vsrc->y;
      vdst->z = vsrc->z;
    }
    /* save surface */
    MRISwrite(mris_small, fname);
  }
  else
  /* save surface */
  {
    MRISwrite(mris_corrected, fname);
  }

  /* smooth original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, dp, 0, parms->match);

  strcat(fname, "m");

  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);

  if (EXTRACT_SMALL_SURFACE) {
    int nvertices;
    /* transfer current coord to 'small' surface */
    for (nvertices = i = 0; i < dp->defect->nvertices + dp->defect->nchull; i++) {
      int vno;
      if (i < dp->defect->nvertices) {
        if (dp->defect->status[i] == DISCARD_VERTEX) {
          continue;
        }
        vno = dp->defect->vertex_trans[dp->defect->vertices[i]];
      }
      else {
        vno = dp->defect->vertex_trans[dp->defect->chull[i - dp->defect->nvertices]];
      }
      if (vno < 0 || vno >= mris_corrected->nvertices) {
        continue;
      }
      vdst = &mris_small->vertices[nvertices++];
      vsrc = &mris_corrected->vertices[vno];
      vdst->x = vsrc->x;
      vdst->y = vsrc->y;
      vdst->z = vsrc->z;
    }
    /* save surface */
    MRISwrite(mris_small, fname);
  }
  else
  /* save surface */
  {
    MRISwrite(mris_corrected, fname);
  }

  if (parms->smooth == 3) {
    DEFECT *defect;
    int *vtrans, vno;
    VERTEX *v;
    defect = dp->defect;
    vtrans = defect->vertex_trans;
    /* write curv for eliminated vertices */
    MRISclearCurvature(mris_corrected);
    for (i = 0; i < defect->nvertices; i++) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      vno = vtrans[defect->vertices[i]];
      if (vno < 0 || vno >= mris_corrected->nvertices) {
        continue;
      }
      v = &mris_corrected->vertices[vno];
      if (v->old_undefval) {
        v->curv = 1;
      }
    }
    strcat(fname, "_c");

    if (EXTRACT_SMALL_SURFACE) {
      /* save surface */
      MRISwriteCurvature(mris_small, fname);
      /* free small surface */
      MRISfree(&mris_small);
    }
    else
    /* save curvature  */
    {
      MRISwriteCurvature(mris_corrected, fname);
    }

    MRISclearCurvature(mris_corrected);
  }

  /* restore the vertex state */
  mrisRestoreVertexState(mris_corrected, dvs);

  /* reset the edges to the unused state
     (unless they were in the original tessellation) */
  for (i = 0; i < dp->nedges; i++) {
    if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION) {
      dp->etable->edges[i].used = NOT_USED;
    }
    if (dp->etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
      dp->etable->edges[i].used = USED_IN_ORIGINAL_TESSELLATION;
    }
  }

  /* free vertices,edges,faces tables */
  TPfree(&dp->tp);
}

static void updateVertexStatistics(
    MRIS *mris, MRIS *mris_corrected, DVS *dvs, RP *rp, DP *dp, int *vertex_trans, float fitness)
{
  DEFECT *defect;
  EDGE_TABLE *etable;
  int i, nedges;
  float total_vertex_fitness = 0.f, new_fitness;
  VERTEX *v;
  static int first_time = 1;

  // TO UPDATE TO BE CHECKED
  if (first_time) {
    first_time = 0;
  };

  fitness = 1.0f;  // to be updated ...

  nedges = dp->nedges;
  etable = dp->etable;
  defect = dp->defect;

  /* first mark the used vertices */
  for (i = 0; i < nedges; i++)
    if (etable->edges[i].used == USED_IN_NEW_TESSELLATION || etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
      mris_corrected->vertices[etable->edges[i].vno1].marked = FINAL_VERTEX;
      mris_corrected->vertices[etable->edges[i].vno2].marked = FINAL_VERTEX;
    }

  /* don't need border vertices */
  for (i = 0; i < defect->nborder; i++) {
    mris_corrected->vertices[vertex_trans[defect->border[i]]].marked = 0;
  }

  /* then compute the total fitness of these used vertices */
  total_vertex_fitness = 0.0f;
  for (i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    v = &mris_corrected->vertices[vertex_trans[defect->vertices[i]]];
    if (v->marked == FINAL_VERTEX) {
      total_vertex_fitness += v->curvbak * fitness;
    }
  }

  if (FZERO(total_vertex_fitness)) {
    total_vertex_fitness = 1.0f;
  }

  total_vertex_fitness /= 100.0f;
  total_vertex_fitness = 1.0f;  // TO BE CHECKED

  /* finally update statistics and reset marks to zero */
  for (i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    v = &mris_corrected->vertices[vertex_trans[defect->vertices[i]]];
    if (v->marked == FINAL_VERTEX) {
      new_fitness = (v->curvbak * fitness / total_vertex_fitness) + (float)rp->nused[i] * rp->vertex_fitness[i];
      rp->vertex_fitness[i] = new_fitness / ((float)rp->nused[i] + 1.0f);
      rp->nused[i]++;
    }
    v->marked = 0;
  }
}

static int deleteWorstVertices(MRIS *mris, RP *rp, DEFECT *defect, int *vertex_trans, float fraction, int count)
{
  int i, nvoxels, niters, init, changed;
  float max;
  int max_i;
  static float threshold = 4.0f;
  nvoxels = 0;

  if (count <= 0) {
    fprintf(WHICH_OUTPUT, "error: count (%d) <= 0 \n", count);
    count = 1;
  }

  /* first kill the non-used vertices */
  for (niters = 0, i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    niters++;
    if (rp->nused[i] == 0) {
      defect->status[i] = DISCARD_VERTEX;
      mris->vertices[vertex_trans[defect->vertices[i]]].ripflag = count;  // TO BE CHECKED
      nvoxels++;
    }
  }

  /* won't kill vertices if less than 10 */
  if (niters < 10) {
    return 0;
  }

  if (fraction > 0.1) {
    threshold /= 2.0f;
  }

  // kill at most 20% of the vertices
  if (niters)  // at least one voxel to be killed!
  {
    niters = MAX(1, (int)(fraction * niters));  // FLO
  }

  // if(nvoxels && fraction <= 0.01) //kill only worst voxel
  //    niters=max(1,niters);

  init = niters;

  changed = 1;
  while (niters && changed) {
    changed = 0;
    // find worst voxel
    max = 0;
    max_i = -1;
    for (i = 0; i < defect->nvertices; i++) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      if (rp->vertex_fitness[i] > max) {
        max = rp->vertex_fitness[i];
        max_i = i;
      }
    }

    if (max_i < threshold && (2 * niters < init)) {
      break;
    }

    if (max_i >= 0) {
      defect->status[max_i] = DISCARD_VERTEX;
      mris->vertices[vertex_trans[defect->vertices[max_i]]].ripflag = count;  // TO BE CHECKED
      nvoxels++;
      changed = 1;
    }
    niters--;
  }
#if 0
  /* reset vertex statistics */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    rp->vertex_fitness[i]=0;
    rp->nused[i]=0;
  }
#endif

  return nvoxels;
}

/* static void defectSurfaceToVolume
   (MRI* mri, float x, float y, float z, int *i, int *j , int *k){ */
/*   (*i)=iVOL(mri,x); */
/*   (*j)=jVOL(mri,x); */
/*   (*k)=kVOL(mri,x); */
/* } */

static MRI *mriDefectVolume(MRIS *mris, EDGE_TABLE *etable, TOPOLOGY_PARMS *parms)
{
  MRI *mri;
  VERTEX *v;
#if 0
  int i,j,p,numu,numv,u,w;
  float px0,px1,py0,py1,pz0,pz1,px,py,pz,d0,d1,d2,dmax,x0,y0,z0;
#endif
  int k, vno1, vno2, l, p, q, n, nmax;
  int width, height, depth;
  float len, dx, dy, dz, x, y, z, x1, y1, z1, x2, y2, z2, xmin, xmax, ymin, ymax, zmin, zmax, scale;

  if (parms->volume_resolution == -1) {
    scale = VOLUME_SCALE;
  }
  else {
    scale = parms->volume_resolution;
  }

  /* find dimension of the volume */
  xmin = ymin = zmin = 1000.0;
  xmax = ymax = zmax = -1000.0;
  for (k = 0; k < etable->nedges; k++) {
    v = &mris->vertices[etable->edges[k].vno1];
    if (v->origx > xmax) {
      xmax = v->origx;
    }
    if (v->origy > ymax) {
      ymax = v->origy;
    }
    if (v->origz > zmax) {
      zmax = v->origz;
    }
    if (v->origx < xmin) {
      xmin = v->origx;
    }
    if (v->origy < ymin) {
      ymin = v->origy;
    }
    if (v->origz < zmin) {
      zmin = v->origz;
    }
    v = &mris->vertices[etable->edges[k].vno2];
    if (v->origx > xmax) {
      xmax = v->origx;
    }
    if (v->origy > ymax) {
      ymax = v->origy;
    }
    if (v->origz > zmax) {
      zmax = v->origz;
    }
    if (v->origx < xmin) {
      xmin = v->origx;
    }
    if (v->origy < ymin) {
      ymin = v->origy;
    }
    if (v->origz < zmin) {
      zmin = v->origz;
    }
  }

  xmin -= 1.0f;
  ymin -= 1.0f;
  zmin -= 1.0f;
  xmax += 1.0f;
  ymax += 1.0f;
  zmax += 1.0f;

  /* allocate the volume */
  width = ceil(scale * (xmax - xmin));
  height = ceil(scale * (ymax - ymin));
  depth = ceil(scale * (zmax - zmin));
  mri = MRIalloc(width, height, depth, MRI_UCHAR);

  if (parms->verbose == VERBOSE_MODE_HIGH)
    fprintf(WHICH_OUTPUT, "      defect volume : %d by %d by %d (scale = %d)\n", width, height, depth, (int)scale);

  mri->xstart = xmin;
  mri->xsize = scale;

  mri->ystart = ymin;
  mri->ysize = scale;

  mri->zstart = zmin;
  mri->zsize = scale;

  /* find all the non-used voxels :
     go through all edges (should be sufficient) */
  for (k = 0; k < etable->nedges; k++) {
    vno1 = etable->edges[k].vno1;
    vno2 = etable->edges[k].vno2;

    /* starting point */
    x1 = mris->vertices[vno1].origx;
    y1 = mris->vertices[vno1].origy;
    z1 = mris->vertices[vno1].origz;

    /* end point */
    x2 = mris->vertices[vno2].origx;
    y2 = mris->vertices[vno2].origy;
    z2 = mris->vertices[vno2].origz;

    if (x2 < x1) {
      /* switch vertices */
      x = x1;
      x1 = x2;
      x2 = x;
      y = y1;
      y1 = y2;
      y2 = y;
      z = z1;
      z1 = z2;
      z2 = z;
    }

    /* length */
    len = scale * sqrt(SQR(x2 - x1) + SQR(y2 - y1) + SQR(z2 - z1));

    if (!FZERO(len)) {
      dx = (x2 - x1) / (2 * len);
      dy = (y2 - y1) / (2 * len);
      dz = (z2 - z1) / (2 * len);
      nmax = ceil(2 * len + 5); /* max number of points */
      for (x = x1, y = y1, z = z1, n = 0; (x < x2) && (n < nmax); x += dx, y += dy, z += dz, n++) {
        l = iVOL(mri, x);
        p = jVOL(mri, y);
        q = kVOL(mri, z);
        if ((l < 0) || (l >= mri->width) || (p < 0) || (p >= mri->height) || (q < 0) || (q >= mri->depth)) {
          continue;
        }
        MRIvox(mri, l, p, q) = 1;
      }
    }
    else {
      dx = dy = dz = 0.0f;
      l = iVOL(mri, x1);
      p = jVOL(mri, y1);
      q = kVOL(mri, z1);
      if ((l < 0) || (l >= mri->width) || (p < 0) || (p >= mri->height) || (q < 0) || (q >= mri->depth)) {
        continue;
      }
      MRIvox(mri, l, p, q) = 1;
    }

    /* last point */
    l = iVOL(mri, x2);
    p = jVOL(mri, y2);
    q = kVOL(mri, z2);
    if ((l < 0) || (l >= mri->width) || (p < 0) || (p >= mri->height) || (q < 0) || (q >= mri->depth)) {
      continue;
    }
    MRIvox(mri, l, p, q) = 1;
  }

#if 0
  /* then, mark the correct surface in the volume */
  for ( k = 0 ; k < mris->nfaces ; k++)
  {
    // calculate three vertices
    //fprintf(WHICH_OUTPUT,
    //"\r  %f            ",100.0*(float)k/(float)mris->nfaces);
    x0 =mris->vertices[mris->faces[k].v[0]].origx;
    y0 =mris->vertices[mris->faces[k].v[0]].origy;
    z0 =mris->vertices[mris->faces[k].v[0]].origz;
    x1 =mris->vertices[mris->faces[k].v[1]].origx;
    y1 =mris->vertices[mris->faces[k].v[1]].origy;
    z1 =mris->vertices[mris->faces[k].v[1]].origz;
    x2 =mris->vertices[mris->faces[k].v[2]].origx;
    y2 =mris->vertices[mris->faces[k].v[2]].origy;
    z2 =mris->vertices[mris->faces[k].v[2]].origz;

    i=iVOL(mri,x0);
    j=iVOL(mri,x1);
    p=iVOL(mri,x2);
    if ((i<0) && (j<0) && (p<0))
    {
      continue;
    }
    if ((i>=mri->width) && (j>=mri->width) && (p>=mri->width))
    {
      continue;
    }

    i=jVOL(mri,y0);
    j=jVOL(mri,y1);
    p=jVOL(mri,y2);
    if ((i<0) && (j<0) && (p<0))
    {
      continue;
    }
    if ((i>=mri->height) && (j>=mri->height) && (p>=mri->height))
    {
      continue;
    }

    i=kVOL(mri,z0);
    j=kVOL(mri,z1);
    p=kVOL(mri,z2);
    if ((i<0) && (j<0) && (p<0))
    {
      continue;
    }
    if ((i>=mri->depth) && (j>=mri->depth) && (p>=mri->depth))
    {
      continue;
    }

    // calculate the sides
    d0 = mri->xsize*sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = mri->ysize*sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = mri->zsize*sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));

    //fprintf(WHICH_OUTPUT,
    // "(%f,%f,%f)and (%f,%f,%f)",x0,x1,x2,d0,d1,mri->xsize);

    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;

    numu = (int)(ceil(2*d0));
    numv = (int)(ceil(2*dmax));

    for ( w = 0 ; w <= numv ; w++ )
    {
      px0 = x0 + (x2-x0)*w/numv;
      py0 = y0 + (y2-y0)*w/numv;
      pz0 = z0 + (z2-z0)*w/numv;
      px1 = x1 + (x2-x1)*w/numv;
      py1 = y1 + (y2-y1)*w/numv;
      pz1 = z1 + (z2-z1)*w/numv;

      for ( u = 0 ; u <= numu ; u++ )
      {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        i=iVOL(mri,px);
        j=jVOL(mri,py);
        p=kVOL(mri,pz);

        if ((i<0)||(i>=mri->width) ||
            (j<0)||(j>=mri->height)||
            (p<0)||(p>=mri->depth))
        {
          continue;
        }

        MRIvox(mri,i,j,p) = 0;
      }
    }
  }
#endif
  return mri;
}

static void defectVolumeWM(MRI *mri, MRI *mri_defect, MRI *mri_wm)
{
  int i, j, k;
  double x, y, z, xv, yv, zv, val;

  if (mri == NULL) {
    return;
  }

  for (k = 0; k < mri_defect->depth; k++)
    for (j = 0; j < mri_defect->height; j++)
      for (i = 0; i < mri_defect->width; i++) {
        /* corresponding surface coords */
        x = xSURF(mri_defect, i);
        y = ySURF(mri_defect, j);
        z = zSURF(mri_defect, k);

#if MATRIX_ALLOCATION
        mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif

        MRIsampleVolume(mri, xv, yv, zv, &val);
        MRIFvox(mri_wm, i, j, k) = val;
      }
}

static void defectVolumeLikelihood(MRI *mri,
                                   MRI *mri_defect,
                                   MRI *mri_white,
                                   MRI *mri_gray,
                                   HISTOGRAM *h_white,
                                   HISTOGRAM *h_gray,
                                   float white_mean,
                                   float gray_mean,
                                   int type,
                                   int contrast)
{
  int i, j, k, n;
  double x, y, z, xv, yv, zv, val, sigma;

  sigma = fabs(white_mean - gray_mean) / 2.0;

  for (n = 0, k = 0; k < mri_defect->depth; k++)
    for (j = 0; j < mri_defect->height; j++)
      for (i = 0; i < mri_defect->width; i++) {
        n++;
        /* corresponding surface coords */
        x = xSURF(mri_defect, i);
        y = ySURF(mri_defect, j);
        z = zSURF(mri_defect, k);

#if MATRIX_ALLOCATION
        mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
        MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif

        MRIsampleVolume(mri, xv, yv, zv, &val);
        MRIvox(mri_defect, i, j, k) = val;

        if (type == 0) {
          MRIFvox(mri_white, i, j, k) = log(h_white->counts[nint(MIN(val, white_mean))]);
          MRIFvox(mri_gray, i, j, k) = log(h_gray->counts[nint(MAX(val, gray_mean))]);
        }
        else {
          // topo_fixer

          // MRIFvox(mri_white,i,j,k) = log(1.0/(1.0+exp(-(val-mu)/sigma)));
          // MRIFvox(mri_gray,i,j,k) = log(1.0/(1.0+exp((val-mu)/sigma)));

          if (contrast == 1) {
            if (val > white_mean) {
              MRIFvox(mri_white, i, j, k) = 0.0;
            }
            else {
              MRIFvox(mri_white, i, j, k) = SQR((val - white_mean) / sigma);
            }
            if (val < gray_mean) {
              MRIFvox(mri_gray, i, j, k) = 0.0;
            }
            else {
              MRIFvox(mri_gray, i, j, k) = SQR((val - gray_mean) / sigma);
            }
          }
          else if (contrast == -1) {
            if (val < white_mean) {
              MRIFvox(mri_white, i, j, k) = 0.0;
            }
            else {
              MRIFvox(mri_white, i, j, k) = SQR((val - white_mean) / sigma);
            }
            if (val > gray_mean) {
              MRIFvox(mri_gray, i, j, k) = 0.0;
            }
            else {
              MRIFvox(mri_gray, i, j, k) = SQR((val - gray_mean) / sigma);
            }
          }
          else {
            MRIFvox(mri_white, i, j, k) = SQR((val - white_mean) / sigma);
            MRIFvox(mri_gray, i, j, k) = SQR((val - gray_mean) / sigma);
          }
        }
      }

  // HISTOplot(h_white,"w.plt");
  //  HISTOplot(h_gray,"g.plt");
  //  fprintf(stderr,"writing out volume");
  // MRIwrite(mri_white,"./mri_w.mgz");
  // MRIwrite(mri_gray,"./mri_g.mgz");
  //  fprintf(stderr,"done!\n");

  // fprintf(WHICH_OUTPUT,"%d voxels out of %d voxels\n",n,mri_defect->width*mri_defect->height*mri_defect->depth);
}


static void computeDefectStatistics(MRI *mri,
                                    MRIS *mris,
                                    DEFECT *defect,
                                    HISTOGRAM *h_white,
                                    HISTOGRAM *h_gray,
                                    MRI *mri_gw,
                                    HISTOGRAM *h_k1,
                                    HISTOGRAM *h_k2,
                                    MRI *mri_k1_k2,
                                    int verbose)
{
  float val, mean, var, mg, mw, vw, vg, total, max, white_val, gray_val, x, y, z, cx, cy, cz, ival, k1, k2, vk1, vk2,
      wv, gv;
  double xv, yv, zv, int_val;
  int i, j, n;
  HISTOGRAM *h;

  cx = cy = cz = ival = 0.0f;
  for (i = 0; i < defect->nvertices; i++) {
    x = mris->vertices[defect->vertices[i]].origx;
    y = mris->vertices[defect->vertices[i]].origy;
    z = mris->vertices[defect->vertices[i]].origz;

#if MATRIX_ALLOCATION
    mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
    MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif
    cx += (float)xv;
    cy += (float)yv;
    cz += (float)zv;
    MRIsampleVolume(mri, xv, yv, zv, &int_val);
    ival += (float)int_val;
  }
  cx /= (float)defect->nvertices;
  cy /= (float)defect->nvertices;
  cz /= (float)defect->nvertices;
  ival /= (float)defect->nvertices;

  if (verbose)
    fprintf(WHICH_OUTPUT,
            "   computing statistics for defect %d: %d vertices\n"
            "   location: [ (%d,%d,%d) - average intensity = %3.3f ]\n",
            defect->defect_number,
            defect->nvertices,
            (int)cx,
            (int)cy,
            (int)cz,
            ival);

  /* computing intensity statistics */

  /* white matter */
  h = h_white;

  // first find max value
  max = 0.0f;
  for (n = 0; n < h->nbins; n++) {
    if (h->counts[n] > max) {
      max = h->counts[n];
    }
  }

  // then, only consider bins with at least value > 1% of max
  total = 0;
  for (mean = 0, var = 0, n = 0; n < h->nbins; n++) {
    val = h->bins[n] - h->bin_size / 2.0;
    if (h->counts[n] < max / 100.0) {
      continue;
    }
    mean += val * h->counts[n];
    var += val * val * h->counts[n];
    total += h->counts[n];
  }

  mean /= total;
  var /= total;
  var = var - SQR(mean);
  if (var < 0) {
    var = 0;
  }

  white_val = mean;

  wv = var;

  /* white matter */
  h = h_gray;

  // first find max value
  max = 0.0f;
  for (n = 0; n < h->nbins; n++) {
    if (h->counts[n] > max) {
      max = h->counts[n];
    }
  }

  // then, only consider bins with at least value > 1% of max
  total = 0;
  for (mean = 0, var = 0, n = 0; n < h->nbins; n++) {
    val = h->bins[n] - h->bin_size / 2.0;
    if (h->counts[n] < max / 100.0) {
      continue;
    }
    mean += (val * h->counts[n]);
    var += (SQR(val) * h->counts[n]);
    total += h->counts[n];
  }

  mean /= total;
  var /= total;

  var = var - SQR(mean);

  if (var < 0) {
    var = 0;
  }

  gray_val = mean;
  gv = var;

  if (verbose)
    fprintf(WHICH_OUTPUT,
            "      -gray ( %2.2f , %2.2f )  -white ( %2.2f , %2.2f )\n",
            gray_val,
            sqrt(gv),
            white_val,
            sqrt(wv));

  for (mw = 0, mg = 0, vw = 0, vg = 0, i = 0; i < 256; i++)
    for (j = 0; j < 256; j++) {
      val = (float)i;
      mw += val * MRIFvox(mri_gw, i, j, 0);
      vw += SQR(val) * MRIFvox(mri_gw, i, j, 0);
      val = (float)j;
      mg += val * MRIFvox(mri_gw, i, j, 0);
      vg += SQR(val) * MRIFvox(mri_gw, i, j, 0);
    }

  vw -= SQR(mw);
  vg -= SQR(mg);

  if (verbose)
    fprintf(WHICH_OUTPUT, "      -gray ( %2.2f , %2.2f )  -white ( %2.2f , %2.2f ) \n", mg, sqrt(vg), mw, sqrt(vw));

  mw = white_val;
  mg = gray_val;

  defect->white_mean = mw;
  defect->gray_mean = mg;

  defect->white_mean_ll = log(h_white->counts[MAX(0, MIN(h_white->nbins - 1, nint(mw)))]);
  defect->gray_mean_ll = log(h_gray->counts[MAX(0, MIN(h_gray->nbins - 1, nint(mg)))]);

  if (verbose)
    fprintf(WHICH_OUTPUT,
            "      -intensity (%f [log = %f ]- %f [log = %f ])\n",
            mg,
            defect->gray_mean_ll,
            mw,
            defect->white_mean_ll);

  /* computing curvature statistics */

  /* principal curvature k1 */
  h = h_k1;

  // first find max value
  max = 0.0f;
  for (n = 0; n < h->nbins; n++) {
    if (h->counts[n] > max) {
      max = h->counts[n];
    }
  }

  // then, only consider bins with at least value > 1% of max
  total = 0;
  for (mean = 0, var = 0, n = 0; n < h->nbins; n++) {
    val = h->bins[n] - h->bin_size / 2.0;
    if (h->counts[n] < max / 100.0) {
      continue;
    }
    mean += val * h->counts[n];
    var += val * val * h->counts[n];
    total += h->counts[n];
  }

  mean /= total;
  var /= total;

  var = var - SQR(mean);
  if (var < 0) {
    var = 0;
  }
  vk1 = var;

  defect->k1_mean = mean;

  /* principal curvature k2 */
  h = h_k2;

  // first find max value
  max = 0.0f;
  for (n = 0; n < h->nbins; n++) {
    if (h->counts[n] > max) {
      max = h->counts[n];
    }
  }

  // then, only consider bins with at least value > 1% of max
  total = 0;
  for (mean = 0, var = 0, n = 0; n < h->nbins; n++) {
    val = h->bins[n] - h->bin_size / 2.0;
    if (h->counts[n] < max / 100.0) {
      continue;
    }
    mean += val * h->counts[n];
    var += val * val * h->counts[n];
    total += h->counts[n];
  }

  mean /= total;
  var /= total;

  var = var - SQR(mean);
  if (var < 0) {
    var = 0;
  }
  vk2 = var;

  defect->k2_mean = mean;

  if (verbose)
    fprintf(WHICH_OUTPUT,
            "      -curv (k1=%3.3f (%3.3f) , "
            "r1 = %3.3f | k2=%3.3f (%3.3f), r2 = %3.3f )\n",
            defect->k1_mean,
            sqrt(vk1),
            -1.0 / defect->k1_mean,
            defect->k2_mean,
            sqrt(vk2),
            -1.0 / defect->k2_mean);

  for (k1 = 0, vk1 = 0, k2 = 0, vk2 = 0, i = 0; i < mri_k1_k2->width; i++)
    for (j = 0; j < mri_k1_k2->height; j++) {
      val = mri_k1_k2->xstart + (float)i * mri_k1_k2->xsize;
      k1 += val * MRIFvox(mri_k1_k2, i, j, 0);
      vk1 += SQR(val) * MRIFvox(mri_k1_k2, i, j, 0);
      val = (float)mri_k1_k2->ystart + (float)j * mri_k1_k2->ysize;
      k2 += val * MRIFvox(mri_k1_k2, i, j, 0);
      vk2 += SQR(val) * MRIFvox(mri_k1_k2, i, j, 0);
    }
  vk1 -= SQR(k1);
  vk2 -= SQR(k2);

  if (verbose)
    fprintf(WHICH_OUTPUT,
            "      -curv (k1=%3.3f (%3.3f) , r1 = %3.3f "
            "| k2=%3.3f (%3.3f), r2 = %3.3f )\n",
            k1,
            sqrt(vk1),
            -1.0 / k1,
            k2,
            sqrt(vk2),
            -1.0 / k2);
}

/* call the right smoothing and matching functions */
static void defectMatch(MRI *mri, MRI_SURFACE *mris, DP *dp, int smooth, int match)
{
  if (smooth) {
    defectSmooth(mris, dp, 25, 0.1, smooth);
  }
  if (match) {
    defectMaximizeLikelihood(mri, mris, dp, 40, 0.5);
  }
}

/* perform a light smoothing of the defect vertices */
static void defectSmooth(MRI_SURFACE *mris, DP *dp, int niter, double alpha, int type)
{
  int i, n;

  float x, y, z;
  float r, F, E, rmin, rmax;
  float dx, dy, dz, sx, sy, sz, sd, sxn, syn, szn, sxt, syt, szt, nc, nx, ny, nz, f;
  int *vertices, nvertices, ninside;
  float mean, var;
  int changed, should_be_smoothed, nstrictlyinside;
  float percentage;

  // ninside=dp->tp.nvertices; //TO BE CHECKED
  ninside = dp->tp.ninside; /* smooth only inside vertices, not border ones */
  if (ninside == 0) {
    return;
  }

  nstrictlyinside = dp->tp.ninside;

  if (type == 0) {
    return;
  }

  switch (type) {
    case 1:
      while (niter--) {
        /* using the tmp vertices */
        for (i = 0; i < ninside; i++) {
          int const vno = dp->tp.vertices[i];
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
          VERTEX                * const v  = &mris->vertices         [vno];

          for (x = 0, y = 0, z = 0, n = 0; n < vt->vnum; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            x += vn->origx;
            y += vn->origy;
            z += vn->origz;
          }
          if (n) {
            x /= (float)n;
            y /= (float)n;
            z /= (float)n;
          }

          v->tx = v->origx + alpha * (x - v->origx);
          v->ty = v->origy + alpha * (y - v->origy);
          v->tz = v->origz + alpha * (z - v->origz);
        }

        for (i = 0; i < ninside; i++) {
          int const vno = dp->tp.vertices[i];
          VERTEX * const v = &mris->vertices[vno];

          MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
        }
      }
      break;
    case 2:

      rmin = -1 / dp->defect->k1_mean;
      rmax = -1 / dp->defect->k2_mean;

      E = (1 / rmin + 1 / rmax) / 2;
      F = 6 / (1 / rmin - 1 / rmax);

      while (niter--) {
        computeDefectFaceNormals(mris, dp, NULL);
        computeDefectVertexNormals(mris, dp);

        /* using the tmp vertices */
        for (i = 0; i < ninside; i++) {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
          VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
          x = v->origx;
          y = v->origy;
          z = v->origz;
          nx = v->nx;
          ny = v->ny;
          nz = v->nz;

          sx = sy = sz = sd = 0;
          n = 0;
          for (n = 0; n < vt->vnum; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];

            sx += dx = vn->origx - x;
            sy += dy = vn->origy - y;
            sz += dz = vn->origz - z;
            sd += sqrt(dx * dx + dy * dy + dz * dz);
            n++;  // ATH is this supposed to be incremented again?
          }
          // mean distance to the neighbors
          sx = sx / (float)n;
          sy = sy / (float)n;
          sz = sz / (float)n;
          sd = sd / (float)n;

          nc = sx * nx + sy * ny + sz * nz;

          // normal component of the mean distance vector
          sxn = nc * nx;
          syn = nc * ny;
          szn = nc * nz;
          // tangential component of the mean distance vector
          sxt = sx - sxn;
          syt = sy - syn;
          szt = sz - szn;

          r = (nc > 0) ? nc : -nc;
          r = SQR(sd) / (2 * r);
          f = (1 + tanh(F * (1 / r - E))) / 2;

          v->tx = v->origx + alpha * (sxt + f * sxn);
          v->ty = v->origy + alpha * (syt + f * syn);
          v->tz = v->origz + alpha * (szt + f * szn);
        }
        for (i = 0; i < ninside; i++) {
          int const vno = dp->tp.vertices[i];
          VERTEX * const v = &mris->vertices[vno];

          MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
        }
      }
      break;
    case 3: /* only smooth the high undefval vertices */
      if (nstrictlyinside == 0) {
        return;
      }

      for (i = 0; i < nstrictlyinside; i++) {
        VERTEX * const v = &mris->vertices[dp->tp.vertices[i]];
        v->old_undefval = 0;
      }

      /* first we need to find which ones need to be smoothed */
      changed = 1;
      percentage = 0.0f;
      while (changed) {
        changed = 0;
        for (i = 0; i < nstrictlyinside; i++) {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
          VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
          if (v->old_undefval) {
            continue; /* already processed */
          }
          /* check if neighboring values are lower or not */
          should_be_smoothed = 1;
          for (n = 0; n < vt->vnum; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            if ((!vn->old_undefval) && vn->undefval < v->undefval) {
              should_be_smoothed = 0;
              break;
            }
          }
          if (should_be_smoothed) {
            changed = 1;
            v->old_undefval = 1;
            percentage += 1.0f;
          }
        }
      }
      percentage = 100.0f * percentage / (float)nstrictlyinside;
      // fprintf(stderr,"%2.3f %% vertices are to be smoothed\n",percentage);

      /* then, we smooth these vertices */
      if (percentage > 50.0f) /* perform a light smoothing only */
      {
        return defectSmooth(mris, dp, niter, alpha, 2);
      }
      niter = 20;
      alpha = 0.5;
      while (niter--) {
        /* using the tmp vertices */
        for (i = 0; i < nstrictlyinside; i++) {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
          VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
          if (v->old_undefval == 0) {
            continue;
          }

          for (x = 0, y = 0, z = 0, n = 0; n < vt->vnum; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            x += vn->origx;
            y += vn->origy;
            z += vn->origz;
          }
          if (n) {
            x /= (float)n;
            y /= (float)n;
            z /= (float)n;
          }

          v->tx = v->origx + alpha * (x - v->origx);
          v->ty = v->origy + alpha * (y - v->origy);
          v->tz = v->origz + alpha * (z - v->origz);
        }

        for (i = 0; i < nstrictlyinside; i++) {
          int const vno = dp->tp.vertices[i];
          VERTEX * const v = &mris->vertices[vno];
          if (v->old_undefval == 0) {
            continue;
          }
          MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
        }
      }
      /* finally we apply a light smoothing of the whole surface */
      defectSmooth(mris, dp, 20, 0.1, 2);

      break;
    default:
      rmin = -1 / dp->defect->k1_mean;
      rmax = -1 / dp->defect->k2_mean;

      E = (1 / rmin + 1 / rmax) / 2;
      F = 6 / (1 / rmin - 1 / rmax);

      // detect high curvatures vertices
      computeDefectFaceNormals(mris, dp, NULL);
      computeDefectVertexNormals(mris, dp);

      mean = var = 0.0;
      for (i = 0; i < ninside; i++) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
        VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
        x = v->origx;
        y = v->origy;
        z = v->origz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;

        sx = sy = sz = sd = 0;
        n = 0;
        for (n = 0; n < vt->vnum; n++) {
          VERTEX const * const vn = &mris->vertices[vt->v[n]];

          sx += dx = vn->origx - x;
          sy += dy = vn->origy - y;
          sz += dz = vn->origz - z;
          sd += sqrt(dx * dx + dy * dy + dz * dz);
          n++;  // ATH is this supposed to be incremented again?
        }
        // mean distance to the neighbors
        sx = sx / (float)n;
        sy = sy / (float)n;
        sz = sz / (float)n;
        sd = sd / (float)n;

        nc = sx * nx + sy * ny + sz * nz;

        // curvature
        // r= (nc>0) ? nc : -nc;
        // r=SQR(sd)/(2*r);
        r = nc;
        r = (2 * r) / SQR(sd);

        mean += r;
        var += (r * r);

        v->curv = r;
      }
      mean = mean / dp->tp.nvertices;
      var = var / dp->tp.nvertices - mean * mean;
      if (var < 0) {
        var = 0.0;
      }
      // fprintf(WHICH_OUTPUT,"mean=%f, var=%f \n",mean,sqrt(var));
      vertices = (int *)malloc(dp->tp.nvertices * sizeof(int));
      nvertices = 0;
      for (i = 0; i < ninside; i++) {
        VERTEX const * const v = &mris->vertices[dp->tp.vertices[i]];
        if (v->curv < mean - 1 * sqrt(var) || v->curv > mean + 1 * sqrt(var)) {
          vertices[nvertices++] = i;
        }
      }

#if 0
    fprintf(WHICH_OUTPUT,"%d (out of %d) vertices with high curvature\n",
            nvertices,dp->tp.ninside);

    MRISclearAnnotations(mris);

    for (i = 0 ; i < nvertices; i++)
    {
      VERTEX * const v = &mris->vertices[dp->tp.vertices[vertices[i]]] ;
      v->annotation=100;
    }
#endif

#if 0
    {
      static int counter=0;
      char fname[100];
      sprintf(fname,"./rh.test_%d",counter);
      MRISwriteAnnotation(mris,fname);
      counter++;
    }
#endif
      // smooth these vertices
      while (niter--) {
        /* using the tmp vertices */
        for (i = 0; i < nvertices; i++) {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[vertices[i]]];
          VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[vertices[i]]];

          for (x = 0, y = 0, z = 0, n = 0; n < vt->vnum; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            x += vn->origx;
            y += vn->origy;
            z += vn->origz;
          }
          if (n) {
            x /= (float)n;
            y /= (float)n;
            z /= (float)n;
          }

          v->tx = v->origx + alpha * (x - v->origx);
          v->ty = v->origy + alpha * (y - v->origy);
          v->tz = v->origz + alpha * (z - v->origz);
        }

        for (i = 0; i < nvertices; i++) {
          int const vno = dp->tp.vertices[vertices[i]];
          VERTEX * const v = &mris->vertices[vno];
          MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
        }
      }
      free(vertices);

      break;
  }
  computeDefectFaceNormals  (mris, dp, NULL);
  computeDefectVertexNormals(mris, dp);
}

static void MRISdefectMaximizeLikelihood(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha, int mode)
{
  float wm, gm, mean;
  int i, n, nvertices, *vertices;
  double x, y, z, xm, ym, zm, nx, ny, nz, dx, dy, dz, g, NRG;
  double xv, yv, zv, white_val, gray_val, val;

  if (mode == 0) {
    return defectMaximizeLikelihood(mri, mris, dp, niter, alpha);
  }

  // matching marked vertices only !

  wm = dp->defect->white_mean;
  gm = dp->defect->gray_mean;

  mean = (wm + gm) / 2.0;

  // find marked vertices
  nvertices = 0;
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked) {
      nvertices++;
    }
  vertices = (int *)malloc(nvertices * sizeof(int));
  if (nvertices == 0) {
    return;
  }
  nvertices = 0;
  for (n = 0; n < mris->nvertices; n++)
    if (mris->vertices[n].marked) {
      vertices[nvertices++] = n;
    }

  while (niter--) {
    /* using the tmp vertices */
    for (NRG = 0, i = 0; i < nvertices; i++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vertices[i]];
      VERTEX                * const v  = &mris->vertices         [vertices[i]];
      x = v->origx;
      y = v->origy;
      z = v->origz;

      /* smoothness term */
      for (xm = 0, ym = 0, zm = 0, n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        xm += vn->origx;
        ym += vn->origy;
        zm += vn->origz;
      }
      if (n) {
        xm /= (double)n;
        ym /= (double)n;
        zm /= (double)n;
      }

      /* image gradient */
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &white_val);

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &gray_val);

      g = (white_val - gray_val) * (val - mean);

      if (fabs(g) > 0.1) {
        g = 0.1 * g / fabs(g);
      }

      dx = 0.5 * (xm - x) + g * nx;
      dy = 0.5 * (ym - y) + g * ny;
      dz = 0.5 * (zm - z) + g * nz;

      NRG += SQR(val - mean);

      v->tx = x + alpha * dx;
      v->ty = y + alpha * dy;
      v->tz = z + alpha * dz;
    }
    NRG = sqrt(NRG / nvertices);

    //    fprintf(WHICH_OUTPUT,"-  NRG = %f  -",NRG);

    /* update orig vertices */
    for (i = 0; i < nvertices; i++) {
      int const vno = vertices[i];
      VERTEX * const v = &mris->vertices[vno];
      MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
    }

    /* recompute normals */
    computeDefectFaceNormals(mris, dp, NULL);
    computeDefectVertexNormals(mris, dp);
  }
  if (vertices) {
    free(vertices);
  }
}

static void defectMaximizeLikelihood_new(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha);
static void defectMaximizeLikelihood_old(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha);

static void defectMaximizeLikelihood(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha) {
    return
#if 1
        defectMaximizeLikelihood_new
#else
        defectMaximizeLikelihood_old
#endif
            (mri, mris, dp, niter, alpha);
}

static void defectMaximizeLikelihood_new(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha)
{
  float const wm   = dp->defect->white_mean;
  float const gm   = dp->defect->gray_mean;
  float const mean = (wm + gm) * 0.5f;

  int const nvertices = dp->tp.ninside; /* matching only for inside vertices */

  DefectFacesCache defectFacesCache;
  initDefectFacesCache(&defectFacesCache);
  while (niter--) {

    /* using the tmp vertices */

    int i;
    for (i = 0; i < nvertices; i++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
      VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
      float const x = v->origx;
      float const y = v->origy;
      float const z = v->origz;

      /* smoothness term */
      double xm = 0, ym = 0, zm = 0;
      if (vt->vnum > 0) {
        int n;
        for (n = 0; n < vt->vnum; n++) {
          VERTEX const * const vn = &mris->vertices[vt->v[n]];
          xm += vn->origx;
          ym += vn->origy;
          zm += vn->origz;
        }

        xm *= 1/(double)n;
        ym *= 1/(double)n;
        zm *= 1/(double)n;
      }
      
      /* image gradient */
      double const nx = v->nx;
      double const ny = v->ny;
      double const nz = v->nz;

      double xv_lo,  yv_lo,  zv_lo;
      double xv_mid, yv_mid, zv_mid;
      double xv_hi,  yv_hi,  zv_hi;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel                  (x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv_lo,  &yv_lo,  &zv_lo);
      mriSurfaceRASToVoxel                  (x,            y,            z,            &xv_mid, &yv_mid, &zv_mid);
      mriSurfaceRASToVoxel                  (x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv_hi,  &yv_hi,  &zv_hi);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv_lo,  &yv_lo,  &zv_lo);
      MRISsurfaceRASToVoxelCached(mris, mri, x,            y,            z,            &xv_mid, &yv_mid, &zv_mid);
      MRISsurfaceRASToVoxelCached(mris, mri, x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv_hi,  &yv_hi,  &zv_hi);
#endif

      double white_val, val, gray_val;
      MRIsampleVolume(mri, xv_lo,  yv_lo,  zv_lo,  &white_val);
      MRIsampleVolume(mri, xv_mid, yv_mid, zv_mid, &val);
      MRIsampleVolume(mri, xv_hi,  yv_hi,  zv_hi,  &gray_val);

      double g = (white_val - gray_val) * (val - mean);
      if (fabs(g) > 0.2) {
        g = SIGN(g)*0.2;
      }

      double const dx = 0.5 * (xm - x) + g * nx;
      double const dy = 0.5 * (ym - y) + g * ny;
      double const dz = 0.5 * (zm - z) + g * nz;

      v->tx = x + alpha * dx;
      v->ty = y + alpha * dy;
      v->tz = z + alpha * dz;
    }

    /* update orig vertices */
    for (i = 0; i < nvertices; i++) {
      int const vno = dp->tp.vertices[i];
      VERTEX* const v = &mris->vertices[vno];
      MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
    }

    /* recompute normals */
    computeDefectFaceNormals  (mris, dp, &defectFacesCache);
    computeDefectVertexNormals(mris, dp);
  } // next iteration
  finiDefectFacesCache(&defectFacesCache);
}

static void defectMaximizeLikelihood_old(MRI *mri, MRI_SURFACE *mris, DP *dp, int niter, double alpha)
{
  float wm, gm, mean;
  int i, n, nvertices;
  double x, y, z, xm, ym, zm, nx, ny, nz, dx, dy, dz, g, NRG;
  double xv, yv, zv, white_val, gray_val, val;

  wm = dp->defect->white_mean;
  gm = dp->defect->gray_mean;

  mean = (wm + gm) / 2.0;

  nvertices = dp->tp.nvertices;
  nvertices = dp->tp.ninside; /* matching only for inside vertices */

  while (niter--) {
    /* using the tmp vertices */
    for (NRG = 0, i = 0; i < nvertices; i++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[dp->tp.vertices[i]];
      VERTEX                * const v  = &mris->vertices         [dp->tp.vertices[i]];
      x = v->origx;
      y = v->origy;
      z = v->origz;

      /* smoothness term */
      for (xm = 0, ym = 0, zm = 0, n = 0; n < vt->vnum; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]];
        xm += vn->origx;
        ym += vn->origy;
        zm += vn->origz;
      }
      if (n) {
        xm /= (double)n;
        ym /= (double)n;
        zm /= (double)n;
      }

      /* image gradient */
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - 0.5 * nx, y - 0.5 * ny, z - 0.5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &white_val);

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);

#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x + 0.5 * nx, y + 0.5 * ny, z + 0.5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &gray_val);

      g = (white_val - gray_val) * (val - mean);

      if (fabs(g) > 0.2) {
        g = 0.2 * g / fabs(g);
      }

      dx = 0.5 * (xm - x) + g * nx;
      dy = 0.5 * (ym - y) + g * ny;
      dz = 0.5 * (zm - z) + g * nz;

      NRG += SQR(val - mean);

      v->tx = x + alpha * dx;
      v->ty = y + alpha * dy;
      v->tz = z + alpha * dz;
    }
    NRG = sqrt(NRG / nvertices);

    //    fprintf(WHICH_OUTPUT,"-  %f  -",NRG);

    /* update orig vertices */
    for (i = 0; i < nvertices; i++) {
      int const vno = dp->tp.vertices[i];
      VERTEX * const v = &mris->vertices[vno];
      MRISsetOriginalXYZ(mris, vno, v->tx, v->ty, v->tz);
    }

    /* recompute normals */
    computeDefectFaceNormals(mris, dp, NULL);
    computeDefectVertexNormals(mris, dp);
  }
}

static void detectDefectFaces(MRIS *mris, DEFECT_PATCH *dp)
{
  int i, n, m, vno, vn1, vn2, nvertices, nthings, *things, nfaces;
  int optimal = dp->defect->optimal_mapping;
  TP *tp;

  /* the tessellated patch */
  tp = &dp->tp;

  /* in 'theory', the euler number of the patch is one */
  nfaces = tp->nedges - dp->tp.ninside + 101; /* add 100 just to make sure */
  things = (int *)malloc(nfaces * sizeof(nfaces));

  /* will use the border flag to mark modified vertices */
  nvertices = tp->nvertices;
#if 0
  for (nthings = i = 0 ; i < nvertices ; i++)
  {
    vno=tp->vertices[i];
    VERTEX * const v = &mris->vertices[vno] ;
    v->border=1;
  }
#endif
  /* detect faces only for modified vertices */
  for (nthings = i = 0; i < nvertices; i++) {
    vno = tp->vertices[i];
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

    for (n = 0; n < v->vnum; n++) {
      vn1 = v->v[n];

      //                        if(!mris->vertices[vn1].border) continue; /* only modified vertices */

      if (optimal && mris->vertices[vn1].fixedval == 0) {
        continue;  // experimental
      }

      for (m = 0; m < v->vnum; m++) {
        vn2 = v->v[m];
        
        // edges should not go between the same vno!
        // but this checks anyway...
        //
        if (vn1 == vn2) {
          continue;
        }

        // if(!mris->vertices[vn2].border) continue; /* only modified vertices */
        if (optimal && mris->vertices[vn2].fixedval == 0) {
          continue;  // experimental
        }

        /* check if this set of vertices could constitue a face */
        if (!vertexNeighbor(mris, vn1, vn2)) {
          continue;
        }

        /* check if this triangle already in the tessellation */
        if (isFace(mris, vno, vn1, vn2)) {
          continue;
        }

        /* check if this potential face contains other vertices */
        if (containsAnotherVertexOnSphere(mris, vno, vn1, vn2, optimal)) {
          continue;
        }

        if (!mrisCanAttachFaceToVertices(mris, vno, vn1, vn2)) {
          static unsigned long count,limit = 1;
          if (count++ >= limit) {
            limit *= 2;
            fs::debug() << "suppressed badly attaching a face, count: " << count;
          }
          continue;
        }

        /* add this new face to the defect faces */
        mrisAddFace(mris, vno, vn1, vn2);
        if (nthings == nfaces) {
          continue;
        }
        things[nthings++] = mris->nfaces - 1;
      }
    }
  }

  /* save the list of new faces */
  if (nthings == nfaces) {
    fprintf(WHICH_OUTPUT, "error in the retessellation \n");
  }
  tp->faces = (int *)malloc(nthings * sizeof(int));
  tp->nfaces = nthings;
  memmove(tp->faces, things, nthings * sizeof(int));
  free(things);
}

#define DEBUG_INFO 0

static int computePatchEulerNumber(MRIS *mris, DP *dp)
{
  int nfaces, nedges, nvertices, euler;

  /* everything has already been computed */
  nfaces = dp->tp.nfaces;
  nedges = dp->tp.nedges;
  nvertices = dp->tp.ninside;
  euler = nvertices - nedges + nfaces;

#if DEBUG_INFO
  fprintf(WHICH_OUTPUT, "euler=%d : (%d,%d,%d)\n", euler, nvertices, nedges, nfaces);
  if (euler != 1) {
    fprintf(WHICH_OUTPUT, "\n\nXXXXXXXXXXXXXXXXXXXXXXXX\n\n");
  }

  if (euler != 1)
    fprintf(WHICH_OUTPUT,
            "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n"
            "euler=%d : (%d,%d,%d)\n"
            "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n",
            euler,
            nvertices,
            nedges,
            nfaces);
#endif
  return euler;
}

static void orientDefectFaces(MRIS *mris, DP *dp)
{
  int n, m, vno0, vno1, fno;
  float dot, cx, cy, cz, a[3], b[3], norm[3];
  TP *tp;
  FACE *face;

  tp = &dp->tp;

  for (n = 0; n < tp->nfaces; n++) {
    fno = tp->faces[n]; /* face index in mris_corrected */

    face = &mris->faces[fno];

    VERTEX const * const v1 = &mris->vertices[face->v[0]];
    VERTEX const * const v2 = &mris->vertices[face->v[1]];
    VERTEX const * const v3 = &mris->vertices[face->v[2]];

    /* compute centroid direction on sphere */
    cx = v1->cx + v2->cx + v3->cx;
    cy = v1->cy + v2->cy + v3->cy;
    cz = v1->cz + v2->cz + v3->cz;

    /* compute normal of face onto sphere*/
    a[0] = v2->cx - v1->cx;
    b[0] = v3->cx - v1->cx;
    a[1] = v2->cy - v1->cy;
    b[1] = v3->cy - v1->cy;
    a[2] = v2->cz - v1->cz;
    b[2] = v3->cz - v1->cz;

    F_CROSS(a, b, norm);
    dot = norm[0] * cx + norm[1] * cy + norm[2] * cz;

    if (dot < 0) {
      /* they disagree - change order of vertices 1 & 2 in face n */
      vno0 = face->v[1];
      vno1 = face->v[2];
      face->v[1] = vno1;
      face->v[2] = vno0;

      /* set vertex face index */
      VERTEX_TOPOLOGY * v = &mris->vertices_topology[vno0]; /* vno0 is now in 2 */
      for (m = 0; m < v->num; m++)
        if (v->f[m] == fno) {
          v->n[m] = 2;
          break;
        }

      v = &mris->vertices_topology[vno1]; /* vno1 is now in 1 */
      for (m = 0; m < v->num; m++)
        if (v->f[m] == fno) {
          v->n[m] = 1;
          break;
        }
    }
  }
}

/* used to temporary rip the faces of the defect so we don't process them many times */
// FLO TO BE CHECKED
#define TEMPORARY_RIPPED_FACE 2

static void computeDefectFaceNormals(MRIS *mris, DP *dp, DefectFacesCache* dfc)
{
  int i, n;
  FACE *face;
  TP *tp;

  /* the tessellated patch */
  tp = &dp->tp;

  if (dfc && dfc->size > 0) {
    int const * const fnos = dfc->fnos;
    for (i = 0; i < dfc->size; i++) {
      computeDefectFaceNormal(mris, fnos[i]);   // measurement says this happens about 40 times per initDefectCache
    }
    return;
  }

  /* compute faces only for modified vertices */
  for (i = 0; i < tp->nvertices; i++) {
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[tp->vertices[i]];
    for (n = 0; n < v->num; n++) {
      face = &mris->faces[v->f[n]];
      if (face->ripflag) {
        continue; /* don't process a face twice */
      }
      computeDefectFaceNormal(mris, v->f[n]);
      if (dfc) insertIntoDefectFacesCache(dfc,v->f[n]); 
      face->ripflag = TEMPORARY_RIPPED_FACE;
    }
  }
  /* unrip faces */
  for (i = 0; i < tp->nvertices; i++) {
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[tp->vertices[i]];
    for (n = 0; n < v->num; n++) {
      face = &mris->faces[v->f[n]];
      if (face->ripflag == TEMPORARY_RIPPED_FACE) /* unrip face */
      {
        face->ripflag = 0;
      }
    }
  }
}

static void computeDefectVertexNormals(MRIS *mris, DP *dp)
{
  TP* const tp = &dp->tp;
  
  /* compute vertex normals only for modified vertices */
  for (int n = 0; n < tp->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[tp->vertices[n]];
    VERTEX                * const v  = &mris->vertices         [tp->vertices[n]];

    /* compute normal at vertex */
    float nx = 0.0f, ny = 0.0f, nz = 0.0f;
    for (int m = 0; m < vt->num; m++) {
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, vt->f[m]);
      nx += fNorm->nx;
      ny += fNorm->ny;
      nz += fNorm->nz;
    }

    /* normalize */
    float len = sqrtf(nx * nx + ny * ny + nz * nz);     // this should have been double precision
                                                        // but fixing might change results

    if (FZERO(len)) {

      // BRB it is possible the patch doesn't attach any faces to a vertex
      //
      if (vt->num > 0) {
        fprintf(WHICH_OUTPUT,
              "normal vector of length zero at vertex %d with %d faces\n",
              tp->vertices[n],
              (int)vt->num);  // TO BE CHECKED
      }

      // if ((int)vt->num == 0) {
      //  ErrorExit(ERROR_BADPARM, "vertex %d has 0 face", tp->vertices[n]);
      // }
      
      len = 1.0f;
    }
    
    v->nx = nx / len;
    v->ny = ny / len;
    v->nz = nz / len;
  }
}

static void printDefectStatistics(DP *dp)
{
  // fprintf(WHICH_OUTPUT,"fitness=%2.4f\n",dp->fitness);
  fprintf(WHICH_OUTPUT,
          "X=%d (v=%d,e=%d,f=%d): %d vertices were discarded\n",
          dp->tp.ninside + dp->tp.nfaces - dp->tp.nedges,
          dp->tp.ninside,
          dp->tp.nedges,
          dp->tp.nfaces,
          dp->tp.ndiscarded);
  fprintf(WHICH_OUTPUT,
          "fll=%2.4f (%2.4f), vll=%2.4f (%2.4f),cll=%2.4f (%2.4f),qcll=%2.4f (%2.4f) umll=%2.4f (%2.4f)\n",
          dp->tp.face_ll,
          (dp->tp.face_ll - dp->defect->initial_face_ll),
          dp->tp.vertex_ll,
          (dp->tp.vertex_ll - dp->defect->initial_vertex_ll),
          dp->tp.curv_ll,
          (dp->tp.curv_ll - dp->defect->initial_curv_ll),
          dp->tp.qcurv_ll,
          (dp->tp.qcurv_ll - dp->defect->initial_qcurv_ll),
          dp->tp.unmri_ll,
          (dp->tp.unmri_ll - dp->defect->initial_unmri_ll));
}

/* check if the vertex vno1 is contained in one of the faces of vno2 */
static int isDiscarded(MRIS *mris, int vno1, int vno2)
{
  int n;
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno2];

  for (n = 0; n < v->num; n++)
    if (isVertexInsideFace(mris, vno1, v->f[n])) /* vno1 is inside face v->f[n]! */
    {
      return 1;
    }

  return 0;
}

static void removeVertex(MRIS *mris, int vno)
{
  int n, m, vnum, *oldlist;
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
  VERTEX          * const v  = &mris->vertices         [vno];

  // to be checked when the ADD_SOME_VERTICES mode is on
  // could be a good "outside" vertex

  if (v->marked == 0) {
    /*          if(parms->verbose>=VERBOSE_MODE_MEDIUM)
                fprintf(stderr,"removeVertex: the vertex %d was not marked - SHOULD NOT HAPPEN\n",vno);
                if(parms->verbose==VERBOSE_MODE_HIGH)*/
    ErrorExit(ERROR_BADPARM, "removeVertex: the vertex %d was not marked - SHOULD NOT HAPPEN\n", vno);
    return;
  }

  for (n = 0; n < vt->vnum; n++) {
    /* remove vno from the list of v->v[n] */
    int const vno2 = vt->v[n];
    VERTEX_TOPOLOGY * const vnt = &mris->vertices_topology[vno2];
    VERTEX          * const vn  = &mris->vertices         [vno2];
    oldlist = vnt->v;
    vnum = vnt->vnum - 1; /* the new # of neighbors */
    if (vnum) {
      vnt->v = (int *)malloc(vnum * sizeof(int));
      for (vnum = m = 0; m < vnt->vnum; m++) {
        if (oldlist[m] == vno) {
          continue;
        }
        vnt->v[vnum++] = oldlist[m];
      }
      free(oldlist);
    }
    else {
      vnt->v = NULL;
    }
    
    modVnum(mris,vno2,vnum,true);
    vnt->nsizeMax = 1;
    
    MRIS_setNsizeCur(mris, vno2, 1);

    /* check if the vertex became singled out */
    if (vnt->vnum == 0) {
      vn->marked = INSIDE_VERTEX;
    }
  }
  
  v->marked = DISCARDED_VERTEX;
  freeAndNULL(vt->v);
  clearVnum(mris,vno);
  vt->v2num = 0;
  vt->v3num = 0;
  vt->vtotal = 0;
}

static int updateVertexTriangle(MRIS *mris, int vno, int fno)
{
  int n, m, vn;
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

  for (m = n = 0; n < v->vnum; n++) {
    /* work with the vertex v->v[n] */
    vn = v->v[n];
    if (mris->vertices[vn].marked == TRIANGLE_VERTEX) {
      continue;
    }
    if (isVertexInsideFace(mris, vn, fno)) {
      /* vertex is inside the new face, so we have to remove it */
      removeVertex(mris, vn);
      n--;
      m++;
    }
  }

  /* we have removed m vertices from vno->vnum */
  return m;
}

/* static void updateTriangle(MRIS *mris,int fno){ */
/*   int n,m,vn; */
/*   VERTEX *v; */

/*   /\* update type *\/ */
/*   for( m = 0 ; m < 3; m++) */
/*     mris->vertices[mris->faces[fno].v[m]].marked=TRIANGLE_VERTEX; */

/*   for( m = 0 ; m < 3; m++){ */
/*     v=&mris->vertices[mris->faces[fno].v[m]]; */
/*     for( n = 0 ; n < v->vnum ; n++){ */
/*       /\* work with the vertex v->v[n] *\/ */
/*       vn=v->v[n]; */
/*       if(mris->vertices[vn].marked==TRIANGLE_VERTEX) continue; */
/*       if(isVertexInsideFace(mris,vn,fno)){ */
/*      /\* vertex is inside the new face, so we have to remove it *\/ */
/*      removeVertex(mris,vn); */
/*      n--; */
/*       } */
/*     } */
/*   } */
/* } */

/* check if the edge vno1 <--> vno2 is the edge of a triangle */
static void possiblyAddNewFaces(MRIS *mris, int vno1, int vno2)
{
  int n, m, vn, fno;

  /* only exterior triangles */

  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno1];
  for (n = 0; n < v->vnum; n++) {
    if (v->v[n] == vno2) {
      continue;
    }
    if (vertexNeighbor(mris, vno2, v->v[n])) {
      /* there is a triangle */
      vn = v->v[n];
      /* make sure this vertex is not outside : could happen if border <-> border edge */
      if (mris->vertices[vn].marked == OUTSIDE_VERTEX) {
        continue; /* outside vertex */
      }

      /* check if face already exists : could happen with border <-> border <-> border triangles */
      if (isFace(mris, vno1, vno2, vn)) {
        continue;
      }

      mrisAddFace(mris, vno1, vno2, vn);
      fno = mris->nfaces - 1;
      /* update type */
      mris->vertices[vno1].marked = TRIANGLE_VERTEX;
      mris->vertices[vno2].marked = TRIANGLE_VERTEX;
      mris->vertices[vn].marked = TRIANGLE_VERTEX;

      m = updateVertexTriangle(mris, vno1, fno);
      updateVertexTriangle(mris, vn, fno);
      updateVertexTriangle(mris, vno2, fno);

      /* we have removed m vertices from the vertex vno1 */
      n -= m;
    }
  }
}

static int possiblyAddEdgesAndFaces(MRIS *mris, int vno1, int vno2, int mode)
{
  int mark1, mark2, tmp;
  VERTEX *v1, *v2;

  if (mode == USE_ALL_VERTICES) {
    mrisAddEdge(mris, vno1, vno2);
    mris->vertices[vno1].marked = TRIANGLE_VERTEX;
    mris->vertices[vno2].marked = TRIANGLE_VERTEX;
    return 1;
  };

  /* the initial marks */
  mark1 = mris->vertices[vno1].marked;
  mark2 = mris->vertices[vno2].marked;

  /* swap vertices if necessary */
  if (mark2 < mark1) {
    tmp = vno1;
    vno1 = vno2;
    vno2 = tmp;
    tmp = mark1;
    mark1 = mark2;
    mark2 = tmp;
  }

  v1 = &mris->vertices[vno1];
  v2 = &mris->vertices[vno2];

  switch (mark1) {
    case INSIDE_VERTEX:
      if (mark2 == TRIANGLE_VERTEX) {
        /* is the first vertex inside or outside the faces of vno2*/
        if (isDiscarded(mris, vno1, vno2)) {
          /* inside one of the faces of vno2 */
          /* remove this vertex that is not connected to anyone */
          v1->marked = DISCARDED_VERTEX;
          return 0;
        }
        else {
          /* outside the faces */
          mrisAddEdge(mris, vno1, vno2);
          /* the first vertex becomes EDGE */
          v1->marked = EDGE_VERTEX;
          return 1;
        }
      }
      else {
        /* the second vertex is either of type INSIDE or EDGE */
        mrisAddEdge(mris, vno1, vno2);
        /* they both become EDGE */
        v1->marked = EDGE_VERTEX;
        v2->marked = EDGE_VERTEX;
        return 1;
      }
      break;
    case EDGE_VERTEX:
      if (mark2 == EDGE_VERTEX) {
        /* add this edge */
        mrisAddEdge(mris, vno1, vno2);
        /* look for new triangles */
        possiblyAddNewFaces(mris, vno1, vno2);
        return 1;
      }
      else {
        /* vno2 is of type TRIANGLE_VERTEX */
        /* is the first vertex inside or outside the faces of vno2 */
        if (isDiscarded(mris, vno1, vno2)) {
          /* inside one of the faces */
          /* remove this vertex and its connections */
          removeVertex(mris, vno1);
          return 0;
        }
        else {
          /* outside */
          mrisAddEdge(mris, vno1, vno2);
          /* check if some new triangles have been formed : if yes, update type*/
          possiblyAddNewFaces(mris, vno1, vno2);
          return 1;
        }
      }
      break;
    case TRIANGLE_VERTEX:
      /* both vertices have the type TRIANGLE_VERTEX */
      mrisAddEdge(mris, vno1, vno2);
      /* check if some new triangles have been formed */
      possiblyAddNewFaces(mris, vno1, vno2);
      return 1;
      break;
  }

  return 0;
};

typedef struct
{
  int intersected;
  int vno1, vno2;
} INTERSECTION_TABLE, IT;

static int retessellateDefect_wkr(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, DVS *dvs, DP *dp);

static int retessellateDefect(MRI_SURFACE *mris, MRI_SURFACE *mris_corrected, DVS *dvs, DP *dp)
{
    int result;
    ROMP_SCOPE_begin
    result =  retessellateDefect_wkr(mris, mris_corrected, dvs, dp);
    ROMP_SCOPE_end
    return result;
}

static int retessellateDefect_wkr(MRIS *mris, MRIS *mris_corrected, DVS *dvs, DP *dp)
{
  mrisCheckVertexFaceTopology(mris);
  mrisCheckVertexFaceTopology(mris_corrected);

  static bool const showStats = false;

  double max_len;
  int i, j, max_i, max_added, nadded, index, ndiscarded;
  int (*intersection_function)(MRI_SURFACE * mris, DEFECT * defect, EDGE * e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2);
  int *vertex_trans;
  DEFECT *defect;
  EDGE *et;
  EDGE_TABLE *etable;
  int nedges, *ordering, n, vno;
  int nthings, *things;
  int modified;
  IT *it;
  int counting;

  /* initialize arrays of tessellated patch to null pointer*/
  TPinit(&dp->tp);

  etable = dp->etable;
  defect = dvs->defect;
  vertex_trans = dvs->vertex_trans;
  et = dp->etable->edges;
  nedges = dp->nedges;
  ordering = dp->ordering;

  /* for the ordering the vertices - using undefval*/
  // useful only if dp->retessellation_mode==USE_ALL_VERTICES
  for (counting = 0, n = 0; n < defect->nvertices; n++) {
    if (defect->status[n] == DISCARD_VERTEX) {
      continue;
    }
    vno = vertex_trans[defect->vertices[n]];
    VERTEX_TOPOLOGY const * const vt = &mris_corrected->vertices_topology[vno];
    VERTEX                * const v  = &mris_corrected->vertices         [vno];
    v->undefval = 0;
    if (vt->vnum) {
      v->undefval = 1;
      counting = 1;
    }
  }

  max_len = 0;
  max_i = 0;
  max_added = 0;
  intersection_function = intersectDefectEdges;

  /* first mark the defective vertices */
  for (ndiscarded = n = 0; n < defect->nvertices; n++) {
    if (defect->status[n] == DISCARD_VERTEX) {
      ndiscarded++;
      continue;
    }
    vno = vertex_trans[defect->vertices[n]];
    mris_corrected->vertices[vno].marked = INSIDE_VERTEX;
  }
  for (n = 0; n < defect->nborder; n++) {
    vno = vertex_trans[defect->border[n]];
    mris_corrected->vertices[vno].marked = BORDER_VERTEX;
  }

  /* count the number of potentially added edges */
  nthings = nadded = 0;

  /* allocate the table of potentially intersected edges */
  it = (IT *)calloc(nedges, sizeof(IT));

  static long stats_count = 0;
  static long stats_limit = 1;
  long stats_modified_loops   = 0;
  long stats_nedges_loops_all = 0, stats_nedges_loops_heavy = 0, stats_intersection_function_calls = 0;
  
  modified = 1;
  while (modified) {
    modified = 0;

    stats_modified_loops++;
    
    IntersectDefectEdgesContext intersectDefectEdgesContext; 
    initIntersectDefectEdgesContext(&intersectDefectEdgesContext, mris);
    
    /* start building the retessellation */
    for (index = 0; index < nedges; index++) {
      stats_nedges_loops_all++;
      
      if (ordering) {
        i = ordering[index];
      }
      else {
        i = index;
      }

      if (it[i].intersected) {
        continue;
      }

      if (et[i].used && et[i].used != USED_IN_ORIGINAL_TESSELLATION) /* already exists in
                                                                        tessellation - don't
                                                                        add it again */
      {
        continue; /* edge status must be USED_IN_TESSELLATION */
      }

      /* check if this edge really exists */
      if (mris_corrected->vertices[et[i].vno1].marked == DISCARDED_VERTEX) {
        continue;
      }
      if (mris_corrected->vertices[et[i].vno2].marked == DISCARDED_VERTEX) {
        continue;
      }

      // TO BE CHECKED (used in RandomRetessellation)
      if (mris_corrected->vertices[et[i].vno1].ripflag) {
        continue;
      }
      if (mris_corrected->vertices[et[i].vno2].ripflag) {
        continue;
      }

      stats_nedges_loops_heavy++;

      if (etable && etable->use_overlap == USE_OVERLAP) /* use pre-computed
                                                           intersection table */
      {
        int intersects = 0;

        for (j = 0; j < etable->noverlap[i]; j++)
          if (et[etable->overlapping_edges[i][j]].used &&
              et[etable->overlapping_edges[i][j]].used != USED_IN_ORIGINAL_TESSELLATION) {
            intersects = 1;
            /* the edge i was refused because of the edge
               etable->overlapping_edges[i][j] */
            it[i].intersected = 2;
            it[i].vno1 = et[etable->overlapping_edges[i][j]].vno1;
            it[i].vno2 = et[etable->overlapping_edges[i][j]].vno2;
            break;
          }
          
        if (intersects) {
          continue;
        }
        if (etable->flags[i] & ET_OVERLAP_LIST_INCOMPLETE) {
          intersection_function = intersectDefectEdges;
        }
        else {
          intersection_function = intersectDefectConvexHullEdges;
        }
      }
      
      stats_intersection_function_calls++;
        // almost always comes this way
        
      if ((*intersection_function)(mris_corrected, defect, &et[i], &intersectDefectEdgesContext, vertex_trans, &it[i].vno1, &it[i].vno2) == 0) {
        /* this edge could potentially be added : no sphere intersection */
        nadded++;

        if (possiblyAddEdgesAndFaces(mris_corrected, et[i].vno1, et[i].vno2, dp->retessellation_mode)) {
          obsoleteIntersectDefectEdgesContext(&intersectDefectEdgesContext);
          nthings++;
          if (et[i].used) /* used in original tessellation */
          {
            et[i].used = USED_IN_BOTH_TEMPORARY_TESSELLATION;
          }
          else {
            et[i].used = USED_IN_NEW_TEMPORARY_TESSELLATION;
          }
          /* useful only if mode = USE_ALL_VERTICES */
          VERTEX * const vertex1 = &mris_corrected->vertices[et[i].vno1];
          VERTEX * const vertex2 = &mris_corrected->vertices[et[i].vno2];
          if (vertex1->undefval == 0 || vertex2->undefval == 0) {
            counting++;
            if (vertex1->undefval == 0) {
              vertex1->undefval = counting;
            }
            if (vertex2->undefval == 0) {
              vertex2->undefval = counting;
            }
          }
        }
        
      }
      else /* intersecting edge with edge e1<-->e2 */
      {
        it[i].intersected = 2;
      }
    }

    finiIntersectDefectEdgesContext(&intersectDefectEdgesContext);

    /* now update the edges */
    for (index = 0; index < nedges; index++) {
      /* keep the same order (not necessary) */
      if (ordering) {
        i = ordering[index];
      }
      else {
        i = index;
      }

      /* already exists in tessellation - don't add it again */
      if (et[i].used == USED_IN_NEW_TEMPORARY_TESSELLATION) {
        et[i].used = NOT_USED;
        if (mris_corrected->vertices[et[i].vno1].marked != TRIANGLE_VERTEX) {
          continue;
        }
        if (mris_corrected->vertices[et[i].vno2].marked != TRIANGLE_VERTEX) {
          continue;
        }
        et[i].used = USED_IN_NEW_TESSELLATION;
      };
      if (et[i].used == USED_IN_BOTH_TEMPORARY_TESSELLATION) {
        et[i].used = USED_IN_ORIGINAL_TESSELLATION;
        if (mris_corrected->vertices[et[i].vno1].marked != TRIANGLE_VERTEX) {
          continue;
        }
        if (mris_corrected->vertices[et[i].vno2].marked != TRIANGLE_VERTEX) {
          continue;
        }
        et[i].used = USED_IN_BOTH_TESSELLATION;
      };

      if (it[i].intersected == 2) {
        if (mris_corrected->vertices[it[i].vno1].marked == DISCARDED_VERTEX ||
            mris_corrected->vertices[it[i].vno2].marked == DISCARDED_VERTEX) {
          modified = 1;
          it[i].intersected = 0;
        }
        else {
          it[i].intersected = 1;
        }
      }
      else {
        it[i].intersected = 1;
      }
    }
    
    mrisCheckVertexFaceTopology(mris_corrected);
  }

  if (++stats_count >= stats_limit) {
    stats_limit *= 2;
    if (showStats) {
      fprintf(stderr, "%s:%d stats_count:%ld "
        "modified_loops:%ld nedges_loops_all:%ld "
        "nedges_loops_heavy:%ld intersection_function_calls:%ld\n",
        __FILE__, __LINE__, stats_count, 
        stats_modified_loops, stats_nedges_loops_all, 
        stats_nedges_loops_heavy, stats_intersection_function_calls);
    }
  }

  /* in this retessellation we have added, at most, nthings edges */
  things = (int *)malloc(nthings * sizeof(int));
  /* then, note the used edges */
  for (nthings = index = 0; index < nedges; index++) {
    /* keep the same order (not necessary) */
    if (ordering) {
      i = ordering[index];
    }
    else {
      i = index;
    }

    /* already exists in tessellation - don't add it again */
    if (et[i].used != USED_IN_NEW_TESSELLATION && et[i].used != USED_IN_BOTH_TESSELLATION) {
      continue;
    }

    things[nthings++] = i;
    if (et[i].len > max_len) {
      max_len = et[i].len;
      max_added = nadded - 1;
      max_i = i;
    }
  }
  
  /* store list of used edges */
  dp->tp.nedges = nthings;
  dp->tp.edges = (int *)malloc(nthings * sizeof(int));
  memmove(dp->tp.edges, things, nthings * sizeof(int));
  free(things);

#if DEBUG_INFO
  fprintf(WHICH_OUTPUT, "%d added edges out of %d potential edges\n", nthings, nadded);
#endif

  /* mark the used vertices */
  for (n = 0; n < dp->tp.nedges; n++) {
    i = dp->tp.edges[n];
    mris_corrected->vertices[et[i].vno1].marked = USED_VERTEX;
    mris_corrected->vertices[et[i].vno2].marked = USED_VERTEX;
  }

  /* in this retessellation we have added, at most, nthings vertices */
  nthings = defect->nvertices + defect->nborder;
  things = (int *)malloc(nthings * sizeof(int));

  /* count the number of added vertices and reset marks to zero*/
  for (nthings = n = 0; n < defect->nvertices; n++) {
    if (defect->status[n] == DISCARD_VERTEX) {
      continue;
    }
    vno = vertex_trans[defect->vertices[n]];
    if (mris_corrected->vertices[vno].marked == USED_VERTEX) {
      things[nthings++] = vno;
    }
    mris_corrected->vertices[vno].marked = 0;
  }
  dp->tp.ninside = nthings; /* we have ninside new vertices */
  for (n = 0; n < defect->nborder; n++) {
    /* every border voxel counts !!!*/
    vno = vertex_trans[defect->border[n]];
    things[nthings++] = vno;
    mris_corrected->vertices[vno].marked = 0;
  }
  dp->tp.nvertices = nthings;
  dp->tp.vertices = (int *)malloc(nthings * sizeof(int));
  memmove(dp->tp.vertices, things, nthings * sizeof(int));
  free(things);

  /* number of discarded vertices */
  dp->tp.ndiscarded = defect->nvertices - ndiscarded - dp->tp.ninside;

#if DEBUG_INFO
  fprintf(WHICH_OUTPUT,
          "%d vertices have been used out of %d(%d): %d were discarded,\n"
          "%d out of %d(%d) are new vertices: %d were discarded,\n"
          "%d out of %d are border vertices : %d were discarded\n",
          dp->tp.nvertices,
          defect->nvertices + defect->nborder - ndiscarded,
          defect->nvertices + defect->nborder,
          defect->nvertices + defect->nborder - ndiscarded - dp->tp.nvertices,
          dp->tp.ninside,
          defect->nvertices - ndiscarded,
          defect->nvertices,
          defect->nvertices - ndiscarded - dp->tp.ninside,
          dp->tp.nvertices - dp->tp.ninside,
          defect->nborder,
          defect->nborder - dp->tp.nvertices + dp->tp.ninside);
#endif

  /* reset the number of original faces in the surface before the retessellation */
  if (dp->retessellation_mode == USE_SOME_VERTICES) {
    mrisRestoreFaceVertexState(mris_corrected, dvs);
    mrisCheckVertexFaceTopology(mris_corrected);
  }

  /* free the allocated memory for the intersection_table */
  free(it);

  return (NO_ERROR);
}

static void computeDisplacement(MRI_SURFACE *mris, DP *dp)
{
  DEFECT *defect;
  int p, vno;
  float val;
  VERTEX *vertex;

  defect = dp->defect;
  for (p = 0; p < defect->nvertices; p++) {
    if (defect->status[p] == DISCARD_VERTEX) {
      continue;
    }
    vno = defect->vertex_trans[defect->vertices[p]];
    vertex = &mris->vertices[vno];
    val = (SQR(vertex->fx - vertex->origx) + SQR(vertex->fy - vertex->origy) + SQR(vertex->fz - vertex->origz));
    vertex->curvbak = sqrt(val);
  }
}

static double mrisDefectPatchFitness(
    ComputeDefectContext* computeDefectContext,
    MRI_SURFACE *mris,
    MRI_SURFACE *mris_corrected,
    MRI *mri,
    DEFECT_PATCH *dp,
    int *vertex_trans,
    DEFECT_VERTEX_STATE *dvs,
    RP *rp,
    HISTOGRAM *h_k1,
    HISTOGRAM *h_k2,
    MRI *mri_k1_k2,
    HISTOGRAM *h_white,
    HISTOGRAM *h_gray,
    HISTOGRAM *h_border,
    HISTOGRAM *h_grad,
    MRI *mri_gray_white,
    HISTOGRAM *h_dot,
    TOPOLOGY_PARMS *parms)
{
  int i, euler;
  VERTEX *v;
  DEFECT *defect = dp->defect;

  defect->vertex_trans = vertex_trans;
  dp->verbose_mode = parms->verbose;

  /* set the arrays to NULL in dp->tp */
  TPinit(&dp->tp);

  retessellateDefect(mris, mris_corrected, dvs, dp);    // BEVIN mris_fix_topology
  mrisCheckVertexFaceTopology(mris_corrected);

  /* detect the new set of faces */
  detectDefectFaces(mris_corrected, dp);

  /* compute the euler number of the patch */
  euler = computePatchEulerNumber(mris_corrected, dp);

  /* orient the patch faces */
  orientDefectFaces(mris_corrected, dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* save original coord into flattened coordinates */
  for (i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    v = &mris_corrected->vertices[vertex_trans[defect->vertices[i]]];
    //          fprintf(WHICH_OUTPUT,"%d-",vertex_trans[defect->vertices[i]]);
    v->fx = v->origx;
    v->fy = v->origy;
    v->fz = v->origz;
  }

  /* smooth and match original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, dp, parms->smooth, parms->match);

  computeDisplacement(mris_corrected, dp);

  /* compute the face normals on the surface*/
  computeDefectFaceNormals(mris_corrected, dp, NULL);

  /* compute vertex normals on original surface */
  // computeDefectVertexNormals(mris_corrected,dp);

  /* compute the patch fitness */
  dp->fitness = mrisComputeDefectLogLikelihood(
      computeDefectContext,
      mris_corrected, mri, dp, h_k1, h_k2, mri_k1_k2, h_white, h_gray, h_border, h_grad, mri_gray_white, h_dot, parms);

  /* update statistics */
  updateVertexStatistics(mris, mris_corrected, dvs, rp, dp, vertex_trans, dp->fitness);

  /* restore the vertex state */
  mrisRestoreVertexState(mris_corrected, dvs);

  /* reset the edges to the unused state (unless they were in the original tessellation) */
  for (i = 0; i < dp->nedges; i++) {
    if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION) {
        dp->etable->edges[i].used  = NOT_USED;
    }
    if (dp->etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
        dp->etable->edges[i].used  = USED_IN_ORIGINAL_TESSELLATION;
    }
  }
  /* free vertices,edges,faces tables */
  TPfree(&dp->tp);

  return (dp->fitness);
}


#define MAX_DEFECT_VERTICES 900000
static long ncross = 0;
static long nmut = 0;
static long nkilled = 0;

static DEFECT_LIST *mrisMergeNeighboringDefects(MRIS *mris, DEFECT_LIST *dl);

static DEFECT_LIST *mrisMergeNeighboringDefects(MRIS *mris, DEFECT_LIST *dl)
{
  int i, j, n, m, *nd, ndefects, merged;
  int vlist[MAX_DEFECT_VERTICES], nadded;
  float len;

  DEFECT *defect;
  DEFECT_LIST *new_dl;

  fprintf(WHICH_OUTPUT, "analyzing neighboring defects...\n");

  ndefects = dl->ndefects;
  nd = (int *)malloc((ndefects + 1) * sizeof(int));

  /* mark inside vertices */
  for (i = 0; i < ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nvertices; n++) {
      mris->vertices[defect->vertices[n]].marked = i + 1;
    }
  }

  /* iteratively find neighboring defects and merge them */
  for (i = 0; i < ndefects; i++) {
    memset(nd, 0, (ndefects + 1) * sizeof(int));
    defect = &dl->defects[i];
    for (n = 0; n < defect->nborder; n++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->border[n]];
      for (m = 0; m < vt->vnum; m++) {
        VERTEX const * const nv = &mris->vertices[vt->v[m]];
        if (nv->marked && nv->marked != i + 1)  // belong to another defect
        {
          nd[nv->marked]++;
        }
      }
    }

    /* find if this segment should be merged with other ones */
    merged = 0;
    for (j = 1; j < ndefects + 1; j++) {
      if (nd[j] >= 2) {
        merged = 1;
        fprintf(WHICH_OUTPUT, "      -merging segment %d into %d\n", j - 1, i);
        /* merging of two segments */

        defect->area = -1;
        /* update inside vertices */
        memmove(vlist, defect->vertices, defect->nvertices * sizeof(int));
        nadded = defect->nvertices;
        memmove(&vlist[nadded], dl->defects[j - 1].vertices, dl->defects[j - 1].nvertices * sizeof(int));
        nadded += dl->defects[j - 1].nvertices;

        free(defect->vertices);
        defect->vertices = (int *)malloc(nadded * sizeof(int));
        memmove(defect->vertices, vlist, nadded * sizeof(int));
        defect->nvertices = nadded;

        free(defect->status);
        defect->status = (char *)malloc(nadded * sizeof(char));
        for (n = 0; n < defect->nvertices; n++) {
          defect->status[n] = KEEP_VERTEX;  // TO BE CHECKED
        }

        /* free inside vertices */
        free(dl->defects[j - 1].vertices);
        free(dl->defects[j - 1].status);
        dl->defects[j - 1].vertices = NULL;
        dl->defects[j - 1].status = NULL;
        dl->defects[j - 1].nvertices = 0;

        /* free border */
        free(dl->defects[j - 1].border);
        dl->defects[j - 1].border = NULL;
        dl->defects[j - 1].nborder = 0;

        /* update defect parameters */
        defect->nx = 0;
        defect->ny = 0;
        defect->nz = 0;
        defect->cx = 0;
        defect->cy = 0;
        defect->cz = 0;
        defect->area = 0;
        for (n = 0; n < defect->nvertices; n++) {
          defect->nx += mris->vertices[defect->vertices[n]].nx;
          defect->ny += mris->vertices[defect->vertices[n]].ny;
          defect->nz += mris->vertices[defect->vertices[n]].nz;
          defect->cx += mris->vertices[defect->vertices[n]].x;
          defect->cy += mris->vertices[defect->vertices[n]].y;
          defect->cz += mris->vertices[defect->vertices[n]].z;
          defect->area += mris->vertices[defect->vertices[n]].origarea;
        }
        defect->nx /= (float)defect->nvertices;
        defect->ny /= (float)defect->nvertices;
        defect->nz /= (float)defect->nvertices;
        defect->cx /= (float)defect->nvertices;
        defect->cy /= (float)defect->nvertices;
        defect->cz /= (float)defect->nvertices;
        defect->area /= (float)defect->nvertices;
      }
    }

    /* update the border vertices and the mark for the inside ones */
    if (merged) {
      /* update marks inside */
      for (n = 0; n < defect->nvertices; n++) {
        mris->vertices[defect->vertices[n]].marked = i + 1;
      }

      /* update border vertices */
      nadded = 0;
      for (n = 0; n < defect->nvertices; n++) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[n]];
        for (m = 0; m < vt->vnum; m++) {
          VERTEX * const nv = &mris->vertices[vt->v[m]];
          if (nv->marked) {
            continue;
          }
          nv->marked = -1; /* border vertex */
          vlist[nadded++] = vt->v[m];
        }
      }

      free(defect->border);
      defect->border = (int *)malloc(nadded * sizeof(int));
      memmove(defect->border, vlist, nadded * sizeof(int));
      defect->nborder = nadded;
      /* set marks to zero */
      for (n = 0; n < defect->nborder; n++) {
        mris->vertices[defect->border[n]].marked = 0;
        ;
      }
      // if merged, we evaluate the defect again...
      i--;
    }
  }

  free(nd);

  /* clean defect from empty spaces */
  new_dl = (DEFECT_LIST *)calloc(1, sizeof(DEFECT_LIST));
  for (i = 0; i < ndefects; i++) {
    defect = &dl->defects[i];
    if (defect->nvertices) {
      /* this defect is not empty */
      memmove(&new_dl->defects[new_dl->ndefects++], defect, sizeof(DEFECT));
    }
  }

  free(dl);

  /* update defect statistics */
  for (i = 0; i < ndefects; i++) {
    defect = &new_dl->defects[i];
    if (defect->area >= 0) {
      continue;
    }
    defect->area = 0;
    defect->cx = 0;
    defect->cy = 0;
    defect->cz = 0;
    for (n = 0; n < defect->nvertices; n++) {
      VERTEX const * const v = &mris->vertices[defect->vertices[n]];
      defect->cx += v->x;
      defect->cy += v->y;
      defect->cz += v->z;
      defect->area += v->origarea;
    }

    defect->cx /= (float)defect->nvertices;
    defect->cy /= (float)defect->nvertices;
    defect->cz /= (float)defect->nvertices;

    defect->nx = defect->ny = defect->nz = 0;
    for (n = 0; n < defect->nborder; n++) {
      VERTEX const * const v = &mris->vertices[defect->border[n]];
      defect->nx += v->nx;
      defect->ny += v->ny;
      defect->nz += v->nz;
    }
    len = sqrt(defect->nx * defect->nx + defect->ny * defect->ny + defect->nz * defect->nz);
    if (FZERO(len)) {
      len = 1.0f;
    }
    defect->nx = defect->nx / len;
    defect->ny = defect->ny / len;
    defect->nz = defect->nz / len;
  }

  fprintf(stderr, "%d defects to be corrected \n", new_dl->ndefects);
  return new_dl;
}

static int defectIdentifyDefectiveVertices(MRI_SURFACE *mris,
                                           DEFECT *defect,
                                           FACE_DEFECT_LIST *fdl,
                                           float area_threshold,
                                           int mark_retained,
                                           int mark_discard,
                                           MHT *mht,
                                           int mode)
{
  int counting, n, p;
  FACE *f;
  VECTOR *v_a, *v_b, *v_n;
  float dot, area;

  v_a = VectorAlloc(3, MATRIX_REAL);
  v_b = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL); /* normal vector */

  /* set marks to zero */
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[n]];
    for (p = 0; p < vt->num; p++) {
      int const fno = vt->f[p];
      f = &mris->faces[fno];
      f->ripflag = 0;
      /* store information */
      setFaceOrigArea(mris, fno, f->area);
    }
  }

  /* store areas and compute spherical area */
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[n]];
    VERTEX          const * const v  = &mris->vertices         [defect->vertices[n]];
    for (p = 0; p < vt->num; p++) {
      f = &mris->faces[vt->f[p]];
      if (f->ripflag) {
        continue;
      }
      f->ripflag = 1;
      /* compute new area with new coord systems */
      VERTEX const * const v0 = &mris->vertices[f->v[0]];
      VERTEX const * const v1 = &mris->vertices[f->v[1]];
      VERTEX const * const v2 = &mris->vertices[f->v[2]];

      VERTEX_CANONICAL_EDGE(v_a, v0, v1);
      VERTEX_CANONICAL_EDGE(v_b, v0, v2);

      /* compute metric properties of first triangle */
      V3_CROSS_PRODUCT(v_a, v_b, v_n);
      area = V3_LEN(v_n) * 0.5f;
      dot = v->cx * V3_X(v_n) + v->cy * V3_Y(v_n) + v->cz * V3_Z(v_n);
      if (dot < 0.0f) /* not in same direction, area < 0 and reverse n */
      {
        f->area = -area;
      }
      else {
        f->area = area;
      }
    }
  }
  /* unrip */
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[n]];
    for (p = 0; p < vt->num; p++) {
      f = &mris->faces[vt->f[p]];
      f->ripflag = 0;
    }
  }

  for (n = 0; n < defect->nvertices; n++) {
    defect->status[n] = 0;
  }

  mrisMarkRetainedPartOfDefect(mris, defect, fdl, area_threshold, mark_retained, mark_discard, mht, mode);

  /* setting back old areas */
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[n]];
    for (p = 0; p < vt->num; p++) {
      FaceNormCacheEntry const * fNorm = getFaceNorm(mris, vt->f[p]);
      f = &mris->faces[vt->f[p]];
      f->area = fNorm->orig_area;
    }
  }

  for (counting = n = 0; n < defect->nvertices; n++)
    if (defect->status[n] == DISCARD_VERTEX) {
      counting++;
    }

  VectorFree(&v_a);
  VectorFree(&v_b);
  VectorFree(&v_n);

  return counting;
}
////////////////////////////////////////////////////////////////////////
//
//
//    This program will generate a topologically
//    correct surface mris_corrected from mris
//
//    canonical vertices = spherical vertices
//    original vertices = original vertices
//    tmp vertices = inflated vertices
//
//    CAREFUL: the fiedls v->curvbak, v->val, v->val2,
//             v->valbak are used by some functions
//
//    Defect output files (label, chull, borders) will be called
//      defectbase_{label,chull,borders}. If left NULL, defectbase
//      will be set to "defect"
//////////////////////////////////////////////////////////////////////

int topology_fixing_exit_after_diag = 0;


MRI_SURFACE *MRIScorrectTopology(
    MRI_SURFACE *mris, MRI *mri, MRI *mri_wm, 
    int nsmooth, TOPOLOGY_PARMS *parms, const char *defectbase)
{
  FACE_DEFECT_LIST *fdl;
  DEFECT_LIST *dl;
  DEFECT *defect;
  int fno, i, n, p, vno, kept_vertices, *face_trans, *vertex_trans, counter = 0, ninitialfaces;
  MHT *mht;
  FACE *f, *fdst;
  HISTOGRAM *h_k1, *h_k2, *h_gray, *h_white, *h_dot, *h_border, *h_grad;
  MRI *mri_gray_white, *mri_k1_k2;
  MRIS *mris_corrected_final;
  char tmpstr[2000];

  if(defectbase == NULL) defectbase = "defect";

#if ADD_EXTRA_VERTICES
  int retessellation_error = -1;
#endif

  OPTIMAL_DEFECT_MAPPING *o_d_m;

  fprintf(WHICH_OUTPUT, "\nCorrection of the Topology\n");

  //    canonical = spherical vertices
  //  original  = original vertices
  //  tmp       = inflated vertices
  //  current   = canonical vertices

  /* saving additional information */
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  // MRISwrite(mris,"orig_uncorrected");
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
// MRISwrite(mris,"inflated_uncorrected");
// current = inflated

#if MATRIX_ALLOCATION
  /* allocation of the transform matrix */
  VoxelFromSRASmatrix = GetSurfaceRASToVoxelMatrix(mri);
#endif

  /* centering the surface using CANONICAL_VERTICES */
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
  //    canonical = spherical vertices
  //  original  = original vertices
  //  tmp       = inflated vertices
  //  current   = spherical vertices
  MRIScenterSphere(mris);

  /* saving into CANONICAL_VERTICES */
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES);
  /* at this point : current = canonical vertices */

  /* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for now on, should avoid reprojecting vertices onto the sphere
     and even more recentering the canonical sphere
     XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  */

  /* marking intersecting edges */
  fdl = MRISmarkAmbiguousVertices(mris, MARK_AMBIGUOUS);
  ninitialfaces = mris->nfaces;
  // if (Gdiag & DIAG_SHOW)
  fprintf(WHICH_OUTPUT, "segmenting defects...\n");
  dl = MRISsegmentDefects(mris, MARK_AMBIGUOUS, MARK_SEGMENTED);

#if ADD_EXTRA_VERTICES
  /* analyze if the segmentation into connected defects is correct */
  dl = mrisDefectAnalysis(mris, dl);
  pushApartDefects(mris, dl);
#endif

  MRISsetVals(mris, 0.0f);
  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    defect->defect_number = i;
    for (n = 0; n < defect->nvertices; n++) {
      mris->vertices[defect->vertices[n]].val = defect->area;
    }
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRISwriteValues(mris, "defect_area");
  }
  // if (Gdiag & DIAG_SHOW)
  fprintf(WHICH_OUTPUT, "%d defects found, arbitrating ambiguous regions...\n", dl->ndefects);

  /* at this point : canonical vertices */
  MRIScomputeMetricProperties(mris);
#if 0
  /* should not modify anything about canonical vertices */
  mrisScaleMaxDimension(mris, FIELD_OF_VIEW*.9f) ;
#endif
  MRISclearMarks(mris);
  MRISclearCurvature(mris);

/* at this point : canonical vertices */
#if MERGE_NEIGHBORING_DEFECTS
  /* to avoid subtle topological errors
     became obsolete with ADD_EXTRA_VERTICES */
  dl = mrisMergeNeighboringDefects(mris, dl);
  MRISclearMarks(mris);
#endif

  mht = MHTcreateFaceTable(mris);

  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    defect->defect_number = i;

    /* Identify some defective vertices */
    mrisMarkRetainedPartOfDefect(mris, defect, fdl, AREA_THRESHOLD, MARK_RETAIN, MARK_DISCARD, mht, parms->keep);
    /* The convex hull is now constituted of the first neighbors only */
    mrisFindDefectConvexHull(mris, defect);
  }

  /* for diagnostic purposes */
  for (i = 0; i < dl->ndefects; i++) {
    int vno2, n2;

    defect = &dl->defects[i];
    defect->defect_number = i;
    for (n = 0; n < defect->nvertices + defect->nborder; n++) {
      if (n < defect->nvertices) {
        if (defect->status[n] == DISCARD_VERTEX) {
          continue;
        }
        vno = defect->vertices[n];
      }
      else {
        vno = defect->border[n - defect->nvertices];
      }
      VERTEX const * const v = &mris->vertices[vno];
      for (n2 = n + 1; n2 < defect->nvertices + defect->nborder; n2++) {
        if (n2 < defect->nvertices) {
          if (defect->status[n2] == DISCARD_VERTEX) {
            continue;
          }
          vno2 = defect->vertices[n2];
        }
        else {
          vno2 = defect->border[n2 - defect->nvertices];
        }
        if (vno == vno2) {
          continue;
        }
        VERTEX const * const vn = &mris->vertices[vno2];
        if (FEQUAL(vn->x, v->x) && FEQUAL(vn->y, v->y) && FEQUAL(vn->z, v->z)) {
          counter++;
          if (Gdiag & DIAG_SHOW) fprintf(WHICH_OUTPUT, "defect %d, vertices %d and %d coincident!\n", i, vno, vno2);
        }
      }
    }
  }
  fprintf(WHICH_OUTPUT, "%d vertices coincident\n", counter);
  if (Gdiag & DIAG_SHOW) {
    fprintf(WHICH_OUTPUT, "\n");
  }
  MHTfree(&mht);

  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  /* at this point : current = original vertices */

  /* v->val = border, v->val2 = white, v->val2bak = gray */
  mrisRipAllDefects(mris, dl, 1);
  mrisFindGrayWhiteBorderMean(mris, mri);
  mrisRipAllDefects(mris, dl, 0);

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  /* at this point : tmp becomes original vertices */
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  MRISaverageVertexPositions(mris, 2);
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
  /* at this point : original vertices become smoothed original vertices */
  // saving out additional information
  // MRISwrite(mris,"orig_smooth_uncorrected");

  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  /* at this point : back to original vertices */

  /* vertex information :
     canonical - canonical
     tmp - original
     original = smoothed
     current = original (non smoothed)
  */

  MRISsetNeighborhoodSizeAndDist(mris, 2);
  h_k1 = HISTOalloc(100);
  h_k2 = HISTOalloc(100);
  mri_k1_k2 = MRIalloc(100, 100, 1, MRI_FLOAT);
  h_dot = HISTOalloc(100);
  mrisComputePrincipalCurvatureDistributions(mris, h_k1, h_k2, mri_k1_k2);
  mrisComputeNormalDotDistribution(mris, h_dot);

  /* now knit each defect together by retessellating the surface,
     using the spherical space for topology (i.e. edge intersection),
     and the original space for geometry (i.e. edge length).
  */
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  /* at this point : 'non-smoothed' original vertices */

  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nvertices; n++) {
      mris->vertices[defect->vertices[n]].curv = defect->status[n] == DISCARD_VERTEX ? -1 : 1;
    }
  }

  if (((Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)) || (Gdiag & DIAG_SAVE_DIAGS)) {
    MRISwriteCurvature(mris, "defect_status");
  }

  /*  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;*/
  MRIScomputeMetricProperties(mris);
  MRISclearCurvature(mris);
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    for (i = 0; i < dl->ndefects; i++) {
      int total_defective_vertices;
      float total_defective_area;
      FILE *fp;
      char fname[STRLEN];

      int req = snprintf(fname, STRLEN, "%s.%s.defect%d.log",
			 mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data(), i);
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      fp = fopen(fname, "wb");
      fprintf(fp, "%d %2.3f\n", dl->defects[i].nvertices, dl->defects[i].area);
      for (n = 0; n < dl->defects[i].nvertices; n++) {
        fprintf(fp, "%d\n", dl->defects[i].vertices[n]);
        mris->vertices[dl->defects[i].vertices[n]].curv = (float)i + 1;
        if (dl->defects[i].vertices[n] == Gdiag_no) {
          DiagBreak();
        }
      }
      fprintf(fp, "\nborder (%d)\n", dl->defects[i].nborder);
      for (n = 0; n < dl->defects[i].nborder; n++) {
        fprintf(fp, "%d\n", dl->defects[i].border[n]);
        if (dl->defects[i].border[n] == Gdiag_no) {
          DiagBreak();
        }
      }

      fprintf(fp, "\nconvex hull (%d)\n", dl->defects[i].nchull);
      for (n = 0; n < dl->defects[i].nchull; n++) {
        fprintf(fp, "%d\n", dl->defects[i].chull[n]);
        if (dl->defects[i].chull[n] == Gdiag_no) {
          DiagBreak();
        }
      }
      fclose(fp);

      req = snprintf(fname, STRLEN, "%s.%s.defects.log",
		     mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data());
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      fp = fopen(fname, "wb");
      for (total_defective_area = 0.0f, total_defective_vertices = i = 0; i < dl->ndefects; i++) {
        total_defective_vertices += dl->defects[i].nvertices;
        total_defective_area += dl->defects[i].area;
      }
      fprintf(fp, "%d %2.3f\n", total_defective_vertices, total_defective_area);
      for (i = 0; i < dl->ndefects; i++) {
        for (n = 0; n < dl->defects[i].nvertices; n++) {
          fprintf(fp, "%d\n", dl->defects[i].vertices[n]);
        }
      }
      fclose(fp);
    }
  }

  // always writing out additional information
  // if(Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
      MRISwrite(mris, "new_orig_uncorrected");
      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
      MRISwrite(mris, "new_qsphere_uncorrected");
      MRISrestoreVertexPositions(mris, TMP_VERTICES);
    }
    MRISclearCurvature(mris);

    for (i = 0; i < dl->ndefects; i++) {
      defect = &dl->defects[i];
      for (n = 0; n < defect->nvertices; n++) {
        mris->vertices[defect->vertices[n]].curv = i + 1;
      }
    }
    // if (DIAG_VERBOSE_ON || (Gdiag & DIAG_SAVE_DIAGS))
    sprintf(tmpstr,"%s_labels",defectbase);
    MRISwriteCurvature(mris, tmpstr);
    for (i = 0; i < dl->ndefects; i++) {
      defect = &dl->defects[i];
      for (n = 0; n < defect->nborder; n++) {
        fflush(stdout);  // nicknote: prevents segfault on Linux PowerPC
        // when -O2 optimization is used w/gcc 3.3.3

        mris->vertices[defect->border[n]].curv = i + 1;
      }
    }
    // if (DIAG_VERBOSE_ON || (Gdiag & DIAG_SAVE_DIAGS))
    sprintf(tmpstr,"%s_borders",defectbase);
    MRISwriteCurvature(mris, tmpstr);
    for (i = 0; i < dl->ndefects; i++) {
      defect = &dl->defects[i];
      for (n = 0; n < defect->nchull; n++) {
        mris->vertices[defect->chull[n]].curv = 2;
      }
      for (n = 0; n < defect->nborder; n++) {
        mris->vertices[defect->border[n]].curv = -1;
      }
      for (n = 0; n < defect->nvertices; n++) {
        fflush(stdout);  // nicknote: prevents segfault on Linux PowerPC
        // when -O2 optimization is used w/gcc 3.3.3

        mris->vertices[defect->vertices[n]].curv = i + 1;
      }
    }
    // if (DIAG_VERBOSE_ON || (Gdiag & DIAG_SAVE_DIAGS))
    sprintf(tmpstr,"%s_chull",defectbase);
    MRISwriteCurvature(mris, tmpstr);
  }
  if (topology_fixing_exit_after_diag) {
    return (NULL);
  }

  /* now start building the target surface */
  MRISclearMarks(mris);
  kept_vertices = mris->nvertices;
  /* mark the defect with the value 1 (only inside vertices not border !) */
  for (i = 0; i < dl->ndefects; i++) {
    mrisMarkDefect(mris, &dl->defects[i], 1);
  }
  if (parms->search_mode == GREEDY_SEARCH) {
    /* don't keep all vertices */
    for (n = 0; n < dl->defects[i].nvertices; n++)
      if (dl->defects[i].status[n] == DISCARD_VERTEX) {
        kept_vertices--;
      }
  }

  face_trans   = (int *)calloc(mris->nfaces,    sizeof(int));
  vertex_trans = (int *)calloc(mris->nvertices, sizeof(int));
  memset(vertex_trans, -1, mris->nvertices    * sizeof(int));
  memset(face_trans,   -1, mris->nfaces       * sizeof(int));
  
  // create a new surface
  MRIS * mris_corrected = MRISoverAlloc(mris->nvertices + 10, 2 * mris->nfaces, kept_vertices, 2 * mris->nfaces);
  
  // keep the extra info into the new one
  mris_corrected->useRealRAS = mris->useRealRAS;
  mris_corrected->hemisphere = mris->hemisphere;
  copyVolGeom(&mris->vg, &mris_corrected->vg);

  mris_corrected->type = MRIS_TRIANGULAR_SURFACE;
  mris_corrected->status         = mris->status;            // this should have been done
  mris_corrected->origxyz_status = mris->origxyz_status;    // this is new 

#if 0
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
#endif
  MRISrestoreVertexPositions(mris, TMP_VERTICES); /* inflated */
  MRIScomputeMetricProperties(mris);

  int newNVertices = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    /* ignore the inside defect vertices but not the bordering ones */
    if (v->marked) {
      continue;
    }
    
    VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[newNVertices];
    VERTEX          * const vdst  = &mris_corrected->vertices         [newNVertices];
    if (mris_corrected->nvertices == Gdiag_no) {
      DiagBreak();
    }

    /* original vertices */
    vdst->x = v->x;
    vdst->y = v->y;
    vdst->z = v->z;
    
    /* smoothed vertices */
    MRISsetOriginalXYZ(mris_corrected, newNVertices, v->origx, v->origy, v->origz);
    
    vdst->tx = v->tx;
    vdst->ty = v->ty;
    vdst->tz = v->tz;
    vdst->nx = v->nx;
    vdst->ny = v->ny;
    vdst->nz = v->nz;
    
    /* canonical vertices */
    vdst->cx = v->cx;
    vdst->cy = v->cy;
    vdst->cz = v->cz;
    
    vdstt->num = vt->num;
    
    vdst->val = v->val;
    vdst->val2 = v->val2;
    vdst->valbak = v->valbak;
    vdst->val2bak = v->val2bak;
    vdst->imag_val = v->imag_val;
    vdst->curv = v->curv;
    vdst->curvbak = v->curvbak;
    vdst->stat = v->stat;
    vdst->mean = v->mean;
    vdst->mean_imag = v->mean_imag;
    vdst->std_error = v->std_error;
    vdst->H = v->H;
    vdst->K = v->K;
    vdst->k1 = v->k1;
    vdst->k2 = v->k2;
    vdst->border = 0;
    vertex_trans[vno] = newNVertices++;
  }
  
  /* now add all the retained vertices in the defects */
  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nvertices; n++) {
      if (defect->vertices[n] == Gdiag_no) {
        DiagBreak();
      }
      /* only add the kept vertices in greedy_search mode */
      if (parms->search_mode != GREEDY_SEARCH || defect->status[n] == KEEP_VERTEX) {
        vno = defect->vertices[n];
        VERTEX          const * const v  = &mris->vertices         [vno];
        if (vno == Gdiag_no) {
          DiagBreak();
        }
        
        int const vno_dst = newNVertices++;
        MRISgrowNVertices(mris_corrected, newNVertices);

        VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[vno_dst];
        VERTEX          * const vdst  = &mris_corrected->vertices         [vno_dst];
        if (mris_corrected->nvertices == Gdiag_no) {
          DiagBreak();
        }

        vdst->x = v->x;
        vdst->y = v->y;
        vdst->z = v->z;
        
        MRISsetOriginalXYZ(mris_corrected, vno_dst, v->origx, v->origy, v->origz);
        
        vdst->tx = v->tx;
        vdst->ty = v->ty;
        vdst->tz = v->tz;
        vdst->cx = v->cx;
        vdst->cy = v->cy;
        vdst->cz = v->cz;
        vdst->nx = v->nx;
        vdst->ny = v->ny;
        vdst->nz = v->nz;
        /* no num*/
        vdst->val = v->val;
        vdst->val2 = v->val2;
        vdst->valbak = v->valbak;
        vdst->val2bak = v->val2bak;
        vdst->imag_val = v->imag_val;
        vdst->curv = v->curv;
        vdst->curvbak = v->curvbak;
        vdst->stat = v->stat;
        vdst->mean = v->mean;
        vdst->mean_imag = v->mean_imag;
        vdst->std_error = v->std_error;
        vdst->H = v->H;
        vdst->K = v->K;
        vdst->k1 = v->k1;
        vdst->k2 = v->k2;
        vdstt->num = 0;
        clearVnum(mris_corrected,vno_dst);
        if (parms->search_mode != GREEDY_SEARCH && defect->status[n] == DISCARD_VERTEX) {
          vdst->ripflag = 1;
        }
        vertex_trans[vno] = vno_dst;
      }
    }
  }

  int newNfaces;
  for (newNfaces = fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    /* don't update triangle with marked vertices */
    if (triangleMarked(mris, fno)) {
      continue;
    }
    /* initialize face */
    fdst = &mris_corrected->faces[newNfaces];
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      fdst->v[n] = vertex_trans[f->v[n]];
    }
    face_trans[fno] = newNfaces++;
  }
  MRISgrowNFaces(mris_corrected, newNfaces);
  
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    FILE *fp;
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.%s.vtrans.log",
		       mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data()); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    fp = fopen(fname, "wb");
    if (!fp) {
      DiagBreak();
    }

    for (vno = 0; vno < mris->nvertices; vno++) {
      fprintf(fp, "%6d --> %6d\n", vno, vertex_trans[vno]);
    }
    fclose(fp);
    req = snprintf(fname, STRLEN, "%s.%s.ftrans.log",
		   mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data());
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fp = fopen(fname, "wb");

    for (vno = 0; vno < mris->nfaces; vno++) {
      fprintf(fp, "%6d --> %6d\n", vno, face_trans[vno]);
    }
    fclose(fp);
  }

  /* now allocate face and neighbor stuff in mris_corrected */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->marked) {
      continue;
    }
    if (vertex_trans[vno] < 0 || vertex_trans[vno] >= mris_corrected->nvertices) {
      continue;
    }
    
    int const vno_dst = vertex_trans[vno];
    VERTEX_TOPOLOGY * const vdstt = &mris_corrected->vertices_topology[vno_dst];

    /* count # of good triangles attached to this vertex */
    for (vdstt->num = n = 0; n < vt->num; n++)
      if (triangleMarked(mris, vt->f[n]) == 0) {
        vdstt->num++;
      }
      
    vdstt->f = (int *)  calloc(vdstt->num, sizeof(int));
    vdstt->n = (uchar *)calloc(vdstt->num, sizeof(uchar));
    for (i = n = 0; n < vt->num; n++) {
      if (triangleMarked(mris, vt->f[n])) {
        continue;
      }
      vdstt->n[i] = vt->n[n];
      vdstt->f[i] = face_trans[vt->f[n]];
      i++;
    }

    /* count # of valid neighbors */
    clearVnum(mris_corrected,vno_dst);
    for (n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].marked == 0) {
        addVnum(mris_corrected,vno_dst,1);
      }

    vdstt->v = (int *)calloc(vdstt->vnum, sizeof(int));
    for (i = n = 0; n < vt->vnum; n++)
      if (mris->vertices[vt->v[n]].marked == 0) {
        vdstt->v[i++] = vertex_trans[vt->v[n]];
      }
      
    vdstt->nsizeMax = 1;
    MRIS_setNsizeCur(mris_corrected, vno_dst, 1);
  }

  mrisCheckVertexFaceTopology(mris);
  
  MRISclearMarks(mris);
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  /* at this point : back to smoothed original vertices for mris */
  /*                 but still original (non smoothed)
                     vertices for mris_corrected */
  MRIScomputeMetricProperties(mris);
  MRISsmoothSurfaceNormals(mris, 10);

  //////////////// watch out ///////////////////////////////////////////
  h_gray = HISTOalloc(256);
  h_white = HISTOalloc(256);
  h_border = HISTOalloc(256);
  h_grad = HISTOalloc(256);
  mri_gray_white = MRIalloc(256, 256, 1, MRI_FLOAT);

  mrisMarkAllDefects(mris, dl, 1);
  // TO BE CHECKED IF MGH PROCESSED BRAIN (T1==110)
  mrisComputeJointGrayWhiteBorderDistributions(mris, mri, mri_gray_white, mri_wm);

  /* compute statistics on original */
  if (parms->search_mode != GREEDY_SEARCH)
    mrisComputeSurfaceStatistics(mris, mri, h_k1, h_k2, mri_k1_k2, mri_gray_white, h_dot);

  mrisMarkAllDefects(mris, dl, 0);
  for (i = 0; i < dl->ndefects; i++) {
    if (parms->correct_defect >= 0 && i != parms->correct_defect) {
      continue;
    }

    defect = &dl->defects[i];
    if (i == Gdiag_no) {
      DiagBreak();
    }
#if 0
    fprintf(WHICH_OUTPUT,
            "\rretessellating defect %d with %d vertices (chull=%d).    ",
            i, defect->nvertices+defect->nborder, defect->nchull) ;
#endif
    mrisMarkAllDefects(mris, dl, 1);
    mrisComputeGrayWhiteBorderDistributions(mris, mri, defect, h_white, h_gray, h_border, h_grad);
    mrisMarkAllDefects(mris, dl, 0);

#if 0
    MRIwrite(mri_gray_white,"mri_gw.mgh");
    MRIwrite(mri_k1_k2,"mri_k1_k2.mgh");
    HISTOplot(h_k1, "k1.plt") ;
    HISTOplot(h_k2, "k2.plt") ;
    HISTOplot(h_white, "hw.plt") ;
    HISTOplot(h_gray, "hg.plt") ;
#endif

#define TESTING_OPTIMAL 1

    if (parms->search_mode != GREEDY_SEARCH && parms->optimal_mapping) {
      // char fname[500];
      int validate, counting;
      MAPPING *mapping;
      FS_VERTEX_INFO *vinfo;
      DVS *dvs;
      // MRIS *mris_small;

      /* in optimal mode, check for self-intersection */
      parms->check_surface_intersection = 1;

      fprintf(stderr, "\nINITIAL MAPPING with %d vertices\n", defect->nvertices);

#if TESTING_OPTIMAL
      /* record vertex state before correction */
      dvs = mrisRecordVertexState(mris_corrected, defect, vertex_trans);

      /* first retessellate using spherical mapping */

      /* identify negative/defective vertices and ignore them */
      counting = defectIdentifyDefectiveVertices(
          mris, defect, fdl, AREA_THRESHOLD, MARK_RETAIN, MARK_DISCARD, mht, parms->keep);

      fprintf(stderr, "%d out of %d vertices were eliminated", counting, defect->nvertices);

      /* retessellate defect */
      mrisTessellateDefect(mris,
                           mris_corrected,
                           defect,
                           vertex_trans,
                           mri,
                           h_k1,
                           h_k2,
                           mri_k1_k2,
                           h_white,
                           h_gray,
                           h_border,
                           h_grad,
                           mri_gray_white,
                           h_dot,
                           parms);

#if 1
      {
        int ne, nv, nf, tt, theoric_euler, euler_nb;
        nf = mris_corrected->nfaces;
        ne = nv = 0;
        for (tt = 0; tt < mris_corrected->nvertices; tt++) {
          if (mris_corrected->vertices[tt].ripflag) {
            continue;
          }
          if (mris_corrected->vertices_topology[tt].vnum == 0) {
            continue;
          }
          ne += mris_corrected->vertices_topology[tt].vnum;
          nv++;
        }
        ne /= 2;
        euler_nb = nv + nf - ne;
        theoric_euler = 3 + defect->defect_number - dl->ndefects;
        fprintf(WHICH_OUTPUT,
                "After retessellation of defect %d, "
                "euler #=%d (%d,%d,%d) : difference with theory (%d) = %d \n",
                i,
                euler_nb,
                nv,
                ne,
                nf,
                theoric_euler,
                theoric_euler - euler_nb);

#if ADD_EXTRA_VERTICES
        if (theoric_euler - euler_nb && retessellation_error < 0) {
          retessellation_error = i;
        }
#endif
      }
#endif

      /* validation : no self-intersection and good retessellation */
      //... to be implemented
      validate = 0;

      if (!validate) {
        /* restore the vertex state */
        mrisRestoreVertexState(mris_corrected, dvs);
#else   // #if TESTING_OPTIMAL
      {
#endif  // #if TESTING_OPTIMAL
        /* generation of optimal mappings */
        defect->optimal_mapping = 1;

        /* modify chull : chull becomes useless */
        free(defect->chull);
        defect->nchull = defect->nborder;
        defect->chull = (int *)malloc(defect->nchull * sizeof(int));
        memmove(defect->chull, defect->border, defect->nchull * sizeof(int));

        /* save border positions */
        for (n = 0; n < defect->nborder; n++) {
          VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->border[n]]];
          vdst->t2x = vdst->cx;
          vdst->t2y = vdst->cy;
          vdst->t2z = vdst->cz;
        }
        /* save inside positions */
        for (n = 0; n < defect->nvertices; n++) {
          VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->vertices[n]]];
          vdst->t2x = vdst->cx;
          vdst->t2y = vdst->cy;
          vdst->t2z = vdst->cz;
        }

        /* generate different mappings (max 10) */
        o_d_m = mrisFindOptimalDefectMapping(mris, defect);

        for (p = 0; p < o_d_m->nmappings; p++) {
          fprintf(stderr, "MAPPING %d\n", p);

          /* select mapping */
          mapping = &o_d_m->mappings[p];

          /* use fixedval to identify good vertices */
          for (vno = 0; vno < mris_corrected->nvertices; vno++) {
            mris_corrected->vertices[vno].fixedval = 0;
          }

          /* use new coordinates */
          for (n = 0; n < defect->nchull; n++) {
            vinfo = &mapping->vertices[o_d_m->vertex_trans[defect->chull[n]]];
            VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->chull[n]]];
            VERTEX * const v = &mris->vertices[defect->chull[n]];
            vdst->cx = vinfo->c_x;
            vdst->cy = vinfo->c_y;
            vdst->cz = vinfo->c_z;
            v->cx = vinfo->c_x;
            v->cy = vinfo->c_y;
            v->cz = vinfo->c_z;
            vdst->ripflag = 0;
            vdst->fixedval = 1;
          }
          for (n = 0; n < defect->nvertices; n++) {
            vinfo = &mapping->vertices[o_d_m->vertex_trans[defect->vertices[n]]];
            VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->vertices[n]]];
            VERTEX * const v = &mris->vertices[defect->vertices[n]];
            vdst->cx = vinfo->c_x;
            vdst->cy = vinfo->c_y;
            vdst->cz = vinfo->c_z;
            v->cx = vinfo->c_x;
            v->cy = vinfo->c_y;
            v->cz = vinfo->c_z;
            vdst->ripflag = 0;
            vdst->fixedval = 1;
            defect->status[n] = 0;
          }
          /* identify negative/defective vertices and ignore them */
          counting = defectIdentifyDefectiveVertices(
              mris, defect, fdl, AREA_THRESHOLD, MARK_RETAIN, MARK_DISCARD, mht, parms->keep);
          fprintf(stderr, "%d out of %d vertices were eliminated", counting, defect->nvertices);

          /* retessellate defect */
          mrisTessellateDefect(mris,
                               mris_corrected,
                               defect,
                               vertex_trans,
                               mri,
                               h_k1,
                               h_k2,
                               mri_k1_k2,
                               h_white,
                               h_gray,
                               h_border,
                               h_grad,
                               mri_gray_white,
                               h_dot,
                               parms);

#if 1
          {
            int ne, nv, nf, tt, theoric_euler, euler_nb;
            nf = mris_corrected->nfaces;
            ne = nv = 0;
            for (tt = 0; tt < mris_corrected->nvertices; tt++) {
              if (mris_corrected->vertices[tt].ripflag) {
                continue;
              }
              if (mris_corrected->vertices_topology[tt].vnum == 0) {
                continue;
              }
              ne += mris_corrected->vertices_topology[tt].vnum;
              nv++;
            }
            ne /= 2;
            euler_nb = nv + nf - ne;
            theoric_euler = 3 + defect->defect_number - dl->ndefects;
            fprintf(WHICH_OUTPUT,
                    "After retessellation of defect %d (v0 =  %d), euler #=%d (%d,%d,%d) : difference with theory (%d) "
                    "= %d \n",
                    i,
                    defect->vertices[0],
                    euler_nb,
                    nv,
                    ne,
                    nf,
                    theoric_euler,
                    theoric_euler - euler_nb);

#if ADD_EXTRA_VERTICES
            if (theoric_euler - euler_nb && retessellation_error < 0) {
              retessellation_error = i;
            }
#endif
          }
#endif
#if 0
          /* save solution */
          sprintf(fname,"./lh.final_%d",p);
          fprintf(stderr,"writting solution into %s",fname);
          MRISrestoreVertexPositions(mris_corrected,ORIGINAL_VERTICES);
#if 0

          MRISwrite(mris_corrected,fname);
#else
          mris_small=extractDefect(mris_corrected,defect);
          MRISwrite(mris_small,fname);
          MRISfree(&mris_small);
#endif
          MRISrestoreVertexPositions(mris_corrected,TMP_VERTICES);
#endif

          /* restore border positions */
          for (n = 0; n < defect->nborder; n++) {
            VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->border[n]]];
            vdst->cx = vdst->t2x;
            vdst->cy = vdst->t2y;
            vdst->cz = vdst->t2z;
          }
          /* save inside positions */
          for (n = 0; n < defect->nvertices; n++) {
            VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->vertices[n]]];
            vdst->cx = vdst->t2x;
            vdst->cy = vdst->t2y;
            vdst->cz = vdst->t2z;
          }
          break;

          fprintf(stderr, "restore\n");
          /* restore the vertex state */
          mrisRestoreVertexState(mris_corrected, dvs);
        }
        /* free the structure o_d_m */
        MRISfree(&o_d_m->mris);
        free(o_d_m->vertex_trans);
        free(o_d_m->face_trans);
        free(o_d_m->orig_mapping.vertices);
        for (n = 0; n < 10; n++) {
          free(o_d_m->mappings[n].vertices);
        }
        free(o_d_m);
      }
/* free the structure dvs */
#if TESTING_OPTIMAL
      mrisFreeDefectVertexState(dvs);
#endif
    }
    else {
      // main part of the routine: retessellation of the defect
      mrisTessellateDefect(mris,
                           mris_corrected,
                           defect,
                           vertex_trans,
                           mri,
                           h_k1,
                           h_k2,
                           mri_k1_k2,
                           h_white,
                           h_gray,
                           h_border,
                           h_grad,
                           mri_gray_white,
                           h_dot,
                           parms);
    }

    /* compute Euler number of surface */
    if (parms->search_mode != GREEDY_SEARCH) {
      int ne, nv, nf, tt, theoric_euler, euler_nb;
      nf = mris_corrected->nfaces;
      ne = nv = 0;
      for (tt = 0; tt < mris_corrected->nvertices; tt++) {
        if (mris_corrected->vertices[tt].ripflag) {
          continue;
        }
        if (mris_corrected->vertices_topology[tt].vnum == 0) {
          continue;
        }
        ne += mris_corrected->vertices_topology[tt].vnum;
        nv++;
      }
      ne /= 2;
      euler_nb = nv + nf - ne;
      theoric_euler = 3 + defect->defect_number - dl->ndefects;
      fprintf(WHICH_OUTPUT,
              "After retessellation of defect %d (v0=%d), "
              "euler #=%d (%d,%d,%d) : "
              "difference with theory (%d) = %d \n",
              i,
              defect->vertices[0],
              euler_nb,
              nv,
              ne,
              nf,
              theoric_euler,
              theoric_euler - euler_nb);
#if ADD_EXTRA_VERTICES
      if (theoric_euler - euler_nb && retessellation_error < 0) {
        retessellation_error = i;
      }
#endif
    }

    if (parms->correct_defect >= 0 && i == parms->correct_defect)
      ErrorExit(ERROR_BADPARM, "TERMINATING PROGRAM AFTER CORRECTED DEFECT\n");
  }
#if ADD_EXTRA_VERTICES
  if (retessellation_error >= 0) {
    fprintf(WHICH_OUTPUT,
            "\nXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\nThe first retessellation error happened for defect %d\n\n",
            retessellation_error);
  }
#endif

  if (Gdiag & DIAG_SAVE_DIAGS) {
    char fname[STRLEN], path[STRLEN];

    MRISclearCurvature(mris_corrected);
    for (i = 0; i < dl->ndefects; i++) {
      defect = &dl->defects[i];
      for (n = 0; n < defect->nvertices; n++) {
        VERTEX * const vdst = &mris_corrected->vertices[vertex_trans[defect->vertices[n]]];
        if (vdst->ripflag == 0) {
          vdst->curv = (i + 1);  // for diagnostics
        }
      }
    }
    FileNamePath(mris->fname, path);
    int req = snprintf(fname, STRLEN, "%s/%s.fixed.defect_labels.mgz",
		       path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    printf("writing corrected defect labels to %s\n", fname);
    MRISwriteCurvature(mris_corrected, fname);
  }
  HISTOfree(&h_white);
  HISTOfree(&h_gray);
  HISTOfree(&h_border);
  HISTOfree(&h_dot);
  HISTOfree(&h_k1);
  HISTOfree(&h_k2);
  HISTOfree(&h_grad);
  MRIfree(&mri_gray_white);
  MRIfree(&mri_k1_k2);

  if (parms->search_mode == GREEDY_SEARCH) {
    mrisAddAllDefectFaces(mris_corrected, dl, vertex_trans);
  }

  mrisCheckSurface(mris_corrected);

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    MHT *mht;

    fprintf(WHICH_OUTPUT, "checking corrected surface for self-intersection...\n");
    MRISsaveVertexPositions(mris_corrected, TMP_VERTICES);
    MRISrestoreVertexPositions(mris_corrected, ORIG_VERTICES);
    mht = MHTcreateFaceTable(mris_corrected);
    MHTfree(&mht);
    MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES);
  }

#if 0
  for (i = 0 ; i < dl->ndefects ; i++)
  {
    defect = &dl->defects[i] ;
    for (n = 0 ; n < defect->nvertices ; n++)
    {
      vno = vertex_trans[defect->vertices[n]] ;
      if (vno == Gdiag_no)
      {
        DiagBreak() ;
      }
      if (vno < 0)
      {
        continue ;
      }
      v = &mris_corrected->vertices[vno] ;
      if (v->vnum < 2)
      {
        fprintf(WHICH_OUTPUT,
                "Warning: vertex %d has only %d neighbors!\n",
                vno, v->vnum) ;
        DiagBreak() ;
      }
    }
  }
  if (Gdiag & DIAG_SHOW)
  {
    fprintf(WHICH_OUTPUT, "\n") ;
  }
#endif

  mrisForgetNeighborhoods(mris_corrected);
  
  mrisCheckVertexFaceTopology(mris_corrected);

  fprintf(WHICH_OUTPUT, "computing original vertex metric properties...\n");
  MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);
  
  if (0) {
    const char * fnm = "./after_MRISrestoreVertexPositions.log";
    fprintf(stdout, "%s:%d Dump shape to %s\n",__FILE__,__LINE__,fnm);
    FILE* file = fopen(fnm,"w");
    mrisDumpShape(file, mris_corrected);
    fclose(file);
  }
  
  /* at this point : smoothed corrected orig vertices */
  MRIScomputeMetricProperties(mris_corrected);
  fprintf(WHICH_OUTPUT, "storing new metric properties...\n");
  /*  MRISstoreMetricProperties(mris_corrected) ;*/
  fprintf(WHICH_OUTPUT, "computing tessellation statistics...\n");
  MRISprintTessellationStats(mris_corrected, stderr);

  // mark everything not in a defect with 1
  MRISsetMarks(mris_corrected, 1);
  for (i = 0; i < dl->ndefects; i++) {
    defect = &dl->defects[i];
    for (n = 0; n < defect->nvertices; n++) {
      vno = vertex_trans[defect->vertices[n]];
      if (vno < 0 || vno >= mris_corrected->nvertices) {
        continue;
      }
      VERTEX * const v = &mris_corrected->vertices[vno];
      v->marked = 0;
    }
    for (n = 0; n < defect->nborder; n++) {
      vno = vertex_trans[defect->border[n]];
      if (vno < 0 || vno >= mris_corrected->nvertices) {
        continue;
      }
      VERTEX * const v = &mris_corrected->vertices[vno];
      v->marked = 0;
    }
  }
  fprintf(WHICH_OUTPUT,
          "performing soap bubble on retessellated vertices for %d "
          "iterations...\n",
          nsmooth);
  /* at this point : smoothed corrected orig vertices */

  MRISsoapBubbleVertexPositions(mris_corrected, nsmooth);
  MRISsaveVertexPositions(mris_corrected, ORIGINAL_VERTICES);
  MRISclearMarks(mris_corrected);

  MRISprintTessellationStats(mris_corrected, stderr);
  MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES);
  /* at this point : back to original vertices */
  fprintf(WHICH_OUTPUT, "tessellation finished, orienting corrected surface...\n");

  if (parms->search_mode == GREEDY_SEARCH) {
    mrisOrientRetessellatedSurface(mris_corrected, dl, vertex_trans);
  }

  if (parms->save_fname) {
    char fname[500];
    if (MRISmarkOrientationChanges(mris_corrected)) {
      MRISclearCurvature(mris);
      for (i = 0; i < mris->nvertices; i++) {
        vno = vertex_trans[i];
        if (vno >= 0 && vno < mris_corrected->nvertices) {
          mris->vertices[i].curv = mris_corrected->vertices[vno].curv;
        }
      }
      sprintf(fname, "%s/orientation_changes", parms->save_fname);
      MRISwriteCurvature(mris, fname);
      MRISrestoreVertexPositions(mris_corrected, CANONICAL_VERTICES);
      sprintf(fname, "%s/new_sphere", parms->save_fname);
      MRISwrite(mris_corrected, fname);
      sprintf(fname, "%s/orientation_changes2", parms->save_fname);
      MRISwriteCurvature(mris_corrected, fname);
      MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);
      MRISclearCurvature(mris_corrected);
      for (i = 0; i < dl->ndefects; i++) {
        defect = &dl->defects[i];
        for (n = 0; n < defect->nborder; n++) {
          vno = vertex_trans[defect->border[n]];
          if (vno < 0 || vno >= mris_corrected->nvertices) {
            continue;
          }
          VERTEX * const v = &mris_corrected->vertices[vno];
          v->curv = -1;
        }
      }
      sprintf(fname, "%s/borders", parms->save_fname);
      MRISwriteCurvature(mris_corrected, fname);
    }
  }

  free(face_trans);
  free(vertex_trans);
  /* free structures */
  for (fno = 0; fno < ninitialfaces; fno++)
    if (fdl->nfaces[fno] > 0) {
      free(fdl->faces[fno]);
    }

  free(fdl->faces);
  free(fdl->nfaces);
  free(fdl);
  for (i = 0; i < dl->ndefects; i++) {
    if (dl->defects[i].vertices) {
      free(dl->defects[i].vertices);
    }
    if (dl->defects[i].status) {
      free(dl->defects[i].status);
    }
    if (dl->defects[i].border) {
      free(dl->defects[i].border);
    }
    if (dl->defects[i].edges) {
      free(dl->defects[i].edges);
    }
  }
  free(dl);

#if MATRIX_ALLOCATION
  if (VoxelFromSRASmatrix) {
    MatrixFree(&VoxelFromSRASmatrix);
  }
#endif

  if (nmut + ncross > 0)
    fprintf(WHICH_OUTPUT,
            "%ld mutations (%2.1f%%), %ld crossovers (%2.1f%%), "
            "%ld vertices were eliminated\n",
            nmut,
            (float)nmut * 100 / (nmut + ncross),
            ncross,
            (float)ncross * 100 / (nmut + ncross),
            nkilled);
  else {
    fprintf(WHICH_OUTPUT, "%ld vertices were eliminated\n", nkilled);
  }

  mris_corrected_final = MRISremoveRippedSurfaceElements(mris_corrected);
  mris_corrected_final->hemisphere = mris->hemisphere;
  strcpy(mris_corrected_final->fname, mris->fname);

  MRISfree(&mris_corrected);
  /* current = orig vertices
     tmp = originial vertices
     orig = smoothed correct vertices = true solution
     canonical = canonical vertices
  */

  return (mris_corrected_final);
} // finished MRIScorrectTopology()

static int mrisMarkAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag)
{
  int j;

  for (j = 0; j < dl->ndefects; j++) {
    mrisMarkDefect(mris, &dl->defects[j], flag);
  }
  return (NO_ERROR);
}

static int mrisRipAllDefects(MRI_SURFACE *mris, DEFECT_LIST *dl, int flag)
{
  int j;

  for (j = 0; j < dl->ndefects; j++) {
    mrisRipDefect(mris, &dl->defects[j], flag);
  }
  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Mark all the vertices in the tessellation that are part
  of faces whose centroid intersects other faces. These are
  regions in which the topology is broken as the spherical
  homeomorphism is non-invertible.
  ------------------------------------------------------*/
FACE_DEFECT_LIST *MRISmarkAmbiguousVertices(MRI_SURFACE *mris, int mark)
{
  FACE *f;
  VERTEX *v;
  int fno, flist[MAX_INT_FACES], i, nfaces, nmarked, n /*, vno, neg*/;
  double area_scale;
#if 0
  double r;
#endif
  MHT *mht;
  FILE *fp = NULL;
  FDL *fdl;

  fdl = (FACE_DEFECT_LIST *)calloc(1, sizeof(FDL));
  if (!fdl) ErrorExit(ERROR_NO_MEMORY, "MRISmarkAmbiguousFaces: could allocate face defect list");
  fdl->nfaces = (int *)calloc(mris->nfaces, sizeof(int));
  if (!fdl->nfaces) ErrorExit(ERROR_NO_MEMORY, "MRISmarkAmbiguousFaces: could allocate face defect list");
  fdl->faces = (int **)calloc(mris->nfaces, sizeof(int *));
  if (!fdl->faces) ErrorExit(ERROR_NO_MEMORY, "MRISmarkAmbiguousFaces: could allocate face defect list");
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.%s.topology.log",
		       mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data());
    if( req >= STRLEN ) {
       std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fp = fopen(fname, "w");
  }

  /* the curvature is for diagnostic purposes so I can write it out */
  MRISclearMarks(mris);
  MRISclearCurvature(mris);
  mrisMarkBadEdgeVertices(mris, mark);

#if 0
  /*
    should not reproject vertices on to sphere
    should not recenter the canonical sphere !
  */
  r = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris, mris, 100.0/r) ;
#endif

  mht = MHTcreateFaceTable(mris);

  area_scale = mris->orig_area / mris->total_area;
  // if (Gdiag & DIAG_SHOW)
  fprintf(WHICH_OUTPUT, "marking ambiguous vertices...\n");

  for (nmarked = fno = 0; fno < mris->nfaces; fno++) {
    if (Gdiag & DIAG_SHOW && !(fno % 25000) && fno)
      fprintf(WHICH_OUTPUT, "%d of %d faces processed, %d ambiguous\n", fno, mris->nfaces - 1, nmarked);
    f = &mris->faces[fno];
    if (fno == Gdiag_no) {
      DiagBreak();
    }
    if (f->ripflag) {
      continue;
    }

    /* only edge-intersection to identify overlapping faces */
    nfaces = mrisFindAllOverlappingFaces(mris, mht, fno, flist);

    /* make sure fno is in list, and add it if it isn't (it should be) */
    for (i = 0; i < nfaces; i++)
      if (flist[i] == fno) {
        break;
      }
    if (i >= nfaces) {
      if (nfaces == 1000) {
        ErrorExit(ERROR_BADPARM, "Too many faces");
      }
      flist[nfaces++] = fno;
    }
#if 1
    if (nfaces > 1)
#else
    if ((nfaces > 1 || area_scale * f->area < 0.001) || ((fno <= 5) && Gdiag & DIAG_SAVE_DIAGS)) /* part of a defect */
#endif
    {
      nmarked++;
#if 0
      if (Gdiag & DIAG_SHOW)
        fprintf(WHICH_OUTPUT, "\r%d of %d faces processed, %d ambiguous",
                fno, mris->nfaces-1, nmarked) ;
#endif

      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        fprintf(WHICH_OUTPUT, "%d faces @ fno %d\n", nfaces, fno);
        for (i = 0; i < nfaces; i++) {
          f = &mris->faces[flist[i]];
          fprintf(WHICH_OUTPUT, "\tface %d area %2.4f (%d, %d, %d)\n", flist[i], f->area, f->v[0], f->v[1], f->v[2]);
        }
        fprintf(WHICH_OUTPUT, "\n");
      }
      if (Gdiag & DIAG_WRITE && fp != NULL && DIAG_VERBOSE_ON) {
        fprintf(fp, "%d faces @ fno %d\n", nfaces, fno);
        for (i = 0; i < nfaces; i++) {
          f = &mris->faces[flist[i]];
          fprintf(fp, "\tface %d area %2.4f (%d, %d, %d)\n", flist[i], f->area, f->v[0], f->v[1], f->v[2]);
        }
        fprintf(fp, "\n");
        fflush(fp);
      }

      fdl->nfaces[fno] = nfaces;
      fdl->faces[fno] = (int *)calloc(nfaces, sizeof(int));
      if (!fdl->faces[fno]) ErrorExit(ERROR_NO_MEMORY, "MRISmarkAmbiguousFaces: could allocate %d defect list", fno);
      for (i = 0; i < nfaces; i++) {
        fdl->faces[fno][i] = flist[i];
        f = &mris->faces[flist[i]];
        for (n = 0; n < VERTICES_PER_FACE; n++) {
          v = &mris->vertices[f->v[n]];
          if (f->v[n] == Gdiag_no) {
            DiagBreak();
          }
          v->marked = mark;
        }
      }
    }
  }
// TO BE CHECKED : ... pbm when "#if 1...
#if 0
  /* expand defective vertices outwards by one */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == mark)
    {
      for (n = 0 ; n < v->vnum ; n++)
      {
        mris->vertices[v->v[n]].marked = mark+1 ;
      }
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == mark+1)
    {
      v->marked = mark ;
    }
  }
#endif

  if (Gdiag & DIAG_SHOW) {
    fprintf(WHICH_OUTPUT, "\n");
  }
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON && fp) {
    fclose(fp);
  }

#if 0
  /* should not reproject vertices on to sphere */
  MRISscaleBrain(mris, mris, r/100.0) ;
#endif

  fprintf(WHICH_OUTPUT, "%d ambiguous faces found in tessellation\n", nmarked);
  MHTfree(&mht);
  return (fdl);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Find all connected components of a defect.
  use undefval,old_undefval, curv and fixedval during segmentation
  into connected components:
  - fixedval is used to mark border vertices
  - curv is used to save the marked flags during divideedge.
  - undefval and old_undefval are used to find the enclosing loop
  ------------------------------------------------------*/
DEFECT_LIST *MRISsegmentDefects(MRI_SURFACE *mris, int mark_ambiguous, int mark_segmented)
{
  DEFECT_LIST *dl;
  int vno, nadded;
  VERTEX *v;
  FILE *fp = NULL;
  DEFECT *defect;

  dl = (DEFECT_LIST *)calloc(1, sizeof(DEFECT_LIST));
  if (!dl) ErrorExit(ERROR_NO_MEMORY, "MRISsegmentDefects: could allocate defect list");

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.%s.topology.log",
		       mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", mris->subject_name.data());
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fp = fopen(fname, "a");
  }

#if ADD_EXTRA_VERTICES /* uses fixedval to mark border vertices */
  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].fixedval = 0;
  }
#endif

  for (nadded = vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->marked != mark_ambiguous) {
      continue;
    }
    if (dl->ndefects == Gdiag_no) {
      DiagBreak();
    }

    defect = &dl->defects[dl->ndefects];
    defect->defect_number = dl->ndefects++;
    //                fprintf(stderr,"DEFECT %d\n",defect->defect_number);

    /* segment defect #defect->defect_number */
    nadded += mrisSegmentDefect(mris, vno, defect, mark_ambiguous, mark_segmented);

#if FIND_ENCLOSING_LOOP
    /* update the defect so it becomes simply connected */
    mrisSimplyConnectedDefect(mris, defect, mark_ambiguous, mark_segmented);
#endif

    if (Gdiag & DIAG_WRITE && fp) {
      int n;
      DEFECT *defect = &dl->defects[dl->ndefects - 1];
      fprintf(fp,
              "defect %d found with %d vertices, area %2.2f\n"
              "\tcentroid (%2.1f,%2.1f,%2.1f)\n",
              dl->ndefects,
              defect->nvertices,
              defect->area,
              defect->cx,
              defect->cy,
              defect->cz);
      for (n = 0; n < defect->nvertices; n++) {
        fprintf(fp, "\t%d\n", defect->vertices[n]);
      }
    }
  }
  if (nadded) fprintf(stderr, "   total of %d vertices have been added to the surface\n", nadded);
  // if (Gdiag & DIAG_SHOW) fprintf(WHICH_OUTPUT, "\n") ;
  if (Gdiag & DIAG_WRITE && fp) {
    fclose(fp);
  }
  return (dl);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Segment the connected region of a defect, starting with vno
  and spreading outwards.
  ------------------------------------------------------*/
static int mrisSegmentDefect(MRI_SURFACE *mris, int vno, DEFECT *defect, int mark_ambiguous, int mark_segmented)
{
  int vlist[MAX_DEFECT_VERTICES], i, j, n, nfilled, nadded, vno1, m;

  float len, nx, ny, nz;

  vno1 = nadded = m = j = 0; /* to avoid compilator warnings */

  if (defect->nvertices + 1 >= MAX_DEFECT_VERTICES)
    ErrorExit(ERROR_NOMEMORY, "mrisSegmentDefect: max number of defective vertices %d exceeded\n", MAX_DEFECT_VERTICES);
  vlist[defect->nvertices++] = vno; /* start the list */

  {
    VERTEX * const v = &mris->vertices[vno];
    v->marked = mark_segmented;
    defect->cx = v->x;
    defect->cy = v->y;
    defect->cz = v->z;
    defect->area = v->origarea;
  }
  
  do {
    nfilled = 0;

    for (i = 0; i < defect->nvertices; i++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vlist[i]];
      for (n = 0; n < vt->vnum; n++) {
        VERTEX * const vn = &mris->vertices[vt->v[n]];
        if (vt->v[n] == Gdiag_no) {
          DiagBreak();
        }
        if (vn->marked == mark_ambiguous) {
          if (defect->nvertices + 1 >= MAX_DEFECT_VERTICES)
            ErrorExit(ERROR_NOMEMORY,
                      "mrisSegmentDefect: max number of defective vertices %d exceeded\n",
                      MAX_DEFECT_VERTICES);
          vlist[defect->nvertices++] = vt->v[n]; /* add it to list */

          vn->marked = mark_segmented;
          defect->cx += vn->x;
          defect->cy += vn->y;
          defect->cz += vn->z;
          defect->area += vn->origarea;
          nfilled++;
        }
      }
    }
  } while (nfilled > 0);

  defect->cx /= (float)defect->nvertices;
  defect->cy /= (float)defect->nvertices;
  defect->cz /= (float)defect->nvertices;
  defect->vertices = (int *)calloc(defect->nvertices, sizeof(int));
  if (!defect->vertices) ErrorExit(ERROR_NO_MEMORY, "mrisSegmentDefect: could allocate defect vertex list");
  defect->status = (char *)calloc(defect->nvertices, sizeof(char));
  if (!defect->status) ErrorExit(ERROR_NO_MEMORY, "mrisSegmentDefect: could allocate defect status list");
  memmove(defect->vertices, vlist, defect->nvertices * sizeof(int));

/* analyze if some extra vertices should be added to the surface */
#if ADD_EXTRA_VERTICES
  for (nadded = i = 0; i < defect->nvertices; i++) {
    v = &mris->vertices[defect->vertices[i]];
    for (n = 0; n < v->vnum; n++) {
      vno1 = v->v[n];
      vn = &mris->vertices[vno1];
      if (vn->marked == 0) {
        /* border vertex */
        /* check if this border vertex is already part of one defect */
        if (vn->fixedval) {
          /* needs to add one extra vertex */
          if (mris->nvertices == mris->max_vertices) {
            ErrorExit(ERROR_BADPARM, "ADD_EXTRA_VERTICES : could not allocate extra vertex\n");
          }
          /* hack to avoid errors : using curv to save marked vertices */
          for (j = 0; j < mris->nvertices; j++) {
            VERTEX *vj;
            vj = &mris->vertices[j];
            vj->curv = vj->marked;
            vj->marked = 0;
          }
          mrisDivideEdge(mris, defect->vertices[i], vno1);
          for (j = 0; j < mris->nvertices - 1; j++) {
            /* reset marks */
            VERTEX *vj;
            vj = &mris->vertices[j];
            vj->marked = (int)vj->curv;
          }

          VERTEX * const vadded = &mris->vertices[mris->nvertices - 1];
          // fprintf(stderr,"adding vertex %d(%d)\n",
          // mris->nvertices-1,vn->fixedval);
          vadded->marked = 0; /* border vertex  */
          vadded->fixedval = 0;

          /* spherical projection */
          sphericalProjection(vadded->cx, vadded->cy, vadded->cz, &vadded->cx, &vadded->cy, &vadded->cz);
          vadded->x = vadded->cx;
          vadded->y = vadded->cy;
          vadded->z = vadded->cz;
          nadded++;
        }
      }
    }
  }
  if (nadded)
    fprintf(stderr, "   defect %d : %d vertices have been added to the surface\n", defect->defect_number, nadded);

#endif

  for (nfilled = i = 0; i < defect->nvertices; i++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[i]];
    VERTEX                * const v  = &mris->vertices         [defect->vertices[i]];
    if (defect->vertices[i] == Gdiag_no) {
      DiagBreak();
    }
    v->val = defect->area;
    defect->status[i] = KEEP_VERTEX;
    for (n = 0; n < vt->vnum; n++) {
      if (vt->v[n] == Gdiag_no) {
        DiagBreak();
      }
      if (mris->vertices[vt->v[n]].marked == 0) /* border vertex */
      {
        mris->vertices[vt->v[n]].marked = 2;
        vlist[nfilled++] = vt->v[n];
      }
    }
  }

  defect->border = (int *)calloc(nfilled, sizeof(int));
  defect->nborder = nfilled;
  memmove(defect->border, vlist, defect->nborder * sizeof(int));
  mrisMarkDefectBorder(mris, defect, 0);

  nx = ny = nz = 0.0f;
  for (n = 0; n < defect->nborder; n++) {
    VERTEX * const v = &mris->vertices[defect->border[n]];
    nx += v->nx;
    ny += v->ny;
    nz += v->nz;
  }
  len = sqrt(nx * nx + ny * ny + nz * nz);
  if (FZERO(len)) {
    len = 1.0f;
  }
  defect->nx = nx / len;
  defect->ny = ny / len;
  defect->nz = nz / len;
  return (nadded);
}
//#if FIND_ENCLOSING_LOOP
/* segment a surface into connected components using undefval to mark
   the different components and avoiding old_undefval */
static int mrisSegmentConnectedComponents(MRIS *mris)
{
  int n, p, vno = 0, seed, found;
  int *vlist, nvertices;
  int *next_vlist, next_nvertices, total_vertices, maxn, maxv;
  VERTEX *v, *vp;

  vlist = (int *)malloc(mris->nvertices * sizeof(int));
  next_vlist = (int *)malloc(mris->nvertices * sizeof(int));

  for (n = 0; n < mris->nvertices; n++) {
    mris->vertices[n].undefval = 0;
  }

  seed = 1; /* first seed */
  maxv = 0;
  maxn = 0;
  while (1) {
    found = 0;
    for (n = 0; n < mris->nvertices; n++) {
      v = &mris->vertices[n];
      if (v->old_undefval) {
        continue;
      }
      if (v->undefval == 0) {
        v->undefval = seed;
        vno = n;
        found = 1;
        break;
      }
    }
    if (found == 0) {
      break;
    }

    /* grow seed point */
    vlist[0] = vno;
    nvertices = 1;
    total_vertices = 1;

    while (nvertices) {
      for (next_nvertices = n = 0; n < nvertices; n++) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vlist[n]];

        for (p = 0; p < vt->vnum; p++) {
          vp = &mris->vertices[vt->v[p]];
          if (vp->old_undefval) {
            continue;
          }
          if (vp->undefval) {
            continue;
          }

          /* new point to be added */
          next_vlist[next_nvertices++] = vt->v[p];
          vp->undefval = seed;
        }
      }
      nvertices = next_nvertices;
      memmove(vlist, next_vlist, mris->nvertices * sizeof(int));
      total_vertices += nvertices;
    }
    if (maxv < total_vertices) {
      maxn = seed;
      maxv = total_vertices;
    }
    // fprintf(stderr,"seed %d has %d vertices \n",seed,total_vertices);

    seed++;
  }
  free(vlist);
  free(next_vlist);

  if (maxn == 0) {
    ErrorExit(ERROR_BADPARM, "0 labels found in mrisSegmentConnectedComponents\n");
  }
  return maxn;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Making sure the defect is simply connected
  region growing method (the largest component is the good one)
  use undefval,old_undefval during the segmentation process
  - old_undefval are forbidden points (borders, inside)
  - undefval are flags for the segmentation into connected components
  ------------------------------------------------------*/
static int mrisSimplyConnectedDefect(MRI_SURFACE *mris, DEFECT *defect, int mark_ambiguous, int mark_segmented)
{
  int n, p, l, label, j, w;
  int nvertices, nedges, inside_face, outside_face;
  int add_vertex, add_edges;
  FACE *f;
  EDGE *edges, *edge, *new_edges;
  int added_edges, vno1, vno2, vno_in, vno_out, isedge;
  int *varray, *barray, vnb, bnb;
  float len, dot, cx, cy, cz;
  VECTOR *v_a, *v_b, *v_n;

  // fprintf(stderr,"\nDEFECT %d\n",defect->defect_number);
  // fprintf(stderr,"before analysis: %d inside vertices, %d border vertices
  // (total=%d)\n",defect->nvertices,defect->nborder,defect->nvertices+defect->nborder);

  v_a = VectorAlloc(3, MATRIX_REAL);
  v_b = VectorAlloc(3, MATRIX_REAL);
  v_n = VectorAlloc(3, MATRIX_REAL);

  /* setting flags to zero */
  for (n = 0; n < mris->nvertices; n++) {
    mris->vertices[n].undefval = 0;
    mris->vertices[n].old_undefval = 0;
  }

  for (n = 0; n < defect->nborder; n++) {
    VERTEX * const v = &mris->vertices[defect->border[n]];
    v->old_undefval = 1; /* forbidden border point */
  }
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX * const v = &mris->vertices[defect->vertices[n]];
    v->old_undefval = 1; /* forbidden inside point */
  }

  /* segment the surface into connected components
     label is the largest component */
  label = mrisSegmentConnectedComponents(mris);

  /* everything not label becomes defective (old_undefval=1)*/
  /* first count them */
  vnb = 0;
  for (n = 0; n < mris->nvertices; n++) {
    /* these first vertices are new ones */
    VERTEX const * const v = &mris->vertices[n];
    if (v->old_undefval) {
      continue;
    }
    if (v->undefval != label) {
      vnb++; /* only counting for now */
    }
  }
  /* allocate the new list of defective vertices
     max size is vnb+defect->nvertices+defect->nborder */
  varray = (int *)malloc((vnb + defect->nvertices + defect->nborder) * sizeof(int));
  /* then add them */
  vnb = 0;
  for (n = 0; n < mris->nvertices; n++) {
    /* this first vertices are new ones */
    VERTEX * const v = &mris->vertices[n];
    if (v->old_undefval) {
      continue;
    }
    if (v->undefval != label) {
      varray[vnb++] = n;
      v->undefval = 1;
    }
    else {
      v->undefval = 0;
    }
  }

  for (n = 0; n < defect->nvertices; n++) {
    VERTEX * const v = &mris->vertices[defect->vertices[n]];
    v->old_undefval = 0;
    v->undefval = 1;
  }

  for (n = 0; n < defect->nborder; n++) {
    /* border vertices become 2 */
    VERTEX * const v = &mris->vertices[defect->border[n]];
    v->old_undefval = 0;
    v->undefval = 2;
  }

  /* max number of edges */
  nedges = defect->nborder * (defect->nborder - 1) / 2;
  edges = (EDGE *)malloc(nedges * sizeof(EDGE));
  added_edges = 0;

  /* now need to find the right border vertices to keep (border = 2)
     right vertices = forms good edge in between inside (1 or 2) and outside 0
     this subtlety is important : the closed loop must be "convexified"
     the right border vertices are marked as 3 */

  for (nvertices = nedges = n = 0; n < defect->nborder; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->border[n]];
    VERTEX                * const v  = &mris->vertices         [defect->border[n]];
    add_vertex = add_edges = 0;
    /* count # of neighbors */
    for (p = 0; p < vt->vnum; p++) {
      if (mris->vertices[vt->v[p]].undefval == 2 || mris->vertices[vt->v[p]].undefval == 3) {
        /* potential edge */
        /* is there an inside and outside triangle in common with this vertex */
        inside_face = 0;
        outside_face = 0;
        vno_in = vno_out = -1;
        for (l = 0; l < vt->num; l++) {
          f = &mris->faces[vt->f[l]];
          /* check if this face has vt->v[p] */
          for (j = 0; j < 3; j++) {
            if (f->v[j] == vt->v[p]) {
              /* check if the last vertex is outside */
              for (w = 0; w < 3; w++) {
                if (mris->vertices[f->v[w]].undefval == 0) {
                  /* outside */
                  outside_face = 1;
                  vno_out = f->v[w];
                  break;
                }
              }
              for (w = 0; w < 3; w++) {
                if (mris->vertices[f->v[w]].undefval == 1 /* inside vertex */
                    || ((mris->vertices[f->v[w]].undefval == 2 || mris->vertices[f->v[w]].undefval == 3) &&
                        f->v[w] != defect->border[n] && f->v[w] != vt->v[p])) {
                  /* border vertex */
                  inside_face = 1;
                  vno_in = f->v[w];
                  break;
                }
              }
            }
          }
        }
        if (inside_face && outside_face) {
          /* good edge -> good vertex */
          add_edges++;
          add_vertex = 1;
          v->undefval = 3;

          /* vertices constituting the edge to be potentially added */
          vno1 = defect->border[n];
          vno2 = vt->v[p];

          /* check is edge already exists */
          for (isedge = w = 0; w < added_edges; w++) {
            edge = &edges[w];
            if ((edge->vno1 == vno1 && edge->vno2 == vno2) || (edge->vno1 == vno2 && edge->vno2 == vno1)) {
              isedge = 1;
              break;
            }
          }
          if (isedge == 0) {
            /* add edge to list */
            /* orient edge first */
            VERTEX          const * const v1 = &mris->vertices[vno1];
            VERTEX          const * const v2 = &mris->vertices[vno2];
            VERTEX                * const vout = &mris->vertices[vno_out];
            VERTEX          const * const vin  = &mris->vertices[vno_in];
            VECTOR_LOAD(v_a, v2->cx - v1->cx, v2->cy - v1->cy, v2->cz - v1->cz);
            VECTOR_LOAD(v_b, vin->cx - vout->cx, vin->cy - vout->cy, vin->cz - vout->cz);
            V3_CROSS_PRODUCT(v_a, v_b, v_n);
            cx = v1->cx + v2->cx;
            cy = v1->cy + v2->cy;
            cz = v1->cz + v2->cz;
            dot = cx * V3_X(v_n) + cy * V3_Y(v_n) + cz * V3_Z(v_n);

            if (dot < 0.0f) {
              /* not in same direction, reverse vno1, vno2 */
              w = vno1;
              vno1 = vno2;
              vno2 = w;
            }

            /* add edge */
            edge = &edges[added_edges++];
            edge->vno1 = vno1;
            edge->vno2 = vno2;
          }
        }
      }
    }
    nedges += add_edges;
    if (add_vertex) {
      nvertices++;
    }
  }
  nedges /= 2;

  VectorFree(&v_a);
  VectorFree(&v_b);
  VectorFree(&v_n);

  //    fprintf(stderr,"Defect %d : Euler Number = %d [ %d , %d  ]
  //    \n",defect->defect_number,nvertices-nedges,nvertices,nedges);

  if (nvertices - nedges) {
    ErrorExit(ERROR_BADPARM, "mrisSimplyConnectedDefect : euler number of loop is not 0!\n");
  }

  /* now update defect */
  /* first list good enclosing edges */
  defect->edges = (EDGE *)malloc(added_edges * sizeof(EDGE));
  defect->nedges = added_edges;
  memmove(defect->edges, edges, added_edges * sizeof(EDGE));
  free(edges);

  /* then keep updating inside and border vertices */
  bnb = 0;
  barray = (int *)malloc((defect->nvertices + defect->nborder) * sizeof(int));
  for (n = 0; n < defect->nvertices; n++) {
    varray[vnb++] = defect->vertices[n];
  }
  for (n = 0; n < defect->nborder; n++) {
    VERTEX const * const v = &mris->vertices[defect->border[n]];
    if (v->undefval == 3) /* good border vertex */
    {
      barray[bnb++] = defect->border[n];
    }
    else /* becomes inside vertex */
    {
      varray[vnb++] = defect->border[n];
    }
  }
  free(defect->vertices);
  free(defect->status);
  free(defect->border);
  defect->nvertices = vnb;
  defect->vertices = (int *)malloc(vnb * sizeof(int));
  defect->status = (char *)calloc(vnb, sizeof(char));
  memmove(defect->vertices, varray, vnb * sizeof(int));
  defect->nborder = bnb;
  defect->border = (int *)malloc(bnb * sizeof(int));
  memmove(defect->border, barray, bnb * sizeof(int));
  free(varray);
  free(barray);

  /* analyze the closed loop and correctly order the edges */
  {
    int init, final, cur, next;
    int *tab, nb = 0;
    tab = (int *)malloc(defect->nedges * sizeof(int));

    new_edges = (EDGE *)malloc(defect->nedges * sizeof(EDGE));

    // fprintf(stderr,"analyze closed loop\n");

    init = defect->edges[0].vno2;
    final = defect->edges[0].vno1;
    /* find other vertex */
    cur = -1;
    for (n = 1; n < defect->nedges; n++) {
      edge = &defect->edges[n];
      if (edge->vno1 == init) {
        // fprintf(stderr,"*");
        cur = edge->vno2;
        break;
      }
      if (edge->vno2 == init) {
        cur = edge->vno1;
        // fprintf(stderr,"!");
        break;
      }
    }
    //          fprintf(stderr,"%d->%d->%d",final,init,cur);
    if (cur < 0) {
      ErrorExit(ERROR_BADPARM, "mrisSimplyConnectedDefect : initialization of the closed loop\n");
    }
    tab[nb++] = init;
    tab[nb++] = cur;
    add_edges = 0;
    new_edges[add_edges].vno1 = final;
    new_edges[add_edges++].vno2 = init;
    new_edges[add_edges].vno1 = init;
    new_edges[add_edges++].vno2 = cur;
    while (cur != final) {
      /*find next vertex */
      next = -1;
      for (n = 0; n < defect->nedges; n++) {
        edge = &defect->edges[n];
        if (edge->vno1 == cur && edge->vno2 != init) {
          // fprintf(stderr,"*");
          next = edge->vno2;
          break;
        }
        if (edge->vno2 == cur && edge->vno1 != init) {
          // fprintf(stderr,".");
          next = edge->vno1;
          break;
        }
      }
      //                        fprintf(stderr,"->%d",next);
      if (next < 0) {
        ErrorExit(ERROR_BADPARM, "mrisSimplyConnectedDefect : propagation of the closed loop\n");
      }
      /* make sure next is not in the loop already */
      for (n = 0; n < nb; n++) {
        if (tab[n] == next) {
          ErrorExit(ERROR_BADPARM, "mrisSimplyConnectedDefect : loop self-intersecting\n");
        }
      }
      tab[nb++] = next;

      init = cur;
      cur = next;
      new_edges[add_edges].vno1 = init;
      new_edges[add_edges++].vno2 = cur;
    }
    free(tab);

    free(defect->edges);
    defect->edges = new_edges;

    for (n = 0; n < defect->nedges; n++) {
      edge = &defect->edges[n];
      // fprintf(stderr,"%d->%d:",edge->vno1,edge->vno2);
    }

    //          fprintf(stderr,"loop of %d=%d vertices\n",nb,defect->nborder);
    if (nb != defect->nborder) {
      ErrorExit(ERROR_BADPARM, "mrisSimplyConnectedDefect : loop smaller than border\n");
    }
  }
  // fprintf(stderr,"\n");

  /* now update defect statistics */
  defect->cx = 0;
  defect->cy = 0;
  defect->cz = 0;
  defect->area = 0;
  for (n = 0; n < defect->nvertices; n++) {
    VERTEX * const v = &mris->vertices[defect->vertices[n]];
    defect->status[n] = KEEP_VERTEX;
    v->marked = mark_segmented;
    defect->cx += v->x;
    defect->cy += v->y;
    defect->cz += v->z;
    defect->area += v->origarea;
  }
  defect->cx /= (float)defect->nvertices;
  defect->cy /= (float)defect->nvertices;
  defect->cz /= (float)defect->nvertices;

  defect->nx = defect->ny = defect->nz = 0.0f;
  for (n = 0; n < defect->nborder; n++) {
    VERTEX * const v = &mris->vertices[defect->border[n]];
    v->marked = 0;
    defect->nx += v->nx;
    defect->ny += v->ny;
    defect->nz += v->nz;
  }
  len = sqrt(defect->nx * defect->nx + defect->ny * defect->ny + defect->nz * defect->nz);
  if (FZERO(len)) {
    len = 1.0f;
  }
  defect->nx /= len;
  defect->ny /= len;
  defect->nz /= len;

// fprintf(stderr,"after analysis: %d inside vertices, %d border vertices
// (total=%d)\n",defect->nvertices,defect->nborder,defect->nvertices+defect->nborder);

#if ADD_EXTRA_VERTICES
  /* mark the border vertices */
  for (n = 0; n < defect->nborder; n++) {
    VERTEX * const v = &mris->vertices[defect->border[n]];
    v->fixedval = defect->defect_number + 1;
  }
#endif

  return NO_ERROR;
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Mark the vertices in a defect as either retained or
  discarded. The decision will be made based on the surface
  area in the defect. If it is above some threshold, then
  fill it by marking the "outside" faces as kept and all others
  as discarded. Otherwise cut it by marking the "inside" face
  as kept and all others as discarded. The inside/outside
  decision will be made by using the dot product of the average
  inflated surface normal with the inflated face centroid.

  By default (genetic search), MRIScorrectTopology will now keep most vertices
  ------------------------------------------------------*/
#define MIN_SPHERE_DIST .01
#define MIN_ORIG_DIST .75

#if 1
int mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris,
                                 DEFECT *defect,
                                 FACE_DEFECT_LIST *fdl,
                                 float area_threshold,
                                 int mark_retain,
                                 int mark_discard,
                                 MHT *mht,
                                 int mode)
{
#if 0
  int      n, i, j, nfaces, fno, flist[100000], n2, vno, retain;
  FACE     *f ;
  VERTEX   *v, *vn ;
  float    dot, x0, y0, z0, x, y, z, dx, dy, dz, dist, fn, dot0, len ;
#endif

  mrisMarkDefect(mris, defect, 0);

  if (!mode) {
    mrisDefectRemoveDegenerateVertices(mris, MIN_SPHERE_DIST, defect);
    mrisDefectRemoveProximalVertices(mris, MIN_ORIG_DIST, defect);
    mrisDefectRemoveNegativeVertices(mris, defect);
  }
  else {
    int i;
    for (i = 0; i < defect->nvertices; i++) {
      defect->status[i] = KEEP_VERTEX;
    }
    mrisDefectRemoveDegenerateVertices(mris, MIN_SPHERE_DIST, defect);
  }

#if 0
  /* throw out anything in a negative face */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
    {
      continue ;
    }
    v = &mris->vertices[defect->vertices[i]] ;
    for (n = 0 ; n < v->num ; n++)
      if (mris->faces[v->f[n]].area < 0.05)
      {
        defect->status[i] = DISCARD_VERTEX ;
      }
  }

  /* compute centroid and average normal of defect using border vertices */
  defect->cx = defect->cy = defect->cz = 0.0f ;
  defect->nx = defect->ny = defect->nz = 0.0f ;
  for (fn = 0.0f, i = 0 ; i < defect->nborder ; i++, fn += 1.0f)
  {
    v = &mris->vertices[defect->border[i]] ;
    defect->nx += v->nx ;
    defect->ny += v->ny ;
    defect->nz += v->nz ;
    defect->cx += v->x  ;
    defect->cy += v->y  ;
    defect->cz += v->z ;
  }
  len = sqrt(SQR(defect->nx) + SQR(defect->ny) + SQR(defect->nz)) ;
  defect->nx /= len ;
  defect->ny /= len ;
  defect->nz /= len ;
  defect->cx /= fn  ;
  defect->cy /= fn  ;
  defect->cz /= fn  ;

  /* discard vertices that are too close to another vertex */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
      {
        continue ;
      }
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
    {
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    }
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
      {
        continue ;
      }
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->origx-v->origx ;
      dy = vn->origy-v->origy ;
      dz = vn->origz-v->origz ;
      dist = sqrt(dx*dx+dy*dy+dz*dz) ;
      if (dist <= 0.75)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(WHICH_OUTPUT, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }

  /* build a list of all faces in this defect */
  for (nfaces = n = 0 ; n < defect->nvertices ; n++)
  {
    if (defect->status[n] == DISCARD_VERTEX)
    {
      continue ;
    }

    vno = defect->vertices[n] ;
    v = &mris->vertices[vno] ;

    /* build a list of faces */
    for (n2 = 0 ; n2 < v->num ; n2++)
    {
      if (mris->faces[v->f[n2]].ripflag == 0)
      {
        if (nfaces == 7127)
        {
          DiagBreak() ;
        }
        if (nfaces == 100000)
        {
          ErrorExit(ERROR_BADPARM, "Too many faces");
        }
        flist[nfaces++] = v->f[n2] ;
        mris->faces[v->f[n2]].ripflag = 1 ;  /* temporary */
      }
    }
  }

  /* for really big defects throw out 'inside' vertices */
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ;
    f = &mris->faces[fno] ;
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0) ;
    dx = x0 - defect->cx ;
    dy = y0 - defect->cy ;
    dz = z0 - defect->cz ;
    dot0 = dx*defect->nx + dy*defect->ny + dz*defect->nz ;

    /* see if there are any faces inside (outside) of this one */
    retain = 1 ;
    for (n2 = 0 ; n2 < fdl->nfaces[fno] ; n2++)
    {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1)
      {
        continue ;
      }
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = x - defect->cx ;
      dy = y - defect->cy ;
      dz = z - defect->cz ;
      dot = dx*defect->nx + dy*defect->ny + dz*defect->nz ;
#define HUGE_DEFECT 10000
#define BIG_DEFECT 5000
      if ((defect->nvertices > HUGE_DEFECT) && (dot > dot0))
      {
        retain = 0 ;
        break ;   /* found a face outside of this one - discard it */
      }
      if (defect->nvertices > BIG_DEFECT &&
          defect->nvertices < HUGE_DEFECT)
      {
        if (dot < dot0)  /* found a face inside this one - keep it */
        {
          retain = 1 ;
          break ;
        }
        else
        {
          retain =  0 ;
        }
      }
    }
    if (!retain)  /* no faces outside of this one */
    {
      for (n2 = 0 ; n2 < VERTICES_PER_FACE ; n2++)
      {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no)
        {
          DiagBreak() ;
        }
        mris->vertices[f->v[n2]].marked = 1 ;
        mris->vertices[f->v[n2]].imag_val = -1 ;
      }
    }
  }

  /* discard all marked vertices */
  for (i = 0 ; i < defect->nvertices ; i++)
  {
    if (defect->status[i] == DISCARD_VERTEX)
    {
      continue ;
    }
    v = &mris->vertices[defect->vertices[i]] ;
    if (v->marked)
    {
      v->marked = 0 ;
      defect->status[i] = DISCARD_VERTEX ;
    }
  }

  /* unmark the faces */
  for (n = 0 ; n < nfaces ; n++)
  {
    fno = flist[n] ;
    f = &mris->faces[fno] ;
    f->ripflag = 0 ;
  }

  /* discard vertices that are too close to another vertex */
  for (i = 0 ; i < defect->nvertices+defect->nborder ; i++)
  {
    float  dx, dy, dz ;

    if (i < defect->nvertices)
    {
      if (defect->status[i] == DISCARD_VERTEX)
      {
        continue ;
      }
      v = &mris->vertices[defect->vertices[i]] ;
    }
    else
    {
      v = &mris->vertices[defect->border[i-defect->nvertices]] ;
    }
    for (j = i+1 ; j < defect->nvertices ; j++)
    {
      if (defect->status[j] == DISCARD_VERTEX)
      {
        continue ;
      }
      vn = &mris->vertices[defect->vertices[j]] ;
      dx = vn->cx-v->cx ;
      dy = vn->cy-v->cy ;
      dz = vn->cz-v->cz ;
      dist = (dx*dx+dy*dy+dz*dz) ;  /* no sqrt */
      if (dist < MIN_SPHERE_DIST)
      {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(WHICH_OUTPUT, "discarding proximal vertex %d\n",
                  defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX ;
        vn->imag_val = -1.0f ;
      }
    }
  }
#endif
  return (NO_ERROR);
}
#else
static int mrisMarkRetainedPartOfDefect(MRI_SURFACE *mris,
                                        DEFECT *defect,
                                        FACE_DEFECT_LIST *fdl,
                                        float area_threshold,
                                        int mark_retain,
                                        int mark_discard,
                                        MHT *mht)
{
#define USING_CUBE 0
#if USING_CUBE || 1
  /* based on faces */
  int n, n2, inside, vno, flist[100000], nfaces, fno, i, j;
  VERTEX *v, *vn;
  FACE *f;
  float dot, x0, y0, z0, x, y, z, dx, dy, dz, fn, len, dist;

  defect->cx = defect->cy = defect->cz = 0.0f;
  inside = defect->area < area_threshold;
  mrisMarkDefect(mris, defect, 1);
  for (fn = 0.0, nfaces = n = 0; n < defect->nvertices; n++) {
    vno = defect->vertices[n];
    v = &mris->vertices[vno];

    /* build a list of faces */
    for (n2 = 0; n2 < v->num; n2++) {
      if (mris->faces[v->f[n2]].ripflag == 0) {
        if (nfaces == 100000) {
          ErrroExit(ERROR_BADPARM, "Too many faces");
        }
        flist[nfaces++] = v->f[n2];
        mris->faces[v->f[n2]].ripflag = 1; /* temporary */
      }
    }
    v->val = inside ? -1.0f : 1.0;
    v->imag_val = 1.0f; /* assume kept until found otherwise */
    defect->cx += v->x;
    defect->cy += v->y;
    defect->cz += v->z;
  }
  defect->nx = defect->ny = defect->nz = 0.0f;
  for (n = 0; n < defect->nborder; n++) {
    vno = defect->border[n];
    v = &mris->vertices[vno];
    defect->nx += v->nx;
    defect->ny += v->ny;
    defect->nz += v->nz;
  }
  mrisMarkDefect(mris, defect, 0);
  len = sqrt(defect->nx * defect->nx + defect->ny * defect->ny + defect->nz * defect->nz);
  if (FZERO(len)) {
    len = 1.0f;
  }
  defect->nx /= len;
  defect->ny /= len;
  defect->nz /= len;
  fn = (float)defect->nvertices;
  defect->cx /= fn;
  defect->cy /= fn;
  defect->cz /= fn;

  /* unrip the faces (used really as a mark, but don't want to add to
  face struct).
  */
  for (n = 0; n < defect->nvertices; n++) {
    v = &mris->vertices[defect->vertices[n]];

    for (n2 = 0; n2 < v->num; n2++) {
      mris->faces[v->f[n2]].ripflag = 0;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(WHICH_OUTPUT, "%d faces found in defect\n", nfaces);
  }
#if USING_CUBE
  for (n = 0; n < nfaces; n++) {
    fno = flist[n];
    f = &mris->faces[fno];
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0);

    /* see if there are any faces inside (outside) of this one */
    for (dot = 0.0f, n2 = 0; n2 < fdl->nfaces[fno]; n2++) {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1) {
        continue;
      }
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = -15.5 - x;
      dot = dx * dx;
      if (dot < ((-15.5 - x0) * (-15.5 - x0))) {
        break;
      }
    }
    if (n2 < fdl->nfaces[fno]) /* found a face inside
     (outside) of this one */
    {
      for (n2 = 0; n2 < VERTICES_PER_FACE; n2++) {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no) {
          DiagBreak();
        }
        mris->vertices[f->v[n2]].marked = 1;
        mris->vertices[f->v[n2]].imag_val = -1;
      }
    }
  }
#else
  for (n = 0; n < nfaces; n++) {
    fno = flist[n];
    f = &mris->faces[fno];
    mrisCalculateFaceCentroid(mris, fno, &x0, &y0, &z0);

    /* see if there are any faces inside (outside) of this one */
    for (dot = 0.0f, n2 = 0; n2 < fdl->nfaces[fno]; n2++) {
      if (triangleNeighbors(mris, fno, fdl->faces[fno][n2]) >= 1) {
        continue;
      }
      mrisCalculateFaceCentroid(mris, fdl->faces[fno][n2], &x, &y, &z);
      dx = x - x0;
      dy = y - y0;
      dz = z - z0;
      dot = dx * defect->nx + dy * defect->ny + dz * defect->nz;
      if ((inside && dot < 0) || (!inside && dot > 0)) {
        break;
      }
    }
    if (n2 < fdl->nfaces[fno]) /* found a face inside
   (outside) of this one */
    {
      for (n2 = 0; n2 < VERTICES_PER_FACE; n2++) {
        if (f->v[n2] == 1245 || f->v[n2] == Gdiag_no) {
          DiagBreak();
        }
        mris->vertices[f->v[n2]].marked = 1;
        mris->vertices[f->v[n2]].imag_val = -1;
      }
    }
  }
#endif
  for (n = 0; n < defect->nvertices; n++) {
    vno = defect->vertices[n];
    v = &mris->vertices[vno];
    defect->status[n] = v->marked ? DISCARD_VERTEX : KEEP_VERTEX;
  }
  mrisMarkDefect(mris, defect, 0);
#else
  int n, n2, inside, vno;
  VERTEX *v, *v2;
  float fn, len;

  defect->cx = defect->cy = defect->cz = defect->nx = defect->ny = defect->nz = 0.0;
  inside = defect->area < area_threshold;
  mrisMarkDefect(mris, defect, 1);
  for (fn = 0.0, n = 0; n < defect->nvertices; n++) {
    vno = defect->vertices[n];
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (inside) {
      v->val = -1.0f;
    }
    else {
      v->val = 1.0f;
    }
    defect->cx += v->x;
    defect->cy += v->y;
    defect->cz += v->z;

    /* use surrounding unmarked vertices as estimate of local normal */
    for (n2 = 0; n2 < v->vnum; n2++) {
      v2 = &mris->vertices[v->v[n2]];
      if (!v2->marked) {
        defect->nx += v2->nx;
        defect->ny += v2->ny;
        defect->nz += v2->nz;
        fn += 1.0;
      }
    }
  }
  mrisMarkDefect(mris, defect, 0);
  defect->nx /= fn;
  defect->ny /= fn;
  defect->nz /= fn;
  len = sqrt(defect->nx * defect->nx + defect->ny * defect->ny + defect->nz * defect->nz);
  if (FZERO(len)) {
    len = 1.0f;
  }
  defect->nx /= len;
  defect->ny /= len;
  defect->nz /= len;
  fn = (float)defect->nvertices;
  defect->cx /= fn;
  defect->cy /= fn;
  defect->cz /= fn;

  for (n = 0; n < defect->nvertices; n++) /* for each vertex in defect */
  {
    vno = defect->vertices[n];
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    defect->status[n] = KEEP_VERTEX;
    v->imag_val = 1.0;
    v->nx = defect->nx;
    v->ny = defect->ny;
    v->nz = defect->nz;
    if (inside) {
      if (mrisFindNextInwardFace(mris, mht, vno, 20.0f) >= 0) {
        defect->status[n] = DISCARD_VERTEX;
        defect->status[n] = KEEP_VERTEX;
        v->imag_val = -1.0;
      }
    }
    else {
      if (mrisFindNextOutwardFace(mris, mht, vno, 20.0f) >= 0) {
        defect->status[n] = DISCARD_VERTEX;
        v->imag_val = -1.0;
      }
    }
  }
#endif

  /* throw out anything in a negative face */
  for (i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    v = &mris->vertices[defect->vertices[i]];
    for (n = 0; n < v->num; n++)
      if (mris->faces[v->f[n]].area < 0) {
        defect->status[i] = DISCARD_VERTEX;
      }
  }

  /* discard vertices that are too close to another vertex */
  for (i = 0; i < defect->nvertices + defect->nborder; i++) {
    float dx, dy, dz;

    if (i < defect->nvertices) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      v = &mris->vertices[defect->vertices[i]];
    }
    else {
      v = &mris->vertices[defect->border[i - defect->nvertices]];
    }
    for (j = i + 1; j < defect->nvertices; j++) {
      if (defect->status[j] == DISCARD_VERTEX) {
        continue;
      }
      vn = &mris->vertices[defect->vertices[j]];
      dx = vn->origx - v->origx;
      dy = vn->origy - v->origy;
      dz = vn->origz - v->origz;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      if (dist <= 0.5) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(WHICH_OUTPUT, "discarding proximal vertex %d\n", defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX;
        vn->imag_val = -1.0f;
      }
    }
  }
  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0; i < defect->nvertices + defect->nborder; i++) {
    float dx, dy, dz;

    if (i < defect->nvertices) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      v = &mris->vertices[defect->vertices[i]];
    }
    else {
      v = &mris->vertices[defect->border[i - defect->nvertices]];
    }
    for (j = i + 1; j < defect->nvertices; j++) {
      if (defect->status[j] == DISCARD_VERTEX) {
        continue;
      }
      vn = &mris->vertices[defect->vertices[j]];
      dx = vn->cx - v->cx;
      dy = vn->cy - v->cy;
      dz = vn->cz - v->cz;
      dist = (dx * dx + dy * dy + dz * dz); /* no sqrt */
      if (dist < MIN_SPHERE_DIST) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(WHICH_OUTPUT, "discarding proximal vertex %d\n", defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX;
        vn->imag_val = -1.0f;
      }
    }
  }
  return (NO_ERROR);
}
#endif

static int mrisRipDefect(MRI_SURFACE *mris, DEFECT *defect, int ripflag)
{
  int n;

  for (n = 0; n < defect->nvertices; n++) {
    mris->vertices[defect->vertices[n]].ripflag = ripflag;
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Mark all the vertices in the given defect.
  ------------------------------------------------------*/
static int mrisMarkDefect(MRI_SURFACE *mris, DEFECT *defect, int mark)
{
  int n;

  for (n = 0; n < defect->nvertices; n++) {
    mris->vertices[defect->vertices[n]].marked = mark;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Mark all the border vertices in the given defect.
  ------------------------------------------------------*/
static int mrisMarkDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect, int mark)
{
  int n;

  for (n = 0; n < defect->nchull; n++) {
    mris->vertices[defect->chull[n]].marked = mark;
  }

  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Mark all the border vertices in the given defect.
  ------------------------------------------------------*/
static int mrisMarkDefectBorder(MRI_SURFACE *mris, DEFECT *defect, int mark)
{
  int n;

  for (n = 0; n < defect->nborder; n++) {
    mris->vertices[defect->border[n]].marked = mark;
  }

  return (NO_ERROR);
}


// compute a signed distance volume from two surfaces ( the original ones + a corrected defect )
void MRIScomputeDistanceVolume(TOPOFIX_PARMS *parms, float distance_to_surface)
{
  int k, i, j, p;
  float x0, x1, x2, y0, y1, y2, z0, z1, z2;
  float x, y, z;

  int fn1;
  int delta;
  int imin, imax, jmin, jmax, kmin, kmax;
  float distance, sign, scale;
  float n_f[3], n_e0[3], n_e1[3], n_e2[3], n_v0[3], n_v1[3], n_v2[3];
  float vec[3], vec0[3], vec1[3], vec2[3], e0[3], e1[3], e2[3], n0[3], n1[3], n2[3];
  float val, valu, val0, val1, val2;

  MRIP *mrip;
  MRIS *mris_defect, *mris_source, *mris;
  MRI *mri_defect, *mri_distance;
  int wsurf;
  int n, n_faces, n_vertices;
  int *ffrom, *fto, *vfrom, *vto;

  int found;

  mrip = parms->mrip;
  mris_defect = parms->mris_defect;  // the defect surface
  mris_source = mrip->mris_source;   // the source surface

  // fprintf(stderr,"INFO:{%d %d} {%d (%d) %d (%d)}
  // \n",mris_defect->nfaces,mris_defect->nvertices,mris_source->nfaces,mris_source->max_faces,mris_source->nvertices,mris_source->max_vertices);

  n_faces = mrip->n_faces;
  n_vertices = mrip->n_vertices;
  ffrom = mrip->ftrans_from;
  fto = mrip->ftrans_to;
  vfrom = mrip->vtrans_from;
  vto = mrip->vtrans_to;

  mri_defect = ((DP *)parms->dp)->mri_defect;
  mri_distance = ((DP *)parms->dp)->mri_defect_sign;

  // marking defect faces of source surface, i.e. mris_source
  // MRISclearFaceMarks(mris_source); //useless
  for (n = 0; n < n_faces; n++) {
    mris_source->faces[fto[n]].marked = 1;
  }
  // marking defect vertices of source surface, i.e. mris_source
  // MRISclearMarks(mris_source); //useless
  for (n = 0; n < n_vertices; n++) {
    mris_source->vertices[vto[n]].marked = 1;
  }
  // unmarking defect faces for defect surface, i.e. mris_defect
  MRISclearFaceMarks(mris_defect);

  /* look at approximately +/- 2mm */
  distance_to_surface = MAX(2.0, distance_to_surface);
  delta = distance_to_surface * mri_defect->xsize;
  delta = 3.0;  // for now, distance_to_surface = 2 * volume_resolution
  scale = mri_defect->xsize;

  /* initialize the signed image */
  for (k = 0; k < mri_distance->depth; k++)
    for (j = 0; j < mri_distance->height; j++)
      for (i = 0; i < mri_distance->width; i++) {
        MRIFvox(mri_distance, i, j, k) = NPY;
      }

  /* find the distance to each surface voxels */
  for (p = 0; p < mris_defect->nfaces + mris_source->nfaces; p++) {
    int fno_init;
    
    if (p >= mris_defect->nfaces) {
      mris = mris_source;
      fno_init = p - mris_defect->nfaces;
      wsurf = 1;  // using source surface
      if (mris->faces[fno_init].marked) {
        continue;
      }
    }
    else {
      fno_init = p;
      mris = mris_defect;
      wsurf = 0;  // using defect surface
    }
    
    int const fno = fno_init;
    FACE * const face = &mris->faces[fno];

    // calculate three vertices
    x0 = mris->vertices[face->v[0]].origx;
    y0 = mris->vertices[face->v[0]].origy;
    z0 = mris->vertices[face->v[0]].origz;
    x1 = mris->vertices[face->v[1]].origx;
    y1 = mris->vertices[face->v[1]].origy;
    z1 = mris->vertices[face->v[1]].origz;
    x2 = mris->vertices[face->v[2]].origx;
    y2 = mris->vertices[face->v[2]].origy;
    z2 = mris->vertices[face->v[2]].origz;

    /* find the bounding box */
    imin = iVOL(mri_defect, MIN3(x0, x1, x2)) - delta;
    imax = iVOL(mri_defect, MAX3(x0, x1, x2)) + delta;

    jmin = jVOL(mri_defect, MIN3(y0, y1, y2)) - delta;
    jmax = jVOL(mri_defect, MAX3(y0, y1, y2)) + delta;

    kmin = kVOL(mri_defect, MIN3(z0, z1, z2)) - delta;
    kmax = kVOL(mri_defect, MAX3(z0, z1, z2)) + delta;

    /* we don't count faces that are outside
       the volume - will not change the sign */
    if (imin > mri_defect->width - 1 || jmin > mri_defect->height - 1 || kmin > mri_defect->depth - 1 || imax < 0 ||
        jmax < 0 || kmax < 0) {
      continue;
    }

    imin = MAX(imin, 0);
    imax = MIN(imax, mri_defect->width - 1);

    jmin = MAX(jmin, 0);
    jmax = MIN(jmax, mri_defect->height - 1);

    kmin = MAX(kmin, 0);
    kmax = MIN(kmax, mri_defect->depth - 1);

    //////////////////////////////////////////////////////////////
    /* generating the pseudo-normals for edges and vertices */
    // normal for the current face
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    n_f[0] = fNorm->nx;
    n_f[1] = fNorm->ny;
    n_f[2] = fNorm->nz;

    /* edge0: x0 <--> x1 */
    {
    
    e0[0] = x1 - x0;
    e0[1] = y1 - y0;
    e0[2] = z1 - z0;
    F_CROSS(n_f, e0, n0);
    fn1 = findOtherEdgeFace(mris, fno, face->v[0], face->v[1]);

    MRIS* oface_mris_init = mris;
    int   oface_fno_init  = fn1;
    if (wsurf) {
      // source surface
      if (mris->faces[fn1].marked) {
        // border face
        // find face in mris_defect
        fn1 = findNonMarkedFace(mris_defect, vfrom[face->v[0]], vfrom[face->v[1]]);
        if (fn1 == -1) {
          // sanity check
          fprintf(stderr, "fn1=-1 in defect- should not happen\n");
          fn1 = fno;
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_defect;
      }
    }
    else {
      // defect surface
      if (fn1 == fno) {
        // border face
        // find face in mris_source
        fn1 = findNonMarkedFace(mris_source, vto[face->v[0]], vto[face->v[1]]);
        if (fn1 == -1) {
          fprintf(stderr, "fn1=-1 in source- should not happen\n");
          fn1 = fno;  // sanity check
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_source;
      }
    }
    MRIS*  const oface_mris = oface_mris_init;
    int    const oface_fno  = oface_fno_init;
    
    FaceNormCacheEntry const * ofNorm = getFaceNorm(oface_mris, oface_fno);
    
    n_e0[0] = fNorm->nx + ofNorm->nx;
    n_e0[1] = fNorm->ny + ofNorm->ny;
    n_e0[2] = fNorm->nz + ofNorm->nz;

    }
    
    /* edge1: x1 <--> x2 */    
    {
    
    e1[0] = x2 - x1;
    e1[1] = y2 - y1;
    e1[2] = z2 - z1;
    F_CROSS(n_f, e1, n1);
    fn1 = findOtherEdgeFace(mris, fno, face->v[1], face->v[2]);
    
    MRIS* oface_mris_init = mris;
    int   oface_fno_init  = fn1;
    if (wsurf) {
      // source surface
      if (mris->faces[fn1].marked) {
        // border face
        // find face in mris_defect
        fn1 = findNonMarkedFace(mris_defect, vfrom[face->v[1]], vfrom[face->v[2]]);
        if (fn1 == -1) {
          // sanity check
          fprintf(stderr, "fn1=-1 in defect- should not happen\n");
          fn1 = fno;
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_defect;
      }
    }
    else {
      // defect surface
      if (fn1 == fno) {
        // border face
        // find face in mris_source
        fn1 = findNonMarkedFace(mris_source, vto[face->v[1]], vto[face->v[2]]);
        if (fn1 == -1) {
          fprintf(stderr, "fn1=-1 in source- should not happen\n");
          fn1 = fno;  // sanity check
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_source;
      }
    }
    MRIS*  const oface_mris = oface_mris_init;
    int    const oface_fno  = oface_fno_init;

    FaceNormCacheEntry const * ofNorm = getFaceNorm(oface_mris, oface_fno);

    n_e1[0] = fNorm->nx + ofNorm->nx;
    n_e1[1] = fNorm->ny + ofNorm->ny;
    n_e1[2] = fNorm->nz + ofNorm->nz;

    }
    
    /* edge2: x2 <--> x0 */
    {
    
    e2[0] = x0 - x2;
    e2[1] = y0 - y2;
    e2[2] = z0 - z2;
    F_CROSS(n_f, e2, n2);
    fn1 = findOtherEdgeFace(mris, fno, face->v[2], face->v[0]);

    MRIS* oface_mris_init = mris;
    int   oface_fno_init  = fn1;

    if (wsurf) {
      // source surface
      if (mris->faces[fn1].marked) {
        // border face
        // find face in mris_defect
        fn1 = findNonMarkedFace(mris_defect, vfrom[face->v[2]], vfrom[face->v[0]]);
        if (fn1 == -1) {
          // sanity check
          fprintf(stderr, "fn1=-1 in defect- should not happen\n");
          fn1 = fno;
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_defect;
      }
    }
    else {
      // defect surface
      if (fn1 == fno) {
        // border face
        // find face in mris_source
        fn1 = findNonMarkedFace(mris_source, vto[face->v[2]], vto[face->v[0]]);
        if (fn1 == -1) {
          fprintf(stderr, "fn1=-1 in source- should not happen\n");
          fn1 = fno;  // sanity check
        }
        oface_fno_init  = fn1;              // note change to fn1 just above
        oface_mris_init = mris_source;
      }
    }

    MRIS*  const oface_mris = oface_mris_init;
    int    const oface_fno  = oface_fno_init;

    FaceNormCacheEntry const * ofNorm = getFaceNorm(oface_mris, oface_fno);
    
    n_e2[0] = fNorm->nx + ofNorm->nx;
    n_e2[1] = fNorm->ny + ofNorm->ny;
    n_e2[2] = fNorm->nz + ofNorm->nz;

    }

    /* vertex pseudo-normals */
    if (wsurf) {
      int vn;
      vn = face->v[0];
      if (mris->vertices[vn].marked) {
        // border
        // test
        if (vfrom[vn] < 0 || vfrom[vn] >= n_vertices) {
          fprintf(stderr, "problem with vfrom vn0\n");
        }
        // fprintf(stderr,"to0(%d-",vn);
        // fprintf(stderr,"%d),",vfrom[vn]);
        vertexPseudoNormal(mris_source, vn, mris_defect, vfrom[vn], n_v0);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v0, 0);
      }
      vn = face->v[1];
      if (mris->vertices[vn].marked) {
        // border
        // test
        if (vfrom[vn] < 0 || vfrom[vn] >= n_vertices) {
          fprintf(stderr, "problem with vfrom vn1\n");
        }
        // fprintf(stderr,"to1(%d-",vn);
        // fprintf(stderr,"%d),",vfrom[vn]);
        vertexPseudoNormal(mris_source, vn, mris_defect, vfrom[vn], n_v1);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v1, 0);
      }
      vn = face->v[2];
      if (mris->vertices[vn].marked) {
        // border
        // test
        if (vfrom[vn] < 0 || vfrom[vn] >= n_vertices) {
          fprintf(stderr, "problem with vfrom vn2\n");
        }
        // fprintf(stderr,"to2(%d-",vn);
        // fprintf(stderr,"%d),",vfrom[vn]);
        vertexPseudoNormal(mris_source, vn, mris_defect, vfrom[vn], n_v2);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v2, 0);
      }
    }
    else {
      int vn = face->v[0];
      // fprintf(stderr,"we have %d and %d \n", mris->vertices_topology[vn].vnum,mris->vertices_topology[vn].num);
      if (mris->vertices_topology[vn].vnum != mris->vertices_topology[vn].num) {
        // border
        // test
        if (vn >= n_vertices) {
          fprintf(stderr, "problem with vto vn0\n");
        }
        // fprintf(stderr,"fo0(%d-",vn);
        // fprintf(stderr,"%d),",vto[vn]);
        vertexPseudoNormal(mris_source, vto[vn], mris_defect, vn, n_v0);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v0, 0);
      }
      vn = face->v[1];
      if (mris->vertices_topology[vn].vnum != mris->vertices_topology[vn].num) {
        // border
        // test
        if (vn >= n_vertices) {
          fprintf(stderr, "problem with vto vn1\n");
        }
        // fprintf(stderr,"fo1(%d-",vn);
        // fprintf(stderr,"%d),",vto[vn]);
        vertexPseudoNormal(mris_source, vto[vn], mris_defect, vn, n_v1);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v1, 0);
      }
      vn = face->v[2];
      if (mris->vertices_topology[vn].vnum != mris->vertices_topology[vn].num) {
        // border
        // test
        if (vn >= n_vertices) {
          fprintf(stderr, "problem with vto vn2\n");
        }
        // fprintf(stderr,"fo2(%d-",vn);
        // fprintf(stderr,"%d),",vto[vn]);
        vertexPseudoNormal(mris_source, vto[vn], mris_defect, vn, n_v2);
      }
      else {
        computeVertexPseudoNormal(mris, vn, n_v2, 0);
      }
    }
    if (std::isnan(n_v0[0]) || std::isnan(n_v1[0]) || std::isnan(n_v2[0])) {
      fprintf(stderr,
              ".%d & %d[%d(%d) %d(%d) %d(%d)][%f %f %f]\n",
              wsurf,
              fno,
              face->v[0],
              mris->vertices[face->v[0]].marked,
              face->v[1],
              mris->vertices[face->v[1]].marked,
              face->v[2],
              mris->vertices[face->v[2]].marked,
              n_v0[0],
              n_v1[0],
              n_v2[0]);
      exit(-1);
    }
    //////////////////////////////////////////////////////////////

    /* finding distance to surface */
    for (k = kmin; k <= kmax; k++)
      for (j = jmin; j <= jmax; j++)
        for (i = imin; i <= imax; i++) {
          x = xSURF(mri_defect, i);
          y = ySURF(mri_defect, j);
          z = zSURF(mri_defect, k);

          vec0[0] = x - x0;
          vec0[1] = y - y0;
          vec0[2] = z - z0;
          vec1[0] = x - x1;
          vec1[1] = y - y1;
          vec1[2] = z - z1;
          vec2[0] = x - x2;
          vec2[1] = y - y2;
          vec2[2] = z - z2;
          vec[0] = (vec0[0] + vec1[0] + vec2[0]) / 3.0;
          vec[1] = (vec0[1] + vec1[1] + vec2[1]) / 3.0;
          vec[2] = (vec0[2] + vec1[2] + vec2[2]) / 3.0;

          /* compute distance to face */
          /* where is the point */
          val0 = F_DOT(vec0, n0);
          val1 = F_DOT(vec1, n1);
          val2 = F_DOT(vec2, n2);

          if ((val0 >= 0) && (val1 >= 0) && (val2 >= 0)) {
            /* the projection of the vertex is inside */
            val = F_DOT(n_f, vec);
            valu = 1;
            sign = val;
            distance = val; /* n_f is already normalized */
          }
          else {
            distance = NPY;
            sign = 0;
            valu = 0;

            if (val0 <= 0) {
              /* compute distance to edge0 */
              val = F_DOT(vec0, e0);
              if (val < 0) {
                /* closer to x0 */
                sign = F_DOT(n_v0, vec0);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
              }
              else if (val < SQR3(e0)) {
                /* closer to edge0 */
                sign = F_DOT(n_e0, vec0);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec0) - SQR(val) / SQR3(e0))));
              }
              else {
                /* closer to x1 */
                sign = F_DOT(n_v1, vec1);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
              }
            };
            if (val1 <= 0) {
              val = F_DOT(vec1, e1);
              if (val < 0) {
                /* closer to x1 */
                sign = F_DOT(n_v1, vec1);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec1));
              }
              else if (val < SQR3(e1)) {
                /* closer to edge1 */
                sign = F_DOT(n_e1, vec1);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec1) - SQR(val) / SQR3(e1))));
              }
              else {
                /* closer to x2 */
                sign = F_DOT(n_v2, vec2);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
              }
            };
            if (val2 <= 0) {
              val = F_DOT(vec2, e2);
              if (val < 0) {
                /* closer to x2 */
                sign = F_DOT(n_v2, vec2);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec2));
              }
              else if (val < SQR3(e2)) {
                /* closer to edge2 */
                sign = F_DOT(n_e2, vec2);
                valu = 3;
                distance = SIGN(sign) * MIN(fabs(distance), sqrt(MAX(0, SQR3(vec2) - SQR(val) / SQR3(e2))));
              }
              else {
                /* closer to x0 */
                sign = F_DOT(n_v0, vec0);
                valu = 2;
                distance = SIGN(sign) * MIN(fabs(distance), NORM3(vec0));
              }
            };
          }

          /* update distance map */
          if (fabs(distance) < fabs(MRIFvox(mri_distance, i, j, k))) {
            MRIFvox(mri_distance, i, j, k) = distance;
          }
        }
  }
  for (n = 0; n < n_faces; n++) {
    mris_source->faces[fto[n]].marked = 0;
  }
  for (n = 0; n < n_vertices; n++) {
    mris_source->vertices[vto[n]].marked = 0;
  }

  // update the volume
  // only keep the points close to the surface (we only care about the sign)
  for (k = 0; k < mri_distance->depth; k++)
    for (j = 0; j < mri_distance->height; j++)
      for (i = 0; i < mri_distance->width; i++) {
        if (i == 0 || j == 0 || k == 0 || i == mri_distance->width - 1 || j == mri_distance->height - 1 ||
            k == mri_distance->depth - 1) {
          MRIFvox(mri_distance, i, j, k) = NPY;
          continue;
        }
        distance = MRIFvox(mri_distance, i, j, k);
        if (distance == NPY) {
          continue;
        }
        if (fabs(distance) * scale > 1.5) {
          MRIFvox(mri_distance, i, j, k) = NPY;
        }
        else if (distance > 0.0) {
          MRIFvox(mri_distance, i, j, k) = MIN(distance * scale, 1.0f);
        }
        else {
          MRIFvox(mri_distance, i, j, k) = MAX(distance * scale, -1.0f);
        }
      }

  found = 1;
  while (found) {
    found = 0;
    for (k = 1; k < mri_distance->depth - 1; k++)
      for (j = 1; j < mri_distance->height - 1; j++)
        for (i = 1; i < mri_distance->width - 1; i++) {
          if (MRIFvox(mri_distance, i, j, k) < NPY) {
            continue;
          }
          
          static int _DX[6] = {-1, 1, 0, 0, 0, 0};
          static int _DY[6] = {0, 0, -1, 1, 0, 0};
          static int _DZ[6] = {0, 0, 0, 0, -1, 1};

          for (n = 0; n < 6; n++) {
            float sign;
            sign = MRIFvox(mri_distance, i + _DX[n], j + _DY[n], k + _DZ[n]);
            if (fabs(sign) < NPY) {
              if (sign > 0) {
                MRIFvox(mri_distance, i, j, k) = 1.0f;
              }
              else {
                MRIFvox(mri_distance, i, j, k) = -1.0f;
              }
              found = 1;
              break;
            }
          }
        }
  }
  // making sure
  for (k = 0; k < mri_distance->depth - 0; k++)
    for (j = 0; j < mri_distance->height - 0; j++)
      for (i = 0; i < mri_distance->width - 0; i++) {
        sign = MRIFvox(mri_distance, i, j, k);
        if (sign > 0.0) {
          MRIFvox(mri_distance, i, j, k) = MIN(1.0f, sign);
        }
        else {
          MRIFvox(mri_distance, i, j, k) = MAX(-1.0f, sign);
        }
      }
}


#define SAVING_SURFACES 0

#define NB_OF_CLUSTERS 2

static int findNumberOfClusters(MRIS *mris) { return NB_OF_CLUSTERS; }
/* using orig vertices to cluster
   fixedval indexes the clusters */
static int clusterDefectVertices(MRIS *mris)
{
  int n, niters, nvertices, k;
  int p, q, nclusters, found;
  int clusters[10], ngroups;
  float dist, ndist, x, y, z, nbv;
  VERTEX *v, *vq;
  int rgb[10], v_rgb;

  ngroups = findNumberOfClusters(mris);
  nclusters = ngroups;

  MRISRGBToAnnot(0, 225, 225, rgb[0]);
  MRISRGBToAnnot(205, 62, 78, rgb[1]);
  MRISRGBToAnnot(120, 62, 78, rgb[2]);
  MRISRGBToAnnot(196, 58, 250, rgb[3]);
  MRISRGBToAnnot(0, 148, 0, rgb[4]);
  MRISRGBToAnnot(220, 248, 164, rgb[5]);
  MRISRGBToAnnot(230, 148, 34, rgb[6]);
  MRISRGBToAnnot(0, 118, 14, rgb[7]);
  MRISRGBToAnnot(12, 48, 255, rgb[8]);
  MRISRGBToAnnot(122, 186, 220, rgb[9]);

  niters = 10;

  nvertices = mris->nvertices;

  /* set marks to zero */
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    v->fixedval = 0;
  }
  /* first cluster = border points */
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->flags == VERTEX_INTERIOR) {
      continue;
    }
    v->fixedval = 1;
  }
  /* initialize other clusters */
  for (p = 1; p < nclusters; p++) {
    found = 0;
    while (!found) {
      //                        fprintf(stderr,".");
      /* draw a random number */
      k = nint(randomNumber(0.0, (double)nvertices - 0.7));
      k = MIN(nvertices, MAX(0, k));
      v = &mris->vertices[k];
      if (v->fixedval) {
        continue;
      }
      // fprintf(stderr,"*%d-%d*",k,p+1);
      found = 1;
      v->fixedval = p + 1;
      clusters[p] = k;
    }
  }
  /* classify vertices */
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    if (v->fixedval) {
      continue;
    }
    /* find the correct cluster */
    /* compute distance with first cluster */
    dist = 1000;
    v->fixedval = 1;
#if 0
    for (q = 0 ; q < mris->nvertices ; q++ )
    {
      vq=&mris->vertices[q];
      if (vq->flags==VERTEX_INTERIOR)
      {
        continue;
      }
      dist=MIN(dist,SQR(vq->origx-v->origx)+
               SQR(vq->origy-v->origy)+SQR(vq->origz-v->origz));
    }
#endif
    /* compute distance for the remaining clusters */
    for (q = 1; q < nclusters; q++) {
      vq = &mris->vertices[clusters[q]];
      ndist = SQR(vq->origx - v->origx) + SQR(vq->origy - v->origy) + SQR(vq->origz - v->origz);
      if (ndist < dist) {
        dist = ndist;
        v->fixedval = q + 1;
      }
    }
  }
  for (q = 0; q < mris->nvertices; q++) {
    vq = &mris->vertices[q];
    vq->curv = vq->fixedval;
  }

  /*    MRISwriteCurvature(mris,"./lh.test.curv1"); */

  while (niters--) {
    for (p = 1; p < nclusters; p++) {
      /* regenerate clusters' center */
      x = y = z = 0;
      nbv = 0;
      for (n = 0; n < mris->nvertices; n++) {
        v = &mris->vertices[n];
        if (v->fixedval != p + 1) {
          continue;
        }
        x += v->origx;
        y += v->origy;
        z += v->origz;
        nbv += 1.0f;
      }
      if (!nbv) {
        fprintf(stderr, "ERROR in clusterDefectVertices");
      }
      else {
        x /= nbv;
        y /= nbv;
        z /= nbv;
      }

      /* update clustering */
      dist = 1000;
      for (n = 0; n < mris->nvertices; n++) {
        v = &mris->vertices[n];
        if (v->fixedval != p + 1) {
          continue;
        }
        ndist = SQR(v->origx - x) + SQR(v->origy - y) + SQR(v->origz - z);
        if (ndist < dist) {
          clusters[p] = n;
          dist = ndist;
        }
      }
    }
    /* re-classify vertices */
    for (n = 0; n < mris->nvertices; n++) {
      v = &mris->vertices[n];
      if (v->flags != VERTEX_INTERIOR) {
        continue;
      }
      /* find the correct cluster */
      /* compute distance with first cluster */
      dist = 1000;
#if 0
      v->fixedval=1;
      for (q=0; q<mris->nvertices; q++)
      {
        vq=&mris->vertices[q];
        if (vq->flags==VERTEX_INTERIOR)
        {
          continue;
        }
        dist=MIN(dist,SQR(vq->origx-v->origx)+
                 SQR(vq->origy-v->origy)+SQR(vq->origz-v->origz));
      }
#endif
      /* compute distance for the remaining clusters */
      for (q = 1; q < nclusters; q++) {
        vq = &mris->vertices[clusters[q]];
        ndist = SQR(vq->origx - v->origx) + SQR(vq->origy - v->origy) + SQR(vq->origz - v->origz);
        if (ndist < dist) {
          dist = ndist;
          v->fixedval = q + 1;
        }
      }
    }
  }
  for (q = 0; q < mris->nvertices; q++) {
    vq = &mris->vertices[q];
    vq->curv = vq->fixedval;
    v_rgb = rgb[(int)(vq->fixedval) % 10];
    vq->annotation = v_rgb;
  }

#if SAVING_SURFACES
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  MRISwrite(mris, "./test.now");
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
  MRISwriteCurvature(mris, "./lh.test.curv2");
  MRISwriteAnnotation(mris, "./lh.ant");
#endif

  return ngroups;
}

static void computeInteriorGradients(MRIS *mris, int option)
{
  int n, p, count;
  float x, y, z;

  for (n = 0; n < mris->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[n];
    VERTEX                * const v  = &mris->vertices         [n];
    
    if (v->flags != VERTEX_INTERIOR) {
      continue;
    }

    x = y = z = 0.0f;

    for (count = p = 0; p < vt->vnum; p++) {
      VERTEX const * const vp = &mris->vertices[vt->v[p]];
      if (v->fixedval != option && vp->fixedval == option) {
        continue;
      }
      x += vp->x;
      y += vp->y;
      // z += vp->z;
      count++;
    }

    if (count) {
      /* compute average */
      x /= (float)count;
      y /= (float)count;
      // z /= (float)count;
      /* generate gradient */
      v->dx += (x - v->x);
      v->dy += (y - v->y);
      // v->dz += (z-v->z);
    }
  }
}

/* use fixedval to exclude some random vertices */
static void computeGradients(MRIS *mris, int option)
{
  int n;
  VERTEX *v;

  /* set gradients to zero */
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    v->dx = v->dy = v->dz = 0.0f;
  }

  /* repulsive term for chull vertices */
  // computeChullGradients(mris);

  /* repulsive term for border vertices */
  // computeBorderGradients(mris);

  /* compute spring term for interior vertices */
  computeInteriorGradients(mris, option);
}

static void circleProjection(float x, float y, float *xd, float *yd, float r)
{
  float dist;

  dist = sqrt(SQR(x) + SQR(y));
  if (dist > r) {
    fprintf(stderr, "!");
    *xd = x * 0.99 * r / dist;
    *yd = y * 0.99 * r / dist;
  }
}

/* generate a place xy mapping */
static void generateDefectMapping(MRIS *mris, int option)
{
  int n, niter = 100;
  double dt = 0.1;
  VERTEX *v;

#if 0
  {
    int n,vlist[8];
    for (n=0; n<mris->nvertices; n++)
    {
      mris->vertices[n].fixedval=0;
    }
    //option=2;
    vlist[0]=508;
    vlist[1]=509;
    vlist[2]=538;
    vlist[3]=539;
    vlist[4]=540;
    vlist[5]=541;
    vlist[6]=569;
    vlist[7]=510;
    for (n=0; n<8; n++)
    {
      mris->vertices[vlist[n]].fixedval=option;
    }
  }
  /*    fprintf(stderr," label is %d or %d or %d :
        option=%d\n",mris->vertices[361].fixedval,
        mris->vertices[392].fixedval,mris->vertices[391].fixedval,option); */
#endif

  while (niter--) {
    if (niter % 20 == 0) {
      fprintf(stderr, ".");
    }

    /* compute gradients */
    computeGradients(mris, option);
    /* apply gradients */
    for (n = 0; n < mris->nvertices; n++) {
      v = &mris->vertices[n];
      if (v->flags == VERTEX_BORDER) {
        continue;
      }
      v->x += dt * v->dx;
      v->y += dt * v->dy;
      circleProjection(v->x, v->y, &v->x, &v->y, mris->radius);
      v->cx = v->x;
      v->cy = v->y;
    }
  }
}

static OPTIMAL_DEFECT_MAPPING *mrisFindOptimalDefectMapping(MRIS *mris_src, DEFECT *defect)
{
  int nvertices, nfaces, nchull;
  int *vertex_trans, *face_trans, *vertex_list;
  int vno, n, m, i, first_inside_vertex, first_border_vertex;
  MAPPING *mapping;
  FS_VERTEX_INFO *vinfo;
  FACE *face, *face_dst;

  MRIS *mris_dst;

  OPTIMAL_DEFECT_MAPPING *o_d_m;
  int option = 0, nclusters;
  EDGE *edge;

  /* first neighbors only */
  MRISsetNeighborhoodSizeAndDist(mris_src, 1);

  nvertices = defect->nvertices + defect->nchull;
  vertex_list = (int *)malloc(nvertices * sizeof(int));

  vertex_trans = (int *)malloc(mris_src->nvertices * sizeof(int));
  face_trans   = (int *)malloc(mris_src->nfaces * sizeof(int));
  memset(vertex_trans, -1, mris_src->nvertices * sizeof(int));
  memset(face_trans,   -1, mris_src->nfaces * sizeof(int));

  nvertices = 0;

  nchull = defect->nchull - defect->nborder; /* number of chull vertices :
                                                should be zero */
  for (n = 0; n < nchull; n++) {
    /* first chull */
    vertex_list[nvertices] = defect->chull[n + defect->nborder];
    vertex_trans[defect->chull[defect->nborder + n]] = nvertices;
    nvertices++;
  }
  first_border_vertex = nvertices;
  for (n = 0; n < defect->nborder; n++) {
    /* then border */
    vertex_list[nvertices] = defect->chull[n];
    vertex_trans[defect->chull[n]] = nvertices;
    nvertices++;
  }
  first_inside_vertex = nvertices;
  for (n = 0; n < defect->nvertices; n++) {
    /* finally inside defect vertices */
    vertex_list[nvertices] = defect->vertices[n];
    vertex_trans[defect->vertices[n]] = nvertices;
    nvertices++;
  }

  /* mark vertices */
  mrisMarkDefectConvexHull(mris_src, defect, 1);
  mrisMarkDefect(mris_src, defect, 1);

  for (nfaces = n = 0; n < mris_src->nfaces; n++) {
    face = &mris_src->faces[n];
    if (mris_src->vertices[face->v[0]].marked && mris_src->vertices[face->v[1]].marked &&
        mris_src->vertices[face->v[2]].marked) {
      face_trans[n] = nfaces;
      nfaces++;
    }
  }

  /* allocate temporary surface */
  mris_dst = MRISalloc(nvertices, nfaces);
  mris_dst->type   = MRIS_TRIANGULAR_SURFACE;
  mris_dst->status = MRIS_SPHERICAL_PATCH;
  mris_dst->origxyz_status = mris_src->origxyz_status;
  
  mris_dst->radius = DEFAULT_RADIUS;

  /* copy faces */
  for (n = 0; n < mris_src->nfaces; n++) {
    if (face_trans[n] < 0) {
      continue;
    }
    face = &mris_src->faces[n];
    face_dst = &mris_dst->faces[face_trans[n]];
    for (m = 0; m < VERTICES_PER_FACE; m++) {
      face_dst->v[m] = vertex_trans[face->v[m]];
    }
  }

  /* copy vertices with their neighbors */
  for (n = 0; n < mris_dst->nvertices; n++) {
    int const vno_dst = n;
    int const vno_src = vertex_list[n];
    
    VERTEX_TOPOLOGY       * const v_dstt = &mris_src->vertices_topology[vno_dst];
    VERTEX                * const v_dst  = &mris_dst->vertices         [vno_dst];
    VERTEX_TOPOLOGY const * const v_srct = &mris_src->vertices_topology[vno_src];
    VERTEX          const * const v_src  = &mris_src->vertices         [vno_src];
    
    /* useless since we reinitialize the locations */
    /* making sure the vertices are in canonical space */
    v_dst->x = v_src->cx;
    v_dst->y = v_src->cy;
    v_dst->z = v_src->cz;
    v_dst->tx = v_src->tx;
    v_dst->ty = v_src->ty;
    v_dst->tz = v_src->tz;
    v_dst->nx = v_src->nx;
    v_dst->ny = v_src->ny;
    v_dst->nz = v_src->nz;
    v_dst->cx = v_src->cx;
    v_dst->cy = v_src->cy;
    v_dst->cz = v_src->cz;
    
    MRISsetOriginalXYZ(mris_dst, vno_dst, v_src->origx, v_src->origy, v_src->origz);
    
    v_dst->ripflag = v_src->ripflag; /* none of them should be ripped */

    if (n < nchull) /* vertex in the convex hull */
    {
      v_dst->flags = VERTEX_CHULL;
    }
    else if (n < defect->nchull) /* vertex in the border */
    {
      v_dst->flags = VERTEX_BORDER;
    }
    else /* vertex inside the defect */
    {
      v_dst->flags = VERTEX_INTERIOR;
    }

    /* count the number of kept neighboring vertices */
    if (n < nchull) {
      /* if n < nchull, we need to watch for the
         right neighboring vertices/faces */
      /* count # of valid neighbors */
      
      int vnum = 0;
      for (m = 0; m < v_srct->vnum; m++)
        if (mris_src->vertices[v_srct->v[m]].marked) {
          vnum++;
        }
        
      mrisVertexReplacingNeighbors(mris_dst, vno_dst, vnum);
      for (i = m = 0; m < v_srct->vnum; m++)
        if (mris_src->vertices[v_srct->v[m]].marked) {
          v_dstt->v[i] = vertex_trans[v_srct->v[m]];
          v_dst->dist[i] = v_src->dist[m];
          v_dst->dist_orig[i] = v_src->dist_orig[m];
          i++;
        }

      /* count # of good triangles attached to this vertex */
      for (v_dstt->num = m = 0; m < v_srct->num; m++) {
        face = &mris_src->faces[v_srct->f[m]];
        if (mris_src->vertices[face->v[0]].marked && mris_src->vertices[face->v[1]].marked &&
            mris_src->vertices[face->v[2]].marked) {
          v_dstt->num++;
        }
      }

      v_dstt->f = (int *)calloc(v_dstt->num, sizeof(int));
      v_dstt->n = (uchar *)calloc(v_dstt->num, sizeof(uchar));
      for (i = m = 0; m < v_srct->num; m++) {
        face = &mris_src->faces[v_srct->f[m]];
        if (mris_src->vertices[face->v[0]].marked && mris_src->vertices[face->v[1]].marked &&
            mris_src->vertices[face->v[2]].marked) {
          v_dstt->n[i] = v_srct->n[m];
          v_dstt->f[i] = face_trans[v_srct->f[m]];
          i++;
        }
      }
    }
    else {
      /* neighboring vertices */
      mrisVertexReplacingNeighbors(mris_dst, vno_dst, v_srct->vnum);
      for (m = 0; m < v_srct->vnum; m++) {
        v_dstt->v[m] = vertex_trans[v_srct->v[m]];
        v_dst->dist[m] = v_src->dist[m];
        v_dst->dist_orig[m] = v_src->dist_orig[m];
      }

      /* neighboring faces */
      v_dstt->num = v_srct->num;
      v_dstt->f = (int *)calloc(v_dstt->num, sizeof(int));
      v_dstt->n = (uchar *)calloc(v_dstt->num, sizeof(uchar));
      for (m = 0; m < v_srct->num; m++) {
        v_dstt->n[m] = v_srct->n[m];
        v_dstt->f[m] = face_trans[v_srct->f[m]];
      }
    }
  }

  mrisCheckVertexFaceTopology(mris_dst);
  
  /* unmark vertices */
  mrisMarkDefectConvexHull(mris_src, defect, 0);
  mrisMarkDefect(mris_src, defect, 0);

  o_d_m = (OPTIMAL_DEFECT_MAPPING *)calloc(1, sizeof(OPTIMAL_DEFECT_MAPPING));
  o_d_m->mris = mris_dst;
  o_d_m->vertex_trans = vertex_trans;
  o_d_m->face_trans = face_trans;
  o_d_m->orig_mapping.vertices = (FS_VERTEX_INFO *)calloc(mris_dst->nvertices, sizeof(FS_VERTEX_INFO));
  o_d_m->orig_mapping.nvertices = mris_dst->nvertices;
  o_d_m->orig_mapping.ninside = first_inside_vertex;
  for (n = 0; n < 10; n++) {
    o_d_m->mappings[n].vertices = (FS_VERTEX_INFO *)calloc(mris_dst->nvertices, sizeof(FS_VERTEX_INFO));
    o_d_m->mappings[n].nvertices = mris_dst->nvertices;
    o_d_m->mappings[n].ninside = first_inside_vertex;
  }

  /* saving original configurations */
  mapping = &o_d_m->orig_mapping;
  for (n = 0; n < defect->nvertices; n++) {
    vinfo = &mapping->vertices[n + mapping->ninside];
    vinfo->status = defect->status[n];
  }
  /* initializing border vertices */
  for (n = 0; n < defect->nedges; n++) {
    /* init border */
    edge = &defect->edges[n];
    vno = vertex_trans[edge->vno1];
    vinfo = &mapping->vertices[vno];
    VERTEX const * const v = &mris_dst->vertices[vno];
    /* circle */
    vinfo->c_x = 20.0 * cos(2 * PI * (float)n / defect->nedges);
    vinfo->c_y = 20.0 * sin(2 * PI * (float)n / defect->nedges);
    vinfo->c_z = 0;
    vinfo->oc_x = v->origx; /* orig coord */
    vinfo->oc_y = v->origy;
    vinfo->oc_z = v->origz;
  }
  /* init radius critical distance */
  mris_dst->radius = 20.0 * cos(PI / (float)defect->nedges);
  // fprintf(stderr,"critical radius is %f\n",mris_dst->radius);

  for (n = 0; n < defect->nvertices; n++) {
    /* init inside */
    vno = vertex_trans[defect->vertices[n]];
    vinfo = &mapping->vertices[vno];
    VERTEX const * const v = &mris_dst->vertices[vno];
    vinfo->c_x = 0;
    vinfo->c_y = 0;
    vinfo->c_z = 0;
    vinfo->oc_x = v->origx; /* orig coord */
    vinfo->oc_y = v->origy;
    vinfo->oc_z = v->origz;
  }

  /* clustering defect vertices into groups */
  nclusters = clusterDefectVertices(mris_dst);

  o_d_m->nmappings = nclusters + 1; /* 1 = smoothed ; then others */

  /////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////
  // Now, we are ready to generate different mappings...

  for (m = 0; m < o_d_m->nmappings; m++) {
    fprintf(stderr, "generating mapping #%d\n", m);

    /* transfer initial coordinates into surface */
    mapping = &o_d_m->orig_mapping;
    for (n = 0; n < mris_dst->nvertices; n++) {
      vinfo = &mapping->vertices[n];
      VERTEX * const v = &mris_dst->vertices[n];
      v->cx = vinfo->c_x; /* plane xy coord */
      v->cy = vinfo->c_y;
      v->cz = vinfo->c_z;
      v->x = v->cx;
      v->y = v->cy;
      v->z = v->cz;
    }
    fprintf(stderr, ".");
#if SAVING_SURFACES
    if (m == 0) {
      MRISwrite(mris_dst, "lh.test_init");
    }
#endif

    option = 0;
    if (m) {
      option = m + 1;
    }

    /* generate the mapping */
    generateDefectMapping(mris_dst, option);

#if SAVING_SURFACES
    {
      char fname[100];
      sprintf(fname, "lh.test_%d", m);
      fprintf(stderr, "writting file into %s\n", fname);
      MRISwrite(mris_dst, fname);
    }
#endif
    /* save the mapping onto the sphere (hack!) */
    mapping = &o_d_m->mappings[m];
    for (n = 0; n < mris_dst->nvertices; n++) {
      vinfo = &mapping->vertices[n];
      VERTEX * const v = &mris_dst->vertices[n];
      vinfo->c_x = v->cx; /* spherical coord */
      vinfo->c_y = v->cy;
      v->cz = sqrt(10000.0 - SQR(v->cx) - SQR(v->cy));
      vinfo->c_z = v->cz;
      v->x = v->cx;
      v->y = v->cy;
      v->z = v->cz;
    }
#if SAVING_SURFACES
    {
      char fname[100];
      sprintf(fname, "lh.tests_%d", m);
      fprintf(stderr, "writting file into %s\n", fname);
      MRISwrite(mris_dst, fname);
    }
#endif
  }

#if 0
  free(e_l_i.border_edges);
  free(e_l_i.inside_edges);
#endif
  free(vertex_list);

  return o_d_m;
}

static int mrisTessellateDefect_wkr(MRI_SURFACE *mris,
                                MRI_SURFACE *mris_corrected,
                                DEFECT *defect,
                                int *vertex_trans,
                                MRI *mri,
                                HISTOGRAM *h_k1,
                                HISTOGRAM *h_k2,
                                MRI *mri_k1_k2,
                                HISTOGRAM *h_white,
                                HISTOGRAM *h_gray,
                                HISTOGRAM *h_border,
                                HISTOGRAM *h_grad,
                                MRI *mri_gray_white,
                                HISTOGRAM *h_dot,
                                TOPOLOGY_PARMS *parms);
				
static int mrisTessellateDefect(MRI_SURFACE *mris,
                                MRI_SURFACE *mris_corrected,
                                DEFECT *defect,
                                int *vertex_trans,
                                MRI *mri,
                                HISTOGRAM *h_k1,
                                HISTOGRAM *h_k2,
                                MRI *mri_k1_k2,
                                HISTOGRAM *h_white,
                                HISTOGRAM *h_gray,
                                HISTOGRAM *h_border,
                                HISTOGRAM *h_grad,
                                MRI *mri_gray_white,
                                HISTOGRAM *h_dot,
                                TOPOLOGY_PARMS *parms) {
  fprintf(stderr,
          "CORRECTING DEFECT %d (vertices=%d, convex hull=%d, v0=%d)\n",
          defect->defect_number,
          defect->nvertices,
          defect->nchull,
          defect->vertices[0]);

  // TIMER_INTERVAL_BEGIN(old);
  
  int result = mrisTessellateDefect_wkr(
    mris,mris_corrected,defect,vertex_trans,mri,h_k1,h_k2,mri_k1_k2,h_white,h_gray,h_border,h_grad,mri_gray_white,h_dot,parms);

  // TIMER_INTERVAL_END(old);
  
  return result;
}

static int compare_edge_length(const void *pe0, const void *pe1)
{
  EDGE *e0, *e1;

  e0 = (EDGE *)pe0;
  e1 = (EDGE *)pe1;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (e0->len > e1->len) {
    return(1);
  }
  else if (e0->len < e1->len) {
    return(-1);
  }
  else if (e0->vno1 > e1->vno1) {
    // if two lengths are the same, force a consistent order by
    // comparing the vno1 values for each edge. Otherwise, returning 0
    // will produce sorting differences on osx
    return(1);
  }
  return(-1);
}
				
static int mrisTessellateDefect_wkr(MRI_SURFACE *mris,
                                MRI_SURFACE *mris_corrected,
                                DEFECT *defect,
                                int *vertex_trans,
                                MRI *mri,
                                HISTOGRAM *h_k1,
                                HISTOGRAM *h_k2,
                                MRI *mri_k1_k2,
                                HISTOGRAM *h_white,
                                HISTOGRAM *h_gray,
                                HISTOGRAM *h_border,
                                HISTOGRAM *h_grad,
                                MRI *mri_gray_white,
                                HISTOGRAM *h_dot,
                                TOPOLOGY_PARMS *parms)
{
  int i, j, vlist[MAX_DEFECT_VERTICES], n, nvertices, nedges, ndiscarded;
  VERTEX *v, *v2;
  EDGE *et;
  /*  double  cx, cy, cz, max_len ;*/
  static int dno = 0;
  double x, y, z, xv, yv, zv, val0, val, total, dx, dy, dz, d, wval, gval, Ix, Iy, Iz;
  float norm1[3], norm2[3], nx, ny, nz;
  int nes; /* number of edges present in original tessellation */
  ES *es;  /* list of edges present in original tessellation */
  /*generate an initial ordering*/
  int *ordering = NULL;


  ROMP_SCOPE_begin
  
  /* first build table of all possible edges among vertices in the defect
     and on its border.
  */
  if (parms->search_mode != GREEDY_SEARCH)
    computeDefectStatistics(mri, mris, defect, h_white, h_gray, mri_gray_white, h_k1, h_k2, mri_k1_k2, 0);

  /* first build table of all possible edges among vertices in the defect
     and on its border.*/
  for (nes = nvertices = i = 0; i < defect->nvertices; i++) {
    if (nvertices >= MAX_DEFECT_VERTICES)
      ErrorExit(ERROR_NOMEMORY, "mrisTessellateDefect: too many vertices in defect (%d)", MAX_DEFECT_VERTICES);
    if (defect->status[i] == KEEP_VERTEX) {
      vlist[nvertices++] = defect->vertices[i];
    }
  }

  for (i = 0; i < defect->nborder; i++) {
    vlist[nvertices++] = defect->border[i];
  }
  //  if (nvertices > 250)  //FLO
  if (DIAG_VERBOSE_ON)
    fprintf(WHICH_OUTPUT,
            "retessellating defect %d with %d vertices (convex hull=%d).\n",
            defect->defect_number,
            nvertices,
            defect->nchull);
  dno++;
  if (nvertices == 0) /* should never happen */
  {
    return (NO_ERROR);
  }

  nedges = (nvertices * (nvertices - 1)) / 2; /* won't be more than this */

  et = (EDGE *)calloc(nedges, sizeof(EDGE));
  if (!et)
    ErrorExit(ERROR_NOMEMORY,
              "Excessive topologic defect encountered: "
              "could not allocate %d edges for retessellation",
              nedges);

  ROMP_SCOPE_end
  
  n = 0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(serial)
#endif
  for (i = 0; i < nvertices; i++) {
    ROMP_PFLB_begin
    
    v = &mris->vertices[vlist[i]];
    if (vlist[i] == Gdiag_no) {
      DiagBreak();
    }
    if (vertex_trans[vlist[i]] == Gdiag_no) {
      DiagBreak();
    }
    for (j = i + 1; j < nvertices; j++, n++) {
      if (vlist[j] == Gdiag_no) {
        DiagBreak();
      }
      if (vlist[j] == Gdiag_no || vlist[i] == Gdiag_no) {
        DiagBreak();
      }
      if (vertex_trans[vlist[j]] == Gdiag_no || vertex_trans[vlist[i]] == Gdiag_no) {
        DiagBreak();
      }
      mrisComputeOrigNormal(mris, vlist[i], norm1);
      mrisComputeOrigNormal(mris, vlist[j], norm2);
      nx = (norm1[0] + norm2[0]) / 2;
      ny = (norm1[1] + norm2[1]) / 2;
      nz = (norm1[2] + norm2[2]) / 2;
      total = sqrt(nx * nx + ny * ny + nz * nz);
      if (FZERO(total)) {
        total = 1;
      }
      nx /= total;
      ny /= total;
      nz /= total;
      v2 = &mris->vertices[vlist[j]];
      x = (v->origx + v2->origx) / 2;
      y = (v->origy + v2->origy) / 2;
      z = (v->origz + v2->origz) / 2;
// MRIworldToVoxel(mri, x, y, z, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x, y, z, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x, y, z, &xv, &yv, &zv);
#endif
      MRIsampleVolumeGradient(mri, xv, yv, zv, &Ix, &Iy, &Iz);
      total = sqrt(Ix * Ix + Iy * Iy + Iz * Iz);
      if (FZERO(total)) {
        total = 1;
      }
      Ix /= total;
      Iy /= total;
      Iz /= total;

      /* assign value for edge table : values in mris_corrected */
      et[n].vno1 = vertex_trans[vlist[i]];
      et[n].vno2 = vertex_trans[vlist[j]];
      if ((et[n].vno1 == 141823 && et[n].vno2 == 141908) || (et[n].vno2 == 141823 && et[n].vno1 == 141908)) {
        DiagBreak();
      }
      if ((vlist[i] == Gdiag_no && vlist[j] == Gx) || (vlist[j] == Gdiag_no && vlist[i] == Gx)) {
        DiagBreak();
      }
      if ((vertex_trans[vlist[i]] == Gdiag_no && vertex_trans[vlist[j]] == Gx) ||
          (vertex_trans[vlist[j]] == Gdiag_no && vertex_trans[vlist[i]] == Gx)) {
        DiagBreak();
      }

      /* sample MR values along line and build estimate of log likelihood
         as distance.
      */
      val0 = (v->val + v2->val) / 2;         /* gray/white border value */
      wval = (v->val2 + v2->val2) / 2;       /* white matter mean */
      gval = (v->val2bak + v2->val2bak) / 2; /* gray matter mean */

      /* sample one end point */
      x = v->origx;
      y = v->origy;
      z = v->origz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);
      total = fabs(val - gval);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);
      total += fabs(val - wval);

      /* sample the other end point */
      x = v2->origx;
      y = v2->origy;
      z = v2->origz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);
      total += fabs(val - gval);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
      mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
      MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
      MRIsampleVolume(mri, xv, yv, zv, &val);
      total += fabs(val - wval);

      dx = v2->origx - v->origx;
      dy = v2->origy - v->origy;
      dz = v2->origz - v->origz;
      for (d = .1; d <= .9; d += .1) {
        /* sample the midpoint end point */
        x = v->origx + d * dx;
        y = v->origy + d * dy;
        z = v->origz + d * dz;
// MRIworldToVoxel(mri, x+.5*nx, y+.5*ny, z+.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
        mriSurfaceRASToVoxel(x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#else
        MRISsurfaceRASToVoxelCached(mris, mri, x + .5 * nx, y + .5 * ny, z + .5 * nz, &xv, &yv, &zv);
#endif
        MRIsampleVolume(mri, xv, yv, zv, &val);
        total += fabs(val - gval);
// MRIworldToVoxel(mri, x-.5*nx, y-.5*ny, z-.5*nz, &xv, &yv, &zv) ;
#if MATRIX_ALLOCATION
        mriSurfaceRASToVoxel(x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#else
        MRISsurfaceRASToVoxelCached(mris, mri, x - .5 * nx, y - .5 * ny, z - .5 * nz, &xv, &yv, &zv);
#endif
        MRIsampleVolume(mri, xv, yv, zv, &val);
        total += fabs(val - wval);
      }

      et[n].len = total / (4.0 + 18.0);
      if (et[n].vno1 == 120811 && et[n].vno2 == 120951) {
        VERTEX *v1, *v2;
        v1 = &mris_corrected->vertices[et[n].vno1];
        v2 = &mris_corrected->vertices[et[n].vno2];
        fprintf(stdout, "v %d (%d) --> %d (%d), len = %2.3f\n", et[n].vno1, vlist[i], et[n].vno2, vlist[j], et[n].len);
        fprintf(stdout,
                "INFLATED:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->tx,
                v1->ty,
                v1->tz,
                v2->tx,
                v2->ty,
                v2->tz,
                sqrt(SQR(v1->tx - v2->tx) + SQR(v1->ty - v2->ty) + SQR(v1->tz - v2->tz)));
        fprintf(stdout,
                "CANON:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->cx,
                v1->cy,
                v1->cz,
                v2->cx,
                v2->cy,
                v2->cz,
                sqrt(SQR(v1->cx - v2->cx) + SQR(v1->cy - v2->cy) + SQR(v1->cz - v2->cz)));
        fprintf(stdout,
                "CURRENT:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->x,
                v1->y,
                v1->z,
                v2->x,
                v2->y,
                v2->z,
                sqrt(SQR(v1->x - v2->x) + SQR(v1->y - v2->y) + SQR(v1->z - v2->z)));
        fprintf(stdout,
                "ORIG:  (%2.1f, %2.1f, %2.1f) --> "
                "(%2.1f, %2.1f, %2.1f), len = %2.2f\n",
                v1->origx,
                v1->origy,
                v1->origz,
                v2->origx,
                v2->origy,
                v2->origz,
                sqrt(SQR(v1->origx - v2->origx) + SQR(v1->origy - v2->origy) + SQR(v1->origz - v2->origz)));
        DiagBreak();
      }
      if (edgeExists(mris_corrected, et[n].vno1, et[n].vno2)) {
        et[n].used = USED_IN_TESSELLATION;
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "excluding existing edge %d <--> %d\n", vlist[i], vlist[j]);
      }

      if (!edgeExists(mris, vlist[i], vlist[j])) {
        /* prioritize edges in original tessellation */
        et[n].len += 100; /*ASCENDING ORDER*/
      }
      else {
        if (et[n].used == USED_IN_TESSELLATION) {
          continue;
        }
        et[n].used = USED_IN_ORIGINAL_TESSELLATION; /* to list
                                                       the original edges */
        nes++;
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

  ROMP_SCOPE_begin
  /* find and discard all edges that intersect one that is already in the
     tessellation.
  */
  for (ndiscarded = i = 0; i < nedges; i++) {
    if (et[i].used != USED_IN_TESSELLATION) {
      continue;
    }

    for (j = i + 1; j < nedges; j++) {
      if (et[j].used == USED_IN_TESSELLATION) {
        continue;
      }

      if (edgesIntersect(mris_corrected, &et[i], &et[j])) {
        ndiscarded++;
        if (j < nedges - 1) {
          memmove(&et[j], &et[j + 1], (nedges - j - 1) * sizeof(EDGE));
        }

        nedges--;
        j--;
      }
    }
  }
  ROMP_SCOPE_end
  

  if (DIAG_VERBOSE_ON) fprintf(WHICH_OUTPUT, "%d of %d overlapping edges discarded\n", ndiscarded, nedges);

  /* sort the edge list by edge length */
  qsort(et, nedges, sizeof(EDGE), compare_edge_length);

  if (!n) /* should never happen */
  {
    return (NO_ERROR);
  }

  ROMP_SCOPE_begin
  
#if 0
  //modify initial ordering, just for fun...
  {
    int k,m,tmp,rd;

    fprintf(WHICH_OUTPUT,"generating random ordering\n");

    ordering=(int*)calloc(nedges,sizeof(int));
    for (j=0; j< nedges; j++)
    {
      ordering[j]=nedges-j-1;
    }

    for (m=0; m<11; m++)
      for (j=0; j< nedges; j++)
      {
        k=nint(randomNumber(0.0, (double)nedges-1)) ;

        tmp=ordering[j];
        ordering[j]=ordering[k];
        ordering[k]=tmp;
      }
  }
#else  // no modification
  ordering = (int *)calloc(nedges, sizeof(int));
  for (j = 0; j < nedges; j++) {
    ordering[j] = j;  // nedges-j-1;
  }
#endif

  /* list the edges used in the original tessellation */
  es = (ES *)malloc(nes * sizeof(ES));
  for (nes = i = 0; i < nedges; i++)
    if (et[i].used == USED_IN_ORIGINAL_TESSELLATION) {
      // et[i].used=0; //reset state
      es[nes].vno1 = et[i].vno1;
      es[nes].vno2 = et[i].vno2;

      es[nes].segment = -1;
      es[nes++].n = i;
    }

  // main part of the routine: the retessellation (using a specific method) !
  if (getenv("USE_GA_TOPOLOGY_CORRECTION") != NULL) {
    parms->search_mode = GENETIC_SEARCH;
  }
  if (getenv("USE_RANDOM_TOPOLOGY_CORRECTION") != NULL) {
    parms->search_mode = RANDOM_SEARCH;
  }

  ROMP_SCOPE_end
  
  switch (parms->search_mode) {
    case GENETIC_SEARCH:
      ROMP_SCOPE_begin
      mrisComputeOptimalRetessellation(mris,
                                       mris_corrected,
                                       mri,
                                       defect,
                                       vertex_trans,
                                       et,
                                       nedges,
                                       es,
                                       nes,
                                       h_k1,
                                       h_k2,
                                       mri_k1_k2,
                                       h_white,
                                       h_gray,
                                       h_border,
                                       h_grad,
                                       mri_gray_white,
                                       h_dot,
                                       parms);
      ROMP_SCOPE_end
      break;
    case RANDOM_SEARCH:
      ROMP_SCOPE_begin
      mrisComputeRandomRetessellation(mris,
                                      mris_corrected,
                                      mri,
                                      defect,
                                      vertex_trans,
                                      et,
                                      nedges,
                                      es,
                                      nes,
                                      h_k1,
                                      h_k2,
                                      mri_k1_k2,
                                      h_white,
                                      h_gray,
                                      h_border,
                                      h_grad,
                                      mri_gray_white,
                                      h_dot,
                                      parms);
      ROMP_SCOPE_end
      break;
    default:
      ROMP_SCOPE_begin
      parms->search_mode = GREEDY_SEARCH;
      mrisRetessellateDefect(mris, mris_corrected, defect, vertex_trans, et, nedges, ordering, NULL);
      ROMP_SCOPE_end
      break;
  }

  free(es);
  if (ordering) {
    free(ordering);
  }
  free(et);
  
  return (NO_ERROR);
}

#define NUM_TO_ADD_FROM_ONE_PARENT 1

static int mrisCrossoverDefectPatches(DEFECT_PATCH *dp1, DEFECT_PATCH *dp2, DEFECT_PATCH *dp_dst, EDGE_TABLE *etable)
{
  int i1, i2, *added, i, isrc, j, nadded;
  double p;
  DEFECT_PATCH *dp_src;

  added = (int *)calloc(dp1->nedges, sizeof(int));
  p = randomNumber(0.0, 1.0);
  if (p < 0.5) /* add from first defect */
  {
    dp_src = dp1;
  }
  else {
    dp_src = dp2;
  }

  for (nadded = isrc = i1 = i2 = i = 0; i < dp_dst->nedges; i++) {
    if (nadded >= NUM_TO_ADD_FROM_ONE_PARENT) {
      nadded = 0;

      if ((dp_src == dp1 && i2 < dp2->nedges) || i1 >= dp1->nedges) /* use dp2 */
      {
        dp_src = dp2;
        isrc = i2;
      }
      else if (i1 < dp1->nedges) {
        dp_src = dp1;
        isrc = i1;
      }
    }
    else /* keep adding from the same parent */
    {
      nadded++;
      if (dp_src == dp1) {
        isrc = i1;
      }
      else {
        isrc = i2;
      }
    }
    if (isrc >= dp_src->nedges) /* shouldn't happen */
    {
      i--;
      continue;
    }

    /* find the next index in the src dp that hasn't been added yet */
    while (added[dp_src->ordering[isrc]])
      if (++isrc >= dp_src->nedges) {
        break;
      }
    if (isrc >= dp_src->nedges) {
      i--; /* process this one again */
      continue;
    }
    if (dp_src->ordering[isrc] == Gdiag_no) {
      DiagBreak();
    }
    if (isrc == Gdiag_no) {
      DiagBreak();
    }

    if (dp_src == dp1) /* update source index for next iteration */
    {
      i1 = isrc + 1;
    }
    else {
      i2 = isrc + 1;
    }

    dp_dst->ordering[i] = dp_src->ordering[isrc];
    if (dp_dst->ordering[i] >= dp_dst->nedges) {
      DiagBreak();
    }

    added[dp_src->ordering[isrc]] = 1; /* make sure every edge
                                          is represented */
  }

  for (i = 0; i < dp_dst->nedges; i++) {
    if ((dp_dst->ordering[i] >= dp_dst->nedges) || (dp_dst->ordering[i] < 0)) {
      DiagBreak();
      dp_dst->ordering[i] = 0;
    }
    if (added[i] == 0) {
      for (j = 0; j < dp_dst->nedges; j++)
        if (dp_dst->ordering[j] < 0) /* nothing in this slot */
        {
          dp_dst->ordering[j] = i;
          added[i] = 1; /* make sure every edge is represented */
          break;
        }
    }
  }

  free(added);

  return (NO_ERROR);
}
#define NTRY 0
#if 1
static int mrisMutateDefectPatch(DEFECT_PATCH *dp, EDGE_TABLE *etable, double pmutation)
{
  int i, j, eti, etj, tmp, *dp_indices, ntry;
  double p;
  EDGE *e;

  //    fprintf(WHICH_OUTPUT,"m");
  dp_indices = (int *)calloc(dp->nedges, sizeof(int));
  for (i = 0; i < dp->nedges; i++) {
    dp_indices[dp->ordering[i]] = i;
  }

  for (i = 0; i < dp->nedges; i++) {
    p = randomNumber(0.0, 1.0);
    eti = dp->ordering[i];

    if (p < pmutation) {
      /* mutation */
      if (etable->use_overlap && etable->noverlap[eti] > 0) {
        if (etable->flags[eti] & ET_OVERLAP_LIST_INCOMPLETE) /* swap any
                                                                two */
        {
          ntry = 0;
          j = (int)randomNumber(0.0, dp->nedges - .1);
          while (ntry < NTRY) {
            e = &etable->edges[dp->ordering[j]]; /*potential new edge */
            if (e->used != USED_IN_ORIGINAL_TESSELLATION) {
              ntry++;
            }
            else {
              break;
            }
            j = (int)randomNumber(0.0, dp->nedges - .1);
          }
          tmp = dp->ordering[i];
          dp->ordering[i] = dp->ordering[j];
          dp->ordering[j] = tmp;
        }
        else /* swap two edges that intersect */
        {
          ntry = 0;
          j = (int)randomNumber(0.0, etable->noverlap[eti] - 0.0001);
          etj = etable->overlapping_edges[eti][j]; /* index of jth
                                                      overlapping edge */
          j = dp_indices[etj];                     /* find where it is in this
                                                      defect patch ordering */
          while (ntry < NTRY) {
            e = &etable->edges[dp->ordering[j]]; /*potential new edge */
            if (e->used != USED_IN_ORIGINAL_TESSELLATION) {
              ntry++;
            }
            else {
              break;
            }
            j = (int)randomNumber(0.0, etable->noverlap[eti] - 0.0001);
            etj = etable->overlapping_edges[eti][j]; /* index of
                                                        jth overlapping
                                                        edge */
            j = dp_indices[etj];                     /* find where it is in this
                                                        defect patch ordering */
          }
          tmp = dp->ordering[i];
          dp->ordering[i] = dp->ordering[j];
          dp->ordering[j] = tmp;
          dp_indices[dp->ordering[i]] = i;
          dp_indices[dp->ordering[j]] = j;
        }
      }
      else {
        /* swap any two */
        ntry = 0;
        j = (int)randomNumber(0.0, dp->nedges - .1);
        while (ntry < NTRY) {
          e = &etable->edges[dp->ordering[j]]; /*potential new edge */
          if (e->used != USED_IN_ORIGINAL_TESSELLATION) {
            ntry++;
          }
          else {
            break;
          }
          j = (int)randomNumber(0.0, dp->nedges - .1);
        }
        tmp = dp->ordering[i];
        dp->ordering[i] = dp->ordering[j];
        dp->ordering[j] = tmp;
      }
    }
  }
  //  fprintf(WHICH_OUTPUT,".");
  free(dp_indices);
  return (NO_ERROR);
}
#else
static int mrisMutateDefectPatch(DEFECT_PATCH *dp, EDGE_TABLE *etable, double pmutation)
{
  int i, j, eti, etj, tmp, *dp_indices;
  double p;

  dp_indices = (int *)calloc(dp->nedges, sizeof(int));
  for (i = 0; i < dp->nedges; i++) {
    dp_indices[dp->ordering[i]] = i;
  }

  for (i = 0; i < dp->nedges; i++) {
    p = randomNumber(0.0, 1.0);
    eti = dp->ordering[i];
    if (p < pmutation && etable->noverlap[eti] > 0) {
      if (etable->flags[eti] & ET_OVERLAP_LIST_INCOMPLETE) /* swap any two */
      {
        j = (int)randomNumber(0.0, dp->nedges - .1);
        tmp = dp->ordering[i];
        dp->ordering[i] = dp->ordering[j];
        dp->ordering[j] = tmp;
      }
      else /* swap two edges that intersect */
      {
        j = (int)randomNumber(0.0, etable->noverlap[eti] - 0.0001);
        etj = etable->overlapping_edges[eti][j]; /* index of jth
                                                    overlapping edge */
        j = dp_indices[etj];                     /* find where it is in this
                                                    defect patch ordering */

        tmp = dp->ordering[i];
        dp->ordering[i] = dp->ordering[j];
        dp->ordering[j] = tmp;
        dp_indices[dp->ordering[i]] = i;
        dp_indices[dp->ordering[j]] = j;
      }
    }
  }

  free(dp_indices);
  return (NO_ERROR);
}
#endif

static int mrisCopyDefectPatch(DEFECT_PATCH *dp_src, DEFECT_PATCH *dp_dst)
{
  int i;

  dp_dst->etable = dp_src->etable;
  dp_dst->defect = dp_src->defect;
  dp_dst->nedges = dp_src->nedges;
  dp_dst->fitness = dp_src->fitness;
  dp_dst->rank = dp_src->rank;
  dp_dst->retessellation_mode = dp_src->retessellation_mode;

  for (i = 0; i < dp_src->nedges; i++) {
    dp_dst->ordering[i] = dp_src->ordering[i];
  }

  return (NO_ERROR);
}

static int defectPatchRank(DEFECT_PATCH *dps, int index, int npatches)
{
  int i, rank = 0;

  for (i = 0; i < npatches; i++) {
    if (i == index) {
      continue;
    }

    if (FEQUAL(dps[i].fitness, dps[index].fitness)) {
      if (index > i) {
        rank++;
      }
    }
    else if (dps[i].fitness > dps[index].fitness) {
      rank++;
    }
  }
  return (rank);
}

#define MAX_EDGES 1000

static int tessellatePatch(MRI *mri,
                           MRI_SURFACE *mris,
                           MRI_SURFACE *mris_corrected,
                           DEFECT *defect,
                           int *vertex_trans,
                           EDGE *et,
                           int nedges,
                           int *ordering,
                           EDGE_TABLE *etable,
                           TOPOLOGY_PARMS *parms)
{
  int i, k, vni;
  DVS *dvs;
  DP dp;
  EDGE_TABLE *new_et = NULL;

  if (DIAG_VERBOSE_ON) {
    fprintf(WHICH_OUTPUT, "tessellating patch....\n");
  }

  /* allocation */
  dvs = mrisRecordVertexState(mris_corrected, defect, vertex_trans);

  if (etable) {
    new_et = etable;
  }
  else {
    new_et = (EDGE_TABLE *)calloc(1, sizeof(EDGE_TABLE));
    new_et->nedges = nedges;
    new_et->edges = et;
    new_et->use_overlap = 0;
  }

  if (parms->retessellation_mode) {
    dp.retessellation_mode = USE_SOME_VERTICES;
  }
  else {
    dp.retessellation_mode = USE_ALL_VERTICES;
  }
  dp.nedges = nedges;
  dp.defect = defect;
  dp.etable = new_et;
  dp.ordering = ordering;
  dp.mri = mri;

  /* compute the final tessellation */
  retessellateDefect(mris, mris_corrected, dvs, &dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* free */
  mrisFreeDefectVertexState(dvs);

  /* detect the new set of faces */
  detectDefectFaces(mris_corrected, &dp);

  /* orient the patch faces */
  orientDefectFaces(mris_corrected, &dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* smooth original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, &dp, parms->smooth, parms->match);

  /* should free the tessellated patch structure */
  TPfree(&dp.tp);

  if (!etable) {
    free(new_et);
  }

  /* discard the vertices that are not used in the
     final tessellation and set marks to zero!!!*/
  for (i = 0; i < dp.nedges; i++)
    if (dp.etable->edges[i].used == USED_IN_NEW_TESSELLATION || dp.etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
      mris_corrected->vertices[dp.etable->edges[i].vno1].marked = FINAL_VERTEX;
      mris_corrected->vertices[dp.etable->edges[i].vno2].marked = FINAL_VERTEX;
    }

  for (k = i = 0; i < dp.defect->nvertices; i++) {
    vni = vertex_trans[dp.defect->vertices[i]];
    if (vni < 0) {
      continue;
    }
    if (mris_corrected->vertices[vni].marked != FINAL_VERTEX) {
      if (dp.defect->status[i] != DISCARD_VERTEX) {
        k++;
      }
      dp.defect->status[i] = DISCARD_VERTEX;
      mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].ripflag = 1;
    }
    mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].marked = 0;
  }

  for (i = 0; i < dp.defect->nborder; i++) {
    mris_corrected->vertices[vertex_trans[dp.defect->border[i]]].marked = 0;
  }

  if (DIAG_VERBOSE_ON) {
    fprintf(WHICH_OUTPUT, "done\n");
  }

  return NO_ERROR;
}

#define USE_SCALING 0

int mrisCountIntersectingFaces(MRIS *mris, int *flist, int nfaces)
{
  int n, m, i, count, intersect, j;
  int vn[3];
  int un[3];
  int vertex_in_common;
  double v0[3], v1[3], v2[3], u0[3], u1[3], u2[3];
  double d0, d1, d2, d, scale, SCALE_FACTOR;
  FACE *f1, *f2;

  SCALE_FACTOR = 100.0f;

  for (count = n = 0; n < nfaces; n++) {
    f1 = &mris->faces[flist[n]];
    /* fill vertices of 1st triangle */
    vn[0] = f1->v[0];
    vn[1] = f1->v[1];
    vn[2] = f1->v[2];
    v0[0] = (double)mris->vertices[f1->v[0]].origx;
    v0[1] = (double)mris->vertices[f1->v[0]].origy;
    v0[2] = (double)mris->vertices[f1->v[0]].origz;
    v1[0] = (double)mris->vertices[f1->v[1]].origx;
    v1[1] = (double)mris->vertices[f1->v[1]].origy;
    v1[2] = (double)mris->vertices[f1->v[1]].origz;
    v2[0] = (double)mris->vertices[f1->v[2]].origx;
    v2[1] = (double)mris->vertices[f1->v[2]].origy;
    v2[2] = (double)mris->vertices[f1->v[2]].origz;
    d0 = SQR(v1[0] - v0[0]) + SQR(v1[1] - v0[1]) + SQR(v1[2] - v0[2]);
    d1 = SQR(v1[0] - v2[0]) + SQR(v1[1] - v2[1]) + SQR(v1[2] - v2[2]);
    d2 = SQR(v0[0] - v2[0]) + SQR(v0[1] - v2[1]) + SQR(v0[2] - v2[2]);
    d = sqrt(MIN(d0, MIN(d1, d2)));
    /* scaling */
    if (FZERO(d)) {
      continue;
    }
    scale = SCALE_FACTOR / d;
#if USE_SCALING
    for (i = 0; i < 3; i++) {
      v0[i] *= scale;
      v1[i] *= scale;
      v2[i] *= scale;
    }
#endif
    for (m = n + 1; m < nfaces; m++) {
      f2 = &mris->faces[flist[m]];
      un[0] = f2->v[0];
      un[1] = f2->v[1];
      un[2] = f2->v[2];
      /* count the number of common vertices */
      vertex_in_common = 0;
      for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
          if (vn[i] == un[j]) {
            vertex_in_common++;
            break;
          }
      }

      if (vertex_in_common > 0) {
        continue;
      }

      // fprintf(WHICH_OUTPUT,"vertex in common = %d \n",vertex_in_common);

      /* fill vertices of 2nd triangle */
      u0[0] = (double)mris->vertices[f2->v[0]].origx;
      u0[1] = (double)mris->vertices[f2->v[0]].origy;
      u0[2] = (double)mris->vertices[f2->v[0]].origz;
      u1[0] = (double)mris->vertices[f2->v[1]].origx;
      u1[1] = (double)mris->vertices[f2->v[1]].origy;
      u1[2] = (double)mris->vertices[f2->v[1]].origz;
      u2[0] = (double)mris->vertices[f2->v[2]].origx;
      u2[1] = (double)mris->vertices[f2->v[2]].origy;
      u2[2] = (double)mris->vertices[f2->v[2]].origz;

      d0 = SQR(u1[0] - u0[0]) + SQR(u1[1] - u0[1]) + SQR(u1[2] - u0[2]);
      d1 = SQR(u1[0] - u2[0]) + SQR(u1[1] - u2[1]) + SQR(u1[2] - u2[2]);
      d2 = SQR(u0[0] - u2[0]) + SQR(u0[1] - u2[1]) + SQR(u0[2] - u2[2]);
      d = sqrt(MIN(d0, MIN(d1, d2)));
      if (FZERO(d)) {
        continue;
      }
/* scaling */
#if USE_SCALING
      for (i = 0; i < 3; i++) {
        u0[i] *= scale;
        u1[i] *= scale;
        u2[i] *= scale;
      }
#endif
      intersect = tri_tri_intersect(v0, v1, v2, u0, u1, u2);

      if (intersect) {
#if 0
        fprintf(stderr,"face %d intersect face %d\n",
                flist[n],flist[m]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[n],f1->v[0],v0[0],v0[1],v0[2]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[n],f1->v[1],v1[0],v1[1],v1[2]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[n],f1->v[2],v2[0],v2[1],v2[2]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[m],f2->v[0],u0[0],u0[1],u0[2]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[m],f2->v[1],u1[0],u1[1],u1[2]);
        fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
                flist[m],f2->v[2],u2[0],u2[1],u2[2]);
        fprintf(stderr,"\nXXXXXXXXXX\n");
#endif
        count++;
        break;
      }
    }
  }
#if 0
  fprintf(stderr,"\n\nYYYYYYYYYYYYYYYYYYYYY LIST FACES \n");
  for (n=0; n < nfaces ; n++)
  {
    f1 = &mris->faces[flist[n]];
    /* fill vertices of 1st triangle */
    v0[0] = (double)mris->vertices[f1->v[0]].origx ;
    v0[1] = (double)mris->vertices[f1->v[0]].origy ;
    v0[2] = (double)mris->vertices[f1->v[0]].origz ;
    v1[0] = (double)mris->vertices[f1->v[1]].origx ;
    v1[1] = (double)mris->vertices[f1->v[1]].origy ;
    v1[2] = (double)mris->vertices[f1->v[1]].origz ;
    v2[0] = (double)mris->vertices[f1->v[2]].origx ;
    v2[1] = (double)mris->vertices[f1->v[2]].origy ;
    v2[2] = (double)mris->vertices[f1->v[2]].origz ;

    fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
            flist[n],f1->v[0],v0[0],v0[1],v0[2]);
    fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
            flist[n],f1->v[1],v1[0],v1[1],v1[2]);
    fprintf(stderr,"face %d : vertex %d [ %f , %f , %f ]\n",
            flist[n],f1->v[2],v2[0],v2[1],v2[2]);
  }
#endif

  return count;
}

#define SAVE_FIT_VALS 0
#if SAVE_FIT_VALS
static float fitness_values[11000];
static float best_values[11000];
#endif

static int mrisComputeOptimalRetessellation_wkr(MRI_SURFACE *mris,
                                            MRI_SURFACE *mris_corrected,
                                            MRI *mri,
                                            DEFECT *defect,
                                            int *vertex_trans,
                                            EDGE *et,
                                            int nedges,
                                            ES *es,
                                            int nes,
                                            HISTOGRAM *h_k1,
                                            HISTOGRAM *h_k2,
                                            MRI *mri_k1_k2,
                                            HISTOGRAM *h_white,
                                            HISTOGRAM *h_gray,
                                            HISTOGRAM *h_border,
                                            HISTOGRAM *h_grad,
                                            MRI *mri_gray_white,
                                            HISTOGRAM *h_dot,
                                            TOPOLOGY_PARMS *parms);


static int mrisComputeOptimalRetessellation(MRI_SURFACE *mris,
                                            MRI_SURFACE *mris_corrected,
                                            MRI *mri,
                                            DEFECT *defect,
                                            int *vertex_trans,
                                            EDGE *et,
                                            int nedges,
                                            ES *es,
                                            int nes,
                                            HISTOGRAM *h_k1,
                                            HISTOGRAM *h_k2,
                                            MRI *mri_k1_k2,
                                            HISTOGRAM *h_white,
                                            HISTOGRAM *h_gray,
                                            HISTOGRAM *h_border,
                                            HISTOGRAM *h_grad,
                                            MRI *mri_gray_white,
                                            HISTOGRAM *h_dot,
                                            TOPOLOGY_PARMS *parms)
{
    int result;
    ROMP_SCOPE_begin
    result = mrisComputeOptimalRetessellation_wkr(mris,
                                            mris_corrected,
                                            mri,
                                            defect,
                                            vertex_trans,
                                            et,
                                            nedges,
                                            es,
                                            nes,
                                            h_k1,
                                            h_k2,
                                            mri_k1_k2,
                                            h_white,
                                            h_gray,
                                            h_border,
                                            h_grad,
                                            mri_gray_white,
                                            h_dot,
                                            parms);
    ROMP_SCOPE_end
    return result;
}

static NOINLINE int mrisComputeOptimalRetessellation_wkr(MRI_SURFACE *mris,
                                            MRI_SURFACE *mris_corrected,
                                            MRI *mri,
                                            DEFECT *defect,
                                            int *vertex_trans,
                                            EDGE *et,
                                            int nedges,
                                            ES *es,
                                            int nes,
                                            HISTOGRAM *h_k1,
                                            HISTOGRAM *h_k2,
                                            MRI *mri_k1_k2,
                                            HISTOGRAM *h_white,
                                            HISTOGRAM *h_gray,
                                            HISTOGRAM *h_border,
                                            HISTOGRAM *h_grad,
                                            MRI *mri_gray_white,
                                            HISTOGRAM *h_dot,
                                            TOPOLOGY_PARMS *parms)
{
  DEFECT_VERTEX_STATE *dvs;
  DEFECT_PATCH dps1[MAX_PATCHES], dps2[MAX_PATCHES], *dps, *dp, *dps_next_generation;
  int i, best_i, j, g, nselected, nreplacements, rank, nunchanged = 0, nelite, ncrossovers, k, l, noverlap;
  int *overlap, ngenerations, nbests, last_euthanasia, nremovedvertices, nfinalvertices;
  double fitness, best_fitness, last_best, fitness_mean, fitness_sigma, fitness_norm, pfitness, two_sigma_sq,
      last_fitness;
  static int dno = 0;     /* for debugging */
  static int nmovies = 1; /* for making movies :
                                 0 is left for the original surface*/
  EDGE_TABLE etable;
  int max_patches = MAX_PATCHES, ranks[MAX_PATCHES], next_gen_index, selected[MAX_PATCHES], nzero, sno = 0, max_edges,
      debug_patch_n = -1, nbest = 0;
  MRI *mri_defect, *mri_defect_white, *mri_defect_gray, *mri_defect_sign;
  char fname[500];
  SEGMENTATION *segmentation;
  RP rp;
  int number_of_patches, nbestpatch;
  int ncross_overs, ntotalcross_overs, ntotalmutations, nmutations;
  int nintersections;
  static int first_time = 1;

  nbestpatch = number_of_patches = 0;
  ncross_overs = nmutations = 0;
  ntotalcross_overs = ntotalmutations = 0;

  if (first_time) {
    char *cp;

    if ((cp = getenv("FS_QCURV")) != NULL) {
      parms->l_qcurv = atof(cp);
      fprintf(WHICH_OUTPUT, "setting qcurv = %2.3f\n", l_qcurv);
    }
    if ((cp = getenv("FS_CURV")) != NULL) {
      parms->l_curv = atof(cp);
      fprintf(WHICH_OUTPUT, "setting curv = %2.3f\n", l_curv);
    }
    if ((cp = getenv("FS_MRI")) != NULL) {
      parms->l_mri = atof(cp);
      fprintf(WHICH_OUTPUT, "setting mri = %2.3f\n", l_mri);
    }
    if ((cp = getenv("FS_UNMRI")) != NULL) {
      parms->l_unmri = atof(cp);
      fprintf(WHICH_OUTPUT, "setting unmri = %2.3f\n", l_unmri);
    }
    first_time = 0;
  }

  max_patches = parms->max_patches;
  max_unchanged = parms->max_unchanged;
  max_edges = MAX_EDGES;

  if (dno == Gdiag_no) {
    DiagBreak();
  }

  if (getenv("FS_DEBUG_PATCH") != NULL) {
    int debug_patch = atoi(getenv("FS_DEBUG_PATCH"));
    if (debug_patch != dno) {
      max_patches = 0;
    }
    else {
      if (getenv("FS_DEBUG_PATCH_N") != NULL) {
        debug_patch_n = atoi(getenv("FS_DEBUG_PATCH_N"));
        fprintf(WHICH_OUTPUT, "terminating after %dth best tessellation\n", debug_patch_n);
      }
    }
  }

#if 0

  dno++ ;  /* for debugging */
  tessellatePatch(mri,mris, mris_corrected, defect,
                  vertex_trans, et, nedges, NULL, NULL,parms);
  MRIScheckIsPolyhedron(mris_corrected, __FILE__,__LINE__);

  return(NO_ERROR) ;

#endif

  if (!max_patches) {
    dno++; /* for debugging */
    // mrisRetessellateDefect(mris, mris_corrected,
    // defect, vertex_trans, et, nedges, NULL, NULL) ;

    tessellatePatch(mri, mris, mris_corrected, defect, vertex_trans, et, nedges, NULL, NULL, parms);

    return (NO_ERROR);
  }

  dno++; /* for debugging */

  if (nedges > 200000) {
#if 1
    printf("An extra large defect has been detected...\n");
    printf("This often happens because cerebellum or dura has not been removed from wm.mgz.\n");
    printf("This may cause recon-all to run very slowly or crash.\n");
    printf("if so, see https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TopologicalDefect_freeview\n");
    max_unchanged = MIN(max_unchanged, 1);
    max_patches = MAX(MIN(max_patches, 3), 1);
    max_edges = MIN(max_edges, 25);
#else
    // add some code here to select a good tessellation (ordering) FLO
    mrisRetessellateDefect(mris, mris_corrected, defect, vertex_trans, et, nedges, NULL, NULL);
    //    tessellatePatch(mri,mris, mris_corrected,
    //                    defect,vertex_trans, et, nedges, NULL, NULL,parms);

    return (NO_ERROR);
#endif
  }
  else if (nedges > 100000) {
    printf("A large defect has been detected...\n");
    printf("This often happens because cerebellum or dura has not been removed from wm.mgz.\n");
    printf("This may cause recon-all to run very slowly or crash.\n");
    printf("if so, see https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/TopologicalDefect_freeview\n");
    max_unchanged = MIN(max_unchanged, 1);
    max_patches = MAX(MIN(max_patches, 10), 1);
    max_edges = MIN(max_edges, 100);
  }
  else if (nedges > 50000) {
    max_patches = MAX(max_patches / 2, 1);
    max_edges = max_edges / 5;
    max_unchanged = max_unchanged / 5;
  }

  ROMP_SCOPE_begin

  etable.use_overlap = parms->edge_table;
  etable.nedges = nedges;
  etable.edges = (EDGE *)calloc(nedges, sizeof(EDGE));
  memmove(etable.edges, et, nedges * sizeof(EDGE));

  if (etable.use_overlap) {
    etable.overlapping_edges = (int **)calloc(nedges, sizeof(int *));
    etable.noverlap = (int *)calloc(nedges, sizeof(int));
    etable.flags = (unsigned char *)calloc(nedges, sizeof(unsigned char));
    overlap = (int *)calloc(nedges, sizeof(int));
    if (!etable.edges || !etable.overlapping_edges || !etable.noverlap || !overlap)
      ErrorExit(ERROR_NOMEMORY,
                "mrisComputeOptimalRetessellation: Excessive "
                "topologic defect encountered: could not allocate %d "
                "edge table",
                nedges);

    for (nzero = i = 0; i < nedges; i++) /* compute overlapping
                                                        for each edge */
    {
      if (nedges > 50000 && !(i % 25000)) {
        fprintf(WHICH_OUTPUT, "%d of %d edges processed\n", i, nedges);
      }
      etable.noverlap[i] = 0;
      for (noverlap = j = 0; j < nedges; j++) {
        if (j == i) {
          continue;
        }
        if (edgesIntersect(mris_corrected, &et[i], &et[j])) {
          overlap[noverlap] = j;
          noverlap++;
        }
        if (noverlap > MAX_EDGES) {
          break;
        }
      }
      if (noverlap > 0) {
        if (noverlap > MAX_EDGES) {
          etable.noverlap[i] = MAX_EDGES;
          etable.flags[i] |= ET_OVERLAP_LIST_INCOMPLETE;
        }
        else {
          etable.noverlap[i] = noverlap;
        }

        etable.overlapping_edges[i] = (int *)calloc(etable.noverlap[i], sizeof(int));
        if (!etable.overlapping_edges[i])
          ErrorExit(ERROR_NOMEMORY,
                    "mrisComputeOptimalRetessellation: Excessive "
                    "topologic defect encountered: could not allocate "
                    "overlap list %d "
                    "with %d elts",
                    i,
                    etable.noverlap[i]);
        memmove(etable.overlapping_edges[i], overlap, etable.noverlap[i] * sizeof(int));
      }
      else {
        nzero++;
      }
    }

    free(overlap);
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  /* allocate the volume constituted by the potential edges */
  mri_defect = mri_defect_white = mri_defect_gray = mri_defect_sign = NULL;
  if (!FZERO(parms->l_unmri)) {
    mri_defect = mriDefectVolume(mris_corrected, &etable, parms);
    mri_defect_white = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_gray = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_sign = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    defectVolumeLikelihood(mri,
                           mri_defect,
                           mri_defect_white,
                           mri_defect_gray,
                           h_white,
                           h_gray,
                           defect->white_mean,
                           defect->gray_mean,
                           0,
                           0);
  };

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  if ((!FZERO(parms->l_unmri)) && parms->save_fname &&
      (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
    sprintf(fname, "%s/white_%d.mgh", parms->save_fname, defect->defect_number);
    MRIwrite(mri_defect_white, fname);
    sprintf(fname, "%s/gray_%d.mgh", parms->save_fname, defect->defect_number);
    MRIwrite(mri_defect_gray, fname);
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  dvs = mrisRecordVertexState(mris_corrected, defect, vertex_trans);
  dps = dps1;

  ngenerations = 0;
  last_euthanasia = -1;
  number_of_patches = 0;
  nremovedvertices = 0;
  nfinalvertices = 0;

  /* generate Random Patch */
  rp.best_ordering = (int *)malloc(nedges * sizeof(int));
  rp.status = (char *)malloc(defect->nvertices * sizeof(char));
  memmove(rp.status, defect->status, defect->nvertices * sizeof(char));
  rp.nused = (int *)calloc(defect->nvertices, sizeof(int));
  rp.vertex_fitness = (float *)calloc(defect->nvertices, sizeof(float));

  nbests = 0;

  ROMP_SCOPE_end

  ComputeDefectContext computeDefectContext;

  ROMP_SCOPE_begin

    constructComputeDefectContext(&computeDefectContext);

  /* generate initial population of patches */
  if (parms->initial_selection) {
    /* segment overlapping edges into clusters */
    segmentation = segmentIntersectingEdges(mris_corrected, defect, vertex_trans, es, &nes);
    if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number)))
      saveSegmentation(mris, mris_corrected, defect, vertex_trans, es, nes, parms->save_fname);

    /* generate initial population of patches */
    best_fitness = -1000000;
    best_i = 0;
    for (i = 0; i < max_patches; i++) {
      dp = &dps2[i];

      if (parms->retessellation_mode) {
        dp->retessellation_mode = USE_SOME_VERTICES;
      }
      else {
        dp->retessellation_mode = USE_ALL_VERTICES;
      }

      dp->nedges = nedges;
      dp->defect = defect;
      dp->etable = &etable;
      dp->ordering = (int *)calloc(nedges, sizeof(int));
      if (!dp->ordering) ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices", i, nedges);
      for (j = 0; j < nedges; j++) {
        dp->ordering[j] = j;
      }

      dp->mri_defect = mri_defect;
      dp->mri_defect_white = mri_defect_white;
      dp->mri_defect_gray = mri_defect_gray;
      dp->mri_defect_sign = mri_defect_sign;

      dp->mri = mri;

      dp = &dps1[i];

      if (parms->retessellation_mode) {
        dp->retessellation_mode = USE_SOME_VERTICES;
      }
      else {
        dp->retessellation_mode = USE_ALL_VERTICES;
      }

      dp->nedges = nedges;
      dp->defect = defect;
      dp->etable = &etable;
      dp->ordering = (int *)calloc(nedges, sizeof(int));
      if (!dp->ordering) ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices", i, nedges);
      for (j = 0; j < nedges; j++) {
        dp->ordering[j] = j;  // nedges-j-1 ;
      }
      /* initial in same order -
         will change later */

      dp->mri_defect = mri_defect;
      dp->mri_defect_white = mri_defect_white;
      dp->mri_defect_gray = mri_defect_gray;
      dp->mri_defect_sign = mri_defect_sign;

      dp->mri = mri;

      /* generate ordering from edge segmentation */
      generateOrdering(dp, segmentation, i);

      fitness = mrisDefectPatchFitness(&computeDefectContext,
                                       mris,
                                       mris_corrected,
                                       mri,
                                       dp,
                                       vertex_trans,
                                       dvs,
                                       &rp,
                                       h_k1,
                                       h_k2,
                                       mri_k1_k2,
                                       h_white,
                                       h_gray,
                                       h_border,
                                       h_grad,
                                       mri_gray_white,
                                       h_dot,
                                       parms);

#if SAVE_FIT_VALS
      fitness_values[number_of_patches] = fitness;
      if (number_of_patches)
        best_values[number_of_patches] = MAX(best_values[number_of_patches - 1], fitness);
      else {
        best_values[number_of_patches] = fitness;
      }
#endif
      number_of_patches++;

      if (parms->verbose == VERBOSE_MODE_LOW)
        fprintf(WHICH_OUTPUT, "for the patch #%d, we have fitness = %f \n", i, fitness);

      /* saving the initial selection */
      if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
        sprintf(fname, "%s/rh.defect_%d_select%d", parms->save_fname, defect->defect_number, i);
        savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        if (parms->movie) {
          sprintf(fname, "%s/rh.defect_%d_movie_%d", parms->save_fname, defect->defect_number, nmovies++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }
      }

      if (!i)  // fisrt patch
      {
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));
        dp->defect->initial_face_ll = dp->tp.face_ll;
        dp->defect->initial_vertex_ll = dp->tp.vertex_ll;
        dp->defect->initial_curv_ll = dp->tp.curv_ll;
        dp->defect->initial_qcurv_ll = dp->tp.qcurv_ll;
        dp->defect->initial_mri_ll = dp->tp.mri_ll;
        dp->defect->initial_unmri_ll = dp->tp.unmri_ll;

        if (parms->verbose == VERBOSE_MODE_LOW) {
          fprintf(WHICH_OUTPUT, "initial defect\n");
          printDefectStatistics(dp);
        }
        best_fitness = fitness;
        best_i = 0;
        // saving first patch
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }

        if (++nbest == debug_patch_n) {
          goto debug_use_this_patch;
        }
      }

      if (fitness > best_fitness) {
        best_fitness = fitness;
        best_i = i;
        if (parms->verbose > VERBOSE_MODE_DEFAULT)
          fprintf(WHICH_OUTPUT, "new optimal fitness found at %d: %2.4f\n", i, fitness);

        nfinalvertices = nremovedvertices;
        nbestpatch = number_of_patches;

        rp.best_fitness = best_fitness;
        /* save ordering*/
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        /* save current status of vertices */
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

        if (parms->verbose == VERBOSE_MODE_LOW) {
          printDefectStatistics(dp);
        }
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname, "%s/rh.defect_%d_best_%d_%d", parms->save_fname, defect->defect_number, ngenerations, i);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }

        if (++nbest == debug_patch_n) {
          goto debug_use_this_patch;
        }
      }
    }
    if (segmentation) {
      SEGMENTATIONfree(&segmentation);
    }
  }
  else {
    /* generate initial population of patches */
    best_fitness = -1000000;
    best_i = 0;
    for (i = 0; i < max_patches; i++) {
      dp = &dps2[i];

      if (parms->retessellation_mode) {
        dp->retessellation_mode = USE_SOME_VERTICES;
      }
      else {
        dp->retessellation_mode = USE_ALL_VERTICES;
      }

      dp->nedges = nedges;
      dp->defect = defect;
      dp->etable = &etable;
      dp->ordering = (int *)calloc(nedges, sizeof(int));
      if (!dp->ordering) ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices", i, nedges);
      for (j = 0; j < nedges; j++) {
        dp->ordering[j] = j;
      } /* initial in same order -
                                               will change later */

      dp->mri_defect = mri_defect;
      dp->mri_defect_white = mri_defect_white;
      dp->mri_defect_gray = mri_defect_gray;
      dp->mri_defect_sign = mri_defect_sign;

      dp->mri = mri;

      dp = &dps1[i];

      if (parms->retessellation_mode) {
        dp->retessellation_mode = USE_SOME_VERTICES;
      }
      else {
        dp->retessellation_mode = USE_ALL_VERTICES;
      }

      dp->nedges = nedges;
      dp->defect = defect;
      dp->etable = &etable;
      dp->ordering = (int *)calloc(nedges, sizeof(int));
      if (!dp->ordering) ErrorExit(ERROR_NOMEMORY, "could not allocate %dth defect patch with %d indices", i, nedges);
      for (j = 0; j < nedges; j++) {
        dp->ordering[j] = j;
      } /* initial in same order -
                                               will change later */

      dp->mri_defect = mri_defect;
      dp->mri_defect_white = mri_defect_white;
      dp->mri_defect_gray = mri_defect_gray;
      dp->mri_defect_sign = mri_defect_sign;

      dp->mri = mri;

      if (i) /* first one is in same order as original edge table */
      {
        mrisMutateDefectPatch(dp, &etable, MUTATION_PCT_INIT);
      }

      fitness = mrisDefectPatchFitness(&computeDefectContext,
                                       mris,
                                       mris_corrected,
                                       mri,
                                       dp,
                                       vertex_trans,
                                       dvs,
                                       &rp,
                                       h_k1,
                                       h_k2,
                                       mri_k1_k2,
                                       h_white,
                                       h_gray,
                                       h_border,
                                       h_grad,
                                       mri_gray_white,
                                       h_dot,
                                       parms);
#if SAVE_FIT_VALS
      fitness_values[number_of_patches] = fitness;
      if (number_of_patches)
        best_values[number_of_patches] = MAX(best_values[number_of_patches - 1], fitness);
      else {
        best_values[number_of_patches] = fitness;
      }
#endif
      number_of_patches++;

      if (i == 0 && Gdiag & 0x1000000) {
        int i;
        char fname[STRLEN];
        int req = snprintf(fname, STRLEN, "%s_defect%d_%03d", mris->fname.data(), dno - 1, sno++); 
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        dp = &dps[best_i];
        mrisRetessellateDefect(
            mris, mris_corrected, dp->defect, vertex_trans, dp->etable->edges, dp->nedges, dp->ordering, dp->etable);
        MRISsaveVertexPositions(mris_corrected, TMP_VERTICES);
        MRISrestoreVertexPositions(mris_corrected, ORIGINAL_VERTICES);
        fprintf(WHICH_OUTPUT, "writing surface snapshow to %s...\n", fname);
        MRISwrite(mris_corrected, fname);
        MRISrestoreVertexPositions(mris_corrected, TMP_VERTICES);
        mrisRestoreVertexState(mris_corrected, dvs);
        /* reset the edges to the unused state
           (unless they were in the original tessellation */
        for (i = 0; i < dp->nedges; i++)
          if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION) {
            dp->etable->edges[i].used = NOT_USED;
          }
      }

      /* saving the initial selection */
      if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
        sprintf(fname, "%s/rh.defect_%d_select%d", parms->save_fname, defect->defect_number, i);
        savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        if (parms->movie) {
          sprintf(fname, "%s/rh.defect_%d_movie_%d", parms->save_fname, defect->defect_number, nmovies++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }
      }

      if (!i) {
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

        dp->defect->initial_face_ll = dp->tp.face_ll;
        dp->defect->initial_vertex_ll = dp->tp.vertex_ll;
        dp->defect->initial_curv_ll = dp->tp.curv_ll;
        dp->defect->initial_qcurv_ll = dp->tp.qcurv_ll;
        dp->defect->initial_mri_ll = dp->tp.mri_ll;
        dp->defect->initial_unmri_ll = dp->tp.unmri_ll;
        if (parms->verbose == VERBOSE_MODE_LOW) {
          fprintf(WHICH_OUTPUT,
                  "defect %d: initial fitness = %2.4e, "
                  "nvertices=%d, nedges=%d, max patches=%d\n",
                  dno - 1,
                  fitness,
                  defect->nvertices,
                  nedges,
                  max_patches);
          printDefectStatistics(dp);
        }
        best_fitness = fitness;
        best_i = 0;

        // saving first patch
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }
        if (++nbest == debug_patch_n) {
          goto debug_use_this_patch;
        }
      }

      if (fitness > best_fitness) {
        best_fitness = fitness;
        best_i = i;
        if (parms->verbose > VERBOSE_MODE_DEFAULT)
          fprintf(WHICH_OUTPUT,
                  "new optimal fitness found at %d: "
                  "%2.4f\n",
                  i,
                  fitness);

        nfinalvertices = nremovedvertices;
        nbestpatch = number_of_patches;

        rp.best_fitness = best_fitness;
        /* save ordering*/
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        /* save current status of vertices */
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

        if (parms->verbose == VERBOSE_MODE_LOW) {
          printDefectStatistics(dp);
        }
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname, "%s/rh.defect_%d_best_%d_%d", parms->save_fname, defect->defect_number, ngenerations, i);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
        }
        if (++nbest == debug_patch_n) {
          goto debug_use_this_patch;
        }
      }
    }
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  /*compute statistics*/
  for (fitness_mean = fitness_sigma = 0.0, i = 0; i < max_patches; i++) {
    dp = &dps[i];
    fitness_mean += dp->fitness;
    fitness_sigma += dp->fitness * dp->fitness;
  }
  fitness_mean /= (float)max_patches;
  fitness_sigma = (fitness_sigma / max_patches - (fitness_mean * fitness_mean));
  if (fitness_sigma < 0) {
    fitness_sigma = 0;
  }
  else {
    fitness_sigma = sqrt(fitness_sigma);
  }

  if (parms->verbose > VERBOSE_MODE_DEFAULT)
    fprintf(WHICH_OUTPUT,
            "Initial population for defect %d:\nbest fitness at %d: %2.4f "
            "(%2.4f +- %2.4f)\n",
            defect->defect_number,
            best_i,
            best_fitness,
            fitness_mean,
            fitness_sigma);

  nelite = nint(ELITISM_PCT * max_patches); /* # to keep in next generation */
  if (nelite < 1) {
    nelite = 1;
  }
  nselected = nint(SELECTION_PCT * max_patches); /* # to allow to crossover */
  if (nselected < 2) {
    nselected = 2;
  }
  nreplacements = nint(REPLACEMENT_PCT * max_patches); /* # to replace with
                                                          mutated versions of
                                                          best */
  if (nreplacements + nelite > max_patches) {
    nreplacements = max_patches - nelite;
  }
  g = 0;
  dps = dps1;

  last_fitness = best_fitness;

  ROMP_SCOPE_end

  ROMP_SCOPE_begin

  while (nunchanged < max_unchanged) {
    if (ngenerations == parms->niters) {
      break;
    }

    ROMP_SCOPE_begin
    
    if (dps == dps1) {
      dps_next_generation = dps2;
    }
    else {
      dps_next_generation = dps1;
    }

    last_best = best_fitness;
    for (i = 0; i < max_patches; i++) {
      dp = &dps[i];
      dp->rank = rank = defectPatchRank(dps, i, max_patches);
      ranks[rank] = i;
    }

    /* first add the 'elite' group that are retained unchanged */
    next_gen_index = 0;
    for (i = 0; i < nelite; i++) mrisCopyDefectPatch(&dps[ranks[i]], &dps_next_generation[next_gen_index++]);

    ROMP_SCOPE_end
    ROMP_SCOPE_begin
    
    /* now replace the worst ones with mutated copies of the best */
    for (i = 0; i < nreplacements; i++) {
      ntotalmutations++;

      dp = &dps_next_generation[next_gen_index++];
      mrisCopyDefectPatch(&dps[ranks[i]], dp);
      mrisMutateDefectPatch(dp, &etable, MUTATION_PCT);
      fitness = mrisDefectPatchFitness(&computeDefectContext,
                                       mris,
                                       mris_corrected,
                                       mri,
                                       dp,
                                       vertex_trans,
                                       dvs,
                                       &rp,
                                       h_k1,
                                       h_k2,
                                       mri_k1_k2,
                                       h_white,
                                       h_gray,
                                       h_border,
                                       h_grad,
                                       mri_gray_white,
                                       h_dot,
                                       parms);
#if SAVE_FIT_VALS
      fitness_values[number_of_patches] = fitness;
      if (number_of_patches)
        best_values[number_of_patches] = MAX(best_values[number_of_patches - 1], fitness);
      else {
        best_values[number_of_patches] = fitness;
      }
#endif
      number_of_patches++;

      if (fitness > best_fitness) {
        nmutations++;
        nunchanged = 0;
        best_fitness = fitness;
        best_i = next_gen_index - 1;

        nfinalvertices = nremovedvertices;
        nbestpatch = number_of_patches;

        rp.best_fitness = best_fitness;
        /* save ordering*/
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        /* save current status of vertices */
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

        if (parms->verbose > VERBOSE_MODE_DEFAULT)
          fprintf(WHICH_OUTPUT,
                  "replacement %d MUTATION: new optimal "
                  "fitness found at %d: %2.4e\n",
                  i,
                  best_i,
                  fitness);
        if (parms->verbose == VERBOSE_MODE_LOW) {
          printDefectStatistics(dp);
        }
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname,
                  "%s/rh.defect_%d_surf_%d_%d",
                  parms->save_fname,
                  defect->defect_number,
                  ngenerations - 1,
                  ranks[i]);
          savePatch(mri, mris, mris_corrected, dvs, &dps[ranks[i]], fname, parms);
          sprintf(
              fname, "%s/rh.defect_%d_best_%d_%dm", parms->save_fname, defect->defect_number, ngenerations, ranks[i]);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          if (parms->movie) {
            sprintf(fname, "%s/rh.defect_%d_movie_%d", parms->save_fname, defect->defect_number, nmovies++);
            savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          }
        }
        nmut++;
        if (++nbest == debug_patch_n) {
          dps = dps_next_generation;
          goto debug_use_this_patch;
        }
      }
    }

    ROMP_SCOPE_end
    ROMP_SCOPE_begin

    for (fitness_mean = fitness_sigma = 0.0, i = 0; i < max_patches; i++) {
      dp = &dps[i];
      fitness_mean += dp->fitness;
      fitness_sigma += dp->fitness * dp->fitness;
    }

    fitness_mean /= (float)max_patches;
    fitness_sigma = (fitness_sigma / max_patches - (fitness_mean * fitness_mean));
    if (fitness_sigma < 0) {
      fitness_sigma = 0;
    }
    else {
      fitness_sigma = sqrt(fitness_sigma);
    }
    if (!std::isfinite(fitness_sigma)) {
      DiagBreak();
    }

    two_sigma_sq = (dps[ranks[0]].fitness - dps[ranks[nselected - 1]].fitness);
    if (FZERO(two_sigma_sq)) {
      two_sigma_sq = 1;
    }
    for (fitness_norm = 0.0, j = 0; j < nselected; j++) {
      i = ranks[j];
      dp = &dps[i];
      fitness_norm += exp((dp->fitness - dps[ranks[0]].fitness) / two_sigma_sq);
      /* make
         them
         positive
         and
         increasing */
    }
    if (FZERO(fitness_norm)) /* something wrong */
    {
      for (i = 0; i < max_patches; i++) {
        dp = &dps[i];
        dp->rank = rank = defectPatchRank(dps, i, max_patches);
        if (dp->fitness >= best_fitness) {
          best_fitness = dp->fitness;
          best_i = i;
        }
        ranks[rank] = i;
      }
      break;
    }

    ROMP_SCOPE_end
    ROMP_SCOPE_begin

    /* selection of chromosomes for cross-over */
    ncrossovers = max_patches - (nelite + nreplacements);
    for (l = k = j = 0; j < nselected; j++) {
      int nadd;

      i = ranks[j];
      dp = &dps[i];
      pfitness = exp((dp->fitness - dps[ranks[0]].fitness) / two_sigma_sq) / fitness_norm;
      nadd = nint(pfitness * ncrossovers);
      if (nadd >= ncrossovers) {
        nadd = ncrossovers - 1;
      }
      else if (nadd == 0) {
        nadd = 1;
      }

      for (k = 0; l < ncrossovers && k < nadd; l++, k++) {
        selected[l] = i;
      }
    }
    for (; l < ncrossovers; l++) /* fill out rest of list */
    {
      double p;
      p = randomNumber(0.0, 1.0);
      for (fitness = 0.0, j = 0; j < nselected; j++) {
        i = ranks[j];
        dp = &dps[i];
        pfitness = exp(dp->fitness / two_sigma_sq) / fitness_norm;
        fitness += pfitness;
        if (fitness > p) {
          break;
        }
      }
      selected[l] = i;
    }

    ROMP_SCOPE_end
    ROMP_SCOPE_begin

    for (i = 0; i < ncrossovers; i++) {
      int p1, p2;
      ntotalcross_overs++;

      p1 = selected[i];
      do /* select second parent at random */
      {
        p2 = selected[(int)randomNumber(0, ncrossovers - .001)];
      } while (p2 == p1);

      ROMP_SCOPE_begin

      dp = &dps_next_generation[next_gen_index++];
      mrisCrossoverDefectPatches(&dps[p1], &dps[p2], dp, &etable);
      fitness = mrisDefectPatchFitness(&computeDefectContext,
                                       mris,
                                       mris_corrected,
                                       mri,
                                       dp,
                                       vertex_trans,
                                       dvs,
                                       &rp,
                                       h_k1,
                                       h_k2,
                                       mri_k1_k2,
                                       h_white,
                                       h_gray,
                                       h_border,
                                       h_grad,
                                       mri_gray_white,
                                       h_dot,
                                       parms);
#if SAVE_FIT_VALS
      fitness_values[number_of_patches] = fitness;
      if (number_of_patches)
        best_values[number_of_patches] = MAX(best_values[number_of_patches - 1], fitness);
      else {
        best_values[number_of_patches] = fitness;
      }
#endif
      number_of_patches++;

      ROMP_SCOPE_end

      if (fitness > best_fitness) {

        ROMP_SCOPE_begin

        ncross_overs++;
        nunchanged = 0;
        best_fitness = fitness;
        best_i = next_gen_index - 1;

        nfinalvertices = nremovedvertices;
        nbestpatch = number_of_patches;

        rp.best_fitness = best_fitness;
        /* save ordering*/
        memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
        /* save current status of vertices */
        memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

        if (parms->verbose > VERBOSE_MODE_DEFAULT)
          fprintf(WHICH_OUTPUT,
                  "CROSSOVER (%d x %d): new optimal fitness "
                  "found at %d: %2.4e\n",
                  dps[p1].rank,
                  dps[p2].rank,
                  best_i,
                  fitness);
        if (parms->verbose == VERBOSE_MODE_LOW) {
          printDefectStatistics(dp);
        }
        if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
          sprintf(fname,
                  "%s/rh.defect_%d_surf_%d_%d",
                  parms->save_fname,
                  defect->defect_number,
                  ngenerations - 1,
                  dps[p1].rank);
          savePatch(mri, mris, mris_corrected, dvs, &dps[p1], fname, parms);
          sprintf(fname,
                  "%s/rh.defect_%d_surf_%d_%d",
                  parms->save_fname,
                  defect->defect_number,
                  ngenerations - 1,
                  dps[p2].rank);
          savePatch(mri, mris, mris_corrected, dvs, &dps[p2], fname, parms);
          sprintf(fname,
                  "%s/rh.defect_%d_best_%d_%dc%d_%d",
                  parms->save_fname,
                  defect->defect_number,
                  ngenerations,
                  best_i,
                  dps[p1].rank,
                  dps[p2].rank);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
          savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          if (parms->movie) {
            sprintf(fname, "%s/rh.defect_%d_movie_%d", parms->save_fname, defect->defect_number, nmovies++);
            savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
          }
        }

        ROMP_SCOPE_end

        ncross++;
        if (++nbest == debug_patch_n) {
          dps = dps_next_generation;
          goto debug_use_this_patch;
        }
      }
      else /* mutate it also */
      {
        ROMP_SCOPE_begin

        mrisMutateDefectPatch(dp, &etable, MUTATION_PCT);
        fitness = mrisDefectPatchFitness(&computeDefectContext,
                                         mris,
                                         mris_corrected,
                                         mri,
                                         dp,
                                         vertex_trans,
                                         dvs,
                                         &rp,
                                         h_k1,
                                         h_k2,
                                         mri_k1_k2,
                                         h_white,
                                         h_gray,
                                         h_border,
                                         h_grad,
                                         mri_gray_white,
                                         h_dot,
                                         parms);
#if SAVE_FIT_VALS
        fitness_values[number_of_patches] = fitness;
        if (number_of_patches)
          best_values[number_of_patches] = MAX(best_values[number_of_patches - 1], fitness);
        else {
          best_values[number_of_patches] = fitness;
        }
#endif
        number_of_patches++;
        ntotalmutations++;

        if (fitness > best_fitness) {
          nmutations++;
          nunchanged = 0;
          best_fitness = fitness;
          best_i = next_gen_index - 1;

          nfinalvertices = nremovedvertices;
          nbestpatch = number_of_patches;

          rp.best_fitness = best_fitness;
          /* save ordering*/
          memmove(rp.best_ordering, dp->ordering, nedges * sizeof(int));
          /* save current status of vertices */
          memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

          if (parms->verbose > VERBOSE_MODE_DEFAULT)
            fprintf(WHICH_OUTPUT,
                    "CROSSOVER (%d x %d) & MUTATION: "
                    "new optimal fitness found at %d: %2.4e\n",
                    dps[p1].rank,
                    dps[p2].rank,
                    best_i,
                    fitness);
          if (parms->verbose == VERBOSE_MODE_LOW) {
            printDefectStatistics(dp);
          }
          if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
            sprintf(fname,
                    "%s/rh.defect_%d_surf_%d_%d",
                    parms->save_fname,
                    defect->defect_number,
                    ngenerations - 1,
                    dps[p1].rank);
            savePatch(mri, mris, mris_corrected, dvs, &dps[p1], fname, parms);
            sprintf(fname,
                    "%s/rh.defect_%d_surf_%d_%d",
                    parms->save_fname,
                    defect->defect_number,
                    ngenerations - 1,
                    dps[p2].rank);
            savePatch(mri, mris, mris_corrected, dvs, &dps[p2], fname, parms);
            sprintf(fname,
                    "%s/rh.defect_%d_best_%d_%dcm%d_%d",
                    parms->save_fname,
                    defect->defect_number,
                    ngenerations,
                    best_i,
                    dps[p1].rank,
                    dps[p2].rank);
            savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);

            sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
            savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
            if (parms->movie) {
              sprintf(fname, "%s/rh.defect_%d_movie_%d", parms->save_fname, defect->defect_number, nmovies++);
              savePatch(mri, mris, mris_corrected, dvs, dp, fname, parms);
            }
          }

          if (++nbest == debug_patch_n) {
            dps = dps_next_generation;
            goto debug_use_this_patch;
          }
          nmut++;
          ncross++;
        }

        ROMP_SCOPE_end
      }
    }

    ROMP_SCOPE_end
    ROMP_SCOPE_begin

    /* make next generation current */
    if (dps == dps1) {
      dps = dps2;
      dps_next_generation = dps1;
    }
    else {
      dps = dps1;
      dps_next_generation = dps2;
    }

    best_fitness = dps[0].fitness;
    best_i = 0;
    for (fitness_mean = fitness_sigma = 0.0, i = 0; i < max_patches; i++) {
      dp = &dps[i];
      dp->rank = rank = defectPatchRank(dps, i, max_patches);
      if (dp->fitness >= best_fitness) {
        best_fitness = dp->fitness;
        best_i = i;
      }
      ranks[rank] = i;
      fitness_mean += dp->fitness;
      fitness_sigma += dp->fitness * dp->fitness;
    }

    fitness_mean /= (float)max_patches;
    fitness_sigma = sqrt(fitness_sigma / max_patches - (fitness_mean * fitness_mean));
    if (parms->verbose > VERBOSE_MODE_DEFAULT)
      fprintf(WHICH_OUTPUT,
              "generation %d complete, optimal fitness = "
              "%2.4e (%2.4e +- %2.4e)\n",
              ++g,
              best_fitness,
              fitness_mean,
              fitness_sigma);
    if (FEQUAL(last_best, best_fitness)) {
      nunchanged++;
    }
    else {
      last_euthanasia = -1;
      nunchanged = 0;
    }

    if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
      sprintf(fname, "%s/rh.defect_%d_generation_%d", parms->save_fname, defect->defect_number, ngenerations);
      savePatch(mri, mris, mris_corrected, dvs, &dps[best_i], fname, parms);
    }

#define NEXT 5

    if (parms->vertex_eliminate) {
      static int count = 1;
      int ndeleted;
      if (nunchanged >= max_unchanged) {
        // will eventually break out
        if (last_euthanasia < 0) {
          last_euthanasia = ngenerations;
          last_fitness = best_fitness;
        } 
	else if (last_euthanasia + NEXT <= ngenerations) {
          break;
        }
        if (last_fitness < best_fitness) {
          last_euthanasia = ngenerations;
        }
        if (nunchanged >= max_unchanged) {
          nunchanged -= NEXT;
        }

        if (parms->verbose == VERBOSE_MODE_LOW) {
          fprintf(WHICH_OUTPUT, "Deleting worst vertices : ");
        }
        ndeleted = deleteWorstVertices(mris_corrected, &rp, defect, vertex_trans, 0.2, count);
        nremovedvertices += ndeleted;
        if (parms->verbose == VERBOSE_MODE_LOW) {
          fprintf(WHICH_OUTPUT, "%d vertices have been deleted\n", ndeleted);
        }
        if (ndeleted == 0) {
          break;
        }
        count++;
      }
      else if (ngenerations >= 10 && (ngenerations % 3 == 0)) {
        ndeleted = deleteWorstVertices(mris_corrected, &rp, defect, vertex_trans, 0.1, count);
        nremovedvertices += ndeleted;
        if (parms->verbose == VERBOSE_MODE_LOW) {
          if (ndeleted == 1) {
            fprintf(WHICH_OUTPUT, "Worst vertex has been deleted\n");
          }
          else if (ndeleted)
            fprintf(WHICH_OUTPUT, "Worst %d vertices have been deleted (corresponds to ~1%%)\n", ndeleted);
        }
      }
    }

    ngenerations++;
    ROMP_SCOPE_end
  }

  ROMP_SCOPE_end

debug_use_this_patch:

  ROMP_SCOPE_begin

  dp = &dps[best_i];

  if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
    /* save eliminated vertices */
    for (i = 0; i < mris->nvertices; i++) {
      mris->vertices[i].curv = 0;
    }

    for (i = 0; i < defect->nvertices; i++) {
      if (vertex_trans[dp->defect->vertices[i]] < 0) {
        continue;
      }
      k = mris_corrected->vertices[vertex_trans[dp->defect->vertices[i]]].ripflag;
      if (k) {
        mris->vertices[defect->vertices[i]].curv = k;
      }
    }
    sprintf(fname, "%s/rh.eliminated_%d", parms->save_fname, defect->defect_number);
    MRISwriteCurvature(mris, fname);

    /* save vertex statistics */
    for (i = 0; i < defect->nvertices; i++) {
      mris->vertices[defect->vertices[i]].curv = 0;
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      mris->vertices[defect->vertices[i]].curv = rp.vertex_fitness[i];
    }
    sprintf(fname, "%s/rh.rp_fitness_%d", parms->save_fname, defect->defect_number);
    MRISwriteCurvature(mris, fname);
  }

  nkilled += nfinalvertices;

  /* use the best ordering to retessellate the defected patch */
  memmove(dp->ordering, rp.best_ordering, nedges * sizeof(int));
  memmove(defect->status, rp.status, defect->nvertices * sizeof(char));

  /* set back to unrip the correct vertices*/
  for (i = 0; i < defect->nvertices; i++) {
    /* careful with vertex_trans */
    int vni;
    vni = vertex_trans[dp->defect->vertices[i]];
    if (vni < 0) {
      continue;
    }
    if (defect->status[i] != DISCARD_VERTEX) {
      mris_corrected->vertices[vni].ripflag = 0;
    }
    else {
      mris_corrected->vertices[vni].ripflag = 1;
    }
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  fitness = mrisDefectPatchFitness(&computeDefectContext,
                                   mris,
                                   mris_corrected,
                                   mri,
                                   dp,
                                   vertex_trans,
                                   dvs,
                                   &rp,
                                   h_k1,
                                   h_k2,
                                   mri_k1_k2,
                                   h_white,
                                   h_gray,
                                   h_border,
                                   h_grad,
                                   mri_gray_white,
                                   h_dot,
                                   parms);

  defect->fitness = fitness; /* saving the fitness of the patch */

  if (fitness != best_fitness)
    fprintf(WHICH_OUTPUT, "Warning - incorrect dp selected!!!!(%f >= %f ) \n", fitness, best_fitness);

  if (parms->verbose == VERBOSE_MODE_LOW) {
    printDefectStatistics(dp);
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin
  
  /* compute the final tessellation */
  retessellateDefect(mris, mris_corrected, dvs, dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  /* detect the new set of faces */
  detectDefectFaces(mris_corrected, dp);

  /* orient the patch faces */
  orientDefectFaces(mris_corrected, dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* smooth original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, dp, parms->smooth, parms->match);

  /* reset the number of original faces in the
     surface before the retessellation */
  // mrisRestoreFaceVertexState(mris_corrected, dvs) ;

  nintersections = 0; /* compiler warnings... */
  if (parms->check_surface_intersection)
    nintersections = mrisCountIntersectingFaces(mris_corrected, dp->tp.faces, dp->tp.nfaces);

  /* saving self-intersection */
  defect->intersect = nintersections;

  if (parms->verbose == VERBOSE_MODE_HIGH) {
    fprintf(WHICH_OUTPUT,
            "PATCH #:%03d:  FITNESS:   %2.2f\n              "
            "MUTATIONS: %d (out of %d)\n              "
            "CROSSOVERS: %d (out of %d)\n",
            defect->defect_number,
            fitness,
            nmutations,
            ntotalmutations,
            ncross_overs,
            ntotalcross_overs);
    fprintf(WHICH_OUTPUT, "              ELIMINATED VERTICES:  %d (out of %d)\n", nfinalvertices, defect->nvertices);
    fprintf(
        WHICH_OUTPUT, "              BEST PATCH #: %d (out of %d generated patches)\n", nbestpatch, number_of_patches);
    if (parms->check_surface_intersection)
      fprintf(
          WHICH_OUTPUT, "              NUMBER OF INTERSECTING FACES: %d (out of %d) \n", nintersections, dp->tp.nfaces);
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  /* should free the tessellated patch structure */
  TPfree(&dp->tp);

  /* discard the vertices that are not used in the
     final tessellation and set marks to zero!!!*/
  for (i = 0; i < dp->nedges; i++)
    if (dp->etable->edges[i].used == USED_IN_NEW_TESSELLATION ||
        dp->etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
      mris_corrected->vertices[dp->etable->edges[i].vno1].marked = FINAL_VERTEX;
      mris_corrected->vertices[dp->etable->edges[i].vno2].marked = FINAL_VERTEX;
    }

  for (k = i = 0; i < dp->defect->nvertices; i++) {
    if (mris_corrected->vertices[vertex_trans[dp->defect->vertices[i]]].marked != FINAL_VERTEX) {
      if (dp->defect->status[i] != DISCARD_VERTEX) {
        k++;
      }
      dp->defect->status[i] = DISCARD_VERTEX;
      mris_corrected->vertices[vertex_trans[dp->defect->vertices[i]]].ripflag = 1;
    }
    mris_corrected->vertices[vertex_trans[dp->defect->vertices[i]]].marked = 0;
  }

  for (i = 0; i < dp->defect->nborder; i++) {
    mris_corrected->vertices[vertex_trans[dp->defect->border[i]]].marked = 0;
  }

  ROMP_SCOPE_end
  ROMP_SCOPE_begin

  /* free everything */
  destructComputeDefectContext(&computeDefectContext);
  mrisFreeDefectVertexState(dvs);

  for (i = 0; i < max_patches; i++) {
    free(dps1[i].ordering);
    free(dps2[i].ordering);
  }

  if (etable.use_overlap) {
    for (i = 0; i < nedges; i++) {
      if (etable.overlapping_edges[i]) {
        free(etable.overlapping_edges[i]);
      }
    }
    free(etable.overlapping_edges);
    free(etable.noverlap);
    free(etable.flags);
  }
  free(etable.edges);

  free(rp.best_ordering);
  free(rp.status);
  free(rp.nused);
  free(rp.vertex_fitness);

  if (mri_defect) {
    MRIfree(&mri_defect);
  }
  if (mri_defect_white) {
    MRIfree(&mri_defect_white);
  }
  if (mri_defect_gray) {
    MRIfree(&mri_defect_gray);
  }
  if (mri_defect_sign) {
    MRIfree(&mri_defect_sign);
  }

#if SAVE_FIT_VALS
  {
    FILE *f;
    int n;
    f = fopen("./optimal1.plt", "w+");
    for (n = 0; n < number_of_patches; n++) {
      fprintf(f, "%d %2.2f\n", n, fitness_values[n]);
    }
    fclose(f);
    f = fopen("./optimal2.plt", "w+");
    for (n = 0; n < number_of_patches; n++) {
      fprintf(f, "%d %2.2f\n", n, best_values[n]);
    }
    fclose(f);
  }
#endif

  ROMP_SCOPE_end

  return (NO_ERROR);
}

static int mrisComputeRandomRetessellation(MRI_SURFACE *mris,
                                           MRI_SURFACE *mris_corrected,
                                           MRI *mri,
                                           DEFECT *defect,
                                           int *vertex_trans,
                                           EDGE *et,
                                           int nedges,
                                           ES *es,
                                           int nes,
                                           HISTOGRAM *h_k1,
                                           HISTOGRAM *h_k2,
                                           MRI *mri_k1_k2,
                                           HISTOGRAM *h_white,
                                           HISTOGRAM *h_gray,
                                           HISTOGRAM *h_border,
                                           HISTOGRAM *h_grad,
                                           MRI *mri_gray_white,
                                           HISTOGRAM *h_dot,
                                           TOPOLOGY_PARMS *parms)
{
  DEFECT_VERTEX_STATE *dvs;
  DEFECT_PATCH dp;
  int niters, m, tmp, best_i, i, j, k, noverlap;
  int *overlap, ngenerations, nbests, last_euthanasia;
  int nremovedvertices, nfinalvertices;
  double fitness, best_fitness;
  static int dno = 0; /* for debugging */
  // static int nmovies=1;
  /* for making movies :
     0 is left for the original surface*/
  EDGE_TABLE etable;
  MRI *mri_defect, *mri_defect_white, *mri_defect_gray, *mri_defect_sign;
  char fname[500];
  RP rp;
  int nbestpatch;
  static int first_time = 1;

  if (first_time) {
    char *cp;

    if ((cp = getenv("FS_QCURV")) != NULL) {
      parms->l_qcurv = atof(cp);
      printf("setting qcurv = %2.3f\n", l_qcurv);
    }
    if ((cp = getenv("FS_CURV")) != NULL) {
      parms->l_curv = atof(cp);
      printf("setting curv = %2.3f\n", l_curv);
    }
    if ((cp = getenv("FS_MRI")) != NULL) {
      parms->l_mri = atof(cp);
      printf("setting mri = %2.3f\n", l_mri);
    }
    if ((cp = getenv("FS_UNMRI")) != NULL) {
      parms->l_unmri = atof(cp);
      printf("setting unmri = %2.3f\n", l_unmri);
    }
    first_time = 0;
  }

  dno++; /* for debugging */

  if (nedges > 200000) {
    // add some code here to select a good tessellation (ordering) FLO
    // mrisRetessellateDefect(mris, mris_corrected, defect,
    // vertex_trans, et, nedges, NULL, NULL) ;
    tessellatePatch(mri, mris, mris_corrected, defect, vertex_trans, et, nedges, NULL, NULL, parms);

    return (NO_ERROR);
  };

  etable.use_overlap = 0;  // XXX parms->edge_table;
  etable.nedges = nedges;
  etable.edges = (EDGE *)calloc(nedges, sizeof(EDGE));
  memmove(etable.edges, et, nedges * sizeof(EDGE));

  if (etable.use_overlap) {
    etable.overlapping_edges = (int **)calloc(nedges, sizeof(int *));
    etable.noverlap = (int *)calloc(nedges, sizeof(int));
    etable.flags = (unsigned char *)calloc(nedges, sizeof(unsigned char));
    overlap = (int *)calloc(nedges, sizeof(int));
    if (!etable.edges || !etable.overlapping_edges || !etable.noverlap || !overlap)
      ErrorExit(ERROR_NOMEMORY,
                "mrisComputeOptimalRetessellation: Excessive "
                "topologic defect encountered: could not allocate "
                "%d edge table",
                nedges);

    for (i = 0; i < nedges; i++) /* compute overlapping for each edge */
    {
      if (nedges > 50000 && !(i % 25000)) {
        printf("%d of %d edges processed\n", i, nedges);
      }
      etable.noverlap[i] = 0;
      for (noverlap = j = 0; j < nedges; j++) {
        if (j == i) {
          continue;
        }
        if (edgesIntersect(mris_corrected, &et[i], &et[j])) {
          overlap[noverlap] = j;
          noverlap++;
        }
        if (noverlap > MAX_EDGES) {
          break;
        }
      }
      if (noverlap > 0) {
        if (noverlap > MAX_EDGES) {
          etable.noverlap[i] = MAX_EDGES;
          etable.flags[i] |= ET_OVERLAP_LIST_INCOMPLETE;
        }
        else {
          etable.noverlap[i] = noverlap;
        }

        etable.overlapping_edges[i] = (int *)calloc(etable.noverlap[i], sizeof(int));
        if (!etable.overlapping_edges[i])
          ErrorExit(ERROR_NOMEMORY,
                    "mrisComputeOptimalRetessellation: Excessive "
                    "topologic defect encountered: could not allocate "
                    "overlap list %d "
                    "with %d elts",
                    i,
                    etable.noverlap[i]);
        memmove(etable.overlapping_edges[i], overlap, etable.noverlap[i] * sizeof(int));
      };
    }
    free(overlap);
  }

  /* allocate the volume constituted by the potential edges */
  mri_defect = mri_defect_white = mri_defect_gray = mri_defect_sign = NULL;
  if (!FZERO(parms->l_unmri)) {
    mri_defect = mriDefectVolume(mris_corrected, &etable, parms);
    mri_defect_white = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_gray = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    mri_defect_sign = MRIalloc(mri_defect->width, mri_defect->height, mri_defect->depth, MRI_FLOAT);
    defectVolumeLikelihood(mri,
                           mri_defect,
                           mri_defect_white,
                           mri_defect_gray,
                           h_white,
                           h_gray,
                           defect->white_mean,
                           defect->gray_mean,
                           0,
                           0);
  };

  if ((!FZERO(parms->l_unmri)) && parms->save_fname &&
      (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
    sprintf(fname, "%s/white_%d.mgh", parms->save_fname, defect->defect_number);
    MRIwrite(mri_defect_white, fname);
    sprintf(fname, "%s/gray_%d.mgh", parms->save_fname, defect->defect_number);
    MRIwrite(mri_defect_gray, fname);
  }

  dvs = mrisRecordVertexState(mris_corrected, defect, vertex_trans);

  ngenerations = 0;
  last_euthanasia = -1;
  nremovedvertices = 0;
  nfinalvertices = 0;
  nbests = 0;
  nbestpatch = 0;
  best_fitness = 0;

  /* generate Random Patch */
  rp.best_ordering = (int *)malloc(nedges * sizeof(int));
  rp.status = (char *)malloc(defect->nvertices * sizeof(char));
  memmove(rp.status, defect->status, defect->nvertices * sizeof(char));
  rp.nused = (int *)calloc(defect->nvertices, sizeof(int));
  rp.vertex_fitness = (float *)calloc(defect->nvertices, sizeof(float));

  if (parms->retessellation_mode) {
    dp.retessellation_mode = USE_SOME_VERTICES;
  }
  else {
    dp.retessellation_mode = USE_ALL_VERTICES;
  }

  dp.nedges = nedges;
  dp.defect = defect;
  dp.etable = &etable;
  dp.ordering = (int *)calloc(nedges, sizeof(int));
  if (!dp.ordering) ErrorExit(ERROR_NOMEMORY, "could not allocate defect patch with %d indices", nedges);
  for (j = 0; j < nedges; j++) {
    dp.ordering[j] = j; /* initial in same order - will change later */
  }

  dp.mri_defect = mri_defect;
  dp.mri_defect_white = mri_defect_white;
  dp.mri_defect_gray = mri_defect_gray;
  dp.mri_defect_sign = mri_defect_sign;
  dp.mri = mri;

  ComputeDefectContext computeDefectContext;
    constructComputeDefectContext(&computeDefectContext);

  niters = 0;
  while (niters <= parms->niters) {
    if (niters % 100 == 0) {
      fprintf(stderr, "iteration %d\n", niters);
    }

    if (niters) {
      /* first one is in same order as
         original edge table -random ordering otherwise */
      for (m = 0; m < 11; m++)
        for (j = 0; j < nedges; j++) {
          k = nint(randomNumber(0.0, (double)nedges - 1));

          tmp = dp.ordering[j];
          dp.ordering[j] = dp.ordering[k];
          dp.ordering[k] = tmp;
        }
    }

    fitness = mrisDefectPatchFitness(&computeDefectContext,
                                     mris,
                                     mris_corrected,
                                     mri,
                                     &dp,
                                     vertex_trans,
                                     dvs,
                                     &rp,
                                     h_k1,
                                     h_k2,
                                     mri_k1_k2,
                                     h_white,
                                     h_gray,
                                     h_border,
                                     h_grad,
                                     mri_gray_white,
                                     h_dot,
                                     parms);

#if SAVE_FIT_VALS
    fitness_values[niters] = fitness;
    if (niters) {
      best_values[niters] = MAX(best_values[niters - 1], fitness);
    }
    else {
      best_values[niters] = fitness;
    }
#endif

    if ((!niters) || (fitness > best_fitness)) {
      best_fitness = fitness;
      best_i = niters;
      printf("new optimal fitness found at %d: %2.4f\n", niters, fitness);

      nfinalvertices = nremovedvertices;
      nbestpatch = niters;

      rp.best_fitness = best_fitness;
      /* save ordering*/
      memmove(rp.best_ordering, dp.ordering, nedges * sizeof(int));
      /* save current status of vertices */
      memmove(rp.status, defect->status, defect->nvertices * sizeof(char));

      if (parms->verbose == VERBOSE_MODE_LOW) {
        printDefectStatistics(&dp);
      }
      if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
        sprintf(fname, "%s/rh.defect_%d_best_%d_%d", parms->save_fname, defect->defect_number, ngenerations, niters);
        savePatch(mri, mris, mris_corrected, dvs, &dp, fname, parms);
        sprintf(fname, "%s/rh.defect_%d_best_%d", parms->save_fname, defect->defect_number, nbests++);
        savePatch(mri, mris, mris_corrected, dvs, &dp, fname, parms);
      }
    }
    niters++;
  }

  if (parms->save_fname && (parms->defect_number < 0 || (parms->defect_number == defect->defect_number))) {
    /* save eliminated vertices */
    for (i = 0; i < mris->nvertices; i++) {
      mris->vertices[i].curv = 0;
    }

    for (i = 0; i < defect->nvertices; i++) {
      if (vertex_trans[dp.defect->vertices[i]] < 0) {
        continue;
      }
      k = mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].ripflag;
      if (k) {
        mris->vertices[defect->vertices[i]].curv = k;
      }
    }
    sprintf(fname, "%s/rh.eliminated_%d", parms->save_fname, defect->defect_number);
    MRISwriteCurvature(mris, fname);

    /* save vertex statistics */
    for (i = 0; i < defect->nvertices; i++) {
      mris->vertices[defect->vertices[i]].curv = 0;
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      mris->vertices[defect->vertices[i]].curv = rp.vertex_fitness[i];
    }
    sprintf(fname, "%s/rh.rp_fitness_%d", parms->save_fname, defect->defect_number);
    MRISwriteCurvature(mris, fname);
  }

  fprintf(WHICH_OUTPUT,
          "%d (out of %d) vertices were eliminated from retessellated patch\n",
          nfinalvertices,
          defect->nvertices);
  fprintf(WHICH_OUTPUT, "best patch # is %d - total number of generated patches %d\n", nbestpatch, niters - 1);
  nkilled += nfinalvertices;

  /* use the best ordering to retessellate the defected patch */
  memmove(dp.ordering, rp.best_ordering, nedges * sizeof(int));
  memmove(defect->status, rp.status, defect->nvertices * sizeof(char));

  /* set back to unrip the correct vertices*/
  for (i = 0; i < defect->nvertices; i++) {
    int vni;
    vni = vertex_trans[dp.defect->vertices[i]];
    if (vni < 0) {
      continue;
    }
    if (defect->status[i] != DISCARD_VERTEX) {
      mris_corrected->vertices[vni].ripflag = 0;
    }
    else {
      mris_corrected->vertices[vni].ripflag = 1;
    }
  }

  fitness = mrisDefectPatchFitness(&computeDefectContext,
                                   mris,
                                   mris_corrected,
                                   mri,
                                   &dp,
                                   vertex_trans,
                                   dvs,
                                   &rp,
                                   h_k1,
                                   h_k2,
                                   mri_k1_k2,
                                   h_white,
                                   h_gray,
                                   h_border,
                                   h_grad,
                                   mri_gray_white,
                                   h_dot,
                                   parms);

  if (fitness != best_fitness)
    fprintf(WHICH_OUTPUT, "Warning - incorrect dp selected!!!!(%f >= %f ) \n", fitness, best_fitness);

  if (parms->verbose == VERBOSE_MODE_LOW) {
    printDefectStatistics(&dp);
  }

  /* compute the final tessellation */
  retessellateDefect(mris, mris_corrected, dvs, &dp);
  mrisCheckVertexFaceTopology(mris_corrected);


  /* detect the new set of faces */
  detectDefectFaces(mris_corrected, &dp);
  /* orient the patch faces */
  orientDefectFaces(mris_corrected, &dp);
  mrisCheckVertexFaceTopology(mris_corrected);

  /* smooth original vertices in the retessellated patch */
  defectMatch(mri, mris_corrected, &dp, parms->smooth, parms->match);

  /* reset the number of original faces in the
     surface before the retessellation */
  // mrisRestoreFaceVertexState(mris_corrected, dvs) ;

  /* should free the tessellated patch structure */
  TPfree(&dp.tp);

  /* discard the vertices that are not used in the
     final tessellation and set marks to zero!!!*/
  for (i = 0; i < dp.nedges; i++)
    if (dp.etable->edges[i].used == USED_IN_NEW_TESSELLATION || dp.etable->edges[i].used == USED_IN_BOTH_TESSELLATION) {
      mris_corrected->vertices[dp.etable->edges[i].vno1].marked = FINAL_VERTEX;
      mris_corrected->vertices[dp.etable->edges[i].vno2].marked = FINAL_VERTEX;
    }

  for (k = i = 0; i < dp.defect->nvertices; i++) {
    if (mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].marked != FINAL_VERTEX) {
      if (dp.defect->status[i] != DISCARD_VERTEX) {
        k++;
      }
      dp.defect->status[i] = DISCARD_VERTEX;
      mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].ripflag = 1;
    }
    mris_corrected->vertices[vertex_trans[dp.defect->vertices[i]]].marked = 0;
  }

  for (i = 0; i < dp.defect->nborder; i++) {
    mris_corrected->vertices[vertex_trans[dp.defect->border[i]]].marked = 0;
  }

  /* free everything */
  destructComputeDefectContext(&computeDefectContext);
  mrisFreeDefectVertexState(dvs);

  free(dp.ordering);

  if (etable.use_overlap) {
    for (i = 0; i < nedges; i++) {
      if (etable.overlapping_edges[i]) {
        free(etable.overlapping_edges[i]);
      }
    }
    free(etable.overlapping_edges);
    free(etable.noverlap);
    free(etable.flags);
  }
  free(etable.edges);

  free(rp.best_ordering);
  free(rp.status);
  free(rp.nused);
  free(rp.vertex_fitness);

  if (mri_defect) {
    MRIfree(&mri_defect);
  }
  if (mri_defect_white) {
    MRIfree(&mri_defect_white);
  }
  if (mri_defect_gray) {
    MRIfree(&mri_defect_gray);
  }
  if (mri_defect_sign) {
    MRIfree(&mri_defect_sign);
  }

#if SAVE_FIT_VALS
  {
    FILE *f;
    int n;
    f = fopen("./random1.plt", "w+");
    for (n = 0; n < niters; n++) {
      fprintf(f, "%d %2.2f\n", n, fitness_values[n]);
    }
    fclose(f);
    f = fopen("./random2.plt", "w+");
    for (n = 0; n < niters; n++) {
      fprintf(f, "%d %2.2f\n", n, best_values[n]);
    }
    fclose(f);
  }
#endif

  return (NO_ERROR);
}

static int mrisRetessellateDefect(MRI_SURFACE *mris,
                                  MRI_SURFACE *mris_corrected,
                                  DEFECT *defect,
                                  int *vertex_trans,
                                  EDGE *et,
                                  int nedges,
                                  int *ordering,
                                  EDGE_TABLE *etable)
{
  double max_len;
  int i, j, max_i, max_added, nadded, index;
  int (*intersection_function)(MRI_SURFACE * mris, DEFECT * defect, EDGE * e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2);

  max_len = 0;
  max_i = 0;
  max_added = nadded = 0;
  intersection_function = intersectDefectEdges;
  for (index = 0; index < nedges; index++) {
    if (ordering) {
      i = ordering[index];
    }
    else {
      i = index;
    }

    if (et[i].used && et[i].used != USED_IN_ORIGINAL_TESSELLATION) /* already exists in
                                                                                    tessellation -
                                                                                    don't add it again */
    {
      continue;
    }

    if (etable && etable->use_overlap) /* use pre-computed
                                                        intersection table */
    {
      int intersects = 0;

      for (j = 0; j < etable->noverlap[i]; j++)
        if (et[etable->overlapping_edges[i][j]].used &&
            et[etable->overlapping_edges[i][j]].used != USED_IN_ORIGINAL_TESSELLATION) {
          intersects = 1;
          break;
        }
      if (intersects) {
        continue;
      }
      if (etable->flags[i] & ET_OVERLAP_LIST_INCOMPLETE) {
        intersection_function = intersectDefectEdges;
      }
      else {
        intersection_function = intersectDefectConvexHullEdges;
      }
    }

    if ((*intersection_function)(mris_corrected, defect, &et[i], NULL, vertex_trans, NULL, NULL) == 0) {
      mrisAddEdge(mris_corrected, et[i].vno1, et[i].vno2);
      if (et[i].used) {
        et[i].used = USED_IN_BOTH_TESSELLATION;
        nadded++;
      }
      else {
        et[i].used = USED_IN_NEW_TESSELLATION;
        nadded++;
      }
      if (et[i].len > max_len) {
        max_len = et[i].len;
        max_added = nadded - 1;
        max_i = i;
      }
    }
  }

  return (NO_ERROR);
}

#if 0
static int
mrisCheckDefectEdges(MRI_SURFACE *mris, DEFECT *defect, int vno,
                     int *vertex_trans)
{
  int    i, n1, n2, n3, vno2, fno1, fno2 ;
  EDGE   edge1, edge2 ;
  VERTEX *v, *vn ;
  FACE   *f ;

  v = &mris->vertices[vno] ;
  edge1.vno1 = vno ;
  for (n1 = 0 ; n1 < v->num ; n1++)
  {
    n2 = v->n[n1] == VERTICES_PER_FACE-1 ? 0 : v->n[n1]+1 ;
    fno1 = v->f[n1] ;
    edge1.vno2 = mris->faces[fno1].v[n2] ;
    for (i = 0 ; i < defect->nborder ; i++)
    {
      vno2 = vertex_trans[defect->border[i]] ;
      if (vno2 < 0)
      {
        continue ;
      }
      vn = &mris->vertices[vno2] ;
      for (n2 = 0 ; n2 < vn->num ; n2++)
      {
        fno2 = vn->f[n2] ;
        f = &mris->faces[fno2] ;
        if (fno2 == Gdiag_no)
        {
          DiagBreak() ;
        }
        if ((fno2 == 103815 && fno1 == 101608) ||
            (fno1 == 103815 && fno2 == 101608))
        {
          DiagBreak() ;
        }
        for (n3 = 0 ; n3 < VERTICES_PER_FACE ; n3++)
        {
          edge2.vno1 = f->v[n3] ;
          edge2.vno2 = f->v[n3 == VERTICES_PER_FACE-1 ? 0 : n3+1] ;
          if (edge2.vno1 == edge1.vno1 || edge2.vno1 == edge1.vno2 ||
              edge2.vno2 == edge1.vno1 || edge2.vno2 == edge1.vno2)
          {
            continue ;
          }
          if (edgesIntersect(mris, &edge1, &edge2))
          {
            fprintf
            (WHICH_OUTPUT,
             "triangles %d (a=%f) and %d (a=%f) intersect!\n",
             fno1, mris->faces[fno1].area, fno2,
             mris->faces[fno2].area) ;
            mrisDumpTriangle(mris, fno1) ;
            mrisDumpTriangle(mris, fno2) ;
          }
        }
      }
    }
  }
  return(NO_ERROR) ;
}
#endif

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the edge e intersects any edges in the defect or it's border.
  ------------------------------------------------------*/

    static bool isHit(MRI_SURFACE * mris, int vno, DEFECT *defect, EDGE* e, int n, const char* whyTracing) {
      bool result = false;

      if (vno == e->vno1 || vno == e->vno2 || vno < 0) {
        if (whyTracing) fprintf(stderr, "isHit false because shared vno when %s\n", whyTracing);
        return false;
      }

      EDGE edge2;
      edge2.vno1 = vno;
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

      if (n >= vt->vnum) {
        if (whyTracing) fprintf(stderr, "isHit vno:%d deleted edge n:%d when %s\n",vno,n, whyTracing);
        return false;
      }
      
      if (whyTracing) fprintf(stderr, "isHit vno:%d v->v[n]:%d e:%p e->vno1:%d e->vno2:%d when %s\n",vno,vt->v[n],e,e->vno1,e->vno2, whyTracing);
      
      if (vt->v[n] == e->vno1 || vt->v[n] == e->vno2) {
        if (whyTracing) fprintf(stderr, "isHit false because shared second vno when %s\n", whyTracing);
        goto Done;
      }
      
      edge2.vno2 = vt->v[n];
      if (defect->optimal_mapping && mris->vertices[vt->v[n]].fixedval == 0) {
        if (whyTracing) fprintf(stderr, "isHit false optimal_mapping\n");
        goto Done;
      }
      if (edgesIntersect(mris, e, &edge2)) {
        result = true;
      }
      if (whyTracing) fprintf(stderr, "isHit edgesIntersect returned %d when %s\n",(int)result, whyTracing);
    Done:
      return result;
    }  

    typedef struct PossiblyIntersectingGreatArcs_callback_context {
      MRI_SURFACE*                  mris;
      DEFECT*                       defect;
      EDGE*                         e;
      IntersectDefectEdgesContext*  ide_ctx;
      
      int                           old_firstHit_vnoLo;     // The same vno pair can be found from several directions
      int                           old_firstHit_vnoHi;     // but any find is ok
      
      bool                          emit_line_py;
      
      int                           tried;                  // count until the first match
                                                            //
      int                           matches;                // total number of matches found, including after the first match
      int                           someMatchingIndexPlus1; // result
    } PossiblyIntersectingGreatArcs_callback_context;

    static bool possiblyIntersectingGreatArcs_callback (void* void_ctx, int indexIntoIDECtxEntries, bool* pIsHit) {
      PossiblyIntersectingGreatArcs_callback_context* ctx = (PossiblyIntersectingGreatArcs_callback_context*)void_ctx;

      if (ctx->someMatchingIndexPlus1 == 0) ctx->tried++;
      
      IntersectDefectEdgesContext_Entry* entry = &ctx->ide_ctx->entries[indexIntoIDECtxEntries];
      *pIsHit = isHit(ctx->mris, entry->vno, ctx->defect, ctx->e, entry->n, NULL);

      if (ctx->emit_line_py) {
        VERTEX_TOPOLOGY const * const v0t = &ctx->mris->vertices_topology[entry->vno]; 
        VERTEX          const * const v0  = &ctx->mris->vertices         [entry->vno]; 
        VERTEX const * const v1 = &ctx->mris->vertices[v0t->v[entry->n]]; 
        fprintf(stderr, " [3, %f, %f, %f, %f, %f, %f, %d], # line.py trial\n", v0->cx,v0->cy,v0->cz, v1->cx,v1->cy,v1->cz, *pIsHit);
      }
      
      int entry_vnoLo = entry->vno;
      int entry_vnoHi = ctx->mris->vertices_topology[entry->vno].v[entry->n];
      if (entry_vnoLo > entry_vnoHi) { int temp = entry_vnoLo; entry_vnoLo = entry_vnoHi; entry_vnoHi = temp; }
      
      if ( ctx->old_firstHit_vnoLo == entry_vnoLo 
        && ctx->old_firstHit_vnoHi == entry_vnoHi
        && !*pIsHit) {
        fprintf(stderr, "callback passed the is_old_firstHit, now it is not hitting\n");
        isHit(ctx->mris, entry->vno, ctx->defect, ctx->e, entry->n, "possiblyIntersectingGreatArcs_callback - showing why not hit");
        *(int*)(-1) = 0;
      }

      if (!*pIsHit)
        return true;    // Keep going when missing, and no need to keep
      
      ctx->matches++;
      
      ctx->someMatchingIndexPlus1 = indexIntoIDECtxEntries + 1;
      
      return false;
    }

static int intersectDefectEdges(MRI_SURFACE *mris, DEFECT *defect, EDGE *e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2)
{
  static long stats_count    = 0;
  static long stats_limit    = 1;
  
  static int once;
  static bool asked_do_old_way,asked_do_new_way,asked_do_stats;
  if (!once++) {
    if (getenv("FREESURFER_intersectDefectEdges_old"))   asked_do_old_way = true;
    if (getenv("FREESURFER_intersectDefectEdges_new"))   asked_do_new_way = true;
    if (getenv("FREESURFER_intersectDefectEdges_stats")) asked_do_stats   = true;
  }
  bool do_old_way = asked_do_old_way;
  bool do_new_way = asked_do_new_way || !asked_do_old_way;
  
  if (!ctx) { do_new_way = false; do_old_way = true; }

  bool showStats = asked_do_stats || (do_old_way && do_new_way);
  if (showStats) {
    stats_count++;
    
    if (stats_count < stats_limit) {
      showStats = false;
    } else {
      if (stats_limit < 20000) stats_limit *= 2; else stats_limit += 20000;
    }
  }

  // This code says 
  //  for each vertex in the defect
  //      for each edge from the vertex 
  //          if the edge shares a vertex with 'e' 
  //          or the edge does not intersect 'e' 
  //          then
  //              ignore it
  //          else 
  //              return it
  //
  // but it only really needs to consider edges that might intersect 'e'
  // so it can be rewritten as
  //    for each vertex in the defect
  //        for each edge from the vertex
  //            add this edge to a list
  //
  //    for each edge in the list
  //        if this edge shares a vertex with 'e' ignore it
  //        else if this edge intersects 'e' return it
  //
  // This second loop can then be replaced with one that only considers edges that MIGHT intersect 'e'
  // The pGreatArcSet code implements this version...
  //
  //  
  int result = 0;
  
  int old_firstHit       = -1;
  int old_firstHit_vnoLo = -1;
  int old_firstHit_vnoHi = -1;
  
  if (do_old_way) {

    static long stats_tried = 0;

    ROMP_SCOPE_begin
    int i;
    for (i = 0; i < defect->nvertices + defect->nchull; i++) {

      int vno;
      if (i < defect->nvertices) {
        vno = vertex_trans[defect->vertices[i]];
      } else {
        vno = vertex_trans[defect->chull[i - defect->nvertices]];
      }

      if (vno == e->vno1 || vno == e->vno2 || vno < 0) {
        continue;
      }

      EDGE edge2;
      edge2.vno1 = vno;
      
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

      int n;
      for (n = 0; n < vt->vnum; n++) {
        if (vt->v[n] == e->vno1 || vt->v[n] == e->vno2) {
          continue;
        }
        edge2.vno2 = vt->v[n];
        if (defect->optimal_mapping && mris->vertices[vt->v[n]].fixedval == 0) {
          continue;  // experimental
        }
        old_firstHit++;

        stats_tried++;

        if (edgesIntersect(mris, e, &edge2)) {
        
          if (v1) {
            (*v1) = edge2.vno1;
          }
          if (v2) {
            (*v2) = edge2.vno2;
          }
          result = 1;
          
          old_firstHit_vnoLo = edge2.vno1;
          old_firstHit_vnoHi = edge2.vno2;
          if (old_firstHit_vnoLo > old_firstHit_vnoHi) { 
            int tmp = old_firstHit_vnoLo; old_firstHit_vnoLo = old_firstHit_vnoHi; old_firstHit_vnoHi = tmp; 
          }
          
          if (stats_count == 48) {
            fprintf(stderr, "%s:%d stats_count:%ld old found edge vno:%d .. vno:%d\n", __FILE__, __LINE__, stats_count,
                old_firstHit_vnoLo, old_firstHit_vnoHi);
          }
          
          goto Done;
        }
      }
    }
  Done:;

    if (showStats) {
        if (stats_limit < 20000) stats_limit *= 2; else stats_limit += 20000;
        fprintf(stderr, "%s:%d"
          " avg_tried:%g, stats_count:%ld\n", 
          __FILE__, __LINE__,
          (float)stats_tried/(float)stats_count, stats_count);
    }

    ROMP_SCOPE_end
  }
  
  if (do_new_way) {
    static long stats_made;                 // the GreatArcSet
    static long stats_reused;
    static long stats_revised;
    static long stats_numberOfGreatArcs;    // #entries in the GreatArcSet
    
    static long stats_addedVertexs, stats_addedArcs, stats_movedArcs, stats_unmovedArcs;
    static long stats_possible;
    static long stats_tried;
    
    GreatArcSet* gas = ctx->greatArcSet;
    
    if (ctx->nvertices_seen < mris->nvertices) {
#ifdef COPE_WITH_VERTEX_MOVEMENT
        int* p = ctx->vnoToFirstNPlus1 = (int*)realloc(vnoToFirstNPlus1,mris->nvertices*sizeof(int));
#else
        int* p = ctx->vnosHighestNSeen = (int*)realloc(ctx->vnosHighestNSeen,mris->nvertices*sizeof(int));
#endif
        int i; for (i = ctx->nvertices_seen; i < mris->nvertices; i++) p[i] = 0;
        stats_addedVertexs += mris->nvertices - ctx->nvertices_seen;
        ctx->nvertices_seen = mris->nvertices;
    }
    
    // Make a GreatArcSet of all the edges that must be considered
    // unless it has been already made
    //
    bool emit_line_py = false;  // showStats || (stats_count < 45);
    if (emit_line_py) fprintf(stderr, "%s:%d emit_line_py set for stats_count:%ld\n", __FILE__, __LINE__, stats_count);
       
    if (gas && !ctx->obsoleted) {
      stats_reused++;
    } else {

      if (!ctx->obsoleted) {
        stats_made++; 
        gas = ctx->greatArcSet = makeGreatArcSet(mris);
        // emit_line_py = true;
        if (emit_line_py) fprintf(stderr, " [0 ], # line.py reset\n");
      } else {
        stats_revised++;
        ctx->obsoleted = false;
      }

      stats_numberOfGreatArcs = 0;
      
      int i;
      for (i = 0; i < defect->nvertices + defect->nchull; i++) {

        int vno;
        if (i < defect->nvertices) {
          vno = vertex_trans[defect->vertices[i]];
        } else {
          vno = vertex_trans[defect->chull[i - defect->nvertices]];
        }

        if (vno < 0) {
          continue;
        }

        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
        VERTEX          const * const v  = &mris->vertices         [vno];
        
#ifdef COPE_WITH_VERTEX_MOVEMENT
        int                                 nextEntryIndexForVno = ctx->vnoToFirstNPlus1[vno] - 1;
        IntersectDefectEdgesContext_Entry*  prevEntry            = NULL;
#else
        int                                 nextEntryIndexForVno = -1;
#endif
        stats_numberOfGreatArcs += vt->vnum;

        int n;
        for (n = 0; n < vt->vnum; n++) {

          VERTEX const *v2 = &mris->vertices[vt->v[n]];

          // During revision, there are four possibilities... This (vno, n) pair
          //    (1) is still at the same location   - no change
          //    (2) has gone                        - the entry will be ignored by the callback (detected by v->vnum <= n)
          //    (3) did not exist before            - add a new one
          //    (4) has moved                       - tell GreatArcSet it has moved     TESTING HAS NOT FOUND ANY OF THESE
          // Need a per-vno chain of entries that will be cheap
          //
          IntersectDefectEdgesContext_Entry* entry;

#ifdef COPE_WITH_VERTEX_MOVEMENT
          if (nextEntryIndexForVno >= 0) {
            entry = &ctx->entries[nextEntryIndexForVno];
            
            if (v ->cx != entry->cx
            ||  v ->cy != entry->cy
            ||  v ->cz != entry->cz
            ||  v2->cx != entry->cx2
            ||  v2->cy != entry->cy2
            ||  v2->cz != entry->cz2
               ) {                                       // true to tell GreatArcSet it has moved

                static int show_moves_count, show_moves_limit = 1;
                if (show_moves_count++ == show_moves_limit) {
                    if (show_moves_limit < 20) show_moves_limit++;
                    else if (show_moves_limit < 1024) show_moves_limit *= 2;
                    else show_moves_limit += 1024;
                    fprintf(stderr, " (x:%6g,y:%6g,z:%6g) moved by \n (  %6g,  %6g,  %6g)\n",
                        entry->cx,         entry->cy,         entry->cz,
                        entry->cx - v->cx, entry->cy - v->cy, entry->cz - v->cz);
                }

                entry->cx  = v ->cx;
                entry->cy  = v ->cy;
                entry->cz  = v ->cz;
                entry->cx2 = v2->cx;
                entry->cy2 = v2->cy;
                entry->cz2 = v2->cz;

                stats_movedArcs++;

                insertGreatArc(gas, nextEntryIndexForVno, true,
                    v->cx,v->cy,v->cz,  v2->cx,v2->cy,v2->cz);
            } else {
                stats_unmovedArcs++;
            }
            
          }
#else
          if (n < ctx->vnosHighestNSeen[vno]) {
            // ones that have been entered have not moved
            ;
          }
#endif
          else 
          {
          
            if (ctx->entriesSize == ctx->entriesCapacity) {
              ctx->entriesCapacity *= 2;
              if (ctx->entriesCapacity == 0) ctx->entriesCapacity = 128;
              ctx->entries = 
                (IntersectDefectEdgesContext_Entry*)realloc(
                  ctx->entries, 
                  ctx->entriesCapacity*sizeof(IntersectDefectEdgesContext_Entry));
            }
          
            nextEntryIndexForVno = ctx->entriesSize++;
            entry = &ctx->entries[nextEntryIndexForVno];
            entry->vno = vno; entry->n = n; 
#ifdef COPE_WITH_VERTEX_MOVEMENT
            entry->next = -1;
            if (prevEntry) prevEntry->next = nextEntryIndexForVno; else ctx->vnoToFirstNPlus1[vno] = nextEntryIndexForVno + 1;
#endif

#ifdef COPE_WITH_VERTEX_MOVEMENT
            entry->cx  = v ->cx;
            entry->cy  = v ->cy;
            entry->cz  = v ->cz;
            entry->cx2 = v2->cx;
            entry->cy2 = v2->cy;
            entry->cz2 = v2->cz;
#endif  
            bool interesting = false;
            if (false && vno == 151634 && vt->v[n] == 487) {
                fprintf(stderr, "%s:%d , interesting greatArc found\n", __FILE__, __LINE__);
                interesting = true;
            }
            
            insertGreatArc(gas, nextEntryIndexForVno, vno, vt->v[n]);
                
            // get some data for visualizing
            // by a separate tool to understand the situation better
            //
            if (emit_line_py) {
              fprintf(stderr, " [1, %f, %f, %f, %f, %f, %f], # line.py defect edge vno:%d..%d\n", 
                v->cx,v->cy,v->cz,  v2->cx,v2->cy,v2->cz,
                vno, vt->v[n]);
            }
            
          }
          
#ifdef COPE_WITH_VERTEX_MOVEMENT
          nextEntryIndexForVno = entry->next;
          prevEntry            = entry;
#endif
        }   // for n

#ifdef COPE_WITH_VERTEX_MOVEMENT
#else
        if (ctx->vnosHighestNSeen[vno] < vt->vnum) {
            if (ctx->vnosHighestNSeen[vno]) stats_addedArcs += vt->vnum - ctx->vnosHighestNSeen[vno];    // added after the first ones
            ctx->vnosHighestNSeen[vno] = vt->vnum;
        }
#endif
      }
    }
    
    // Get the subset that need to be examined
    //
    VERTEX const * const ev  = &mris->vertices[e->vno1];
    VERTEX const * const ev2 = &mris->vertices[e->vno2];
    
    PossiblyIntersectingGreatArcs_callback_context callback_context;
    {
      bzero(&callback_context, sizeof(callback_context));
      callback_context.mris    = mris;
      callback_context.defect  = defect;
      callback_context.e       = e;
      callback_context.ide_ctx = ctx;
      callback_context.old_firstHit_vnoLo = old_firstHit_vnoLo;
      callback_context.old_firstHit_vnoHi = old_firstHit_vnoHi;
      
      callback_context.emit_line_py = emit_line_py;

      if (emit_line_py) fprintf(stderr, " [2, %f, %f, %f, %f, %f, %f], # line.py target\n", ev->cx,ev->cy,ev->cz, ev2->cx,ev2->cy,ev2->cz);

      possiblyIntersectingGreatArcs(
        gas,
        &callback_context,
        possiblyIntersectingGreatArcs_callback,
        e->vno1,              e->vno2,
        ev->cx,ev->cy,ev->cz, ev2->cx,ev2->cy,ev2->cz, 
        stats_count == stats_limit-1);                    // tracing

    }
        
    stats_possible += stats_numberOfGreatArcs;
    stats_tried    += callback_context.tried;
    
    // Process the found one
    //
    if (callback_context.someMatchingIndexPlus1 > 0) {

      IntersectDefectEdgesContext_Entry* entry = &ctx->entries[callback_context.someMatchingIndexPlus1 - 1];
      int vno = entry->vno;
      int n   = entry->n;
      VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

      if (v1) {
        (*v1) = vno;
      }
      if (v2) {
        (*v2) = v->v[n];
      }
      
      if (do_old_way && !result) { 
        fprintf(stderr, "%s:%d old stats_count:%ld didn't find, new did\n", __FILE__, __LINE__, stats_count);
        *(int*)(-1) = 0; 
      }
      
      result = 1;
      goto Done2;
    }

    if (do_old_way && result) { 
      
      fprintf(stderr, "%s:%d stats_count:%ld old found vno0:%d vno1:%d, new did not during \n", __FILE__, __LINE__, stats_count, 
        old_firstHit_vnoLo, old_firstHit_vnoHi); 
        
      possiblyIntersectingGreatArcs_Debug(                              // show how vno0..vno1 interacts with the arc
        gas,
        ev->cx,ev->cy,ev->cz, ev2->cx,ev2->cy,ev2->cz,                  // the arc
        old_firstHit_vnoLo, old_firstHit_vnoHi);                        // vno0..vno1

      *(int*)(-1) = 0; 
    }

Done2:;
    
    if (showStats) {
        if (stats_limit < 20000) stats_limit *= 2; else stats_limit += 20000;
        fprintf(stderr, "%s:%d"
          " made:%ld revised:%ld reused:%ld"
          " addedVertexs:%g, addedArcs:%g movedArcs:%g movedArcs/revision:%g "
          " unmovedArcs/revision:%g "
          " avg_possible:%g  avg_tried:%g  stats_count:%ld\n", 
          __FILE__, __LINE__,
          stats_made, stats_revised, stats_reused,
          (float)stats_addedVertexs,
          (float)stats_addedArcs,
          (float)stats_movedArcs,
          (float)stats_movedArcs/(float)stats_revised,
          (float)stats_unmovedArcs/(float)stats_revised,
          (float)stats_possible/(float)stats_count, 
          (float)stats_tried   /(float)stats_count,
          stats_count);
    }
  }

  return result;
}

static int intersectDefectConvexHullEdges(
    MRI_SURFACE *mris, DEFECT *defect, EDGE *e, IntersectDefectEdgesContext* ctx, int *vertex_trans, int *v1, int *v2)
{
  int i, n, vno;
  EDGE edge2;

  for (i = defect->nborder; i < defect->nchull; i++) {
    vno = vertex_trans[defect->chull[i]];
    if (vno == e->vno1 || vno == e->vno2 || vno < 0) {
      continue;
    }

    edge2.vno1 = vno;
    VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];
    for (n = 0; n < v->vnum; n++) {
      if (v->v[n] == e->vno1 || v->v[n] == e->vno2) {
        continue;
      }
      edge2.vno2 = v->v[n];
      if (defect->optimal_mapping && mris->vertices[v->v[n]].fixedval == 0) {
        continue;  // experimental
      }
      if (edgesIntersect(mris, e, &edge2)) {
        if (v1) {
          (*v1) = edge2.vno1;
        }
        if (v2) {
          (*v2) = edge2.vno2;
        }
        return (1);
      }
    }
  }

  return (0);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Add the appropriate face to the tessellation.
  ------------------------------------------------------*/
static int mrisAddAllDefectFaces(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vertex_trans)
{
  DEFECT *defect;
  int i, vno1, vno2, nfaces, n, m, dno;

#if (!SPHERE_INTERSECTION)
  double origin[3], e0[3], e1[3];
#endif
  EDGE edge;

  for (dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    for (i = 0; i < defect->nborder + defect->nvertices; i++) {
      if (i < defect->nvertices) {
        if (defect->status[i] == DISCARD_VERTEX) {
          continue;
        }
        vno1 = vertex_trans[defect->vertices[i]];
      }
      else {
        vno1 = vertex_trans[defect->border[i - defect->nvertices]];
      }
      if (vno1 < 0) {
        continue;
      }

      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno1];
      edge.vno1 = vno1;
      for (m = 0; m < vt->vnum; m++) {
        edge.vno2 = vno2 = vt->v[m];
#if (!SPHERE_INTERSECTION)
        mrisComputeCanonicalEdgeBasis(mris, &edge, &edge, origin, e0, e1);
#endif
        if (vno1 == 108332 && vno2 == 109240) {
          DiagBreak();
        }

        /*
          for every vertex which is a neighbor of both of these,
          add 1 triangle */
        for (nfaces = n = 0; n < vt->vnum; n++) {
          if (vt->v[n] == vno2) {
            continue;
          }
          if (vertexNeighbor(mris, vno2, vt->v[n]) && !isFace(mris, vno1, vno2, vt->v[n]) &&
#if SPHERE_INTERSECTION
              !containsAnotherVertexOnSphere(mris, vno1, vno2, vt->v[n], 0)) {
#else
              !containsAnotherVertex(mris, vno1, vno2, vt->v[n], e0, e1, origin)) {
#endif
            if (nfaces++ > 1) {
              DiagBreak();
            }
            if (mris->nfaces == Gdiag_no) {
              DiagBreak();
            }
            mrisAddFace(mris, vno1, vno2, vt->v[n]);
          }
        }
      }
    }
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Orient the faces of the tessellation so that the
  point outward (i.e. in the same direction as the original
  surface).
  ------------------------------------------------------*/
#if 1
static int mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vtrans)
{
  int vno, n, fno, m, vno0, vno1; /*n0, n1;*/
  FACE *f;
  float dot, norm[3], cx, cy, cz, a[3], b[3];

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);

  MRIScomputeMetricProperties(mris);

  for (fno = 0; fno < mris->nfaces; fno++) {
    if (fno == Gdiag_no) {
      DiagBreak();
    }
    f = &mris->faces[fno];
    VERTEX const * const v1 = &mris->vertices[f->v[0]];
    VERTEX const * const v2 = &mris->vertices[f->v[1]];
    VERTEX const * const v3 = &mris->vertices[f->v[2]];
    cx = v1->cx + v2->cx + v3->cx;
    cy = v1->cy + v2->cy + v3->cy;
    cz = v1->cz + v2->cz + v3->cz;

    a[0] = v2->cx - v1->cx;
    b[0] = v3->cx - v1->cx;
    a[1] = v2->cy - v1->cy;
    b[1] = v3->cy - v1->cy;
    a[2] = v2->cz - v1->cz;
    b[2] = v3->cz - v1->cz;

    F_CROSS(a, b, norm);

    dot = norm[0] * cx + norm[1] * cy + norm[2] * cz;
    if (dot < 0) /* they disagree - change order of vertices in face */
    {
      vno0 = f->v[1];
      vno1 = f->v[2];
      f->v[1] = vno1;
      f->v[2] = vno0;
      mrisSetVertexFaceIndex(mris, vno0, fno);
      mrisSetVertexFaceIndex(mris, vno1, fno);
    }
  }

  mrisCheckVertexFaceTopology(mris);
  
  MRIScomputeMetricProperties(mris);

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    dot = v->nx * v->x + v->ny * v->y + v->nz * v->z;
    if (dot < 0) {
      fprintf(stdout, "vertex %d seems to have inverted normal!\n", vno);
      DiagBreak();
    }

    for (m = 0; m < vt->num; m++) {
      fno = vt->f[m];
      f = &mris->faces[fno];
      if (fno == Gdiag_no) {
        DiagBreak();
      }
      if (f->area < 0) {
        fprintf(stdout, "face %d seems to have negative area!\n", fno);
        DiagBreak();
      }
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mrisNormalFace(mris, fno, n, norm); /* compute face normal */
        dot = norm[0] * v->x + norm[1] * v->y + norm[2] * v->z;

        if (dot < 0) /* they disagree - change order
                                                  of vertices in face */
        {
          fprintf(stdout, "face %d seems to have inverted normal!\n", fno);
          DiagBreak();
        }
      }
    }
  }

  MRISclearMarks(mris);
  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);
  return (NO_ERROR);
}
#else
static int mrisOrientRetessellatedSurface(MRI_SURFACE *mris, DEFECT_LIST *dl, int *vtrans)
{
#if 1
  int dno, fno, vno, i, n, m, vno0, vno1, n0, n1, oriented, nreversed;
  VERTEX *v, *vn;
  FACE *f;
  DEFECT *defect;
  float dot, norm[3], len, nx, ny, nz;

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRISrestoreVertexPositions(mris, ORIG_VERTICES);
  MRISaverageVertexPositions(mris, 200);
  MRIScomputeMetricProperties(mris);

  /* compute average boundary normal in the smoothed space */
  for (dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    nx = ny = nz = 0.0f;
    for (i = 0; i < defect->nborder; i++) {
      vno = vtrans[defect->border[i]];
      if (vno < 0) /* not in new tessellation */
      {
        continue;
      }
      v = &mris->vertices[vno];
      nx += v->nx;
      ny += v->ny;
      nz += v->nz;
    }
    len = sqrt(nx * nx + ny * ny + nz * nz);
    if (FZERO(len)) {
      len = 1.0f;
    }
    defect->nx = nx / len;
    defect->ny = ny / len;
    defect->nz = nz / len;
  }

  /* orient defect faces so that the outward direction agrees with boundary */
  for (oriented = dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    for (i = 0; i < defect->nvertices; i++) {
      vno = vtrans[defect->vertices[i]];
      if (vno < 0) /* not in new tessellation */
      {
        continue;
      }
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }

      /* go through each face and orient it to agree with defect normal */
      for (m = 0; m < v->num; m++) {
        fno = v->f[m];
        f = &mris->faces[fno];
        if (fno == Gdiag_no) {
          DiagBreak();
        }
        for (n = 0; n < VERTICES_PER_FACE; n++) {
          mrisNormalFace(mris, fno, n, norm); /*compute face normal */
          dot = norm[0] * defect->nx + norm[1] * defect->ny + norm[2] * defect->nz;
          if (dot < 0) /* they disagree - change
          order of vertices in face */
          {
            oriented++;
            n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
            n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
            vno0 = f->v[n0];
            vno1 = f->v[n1];
            f->v[n0] = vno1;
            f->v[n1] = vno0;
            mrisSetVertexFaceIndex(mris, vno0, fno);
            mrisSetVertexFaceIndex(mris, vno1, fno);
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stdout, "reversing face %d orientation\n", fno);
          }
        }
      }
    }
  }

  MRISsetNeighborhoodSizeAndDist(mris, 2);

  i = 0;
  do {
    MRIScomputeMetricProperties(mris);
    nreversed = 0;
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      for (nx = ny = nz = 0.0, n = 0; n < v->vtotal; n++) {
        vn = &mris->vertices[v->v[n]];
        dot = vn->nx * v->nx + vn->ny * v->ny + vn->nz * v->nz;
        if (dot < 0) {
          DiagBreak();
        }
        nx += vn->nx;
        ny += vn->ny;
        nz += vn->nz;
      }
      dot = nx * v->nx + ny * v->ny + nz * v->nz;
      if (dot < 0) {
        v->nx *= -1.0;
        v->ny *= -1.0;
        v->nz *= -1.0;
        DiagBreak();
      }

      for (m = 0; m < v->num; m++) {
        fno = v->f[m];
        f = &mris->faces[fno];
        if (fno == Gdiag_no) {
          DiagBreak();
        }
        dot = nx * f->nx + ny * f->ny + nz * f->nz;
        if (dot < 0) /* they disagree - change order
        of vertices in face */
        {
          DiagBreak();
        }
        for (n = 0; n < VERTICES_PER_FACE; n++) {
          mrisNormalFace(mris, fno, n, norm); /* compute
         face normal */
          dot = norm[0] * nx + norm[1] * ny + norm[2] * nz;
          if (dot < 0) /* they disagree - change
          order of vertices in face */
          {
            nreversed++;
            n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
            n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
            vno0 = f->v[n0];
            vno1 = f->v[n1];
            f->v[n0] = vno1;
            f->v[n1] = vno0;
            mrisSetVertexFaceIndex(mris, vno0, fno);
            mrisSetVertexFaceIndex(mris, vno1, fno);
            if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) fprintf(stdout, "reversing face %d orientation\n", fno);
          }
        }
      }
    }
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "\rpass %d: %d oriented   ", i + 1, nreversed);
    }

    oriented += nreversed;
    if (++i > 20) /* shouldn't happen, but... */
    {
      break;
    }
  } while (nreversed > 0);

  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "orientation complete in %d passes\n", i);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (m = 0; m < v->num; m++) {
      fno = v->f[m];
      f = &mris->faces[fno];
      if (fno == Gdiag_no) {
        DiagBreak();
      }
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mrisNormalFace(mris, fno, n, norm); /* compute face normal */
        dot = norm[0] * v->nx + norm[1] * v->ny + norm[2] * v->nz;
        if (dot < 0) /* they disagree - change order
        of vertices in face */
        {
          DiagBreak();
        }
        dot = norm[0] * f->nx + norm[1] * f->ny + norm[2] * f->nz;
        if (dot < 0) /* they disagree - change order
        of vertices in face */
        {
          DiagBreak();
        }
      }
    }
  }

#if 0
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  /* back to inflated */
#endif
  MRIScomputeMetricProperties(mris);
  return (oriented);
#else
  int vno, vno0, vno1, i, n, fno, oriented, num, dno, m, blist[MAX_DEFECT_VERTICES], nb, tmp[MAX_DEFECT_VERTICES],
      nbnew, n0, n1;
  FACE *f;
  VERTEX *v, *vn;
  float norm[3], dot, len;
  DEFECT *defect;

  oriented = 0;
  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);
  MRISsetNeighborhoodSizeAndDist(mris, 2);

  /* first orient vertices */

  /* mark all defective vertices */
  for (nb = dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    for (n = 0; n < defect->nvertices; n++) {
      vno = vtrans[defect->vertices[n]];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (vno < 0) {
        continue;
      }
      mris->vertices[vno].marked = 1;
    }
    for (n = 0; n < defect->nborder; n++) {
      mris->vertices[vtrans[defect->border[n]]].marked = 1;
    }
    memmove(blist + nb, defect->border, defect->nborder * sizeof(int));
    nb += defect->nborder;
  }

  /* starts out as a list of border vertices and will grow inwards */
  for (i = 0; i < nb; i++) {
    blist[i] = vtrans[blist[i]];
  }

  do /* grow border inwards one edge length at each iteration */
  {
    for (nbnew = i = 0; i < nb; i++) {
      vno = blist[i];
      if (vno < 0) /* not in new tessellation - shouldn't happen */
      {
        continue;
      }
      v = &mris->vertices[vno];
      if (vno == Gdiag_no || vno == 135681) {
        DiagBreak();
      }
      v->nx = v->ny = v->nz = 0;
      for (dot = 0.0, num = n = 0; n < v->vtotal; n++) {
        vn = &mris->vertices[v->v[n]];

        /* only use already oriented (or non-defective) vertices */
        if (vn->marked) {
          continue;
        }
        v->nx += vn->nx;
        v->ny += vn->ny;
        v->nz += vn->nz;
        num++;
      }
      if (!num) /* surrounded by unoriented defective vertices */
      {
        tmp[nbnew++] = vno; /* so it will be processed
   next time again */
        continue;
      }
      len = sqrt(v->nx * v->nx + v->ny * v->ny + v->nz * v->nz);
      if (FZERO(len)) {
        len = 1.0f;
      }
      v->nx /= len;
      v->ny /= len;
      v->nz /= len;
      v->marked = 3; /* it's proper orientation has been established */
    }
    for (i = 0; i < nb; i++) {
      vno = blist[i];
      if (vno < 0) /* not in new tessellation - shouldn't happen */
      {
        continue;
      }
      v = &mris->vertices[vno];
      if (v->marked == 3) /* was oriented properly */
      {
        v->marked = 0;
      }
    }

    /* now build new list of vertices, moving inward by one vertex */
    for (i = 0; i < nb; i++) {
      vno = blist[i];
      v = &mris->vertices[vno];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      for (n = 0; n < v->vnum; n++) {
        vn = &mris->vertices[v->v[n]];

        /* only use already oriented (or non-defective) vertices */
        if (vn->marked == 1) {
          vn->marked = 2; /* don't add it more than once */
          tmp[nbnew++] = v->v[n];
        }
      }
    }
    nb = nbnew;
    if (nb > 0) {
      memmove(blist, tmp, nb * sizeof(int));
    }
  } while (nb > 0);

  MRISclearMarks(mris);

  /* now orient faces to agree with their vertices */
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (fno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      mrisNormalFace(mris, fno, n, norm); /* how about vertex 2 ???? */
      vno = f->v[n];
      v = &mris->vertices[vno];
      dot = norm[0] * v->nx + norm[1] * v->ny + norm[2] * v->nz;
      if (dot < 0) /* change order of vertices in face */
      {
        oriented++;
        n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
        n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
        vno0 = f->v[n0];
        vno1 = f->v[n1];
        f->v[n0] = vno1;
        f->v[n1] = vno0;
        mrisSetVertexFaceIndex(mris, vno0, fno);
        mrisSetVertexFaceIndex(mris, vno1, fno);
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stdout, "reversing face %d orientation\n", fno);
        }
      }
    }
  }

  MRIScomputeMetricProperties(mris);
  MRISclearMarks(mris);

  /* mark all vertices that have a normal which disagrees with any neighbor */
  for (dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    for (i = 0; i < defect->nvertices; i++) {
      vno = vtrans[defect->vertices[i]];
      if (vno < 0) {
        continue;
      }
      v = &mris->vertices[vno];
      for (n = 0; n < v->vnum; n++) {
        vn = &mris->vertices[v->v[n]];
        dot = vn->nx * v->nx + vn->ny * v->ny + vn->nz * v->nz;
        if (dot < 0) {
          v->marked = 1;
          if (!vn->marked) {
            vn->marked = 1;
          }
          break;
        }
      }
    }
  }

  /* go back and orient the ambiguous vertices based on the normal of
  the defect.
  */
  for (dno = 0; dno < dl->ndefects; dno++) {
    defect = &dl->defects[dno];
    for (i = 0; i < defect->nvertices; i++) {
      vno = vtrans[defect->vertices[i]];
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      if (vno < 0) {
        continue;
      }
      v = &mris->vertices[vno];
      if (!v->marked) {
        continue;
      }
      dot = defect->nx * v->nx + defect->ny * v->ny + defect->nz * v->nz;
      if (dot < 0) {
        v->nx *= -1;
        v->ny *= -1;
        v->nz *= -1;
      }
    }
  }
  MRISclearMarks(mris);

  /* now orient faces to agree with their vertices */
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    if (fno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      mrisNormalFace(mris, fno, n, norm); /* how about vertex 2 ???? */
      vno = f->v[n];
      v = &mris->vertices[vno];
      dot = norm[0] * v->nx + norm[1] * v->ny + norm[2] * v->nz;
      if (dot < 0) /* change order of vertices in face */
      {
        oriented++;
        m = (n + 1) >= VERTICES_PER_FACE ? 0 : n + 1;
        vno0 = f->v[n];
        vno1 = f->v[m];
        f->v[n] = vno1;
        f->v[m] = vno0;
        mrisSetVertexFaceIndex(mris, vno0, fno);
        mrisSetVertexFaceIndex(mris, vno1, fno);
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stdout, "reversing face %d orientation\n", fno);
        }
        mrisNormalFace(mris, fno, n, norm); /* how about
    vertex 2 ???? */
        dot = norm[0] * v->nx + norm[1] * v->ny + norm[2] * v->nz;
      }
    }
  }

  MRIScomputeMetricProperties(mris);

  if (Gdiag_no >= 0) {
    v = &mris->vertices[Gdiag_no];
    for (dot = 0.0, n = 0; n < v->vnum; n++) {
      vn = &mris->vertices[v->v[n]];

      /* only use already oriented (or non-defective) vertices */
      dot += v->nx * vn->nx + v->ny * vn->ny + v->nz * vn->nz;
    }
    if (dot < 0) {
      DiagBreak();
    }
    for (dot = 0.0, n = 0; n < v->num; n++) {
      fno = v->f[n];
      f = &mris->faces[fno];

      /* only use already oriented (or non-defective) vertices */
      dot += v->nx * f->nx + v->ny * f->ny + v->nz * f->nz;
    }
    if (dot < 0) {
      DiagBreak();
    }
    MRISrestoreVertexPositions(mris, ORIG_VERTICES);
    MRIScomputeMetricProperties(mris);
    v = &mris->vertices[Gdiag_no];
    for (dot = 0.0, n = 0; n < v->vnum; n++) {
      vn = &mris->vertices[v->v[n]];

      /* only use already oriented (or non-defective) vertices */
      dot += v->nx * vn->nx + v->ny * vn->ny + v->nz * vn->nz;
    }
    if (dot < 0) {
      DiagBreak();
    }
  }

#if 0
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
#endif
  return (oriented);
#endif
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  find the convex hull of the defect on the sphere (actually,
  just a disk, but at least it's convex....)
  ------------------------------------------------------*/
static int mrisFindDefectConvexHull(MRI_SURFACE *mris, DEFECT *defect)
{
#if SMALL_CONVEX_HULL
  int chull[MAX_DEFECT_VERTICES], nfound, n, i, vno;

  defect->chull = chull;
  defect->nchull = defect->nborder; /* include border vertices first*/
  memmove(chull, defect->border, defect->nborder * sizeof(int));

  MRISclearMarks(mris);
  mrisMarkDefectConvexHull(mris, defect, 1);
  mrisMarkDefect(mris, defect, 1);

  for (nfound = i = 0; i < defect->nborder; i++) {
    vno = defect->border[i];
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

    /* vertex inside convex hull - add all its nbrs */
    for (n = 0; n < vt->vnum; n++) {
      VERTEX * const vn = &mris->vertices[vt->v[n]];
      if (vn->marked) /* already in defect or convex hull */
      {
        continue;
      }
      chull[defect->nchull + nfound++] = vt->v[n]; /* first neighbors only! */
      vn->marked = 1;
    }
  }
  defect->nchull += nfound;

  MRISclearMarks(mris);
  defect->chull = (int *)calloc(defect->nchull, sizeof(int));
  if (!defect->chull) ErrorExit(ERROR_NO_MEMORY, "mrisFindConvexHull: could not allocate %d vlist\n", defect->nchull);
  memmove(defect->chull, chull, defect->nchull * sizeof(int));

  return (NO_ERROR);
#else
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int chull[MAX_DEFECT_VERTICES], nfound, n, i, vno;

  xmin = ymin = zmin = 100000;
  xmax = ymax = zmax = 0.0f;

  /* now compute max radius on surface of sphere */
  for (i = 0; i < defect->nvertices + defect->nborder; i++) {
    if (i < defect->nvertices) {
      vno = defect->vertices[i];
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
    }
    else {
      vno = defect->border[i - defect->nvertices];
    }
    VERTEX const * const v  = &mris->vertices[vno];
    if (v->cx < xmin) {
      xmin = v->cx;
    }
    if (v->cy < ymin) {
      ymin = v->cy;
    }
    if (v->cz < zmin) {
      zmin = v->cz;
    }
    if (v->cx > xmax) {
      xmax = v->cx;
    }
    if (v->cy > ymax) {
      ymax = v->cy;
    }
    if (v->cz > zmax) {
      zmax = v->cz;
    }
  }
  defect->chull = chull;
  defect->nchull = defect->nborder;
  memmove(chull, defect->border, defect->nborder * sizeof(int));

  MRISclearMarks(mris);
  mrisMarkDefectConvexHull(mris, defect, 1);
  mrisMarkDefect(mris, defect, 1);

  do {
    nfound = 0;
    for (i = 0; i < defect->nchull; i++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->chull[i]];
      VERTEX          const * const v  = &mris->vertices         [defect->chull[i]];
      if (defect->chull[i] == Gdiag_no) {
        DiagBreak();
      }
      if ((v->cx >= xmin && v->cx <= xmax) && (v->cy >= ymin && v->cy <= ymax) && (v->cz >= zmin && v->cz <= zmax)) {
        /* vertex inside convex hull - add all its nbrs */
        for (n = 0; n < vt->vnum; n++) {
          VERTEX * const vn = &mris->vertices[vt->v[n]];
          if (vn->marked) /* already in defect or convex hull */
          {
            continue;
          }
          chull[defect->nchull + nfound++] = vt->v[n];
          vn->marked = 1;
        }
      }
    }
    defect->nchull += nfound;
  } while (nfound > 0);

  MRISclearMarks(mris);

  defect->chull = (int *)calloc(defect->nchull, sizeof(int));
  if (!defect->chull) ErrorExit(ERROR_NO_MEMORY, "mrisFindConvexHull: could not allocate %d vlist\n", defect->nchull);
  memmove(defect->chull, chull, defect->nchull * sizeof(int));

  return (NO_ERROR);
#endif
}


static int mrisDefectRemoveNegativeVertices(MRI_SURFACE *mris, DEFECT *defect)
{
  int i, n;
  for (i = 0; i < defect->nvertices; i++) {
    if (defect->status[i] == DISCARD_VERTEX) {
      continue;
    }
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[defect->vertices[i]];
    for (n = 0; n < vt->num; n++)
      if (mris->faces[vt->f[n]].area < 0.0) {
        defect->status[i] = DISCARD_VERTEX;
      }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisDefectRemoveDegenerateVertices(MRI_SURFACE *mris, float min_sphere_dist, DEFECT *defect)
{
  float dx, dy, dz, dist;
  int i, j;
  VERTEX *v, *vn;

  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0; i < defect->nvertices + defect->nborder; i++) {
    if (i < defect->nvertices) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      v = &mris->vertices[defect->vertices[i]];
      if (defect->vertices[i] == Gdiag_no) {
        DiagBreak();
      }
    }
    else {
      v = &mris->vertices[defect->border[i - defect->nvertices]];
      if (defect->border[i - defect->nvertices] == Gdiag_no) {
        DiagBreak();
      }
    }
    for (j = i + 1; j < defect->nvertices; j++) {
      if (defect->status[j] == DISCARD_VERTEX) {
        continue;
      }
      vn = &mris->vertices[defect->vertices[j]];
      dx = vn->cx - v->cx;
      dy = vn->cy - v->cy;
      dz = vn->cz - v->cz;
      dist = (dx * dx + dy * dy + dz * dz); /* no sqrt */
      if (dist < min_sphere_dist) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n", defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX;
      }
    }
  }
  return (NO_ERROR);
}
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int mrisDefectRemoveProximalVertices(MRI_SURFACE *mris, float min_orig_dist, DEFECT *defect)
{
  float dx, dy, dz, dist;
  int i, j;
  VERTEX *v, *vn;

  /* discard vertices that are too close to another vertex on sphere */
  for (i = 0; i < defect->nvertices + defect->nborder; i++) {
    if (i < defect->nvertices) {
      if (defect->status[i] == DISCARD_VERTEX) {
        continue;
      }
      v = &mris->vertices[defect->vertices[i]];
      if (defect->vertices[i] == Gdiag_no) {
        DiagBreak();
      }
    }
    else {
      v = &mris->vertices[defect->border[i - defect->nvertices]];
      if (defect->border[i - defect->nvertices] == Gdiag_no) {
        DiagBreak();
      }
    }

    for (j = i + 1; j < defect->nvertices; j++) {
      if (defect->status[j] == DISCARD_VERTEX) {
        continue;
      }
      vn = &mris->vertices[defect->vertices[j]];
      dx = vn->origx - v->origx;
      dy = vn->origy - v->origy;
      dz = vn->origz - v->origz;
      dist = (dx * dx + dy * dy + dz * dz); /* no sqrt */
      if (dist < min_orig_dist) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "discarding proximal vertex %d\n", defect->vertices[j]);
        defect->status[j] = DISCARD_VERTEX;
      }
    }
  }
  return (NO_ERROR);
}

/*!
  \fn int MRISsegmentMarked(MRI_SURFACE *mris, LABEL ***plabel_array, int *pnlabels, float min_label_area)
  \brief Appears to create a label for each connected component
  defined by v->marked=1; the surface area of the label must be >
  min_label_area.
*/
int MRISsegmentMarked(MRI_SURFACE *mris, LABEL ***plabel_array, int *pnlabels, float min_label_area)
{
  int vno, nfound, n, nlabels, *marks;
  VERTEX *v;
  LABEL *area = NULL, **tmp, **label_array;

  marks = (int *)calloc(mris->nvertices, sizeof(int));
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *));
  if (!label_array || !marks)
    ErrorExit(ERROR_NOMEMORY, "%s: MRISsegmentMarked could not allocate tmp storage", Progname);

  /* save current marks */
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    marks[vno] = v->marked;
    if (v->marked != 0) {
      v->marked = 1;
    }
  }

  nlabels = 0;
  do {
    nfound = 0;

    v = &mris->vertices[0];

    /* find a marked vertex */
    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag || v->marked != 1) {
        continue;
      }
      break;
    }

    if (vno < mris->nvertices) {
      area = LabelAlloc(mris->nvertices, NULL, NULL);
      area->n_points = 1;
      area->lv[0].x = v->x;
      area->lv[0].y = v->y;
      area->lv[0].z = v->z;
      area->lv[0].vno = vno;
      LabelFillMarked(area, mris);
      if (LabelArea(area, mris) >= min_label_area) {
        label_array[nlabels++] = LabelCopy(area, NULL);
      }
      LabelFree(&area);
      nfound = 1;
    }
    else {
      nfound = 0;
    }

  } while (nfound > 0);

  /* restore original marks */
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->marked = marks[vno];
  }

  free(marks);

  /* crunch label array down to a reasonable size */
  tmp = label_array;
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *));
  if (!label_array) ErrorExit(ERROR_NOMEMORY, "%s: MRISsegmentMarked could not allocate tmp storage", Progname);
  for (n = 0; n < nlabels; n++) {
    label_array[n] = tmp[n];
  }
  free(tmp);
  *plabel_array = label_array;
  *pnlabels = nlabels;
  return (NO_ERROR);
}

int MRISsegmentAnnotated(MRI_SURFACE *mris, LABEL ***plabel_array, int *pnlabels, float min_label_area)
{
  int vno, nfound, n, nlabels, last_vno;
  VERTEX *v;
  LABEL *area = NULL, **tmp, **label_array;

  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *));
  if (!label_array) ErrorExit(ERROR_NOMEMORY, "%s: MRISsegmentAnnotated could not allocate tmp storage", Progname);

  MRISclearMarks(mris);

  nlabels = 0;
  last_vno = -1;
  do {
    nfound = 0;

    v = &mris->vertices[0];

    /* find an un-marked vertex */
    for (vno = last_vno + 1; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag || v->annotation == 0 || v->marked) {
        continue;
      }
      break;
    }
    if (vno < mris->nvertices) {
      area = LabelAlloc(mris->nvertices, NULL, NULL);
      area->n_points = 1;
      area->lv[0].x = v->x;
      area->lv[0].y = v->y;
      area->lv[0].z = v->z;
      area->lv[0].vno = vno;
      LabelFillAnnotated(area, mris);
      if (LabelArea(area, mris) >= min_label_area) {
        label_array[nlabels++] = LabelCopy(area, NULL);
      }
      LabelMarkSurface(area, mris);
      LabelFree(&area);
      nfound = 1;
      last_vno = vno;
    }
    else {
      nfound = 0;
    }

  } while (nfound > 0);

  /* crunch label array down to a reasonable size */
  tmp = label_array;
  label_array = (LABEL **)calloc(mris->nvertices, sizeof(LABEL *));
  if (!label_array) ErrorExit(ERROR_NOMEMORY, "%s: MRISsegmentAnnotated could not allocate tmp storage", Progname);
  for (n = 0; n < nlabels; n++) {
    label_array[n] = tmp[n];
  }
  free(tmp);
  *plabel_array = label_array;
  *pnlabels = nlabels;
  return (NO_ERROR);
}

static int mrisComputeGrayWhiteBorderDistributions(MRI_SURFACE *mris,
                                                   MRI *mri,
                                                   DEFECT *defect,
                                                   HISTOGRAM *h_white,
                                                   HISTOGRAM *h_gray,
                                                   HISTOGRAM *h_border,
                                                   HISTOGRAM *h_grad)
{
  int vno, n2, n, i, nvertices, bin;
  HISTOGRAM *h_white_raw, *h_gray_raw, *h_border_raw, *h_grad_raw;
  double grad, min_grad, max_grad, bin_val, bin_size;

  HISTOclear(h_gray, h_gray);
  HISTOclear(h_white, h_white);
  HISTOclear(h_border, h_border);
  HISTOclear(h_grad, h_grad);

  mrisMarkDefect(mris, defect, 1); /* avoid vertices in the defect */
  for (bin = 0; bin < h_gray->nbins; bin++) {
    h_gray->bins[bin] = h_white->bins[bin] = h_border->bins[bin] = bin;
  }

  min_grad = 100000;
  max_grad = -min_grad;

  for (nvertices = i = 0; i < defect->nchull; i++) {
    vno = defect->chull[i];
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      for (n2 = 0; n2 < vnt->vtotal; n2++) {
        VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
        if (vn2->marked) /* already processed */
        {
          continue;
        }
        grad = vn2->val2 - vn2->val2bak;
        if (grad < min_grad) {
          min_grad = grad;
        }
        if (grad > max_grad) {
          max_grad = grad;
        }
      }
    }
  }

  /* add one bin at either end for almost zero probability events */
  bin_size = (max_grad - min_grad) / (h_grad->nbins - 2);
  h_grad->bin_size = bin_size;
  for (bin_val = min_grad - bin_size, bin = 0; bin < h_grad->nbins; bin++, bin_val += bin_size) {
    h_grad->bins[bin] = bin_val;
  }

  min_grad = h_grad->bins[0];

  for (nvertices = i = 0; i < defect->nchull; i++) {
    vno = defect->chull[i];
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      for (n2 = 0; n2 < vnt->vtotal; n2++) {
        VERTEX * const vn2 = &mris->vertices[vnt->v[n2]];
        if (vn2->marked) /* already processed */
        {
          continue;
        }
        nvertices++;

        if (vn2->val2bak < 70) {
          DiagBreak();
        }

        bin = nint(vn2->val2);
        bin = MIN(h_white->nbins, MAX(0, bin));
        h_white->counts[bin]++; /* wm value */
        bin = nint(vn2->val2bak);
        bin = MIN(h_gray->nbins, MAX(0, bin));
        h_gray->counts[bin]++; /* gray value */
        bin = nint(vn2->val);
        bin = MIN(h_border->nbins, MAX(0, bin));
        h_border->counts[bin]++; /* border value */
        grad = vn2->val2 - vn2->val2bak;
        bin = (int)((grad - min_grad) / bin_size);
        bin = MIN(h_grad->nbins, MAX(0, bin));
        h_grad->counts[bin]++;

        vn2->marked = 1; /* don't process it twice */
      }
    }
  }

  /* unmark them all */
  for (i = 0; i < defect->nchull; i++) {
    vno = defect->chull[i];
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      for (n2 = 0; n2 < vnt->vtotal; n2++) {
        VERTEX * const vn2 = &mris->vertices[vnt->v[n2]];
        vn2->marked = 0;
      }
    }
  }

  for (bin = 0; bin < h_gray->nbins; bin++) {
    if (h_gray->counts[bin] == 0) {
      h_gray->counts[bin] = 0.1;
    }
    if (h_white->counts[bin] == 0) {
      h_white->counts[bin] = 0.1;
    }
    if (h_border->counts[bin] == 0) {
      h_border->counts[bin] = 0.1;
    }
    if (h_grad->counts[bin] == 0) {
      h_grad->counts[bin] = 0.1;
    }
    h_grad->counts[bin] /= (float)nvertices;
    h_gray->counts[bin] /= (float)nvertices;
    h_white->counts[bin] /= (float)nvertices;
    h_border->counts[bin] /= (float)nvertices;
  }

  h_grad_raw = HISTOcopy(h_grad, NULL);
  h_gray_raw = HISTOcopy(h_gray, NULL);
  h_white_raw = HISTOcopy(h_white, NULL);
  h_border_raw = HISTOcopy(h_border, NULL);

  // to correct the bug in HISTOsmooth
  h_grad_raw->bin_size = h_grad->bin_size;
  h_gray_raw->bin_size = h_gray->bin_size;
  h_white_raw->bin_size = h_white->bin_size;
  h_border_raw->bin_size = h_border->bin_size;

  HISTOsmooth(h_grad_raw, h_grad, 2);
  HISTOsmooth(h_gray_raw, h_gray, 2);
  HISTOsmooth(h_white_raw, h_white, 2);
  HISTOsmooth(h_border_raw, h_border, 2);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    HISTOplot(h_gray, "g.plt");
    HISTOplot(h_white, "w.plt");
    HISTOplot(h_border, "b.plt");
    HISTOplot(h_gray_raw, "gr.plt");
    HISTOplot(h_white_raw, "wr.plt");
    HISTOplot(h_border_raw, "br.plt");
    HISTOplot(h_grad, "d.plt");
    HISTOplot(h_grad_raw, "dr.plt");
  }
  mrisMarkDefect(mris, defect, 0);
  HISTOfree(&h_gray_raw);
  HISTOfree(&h_white_raw);
  HISTOfree(&h_border_raw);
  HISTOfree(&h_grad_raw);
  return (NO_ERROR);
}


static void mrisComputeSurfaceStatistics(
    MRIS *mris, MRI *mri, HISTOGRAM *h_k1, HISTOGRAM *h_k2, MRI *mri_k1_k2, MRI *mri_gray_white, HISTOGRAM *h_dot)
{
  int n, nvertices, nfaces;
  TP tp;
  float total_ll;

  TPinit(&tp);

  fprintf(WHICH_OUTPUT, "Computing Initial Surface Statistics\n");

  MRISsaveVertexPositions(mris, TMP_VERTICES);
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
  MRIScomputeMetricProperties(mris);
  MRIScomputeSecondFundamentalForm(mris);

  /* initializing the table of vertices */
  nvertices = mris->nvertices;
  tp.vertices = (int *)malloc(nvertices * sizeof(int));
  for (nvertices = n = 0; n < mris->nvertices; n++) {
    if (mris->vertices[n].marked) {
      continue;
    }
    tp.vertices[nvertices++] = n;
  }
  tp.nvertices = nvertices;

  /* initializing the table of faces */
  nfaces = mris->nfaces;
  tp.faces = (int *)malloc(nfaces * sizeof(int));
  for (nfaces = n = 0; n < mris->nfaces; n++) {
    if (triangleMarked(mris, n)) {
      continue;
    }
    tp.faces[nfaces++] = n;
  }
  tp.nfaces = nfaces;

  mrisComputeDefectCurvatureLogLikelihood(mris, &tp, h_k1, h_k2, mri_k1_k2);
  mrisComputeDefectNormalDotLogLikelihood(mris, &tp, h_dot);
  mrisDefectFaceMRILogLikelihood(mris, mri, &tp, NULL, NULL, NULL, mri_gray_white);
  mrisDefectVertexMRILogLikelihood(mris, mri, &tp, NULL, NULL, NULL, mri_gray_white);
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);

  total_ll = (tp.face_ll + tp.vertex_ll + tp.curv_ll + tp.qcurv_ll);

  fprintf(WHICH_OUTPUT, "      -face       loglikelihood: %2.4f  (%2.4f)\n", tp.face_ll, tp.face_ll / 2.0);
  fprintf(WHICH_OUTPUT, "      -vertex     loglikelihood: %2.4f  (%2.4f)\n", tp.vertex_ll, tp.vertex_ll / 2.0);
  fprintf(WHICH_OUTPUT, "      -normal dot loglikelihood: %2.4f  (%2.4f)\n", tp.curv_ll, tp.curv_ll);
  fprintf(WHICH_OUTPUT, "      -quad curv  loglikelihood: %2.4f  (%2.4f)\n", tp.qcurv_ll, tp.qcurv_ll / 2.0);
  fprintf(WHICH_OUTPUT, "      Total Loglikelihood : %2.4f\n", total_ll);

  MRISrestoreVertexPositions(mris, TMP_VERTICES);
  MRIScomputeMetricProperties(mris);

  /* free arrays */
  TPfree(&tp);
}


/*!
  \fn int MRISdefects2Seg(MRIS *surf, MRI *defects, int offset, MRI *vol)
  \brief Sample the defect numbers into the volume to create a
  segmentation. Works by going through all the vertices and finding
  the ones with non-zero defectno.  For each neighboring face, the
  voxel above and below is set to defectno+offset. vol must already
  exist and have the geometry of orig.mgz; it should be ready (eg,
  zeroed) to be filled with the defectno. surf should be the
  ?h.orig.nofix.  defects should be an MRI struct with each
  voxel/vertex indicating the defectno (eg, ?h.defect_labels).  A
  segmentation color table is imbedded in the vol.  It might be nice
  to be able to dilate the defects.
 */
int MRISdefects2Seg(MRIS *surf, MRI *defects, int offset, MRI *vol)
{
  int n, defectno,defectnomax;
  VERTEX *v, *vf;
  double delta, projsign;
  double x,y,z, c,r,s, cx,cy,cz;
  int ic,ir,is,nthface,faceno,vno,k,m,oob;
  float snorm[3];

  delta = vol->xsize/5.0; // could be smarter

  MRIS_SurfRAS2VoxelMap* sras2v_map = MRIS_makeRAS2VoxelMap(vol, surf);

  if(surf->nvertices != defects->width){
    printf("ERROR: MRISdefectNo2Vol(): dimension mismatch surf=%d, defects=%d\n",
	   surf->nvertices,defects->width);
    fflush(stdout);
    return(1);
  }

  // Get the maximum number of defects
  defectnomax = 0;
  for(n=0; n < surf->nvertices; n++){
    defectno = MRIgetVoxVal(defects,n,0,0,0);
    if(defectno > defectnomax) defectnomax = defectno;
  }
  int nentries = defectnomax+offset+2;
  if(vol->ct != NULL) nentries = MAX(nentries,vol->ct->nentries);
  COLOR_TABLE *ctab = CTABalloc(nentries);
  CTABunique(ctab, 10);
  for(n=0; n <= defectnomax; n++){
    sprintf(ctab->entries[n+offset]->name,"Defect-%03d",n+offset);
  }
  if(vol->ct == NULL) {
    vol->ct = ctab;
  }
  else {
    CTABmerge(ctab,vol->ct);
    CTABfree(&vol->ct);
    vol->ct=ctab;
  }

  for(n=0; n < surf->nvertices; n++){
    defectno = MRIgetVoxVal(defects,n,0,0,0);
    if(defectno==0) continue; // This vertex is not part of a defect
    v = &(surf->vertices[n]);
    // Go through all the faces associated with this vertex
    VERTEX_TOPOLOGY *vt = &(surf->vertices_topology[n]);
    for(nthface=0; nthface <  vt->num; nthface++){
      faceno = vt->f[nthface];
      // Compute the centroid of the face
      FACE *face = &(surf->faces[faceno]);
      cx = 0;
      cy = 0;
      cz = 0;
      for(k=0; k<3; k++){
	m = k + 1;
	if(m>2) m = 0;
	vno = face->v[k];
	vf = &(surf->vertices[vno]);
	cx += vf->x;
	cy += vf->y;
	cz += vf->z;
      }
      cx /= 3.0;
      cy /= 3.0;
      cz /= 3.0;
      // Compute the normal of the face
      mrisNormalFace(surf, faceno, 0, snorm);
      // Sample the volume at a voxel just above and just below the face
      for(projsign = -1; projsign <= +1; projsign+=2){
	x = cx + (projsign*delta)*snorm[0];
	y = cy + (projsign*delta)*snorm[1];
	z = cz + (projsign*delta)*snorm[2];
	MRIS_useRAS2VoxelMap(sras2v_map, vol, x, y, z, &c, &r, &s);
	ic = nint(c);
	ir = nint(r);
	is = nint(s);
	oob = MRIindexNotInVolume(vol,ic,ir,is);
	if(oob) continue;
	// Assign the voxel to the given vertex no
	MRIsetVoxVal(vol,ic,ir,is,0,defectno+offset);
      }
    }
  }
  MRIS_freeRAS2VoxelMap(&sras2v_map);

  return(0);
}

