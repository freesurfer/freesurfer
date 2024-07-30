#pragma once
/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
// The positions of the vertices and the immediate consequences of that 
//
#include "mrisurf_topology.h"
#include "topo_parms.h"
#include "realm.h"

int  mrisCheckSurface( MRIS       * mris);
bool mrisCheckDist    (MRIS const * mris);
bool mrisCheckDistOrig(MRIS const * mris);

void MRISmakeDist(MRIS *mris, int vno);




typedef struct PerThreadMRIDistance {
  MRI const * mri_distance;
  int heightTimesDepth,depth;
  float* elts;
} PerThreadMRIDistance;
PerThreadMRIDistance* makePerThreadMRIDistance(MRI const * const mri_distance);
void freePerThreadMRIDistance(PerThreadMRIDistance** ptdp);
float* perThreadMRIDistanceElt(PerThreadMRIDistance* ptd, int i, int j, int k);
void updateDistanceElt(volatile float* f, float distance, bool lockNeeded);
void updateDistanceEltFromSignArgAndSquareLockNeeded(volatile float* f, float distanceSignArg, float distanceSquared);
void updateDistanceEltFromSignArgAndSquareNoLockNeeded(
    volatile float* f, 
    float distanceSignArg, 
    float distanceSquared,
    float sqrtfDistanceSquared);

int mrisFindAllOverlappingFaces(MRIS *mris, MHT *mht, int fno, int *flist);
    
int mrisCalculateFaceCentroid(MRIS *mris, int fno, float *px, float *py, float *pz);

#define NOT_PROCESSED_YET 1000
#define NPY               NOT_PROCESSED_YET

/* volume floats */
#define xVOL(mri, x) (mri->xsize * (x - mri->xstart))
#define yVOL(mri, y) (mri->ysize * (y - mri->ystart))
#define zVOL(mri, z) (mri->zsize * (z - mri->zstart))
/* volume integers */
#define iVOL(mri, x) ((int)(xVOL(mri, x) + 0.5))
#define jVOL(mri, y) ((int)(yVOL(mri, y) + 0.5))
#define kVOL(mri, z) ((int)(zVOL(mri, z) + 0.5))
/* surface floats */
#define xSURF(mri, x) (mri->xstart + (float)x / mri->xsize)
#define ySURF(mri, y) (mri->ystart + (float)y / mri->ysize)
#define zSURF(mri, z) (mri->zstart + (float)z / mri->zsize)


int mrisComputeTangentPlanes(MRIS *mris);

extern double l_mri;
extern double l_unmri;
extern double l_curv;
extern double l_qcurv;

extern double l_vol;
extern double l_surf;
extern double l_wm;

void computeVertexPseudoNormal(MRIS const *mris, int vno, float norm[3], int verbose);

void  setFaceNorm    (MRIS const * const mris, int fno, float nx, float ny, float nz);
void  setFaceOrigArea(MRIS const * const mris, int fno, float orig_area);
float getFaceOrigArea(MRIS const * const mris, int fno);

int mrisChooseFace(MRIS *mris, MHT *mht, VERTEX *v);

void computeDefectFaceNormal(MRIS const * const mris, int const fno);

void mrisurf_deferSetFaceNorms  (MRIS* mris);
void mrisurf_recomputeFaceNorms (MRIS* mris);
void mrisurf_undeferSetFaceNorms(MRIS* mris);

int mrisMarkIntersections(MRIS *mris, int FillHoles);
int MRISmarkEdge(MRIS *surf, MRI *mask, int metricid, double thresh, int FillHoles);


#define OUTSIDE_VERTEX 0
#define INSIDE_VERTEX 1               /* not yet used */
#define EDGE_VERTEX 2                 /* part of an edge */
#define TRIANGLE_VERTEX 3             /* part of a triangle */
#define DISCARDED_VERTEX 4            /* excluded from the current tessellation */
#define BORDER_VERTEX TRIANGLE_VERTEX /* part of a triangle! */
#define USED_VERTEX 5                 /* used in the final tessellation */

void mrisDumpFace(MRIS const * mris, int fno, FILE *fp);


int    MRIScomputeAllDistances           (MRIS *mris);
void   MRIScomputeAvgInterVertexDist     (MRIS *Surf, double *StdDev);
void   mrisSetAvgInterVertexDist         (MRIS *Surf, double to);
int    mrisTrackTotalDistance            (MRIS *mris);
int    mrisTrackTotalDistanceNew         (MRIS *mris);

float  mrisComputeArea                   (MRIS *mris, int fac, int n);
float  MRIScomputeOrigArea               (MRIS* mris);
void   MRISsetOrigArea                   (MRIS* mris);


FaceNormCacheEntry const * getFaceNorm(MRIS    const * const mris, int fno);
FaceNormCacheEntry const * getFaceNorm(MRIS_MP const * const mris, int fno);
void setFaceNorm(MRIS    const * const mris, int fno, float nx, float ny, float nz);
void setFaceNorm(MRIS_MP const * const mris, int fno, float nx, float ny, float nz);

int mrisComputeBoundaryNormals(MRIS *mris);

int mrisComputeCurvatureMinMax(MRIS *mris);

int mrisComputeNormalDotDistribution(MRIS *mris, HISTOGRAM *h_dot);

int mrisComputePrincipalCurvatureDistributions(MRIS *mris,
                                                      HISTOGRAM *h_k1,
                                                      HISTOGRAM *h_k2,
                                                      MRI *mri_k1_k2);

float mrisDefectFaceMRILogLikelihood(
    MRIS *mris, MRI *mri, TP *tp, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white);

float mrisDefectVertexMRILogLikelihood(
    MRIS *mris, MRI *mri, TP *tp, HISTOGRAM *h_white, HISTOGRAM *h_gray, HISTOGRAM *h_grad, MRI *mri_gray_white);

float  mrisSampleAshburnerTriangleEnergy    (MRIS * const mris, int const vno, INTEGRATION_PARMS * const parms, float cx, float cy, float cz);
float  mrisSampleMinimizationEnergy         (MRIS *mris, int const vno,     INTEGRATION_PARMS *parms, float cx, float cy, float cz);
float  mrisSampleParallelEnergyAtVertex     (MRIS *mris, int const vno,     INTEGRATION_PARMS *parms);
float  mrisSampleParallelEnergy             (MRIS *mris, int const vno,     INTEGRATION_PARMS *parms, float cx, float cy, float cz);
float  mrisSampleNormalEnergy               (MRIS *mris, int const vno,     INTEGRATION_PARMS *parms, float cx, float cy, float cz);
float  mrisSampleSpringEnergy               (MRIS *mris, int const vno, float cx, float cy, float cz, INTEGRATION_PARMS *parms);

int mrisComputeOrigNormal (MRIS *mris, int vno, float norm[]);
int mrisComputeWhiteNormal(MRIS *mris, int vno, float norm[]);
int mrisComputePialNormal (MRIS *mris, int vno, float norm[]);

int mrisFindUnambiguousFace(MRIS *mris, MHT *mht, VERTEX *v, int *pnfound);


typedef struct ComputeDefectContext {
    RealmTree* realmTree;
    MRIS*      mris_deferred_norms;
} ComputeDefectContext;

static void constructComputeDefectContext(ComputeDefectContext* computeDefectContext) {
    bzero(computeDefectContext, sizeof(*computeDefectContext));
}

#define VERTEX_EDGE(vec, v0, v1)           VECTOR_LOAD(vec, v1->x     - v0->x,     v1->y     - v0->y,     v1->z     - v0->z)
#define VERTEX_ORIG_EDGE(vec, v0, v1)      VECTOR_LOAD(vec, v1->origx - v0->origx, v1->origy - v0->origy, v1->origz - v0->origz)
#define VERTEX_CANONICAL_EDGE(vec, v0, v1) VECTOR_LOAD(vec, v1->cx    - v0->cx,    v1->cy    - v0->cy,    v1->cz    - v0->cz)


void MRIScomputeMetricProperties(MRIS_MP* mris);
void MRIScomputeMetricPropertiesFaster(MRIS *mris);
