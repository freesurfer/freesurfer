/*
 * @file utilities common to mrisurf*.c but not used outside them
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
#include "mrisurf_base.h"

// This is for backwards compatibility for when we don't want to fix the vertex area. 
// Default is to fix, but this can be changed by setting to 0.
//
int fix_vertex_area = 1;

int MRISsetFixVertexAreaValue(int value)
{
  fix_vertex_area = value;
  return (fix_vertex_area);
}

int MRISgetFixVertexAreaValue(void)
{
  return (fix_vertex_area);
}

int MRISsetInflatedFileName(char *inflated_name)
{
  char fname[STRLEN];

  mrisurf_surface_names[0] = inflated_name;
  sprintf(fname, "%s.H", inflated_name);
  curvature_names[0] = (char *)calloc(strlen(fname) + 1, sizeof(char));
  strcpy(curvature_names[0], fname);
  return (NO_ERROR);
}

int MRISsetSulcFileName(const char *sulc_name)
{
  curvature_names[1] = (char *)calloc(strlen(sulc_name) + 1, sizeof(char));
  strcpy(curvature_names[1], sulc_name);
  return (NO_ERROR);
}

int MRISsetOriginalFileName(char *orig_name)
{
  mrisurf_surface_names[1] = mrisurf_surface_names[2] = orig_name;
  return (NO_ERROR);
}



double NEG_AREA_K = 10.0; /* was 200 */



// An accelerator for a hot function
//
typedef struct ActiveRealmTree {
    RealmTree*          realmTree;
    MRIS const*         mris;
    GetXYZ_FunctionType getXYZ;
} ActiveRealmTree;



#ifdef HAVE_OPENMP
static volatile bool    activeRealmTreesLockInited = false;
static omp_lock_t       activeRealmTreesLock;
#endif

static void acquireActiveRealmTrees() {
#ifdef HAVE_OPENMP
    if (!activeRealmTreesLockInited) {	// avoid critical if possible
        #pragma omp critical
    	if (!activeRealmTreesLockInited) {
	    omp_init_lock(&activeRealmTreesLock);
	    activeRealmTreesLockInited = true;
    	}
    }
    omp_set_lock(&activeRealmTreesLock);
#endif
}

static void releaseActiveRealmTrees() {
#ifdef HAVE_OPENMP
    omp_unset_lock(&activeRealmTreesLock);
#endif
}

static int              activeRealmTreesCapacity;
static ActiveRealmTree* activeRealmTrees;

void insertActiveRealmTree(MRIS const * const mris, RealmTree* realmTree, GetXYZ_FunctionType getXYZ)
{
    acquireActiveRealmTrees();
    if (mrisurf_activeRealmTreesSize == activeRealmTreesCapacity) {
        activeRealmTreesCapacity++;
        activeRealmTrees = 
            (ActiveRealmTree*)realloc(
                activeRealmTrees, activeRealmTreesCapacity*sizeof(ActiveRealmTree));
    }

    ActiveRealmTree* art = &activeRealmTrees[mrisurf_activeRealmTreesSize++];
    art->realmTree = realmTree;
    art->mris      = mris;
    art->getXYZ    = getXYZ;
    releaseActiveRealmTrees();
}

void removeActiveRealmTree(RealmTree* realmTree)
{
    acquireActiveRealmTrees();
    int i;
    for (i = 0; activeRealmTrees[i].realmTree != realmTree; i++) ;
    i++;
    for (; i != mrisurf_activeRealmTreesSize; i++) activeRealmTrees[i-1] = activeRealmTrees[i];
    mrisurf_activeRealmTreesSize--;
    releaseActiveRealmTrees();
}

int noteVnoMovedInActiveRealmTreesCount;
void noteVnoMovedInActiveRealmTrees(MRIS const * const mris, int vno) {
    acquireActiveRealmTrees();
    int i;
    noteVnoMovedInActiveRealmTreesCount++;
    for (i = 0; i < mrisurf_activeRealmTreesSize; i++) {
        if (activeRealmTrees[i].mris != mris) continue;
        if (0)
            fprintf(stderr,"Thread:%d updating realmTree:%p vno:%d\n", 
                omp_get_thread_num(), activeRealmTrees[i].realmTree, vno);
        noteIfXYZChangedRealmTree(
            activeRealmTrees[i].realmTree, 
            activeRealmTrees[i].mris,
            activeRealmTrees[i].getXYZ, 
            vno);
    }
    releaseActiveRealmTrees();
}

void notifyActiveRealmTreesChangedNFacesNVertices(MRIS const * const mris) {
    acquireActiveRealmTrees();
    int i;
    for (i = 0; i < mrisurf_activeRealmTreesSize; i++) {
        if (activeRealmTrees[i].mris != mris) continue;
        if (0)
            fprintf(stderr,"Thread:%d updating realmTree:%p\n", 
                omp_get_thread_num(), activeRealmTrees[i].realmTree);
        updateRealmTree(
            activeRealmTrees[i].realmTree, 
            activeRealmTrees[i].mris,
            activeRealmTrees[i].getXYZ);
    }
    releaseActiveRealmTrees();
}



/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/

static void MRISchangedNFacesNVertices(MRI_SURFACE * mris, bool scrambled) {
    // useful for debugging
}

bool MRISreallocVertices(MRI_SURFACE * mris, int max_vertices, int nvertices) {
  cheapAssert(nvertices >= 0);
  cheapAssert(max_vertices >= nvertices);

  mris->vertices = (VERTEX *)realloc(mris->vertices, max_vertices*sizeof(VERTEX));
  if (!mris->vertices) return false;
  #ifndef SEPARATE_VERTEX_TOPOLOGY
    mris->vertices_topology = mris->vertices;
  #else
    mris->vertices_topology = (VERTEX_TOPOLOGY *)realloc(mris->vertices_topology, max_vertices*sizeof(VERTEX_TOPOLOGY));
  #endif
  
  int change = max_vertices - mris->nvertices;
  if (change > 0) { // all above nvertices must be zero'ed, since MRISgrowNVertices does not...
      bzero(mris->vertices          + mris->nvertices, change*sizeof(VERTEX));
#ifdef SEPARATE_VERTEX_TOPOLOGY
      bzero(mris->vertices_topology + mris->nvertices, change*sizeof(VERTEX_TOPOLOGY));
#endif
  }

  *(int*)(&mris->max_vertices) = max_vertices;    // get around const
  *(int*)(&mris->nvertices)    = nvertices;       // get around const
    
  notifyActiveRealmTreesChangedNFacesNVertices(mris);
  
  return true;
}

void MRISgrowNVertices(MRI_SURFACE * mris, int nvertices) {
  if (nvertices > mris->max_vertices) {
    ErrorExit(ERROR_NOMEMORY, "MRISgrowNVertices: max vertices reached");
  }
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}

void MRIStruncateNVertices(MRI_SURFACE * mris, int nvertices) {
  cheapAssert(mris->nvertices >= nvertices);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}

void MRISremovedVertices(MRI_SURFACE * mris, int nvertices) {
  cheapAssert(mris->nvertices >= nvertices);
  MRISchangedNFacesNVertices(mris, true);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}

bool MRISreallocFaces(MRI_SURFACE * mris, int max_faces, int nfaces) {
  cheapAssert(nfaces >= 0);
  cheapAssert(max_faces >= nfaces);
  
  mris->faces  =
    (FACE *)realloc(mris->faces, max_faces*sizeof(FACE));
  if (!mris->faces) return false;

  mris->faceNormCacheEntries =
    (FaceNormCacheEntry*)realloc(mris->faceNormCacheEntries, max_faces*sizeof(FaceNormCacheEntry));
  if (!mris->faceNormCacheEntries) return false;
 
  mris->faceNormDeferredEntries =
    (FaceNormDeferredEntry*)realloc(mris->faceNormDeferredEntries, max_faces*sizeof(FaceNormDeferredEntry));
  if (!mris->faceNormDeferredEntries) return false;
  
  
  int change = max_faces - mris->nfaces;
  if (change > 0) { // all above nfaces must be zero'ed, since MRISgrowNFaces does not...
      bzero(mris->faces                   + mris->nfaces, change*sizeof(FACE));
      bzero(mris->faceNormCacheEntries    + mris->nfaces, change*sizeof(FaceNormCacheEntry));
      bzero(mris->faceNormDeferredEntries + mris->nfaces, change*sizeof(FaceNormDeferredEntry));
  }

  *(int*)(&mris->max_faces) = max_faces;    // get around const
  *(int*)(&mris->nfaces)    = nfaces;       // get around const
    
  return true;
}

bool MRISallocateFaces(MRI_SURFACE * mris, int nfaces) {
    return MRISreallocFaces(mris, nfaces, nfaces);
}

void MRISgrowNFaces(MRI_SURFACE * mris, int nfaces) {
  if (nfaces > mris->max_faces) {
    ErrorExit(ERROR_NOMEMORY, "mrisAddFace: max faces reached");
  }
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}

void MRIStruncateNFaces(MRI_SURFACE * mris, int nfaces) {
  cheapAssert(mris->nfaces >= nfaces);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}

void MRISremovedFaces(MRI_SURFACE * mris, int nfaces) {
  cheapAssert(mris->nfaces >= nfaces);
  MRISchangedNFacesNVertices(mris, true);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}


void MRISoverAllocVerticesAndFaces(MRI_SURFACE* mris, int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRISchangedNFacesNVertices(mris, false);
  if (nvertices < 0) ErrorExit(ERROR_BADPARM, "ERROR: MRISalloc: nvertices=%d < 0\n", nvertices);
  if (nfaces    < 0) ErrorExit(ERROR_BADPARM, "ERROR: MRISalloc: nfaces=%d < 0\n", nfaces);

  if (max_vertices <= nvertices) max_vertices = nvertices;
  if (max_faces    <= nfaces   ) max_faces    = nfaces;

  if (!MRISreallocVertices(mris, max_vertices, nvertices)) ErrorExit(ERROR_NO_MEMORY, 
    "MRISalloc(%d, %d): could not allocate vertices", max_vertices, sizeof(VERTEX));

  if (!MRISreallocFaces(mris, max_faces, nfaces)) ErrorExit(ERROR_NO_MEMORY, 
    "MRISalloc(%d, %d): could not allocate faces", nfaces,
    sizeof(FACE)+sizeof(FaceNormCacheEntry));
}


void MRISreallocVerticesAndFaces(MRI_SURFACE *mris, int nvertices, int nfaces) {
  MRISoverAllocVerticesAndFaces(mris, nvertices, nfaces, nvertices, nfaces);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISoverAlloc(int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRI_SURFACE* mris = (MRI_SURFACE *)calloc(1, sizeof(MRI_SURFACE));
  if (!mris) ErrorExit(ERROR_NO_MEMORY, 
                "MRISalloc(%d, %d): could not allocate mris structure", max_vertices, max_faces);

  mris->nsize = 1;      // only 1-connected neighbors initially

  mris->useRealRAS = 0; // just initialize
  mris->vg.valid = 0;   // mark as invalid

  MRISoverAllocVerticesAndFaces(mris, max_vertices, max_faces, nvertices, nfaces);

  return (mris);
}


MRI_SURFACE* MRISalloc(int nvertices, int nfaces)
{
  return MRISoverAlloc(nvertices, nfaces, nvertices, nfaces);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISfree(MRI_SURFACE **pmris)
{
  MRI_SURFACE *mris;
  int vno,e,k,faceno;

  mris = *pmris;
  *pmris = NULL;

  if (mris->dx2) {
    free(mris->dx2);
  }
  if (mris->dy2) {
    free(mris->dy2);
  }
  if (mris->dz2) {
    free(mris->dz2);
  }
  if (mris->labels) {
    free(mris->labels);
  }
  if (mris->ct) {
    CTABfree(&mris->ct);
  }
  for (vno = 0; vno < mris->nvertices; vno++) {
    if (mris->vertices_topology[vno].f) {
      free(mris->vertices_topology[vno].f);
      mris->vertices_topology[vno].f = NULL;
    }
    if (mris->vertices_topology[vno].n) {
      free(mris->vertices_topology[vno].n);
      mris->vertices_topology[vno].n = NULL;
    }
    if (mris->vertices[vno].dist) {
      free(mris->vertices[vno].dist);
      mris->vertices[vno].dist = NULL;
    }
    if (mris->vertices[vno].dist_orig) {
      free(mris->vertices[vno].dist_orig);
      mris->vertices[vno].dist_orig = NULL;
    }
    if (mris->vertices_topology[vno].v) {
      free(mris->vertices_topology[vno].v);
      mris->vertices_topology[vno].v = NULL;
    }
  }

  free(mris->faceNormDeferredEntries);
  free(mris->faceNormCacheEntries);
  
  if (mris->vertices) {
    free(mris->vertices);
  }
  if(mris->faces) {
    for(faceno = 0; faceno < mris->nfaces; faceno++){
      if(mris->faces[faceno].norm) DMatrixFree(&mris->faces[faceno].norm);
      for(k=0; k < 3; k++){
	if(mris->faces[faceno].gradNorm[k]) DMatrixFree(&mris->faces[faceno].gradNorm[k]);
      }
    }
    free(mris->faces);
  }
  if(mris->edges) {
    for(e=0; e < mris->nedges; e++){
      for(k=0; k<4; k++){
	if(mris->edges[e].gradDot[k]) DMatrixFree(&mris->edges[e].gradDot[k]);
      }
    }
    free(mris->edges);
  }
  if (mris->free_transform) {
    if (mris->SRASToTalSRAS_) {
      MatrixFree(&(mris->SRASToTalSRAS_));
    }
    if (mris->TalSRASToSRAS_) {
      MatrixFree(&(mris->TalSRASToSRAS_));
    }
    if (mris->lta) {
      LTAfree(&(mris->lta));
    }
  }

  {
    int i;
    for (i = 0; i < mris->ncmds; i++) {
      free(mris->cmdlines[i]);
    }
  }
  if (mris->m_sras2vox) {
    MatrixFree(&mris->m_sras2vox);
  }

  free(mris);
  return (NO_ERROR);
}


int MRISfreeDists(MRI_SURFACE *mris)
{
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    if (mris->vertices[vno].dist) {
      free(mris->vertices[vno].dist);
    }
    if (mris->vertices[vno].dist_orig) {
      free(mris->vertices[vno].dist_orig);
    }
    mris->vertices         [vno].dist = mris->vertices[vno].dist_orig = NULL;
    mris->vertices_topology[vno].vtotal = 0;
  }

  return (NO_ERROR);
}


char *mrisurf_surface_names[3] = {"inflated", "smoothwm", "smoothwm"};
char *curvature_names[3] = {"inflated.H", "sulc", NULL};
int MRISsetCurvatureName(int nth, char *name)
{
  if (nth > 2) {
    printf("ERROR: MRISsetCurvatureName() nth=%d > 2\n", nth);
    return (1);
  }
  curvature_names[nth] = strcpyalloc(name);
  return (0);
}
int MRISprintCurvatureNames(FILE *fp)
{
  int k;
  for (k = 0; k < sizeof(curvature_names) / sizeof(curvature_names[0]); k++) {
    if (curvature_names[k])
      printf("%d %s\n", k, curvature_names[k]);
    else if (mrisurf_surface_names[k])
      printf("%d %s (computed)\n", k, mrisurf_surface_names[k]);
  }
  return (0);
}


static const float sigmas_default[] = {
    //    16.00, 4.0f, 2.0f, 1.0f, 0.5f
    //    8.00f, 4.00f, 2.0f, 0.5f
    4.00f,
    2.0f,
    1.0f,
    0.5f};
const float * sigmas = sigmas_default;      // can be changed by caller
double nsigmas = (sizeof(sigmas_default) / sizeof(sigmas_default[0]));

int MRISsetRegistrationSigmas(float *sigmas_local, int nsigmas_local)
{
  if (sigmas_local == NULL) {
    nsigmas = (sizeof(sigmas_default) / sizeof(sigmas_default[0]));
    sigmas  = sigmas_default;
  }
  else {
    nsigmas = nsigmas_local;
    sigmas  = sigmas_local;
  }
  return (NO_ERROR);
}


VOXEL_LIST **vlst_alloc(MRI_SURFACE *mris, int max_vox)
{
  int vno;
  VOXEL_LIST **vl;

  vl = (VOXEL_LIST **)calloc(mris->nvertices, sizeof(VOXEL_LIST *));
  if (vl == NULL) ErrorExit(ERROR_NOMEMORY, "vlst_alloc(%d): could not allocate vl array", max_vox);
  for (vno = 0; vno < mris->nvertices; vno++) {
    vl[vno] = VLSTalloc(max_vox);
    if (vl[vno] == NULL) ErrorExit(ERROR_NOMEMORY, "vlst_alloc(%d): could not allocate list for vno %d", max_vox, vno);
    vl[vno]->nvox = 0;
  }
  return (vl);
}

int vlst_free(MRI_SURFACE *mris, VOXEL_LIST ***pvl)
{
  VOXEL_LIST **vl = *pvl;
  int vno;

  *pvl = NULL;
  for (vno = 0; vno < mris->nvertices; vno++) VLSTfree(&vl[vno]);
  free(vl);
  return (NO_ERROR);
}

int vlst_enough_data(MRI_SURFACE *mris, int vno, VOXEL_LIST *vl, double displacement)
{
  double dot, dx, dy, dz;
  int i, nin, nout;
  VERTEX *v;
  double xs, ys, zs;

  v = &mris->vertices[vno];
  xs = v->x + displacement * v->nx;
  ys = v->y + displacement * v->ny;
  zs = v->z + displacement * v->nz;
  for (nin = nout = i = 0; i < vl->nvox; i++) {
    dx = vl->xd[i] - xs;
    dy = vl->yd[i] - ys;
    dz = vl->zd[i] - zs;
    dot = dx * v->nx + dy * v->ny + dz * v->nz;
    if (dot > 0)  // outside surface
      nout++;
    else  // inside the surface
      nin++;
  }
  return ((nout >= 2 && nin >= 2) || (nout >= 10 && nin >= 1) || (nout + nin > 4));
}

int vlst_add_to_list(VOXEL_LIST *vl_src, VOXEL_LIST *vl_dst)
{
  int i;

  for (i = 0; i < vl_src->nvox; i++)
    VLSTaddUnique(vl_dst, vl_src->xi[i], vl_src->yi[i], vl_src->zi[i], vl_src->xd[i], vl_src->yd[i], vl_src->zd[i]);

  return (NO_ERROR);
}



/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static int edgesIntersectStable(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2);

#if SPHERE_INTERSECTION

/* check for intersection on the sphere */
int edgesIntersect(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2)
{
  VERTEX *v1, *v2;
  double n0[3], n1[3], n2[3], u0[3], u1[3], u2[3], u3[3], u[3];
  double a0, a1, a2, a3, a, b;

  // first test if some vertices are the same
  if (edge1->vno1 == edge2->vno1 || edge1->vno1 == edge2->vno2 || edge1->vno2 == edge2->vno1 ||
      edge1->vno2 == edge2->vno2) {
    return (0);
  }

  // INTERSECTION ONTO THE SPHERE

  // compute normals
  v1 = &mris->vertices[edge1->vno1];
  v2 = &mris->vertices[edge1->vno2];
  u0[0] = v1->cx;
  u0[1] = v1->cy;
  u0[2] = v1->cz;
  u1[0] = v2->cx;
  u1[1] = v2->cy;
  u1[2] = v2->cz;
  F_CROSS(u0, u1, n0);

  v1 = &mris->vertices[edge2->vno1];
  v2 = &mris->vertices[edge2->vno2];
  u2[0] = v1->cx;
  u2[1] = v1->cy;
  u2[2] = v1->cz;
  u3[0] = v2->cx;
  u3[1] = v2->cy;
  u3[2] = v2->cz;
  F_CROSS(u2, u3, n1);

  a0 = F_DOT(u0, n1);
  a1 = F_DOT(u1, n1);
  a2 = F_DOT(u2, n0);
  a3 = F_DOT(u3, n0);

// if a0,a1, a2, or a3 are very small,
// should definitely use edgesIntersectStable instead
#define SMALL_VALUES 0.0000001
  if (fabs(a0) < SMALL_VALUES || fabs(a1) < SMALL_VALUES || fabs(a2) < SMALL_VALUES || fabs(a3) < SMALL_VALUES) {
    return edgesIntersectStable(mris, edge1, edge2);
  }

  a = a0 * a1;
  b = a2 * a3;

  if (a > 0) {
    return (0);
  }
  else if (a < 0) {
    //      fprintf(stderr,"-");
    if (b > 0) {
      return (0);
    }
    else {
      return (1);
    }
  }
  else {
    //    fprintf(stderr,"+");
    if (b > 0) {
      return (0);
    }
    else if (b < 0) {
      return (1);
    }
    else {
      // special case again! The points are exactly aligned...
      double c0, c1, c2, c3, x1_min, x1_max, x2_min, x2_max;

#if 0
      fprintf(stderr,"^");
      {

        int r1,r2;
        r1=sphereEdgesIntersect(mris,edge1,edge2);
        r2=newEdgesIntersect3(mris,edge1,edge2);

        fprintf(stderr,"%d-%d",r1,r2);

        if (r1!=r2)
        {
          debugEdgesIntersect(mris,edge1,edge2);
        }

      }
#endif

      u[0] = (u0[0] + u1[0] + u2[0] + u3[0]);
      u[1] = (u0[1] + u1[1] + u2[1] + u3[1]);
      u[2] = (u0[2] + u1[2] + u2[2] + u3[2]);

      F_CROSS(u, n0, n2);

      c0 = F_DOT(u0, n2);
      c1 = F_DOT(u1, n2);
      c2 = F_DOT(u2, n2);
      c3 = F_DOT(u3, n2);

      x1_min = MIN(c0, c1);
      x1_max = MAX(c0, c1);
      x2_min = MIN(c2, c3);
      x2_max = MAX(c2, c3);

      return (x2_max >= x1_min && x2_min <= x1_max);
    }
  }
  return (0);
}

/* it should be a more stable implementation but a bit
   slower than the previous one */
static int edgesIntersectStable(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2)
{
  VERTEX *v1, *v2;
  double n0[3], n1[3], n2[3], u0[3], u1[3], u2[3], u3[3], u[3];
  double og0[3], og1[3], v_0[3], v_1[3];
  double a0, a1, a2, a3, a, b;

  // first test if some vertices are the same
  if (edge1->vno1 == edge2->vno1 || edge1->vno1 == edge2->vno2 || edge1->vno2 == edge2->vno1 ||
      edge1->vno2 == edge2->vno2) {
    return (0);
  }

  // INTERSECTION ONTO THE SPHERE
  v1 = &mris->vertices[edge1->vno1];
  v2 = &mris->vertices[edge1->vno2];
  og0[0] = (v1->cx + v2->cx) / 2.0;
  og0[1] = (v1->cy + v2->cy) / 2.0;
  og0[2] = (v1->cz + v2->cz) / 2.0;

  // compute tangent vector
  v_0[0] = (v2->cx - v1->cx);
  v_0[1] = (v2->cy - v1->cy);
  v_0[2] = (v2->cz - v1->cz);

  // compute normal
  F_CROSS(og0, v_0, n0);

  v1 = &mris->vertices[edge2->vno1];
  v2 = &mris->vertices[edge2->vno2];
  og1[0] = (v1->cx + v2->cx) / 2.0;
  og1[1] = (v1->cy + v2->cy) / 2.0;
  og1[2] = (v1->cz + v2->cz) / 2.0;

  // compute tangent vector
  v_1[0] = (v2->cx - v1->cx);
  v_1[1] = (v2->cy - v1->cy);
  v_1[2] = (v2->cz - v1->cz);

  // compute normal
  F_CROSS(og1, v_1, n1);

  // compute vectors
  v1 = &mris->vertices[edge1->vno1];
  v2 = &mris->vertices[edge1->vno2];
  u0[0] = v1->cx - og1[0];
  u0[1] = v1->cy - og1[1];
  u0[2] = v1->cz - og1[2];
  u1[0] = v2->cx - og1[0];
  u1[1] = v2->cy - og1[1];
  u1[2] = v2->cz - og1[2];
  v1 = &mris->vertices[edge2->vno1];
  v2 = &mris->vertices[edge2->vno2];
  u2[0] = v1->cx - og0[0];
  u2[1] = v1->cy - og0[1];
  u2[2] = v1->cz - og0[2];
  u3[0] = v2->cx - og0[0];
  u3[1] = v2->cy - og0[1];
  u3[2] = v2->cz - og0[2];

  a0 = F_DOT(u0, n1);
  a1 = F_DOT(u1, n1);
  a2 = F_DOT(u2, n0);
  a3 = F_DOT(u3, n0);

  a = a0 * a1;
  b = a2 * a3;

  if (a > 0) {
    return (0);
  }
  else if (a < 0) {
    if (b > 0) {
      return (0);
    }
    else {
      return (1);
    }
  }
  else {
    //    fprintf(stderr,"+");
    if (b > 0) {
      return (0);
    }
    else if (b < 0) {
      return (1);
    }
    else {
      // special case again! The points are exactly aligned...
      double c0, c1, c2, c3, x1_min, x1_max, x2_min, x2_max;

      // fprintf(stderr,"-");

      u[0] = (u0[0] + u1[0] + u2[0] + u3[0]);
      u[1] = (u0[1] + u1[1] + u2[1] + u3[1]);
      u[2] = (u0[2] + u1[2] + u2[2] + u3[2]);

      F_CROSS(u, n0, n2);

      c0 = F_DOT(u0, n2);
      c1 = F_DOT(u1, n2);
      c2 = F_DOT(u2, n2);
      c3 = F_DOT(u3, n2);

      x1_min = MIN(c0, c1);
      x1_max = MAX(c0, c1);
      x2_min = MIN(c2, c3);
      x2_max = MAX(c2, c3);

      return (x2_max >= x1_min && x2_min <= x1_max);
    }
  }
  return (0);
}

#else
/* original version */
int edgesIntersect(MRI_SURFACE *mris, EDGE *edge1, EDGE *edge2)
{
  VERTEX *v1, *v2;
  double b1, b2, m1, m2, x1_start, x1_end, y1_start, y1_end, x2_start, x2_end, y2_start, y2_end, x, y, x1min, x1max,
      y1min, y1max, x2min, x2max, y2min, y2max, cx, cy, cz;
  double origin[3], e0[3], e1[3], e0_0, e0_1, e0_2, e1_0, e1_1, e1_2;

  if (edge1->vno1 == edge2->vno1 || edge1->vno1 == edge2->vno2 || edge1->vno2 == edge2->vno1 ||
      edge1->vno2 == edge2->vno2) {
    return (0);
  }

  mrisComputeCanonicalEdgeBasis(mris, edge1, edge2, origin, e0, e1);
  e0_0 = e0[0];
  e0_1 = e0[1];
  e0_2 = e0[2];
  e1_0 = e1[0];
  e1_1 = e1[1];
  e1_2 = e1[2];
  if (edge1->vno1 == Gdiag_no || edge2->vno1 == Gdiag_no || edge1->vno2 == Gdiag_no || edge2->vno2 == Gdiag_no) {
    DiagBreak();
  }
  v1 = &mris->vertices[edge1->vno1];
  v2 = &mris->vertices[edge1->vno2];
  cx = v1->cx - origin[0];
  cy = v1->cy - origin[1];
  cz = v1->cz - origin[2];
  x1_start = cx * e0_0 + cy * e0_1 + cz * e0_2;
  y1_start = cx * e1_0 + cy * e1_1 + cz * e1_2;
  cx = v2->cx - origin[0];
  cy = v2->cy - origin[1];
  cz = v2->cz - origin[2];
  x1_end = cx * e0_0 + cy * e0_1 + cz * e0_2;
  y1_end = cx * e1_0 + cy * e1_1 + cz * e1_2;
  x1min = MIN(x1_start, x1_end);
  x1max = MAX(x1_start, x1_end);
  y1min = MIN(y1_start, y1_end);
  y1max = MAX(y1_start, y1_end);
  if (!FEQUAL(x1_start, x1_end)) {
    m1 = (y1_end - y1_start) / (x1_end - x1_start);
    b1 = y1_end - m1 * x1_end;
  }
  else {
    m1 = b1 = 0; /* will be handled differently */
  }

  v1 = &mris->vertices[edge2->vno1];
  v2 = &mris->vertices[edge2->vno2];
  cx = v1->cx - origin[0];
  cy = v1->cy - origin[1];
  cz = v1->cz - origin[2];
  x2_start = cx * e0_0 + cy * e0_1 + cz * e0_2;
  y2_start = cx * e1_0 + cy * e1_1 + cz * e1_2;
  cx = v2->cx - origin[0];
  cy = v2->cy - origin[1];
  cz = v2->cz - origin[2];
  x2_end = cx * e0_0 + cy * e0_1 + cz * e0_2;
  y2_end = cx * e1_0 + cy * e1_1 + cz * e1_2;
  x2min = MIN(x2_start, x2_end);
  x2max = MAX(x2_start, x2_end);
  y2min = MIN(y2_start, y2_end);
  y2max = MAX(y2_start, y2_end);
#if 0
#define INT_EPSILON 0.000
  x1max += INT_EPSILON ;
  y1max += INT_EPSILON ;
  x2max += INT_EPSILON ;
  y2max += INT_EPSILON ;
  x1min -= INT_EPSILON ;
  y1min -= INT_EPSILON ;
  x2min -= INT_EPSILON ;
  y2min -= INT_EPSILON ;
#endif

#if 1
  /* non-overlapping intervals can't intersect */
  if (x1min > x2max || x1max < x2min || y1min > y2max || y1max < y2min) {
    return (0);
  }
#endif

  /* handle special cases first */
  if (FEQUAL(x1_start, x1_end)) /* l1 is vertical */
  {
    if (FEQUAL(x2_start, x2_end)) /* both vertical */
    {
      return (FEQUAL(x1_start, x2_start)); /* same line */
    }
    m2 = (y2_end - y2_start) / (x2_end - x2_start);
    b2 = y2_end - m2 * x2_end;
    y = m2 * x1_start + b2; /* y coordinate of intersection */
    if (y >= y1min && y <= y1max) {
      return (1);
    }
  }
  else if (FEQUAL(x2_start, x2_end)) /* l2 is vertical */
  {
    if (FEQUAL(x1_start, x1_end)) /* both vertical */
    {
      return (FEQUAL(x1_start, x2_start)); /* same line */
    }
    y = m1 * x2_start + b1; /* y coord of intersection */
    if (y >= y2min && y <= y2max) {
      return (1);
    }
  }
  else /* neither line is vertical, compute intersection point */
  {
    m2 = (y2_end - y2_start) / (x2_end - x2_start);
    b2 = y2_end - m2 * x2_end;
    if (FEQUAL(m1, m2)) /* parallel lines */
    {
      if (FEQUAL(b1, b2)) /* same line, see if segments overlap */
      {
        return (x2max >= x1min && x2min <= x1max);
      }
      return (0);
    }
    x = (b2 - b1) / (m1 - m2);        /* intersection point */
    if ((x <= x1max && x >= x1min) && /* is it in l1 interval? */
        (x <= x2max && x >= x2min))   /* is it in l2 interval? */
    {
      return (1); /* both true - intersection */
    }
  }
  return (0);
}
#endif


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the line in the normal direction passing
  through vertex v intersects any of the triangles at the
  given location.
  ------------------------------------------------------*/
int mrisDirectionTriangleIntersection(
    MRI_SURFACE *mris, float x0, float y0, float z0, float nx, float ny, float nz, MHT *mht, double *pdist, int vno)
{
  double dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3], dot;
  float x, y, z, dx, dy, dz;
  FACE *face;
  int i, found, fno, ret;
  static MHBT *prev_bucket = NULL;

  dist = *pdist;
  dir[0] = nx;
  dir[1] = ny;
  dir[2] = nz;
  pt[0] = x0;
  pt[1] = y0;
  pt[2] = z0;
  x = x0 + nx * dist;
  y = y0 + ny * dist;
  z = z0 + nz * dist;

  min_dist = 10000.0f;

  MHBT *bucket = MHTacqBucket(mht, x, y, z);
  if (bucket == NULL) {
    return (0);
  }

  prev_bucket = bucket;

  MHB *bin;
  for (bin = bucket->bins, found = i = 0; i < bucket->nused; i++, bin++) {
    fno = bin->fno;
    face = &mris->faces[fno];
    if (face->v[0] == vno || face->v[1] == vno || face->v[2] == vno) continue;
    if (fno == 1287 || fno == 5038) {
      DiagBreak();
    }

    load_triangle_vertices(mris, fno, U0, U1, U2, CURRENT_VERTICES);
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt);
    if (ret) {
      dx = int_pt[0] - x0;
      dy = int_pt[1] - y0;
      dz = int_pt[2] - z0;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      dot = dx * nx + dy * ny + dz * nz;
      if (dot >= 0 && dist < min_dist) {
        found = 1;
        *pdist = min_dist = dist;
      }
    }
  }
  
  MHTrelBucket(&bucket);
  return (found);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  See if the line in the normal direction passing
  through vertex v intersects any of the triangles at the
  given location.
  ------------------------------------------------------*/
int mrisAllNormalDirectionCurrentTriangleIntersections(
    MRI_SURFACE *mris, VERTEX *v, MHT *mht, double *pdist, int *flist)
{
  double dist, min_dist, U0[3], U1[3], U2[3], pt[3], dir[3], int_pt[3];
  float nx, ny, nz, x, y, z, dx, dy, dz, dot;
  int i, found, fno, ret;
  static MHBT *last_bucket = NULL;
  static VERTEX *last_v = NULL;

  dist = *pdist;
  nx = v->nx;
  ny = v->ny;
  nz = v->nz;
  dir[0] = v->nx;
  dir[1] = v->ny;
  dir[2] = v->nz;
  pt[0] = v->x;
  pt[1] = v->y;
  pt[2] = v->z;
  x = v->x + nx * dist;
  y = v->y + ny * dist;
  z = v->z + nz * dist;

  MHBT *bucket = MHTacqBucket(mht, x, y, z);
  if (bucket == NULL) {
    return (-1);
  }

  if (last_v == v && bucket == last_bucket) {
    MHTrelBucket(&bucket);
    return (-2);
  }

  last_v = v;
  last_bucket = bucket;

  min_dist = 10000.0f;
  MHB *bin;
  for (bin = bucket->bins, found = i = 0; i < bucket->nused; i++, bin++) {
    fno = bin->fno;

    load_triangle_vertices(mris, fno, U0, U1, U2, CURRENT_VERTICES);
    ret = triangle_ray_intersect(pt, dir, U0, U1, U2, int_pt);
    if (ret) {
      dx = int_pt[0] - v->x;
      dy = int_pt[1] - v->y;
      dz = int_pt[2] - v->z;
      dot = dx * nx + dy * ny + dz * nz;
      if (dot < 0) /* in direciton antiparallel to normal direction */
      {
        continue;
      }
      dist = sqrt(dot);
      if (flist) {
        flist[found] = fno;
      }
      found++;
      if (found == 1000)
        ErrorExit(ERROR_BADPARM,
                  "too many intersecting faces.  "
                  "check filled volume for correctness\n");
      if (dist < min_dist) {
        *pdist = min_dist = dist;
      }
    }
  }

  MHTrelBucket(&bucket);

  return (found);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISsampleAtEachDistance(MRI_SURFACE *mris, int nbhd_size, int nbrs_per_distance)
{
  int n, nbrs_array[MAX_NBHD_SIZE];

  if (!nbhd_size) {
    return (NO_ERROR);
  }

  if (Gdiag & (DIAG_HEARTBEAT | DIAG_SHOW)) {
    fprintf(stdout, "sampling long-range distances");
  }
  for (n = 0; n <= nbhd_size; n++) {
    nbrs_array[n] = nbrs_per_distance;
  }
  return MRISsampleDistances(mris, nbrs_array, nbhd_size);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Expand the list of neighbors of each vertex, reallocating
  the v->v array to hold the expanded list.
  ------------------------------------------------------*/
#define MAX_V 5000                                         /* max for any one node, actually way too big */
#define TRIANGLE_DISTANCE_CORRECTION 1.09f                 /*1.1f*/
/*1.066f*/ /*1.12578*/ /* 1.13105f*/                       /*1.1501f  (1.1364f)*/
#define QUADRANGLE_DISTANCE_CORRECTION ((1 + sqrt(2)) / 2) /* 1.2071  */
int MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd)
{
  int i, n, vno, vnum, old_vnum, total_nbrs, max_possible, max_v, vtotal;
  int *vnbrs, *vall, *vnb, found, n2, vnbrs_num, vall_num, nbhd_size, done, checks = 0;
  float xd, yd, zd, min_dist, dist, dist_scale, *old_dist;
  float min_angle, angle;
  VECTOR *v1, *v2;
  float c[100], *orig_old_dist;
  int nc[100], *old_v;
  int diag_vno1, diag_vno2;
  char *cp;

  orig_old_dist = old_dist = (float *)calloc(MAX_V, sizeof(float));
  old_v = (int *)calloc(MAX_V, sizeof(int));
  if (old_dist == NULL || old_v == NULL)
    ErrorExit(ERROR_NOMEMORY, "MRISsampleDistances: could not allocated old_dist buffer");

  if ((cp = getenv("VDIAG1")) != NULL) {
    diag_vno1 = atoi(cp);
  }
  else {
    diag_vno1 = -1;
  }
  if ((cp = getenv("VDIAG2")) != NULL) {
    diag_vno2 = atoi(cp);
  }
  else {
    diag_vno2 = -1;
  }
  if (diag_vno1 >= 0) {
    printf("\nlooking for vertex pair %d, %d\n", diag_vno1, diag_vno2);
  }

  memset(c, 0, 100 * sizeof(float));
  memset(nc, 0, 100 * sizeof(int));
  v1 = VectorAlloc(3, MATRIX_REAL);
  v2 = VectorAlloc(3, MATRIX_REAL);

  /* adjust for Manhattan distance */
  if (IS_QUADRANGULAR(mris)) {
    dist_scale = (1.0 + sqrt(2.0)) / 2.0f;
  }
  else {
    dist_scale = TRIANGLE_DISTANCE_CORRECTION;
  }

  vnbrs = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vall = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vnb = (int *)calloc(MAX_NBHD_VERTICES, sizeof(int));
  vtotal = total_nbrs = 0;
  for (vtotal = max_possible = 0, n = 1; n <= max_nbhd; n++) {
    max_possible += nbrs[n];
    if (n > mris->nsize) {
      vtotal += nbrs[n];
    }
  }
  if (Gdiag & DIAG_HEARTBEAT)
    fprintf(stdout,
            "\nsampling %d dists/vertex (%2.1f at each dist) = %2.1fMB\n",
            vtotal,
            (float)vtotal / ((float)max_nbhd - (float)mris->nsize),
            (float)vtotal * MRISvalidVertices(mris) * sizeof(float) * 3.0f / (1024.0f * 1024.0f));

  for (vno = 0; vno < mris->nvertices; vno++) {
  
    if ((Gdiag & DIAG_HEARTBEAT) && (!(vno % (mris->nvertices / 10))))
      fprintf(stdout, "%%%1.0f done\n", 100.0f * (float)vno / (float)mris->nvertices);
    if ((vno > 139000 || (!(vno % 100))) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
      if (mrisCheckSurfaceNbrs(mris) != 1) {
        printf("surface bad at vertex %d!\n", vno);
        DiagBreak();
      }
    }
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
    VERTEX          * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    if (v->ripflag) {
      continue;
    }

    /* small neighborhood is always fixed, don't overwrite them */
    vtotal = vt->vtotal;
    if (vt->nsize == 3) {
      vt->vtotal = vt->v3num;
    }
    else if (vt->nsize == 2) {
      vt->vtotal = vt->v2num;
    }
    else {
      vt->vtotal = vt->vnum;
    }

    max_v = vt->vtotal + max_possible;
    if (vtotal < max_v) /* won't fit in current allocation,
                         reallocate stuff */
    {
      /* save and restore neighbor list */
      memmove(old_v, vt->v, vt->vtotal * sizeof(vt->v[0]));
      free(vt->v);
      vt->v = (int *)calloc(max_v, sizeof(int));
      if (!vt->v)
        ErrorExit(ERROR_NO_MEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "nbrs at v=%d",
                  max_v,
                  vno);
      memmove(vt->v, old_v, vt->vtotal * sizeof(vt->v[0]));

      /* save and restore distance vector */
      if (vt->vtotal >= MAX_V || old_dist != orig_old_dist)
        printf("!!!!!!!!!!!!! v %d has too many (%d) vertices to save distances for !!!!!!!!!!!!!!!!!\n", vno, vtotal);
      memmove(old_dist, v->dist, vt->vtotal * sizeof(v->dist[0]));
      free(v->dist);
      v->dist = (float *)calloc(max_v, sizeof(float));
      if (!v->dist)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d",
                  max_v,
                  vno);
      memmove(v->dist, old_dist, vt->vtotal * sizeof(v->dist[0]));

      /* save and restore original distance vector */
      {
        int k;
        for (k = 0; k < vt->vtotal; k++) {
          old_dist[k] = 0;
          v->dist_orig[k] *= 1;
          old_dist[k] = v->dist_orig[k];
        }
        if (vno == Gdiag_no) DiagBreak();
      }
      memmove(old_dist, v->dist_orig, vt->vtotal * sizeof(v->dist_orig[0]));
      free(v->dist_orig);
      v->dist_orig = (float *)calloc(max_v, sizeof(float));
      if (!v->dist_orig)
        ErrorExit(ERROR_NOMEMORY,
                  "MRISsampleDistances: could not allocate list of %d "
                  "dists at v=%d",
                  max_v,
                  vno);
      memmove(v->dist_orig, old_dist, vt->vtotal * sizeof(v->dist_orig[0]));
    }

    if ((vno > 139000 || !(vno % 100)) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
      if (mrisCheckSurfaceNbrs(mris) != 1) {
        printf("surface bad at vertex %d!\n", vno);
        DiagBreak();
      }
    }
    /*
     find all the neighbors at each extent (i.e. 1-neighbors, then
     2-neighbors, etc..., marking their corrected edge-length distances
     as you go.
    */
    vall[0] = vno;
    vall_num = 1;
    old_vnum = 0;
    v->marked = 1; /* a hack - it is a zero neighbor */
    for (nbhd_size = 1; vall_num < MAX_NBHD_VERTICES && nbhd_size <= max_nbhd; nbhd_size++) {

      /* expand neighborhood outward by a ring of vertices */
      vnbrs_num = 0; /* will count neighbors in this ring */
      vnum = vall_num;
      for (found = 0, n = old_vnum; vall_num < MAX_NBHD_VERTICES && n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX          const * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }

        /* search through vn's neighbors to find an unmarked vertex */
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag || vn2->marked) {
            continue;
          }

          /* found one, mark it and put it in the vall list */
          found++;
          vn2->marked = nbhd_size;
          vnb[vnum] = vall[n];
          vall[vnum++] = vnt->v[n2];
          if (nbrs[nbhd_size] > 0) /* want to store this distance */
          {
            vnbrs[vnbrs_num++] = vnt->v[n2];
          }
        }
      } /* done with all neighbors at previous distance */

      /* found all neighbors at this extent - calculate distances */
      old_vnum = vall_num; /* old_vnum is index of 1st
                          nbr at this distance*/
      vall_num += found;   /* vall_num is total # of nbrs */
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX                * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }
        for (min_dist = UNFOUND_DIST, n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (vn2->ripflag) {
            continue;
          }
          if (!vn2->marked || vn2->marked == nbhd_size) {
            continue;
          }
          xd = vn2->x - vn->x;
          yd = vn2->y - vn->y;
          zd = vn2->z - vn->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1) {
            dist /= dist_scale;
          }
          if (vn2->d + dist < min_dist) {
            min_dist = vn2->d + dist;
          }
#else
          dist = (dist + vn2->d * distance_scale[vn2->marked]) / distance_scale[vn2->marked + 1];
          if (dist < min_dist) {
            min_dist = dist;
          }
#endif
        }
        vn->d = min_dist;
        if (nbhd_size <= 2) {
          xd = vn->x - v->x;
          yd = vn->y - v->y;
          zd = vn->z - v->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
          vn->d = dist;
        }
        if (vn->d >= UNFOUND_DIST / 2) {
          printf(
              "***** WARNING - surface distance not found at "
              "vno %d, vall[%d] = %d (vnb[%d] = %d ******",
              vno,
              n,
              vall[n],
              n,
              vnb[n]);
          mrisCheckSurfaceNbrs(mris);
          DiagBreak();
          exit(1);
        }
        if ((vall[n] == diag_vno1 && vno == diag_vno2) || (vall[n] == diag_vno2 && vno == diag_vno1))
          printf("vn %d is %2.3f mm from v %d at ring %d\n", vall[n], vn->d, vno, nbhd_size);
      }

      /*
       now check each to see if a neighbor at the same 'distance'
       is actually closer than neighbors which are 'nearer' (i.e. maybe
       the path through a 3-neighbor is shorter than that through any
       of the 2-neighbors.
      */
      for (n = old_vnum; n < vall_num; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vall[n]];
        VERTEX                * const vn  = &mris->vertices         [vall[n]];
        if (vn->ripflag) {
          continue;
        }
        min_dist = vn->d;
        for (n2 = 0; n2 < vnt->vnum; n2++) {
          VERTEX const * const vn2 = &mris->vertices[vnt->v[n2]];
          if (!vn2->marked || vn2->marked != nbhd_size || vn2->ripflag) {
            continue;
          }
          xd = vn2->x - vn->x;
          yd = vn2->y - vn->y;
          zd = vn2->z - vn->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
#if !MULTI_DIST_SCALING
          if (nbhd_size > 1) {
            dist /= dist_scale;
          }
          if (vn2->d + dist < min_dist) {
            min_dist = vn2->d + dist;
          }
#else
          dist = (dist + vn2->d * distance_scale[vn2->marked]) / distance_scale[vn2->marked + 1];
          if (dist < min_dist) {
            min_dist = dist;
          }
#endif
        }
        vn->d = min_dist;
        if ((min_dist > 5 * nbhd_size) || min_dist > 60) {
          DiagBreak();
        }
        {
          xd = vn->x - v->x;
          yd = vn->y - v->y;
          zd = vn->z - v->z;
          dist = sqrt(xd * xd + yd * yd + zd * zd);
          c[nbhd_size] += vn->d / dist;
          nc[nbhd_size]++;
        }
      }

      /* if this set of neighbors are to be stored, sample from them */
      if (nbrs[nbhd_size] <= 0) {
        continue;
      }

      /* make sure the points are not too close together */
      min_angle = 0.9 * 2.0 * M_PI / (float)nbrs[nbhd_size];

      /*
       at this point, the vall list contains all the
       neighbors currently found
       at ALL distances, while the vnbrs list contains ONLY the
       nbhd_size-neighbors.
      */
      if (found <= nbrs[nbhd_size]) /* just copy them all in */
      {
        for (n = 0; n < found; n++, vt->vtotal++) {
          vt->v[vt->vtotal] = vnbrs[n];
          v->dist_orig[vt->vtotal] = mris->vertices[vnbrs[n]].d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
        }
      }
      else /* randomly sample from them */
      {
        int vstart = vt->vtotal;
        for (n = 0; n < nbrs[nbhd_size]; n++, vt->vtotal++) {
          int j, niter = 0;
          do {
            do {
              i = nint(randomNumber(0.0, (double)found - 1));
            } while (vnbrs[i] < 0);
            /*
             now check to make sure that the angle between this
             point and the others already selected is not too
             small to make sure the points are not bunched.
            */
            VERTEX const * const vn = &mris->vertices[vnbrs[i]];
            VECTOR_LOAD(v1, vn->x - v->x, vn->y - v->y, vn->z - v->z);
            done = 1;
            for (j = vstart; done && j < vt->vtotal; j++) {
              VERTEX const * const vn2 = &mris->vertices[vt->v[j]];
              VECTOR_LOAD(v2, vn2->x - v->x, vn2->y - v->y, vn2->z - v->z);
              angle = Vector3Angle(v1, v2);
              if (angle < min_angle) {
                done = 0;
              }
            }
            if (++niter > found) /* couldn't find enough at this difference */
            {
              min_angle *= 0.75f; /* be more liberal */
              niter = 0;
            }
          } while (!done && !FZERO(min_angle));
          VERTEX const * const vn = &mris->vertices[vnbrs[i]];
          vt->v[vt->vtotal] = vnbrs[i];
          v->dist_orig[vt->vtotal] = vn->d;
          if (v->dist_orig[vt->vtotal] > 60) {
            DiagBreak();
          }
          if (FZERO(vn->d)) {
            DiagBreak();
          }
          vnbrs[i] = -1;
        }
      }
    }

    if ((vno > 9.0 * mris->nvertices / 10.0) && 0) {
      if (checks++ == 0) {
        printf("checking surface at vno %d\n", vno);
      }
      if (mrisCheckSurfaceNbrs(mris) != 1) {
        printf("surface bad at vertex %d!\n", vno);
        DiagBreak();
      }
    }

    if ((Gdiag_no == vno) && DIAG_VERBOSE_ON) {
      FILE *fp;
      char fname[STRLEN];

      sprintf(fname, "v%d", vno);
      fp = fopen(fname, "w");
      fprintf(fp, "%d\n", vall_num);
      for (n = 0; n < vall_num; n++) {
        fprintf(fp, "%d\n", vall[n]);
      }
      fclose(fp);

      sprintf(fname, "vn%d", vno);
      fp = fopen(fname, "w");
      fprintf(fp, "%d\n", vt->vtotal);
      for (n = 0; n < vt->vtotal; n++) {
        fprintf(fp, "%d\n", vt->v[n]);
      }
      fclose(fp);
      for (n = 0; n < mris->nvertices; n++) {
        VERTEX * const vn = &mris->vertices[n];
#if 0
        if (vn->ripflag)
        {
          continue ;
        }
#endif
        vn->curv = vn->d;
      }
      sprintf(fname, "%s.dist", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh");
      MRISwriteCurvature(mris, fname);
    }

    /*
     done building arrays - allocate distance vectors and
     sample from the found neighbors list.
    */
    /* now unmark them all */
    for (n = 0; n < vall_num; n++) {
      mris->vertices[vall[n]].marked = 0;
      mris->vertices[vall[n]].d = 0.0;
    }

    total_nbrs += vt->vtotal;
  }

  /* now fill in immediate neighborhood(Euclidean) distances */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    if (vt->nsize == 3) {
      vtotal = vt->v3num;
    }
    else if (vt->nsize == 2) {
      vtotal = vt->v2num;
    }
    else {
      vtotal = vt->vnum;
    }
    for (n = 0; n < vtotal; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      xd = v->x - vn->x;
      yd = v->y - vn->y;
      zd = v->z - vn->z;
      v->dist_orig[n] = sqrt(xd * xd + yd * yd + zd * zd);
    }
  }

  // make sure distances are symmetric
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    for (n = 0; n < vt->vtotal; n++) {
      VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vt->v[n]];
      VERTEX                * const vn  = &mris->vertices         [vt->v[n]];
      if (vn->ripflag) {
        continue;
      }
      for (i = 0; i < vnt->vtotal; i++) {
        if (vnt->v[i] == vno)  // distance in both lists - make it the average
        {
          double dist;
          dist = (vn->dist_orig[i] + v->dist_orig[n]) / 2;
          vn->dist_orig[i] = dist;
          v->dist_orig[n] = dist;
          break;
        }
      }
    }
  }

  /* check reasonableness of distances */
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    for (n = 0; n < vt->vtotal; n++) {
      if (DZERO(v->dist_orig[n])) fprintf(stderr, "zero distance at v %d, n %d (vn = %d)\n", vno, n, vt->v[n]);
    }
  }

  if (Gdiag_no >= 0) {
    FILE *fp;
    char fname[STRLEN];
    int i;

    sprintf(fname, "v%d.log", Gdiag_no);
    fp = fopen(fname, "w");

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[Gdiag_no];    
    VERTEX          const * const v  = &mris->vertices         [Gdiag_no];
    for (i = 0; i < vt->vtotal; i++) {
      fprintf(fp, "%d: %d, %f\n", i, vt->v[i], v->dist_orig[i]);
    }
    fclose(fp);
  }

  mris->avg_nbrs = (float)total_nbrs / (float)MRISvalidVertices(mris);
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "avg_nbrs = %2.1f\n", mris->avg_nbrs);
  }

#if MULTI_DIST_SCALING
  if (Gdiag & DIAG_SHOW) {
    for (n = 0; n <= max_nbhd; n++) {
      if (nc[n]) {
        c[n] /= (float)nc[n];
      }
      fprintf(stdout, "c[%d] = %2.5f (%d samples)\n", n, c[n], nc[n]);
    }
    fprintf(stdout, "c[] = { ");
    for (n = 0; n <= max_nbhd; n++) {
      fprintf(stdout, "%2.5f", c[n]);
      if (n < max_nbhd) {
        fprintf(stdout, ", ");
      }
    }
  }
#endif
  free(vnbrs);
  free(vall);
  free(vnb);
  free(old_dist);
  free(old_v);
  VectorFree(&v1);
  VectorFree(&v2);
  if (Gdiag & DIAG_HEARTBEAT) {
    fprintf(stdout, " done.\n");
  }

  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRIScanonicalToWorld(MRI_SURFACE *mris, double phi, double theta, double *pxw, double *pyw, double *pzw)
{
  double x, y, z, radius;

  radius = mris->radius;
  *pxw = x = radius * sin(phi) * cos(theta);
  *pyw = y = radius * sin(phi) * sin(theta);
  *pzw = z = radius * cos(phi);
  return (NO_ERROR);
}


int MRISallocExtraGradients(MRI_SURFACE *mris)
{
  if (mris->dx2) {
    return (NO_ERROR); /* already allocated */
  }

  mris->dx2 = (float *)calloc(mris->nvertices, sizeof(float));
  mris->dy2 = (float *)calloc(mris->nvertices, sizeof(float));
  mris->dz2 = (float *)calloc(mris->nvertices, sizeof(float));
  if (!mris->dx2 || !mris->dy2 || !mris->dz2)
    ErrorExit(ERROR_NO_MEMORY, "MRISallocExtraGradients: could allocate gradient vectors");
  return (NO_ERROR);
}


int slprints(char *apch_txt)
{
  //
  // PRECONDITIONS
  //      o String <apch_txt> to print to stdout
  //
  // POSTCONDITIONS
  //      o Pre-pends a syslogd style prefix to <apch_txt>:
  //              <date> <hostname> <apch_txt>
  //

  struct tm *ptm_local;
  time_t t;
  char pch_hostname[255];
  size_t len_hostname;
  char *pch_timeMon = NULL;
  char *pch_time = NULL;
  char pch_output[65536];

  t = time(NULL);
  len_hostname = 255;
  gethostname(pch_hostname, len_hostname);
  strcpy(pch_hostname, strtok(pch_hostname, "."));
  strcpy(pch_output, "");
  ptm_local = localtime(&t);
  pch_timeMon = strtok(asctime(ptm_local), "\n");
  pch_time = strdup(pch_timeMon + 4);
  sprintf(pch_output, "%s %s", pch_time, pch_hostname);
  sprintf(pch_output, "%s %s", pch_output, apch_txt);
  printf("%s", pch_output);
  return strlen(pch_output);
}

void cprints(char *apch_left, char *apch_right)
{
  //
  // PRECONDITIONS
  //  o The length of each text string should be such that
  //    the string will "fit" into its column.
  //
  // POSTCONDITIONS
  //  o Column prints the left (action) and right (status)
  //    text strings in a formatted manner.
  //

  static char pch_right[16384];

  sprintf(pch_right, " [ %s ]\n", apch_right);
  if (strlen(apch_left)) {
    fprintf(stderr, "%*s", G_LC, apch_left);
  }
  if (strlen(apch_right)) {
    fprintf(stderr, "%*s", G_RC, pch_right);
  }
  fflush(stderr);
}

void cprintd(char *apch_left, int a_right)
{
  //
  // PRECONDITIONS
  //  o The length of each text string should be such that
  //    the string will "fit" into its column.
  //
  // POSTCONDITIONS
  //  o Column prints the left (action) and right (status)
  //    text strings in a formatted manner.
  //

  static char pch_right[16384];

  sprintf(pch_right, " [ %d ]\n", a_right);
  if (strlen(apch_left)) {
    fprintf(stderr, "%*s", G_LC, apch_left);
  }
  else {
    fprintf(stderr, "%*s", G_LC, " ");
  }
  fprintf(stderr, "%*s", G_RC, pch_right);
  fflush(stderr);
}

void cprintf(char *apch_left, float af_right)
{
  //
  // PRECONDITIONS
  //  o The length of each text string should be such that
  //    the string will "fit" into its column.
  //
  // POSTCONDITIONS
  //  o Column prints the left (action) and right (status)
  //    text strings in a formatted manner.
  //

  static char pch_right[16384];

  sprintf(pch_right, " [ %f ]\n", af_right);
  if (strlen(apch_left)) {
    fprintf(stderr, "%*s", G_LC, apch_left);
  }
  else {
    fprintf(stderr, "%*s", G_LC, " ");
  }
  fprintf(stderr, "%*s", G_RC, pch_right);
  fflush(stderr);
}

short VECTOR_elementIndex_findNotEqual(VECTOR *apV, float af_searchTerm)
{
  //
  // PRECONDITIONS
  //  o The <apV> is a column vector.
  //
  // POSTCONDITIONS
  //  o The index of the first element in <apV> that is not equal to
  //    the <af_searchTerm> is returned in the function name. If no
  //    hits are found, a '0' is returned.
  //  o NOTE that the first element index in the vector is '1', not
  //    zero (i.e. MatLAB convention).
  //

  int i = 0;
  short b_ret = 0;
  for (i = 0; i < apV->rows; i++)
    if (VECTOR_ELT(apV, i + 1) != af_searchTerm) {
      b_ret = i;
      break;
    }
  return b_ret;
}

short VECTOR_elementIndex_find(VECTOR *apV, float af_searchTerm)
{
  //
  // PRECONDITIONS
  //  o The <apV> is a column vector.
  //
  // POSTCONDITIONS
  //  o The index of the first element in <apV> that is equal to
  //    the <af_searchTerm> is returned in the function name. If no
  //    hits are found, a '0' is returned.
  //  o NOTE that the first element index in the vector is '1', not
  //    zero (i.e. MatLAB convention).
  //

  int i = 0;
  short b_ret = 0;
  for (i = 1; i <= apV->rows; i++)
    if (VECTOR_ELT(apV, i) == af_searchTerm) {
      b_ret = i;
      break;
    }
  return b_ret;
}


short MRIS_vertexProgress_print(MRIS *apmris, int avertex, char *apch_message)
{
  //
  // PRECONDITIONS
  //  o <avertex> is the current vertex being processed in a stream.
  //  o If <apch_message> is non-NULL, then prefix the progress bar
  //    with <apch_message> (and terminate progress bar with [ ok ]).
  //
  // POSTCONDITIONS
  //  o For every 5% of processed vertices a "#" is written to stderr
  //

  static int totalVertices = 0;
  static int fivePerc = 0;

  totalVertices = apmris->nvertices;
  fivePerc = 0.05 * totalVertices;

  if (!avertex) {
    if (apch_message != NULL) {
      fprintf(stderr, "%*s", G_LC, apch_message);
    }
    fprintf(stderr, " [");
  }
  if (avertex % fivePerc == fivePerc - 1) {
    fprintf(stderr, "#");
  }
  if (avertex == apmris->nvertices - 1) {
    fprintf(stderr, "] ");
    if (apch_message != NULL) {
      fprintf(stderr, "%*s\n", 1, "[ ok ]");
    }
  }
  return 1;
}


int signum_eval(float af)
{
  //
  // PRECONDITIONS
  //  o <af> is an input float.
  //
  // POSTCONDITIONS
  //  o if <af> < 0, a -1 is returned, else +1 is returned
  //

  return (af < 0) ? -1 : 1;
}


/*!
  \fn int CompareVertexCoords(const void *v1, const void *v2)
  \brief Compares the xyz coords from two vertices in the VERTEX_SORT structure. This
  if a function that can be used in qsort to sort vertices by xyz.
 */
int CompareVertexCoords(const void *v1, const void *v2)
{
  VERTEX_SORT *vtx1, *vtx2;
  vtx1 = (VERTEX_SORT*) v1;
  vtx2 = (VERTEX_SORT*) v2;
  // Traditionally, vertices tended to be sorted in the ydir
  if(vtx1->y < vtx2->y) return(-1);
  if(vtx1->y > vtx2->y) return(+1);
  if(vtx1->x < vtx2->x) return(-1);
  if(vtx1->x > vtx2->x) return(+1);
  if(vtx1->z < vtx2->z) return(-1);
  return(+1);
}

/*!
  \fn int CompareFaceVertices(const void *f1, const void *f2);
  \brief Compares the vertex indices from two faces in the FACE_SORT
  structure. This if a function that can be used in qsort to sort
  faces by vertex index.
 */
int CompareFaceVertices(const void *vf1, const void *vf2)
{
  FACE_SORT *f1, *f2;
  f1 = (FACE_SORT*) vf1;
  f2 = (FACE_SORT*) vf2;
  if(f1->v0 < f2->v0) return(-1);
  if(f1->v0 > f2->v0) return(+1);
  if(f1->v1 < f2->v1) return(-1);
  if(f1->v1 > f2->v1) return(+1);
  if(f1->v2 < f2->v2) return(-1);
  return(+1);
}


#ifdef FS_CUDA_TIMINGS
/* this stuff is needed for some of the benchmarking I have been doing */
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
static int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}
#endif  // FS_CUDA_TIMINGS
