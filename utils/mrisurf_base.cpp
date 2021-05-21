#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
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
#include "mrisurf_base.h"
#include "mrisurf_metricProperties.h"


int mris_sort_compare_float(const void *pc1, const void *pc2)
{
  float c1, c2;

  c1 = *(float *)pc1;
  c2 = *(float *)pc2;

  /*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (c1 > c2) {
    return (1);
  }
  else if (c1 < c2) {
    return (-1);
  }

  return (0);
}


const char* MRIS_Status_text(MRIS_Status s1)
{
    switch (s1) {
#define SEP
#define ELT(S,D) case S : return #S;
    MRIS_Status_ELTS
#undef ELT
#undef SEP
    default: cheapAssert(false); return "<Bad MRIS_Status>";
    }
}


MRIS_Status_DistanceFormula MRIS_Status_distanceFormula(MRIS_Status s1)
{
    switch (s1) {
#define SEP
#define ELT(S,D) case S : return MRIS_Status_DistanceFormula_##D;
    MRIS_Status_ELTS
#undef ELT
#undef SEP
    default: cheapAssert(false); return MRIS_Status_DistanceFormula_0;
    }
}


bool areCompatible(MRIS_Status s1, MRIS_Status s2)
{
  return MRIS_Status_distanceFormula(s1) == MRIS_Status_distanceFormula(s2);
}


void checkOrigXYZCompatibleWkr(MRIS_Status s1, MRIS_Status s2, const char* file, int line) {
  if (areCompatible(s1,s2)) return;
  fs::debug() << "using incompatible " << MRIS_Status_text(s1) << " and " << MRIS_Status_text(s2);
}


int  MRIS_acquireTemp(MRIS* mris, MRIS_TempAssigned temp) {
  int const bits = 1 << temp;
  int const * tc = &mris->tempsAssigned;
  cheapAssert(!(bits & *tc));
  int* tv = (int*)tc;
  *tv |= bits;
  return 123;   // hack
}

void MRIS_checkAcquiredTemp(MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result) {
  int const bits = 1 << temp;
  int const * tc = &mris->tempsAssigned;
  cheapAssert((bits & *tc));
}


void MRIS_releaseTemp(MRIS* mris, MRIS_TempAssigned temp, int MRIS_acquireTemp_result) {
  int const bits = 1 << temp;
  int const * tc = &mris->tempsAssigned;
  cheapAssert((bits & *tc));
  int* tv = (int*)tc;
  *tv &= ~bits;
}


// Create temps, and don't let the nvertices change until it is freed
//
float* MRISmakeFloatPerVertex(MRIS *mris) {
  MRISacquireNverticesFrozen(mris);
  float* p = (float*)malloc(mris->nvertices*sizeof(float));
  return p;  
}

void MRISfreeFloatPerVertex(MRIS *mris, float** pp) {
  freeAndNULL(*pp);
  MRISreleaseNverticesFrozen(mris);
}


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
  char * copy = (char *)calloc(strlen(fname) + 1, sizeof(char));
  strcpy(copy, fname);
  curvature_names[0] = copy;
  return (NO_ERROR);
}

int MRISsetSulcFileName(const char *sulc_name)
{
  char * copy = (char *)calloc(strlen(sulc_name) + 1, sizeof(char));
  strcpy(copy, sulc_name);
  curvature_names[1] = copy;
  return (NO_ERROR);
}

int MRISsetOriginalFileName(const char *orig_name)
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


// dist and dist_orig must be managed separately because 
//      changing xyz        invalidates dist 
//      changing origxyz    invalidates dist_orig
//      changing v          invalidates both
//
// By keeping their storage allocation independent of the VERTEX
// it is very fast to change VERTEX::dist and VERTEX::dist_orig between NULL and !NULL
// so the code can make the ->dist or ->dist_orig become NULL, making the bad values immediately inaccessible
// and then quickly having them available again when they are computed
//     
static void freeDistOrDistOrig(bool doOrig, VERTEX* v) {
  if (doOrig) { *(void**)(&v->dist_orig) = NULL; v->dist_orig_capacity = 0; }
  else        { *(void**)(&v->dist     ) = NULL; v->dist_capacity      = 0; }
}

static void changeDistOrDistOrig(bool doOrig, MRIS *mris, int vno, int oldSize, int neededCapacity) 
{
  VERTEX const * const v = &mris->vertices[vno];

  const int    * pcCap = doOrig ? &v->dist_orig_capacity : &v->dist_capacity; 
  float* const * pc    = doOrig ? &v->dist_orig          : &v->dist;

  if (*pcCap < neededCapacity) {
  
    void** p_storage = doOrig ? &mris->dist_orig_storage[vno] : &mris->dist_storage[vno];
    float* alloced = (float*)realloc(*p_storage, neededCapacity*sizeof(float));
    *p_storage = (void*)alloced;

    if (!*pc) {
      char const flag = (char)(1)<<doOrig;
      mris->dist_alloced_flags |= flag;
    }
    
    *(float**)pc    = alloced;
    *(int   *)pcCap = neededCapacity;
  }
  
  if (neededCapacity > oldSize) {
    bzero((*pc) + oldSize, (neededCapacity-oldSize)*sizeof(float));
  }
}

float* mrisStealDistStore(MRIS* mris, int vno, int capacity) {

  // mrisSetDist may be called later, so stealing it now (if it is not being used) 
  // will avoid the need to malloc some now and to free later
  //
  VERTEX const * const v = &mris->vertices[vno];
  
  float* p = NULL;
  if (v->dist == NULL) { p = (float *)mris->dist_storage[vno]; mris->dist_storage[vno] = NULL; }
  
  return (float*)realloc(p, capacity*sizeof(float));
}

void mrisSetDist(MRIS* mris, int vno, float* dist, int newCapacity) {
  // Some code has already computed the dist into a vector.
  // Here is a quick way of putting it into the vertex
  // The capacity must not be smaller than the current capacity
  //
  bool const doOrig = false;
  char const flag = (char)(1)<<doOrig;
  mris->dist_alloced_flags |= flag;
    
  VERTEX const * const v = &mris->vertices[vno];

  cheapAssert(newCapacity >= v->dist_capacity);

  float * const * pc = &v->dist;
  float *       * p  = (float**)pc;
  
  if (mris->dist_storage[vno]) free(mris->dist_storage[vno]);   // very frequently avoid the lock
  mris->dist_storage[vno] = dist;
  *p = dist;
  
  const int * pcCap = &v->dist_capacity; *(int*)pcCap = newCapacity;
}

static void growDistOrDistOrig(bool doOrig, MRIS *mris, int vno, int minimumCapacity)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  int oldSize = doOrig ? v->dist_orig_capacity : vt->vtotal;
  int vSize = mrisVertexVSize(mris, vno);
  if (vSize < minimumCapacity) vSize = minimumCapacity;
  changeDistOrDistOrig(doOrig, mris, vno, oldSize, vSize);
}

void MRISgrowDist(MRIS *mris, int vno, int minimumCapacity)
{
  growDistOrDistOrig(false, mris, vno, minimumCapacity);
}

void MRISgrowDistOrig(MRIS *mris, int vno, int minimumCapacity)
{
  growDistOrDistOrig(true, mris, vno, minimumCapacity);
}

static void makeDistOrDistOrig(bool doOrig, MRIS *mris, int vno) 
{
  int vSize = mrisVertexVSize(mris, vno);
  changeDistOrDistOrig(doOrig, mris, vno, 0, vSize);
}

void MRISmakeDist(MRIS *mris, int vno) 
{
  makeDistOrDistOrig(false, mris, vno);
}

void MRISmakeDistOrig(MRIS *mris, int vno) 
{
  makeDistOrDistOrig(true, mris, vno);
}


// Freeing
//
static void freeDistsOrDistOrigs(bool doOrig, MRIS *mris)
{
  char const flag = (char)(1)<<doOrig;

  if (!(mris->dist_alloced_flags & flag)) return;
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    freeDistOrDistOrig(doOrig, &mris->vertices[vno]);
  }
  
  mris->dist_alloced_flags ^= flag;
}

void MRISfreeDistsButNotOrig(MRIS* mris)
  // Maybe should not be here...
{
  mris->dist_nsize = 0;
  freeDistsOrDistOrigs(false,mris);
}

void MRISfreeDistOrigs(MRIS* mris)
  // Maybe should not be here...
{
  freeDistsOrDistOrigs(true,mris);
}

void MRISfreeDistsButNotOrig(MRIS_MP* mris)
{
  cheapAssert(!"MRISfreeDistsButNotOrig(MRIS_MP* mris) NYI");  
}
void MRISfreeDistsButNotOrig(MRISPV* mris)
{
  cheapAssert(!"MRISfreeDistsButNotOrig(MRISPV* mris) NYI");  
}


// VERTEX, FACE, EDGE primitives
//
static void VERTEX_TOPOLOGYdtr(VERTEX_TOPOLOGY* v) {
  freeAndNULL(v->v);
  freeAndNULL(v->f);
  freeAndNULL(v->n);
  if (v->e) freeAndNULL(v->e);
}

static void VERTEXdtr(VERTEX* v) {
  freeDistOrDistOrig(false, v);
  freeDistOrDistOrig(true,  v);
}

static void FACEdtr(FACE* f) {
  if (f->norm) DMatrixFree(&f->norm);
  int k;
  for (k=0; k < 3; k++){
    if (f->gradNorm[k]) DMatrixFree(&f->gradNorm[k]);
  }
}

static void MRIS_EDGEdtr(MRI_EDGE* e) {
  for (int k=0; k<4; k++){
    if (e->gradDot[k]) DMatrixFree(&e->gradDot[k]);
  }
}

static void MRIS_CORNERdtr(MRI_CORNER* c) {
  for (int k=0; k<3; k++){
    if (c->gradDot[k]) DMatrixFree(&c->gradDot[k]);
  }
}

// MRIS 
//
void mrisDumpVertex(FILE* file, MRIS const * mris, int vno) {
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];
  
  fprintf(file, "vno:%d ripflag:%d marked:%d nsizeCur:%d nsizeMax:%d vnum:%d v2num:%d v3num:%d vtotal:%d nsizeCur:%d nsizeMax:%d x:%f y:%f z:%f",
    vno, v->ripflag, v->marked, vt->nsizeCur, vt->nsizeMax, vt->vnum, vt->v2num, vt->v3num, vt->vtotal, vt->nsizeCur, vt->nsizeMax, v->x, v->y, v->z);
  
  int const vsize = mrisVertexVSize(mris, vno);
  int i;
  for (i = 0; i < vsize; i++) {
    fprintf(file, " [%d] v:%d dist:%f dist_orig:%f",
      i, vt->v[i], v->dist ? v->dist[i] : -1.0f, v->dist_orig ? v->dist_orig[i] : -1.0f);
  } fprintf(file, "\n");

  for (i = 0; i < vt->num; i++) {
    fprintf(file, " [%d] f:%d n:%d",
      i, vt->f[i], vt->n[i]);
  } fprintf(file, "\n");
  
}

void mrisDumpFace(MRIS const * mris, int fno, FILE* file) {
  FACE const * const f = &mris->faces[fno];
  fprintf(file, "fno:%d ripflag:%d marked:%d\n",
    fno, f->ripflag, f->marked);
  int i;
  for (i = 0; i < VERTICES_PER_FACE; i++) {
    fprintf(file, " [%d] v:%d",
      i, f->v[i]);
  } fprintf(file, "\n");
}

void mrisDumpShape(FILE* file, MRIS const * mris) {
  fprintf(file, "mrisDumpShape {\n");
  fprintf(file, "nvertices:%d nfaces:%d max_nsize:%d nsize:%d\n",
    mris->nvertices, mris->nfaces, mris->max_nsize, mris->nsize);
  int vno;
  for (vno = 0; vno < MIN(10,mris->nvertices); vno++) mrisDumpVertex(file, mris, vno);
  int fno;
  for (fno = 0; fno < MIN(10,mris->nfaces); fno++) mrisDumpFace(mris, fno, file);
  fprintf(file, "} // mrisDumpShape\n");
}

static void MRISchangedNFacesNVertices(MRIS * mris, bool scrambled) {
  // useful for debugging
}

void MRISacquireNverticesFrozen(MRIS *mris) {
  #pragma omp atomic
  (*(int*)&mris->nverticesFrozen)++;
}
void MRISreleaseNverticesFrozen(MRIS *mris) {
  #pragma omp atomic
  (*(int*)&mris->nverticesFrozen)--;
}

bool MRISreallocVertices(MRIS * mris, int max_vertices, int nvertices) {

  // mris->max_vertices should only ever grow
  // because this code does not know how to free existing ones
  //
  if (max_vertices < mris->max_vertices) max_vertices = mris->max_vertices;
  
  cheapAssert(nvertices >= 0);
  cheapAssert(max_vertices >= nvertices);
  cheapAssert(!mris->nverticesFrozen);
  
  mris->vertices = (VERTEX *)realloc(mris->vertices, max_vertices*sizeof(VERTEX));
  if (!mris->vertices) return false;
  #ifndef SEPARATE_VERTEX_TOPOLOGY
    mris->vertices_topology = mris->vertices;
  #else
    mris->vertices_topology = (VERTEX_TOPOLOGY *)realloc(mris->vertices_topology, max_vertices*sizeof(VERTEX_TOPOLOGY));
  #endif
  
  mris->dist_storage      = (void**)realloc(mris->dist_storage,      max_vertices*sizeof(void*));
  mris->dist_orig_storage = (void**)realloc(mris->dist_orig_storage, max_vertices*sizeof(void*));
  
  int change = max_vertices - mris->nvertices;
  if (change > 0) { // all above nvertices must be zero'ed, since MRISgrowNVertices does not...
    bzero(mris->vertices          + mris->nvertices,     change*sizeof(VERTEX));
#ifdef SEPARATE_VERTEX_TOPOLOGY
    bzero(mris->vertices_topology + mris->nvertices,     change*sizeof(VERTEX_TOPOLOGY));
#endif
  }
  
  change = max_vertices - mris->max_vertices;
  if (change > 0) {
    bzero(mris->dist_storage      + mris->max_vertices, change*sizeof(void*));
    bzero(mris->dist_orig_storage + mris->max_vertices, change*sizeof(void*));
  }

  *(int*)(&mris->max_vertices) = max_vertices;    // get around const
  *(int*)(&mris->nvertices)    = nvertices;       // get around const
    
  notifyActiveRealmTreesChangedNFacesNVertices(mris);
  
  return true;
}

void MRISgrowNVertices(MRIS * mris, int nvertices) {
  cheapAssert(nvertices <= mris->max_vertices);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}

void MRIStruncateNVertices(MRIS * mris, int nvertices) {
  cheapAssert(mris->nvertices >= nvertices);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}

void MRISremovedVertices(MRIS * mris, int nvertices) {
  cheapAssert(mris->nvertices >= nvertices);
  MRISchangedNFacesNVertices(mris, true);
  *(int*)(&mris->nvertices) = nvertices;  // get around const
}



bool MRISreallocFaces(MRIS * mris, int max_faces, int nfaces) {
  cheapAssert(nfaces >= 0);
  cheapAssert(max_faces >= nfaces);
  
  mris->faces  =
    (FACE *)realloc(mris->faces, max_faces*sizeof(FACE));
  cheapAssert(mris->faces);

  mris->faceNormCacheEntries =
    (FaceNormCacheEntry*)realloc(mris->faceNormCacheEntries, max_faces*sizeof(FaceNormCacheEntry));
  cheapAssert(mris->faceNormCacheEntries);
 
  mris->faceNormDeferredEntries =
    (FaceNormDeferredEntry*)realloc(mris->faceNormDeferredEntries, max_faces*sizeof(FaceNormDeferredEntry));
  cheapAssert(mris->faceNormDeferredEntries);
  
  
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

bool MRISallocateFaces(MRIS * mris, int nfaces) {
    return MRISreallocFaces(mris, nfaces, nfaces);
}

void MRISgrowNFaces(MRIS * mris, int nfaces) {
  cheapAssert(nfaces <= mris->max_faces);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}

void MRIStruncateNFaces(MRIS * mris, int nfaces) {
  cheapAssert(mris->nfaces >= nfaces);
  MRISchangedNFacesNVertices(mris, false);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}

void MRISremovedFaces(MRIS * mris, int nfaces) {
  cheapAssert(mris->nfaces >= nfaces);
  MRISchangedNFacesNVertices(mris, true);
  *(int*)(&mris->nfaces) = nfaces;  // get around const
}


void MRISoverAllocVerticesAndFaces(MRIS* mris, int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRISchangedNFacesNVertices(mris, false);
  cheapAssert(nvertices >= 0);
  cheapAssert(nfaces    >= 0);

  if (max_vertices <= nvertices) max_vertices = nvertices;
  if (max_faces    <= nfaces   ) max_faces    = nfaces;

  MRISreallocVertices(mris, max_vertices, nvertices);
  MRISreallocFaces(mris, max_faces, nfaces);
}


void MRISreallocVerticesAndFaces(MRIS *mris, int nvertices, int nfaces) {
  MRISoverAllocVerticesAndFaces(mris, nvertices, nfaces, nvertices, nfaces);
}


MRIS* MRIScopyMetadata(MRIS const * source, MRIS * target)
{
  if (!target) target = MRISalloc(source->nvertices, source->nfaces);

  target->type = source->type;
  target->status = source->status;
  target->hemisphere = source->hemisphere;

  copyVolGeom(&source->vg, &target->vg);

  if (source->lta) target->lta = LTAcopy(source->lta, nullptr);
  if (source->SRASToTalSRAS_) target->SRASToTalSRAS_ = MatrixCopy(source->SRASToTalSRAS_, nullptr);
  if (source->TalSRASToSRAS_) target->TalSRASToSRAS_ = MatrixCopy(source->TalSRASToSRAS_, nullptr);

  return target;
}


void MRISctr(MRIS *mris, int max_vertices, int max_faces, int nvertices, int nfaces) {
  // This should be the ONLY place where an MRIS is constructed
  //
  bzero(mris, sizeof(MRIS));
  MRISoverAllocVerticesAndFaces(mris, max_vertices, max_faces, nvertices, nfaces);
  
  mris->max_nsize = mris->nsize = 1;      // only 1-connected neighbors initially
}

void MRISdtr(MRIS *mris) {
  freeAndNULL(mris->dx2);
  freeAndNULL(mris->dy2);
  freeAndNULL(mris->dz2);
  freeAndNULL(mris->labels);
  
  if (mris->ct) CTABfree(&mris->ct);

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGYdtr(&mris->vertices_topology[vno]);
    VERTEXdtr         (&mris->vertices         [vno]);
  }
  freeAndNULL(mris->vertices);

#ifdef SEPARATE_VERTEX_TOPOLOGY
  freeAndNULL(mris->vertices_topology);
#endif

  freeAndNULL(mris->faceNormDeferredEntries);
  freeAndNULL(mris->faceNormCacheEntries);
  
  int faceno;
  for(faceno = 0; faceno < mris->nfaces; faceno++){
    FACEdtr(&mris->faces[faceno]);
  }
  freeAndNULL(mris->faces);

  int edgeNo;
  for (edgeNo=0; edgeNo < mris->nedges; edgeNo++) {
    MRIS_EDGEdtr(&mris->edges[edgeNo]);
  }
  freeAndNULL(mris->edges);

  int cornerNo;
  for (cornerNo=0; cornerNo < mris->ncorners; cornerNo++) {
    MRIS_CORNERdtr(&mris->corners[cornerNo]);
  }
  freeAndNULL(mris->corners);

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

  int i;
  for (i = 0; i < mris->ncmds; i++) {
    freeAndNULL(mris->cmdlines[i]);
  }

  if (mris->m_sras2vox) {
    MatrixFree(&mris->m_sras2vox);
  }

  for (vno = 0; vno < mris->max_vertices; vno++) {
    freeAndNULL(mris->dist_storage     [vno]);
    freeAndNULL(mris->dist_orig_storage[vno]);
  }

  freeAndNULL(mris->dist_storage);
  freeAndNULL(mris->dist_orig_storage);
}


// Heap allocators are simply malloc and ctr the object, freeing is the opposite
//
MRIS *MRISoverAlloc(int max_vertices, int max_faces, int nvertices, int nfaces)
{
  MRIS* mris = (MRIS *)malloc(sizeof(MRIS));
  if (!mris) ErrorExit(ERROR_NO_MEMORY, 
                "MRISalloc(%d, %d): could not allocate mris structure", max_vertices, max_faces);

  MRISctr(mris, max_vertices, max_faces, nvertices, nfaces);

  return (mris);
}

MRIS* MRISalloc(int nvertices, int nfaces)
{
  return MRISoverAlloc(nvertices, nfaces, nvertices, nfaces);
}

void MRISfree(MRIS **pmris)
{
  MRIS *mris = *pmris;
  *pmris = NULL;
  MRISdtr(mris);
  free(mris);
}


//======
// ripped
//
char* MRISexportFaceRipflags(MRIS* mris) {
  char* flags = (char*)malloc(mris->nfaces * sizeof(char));
  int fno;
  for (fno = 0 ; fno < mris->nfaces ; fno++){
    flags[fno] = mris->faces[fno].ripflag;
  }
  return flags;
}

void  MRISimportFaceRipflags(MRIS* mris, const char* flags) {
  int fno;
  for (fno = 0 ; fno < mris->nfaces ; fno++){
    mris->faces[fno].ripflag = flags[fno];
  }
}


char* MRISexportVertexRipflags(MRIS* mris) {
  char* flags = (char*)malloc(mris->nvertices * sizeof(char));
  int vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    flags[vno] = mris->vertices[vno].ripflag;
  }
  return flags;
}

void MRISimportVertexRipflags(MRIS* mris, const char* flags) {
  int vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    mris->vertices[vno].ripflag = flags[vno];
  }
} 

int MRIScountRipped(MRIS *mris)
{
  int nripped=0, vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++){
    if(mris->vertices[vno].ripflag) nripped++;
  }
  return(nripped);
}


int MRISstoreRipFlags(MRIS *mris)
{
  int vno, fno;
  VERTEX *v;
  FACE *f;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->oripflag = v->ripflag;
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    f->oripflag = f->ripflag;
  }
  return (NO_ERROR);
}

int MRISrestoreRipFlags(MRIS *mris)
{
  int vno, fno;
  VERTEX *v;
  FACE *f;

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    v->ripflag = v->oripflag;
  }

  for (fno = 0; fno < mris->nfaces; fno++) {
    f = &mris->faces[fno];
    f->ripflag = f->oripflag;
  }

  return (NO_ERROR);
}

int MRISunrip(MRIS *mris)
{
  int vno, fno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    mris->vertices[vno].ripflag = 0;
  }
  for (fno = 0; fno < mris->nfaces; fno++) {
    mris->faces[fno].ripflag = 0;
  }
  return (NO_ERROR);
}


/*
  Returns an (nvertices x 3) array of the surface's vertex positions.
*/
float* MRISgetVertexArray(MRIS *mris)
{
  float * const buffer = new float[mris->nvertices * 3];
  float * ptr = buffer;
  for (int n = 0 ; n < mris->nvertices ; n++) {
    *ptr++ = mris->vertices[n].x;
    *ptr++ = mris->vertices[n].y;
    *ptr++ = mris->vertices[n].z;
  }
  return buffer;
}


/*
  Returns an (nvertices x 3) array of the surface's vertex normals.
*/
float* MRISgetVertexNormalArray(MRIS *mris)
{
  float * const buffer = new float[mris->nvertices * 3];
  float * ptr = buffer;
  for (int n = 0 ; n < mris->nvertices ; n++) {
    *ptr++ = mris->vertices[n].nx;
    *ptr++ = mris->vertices[n].ny;
    *ptr++ = mris->vertices[n].nz;
  }
  return buffer;
}


/*
  Returns an (nfaces x 3) array of the surface's face indices.
*/
int* MRISgetFaceArray(MRIS *mris)
{
  int * const buffer = new int[mris->nfaces * 3];
  int * ptr = buffer;
  for (int n = 0 ; n < mris->nfaces ; n++) {
    *ptr++ = mris->faces[n].v[0];
    *ptr++ = mris->faces[n].v[1];
    *ptr++ = mris->faces[n].v[2];
  }
  return buffer;
}


/*
  Returns an (nfaces x 3) array of the surface's face normals.
*/
float* MRISgetFaceNormalArray(MRIS *mris)
{
  float * const buffer = new float[mris->nfaces * 3];
  float * ptr = buffer;
  for (int n = 0 ; n < mris->nfaces ; n++) {
    FaceNormCacheEntry const * norm = getFaceNorm(mris, n);
    *ptr++ = norm->nx;
    *ptr++ = norm->ny;
    *ptr++ = norm->nz;
  }
  return buffer;
}


char const * mrisurf_surface_names[3] = {"inflated", "smoothwm", "smoothwm"};
char const * curvature_names      [3] = {"inflated.H", "sulc", NULL};

int MRISsetCurvatureName(int nth, char const *name)
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
  for (unsigned int k = 0; k < sizeof(curvature_names) / sizeof(curvature_names[0]); k++) {
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


VOXEL_LIST **vlst_alloc(MRIS *mris, int max_vox)
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

int vlst_free(MRIS *mris, VOXEL_LIST ***pvl)
{
  VOXEL_LIST **vl = *pvl;
  int vno;

  *pvl = NULL;
  for (vno = 0; vno < mris->nvertices; vno++) VLSTfree(&vl[vno]);
  free(vl);
  return (NO_ERROR);
}

int vlst_enough_data(MRIS *mris, int vno, VOXEL_LIST *vl, double displacement)
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
static int edgesIntersectStable(MRIS *mris, EDGE *edge1, EDGE *edge2);

#if SPHERE_INTERSECTION

/* check for intersection on the sphere */
int edgesIntersect(MRIS *mris, EDGE *edge1, EDGE *edge2)
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
static int edgesIntersectStable(MRIS *mris, EDGE *edge1, EDGE *edge2)
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
int edgesIntersect(MRIS *mris, EDGE *edge1, EDGE *edge2)
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
    MRIS *mris, float x0, float y0, float z0, float nx, float ny, float nz, MHT *mht, double *pdist, int vno)
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
    MRIS *mris, VERTEX *v, MHT *mht, double *pdist, int *flist)
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

int MRISallocExtraGradients(MRIS *mris)
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

static float checkNotNanf(float f) {
  cheapAssert(!devIsnan(f));
  return f;
}

void MRIScheckForNans(MRIS *mris)
{
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
    checkNotNanf(v->odx); checkNotNanf(v->ody); checkNotNanf(v->odz);
    checkNotNanf(v-> dx); checkNotNanf(v-> dy); checkNotNanf(v-> dz);
  }  
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
  std::string pch_output;

  t = time(NULL);
  len_hostname = 255;
  gethostname(pch_hostname, len_hostname);
  strcpy(pch_hostname, strtok(pch_hostname, "."));
  ptm_local = localtime(&t);
  pch_timeMon = strtok(asctime(ptm_local), "\n");
  pch_time = strdup(pch_timeMon + 4);
  pch_output = std::string(pch_time) + std::string(pch_hostname) + apch_txt;
  printf("%s", pch_output.c_str());
  return pch_output.size();
}

void cprints(const char *apch_left, const char *apch_right)
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

void cprintd(const char *apch_left, int a_right)
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

void cprintf(const char *apch_left, float af_right)
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


short MRIS_vertexProgress_print(MRIS *apmris, int avertex, const char *apch_message)
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
