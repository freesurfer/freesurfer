#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
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
#include "mrisurf_timeStep.h"

#include "mrisurf_sseTerms.h"
#include "mrisurf_compute_dxyz.h"
#include "mrinorm.h"

#include "mrisurf_base.h"


int (*gMRISexternalTimestep)(MRI_SURFACE *mris, INTEGRATION_PARMS *parms) = NULL;


struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t {
  float xLo, xHi, yLo, yHi, zLo, zHi;
  int nextVnoPlus1;   // Plus1 so 0 can act as NULL
};
  
typedef struct MRISAsynchronousTimeStep_optionalDxDyDzUpdate_Context {
  float xLo, xHi, yLo, yHi, zLo, zHi;
  float xSubvolLen,ySubvolLen,zSubvolLen;
  float xSubvolVerge,ySubvolVerge,zSubvolVerge;
  int*   vnoToSvi; 
  struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t* vertexInfos;
} MRISAsynchronousTimeStep_optionalDxDyDzUpdate_Context;

typedef struct MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context {
  MRISAsynchronousTimeStep_optionalDxDyDzUpdate_Context* allVertexsContext;
  unsigned long hash;
  size_t count;
  size_t limit;
  bool   trace;
  int    svi;
} MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context;

#define MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerEdge 4
static const size_t MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerThread = 
    (MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerEdge * 
     MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerEdge * 
     MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerEdge + 1);
    // one extra for the faces that intersect more than one subvol

static int MRISAsynchronousTimeStep_optionalDxDyDzUpdate_svi(
  MRISAsynchronousTimeStep_optionalDxDyDzUpdate_Context* ctx,
  float x, float y, float z) 
{
  const size_t numSubvolsPerEdge   = MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerEdge;
  int const svxLo = MIN((int)numSubvolsPerEdge-1,(int)((x - ctx->xSubvolVerge - ctx->xLo)/ctx->xSubvolLen));
  int const svyLo = MIN((int)numSubvolsPerEdge-1,(int)((y - ctx->ySubvolVerge - ctx->yLo)/ctx->ySubvolLen));
  int const svzLo = MIN((int)numSubvolsPerEdge-1,(int)((z - ctx->zSubvolVerge - ctx->zLo)/ctx->zSubvolLen));

  int const svxHi = MIN((int)numSubvolsPerEdge-1,(int)((x + ctx->xSubvolVerge - ctx->xLo)/ctx->xSubvolLen));
  int const svyHi = MIN((int)numSubvolsPerEdge-1,(int)((y + ctx->ySubvolVerge - ctx->yLo)/ctx->ySubvolLen));
  int const svzHi = MIN((int)numSubvolsPerEdge-1,(int)((z + ctx->zSubvolVerge - ctx->zLo)/ctx->zSubvolLen));

  if (svxLo != svxHi || svyLo != svyHi || svzLo != svzHi) 
    return MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerThread - 1;       // the svi for any items not deep inside a subvolume

  return svxLo*numSubvolsPerEdge*numSubvolsPerEdge + svyLo*numSubvolsPerEdge + svzLo;
}



static bool proposed_ODXYZ_valid(
    MRI_SURFACE *mris,
    MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx, 
    int vno,
    float odx, float ody, float odz) 
{
  struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t* vertexInfos
    = !ctx ? NULL : ctx->allVertexsContext->vertexInfos;
  
  if (!vertexInfos) return true;            // The lack of vertexInfo means all vertexs can move anywhere
  
  VERTEX const * v = &mris->vertices[vno];
  float const x = v->x + odx;
  float const y = v->y + ody;
  float const z = v->z + odz;

  int current_svi = ctx->allVertexsContext->vnoToSvi[vno];
  if (current_svi == MRISAsynchronousTimeStep_optionalDxDyDzUpdate_numSubvolsPerThread - 1) return true;
    // already not deep in any subvolume
    
  int proposed_svi = 
    MRISAsynchronousTimeStep_optionalDxDyDzUpdate_svi(ctx->allVertexsContext, 
      x,y,z);

  return (current_svi == proposed_svi);
    // hasn't changed subvolumes, and hence can be done as before
}

static void checkODXYZ_valid(
    MRI_SURFACE *mris,
    MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx, 
    int vno) 
{
  VERTEX const * v = &mris->vertices[vno];
  if (proposed_ODXYZ_valid(mris, ctx, vno, v->odx, v->ody, v->odz)) return; 
  struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t* vertexInfos
    = !ctx ? NULL : ctx->allVertexsContext->vertexInfos;
  float x = v->x + v->odx;
  float y = v->y + v->ody;
  float z = v->z + v->odz;
  struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t * vi = &vertexInfos[vno];
  fprintf(stdout, "%s:%d proposed to move vno:%d further than originally planned"
    "xLo:%g x:%g xHi:%g\n"
    "yLo:%g y:%g yHi:%g\n"
    "zLo:%g z:%g zHi:%g\n", 
    __FILE__, __LINE__, vno,
     vi->xLo, x, vi->xHi,
     vi->yLo, y, vi->yHi,
     vi->zLo, z, vi->zHi);
  *(int*)-1 = 0;
}

static int vertexOdFrozen = 0;
static void checkVertexOdNotFrozen(const char* file, int line) {
    if (vertexOdFrozen) {
        fprintf(stdout, "%s:%d vertexOdFrozen\n",file,line);
        *(int*)-1 = 0;
    }
}
#define CHANGES_ODXYZ checkVertexOdNotFrozen(__FILE__,__LINE__);
#define CHANGES_ODXYZ_OKAY


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MIN_NBR_DIST (0.01)

static float minNeighborDistance(MRI_SURFACE *mris) {
    return MAX(MIN_NBR_DIST, MIN_NBR_DIST * (mris->vg.xsize + mris->vg.ysize + mris->vg.zsize) / 3);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
#define MIN_MM 0.001

static bool mrisRemoveNeighborGradientComponent(MRI_SURFACE *mris, int vno, 
  MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx)
{
  float const min_nbr_dist = minNeighborDistance(mris);

  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];
  
  if (v->ripflag) {
    return true;
  }

  float const x = v->x;
  float const y = v->y;
  float const z = v->z;

  int n;
  for (n = 0; n < vt->vnum; n++) {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];

    float dx = vn->x - x;
    float dy = vn->y - y;
    float dz = vn->z - z;
    float dist = sqrt(dx * dx + dy * dy + dz * dz);

    /* too close - take out gradient component in this dir. */
    if (dist <= min_nbr_dist) {
      dx /= dist;
      dy /= dist;
      dz /= dist;
      float dot = dx * v->odx + dy * v->ody + dz * v->odz;
      if (dot > 0.0) {
        if (vno == Gdiag_no)
          printf("v %d: removing neighbor gradient dist (%2.2f, %2.2f, %2.2f) --> ", vno, v->odx, v->ody, v->odz);

        float odx = v->odx - dot * dx;
        float ody = v->ody - dot * dy;
        float odz = v->odz - dot * dz;

        if (!proposed_ODXYZ_valid(mris, ctx, vno, odx, ody, odz)) return false;
            // the proposed odxyz is beyond the range planned the parallelism can cope with

        v->odx = odx;
        v->ody = ody;
        v->odz = odz;
        
        if (vno == Gdiag_no) printf(" (%2.2f, %2.2f, %2.2f)\n", v->odx, v->ody, v->odz);
      }
    }
  }

  return true;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:    true if and only if any required movement was done.

  Description
  ------------------------------------------------------*/
static bool mrisLimitGradientDistance(
    MRI_SURFACE *mris, MHT const *mht, int vno,
    MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context * ctx)
{
  VERTEX *v = &mris->vertices[vno];

  if (v->ripflag) {
    return true;
  }

  if (!mrisRemoveNeighborGradientComponent(mris, vno, ctx)) return false;
  
  if (MHTisVectorFilled(mht, vno, v->odx, v->ody, v->odz)) {
    v->odx = v->ody = v->odz = 0.0;
    if (vno == Gdiag_no) printf("(%2.2f, %2.2f, %2.2f)\n", v->odx, v->ody, v->odz);
    v->cropped++;
    return true;
  }

  v->cropped = 0;
  return true;    
}


double mrisAdaptiveTimeStep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  double delta_t, sse, starting_sse;

  starting_sse = MRIScomputeSSE(mris, parms);

  MRISstoreCurrentPositions(mris);
  delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, parms->tol, parms->n_averages);

  sse = MRIScomputeSSE(mris, parms);

  if (sse > starting_sse) /* error increased - turn off momentum */
  {
    mrisClearMomentum(mris);
    parms->dt *= parms->dt_decrease;
    if (parms->dt <= parms->base_dt) {
      parms->dt = parms->base_dt;
    }

    if (sse / starting_sse > parms->error_ratio) /* undo time step */
    {
      if (Gdiag & DIAG_SHOW)
        fprintf(stdout, "sse increased by %2.0f%%, undoing time step...\n", (float)sse / starting_sse * 100.0f);
      if (parms->dt > parms->base_dt) /* reset step size */
      {
        parms->dt = parms->base_dt;
      }

      /* undo the time step */
      MRISrestoreOldPositions(mris);
      mrisProjectSurface(mris);
    }
  }
  else /* error decreased */
  {
    parms->dt *= parms->dt_increase;
  }

  return (delta_t);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static void mrisAsynchronousTimeStep_update_odxyz(
    MRI_SURFACE * const mris, 
    float         const momentum, 
    float         const delta_t, 
    float         const max_mag,
    VERTEX *      const v)
{
    v->odx = delta_t * v->dx + momentum * v->odx;
    v->ody = delta_t * v->dy + momentum * v->ody;
    v->odz = delta_t * v->dz + momentum * v->odz;

    double mag = sqrt(v->odx * v->odx + v->ody * v->ody + v->odz * v->odz);
    
    if (mag > max_mag) {        /* don't let step get too big */
        mag = max_mag / mag;
        v->odx *= mag;
        v->ody *= mag;
        v->odz *= mag;
    }
}


static bool mrisAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex(    // returns false if tries to move outside available range
    MRI_SURFACE * const mris, 
    MHT *         const mht, 
    bool          const updateDxDyDz,
    int           const vno,
    MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context* ctx) 
{
  if (vno == Gdiag_no) {
      DiagBreak();
  }

  //VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];

  /* erase the faces this vertex is part of */

  // This will be a challenge to parallelize
  //
  if (mht) {
    MHTremoveAllFaces(mht, mris, vno);
  }

  bool canMove = true;
  if (mht) {
    canMove = mrisLimitGradientDistance(mris, mht, vno, ctx);
  }

  if (canMove) {

    MRISsetXYZ(mris, vno,
      v->x + v->odx,
      v->y + v->ody,
      v->z + v->odz);

    if ((fabs(v->x) > 128.0f) || (fabs(v->y) > 128.0f) || (fabs(v->z) > 128.0f)) {
      DiagBreak();
    }

    // If this is done here, then the undo step needs to cope 
    // Perhaps this is why mrisAsynchronousTimeStepNew does not do it
    //
    if (updateDxDyDz) {
      v->dx = v->odx; /* for mrisTrackTotalDistances */
      v->dy = v->ody;
      v->dz = v->odz;
    }
  }

  if (mht) {
    MHTaddAllFaces(mht, mris, vno);
  }

  return canMove;
}

static void mrisAsynchronousTimeStep_optionalDxDyDzUpdate( // BEVIN mris_make_surfaces 1
    MRI_SURFACE * const mris, 
    float         const momentum, 
    float         const delta_t, 
    MHT *         const mht, 
    float         const max_mag,
    int*          const directionPtr,
    bool          const updateDxDyDz)
{    
  *directionPtr *= -1;

  /* take a step in the gradient direction modulated by momentum */

  // Rigid bodies are simple
  //
  if (mris->status == MRIS_RIGID_BODY) {
    mris->da = delta_t * mris->alpha + momentum * mris->da;
    mris->db = delta_t * mris->beta  + momentum * mris->db;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg);
    return;
  }

  // The following algorithm moves each face in turn, but the movement is constrained by the placement of the other faces.
  // This makes it hard to parallelize, because deciding that two or more faces can move and then moving them could cause
  // collisions, each moving into the same previously vacant space, not realizing the other was moving into it.
  //
  // The parallelization is made harder by the shared MHT cache.
  //
  // The direction state stops bias in one direction during the expansion, but really indicates that any order of the 
  // vertices is acceptable.  Indeed random may be preferable.
  //
  // Parallelization is done by 
  //    partitioning the volume into subvolumes, 
  //    moving those vertices whose faces don't cross a subvolume wall
  //        the subvolumes can thus be done in parallel
  //            and, if they don't share portions of the mht, sharing in it is not a problem
  //    moving those vertices that cross the subvolume wall
  //        if there is only a few of them, this can be done serial, otherwise we could use a different partitioning
  //
#ifdef HAVE_OPENMP
  if (0)
#endif
  {

    // This is the origin, serial, algorithm
    //
    int i;
    for (i = 0; i < mris->nvertices; i++) {

            // BEVIN   <   WAS THE OLD COMPARISON
            // REVERSING IT TO SEE THE EFFECT
            //
      int const vno =                 // this strongly suggests it is NOT a parallelizable algorithm
        (*directionPtr > 0)           // since the direction through the nodes should not make a difference!
        ? (mris->nvertices - i - 1)
        : i;

      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        ROMP_PF_continue;
      }

      mrisAsynchronousTimeStep_update_odxyz(
        mris, momentum, delta_t, max_mag, v);
    
      mrisAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex(
        mris, mht, updateDxDyDz, vno, NULL);

    }
    
    return;
  } 
  
#ifdef HAVE_OPENMP

  // This is the parallel algorithm
  bool   const debug      = false;
  size_t const numThreads = omp_get_max_threads();

  // Calculate the bounding boxes
  //
  float xLo = 1e8, xHi=-xLo, yLo=xLo, yHi=xHi, zLo=xLo, zHi=yHi;
    // Box that will contains all the vertices after their largest possible movement

  typedef struct MrisAsynchronousTimeStep_optionalDxDyDzUpdate_PerVertexInfo_t PerVertexInfo;

  PerVertexInfo* vertexInfos = (PerVertexInfo*)calloc(mris->nvertices, sizeof(PerVertexInfo));

    // Compute the limits of each vertex's possible movement
    { 
      // Older versions of GCC don't support max and min reductions
      // so this codes it explicitly
      //
      float xLos[_MAX_FS_THREADS], 
            xHis[_MAX_FS_THREADS],
	    yLos[_MAX_FS_THREADS],
	    yHis[_MAX_FS_THREADS],
	    zLos[_MAX_FS_THREADS],
	    zHis[_MAX_FS_THREADS];
     
      { int i;
        for (i = 0; i < _MAX_FS_THREADS; i++) {
	  xLos[i] = xLo; xHis[i] = xHi;
	  yLos[i] = yLo; yHis[i] = yHi;
	  zLos[i] = zLo; zHis[i] = zHi;
	}
      }
      
      if (debugNonDeterminism) {
        fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
        mris_print_hash(stdout, mris, "mris ", "\n");
      }

      int vnoStep = (mris->nvertices + _MAX_FS_THREADS - 1)/_MAX_FS_THREADS;
      int vnoLo; 
      ROMP_PF_begin
#if defined(HAVE_OPENMP)
      #pragma omp parallel for if_ROMP(assume_reproducible) /* reduction(max:xHi,yHi,zHi) reduction(min:xLo,yLo,zLo) */
#endif
      for (vnoLo = 0; vnoLo < mris->nvertices; vnoLo += vnoStep) {
        ROMP_PFLB_begin
	int const tid =
#if defined(HAVE_OPENMP)
            omp_get_thread_num();
#else
	    0;
#endif	
	int vnoHi = MIN(mris->nvertices, vnoLo + vnoStep);
        int vno;
	for (vno = vnoLo; vno < vnoHi; vno++) { 
          VERTEX * const v = &mris->vertices[vno];

          mrisAsynchronousTimeStep_update_odxyz(
            mris, momentum, delta_t, max_mag, v);

          PerVertexInfo* pvi = vertexInfos + vno;
          float tLo, tHi;
          tLo = v->x; tHi = v->x + v->odx; if (tLo > tHi) { float temp = tLo; tLo = tHi; tHi = temp; }
          xLos[tid] = MIN(xLos[tid], pvi->xLo = tLo); xHis[tid] = MAX(xHis[tid], pvi->xHi = tHi);
          tLo = v->y; tHi = v->y + v->ody; if (tLo > tHi) { float temp = tLo; tLo = tHi; tHi = temp; }
          yLos[tid] = MIN(yLos[tid], pvi->yLo = tLo); yHis[tid] = MAX(yHis[tid], pvi->yHi = tHi);
          tLo = v->z; tHi = v->z + v->odz; if (tLo > tHi) { float temp = tLo; tLo = tHi; tHi = temp; }
          zLos[tid] = MIN(zLos[tid], pvi->zLo = tLo); zHis[tid] = MAX(zHis[tid], pvi->zHi = tHi);
	}
	
        ROMP_PFLB_end
      }
      ROMP_PF_end

      { int i;
        for (i = 0; i < _MAX_FS_THREADS; i++) {
	  xLo = MIN(xLo,xLos[i]); xHi = MAX(xHi, xHis[i]);
	  yLo = MIN(yLo,yLos[i]); yHi = MAX(yHi, yHis[i]);
	  zLo = MIN(zLo,zLos[i]); zHi = MAX(zHi, zHis[i]);
	}
      }
    }

  if (debug) {
    fprintf(stderr, "%s:%d vertex limits x:%g..%g y:%g..%g z:%g..%g\n", __FILE__, __LINE__, xLo,xHi,yLo,yHi,zLo,zHi);
  }
      
  // Create the subvolumes
  //
  const size_t numSubvolsPerEdge = 4;
  const size_t numSubvolsPerThread = (numSubvolsPerEdge * numSubvolsPerEdge * numSubvolsPerEdge + 1);
  // Info that must be shared with called functions
  //  
  MRISAsynchronousTimeStep_optionalDxDyDzUpdate_Context allVertexsContext;
  bzero(&allVertexsContext, sizeof(allVertexsContext));
  allVertexsContext.xLo = xLo; allVertexsContext.xHi = xHi; 
  allVertexsContext.yLo = yLo; allVertexsContext.yHi = yHi; 
  allVertexsContext.zLo = zLo; allVertexsContext.zHi = zHi;

  // Decide the subvolume geometry
  //
  {
    float const xSubvolLen = (xHi - xLo) / numSubvolsPerEdge; float xSubvolVerge = xSubvolLen*0.02;
    float const ySubvolLen = (yHi - yLo) / numSubvolsPerEdge; float ySubvolVerge = ySubvolLen*0.02;
    float const zSubvolLen = (zHi - zLo) / numSubvolsPerEdge; float zSubvolVerge = zSubvolLen*0.02;
      // The verge is a safety margin.  
      // Any face that gets this close to the surface of the subvolume is assumed might influence across it due to rounding errors.

    float const min_nbr_dist = minNeighborDistance(mris);
    {
      // Even things this close to the subvolume surface can affect the adjacent subvolumes
      //
      xSubvolVerge += min_nbr_dist;
      ySubvolVerge += min_nbr_dist;
      zSubvolVerge += min_nbr_dist;
    }
    
    allVertexsContext.xSubvolLen   = xSubvolLen;
    allVertexsContext.ySubvolLen   = ySubvolLen;
    allVertexsContext.zSubvolLen   = zSubvolLen;
    allVertexsContext.xSubvolVerge = xSubvolVerge;
    allVertexsContext.ySubvolVerge = ySubvolVerge;
    allVertexsContext.zSubvolVerge = zSubvolVerge;
  }

  // Assign the faces to the subvolumes
  //
  typedef struct PerFaceInfo_t {
    int svi;            // the subvolume this face is assigned to
  } PerFaceInfo;
  PerFaceInfo* faceInfos = (PerFaceInfo*)calloc(mris->nfaces, sizeof(PerFaceInfo));

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
  }

  { int fno;
    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(assume_reproducible)
    for (fno = 0; fno < mris->nfaces; fno++) {
      ROMP_PFLB_begin
      
      FACE const * const face = &mris->faces[fno];

      // Compute the box that will contain the face as each vertex moves through its possible values
      //
      int fi = 0;
      int vno = face->v[fi];
      PerVertexInfo* pvi = vertexInfos + vno; 

      float fxLo = pvi->xLo, fxHi = pvi->xHi, fyLo = pvi->yLo, fyHi = pvi->yHi, fzLo = pvi->zLo, fzHi = pvi->zHi;

      if (debug && fno < 10) 
      #pragma omp critical
      {
        fprintf(stderr, "%s:%d fno:%d has vertices\n", __FILE__, __LINE__, fno);
        fprintf(stderr, " vno:%d x:%g..%g y:%g..%g z:%g..%g\n", vno, pvi->xLo,pvi->xHi,pvi->yLo,pvi->yHi,pvi->zLo,pvi->zHi);
      }

      for (fi = 1; fi < VERTICES_PER_FACE; fi++) {
        vno = face->v[fi];
        pvi = vertexInfos + vno; 
        if (debug && fno < 10) 
        #pragma omp critical
        {
          fprintf(stderr, " vno:%d x:%g..%g y:%g..%g z:%g..%g\n", vno, pvi->xLo,pvi->xHi,pvi->yLo,pvi->yHi,pvi->zLo,pvi->zHi);
        }
        fxLo = MIN(fxLo,pvi->xLo); fyLo = MIN(fyLo,pvi->yLo); fzLo = MIN(fzLo,pvi->zLo);
        fxHi = MAX(fxHi,pvi->xHi); fyHi = MAX(fyHi,pvi->yHi); fzHi = MAX(fzHi,pvi->zHi);
      }

      // Compute the subvolumes that this face intersects
      //
      int const sviLo = MRISAsynchronousTimeStep_optionalDxDyDzUpdate_svi(&allVertexsContext, fxLo, fyLo, fzLo);
      int const sviHi = MRISAsynchronousTimeStep_optionalDxDyDzUpdate_svi(&allVertexsContext, fxHi, fyHi, fzHi);
      
      // Choose the subvolume to put this face into
      //
      int const svi = (sviLo==sviHi) ? sviLo : numSubvolsPerThread-1;   // if not same subvol, into the shared subvol

      // Assign the face to the subvolume
      //
      PerFaceInfo* pfi = faceInfos + fno;
      pfi->svi = svi;
      
      if (debug && fno < 10) 
      #pragma omp critical
      {
        fprintf(stderr, "%s:%d fno:%d assigned to svi:%d\n", __FILE__, __LINE__, fno, svi);
      }
      
      ROMP_PFLB_end
    }
    ROMP_PF_end
  }

  // In parallel, assign each vertex to either a subvolume or to the cross-subvolumes subvolume
  //
  const size_t numSubvols = numSubvolsPerThread * numThreads;
    // allocated per thread to avoid need to lock below

  typedef struct SubvolInfo_t {
    int firstVnoPlus1;   // plus1 so that 0 can act as NULL
    int lastVnoPlus1;    // plus1 so that 0 can act as NULL
  } SubvolInfo;
  SubvolInfo* subvols = (SubvolInfo*)calloc(numSubvols, sizeof(SubvolInfo));
  
  int* vnoToSvi = (int*)calloc(mris->nvertices, sizeof(int));
  allVertexsContext.vnoToSvi = vnoToSvi;
  
  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
  }

  { int i;
    ROMP_PF_begin
    #pragma omp parallel for if_ROMP(assume_reproducible)
    for (i = 0; i < mris->nvertices; i++) {
      ROMP_PFLB_begin
      
      int const vno =                       // I think this is to avoid some bias
        (*directionPtr < 0)                 // while maintaining some semblance of good cache behavior
        ? (mris->nvertices - i - 1)
        : i;
      
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];
      if (v->ripflag || vt->num == 0) {
        vnoToSvi[vno] = -1;
        ROMP_PF_continue;
      }
      
      // Choose the same svi as all its faces have, 
      // but if they disagree, use the shared subvol
      //
      int fi  = 0;
      int fno = vt->f[fi];
      PerFaceInfo* pfi = faceInfos + fno;
      int svi = pfi->svi;
      
      for (fi = 1; fi < vt->num; fi++) {
        int fno = vt->f[fi];
        pfi = faceInfos + fno;
        if (svi != pfi->svi) { svi = numSubvolsPerThread - 1; break; }
      }

      // Place into the subvol
      //
      const int tid = omp_get_thread_num();
      SubvolInfo* subvol = subvols + tid*numSubvolsPerThread + svi;
      PerVertexInfo* pvi = vertexInfos + vno; 

      pvi->nextVnoPlus1 = subvol->firstVnoPlus1;    
      subvol->firstVnoPlus1 = vno + 1;
      if (!subvol->lastVnoPlus1) subvol->lastVnoPlus1 = vno + 1;
  
      vnoToSvi[vno] = svi;
          
      ROMP_PFLB_end
    }
    ROMP_PF_end 
  }

  // Merge the per thread subvolumes into the tid==0 subvolume
  // It is important that the order be independent of the number of threads
  // so sort each subvols list...
  //
  { int * temp = (int*)malloc(sizeof(int)*mris->nvertices);

    for (unsigned int svi = 0; svi < numSubvolsPerThread; svi++) {
    
      // build the list in the temp
      //
      int tempSize = 0;

      for (unsigned int tid = 0; tid < numThreads; tid++) {
        SubvolInfo* subvol = subvols + tid*numSubvolsPerThread + svi;
        int vno = subvol->firstVnoPlus1 - 1;
        while (vno >= 0) {
          temp[tempSize++] = vno;
          PerVertexInfo* pvi = vertexInfos + vno;
          vno = pvi->nextVnoPlus1 - 1;
        }
      }
    
      // sort it, so it is thread count independent
      //
#if 0
      qsort(temp, tempSize, sizeof(int), int_compare);
#else
      sort_int(temp, tempSize, true);
#endif

      if (false && debugNonDeterminism) {
        unsigned long hash = fnv_add(fnv_init(), (const unsigned char*)temp, tempSize*sizeof(int));
        fprintf(stdout, "%s:%d stdout hash of the allocation into the %d subvolume:%ld\n",__FILE__,__LINE__,svi,hash);
      }


      // put it into the thread 0 subvolume for this svi
      //
      SubvolInfo* subvol0 = subvols + svi;
      subvol0->firstVnoPlus1 = 0;
      subvol0->lastVnoPlus1  = 0;
      if (tempSize > 0) {
        subvol0->firstVnoPlus1 = temp[         0] + 1;
        subvol0->lastVnoPlus1  = temp[tempSize-1] + 1;
        int i;
        for (i = 0; i < tempSize-1; i++) {
          int vno     = temp[i];
          int nextVno = temp[i+1];
          PerVertexInfo* pvi = vertexInfos + vno;
          pvi->nextVnoPlus1  = nextVno + 1;
        }
        int vno = temp[tempSize-1];
        PerVertexInfo* pvi = vertexInfos + vno;
        pvi->nextVnoPlus1  = 0;
      }
    }
  }

  // Make more info available to functions called below
  //  
  // Pass 0: In parallel, process each subvolume
  // Pass 1: In serial, process the cross-subvolume (parallel but only one hence serial)
  { 
    MHT_maybeParallel_begin();
    
    int pass;
    for (pass=0; pass<2; pass++) {
      int const sviLo = (pass==0) ? 0                     : numSubvolsPerThread-1;
      int const sviHi = (pass==0) ? numSubvolsPerThread-1 : numSubvolsPerThread;
      int svi;

      allVertexsContext.vertexInfos = (pass == 1) ? NULL : vertexInfos;   // on the second pass, the vertexs can move anywhere

      if (debugNonDeterminism) {
        fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);             // the results differ here on the 2nd pass!
        mris_print_hash(stdout, mris, "mris ", "\n");
      }
      ROMP_PF_begin
      #pragma omp parallel for if_ROMP(assume_reproducible)
      for (svi = sviLo; svi < sviHi; svi++) {
        ROMP_PFLB_begin
        MRISAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex_Context ctx;
        ctx.allVertexsContext = &allVertexsContext;
        ctx.hash  = fnv_init();
        ctx.count = 0;
        ctx.limit = 3035;
        ctx.trace = false;
        ctx.svi   = svi;
        SubvolInfo* subvol = subvols + svi;
        int vno = subvol->firstVnoPlus1 - 1;
        while (vno >= 0) {
          int vnoToRetry = -1;
          if (!mrisAsynchronousTimeStep_optionalDxDyDzUpdate_oneVertex(
            mris, mht, updateDxDyDz, vno,
            &ctx)) {
            vnoToRetry = vno;
          }
          PerVertexInfo* pvi = vertexInfos + vno;
          vno = pvi->nextVnoPlus1 - 1;
          if (vnoToRetry >= 0) 
#ifdef HAVE_OPENMP
          #pragma omp critical
#endif
          {
            if (pass == 1) *(int*)-1 = 0;                     // on the second pass, the vertexs can move anywhere, so this should not happen
            int shared_svi = numSubvolsPerThread - 1;
            SubvolInfo* subvol = subvols + shared_svi;
            pvi->nextVnoPlus1 = subvol->firstVnoPlus1;    
            subvol->firstVnoPlus1 = vno + 1;
            if (!subvol->lastVnoPlus1) subvol->lastVnoPlus1 = vno + 1;
		  }
        }
        ROMP_PFLB_end
      }
      ROMP_PF_end
    }

    MHT_maybeParallel_end();

  }

  // Free the temporary data
  free(vnoToSvi);
  free(faceInfos);
  free(subvols);
  free(vertexInfos);
#endif
}


// These functions were almost identical, but used separate "static int direction"
// Before parallelizing the above, I have merged them
//
double mrisAsynchronousTimeStep(
    MRI_SURFACE * const mris, 
    float         const momentum, 
    float         const delta_t, 
    MHT *         const mht, 
    float         const max_mag)
{
    MRISfreeDistsButNotOrig(mris);

    static int direction = -1;
    mrisAsynchronousTimeStep_optionalDxDyDzUpdate(
        mris, momentum, delta_t, mht, max_mag, &direction, true);

    return delta_t;
}


double mrisAsynchronousTimeStepNew(MRI_SURFACE *mris, float momentum, float delta_t, MHT *mht, float max_mag)
{
    MRISfreeDistsButNotOrig(mris);

    static int direction = -1;
    mrisAsynchronousTimeStep_optionalDxDyDzUpdate(
        mris, momentum, delta_t, mht, max_mag, &direction, false);

    return delta_t;
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double MRISmomentumTimeStep(MRI_SURFACE *mris, float momentum, float dt, float tol, float n_averages)
{
  double const delta_t = dt * sqrt((double)n_averages + 1.0);

  /* take a step in the gradient direction modulated by momentum */

  // Rigid bodies are simple
  //
  if (mris->status == MRIS_RIGID_BODY) {
    mris->da = delta_t * mris->alpha + momentum * mris->da;
    mris->db = delta_t * mris->beta  + momentum * mris->db;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg);
    return (delta_t);
  };
  
  if (mris->status == MRIS_SPHERICAL_PATCH) {
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      v->odx = delta_t * v->dx + momentum * v->odx;
      v->ody = delta_t * v->dy + momentum * v->ody;
      v->odz = delta_t * v->dz + momentum * v->odz;
      double mag = sqrt(v->odx * v->odx + v->ody * v->ody + v->odz * v->odz);
      if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
      {
        mag = MAX_MOMENTUM_MM / mag;
        v->odx *= mag;
        v->ody *= mag;
        v->odz *= mag;
      }
    }
    if (mris->patch != 2) {
      mrisApplyTopologyPreservingGradient(mris, 0, 1);
    }
  }
  else
  {
    MRISfreeDistsButNotOrig(mris);
        // MRISsetXYZ will invalidate all of these,
        // so make sure they are recomputed before being used again!
    
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      v->odx = delta_t * v->dx + momentum * v->odx;
      v->ody = delta_t * v->dy + momentum * v->ody;
      v->odz = delta_t * v->dz + momentum * v->odz;
      double mag = sqrt(v->odx * v->odx + v->ody * v->ody + v->odz * v->odz);
      if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
      {
        mag = MAX_MOMENTUM_MM / mag;
        v->odx *= mag;
        v->ody *= mag;
        v->odz *= mag;
      }
      if (vno == Gdiag_no) {
        float dx = v->x - v->origx;
        float dy = v->y - v->origy;
        float dz = v->z - v->origz;
        float dist = sqrt(dx * dx + dy * dy + dz * dz);
        float dot = dx * v->nx + dy * v->ny + dz * v->nz;
        fprintf(stdout,
                "moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
                "(%2.1f, %2.1f, %2.1f)\n",
                vno,
                v->odx,
                v->ody,
                v->odz,
                v->odx * v->nx + v->ody * v->ny + v->odz * v->nz,
                v->x,
                v->y,
                v->z);
        fprintf(
            stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n", v->nx, v->ny, v->nz, dist, dot);
      }
      
      MRISsetXYZ(mris, vno,
        v->x + v->odx,
        v->y + v->ody,
        v->z + v->odz);
    }
  }
  
  if (mris->status != MRIS_SPHERICAL_PATCH) {
    mrisProjectSurface(mris);
  }

#if DEBUG_HOMEOMORPHISM
  if (mris->patch == 2) /* bad trick  for debugging  */
  {
    mrisProjectSurface(mris);
  }
#endif

  return (delta_t);
}


MRI *MRISsolveLaplaceEquation(MRI_SURFACE *mris, MRI *mri, double res)
{
  MRI *mri_white, *mri_pial, *mri_laplace, *mri_control, *mri_tmp = NULL;
  int x, y, z, ncontrol, nribbon, v, i, xm1, xp1, ym1, yp1, zm1, zp1;
  VOXLIST *vl;
  float wval, pval, max_change, change, val, oval;

  MRISrestoreVertexPositions(mris, PIAL_VERTICES);
  mri_pial = MRISfillInterior(mris, res, NULL);
  mri_white = MRIclone(mri_pial, NULL);
  MRISrestoreVertexPositions(mris, WHITE_VERTICES);
  MRISfillInterior(mris, res, mri_white);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "pi.%2.2f.mgz", res);
    MRIwrite(mri_pial, fname);
    sprintf(fname, "wi.%2.2f.mgz", res);
    MRIwrite(mri_white, fname);
  }
  mri_laplace = MRIcloneDifferentType(mri_white, MRI_FLOAT);
  mri_control = MRIcloneDifferentType(mri_white, MRI_UCHAR);
  ncontrol = nribbon = 0;
  for (x = 0; x < mri_white->width; x++)
    for (y = 0; y < mri_white->height; y++)
      for (z = 0; z < mri_white->depth; z++) {
        wval = MRIgetVoxVal(mri_white, x, y, z, 0);
        pval = MRIgetVoxVal(mri_pial, x, y, z, 0);
        if (wval) {
          MRIsetVoxVal(mri_control, x, y, z, 0, CONTROL_MARKED);
          MRIsetVoxVal(mri_laplace, x, y, z, 0, 0.0);
          ncontrol++;
        }
        else if (FZERO(pval))  // outside pial surface
        {
          MRIsetVoxVal(mri_control, x, y, z, 0, CONTROL_MARKED);
          MRIsetVoxVal(mri_laplace, x, y, z, 0, 1.0);
          ncontrol++;
        }
        else
          nribbon++;
      }

  vl = VLSTalloc(nribbon);
  vl->mri = mri_laplace;
  nribbon = 0;
  for (x = 0; x < mri_white->width; x++)
    for (y = 0; y < mri_white->height; y++)
      for (z = 0; z < mri_white->depth; z++) {
        wval = MRIgetVoxVal(mri_white, x, y, z, 0);
        pval = MRIgetVoxVal(mri_pial, x, y, z, 0);
        if (FZERO(MRIgetVoxVal(mri_control, x, y, z, 0))) {
          vl->xi[nribbon] = x;
          vl->yi[nribbon] = y;
          vl->zi[nribbon] = z;
          nribbon++;
        }
      }

  i = 0;
  do {
    max_change = 0.0;
    mri_tmp = MRIcopy(mri_laplace, mri_tmp);
    for (v = 0; v < vl->nvox; v++) {
      x = vl->xi[v];
      y = vl->yi[v];
      z = vl->zi[v];
      xm1 = mri_laplace->xi[x - 1];
      xp1 = mri_laplace->xi[x + 1];
      ym1 = mri_laplace->yi[y - 1];
      yp1 = mri_laplace->yi[y + 1];
      zm1 = mri_laplace->zi[z - 1];
      zp1 = mri_laplace->zi[z + 1];
      oval = MRIgetVoxVal(mri_laplace, x, y, z, 0);
      val = (MRIgetVoxVal(mri_laplace, xm1, y, z, 0) + MRIgetVoxVal(mri_laplace, xp1, y, z, 0) +
             MRIgetVoxVal(mri_laplace, x, ym1, z, 0) + MRIgetVoxVal(mri_laplace, x, yp1, z, 0) +
             MRIgetVoxVal(mri_laplace, x, y, zm1, 0) + MRIgetVoxVal(mri_laplace, x, y, zp1, 0)) *
            1.0 / 6.0;
      change = fabs(val - oval);
      if (change > max_change) max_change = change;
      MRIsetVoxVal(mri_tmp, x, y, z, 0, val);
    }
    MRIcopy(mri_tmp, mri_laplace);
    i++;
    if (i % 10 == 0) printf("iter %d complete, max change %f\n", i, max_change);
  } while (max_change > 1e-3);

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    MRIwrite(mri_white, "w.mgz");
    MRIwrite(mri_pial, "p.mgz");
  }
  {
    char fname[STRLEN];
    sprintf(fname, "laplace.%2.2f.mgz", mri_laplace->xsize);
    MRIwrite(mri_laplace, fname);
  }
  MRIfree(&mri_white);
  MRIfree(&mri_pial);
  VLSTfree(&vl);
  return (mri_laplace);
}

int MRISmeasureLaplaceStreamlines(MRI_SURFACE *mris, MRI *mri_laplace, MRI *mri_intensity, MRI *mri_profiles)
{
  MRI * mri_mag  = MRIclone(mri_laplace, NULL);
  MRI * mri_grad = MRIsobel(mri_laplace, NULL, mri_mag);

  MRISrestoreVertexPositions(mris, PIAL_VERTICES);

  // normalize the gradient to be unit vectors
  {
    int z,y,x;
    for (z = 0; z < mri_grad->depth; z++)
      for (y = 0; y < mri_grad->height; y++)
        for (x = 0; x < mri_grad->width; x++) {
          double norm = MRIgetVoxVal(mri_mag, x, y, z, 0);
          if (FZERO(norm)) continue;
          int f;
          for (f = 0; f < mri_grad->nframes; f++) {
            double val = MRIgetVoxVal(mri_grad, x, y, z, f);
            MRIsetVoxVal(mri_grad, x, y, z, f, val / norm);
          }
        }
  }

  double const dt = 0.1;  // dt is in voxels
     
  double const voxsize = (mri_laplace->xsize + mri_laplace->ysize + mri_laplace->zsize) / 3;
  
  MRI * mri_surf_lap_grad = NULL;
  
  int nmissing = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    char const * const cp = getenv("USE_LAPLACE_GRAD");
    if (cp) {
      if (vno == 0) {
        printf("reading label stat from surface overlay %s\n", cp);
        mri_surf_lap_grad = MRIread(cp);
        if (mri_surf_lap_grad == NULL) ErrorExit(ERROR_NOFILE, "could not read USE_LAPLACE_GRAD file %s", cp);
      }
    }

    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
    VERTEX          * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    double xv, yv, zv;
    // check to see if this location doesn't have the resolution to represent the pial surface
    MRISsurfaceRASToVoxel(mris, mri_laplace, v->pialx, v->pialy, v->pialz, &xv, &yv, &zv);
    
    double val;
    MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
    if (val < 0.5) nmissing++;

    MRISsurfaceRASToVoxel(mris, mri_laplace, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv);
    double dist = 0.0;
    int npoints = 0;
    do {
      double dx, dy, dz;
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 0, &dx);
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 1, &dy);
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 2, &dz);
      double norm = sqrt(dx * dx + dy * dy + dz * dz);
      npoints++;
      if (FZERO(norm) || dist > 10) {
        if (val < .9) DiagBreak();
        if (dist > 10) DiagBreak();
        break;
      }
      dx /= norm;
      dy /= norm;
      dz /= norm;
      xv += dx * dt;
      yv += dy * dt;
      zv += dz * dt;
      dist += dt * voxsize;
      MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
      if (vno == Gdiag_no) printf("v %d:   (%2.2f %2.2f %2.2f): dist=%2.2f, val=%2.2f\n", vno, xv, yv, zv, dist, val);
    } while (val < 1);

    int nbr_is_diag = 0;
    if (Gdiag_no >= 0) {
      int n;
      for (n = 0; n < vt->vtotal; n++)
        if (vt->v[n] == Gdiag_no) {
          nbr_is_diag = 1;
          break;
        }
    }
    
    if (vno == Gdiag_no || nbr_is_diag) {
      LABEL *area = LabelAlloc(npoints, NULL, NULL);
      MRISsurfaceRASToVoxel(mris, mri_laplace, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv);
      dist = 0.0;
      int i = 0;
      do {
        MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
        
        double xv2, yv2, zv2;
        MRIvoxelToVoxel(mri_laplace, mri_intensity, xv, yv, zv, &xv2, &yv2, &zv2);
        
        double ival;
        MRIsampleVolumeFrame(mri_intensity, xv2, yv2, zv2, 0, &ival);
        
        double xs,ys,zs;
        MRISsurfaceRASFromVoxel(mris, mri_intensity, xv2, yv2, zv2, &xs, &ys, &zs);
        
        area->lv[i].x = xs;
        area->lv[i].y = ys;
        area->lv[i].z = zs;
        area->lv[i].stat = ival;
        if (mri_surf_lap_grad) area->lv[i].stat = MRIgetVoxVal(mri_surf_lap_grad, vno, 0, 0, 0);
        area->lv[i].vno = vno;
        area->n_points++;
        
        double dx,dy,dz;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 0, &dx);
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 1, &dy);
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 2, &dz);
        
        double norm = sqrt(dx * dx + dy * dy + dz * dz);
        if (FZERO(norm) || dist > 10) {
          if (val < .9) DiagBreak();
          if (dist > 10) DiagBreak();
          break;
        }
        dx /= norm;
        dy /= norm;
        dz /= norm;
        xv += dx * dt;
        yv += dy * dt;
        zv += dz * dt;
        MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
        dist += dt * voxsize;
        i++;
      } while (val < 1);

      char fname[STRLEN];
      sprintf(fname, "vno%d.label", vno);
      printf("writing label %s\n", fname);
      LabelWrite(area, fname);
      LabelFree(&area);
    }

    double const points_per_frame = (double)npoints / (double)mri_profiles->nframes;

    MRISsurfaceRASToVoxel(mris, mri_laplace, v->whitex, v->whitey, v->whitez, &xv, &yv, &zv);
    int i = 0, num = 0, frame = 0;
    dist = 0.0;
    do {
      MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
      
      double xv2, yv2, zv2;
      MRIvoxelToVoxel(mri_laplace, mri_intensity, xv, yv, zv, &xv2, &yv2, &zv2);
      
      double ival;
      MRIsampleVolumeFrame(mri_intensity, xv2, yv2, zv2, 0, &ival);
      
      double dx,dy,dz;
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 0, &dx);
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 1, &dy);
      MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 2, &dz);
      
      double norm = sqrt(dx * dx + dy * dy + dz * dz);
      double val = MRIgetVoxVal(mri_profiles, vno, 0, 0, frame);
      val = val + ival;
      if (FZERO(norm) || dist > 10) {
        if (num > 0) val /= num;
        MRIsetVoxVal(mri_profiles, vno, 0, 0, frame, val);
        if (val < .9) DiagBreak();
        if (dist > 10) DiagBreak();
        break;
      }
      dx /= norm;
      dy /= norm;
      dz /= norm;
      xv += dx * dt;
      yv += dy * dt;
      zv += dz * dt;
      i++;
      num++;
      int new_frame = (int)(i / points_per_frame);
      if (new_frame > frame) {
        if (num > 0) val /= num;
        num = 0;
        MRIsetVoxVal(mri_profiles, vno, 0, 0, frame, val);
        frame = new_frame;
        if (frame >= mri_profiles->nframes) break;
      }
      else
        MRIsetVoxVal(mri_profiles, vno, 0, 0, frame, val);

      MRIsampleVolumeFrame(mri_laplace, xv, yv, zv, 0, &val);
      dist += dt * voxsize;
    } while (val < 1);

    v->curv = dist;
  }

  //  printf("%d of %d pial surface nodes not resolved - %2.3f %%\n",
  //         nmissing, mris->nvertices, 100.0*nmissing/mris->nvertices) ;

  if (mri_surf_lap_grad) MRIfree(&mri_surf_lap_grad);
  MRIfree(&mri_mag);
  MRIfree(&mri_grad);
  return (NO_ERROR);
}


/*
  \fn MRIS *MRIStessellate(MRI *mri,  int value, int all_flag)
  \brief Creates a surface from a volume by tiling the outside of all
  the voxels that match value. The all_flag should be set if
  remove_non_hippo_voxels() has been run in mri. This function was
  derived from mri_tessellate.c and produces the same exact result. .
  It calls the TESSxxx() functions below.  The code is pretty
  horrific. Don't blame me, I did not write it, I just copied the from
  mri_tessellate.
 */
MRIS *MRIStessellate(MRI *mri, int value, int all_flag)
{
  int imnr, i, j, f_pack, v_ind, f;
  int xnum, ynum, numimg;
  int face_index, vertex_index;
  int *face_index_table0;
  int *face_index_table1;
  tface_type *face;
  tvertex_type *vertex;
  int *vertex_index_table;
  int k, n;
  MATRIX *m;
  VECTOR *vw, *vv;
  VERTEX *vtx;
  int useRealRAS = 0;
  MRIS *surf;
  int which, fno, vertices[VERTICES_PER_FACE + 1];
  FACE *sface;

  xnum = mri->width;
  ynum = mri->height;
  numimg = mri->depth;
  face_index = 0;
  vertex_index = 0;

  face = (tface_type *)lcalloc(MAXFACES, sizeof(tface_type));

  face_index_table0 = (int *)lcalloc(600 * ynum * xnum, sizeof(int));
  face_index_table1 = (int *)lcalloc(600 * ynum * xnum, sizeof(int));

  vertex = (tvertex_type *)lcalloc(MAXVERTICES, sizeof(tvertex_type));
  vertex_index_table = (int *)lcalloc(800 * ynum * xnum, sizeof(int));

  for (imnr = 0; imnr <= numimg; imnr++) {
    if (Gdiag_no > 0) {
      if ((vertex_index || face_index) && !(imnr % 10))
        printf("slice %d: %d vertices, %d faces\n", imnr, vertex_index, face_index);
    }
    // i is for width
    for (i = 0; i <= ynum; i++)
      for (j = 0; j <= xnum; j++) {
        //              z, y,  x,     z,   y,   x
        if (TESSfacep(mri, imnr, i - 1, j - 1, imnr - 1, i - 1, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr, i - 1, j, imnr - 1, i - 1, j, value, all_flag) ||
            TESSfacep(mri, imnr, i, j, imnr - 1, i, j, value, all_flag) ||
            TESSfacep(mri, imnr, i, j - 1, imnr - 1, i, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr - 1, i, j - 1, imnr - 1, i - 1, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr - 1, i, j, imnr - 1, i - 1, j, value, all_flag) ||
            TESSfacep(mri, imnr, i, j, imnr, i - 1, j, value, all_flag) ||
            TESSfacep(mri, imnr, i, j - 1, imnr, i - 1, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr - 1, i - 1, j, imnr - 1, i - 1, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr - 1, i, j, imnr - 1, i, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr, i, j, imnr, i, j - 1, value, all_flag) ||
            TESSfacep(mri, imnr, i - 1, j, imnr, i - 1, j - 1, value, all_flag)) {
          v_ind = TESSaddVertex(mri, imnr, i, j, &vertex_index, vertex_index_table, vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j - 1,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        0,
                        2,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j,
                        imnr - 1,
                        i - 1,
                        j,
                        0,
                        3,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j,
                        imnr - 1,
                        i,
                        j,
                        0,
                        0,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j - 1,
                        imnr - 1,
                        i,
                        j - 1,
                        0,
                        1,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j - 1,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        2,
                        2,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j,
                        imnr - 1,
                        i - 1,
                        j,
                        2,
                        1,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j,
                        imnr,
                        i - 1,
                        j,
                        2,
                        0,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j - 1,
                        imnr,
                        i - 1,
                        j - 1,
                        2,
                        3,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        4,
                        2,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j,
                        imnr - 1,
                        i,
                        j - 1,
                        4,
                        3,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j,
                        imnr,
                        i,
                        j - 1,
                        4,
                        0,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j,
                        imnr,
                        i - 1,
                        j - 1,
                        4,
                        1,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);

          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        imnr,
                        i - 1,
                        j - 1,
                        1,
                        2,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j,
                        imnr,
                        i - 1,
                        j,
                        1,
                        1,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j,
                        imnr,
                        i,
                        j,
                        1,
                        0,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j - 1,
                        imnr,
                        i,
                        j - 1,
                        1,
                        3,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        imnr - 1,
                        i,
                        j - 1,
                        3,
                        2,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j,
                        imnr - 1,
                        i,
                        j,
                        3,
                        3,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j,
                        imnr,
                        i,
                        j,
                        3,
                        0,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j - 1,
                        imnr,
                        i,
                        j - 1,
                        3,
                        1,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i - 1,
                        j - 1,
                        imnr - 1,
                        i - 1,
                        j,
                        5,
                        2,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr - 1,
                        i,
                        j - 1,
                        imnr - 1,
                        i,
                        j,
                        5,
                        1,
                        v_ind,
                        1,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i,
                        j - 1,
                        imnr,
                        i,
                        j,
                        5,
                        0,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
          TESScheckFace(mri,
                        imnr,
                        i - 1,
                        j - 1,
                        imnr,
                        i - 1,
                        j,
                        5,
                        3,
                        v_ind,
                        0,
                        all_flag,
                        value,
                        face,
                        &face_index,
                        face_index_table0,
                        face_index_table1,
                        vertex);
        }
      }
    for (i = 0; i < ynum; i++) {
      for (j = 0; j < xnum; j++) {
        for (f = 0; f < 6; f++) {
          f_pack = f * ynum * xnum + i * xnum + j;
          face_index_table0[f_pack] = face_index_table1[f_pack];
        }
      }
    }
  }

  printf("MRIStessellate: nvertices = %d, nfaces = %d\n", vertex_index, face_index);
  surf = MRISalloc(vertex_index, 2 * face_index);
  surf->type = MRIS_TRIANGULAR_SURFACE;

  if (useRealRAS == 1)
    m = extract_i_to_r(mri);
  else
    m = surfaceRASFromVoxel_(mri);
  if (Gdiag_no > 0) {
    printf("MRIStessellate: using vox2ras matrix:\n");
    MatrixPrint(stdout, m);
  }

  vv = VectorAlloc(4, MATRIX_REAL);
  vv->rptr[4][1] = 1;
  vw = VectorAlloc(4, MATRIX_REAL);
  for (k = 0; k < vertex_index; k++) {
    vtx = &surf->vertices[k];

    // V4_LOAD(vv, vertex[k].j-0.5, vertex[k].i-0.5, vertex[k].imnr-0.5, 1);
    vv->rptr[1][1] = vertex[k].j - 0.5;
    vv->rptr[2][1] = vertex[k].i - 0.5;
    vv->rptr[3][1] = vertex[k].imnr - 0.5;

    MatrixMultiply(m, vv, vw);
    // we are doing the same thing as the following, but we save time in
    // calculating the matrix at every point
    // if (useRealRAS == 1)  // use the physical RAS as the vertex point
    //   MRIvoxelToWorld(mri,
    //                   vertex[k].j-0.5,
    //                   vertex[k].i-0.5,
    //                   vertex[k].imnr-0.5,
    //                   &x, &y, &z);
    // else
    //   MRIvoxelToSurfaceRAS(mri,
    //                        vertex[k].j-0.5,
    //                        vertex[k].i-0.5,
    //                        vertex[k].imnr-0.5,
    //                        &x, &y, &z);
    
    MRISsetXYZ(surf, k,
      V3_X(vw),
      V3_Y(vw),
      V3_Z(vw));
  }

  k = -1;
  for (fno = 0; fno < surf->nfaces; fno += 2) {
    k++;

    /* quan
    drangular face */
    for (n = 0; n < 4; n++) vertices[n] = face[k].v[n];

    /* if we're going to be arbitrary, we might as well be really arbitrary */
    /*  NOTE: for this to work properly in the write, the first two
        vertices in the first face (EVEN and ODD) must be 0 and 1. */
    which = WHICH_FACE_SPLIT(vertices[0], vertices[1]);

    /* 1st triangle */
    if (EVEN(which)) {
      surf->faces[fno].v[0] = vertices[0];
      surf->faces[fno].v[1] = vertices[1];
      surf->faces[fno].v[2] = vertices[3];

      /* 2nd triangle */
      surf->faces[fno + 1].v[0] = vertices[2];
      surf->faces[fno + 1].v[1] = vertices[3];
      surf->faces[fno + 1].v[2] = vertices[1];
    }
    else {
      surf->faces[fno].v[0] = vertices[0];
      surf->faces[fno].v[1] = vertices[1];
      surf->faces[fno].v[2] = vertices[2];

      /* 2nd triangle */
      surf->faces[fno + 1].v[0] = vertices[0];
      surf->faces[fno + 1].v[1] = vertices[2];
      surf->faces[fno + 1].v[2] = vertices[3];
    }
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      surf->vertices_topology[surf->faces[fno    ].v[n]].num++;
      surf->vertices_topology[surf->faces[fno + 1].v[n]].num++;
    }
  }

  for (k = 0; k < vertex_index; k++) {
    surf->vertices_topology[k].f = (int   *)calloc(surf->vertices_topology[k].num, sizeof(int));
    surf->vertices_topology[k].n = (uchar *)calloc(surf->vertices_topology[k].num, sizeof(uchar));
    surf->vertices_topology[k].num = 0;
  }
  for (k = 0; k < surf->nfaces; k++) {
    sface = &(surf->faces[k]);
    for (n = 0; n < VERTICES_PER_FACE; n++) surf->vertices_topology[sface->v[n]].f[surf->vertices_topology[sface->v[n]].num++] = k;
  }

  for (k = 0; k < vertex_index; k++) {
    for (n = 0; n < surf->vertices_topology[k].num; n++) {
      for (i = 0; i < VERTICES_PER_FACE; i++) {
        if (surf->faces[surf->vertices_topology[k].f[n]].v[i] == k) surf->vertices_topology[k].n[n] = i;
      }
    }
  }

  mrisCheckVertexFaceTopology(surf);
  
  getVolGeom(mri, &(surf->vg));
  mrisCompleteTopology(surf);
  MRIScomputeMetricProperties(surf);

  return (surf);
}

void TESSaddFace(MRI *mri,
                 int imnr,
                 int i,
                 int j,
                 int f,
                 int prev_flag,
                 int *pface_index,
                 tface_type *face,
                 int *face_index_table0,
                 int *face_index_table1)
{
  int xnum, ynum, pack;

  xnum = mri->width;
  ynum = mri->height;
  pack = f * ynum * xnum + i * xnum + j;

  if (*pface_index >= MAXFACES - 1) ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded", Progname, MAXFACES);
  if (prev_flag) {
    face_index_table0[pack] = *pface_index;
  }
  else {
    face_index_table1[pack] = *pface_index;
  }
  face[*pface_index].imnr = imnr;  // z
  face[*pface_index].i = i;        // y
  face[*pface_index].j = j;        // x
  face[*pface_index].f = f;
  face[*pface_index].num = 0;
  (*pface_index)++;
}

void TESScheckFace(MRI *mri,
                   int im0,
                   int i0,
                   int j0,
                   int im1,
                   int i1,
                   int j1,
                   int f,
                   int n,
                   int v_ind,
                   int prev_flag,
                   int all_flag,
                   int value,
                   tface_type *face,
                   int *pface_index,
                   int *face_index_table0,
                   int *face_index_table1,
                   tvertex_type *vertex)
{
  int xnum, ynum, numimg, f_pack, f_ind;
  int imax, imin, jmax, jmin;
  xnum = mri->width;
  ynum = mri->height;
  numimg = mri->depth;
  f_pack = f * ynum * xnum + i0 * xnum + j0;  // f= 0, 1, 2, 3, 4, 5

  jmax = mri->width;
  jmin = 0;
  imax = mri->height;
  imin = 0;

  if ((im0 >= 0 && im0 < numimg && i0 >= imin && i0 < imax && j0 >= jmin && j0 < jmax && im1 >= 0 && im1 < numimg &&
       i1 >= imin && i1 < imax && j1 >= jmin && j1 < jmax)) {
    if ((all_flag && ((MRIgetVoxVal(mri, j0, i0, im0, 0) !=
                       MRIgetVoxVal(mri, j1, i1, im1, 0)) /* && (MRIvox(mri, j1, i1, im1) == 0)*/)) ||
        (((MRIgetVoxVal(mri, j0, i0, im0, 0) == value) && (MRIgetVoxVal(mri, j1, i1, im1, 0) != value)))) {
      if (n == 0) TESSaddFace(mri, im0, i0, j0, f, prev_flag, pface_index, face, face_index_table0, face_index_table1);
      if (prev_flag)
        f_ind = face_index_table0[f_pack];
      else
        f_ind = face_index_table1[f_pack];
      face[f_ind].v[n] = v_ind;
      if (vertex[v_ind].num < 9) vertex[v_ind].f[vertex[v_ind].num++] = f_ind;
    }
  }
}

int TESSaddVertex(MRI *mri, int imnr, int i, int j, int *pvertex_index, int *vertex_index_table, tvertex_type *vertex)
{
  int xnum = mri->width;

  int pack = i * (xnum + 1) + j;

  if (*pvertex_index >= MAXVERTICES - 1)
    ErrorExit(ERROR_NOMEMORY, "%s: max vertices %d exceeded", Progname, MAXVERTICES);
  vertex_index_table[pack] = *pvertex_index;
  vertex[*pvertex_index].imnr = imnr;  // z
  vertex[*pvertex_index].i = i;        // y
  vertex[*pvertex_index].j = j;        // x
  vertex[*pvertex_index].num = 0;

  return ((*pvertex_index)++);
}

int TESSfacep(MRI *mri, int im0, int i0, int j0, int im1, int i1, int j1, int value, int all_flag)
{
  int numimg, imax, imin, jmax, jmin;
  numimg = mri->depth;
  // it is so confusing this guy uses j for width and i for height
  jmax = mri->width;
  jmin = 0;
  imax = mri->height;
  imin = 0;
  return (im0 >= 0 && im0 < numimg && i0 >= imin && i0 < imax && j0 >= jmin && j0 < jmax && im1 >= 0 && im1 < numimg &&
          i1 >= imin && i1 < imax && j1 >= jmin && j1 < jmax &&
          MRIgetVoxVal(mri, j0, i0, im0, 0) != MRIgetVoxVal(mri, j1, i1, im1, 0) &&
          ((MRIgetVoxVal(mri, j0, i0, im0, 0) == value || MRIgetVoxVal(mri, j1, i1, im1, 0) == value) || all_flag));
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
    incrementally stores distance from white surface in v->val field -
    useful for calling iteratively
  ------------------------------------------------------*/
#define MAX_EXP_MM 0.1
int MRISexpandSurface(MRI_SURFACE *mris, float distance, INTEGRATION_PARMS *parms, int use_thick, int nsurfaces)
{
  int vno, n, niter, avgs, nrounds, surf_no, orig_start_t = parms->start_t;
  VERTEX *v;
  double dist, dx = 0., dy = 0., dz = 0.0, dtotal, *pial_x, *pial_y, *pial_z, l_spring_orig;
  char fname[STRLEN], *cp;
  INTEGRATION_PARMS thick_parms;
  MRI_SURFACE *mris_ico;
  MHT *mht = NULL;

  l_spring_orig = parms->l_spring;
  if (Gdiag & DIAG_SHOW) {
    mrisLogIntegrationParms(stderr, mris, parms);
  }
  if (Gdiag & DIAG_WRITE) {
    mrisLogIntegrationParms(parms->fp, mris, parms);
    mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
  }
  mrisLogStatus(mris, parms, stdout, 0.0f, -1);

  pial_x = (double *)calloc(mris->nvertices, sizeof(double));
  pial_y = (double *)calloc(mris->nvertices, sizeof(double));
  pial_z = (double *)calloc(mris->nvertices, sizeof(double));

  const char *hemi = mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh";

  if (pial_x == NULL || pial_y == NULL || pial_z == NULL) {
    ErrorExit(ERROR_NOMEMORY, "MRISexpandSurface: could not allocaet %d element vertex array", mris->nvertices);
  }

  if (nsurfaces > 1) {
    if (parms->smooth_averages > 0) {
      MRISrestoreVertexPositions(mris, PIAL_VERTICES);
      MRISaverageVertexPositions(mris, parms->smooth_averages);
      MRISsaveVertexPositions(mris, PIAL_VERTICES);
      printf("smoothing vertex locations %d times\n", parms->smooth_averages);
      MRISrestoreVertexPositions(mris, WHITE_VERTICES);
      MRISaverageVertexPositions(mris, parms->smooth_averages);
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
      MRISsaveVertexPositions(mris, WHITE_VERTICES);
    }
    sprintf(fname, "%s.%s%3.3d", hemi, parms->base_name, 0);
    printf("writing expanded surface to %s...\n", fname);
    MRISwrite(mris, fname);
  }

  if (use_thick) {
    MRISripZeroThicknessRegions(mris);
  }
  
  if (parms == NULL) {

    MRISfreeDistsButNotOrig(mris);
        // MRISsetXYZ will invalidate all of these,
        // so make sure they are recomputed before being used again!

    for (vno = 0; vno < mris->nvertices; vno++) {
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      MRISsetXYZ(mris, vno,
        v->x + distance * v->nx,
        v->y + distance * v->ny,
        v->z + distance * v->nz);
    }
    
  }
  else {
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
    if (use_thick) {
      memset(&thick_parms, 0, sizeof(thick_parms));
      thick_parms.dt = 0.2;
      thick_parms.momentum = .5;
      thick_parms.l_nlarea = 1;
      thick_parms.l_thick_min = 0;
      thick_parms.tol = 1e-1;
      thick_parms.l_thick_parallel = 1;
      thick_parms.remove_neg = 1;
      cp = getenv("FREESURFER_HOME");
      if (cp == NULL) {
        ErrorExit(ERROR_BADPARM, "%s: FREESURFER_HOME not defined in environment", cp);
      }
      sprintf(fname, "%s/lib/bem/ic7.tri", cp);
      mris_ico = MRISread(fname);
      if (!mris_ico) {
        ErrorExit(ERROR_NOFILE, "%s: could not open surface file %s", Progname, fname);
      }
      MRISscaleBrain(mris_ico, mris_ico, mris->radius / mris_ico->radius);
      MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES);
      MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
      MRISminimizeThicknessFunctional(mris, &thick_parms, 10.0);

      // compute locations of pial surface vertices nearest to each orig one
      for (vno = 0; vno < mris->nvertices; vno++) {
        v = &mris->vertices[vno];
        if (v->ripflag) {
          continue;
        }
        if (vno == Gdiag_no) {
          DiagBreak();
        }
        pial_x[vno] = v->tx;
        pial_y[vno] = v->ty;
        pial_z[vno] = v->tz;
      }
      MRISrestoreVertexPositions(mris, WHITE_VERTICES);
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
    }
    if ((parms->write_iterations > 0) && (Gdiag & DIAG_WRITE) && !parms->start_t) {
      mrisWriteSnapshot(mris, parms, 0);
    }
    mrisClearMomentum(mris);
    nrounds = log2(parms->n_averages) - log2(parms->min_averages + 1) + 2;
    printf("nrounds = %d\n", nrounds);
    if (use_thick)  // distance is a % of the total thickness
    {
      niter = nint(fabs(distance) * 2 / (nsurfaces * parms->dt * MAX_EXP_MM));
    }
    else {
      niter = nint(2 * fabs(distance) / (parms->l_location * nrounds * nsurfaces * parms->dt * MAX_EXP_MM));
    }
    if (Gdiag_no >= 0) {
      v = &mris->vertices[Gdiag_no];
      if (use_thick)
        printf("v %d: thickness=%2.2f, moving outwards %2.2fmm\n", Gdiag_no, v->curv, v->curv * distance);
      else {
        printf("v %d: moving outwards %2.2fmm\n", Gdiag_no, distance);
      }
    }

    for (surf_no = 0; surf_no < nsurfaces; surf_no++) {
      // compute target locations for each vertex
      ROMP_PF_begin
      // ifdef HAVE_OPENMP
      // pragma omp parallel for if_ROMP(experimental)
      // endif
      for (vno = 0; vno < mris->nvertices; vno++) {
        ROMP_PFLB_begin
	
        v = &mris->vertices[vno];
        if (v->ripflag) {
          ROMP_PF_continue;
        }
        if (vno == Gdiag_no) {
          DiagBreak();
        }
        if (use_thick) {
          dx = pial_x[vno] - v->origx;
          dy = pial_y[vno] - v->origy;
          dz = pial_z[vno] - v->origz;
          dtotal = sqrt(dx * dx + dy * dy + dz * dz);
          if (FZERO(dtotal)) {
            v->ripflag = 1;
          }
          dist = (surf_no + 1) * distance * dtotal / ((float)nsurfaces);
          dx /= dtotal;
          dy /= dtotal;
          dz /= dtotal;
          v->targx = v->origx + dx * dist;
          v->targy = v->origy + dy * dist;
          v->targz = v->origz + dz * dist;
        }
        else  // just move outwards along surface normal
        {
          dist = distance;
          v->targx = v->origx + v->nx * distance;
          v->targy = v->origy + v->ny * distance;
          v->targz = v->origz + v->nz * distance;
        }
        if (parms->mri_brain && parms->target_intensity >= 0 && 0) {
          double step_size, xw, yw, zw, xv, yv, zv, val, val0, d;
          MRI *mri;

          mri = parms->mri_brain;
          step_size = (mri->xsize + mri->ysize + mri->zsize) / 6.0;  // shannon

          MRISsurfaceRASToVoxelCached(mris, mri, v->origx, v->origy, v->origz, &xv, &yv, &zv);
          MRIsampleVolume(mri, xv, yv, zv, &val0);
          v->val2 = val0;
          for (d = 0; d < dist; d += step_size) {
            xw = v->origx + v->nx * d;
            yw = v->origy + v->ny * d;
            zw = v->origz + v->nz * d;

            MRISsurfaceRASToVoxelCached(mris, mri, xw, yw, zw, &xv, &yv, &zv);
            MRIsampleVolume(mri, xv, yv, zv, &val);
            if (vno == Gdiag_no) printf("v %d: dist %2.2f, (%2.1f, %2.1f, %2.1f) = %2.1f\n", vno, d, xv, yv, zv, val);
            if ((val0 > parms->target_intensity && val < parms->target_intensity) ||
                (val0 < parms->target_intensity && val > parms->target_intensity) ||
                MRIindexNotInVolume(mri, xv, yv, zv))
              break;
          }
          d -= STEP_SIZE;
          v->targx = v->origx + v->nx * d;
          v->targy = v->origy + v->ny * d;
          v->targz = v->origz + v->nz * d;
        }
	ROMP_PFLB_end
      }
      ROMP_PF_end
      
      if (Gdiag_no >= 0)
        printf("v%d, target - %2.1f %2.1f %2.1f\n",
               Gdiag_no,
               mris->vertices[Gdiag_no].targx,
               mris->vertices[Gdiag_no].targy,
               mris->vertices[Gdiag_no].targz);
      if (Gdiag & DIAG_WRITE) {
        char fname[STRLEN];
        int vno, max_vno;
        VERTEX *v;
        double dist, max_dist;

        max_dist = 0;
        max_vno = 0;
        for (vno = 0; vno < mris->nvertices; vno++) {
          v = &mris->vertices[vno];
          dist = sqrt(SQR(v->x - v->targx) + SQR(v->y - v->targy) + SQR(v->z - v->targz));
          v->d = dist;
          if (dist > max_dist) {
            max_dist = dist;
            max_vno = vno;
          }
        }

        sprintf(fname, "%s.target_dist.%d.mgz", parms->base_name, surf_no);
        printf("writing target surface distances to %s\n", fname);
        MRISwriteD(mris, fname);
        sprintf(fname, "%s.targets.%d", parms->base_name, surf_no);
        printf("writing target surface to %s\n", fname);
        MRISsaveVertexPositions(mris, TMP_VERTICES);
        MRISrestoreVertexPositions(mris, TARGET_VERTICES);
        MRISwrite(mris, fname);
        MRISrestoreVertexPositions(mris, TMP_VERTICES);
      }

      for (avgs = parms->n_averages; avgs >= parms->min_averages; avgs /= 2) {
        parms->l_spring = l_spring_orig * sqrt((float)parms->min_averages / (float)parms->n_averages);
        printf(
            "***************** integrating with averages = %d, niter = %d, l_spring = %2.3f ***********************\n",
            avgs,
            niter,
            parms->l_spring);
        for (n = parms->start_t; n < parms->start_t + niter; n++) {
          printf("\rstep %d of %d     ", n + 1 - (surf_no * niter * nrounds), orig_start_t + nrounds * niter);
          if (Gdiag_no >= 0) printf("\n");
          fflush(stdout);
          MRIScomputeMetricProperties(mris);
          if (!(parms->flags & IPFLAG_NO_SELF_INT_TEST)) {
            double vmean, vsigma;
            float voxel_res = (mris->vg.xsize + mris->vg.ysize + mris->vg.zsize) / 3;

            vmean = MRIScomputeTotalVertexSpacingStats(mris, &vsigma, NULL, NULL, NULL, NULL);
            if (FZERO(voxel_res)) voxel_res = 1;
            MHTfree(&mht); mht = MHTcreateFaceTable_Resolution(mris, CURRENT_VERTICES, vmean);
          }

          if (parms->mri_brain && parms->target_intensity >= 0)
            mrisUpdateTargetLocations(mris, parms->mri_brain, parms->target_intensity);
          //      MRISsmoothSurfaceNormals(mris, avgs) ;
          MRISclearGradient(mris);
          mrisComputeTargetLocationTerm(mris, parms->l_location, parms);
          {
            int vno;
            double intensity_dif;
            VERTEX *v;
            for (vno = 0; vno < mris->nvertices; vno++) {
              v = &mris->vertices[vno];
              if (v->marked == 0 || v->ripflag) continue;
              intensity_dif = abs(v->val - parms->target_intensity) / 10;
              if (intensity_dif > 1) {
                intensity_dif = MIN(intensity_dif, 10);
                if (vno == Gdiag_no)
                  printf(
                      "v %d, increasing weight of intensity term by %2.1f (val=%2.1f)\n", vno, intensity_dif, v->val);
                v->dx *= intensity_dif;
                v->dy *= intensity_dif;
                v->dz *= intensity_dif;
              }
            }
          }
          mrisComputeSpringTerm(mris, parms->l_spring);
          mrisComputeConvexityTerm(mris, parms->l_convex);
          mrisComputeLaplacianTerm(mris, parms->l_lap);
          mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm);
          mrisComputeThicknessSmoothnessTerm(mris, parms->l_tsmooth, parms);
          mrisComputeThicknessMinimizationTerm(mris, parms->l_thick_min, parms);
          mrisComputeThicknessParallelTerm(mris, parms->l_thick_parallel, parms);
          mrisComputeNormalSpringTerm(mris, parms->l_nspring);
          mrisComputeSurfaceNormalIntersectionTerm(mris, mht, parms->l_norm, 0.1);
          mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv);
          mrisComputeMaxSpringTerm(mris, parms->l_max_spring);
          mrisComputeNonlinearSpringTerm(mris, parms->l_nlspring, parms);
          mrisComputeTangentialSpringTerm(mris, parms->l_tspring);
          mrisComputeNonlinearTangentialSpringTerm(mris, parms->l_nltspring, parms->min_dist);
          mrisComputeAngleAreaTerms(mris, parms);
          MRISaverageGradients(mris, avgs);
          
          mrisAsynchronousTimeStep(mris, parms->momentum, parms->dt, mht, MAX_EXP_MM);

          if ((parms->write_iterations > 0) && !((n + 1) % parms->write_iterations) && (Gdiag & DIAG_WRITE)) {
            mrisWriteSnapshot(mris, parms, n + 1);
          }
        }
        parms->start_t += niter;
        if (avgs == 0) {
          break;
        }
      }
      if (parms->smooth_averages > 0) {
        printf("\nsmoothing vertex locations %d times\n", parms->smooth_averages);
        MRISaverageVertexPositions(mris, parms->smooth_averages);
      }
      else {
        printf("\n");
      }
      sprintf(fname, "%s.%s%3.3d", hemi, parms->base_name, surf_no + 1);
      if (nsurfaces > 1) {
        printf("writing expanded surface to %s...\n", fname);
        MRISwrite(mris, fname);
      }
    }
  }

  if (mht) {
    MHTfree(&mht);
  }
  free(pial_x);
  free(pial_y);
  free(pial_z);
  return (NO_ERROR);
}


static int mrisSmoothingTimeStep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

int MRISremoveOverlapWithSmoothing(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int negative, old_neg, same = 0, min_neg, min_neg_iter, last_expand;

  parms->dt = .99;
  parms->max_nbrs = 0;
  min_neg_iter = 0;
  last_expand = 0;
  parms->t = parms->start_t;
  min_neg = negative = MRIScountNegativeTriangles(mris);

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    if (!parms->fp) {
      sprintf(fname, "%s.%s.out", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);
      INTEGRATION_PARMS_openFp(parms, fname, "a");
      if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    }
  }

  printf("%03d: dt=%2.4f, %3d negative triangles  VmPeak %d\n", -1, 0.0, negative,GetVmPeak());
  while (negative > 0) {
    old_neg = negative;

    if (parms->fp && parms->t % 100 == 0)
      fprintf(parms->fp, "%03d: dt=%2.4f, %d negative triangles\n", parms->t, parms->dt, negative);

    printf("%03d: dt=%2.4f, %3d negative triangles\n", parms->t, parms->dt, negative);
    mrisSmoothingTimeStep(mris, parms);
    parms->t++;  // advance time-step counter
    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    negative = MRIScountNegativeTriangles(mris);
    if (negative < min_neg) {
      min_neg = negative;
      min_neg_iter = parms->t;
    }
    else if ((((parms->t - min_neg_iter) % 10) == 0) && parms->t > min_neg_iter) {
      if (parms->dt > 0.01) {
        parms->dt *= 0.95;
      }
    }
    else if ((parms->t > min_neg_iter + 50) && parms->t > last_expand + 25) {
      if (parms->max_nbrs < 1) {
        parms->dt = 0.99;
        parms->max_nbrs++;
        printf("expanding nbhd size to %d\n", parms->max_nbrs);
      }
      last_expand = parms->t;
      same = 0;
    }
    if (parms->t > min_neg_iter + 1000) {
      printf("terminating loop due to lack of progress\n");
      break;
    }
    if (old_neg == negative) {
      if (same++ > 25 && parms->max_nbrs < 1) {
        parms->max_nbrs++;
        parms->dt = 0.99;
        printf("expanding nbhd size to %d\n", parms->max_nbrs);
        last_expand = parms->t;
        same = 0;
      }
    }
    else {
      same = 0;
    }
    if (parms->t - parms->start_t > parms->niterations) {
      break;
    }

    fflush(stdout);
  }

  if (negative > 0) {
    printf("%03d: %d negative triangles\n", parms->t, negative);
  }
  if (parms->fp) {
    INTEGRATION_PARMS_closeFp(parms);
  }
  parms->t += parms->start_t;
  parms->start_t = parms->t;
  return (NO_ERROR);
}

static int mrisSmoothingTimeStep(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int vno, n, m, fno;
  FACE *face;
  double dx, dy, dz, x, y, z, max_dx, max_dy, max_dz;

  MRIScomputeMetricProperties(mris);
  MRISclearMarks(mris);
  for (fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->area < 0) {
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        VERTEX * const v = &mris->vertices[face->v[n]];
        v->area = -1;
        v->marked = 1;
      }
    }
  }

  MRISdilateMarked(mris, parms->max_nbrs);

  max_dx = max_dy = max_dz = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->marked == 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    x = v->x;
    y = v->y;
    z = v->z;

    dx = dy = dz = 0.0;
    n = 0;
    for (m = 0; m < vt->vnum; m++) {
      VERTEX const * vn = &mris->vertices[vt->v[m]];
      if (!vn->ripflag) {
        if (vno == Gdiag_no)
          fprintf(stdout,
                  "v %d --> %d spring term:         (%2.3f, %2.3f, %2.3f)\n",
                  vno,
                  vt->v[m],
                  vn->x - x,
                  vn->y - y,
                  vn->z - z);
        dx += vn->x - x;
        dy += vn->y - y;
        dz += vn->z - z;
        n++;
      }
    }
    if (n > 0) {
      dx = dx / n;
      dy = dy / n;
      dz = dz / n;
    }

    v->dx = dx;
    v->dy = dy;
    v->dz = dz;
    if (fabs(dx) > fabs(max_dx)) {
      max_dx = dx;
    }
    if (fabs(dy) > fabs(max_dy)) {
      max_dy = dy;
    }
    if (fabs(dz) > fabs(max_dz)) {
      max_dz = dz;
    }

    if (vno == Gdiag_no) fprintf(stdout, "v %d spring term:         (%2.3f, %2.3f, %2.3f)\n", vno, dx, dy, dz);
  }
  
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("max delta = (%2.4f, %2.4f, %2.4f) [%2.3f]\n",
           max_dx,
           max_dy,
           max_dz,
           sqrt(SQR(max_dx) + SQR(max_dy) + SQR(max_dz)));

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag || v->marked == 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    MRISsetXYZ(mris, vno,
      v->x + v->dx * parms->dt,
      v->y + v->dy * parms->dt,
      v->z + v->dz * parms->dt);
  }
  
  MRISclearMarks(mris);

  return (NO_ERROR);
}


int mrisScaleTimeStepByCurvature(MRI_SURFACE *mris)
{
  int vno;
  VERTEX *v;
  float scale;

  return (0);  // disabled
  MRIScomputeSecondFundamentalForm(mris);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    if (v->H < 1) {
      continue;
    }
    scale = 1.0 + (v->H - 1.0) / 10.0;
    if (vno == Gdiag_no) {
      printf("scaling vertex %d gradient by %2.2f\n", vno, scale);
    }
    v->odx *= scale;
    v->ody *= scale;
    v->odz *= scale;
  }

  return (NO_ERROR);
}


int MRISremoveCompressedRegions(MRI_SURFACE *mris, double min_dist)
{
  int compressed, iters = 0, old_compressed, n_averages, nsize;
  MHT *mht = NULL;
  double delta_t, l_spring, l_convex, l_max_spring;
  ;
  char fname[STRLEN];
  INTEGRATION_PARMS parms;

  {
    int fno;
    for (fno = 0; fno < mris->nfaces; fno++) {
      setFaceOrigArea(mris,fno,0.5f);
    }
  }

  memset(&parms, 0, sizeof(parms));

  parms.l_parea = .002;
  l_spring = .01;
  l_convex = 0;
  l_max_spring = .1;
  n_averages = 256;

  compressed = old_compressed = mrisCountCompressed(mris, min_dist);
  if (compressed > 0) {
    printf("removing compressed regions in tessellation\n");
  }
  nsize = mris->nsize;
  MRISsetNeighborhoodSizeAndDist(mris, 3);
  MRIScomputeSecondFundamentalForm(mris);
  sprintf(fname, "uncompress%4.4d", iters);
  MRISwrite(mris, fname);
  MRISclearGradient(mris);
  parms.flags |= IPFLAG_NO_SELF_INT_TEST;
  while ((compressed > 0) && (iters++ < 500)) {
    MRISclearGradient(mris);
    if ((parms.flags & IPFLAG_NO_SELF_INT_TEST) == 0) {
      MHTfree(&mht);
      mht = MHTcreateFaceTable(mris);
    }
    else {
      mht = NULL;
    }
    if (Gdiag_no >= 0 && DIAG_VERBOSE_ON) {
      MRISprintVertexStats(mris, Gdiag_no, Gstdout, CURRENT_VERTICES);
    }

    //    mrisComputeNonlinearTangentialSpringTerm(mris, 1, min_dist) ;
    //    mrisComputeTangentialSpringTerm(mris, 1) ;
    //    mrisComputeNormalizedSpringTerm(mris, 1) ;
    if (0)  // disable
    {
      mrisComputePlaneTerm(mris, 0, 0);
    }
    mrisComputeSpringTerm(mris, l_spring);
    mrisComputeBorderTerm(mris, 1 - l_spring);
    mrisComputeConvexityTerm(mris, l_convex);
    MRISaverageGradients(mris, n_averages);
    mrisComputeMaxSpringTerm(mris, l_max_spring);
    mrisComputeAngleAreaTerms(mris, &parms);
    
    delta_t = mrisAsynchronousTimeStep(mris, 0, .2, mht, min_dist);

    MRIScomputeMetricProperties(mris);
    MRIScomputeSecondFundamentalForm(mris);
    old_compressed = compressed;
    compressed = mrisCountCompressed(mris, min_dist);
    sprintf(fname, "uncompress%4.4d", iters);
    printf("ITER %d: compressed = %d, delta = %d, writing %s\n", iters, compressed, compressed - old_compressed, fname);
    MRISwrite(mris, fname);
  }

  if (mht) {
    MHTfree(&mht);
  }
  MRISsetNeighborhoodSizeAndDist(mris, nsize);
  return (NO_ERROR);
}

