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
#define COMPILING_MRISURF_TOPOLOGY
 
#include "mrisurf_topology.h"

#include "mrisurf_base.h"


// MRIS code dealing with the existence and connectedness of the vertices, edges, and faces
//                   and with their partitioning into sets (ripped, marked, ...)
//                   but not with their placement in the xyz coordinate space


//=============================================================================
// Support for consistency checking
//
static bool shouldReportWkr(int line) {
  bool wasReported;
  if (copeWithLogicProblem2(&wasReported, 
            "FREESURFER_crash_on_bad_topology", "Bad vertex or face found", 
            __FILE__, line, "") == LogicProblemResponse_fix)
    cheapAssert(false);
  return wasReported;
}


//=============================================================================
// Vertexs and edges
//
bool mrisCheckVertexVertexTopologyWkr(const char* file, int line, MRIS const *mris, bool always)
{
  { static bool laterTime, forceAlways;
    if (!laterTime) { 
      laterTime = true; 
      forceAlways = !!getenv("FREESURFER_always_check_vv"); 
      if (forceAlways) fprintf(stdout,"%s:%d FREESURFER_always_check_vv set\n", __FILE__, __LINE__);
    }
    always |= forceAlways;
  }
  if (!always && !spendTimeCheckingForLogicProblem(file, line)) return true;

  int reported = 0;
  #define shouldReport (shouldReportWkr(__LINE__) && (reported=1))    // sets reported only if should report
  
  if (mris->nsize > mris->max_nsize) {
    if (shouldReport) {
      fprintf(stdout, "mris->nsize:%d > mris->max_nsize:%d\n", mris->nsize, mris->max_nsize);
    }
  }
  
  int vno1;
  for (vno1 = 0; vno1 < mris->nvertices; vno1++) {
    VERTEX_TOPOLOGY const * const v         = &mris->vertices_topology[vno1];
    VERTEX          const * const v_nontopo = &mris->vertices         [vno1];

#if 0
    // HACK TO HELP FIND A SPECIFIC BUG
    //
    if (vno1 == 0) {
      typedef struct Debugging { bool dist_orig_is_0; } Debugging;        
      if (!v_nontopo->debugging) {
        mris->vertices[vno1].debugging = calloc(1,sizeof(Debugging));
      }
      Debugging* debugging = (Debugging*)(v_nontopo->debugging);
      bool dist_orig_is_0 = !v_nontopo->dist_orig || !v_nontopo->dist_orig[0];
      if (debugging->dist_orig_is_0 != dist_orig_is_0) {
        fprintf(stdout,"%s:%d dist_orig_is_0:%d\n",file,line,dist_orig_is_0);
        debugging->dist_orig_is_0 = dist_orig_is_0;
      }
    }
#endif
    
    if ( (v_nontopo->dist      && !(mris->dist_alloced_flags&1))
      || (v_nontopo->dist_orig && !(mris->dist_alloced_flags&2))) {
      if (shouldReport) {
        fprintf(stdout, "dist or dist_orig non-null when !mris->dist_alloced\n");
      }
    }

    if (mris->vertices[vno1].ripflag) continue;
      
    int vSize = mrisVertexVSize(mris, vno1);

    if (mris->nsize > 0 
     && mris->nsize != v->nsizeCur 
     && vSize > 0) {
      if (shouldReport) {
        fprintf(stdout, "[vno1:%d].nsizeCur:%d != mris->nsize:%d vSize:%d vnum:%d vtotal:%d\n", 
          vno1, v->nsizeCur, mris->nsize, vSize, v->vnum, v->vtotal);
      }
    }

    if (v->nsizeMax > 0 &&
        v->nsizeCur > v->nsizeMax ) {
      if (shouldReport) {
        fprintf(stdout, "[vno1:%d].nsizeCur:%d exceeds nsizeMax:%d\n", vno1, v->nsizeCur, v->nsizeMax);
      }
    }

    int const vtotalExpected = !v->nsizeCur ? 0 : VERTEXvnum(v,v->nsizeCur);
    
    if (v->nsizeCur > 0 
     && (mris->vtotalsMightBeTooBig ? (v->vtotal < vtotalExpected) : (v->vtotal != vtotalExpected))) {
      //
      // MRISsampleDistances sets vtotal beyond vtotalExpected
      //
      if (shouldReport) {
        fprintf(stdout, "[vno1:%d].vtotal:%d differs from expected:%d for nsize:%d, ripflag:%d\n", vno1, v->vtotal, vtotalExpected, v->nsizeCur, 0);
      }
    }

    int n;
    for (n = 0; n < vSize; n++) {
      int vno2 = v->v[n];

      if (vno2 < 0 || mris->nvertices <= vno2) {
        if (shouldReport) {
          fprintf(stdout, "[vno1:%d].v[%d] is bad vno2:%d\n", vno1, n, vno2);
        }
      }

      if (v->vnum <= n) continue;
      
      // immediate neighborlyness is commutative
      if (!mrisVerticesAreNeighbors(mris, vno2, vno1)) {
        if (shouldReport) {
          fprintf(stdout, "[vno1:%d].v[%d] not found in [vno2:%d].v[*]\n", vno1, n, vno2);
        }
      }
      
      // neighbors should only appear once
      int i;
      for (i = 0; i < n; i++) {
        if (vno2 == v->v[i]) {
          if (shouldReport) {
            fprintf(stdout, "[vno1:%d].v[%d]:%d same as [vno1:%d].v[%d]\n", vno1, n, vno2, i, v->v[i]);
          }
        }
      }
    }
  }
  
  static int laterTime, checkDist;
  if (!laterTime) { 
    laterTime++; 
    checkDist = !!getenv("FREESURFER_checkDist"); 
    if (checkDist) fprintf(stdout, "%s:%d checking dist[] contents\n", __FILE__, __LINE__);
  }
  if (checkDist) {
    extern bool mrisCheckDist    (MRIS const * mris);               // gross hack for now
    extern bool mrisCheckDistOrig(MRIS const * mris);               // gross hack for now
    if (!mrisCheckDist    (mris)) reported = 1;
    if (!mrisCheckDistOrig(mris)) reported = 1;
  }

#undef shouldReport

  return reported == 0;
}


static int  MRIScheckIsPolyhedronCalls;

void MRIScheckIsPolyhedron(MRIS *mris, const char* file, int line) {
  MRIScheckIsPolyhedronCalls++;
  
  return ;  // HACK
  
  bool const tearsOk = true;
  
  static const char* prevFile = "<none>";
  static        int  prevLine = 0;
  
  int const nvertices = mris->nvertices;
  int const nfaces    = mris->nfaces;

  typedef struct Edges {
    size_t hiVnosCount;
    size_t hiVnosBegin;
  } Edges;
  Edges* edges  = (Edges*)malloc(nvertices*sizeof(Edges));
   
  size_t hiVnosCapacity = 0;
  int*   hiVnos         = NULL;
  size_t hiVnosSize     = 0;

  // Record all the edges - defined by 1-hop neighbours
  // and check that they are commutative and unique
  //
  size_t maxNeighbours;
  for (maxNeighbours = 16; maxNeighbours < 1024; maxNeighbours *= 2) {
    hiVnosCapacity = nvertices*maxNeighbours;
    hiVnos = (int*)malloc(hiVnosCapacity * sizeof(int));
    
    int vno;
    for (vno = 0; vno < nvertices; vno++) {
      size_t const hiVnosBegin = hiVnosSize;
      
      if (mris->vertices[vno].ripflag) continue;
      
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      
      int n;
      for (n = 0; n < vt->vnum; n++) {
        int const vno2 = vt->v[n];
        if (mris->vertices[vno2].ripflag) continue;
        if (vno < vno2) {
          if (hiVnosSize == hiVnosCapacity) goto TryLargerCapacity;
          hiVnos[hiVnosSize++] = vno2;
        } else {
          Edges const * const edgesForVno2 = &edges[vno2];
          int count = 0;
          for (unsigned int i = 0; i < edgesForVno2->hiVnosCount; i++) {
            if (hiVnos[edgesForVno2->hiVnosBegin + i] == vno) count++;
          }
          cheapAssert(count == 1);  // each should be entered exactly once
        }
      }

      Edges* const edgesForVno = &edges[vno];
      edgesForVno->hiVnosCount = hiVnosSize - hiVnosBegin;
      edgesForVno->hiVnosBegin = hiVnosBegin;
    }
    break;
    
  TryLargerCapacity:
    freeAndNULL(hiVnos);
  }
  
  // Count all the faces contributing to edges
  //
  char* nFacesPerEdge = (char*)calloc(hiVnosSize, sizeof(char));
  bool badCount = false;
  
  do {
  
    int reports1 = 0, reportsNot1 = 0;
    int fno;
    for (fno = 0; fno < nfaces; fno++) {
      FACE const * const face = &mris->faces[fno];
      int n1;
      for (n1 = 0; n1 < VERTICES_PER_FACE; n1++) {
        int const n2 = (n1 > 0) ? n1-1 : VERTICES_PER_FACE-1;
        int vno1 = face->v[n1];
        int vno2 = face->v[n2];
        if (vno1 > vno2) { vno1 = vno2; vno2 = face->v[n1]; } else cheapAssert(vno1 != vno2);
        Edges const * const edgesForVno = &edges[vno1];
        for (unsigned int i = 0; i < edgesForVno->hiVnosCount; i++) {
          if (hiVnos   [edgesForVno->hiVnosBegin + i] != vno2) continue;
          // found the entry for an edge
          if (!badCount) {
            nFacesPerEdge[edgesForVno->hiVnosBegin + i]++;  // yeah, might wrap - unlikely
          } else {
            int const count = nFacesPerEdge[edgesForVno->hiVnosBegin + i];
            if (count != 2) {
              if ((count == 1) ? (reports1++    < 10)
                               : (reportsNot1++ < 10)) {
                fprintf(stdout, "Fno %d contributes edge %d..%d towards a count of %d\n", 
                  fno, vno1, vno2, count);
              }
            }
          }
        }
      }
    }
    
    cheapAssert(!reportsNot1);
    
    if (badCount) break;

    for (unsigned int edgeNo = 0; edgeNo < hiVnosSize; edgeNo++) {
      int count = nFacesPerEdge[edgeNo];
      if (count == 2) continue;
      if (count <  2 && tearsOk) continue;
      badCount = true;
      fprintf(stdout, "MRIScheckIsPolyhedron failed during call %d at %s:%d, previously good at %s:%d because count:%d for edgeNo:%d\n", 
        MRIScheckIsPolyhedronCalls, file, line, prevFile, prevLine, count, edgeNo);
      break;
    }
    
  } while (badCount);
    
  freeAndNULL(nFacesPerEdge);
  freeAndNULL(hiVnos);
  freeAndNULL(edges);

  prevFile = file;
  prevLine = line;
}
int MRISvalidVertices(MRIS *mris)
{
  int vno, nvertices, nvalid;

  nvertices = mris->nvertices;
  for (vno = nvalid = 0; vno < nvertices; vno++)
    if (!mris->vertices[vno].ripflag) {
      nvalid++;
    }

  return (nvalid);
}


int edgeExists(MRI_SURFACE *mris, int vno1, int vno2)
{
  int n;
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno1];
  for (n = 0; n < v->vnum; n++)
    if (v->v[n] == vno2) {
      return (1);
    }
  return (0);
}


/*-----------------------------------------------------
  Go through all the (valid) faces in vno1 and see how many
  valid links to vno2 exists.
  ------------------------------------------------------*/
int mrisCountValidLinks(MRIS *mris, int vno1, int vno2)
{
  int nvalid, fno, vno;
  FACE *face;
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno1];
  for (nvalid = fno = 0; fno < vt->num; fno++) {
    face = &mris->faces[vt->f[fno]];
    if (face->ripflag) {
      continue;
    }
    for (vno = 0; vno < VERTICES_PER_FACE; vno++)
      if (face->v[vno] == vno1) {
        if ((vno == VERTICES_PER_FACE - 1) && (face->v[0] == vno2)) {
          nvalid++;
        }
        else if ((vno < VERTICES_PER_FACE - 1) && (face->v[vno + 1] == vno2)) {
          nvalid++;
        }
      }
  }
  return (nvalid);
}


short modVnum(MRIS const *mris, int vno, short add, bool clear) 
{
  VERTEX_TOPOLOGY* vt = &mris->vertices_topology[vno];
  short * p = const_cast<short*>(&vt->vnum);
  if (clear) *p = 0;
  *p += add;
  static short maxVnumSeen = 30;
  if (*p > maxVnumSeen) {
    maxVnumSeen = *p;
    fs::debug() << "modVnum: vertex " << vno << " has " << *p << " immediate neighbours";
  }
  return *p;
}


static void mrisAddEdgeWkr(MRIS *mris, int vno1, int vno2) 
{
  if (vno1 < 0 || vno2 < 0) {
    DiagBreak();
  }
  if (vno1 == Gdiag_no || vno2 == Gdiag_no) {
    DiagBreak();
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "adding edge %d <--> %d\n", vno1, vno2);
  }

  auto adjustV = [&](int vnoA, int vnoB) {  
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vnoA];
    v->v = (int*)realloc(v->v, (v->vnum+1)*sizeof(int));
    if (!v->v) ErrorExit(ERROR_NO_MEMORY, "mrisAddEdge(%d, %d): could not allocate %d len vlist", v->vnum);

    v->v[vnumAdd(mris,vnoA,1)] = vnoB;
    v->vtotal = v->vnum; v->nsizeCur = v->nsizeMax = 1;

    static int maxVnum = 20;
    if (maxVnum < v->vnum) {
      maxVnum = v->vnum;
      fs::debug() << "vertex " << vnoA << " has " << v->vnum << " immediate neighbours";
    }
  };
  
  adjustV(vno1,vno2);
  adjustV(vno2,vno1);
}

void mrisAddEdge(MRIS *mris, int vno1, int vno2)
{
  costlyAssert(!mrisVerticesAreNeighbors(mris, vno1, vno2));
  mrisAddEdgeWkr(mris, vno1, vno2);
}


void mrisRemoveEdge(MRIS *mris, int vno1, int vno2)
{
  // BUG doesn't adjust v2num etc.
  
  VERTEX_TOPOLOGY * const v1 = &mris->vertices_topology[vno1];

  int i;
  for (i = 0; i < v1->vnum; i++)
    if (v1->v[i] == vno2) {
      addVnum(mris,vno1,-1);
      if (i < v1->vnum) /* not the (previous) last index */
      {
        memmove(&v1->v[i], &v1->v[i + 1], (v1->vnum - i) * sizeof(v1->v[0]));
      }
      return ;
    }
  // TODO what about any faces that use this edge?
}


int mrisRemoveVertexLink(MRIS *mris, int vno1, int vno2)
{
  // BUG doesn't adjust v2num etc.
  
  int n;
  VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno1];
  for (n = 0; n < v->vnum; n++)
    if (v->v[n] == vno2) {
      break;
    }

  if (n < v->vnum) {
    memmove(v->v + n, v->v + n + 1, (v->vtotal - (n + 1)) * sizeof(int));
    addVnum(mris,vno1,-1);
    v->vtotal--;
  }
  return (NO_ERROR);
}


int MRISremoveTriangleLinks(MRIS *mris)
{
  int fno, which;
  FACE *f;

  if (!IS_QUADRANGULAR(mris)) {
    return (NO_ERROR);
  }
  if (mris->triangle_links_removed) {
    return (NO_ERROR);
  }

  mris->triangle_links_removed = 1;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "removing non-quadrangular links.\n");
  }

  for (fno = 0; fno < mris->nfaces; fno += 2) {
    f = &mris->faces[fno];
    if (f->ripflag) {
      continue;
    }
    which = WHICH_FACE_SPLIT(f->v[0], f->v[1]);
    if (EVEN(which)) {
      mrisRemoveVertexLink(mris, f->v[1], f->v[2]);
      mrisRemoveVertexLink(mris, f->v[2], f->v[1]);
    }
    else {
      mrisRemoveVertexLink(mris, f->v[0], f->v[2]);
      mrisRemoveVertexLink(mris, f->v[2], f->v[0]);
    }
  }
  return (NO_ERROR);
}


int MRIScountEdges(MRIS* mris)
{
  int nedges = 0;

  int vno;
  for (vno=0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    int n;
    for (n = 0; n < vt->vnum; n++) {
      if (vt->v[n] < vno) continue;
      nedges++;
    }
  }

  return nedges;
}


int mrisCountTotalNeighbors(MRIS* mris)
{
  int total = 0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;
    total += vt->vtotal + 1; /* include this vertex in count */
  }
  return (total);
}



//=============================================================================
// EulerNumber and similar calculations
//
int MRIScomputeEulerNumber(MRIS *mris, int *pnvertices, int *pnfaces, int *pnedges)
{
  int eno, nfaces, nedges, nvertices, vno, fno, vnb, i;

  for (nfaces = fno = 0; fno < mris->nfaces; fno++)
    if (!mris->faces[fno].ripflag) {
      nfaces++;
    }

  for (nvertices = vno = 0; vno < mris->nvertices; vno++)
    if (!mris->vertices[vno].ripflag) {
      nvertices++;
    }

  for (nedges = vno = 0; vno < mris->nvertices; vno++)
    if (!mris->vertices[vno].ripflag) {
      VERTEX_TOPOLOGY const * const v1 = &mris->vertices_topology[vno];
      for (i = 0; i < v1->vnum; i++) {
        vnb = v1->v[i];
        /* already counted */
        if ((vnb > vno) && !mris->vertices[vnb].ripflag) {
          nedges++;
        }
      }
    }

  *pnfaces = nfaces;
  *pnvertices = nvertices;
  *pnedges = nedges;
  eno = nvertices - nedges + nfaces;
  return (eno);
}


int MRIStopologicalDefectIndex(MRIS *mris)
{
  int eno, nfaces, nedges, nvertices, vno, fno, vnb, i, dno;

  for (nfaces = fno = 0; fno < mris->nfaces; fno++)
    if (!mris->faces[fno].ripflag) {
      nfaces++;
    }

  for (nvertices = vno = 0; vno < mris->nvertices; vno++)
    if (!mris->vertices[vno].ripflag) {
      nvertices++;
    }

  for (nedges = vno = 0; vno < mris->nvertices; vno++)
    if (!mris->vertices[vno].ripflag) {
      VERTEX_TOPOLOGY const * const v1 = &mris->vertices_topology[vno];
      for (i = 0; i < v1->vnum; i++) {
        vnb = v1->v[i];
        /* already counted */
        if ((vnb > vno) && !mris->vertices[vnb].ripflag) {
          nedges++;
        }
      }
    }
#if 0
  nedges += nfaces ;        /* one additional edge added for each triangle */
  nfaces *= 2 ;             /* two triangular faces per face */
#endif
  eno = nvertices - nedges + nfaces;

  dno = abs(2 - eno) + abs(2 * nedges - 3 * nfaces);
  return (dno);
}


//=============================================================================
// MRI_EDGE support

/*
  Allocates and assigns edge structures. 
  I think the corner values are not correct.
*/
int MRISedges(MRIS *surf)
{
  int edgeno = 0;
  int vtxno0, vtxno1;

  surf->nedges = MRIScountEdges(surf);
  surf->edges  = (MRI_EDGE *)calloc(surf->nedges,sizeof(MRI_EDGE));
  //printf("MRISedges(): nv=%d, nf=%d, ne=%d\n",surf->nvertices,surf->nfaces,surf->nedges);

  // This is not thread safe and cannot be made thread safe
  for(vtxno0=0; vtxno0 < surf->nvertices; vtxno0++){
    VERTEX_TOPOLOGY const * const v0t = &surf->vertices_topology[vtxno0];
    
    int nthnbr;
    for(nthnbr = 0; nthnbr < v0t->vnum; nthnbr++){
      int const vtxno1 = v0t->v[nthnbr];
      if(vtxno1 < vtxno0) continue;
      
      
      VERTEX_TOPOLOGY const * const v1t = &surf->vertices_topology[vtxno1];
      MRI_EDGE * const edge = &(surf->edges[edgeno]);
      edge->vtxno[0] = vtxno0;
      edge->vtxno[1] = vtxno1;

      // Find the two faces that adjoin this edge
      {
        int k = 0; // k is the nthface out of 2
        int n;
        for(n=0; n < v0t->num; n++){ // go thru faces of vtx1
          int m;
	  for(m=0; m < v1t->num; m++){ // go thru faces of vtx2
	    if(v0t->f[n] == v1t->f[m]){ // same face
	      if(k>1){
	        printf("ERROR: MRISedges(): too many faces: %d %d n=%d, m=%d, k=%d\n",vtxno0,vtxno1,n,m,k);
	        return(1);
	      }
	      FACE const * const f = &(surf->faces[v0t->f[n]]);
	      // get the ordinal corner numbers for the two vertices
              int c;
	      for(c=0; c<2; c++){
	        if(f->v[c] == vtxno0) edge->corner[k][0] = c;
	        if(f->v[c] == vtxno1) edge->corner[k][1] = c;
	      }
	      edge->faceno[k] = v0t->f[n];
	      k++;
	    }
	  }
        }
        if(k != 2){
	  printf("ERROR: MRISedges(): not enough faces %d %d k=%d\n",vtxno0,vtxno1,k);
	  return(1);
        }
      }
      
      // Now find the two vertices of the faces that are not common between the faces;
      // these are the 3rd corner for the two triangles
      {
        int k;
        for(k=0; k<2; k++){ // go thru the two faces
	  FACE const * const f = &(surf->faces[edge->faceno[k]]); // kth face
          int c;
	  for(c=0; c<3; c++){ // go thru the vertices of the kth face
	    if(f->v[c] != vtxno0 && f->v[c] != vtxno1){
	      // this is the 3rd corner if it is not vtx1 or vtx2
	      edge->vtxno[k+2] = f->v[c];
	      break;
	    }
	  }
        }
      }
      
      // Fill the corner matrix
      {
        int n;
        for(n=0; n<4; n++){ // go thru the four edge vertices
	  int const m = edge->vtxno[n];
          int k;
	  for(k=0; k<2; k++){ // go thru the two faces
	    FACE const * const f = &(surf->faces[edge->faceno[k]]); // kth face
            int c;
	    for(c=0; c<3; c++){ // go thru the vertices of the kth face
	      // Compare the edge vertex no against the face vertex no
	      if(f->v[c] == m){
	        edge->corner[n][k] = c;
	        break;
	      }
	    }
	  }
        }
      }
      
      // Note: corner[2][1] and corner[3][0] are not used
      edge->corner[2][1] = 9;
      edge->corner[3][0] = 9;

      edgeno++;
    }
  }

  // Allocate the edge number array in each vertex
  VERTEX_TOPOLOGY *vt;
  for(vtxno0 = 0; vtxno0 < surf->nvertices; vtxno0++){
    vt = &(surf->vertices_topology[vtxno0]);
    if(vt->e) free(vt->e);
    vt->e = (int*)calloc(sizeof(int),vt->vtotal);
  }

  // Assign edge numbers to vertices
  // could add directionality here too
  int k;
  for(edgeno = 0; edgeno < surf->nedges; edgeno++){
    vtxno0 = surf->edges[edgeno].vtxno[0];
    vtxno1 = surf->edges[edgeno].vtxno[1];
    vt = &(surf->vertices_topology[vtxno0]);
    for(k=0; k < vt->vtotal; k++)
      if(vt->v[k] == vtxno1) vt->e[k] = edgeno;
    vt = &(surf->vertices_topology[vtxno1]);
    for(k=0; k < vt->vtotal; k++)
      if(vt->v[k] == vtxno0) vt->e[k] = edgeno;
  }

  return(0);
}


/*!
  \fn int MRIScorners(MRIS *surf)
  \brief Create triangle corner topology. Will return immediately with
  0 if already done. If edge topology has not been built, then that
  will be done as well. A corner is an angle of a triangle.
 */
int MRIScorners(MRIS *surf)
{
  int faceno, cornerno, k,m,n;
  FACE *face;
  MRI_CORNER *c;

  if(surf->corners) return(0);

  //printf("Building triangle corner toplology\n");
  if(!surf->edges) MRISedges(surf);

  surf->ncorners = 3*surf->nfaces;
  surf->corners = (MRI_CORNER*) calloc(sizeof(MRI_CORNER),surf->ncorners);
  //printf("MRIScorners(): nv=%d, nf=%d, ne=%d, nc=%d\n",
  //	 surf->nvertices,surf->nfaces,surf->nedges,surf->ncorners);

  // First assign vertices to each corner
  cornerno = 0;
  for(faceno=0; faceno < surf->nfaces; faceno++){
    face = &(surf->faces[faceno]);
    for(k=0; k < 3; k++){
      c = &(surf->corners[cornerno]);
      c->cornerno = cornerno;
      c->faceno = faceno;
      m = k + 1;
      if(m>2) m -= 3;
      n = k + 2;
      if(n>2) n -= 3;
      c->vtxno[0] = face->v[k];
      c->vtxno[1] = face->v[m];
      c->vtxno[2] = face->v[n];
      cornerno++;
    }
  }

  // Now assign edges and edge direction to each corner
  for(cornerno = 0; cornerno < surf->ncorners; cornerno++){
    c = &(surf->corners[cornerno]);
    VERTEX_TOPOLOGY *v = &(surf->vertices_topology[c->vtxno[0]]);
    // Go through the two cornder vertex neighbors (ie, edges) for this corner
    int nthedge;
    for(nthedge=0; nthedge < 2; nthedge++){
      // Go through all the vertex neighbors of the corner vertex
      int hit = 0;
      for(k=0; k < v->vtotal; k++){
	// If this vertex neighbor is the same vertex as this corner neighbor
	if(v->v[k] == c->vtxno[nthedge+1]){
	  // Set the edgeno for this edge to the edgeno connecting central vertex with neighbor
	  c->edgeno[nthedge] = v->e[k];
	  // Determine the direction of the edge. If the first vertex
	  // of the edge is the same as the center vertex of the
	  // corner, then they are pointing the in same direction
	  MRI_EDGE *edge = &(surf->edges[c->edgeno[nthedge]]);
	  if(edge->vtxno[0] == c->vtxno[0]) c->edgedir[nthedge] = +1;
	  else                              c->edgedir[nthedge] = -1;
	  hit = 1;
	  break;
	}
      } // loop over neighbors
      if(!hit){
	printf("MRIScorners(): ERROR: could not find an edge for corner %d, edge %d, %d->%d\n",
	       cornerno,nthedge,c->vtxno[0],c->vtxno[nthedge+1]);
	fflush(stdout);
      }
    } // loop over edges
  }

  return(0);
}



//=============================================================================
// Neighbourhoods
//
// These are the vertexs that can be reached by following vno = mris->vertex_topology[vno]->v[<all>] nsize hops.
// They are stored, sorted, in the v->v vector with the vnum, v2num, and v3num storing where the hop count changes.
//
// Obviously adding or removing edges invalids v2num and v3num.
//      mris->nsizeMaxClock changes to help detect this bug, but it wraps so is not a guarantee.
//      Here it is checked to assert vtotal is valid.
//
void MRIS_setNsizeCur(MRIS *mris, int vno, int nsize) {
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    

  cheapAssert(nsize <= vt->nsizeMax);
  
  switch (nsize) {
  case 0: vt->vtotal = 0;         break;    // seen to happen
  case 1: vt->vtotal = vt->vnum;  break;
  case 2: vt->vtotal = vt->v2num; break;
  case 3: vt->vtotal = vt->v3num; break;
  default: cheapAssert(false);
  }

  uchar const * pc = &vt->nsizeCur;
  uchar       * p  = (uchar*)pc;
  *p = nsize;
}

int mrisStoreVtotalInV3num(MRIS *mris)
{
  // WEIRD find out who is doing this and why
  
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno];
    v->v3num = v->vtotal;
  }
  
  return (NO_ERROR);
}


static void resizeVertexV(MRIS* mris, int vno, int newSize, int oldSize) {

  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno]; 
  VERTEX          * const v  = &mris->vertices         [vno];

  // allocating zero is a free: keep the pointers around to optimize growing them again
  //           non-zero:        change to the new size
  
  if (newSize > 0) {
    int const intSize = newSize*sizeof(int);
    vt->v = (int*)realloc(vt->v, intSize);
  }
      
  // Zero the added storage, if any
  //
  if (oldSize < newSize) {
    int const intSizeChange = (newSize - oldSize)*sizeof(int);
    bzero(vt->v + oldSize, intSizeChange);
  }
  
  // Follow the tradition of keep ::dist at least as big as ::v
  // but don't make it if not already made
  //
  if (v->dist) MRISgrowDist(mris, vno, newSize);
}


void mrisVertexReplacingNeighbors(MRIS * const mris, int const vno, int const vnum)
{
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];

  mris->max_nsize = vt->nsizeMax = vt->nsizeCur = 1; 
  if (mris->nsize != 1) mris->nsize = -1;               // vertices are now inconsistent
  MRISfreeDistsButNotOrig(mris);                        // distances are wrong
  
  int const old_vt_vnum = vt->vnum;

  resizeVertexV(mris,vno, vnum, old_vt_vnum);

  modVnum(mris,vno,vnum,true);
  vt->v2num = 0;
  vt->v3num = 0;
  
  vt->vtotal = vt->vnum;
}


void mrisForgetNeighborhoods(MRIS * const mris) {
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    mrisVertexReplacingNeighbors(mris, vno, mris->vertices_topology[vno].vnum);
  }
  mris->nsize = 1;                                     // vertices are now consistent again
}


int mrisVertexVSize(MRIS const * mris, int vno) {
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];
  int c = 0;
  switch (v->nsizeMax) {
  case 1: c = v->vnum;  break;
  case 2: c = v->v2num; break;
  case 3: c = v->v3num; break;
  default: break;
  }
  if (c < v->vtotal) c = v->vtotal;
  return c;
}


int MRIScountTotalNeighbors(MRIS *mris, int nsize)
{
  int vno, total_nbrs;

  for (total_nbrs = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    cheapAssert(nsize <= vt->nsizeMax);
    switch (nsize) {
      case 0:
        total_nbrs++;               // weird, because selves aren't counted in the other cases
        break;
      case 1:
        total_nbrs += vt->vnum;
        break;
      case 2:
        total_nbrs += vt->v2num;
        break;
      case 3:
        total_nbrs += vt->v3num;
        break;
      default:
        ErrorExit(ERROR_UNSUPPORTED, "MRIScountNeighbors(%d): max nbhd size = 3", nsize);
    }
  }
  return (total_nbrs);
}


int MRISresetNeighborhoodSize(MRI_SURFACE *mris, int nsize)
{
  int new_mris_nsize = nsize;
  int ntotal         = 0;
  int vtotal         = 0;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
    VERTEX          * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      continue;
    }
    switch (nsize) {
      default: /* reset back to original */
        cheapAssert(vt->nsizeMax > 0);
        MRIS_setNsizeCur(mris, vno, vt->nsizeMax);
        if (new_mris_nsize < 0) new_mris_nsize = vt->nsizeMax;
        else                    cheapAssert(new_mris_nsize == vt->nsizeMax);
        break;
      case 1:
      case 2:
      case 3:
        MRIS_setNsizeCur(mris, vno, nsize);
        break;
    }

    vtotal += vt->vtotal;
    ntotal++;
  }
  mris->nsize = new_mris_nsize;
  mris->avg_nbrs = (float)vtotal / (float)ntotal;

  mrisCheckVertexFaceTopology(mris);

  return (NO_ERROR);
}


void MRISgetNeighborsBeginEnd(
  MRIS const * mris, 
  int     vno, 
  size_t  inner_nbhd_size, 
  size_t  outer_nbhd_size, 
  size_t* neighborsIndexBegin,    // set so VERTEX v[*neighborsIndexBegin] is the first in this list with inner <= links <= outer 
  size_t* neighborsIndexEnd) {    // set so VERTEX v[*neighborsIndexEnd]   is the first in this list with outer < links
    
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    
  size_t b = 0, e = 0;
    
  if (inner_nbhd_size <= outer_nbhd_size) {

    cheapAssert(mris->nsizeMaxClock == vt->nsizeMaxClock || outer_nbhd_size <= 1);

    switch (inner_nbhd_size) {
    case 1: b = 0;         break;
    case 2: b = vt->vnum;  break;
    case 3: b = vt->v2num; break;
    default: cheapAssert(false);
    }
    switch (outer_nbhd_size) {
    case 1: e = vt->vnum;  break;
    case 2: e = vt->v2num; break;
    case 3: e = vt->v3num; break;
    default: cheapAssert(false);
    }
  }
    
  *neighborsIndexBegin = b;
  *neighborsIndexEnd   = e;
}


void MRIS_VertexNeighbourInfo_check(
    MRIS_VertexNeighbourInfo* lhs,  
    MRIS_VertexNeighbourInfo* rhs) {
    
  cheapAssert(lhs->hops == rhs->hops);

  for (unsigned int i = 0; i <= lhs->hops; i++)
    cheapAssert(lhs->vnum[i] == rhs->vnum[i]); 
    
  // Even test the algorithms put them in the same order!
  //
  for (int i = 0; i < lhs->vnum[lhs->hops]; i++) 
    cheapAssert(lhs->v[i] == rhs->v[i]);
}


void MRIS_VertexNeighbourInfo_load_from_VERTEX (MRIS_VertexNeighbourInfo* info, MRIS* mris, int vno) {
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  info->hops = vt->nsizeMax;
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#endif
  switch (info->hops) {
  default: cheapAssert(false);
  case 3: info->vnum[3] = vt->v3num;
  case 2: info->vnum[2] = vt->v2num;
  case 1: info->vnum[1] = vt->vnum;
  case 0: info->vnum[0] = 1;
  }
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
  int i;
  for (i = 0; i < info->vnum[info->hops]; i++) info->v[i] = vt->v[i];
}

void MRIS_VertexNeighbourInfo_load_from_vlist (MRIS_VertexNeighbourInfo* info, MRIS* mris, int vno, size_t listSize, int* vlist, int* hops) {
  info->vnum[0] = 1;
  for (unsigned int i = 0; i < listSize; i++) {
    info->v[i] = vlist[i];
    info->hops = hops[i];
    // cheapAssert(info->hops >= 0); Always true for unsigned
    cheapAssert(info->hops <= 3);
    info->vnum[info->hops] = i+1;
  }
}

void MRIS_VertexNeighbourInfo_load_by_algorithm(MRIS_VertexNeighbourInfo* info, MRIS* mris, int vno) {
  // This algorithm is deliberately simple since it needs to be definitive and is not performance critical
}


// Fills the vlist parameter with the indices of the vertices up to and include nlinks hops along edges.
//
    // assumes mris->vertices[*].mark are all zero
    // leaves them zero

static int MRISfindNeighborsAtVertex_new(MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops, bool noCache);
static int MRISfindNeighborsAtVertex_old(MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops);

int MRISfindNeighborsAtVertex(MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops)
{
  static int  laterTime, interestingLaterTime;
  static bool use_new, use_old;
  
  if (!laterTime++) {
    use_new = getenv("MRISfindNeighborsAtVertex_new");
    use_old = getenv("MRISfindNeighborsAtVertex_old");
    if (!(use_new || use_old)) use_new = true;
  }
  bool const use_both = use_new && use_old;

  int const breakLine = __LINE__ + 2;
  if (laterTime == interestingLaterTime)
    fprintf(stdout, "%s:%d laterTime:%d\n", __FILE__, __LINE__, laterTime);  

  int result_old = use_old
    ? MRISfindNeighborsAtVertex_old(mris, vno, nlinks, listCapacity, vlist, hops)
    : 0;

  int * vlistTmp = NULL, *hopsTmp = NULL;
  if (use_old && use_new) {
    vlistTmp = (int*)malloc(MAX_NEIGHBORS*sizeof(int));
    hopsTmp  = (int*)malloc(MAX_NEIGHBORS*sizeof(int));
  }
   
  int result_new = use_new
    //
    // when testing, concatenates the two lists
    //
    ? MRISfindNeighborsAtVertex_new(mris, vno, nlinks, 
        use_both ? MAX_NEIGHBORS : listCapacity, 
        use_both ? vlistTmp      : vlist, 
        use_both ? hopsTmp       : hops,
        false)
    : 0;

  if (true && nlinks < 4) {
    int i;
    for (i = 0; i < result_new; i++) {
      if (mris->vertices_topology[vno].v[i] < 0 
      ||  mris->vertices_topology[vno].v[i] >= mris->nvertices) {
        fprintf(stdout, "%s:%d Bad mris->vertices_topology[%d].v[%d]\n", __FILE__, __LINE__, vno, i);
        cheapAssert(false);
      }
    }
  }
  if (use_both) {
    bool good = true;
    static bool shownTesting;
    if (!shownTesting) {
      shownTesting = true;
      fprintf(stdout, "%s:%dTesting MRISfindNeighborsAtVertex\n", __FILE__, __LINE__);
    }
    if (result_old != result_new) {
      fprintf(stdout, " result_old:%d != result_new:%d\n",result_old,result_new);
      good = false;
    }
    int i;
    for (i = 0; i < result_old && i < result_new;  i++) {
      if (vlist[i] != vlistTmp[i]
      ||  hops [i] != hopsTmp[i]) {
        fprintf(stdout, " vlist[%d] old:%d != new:%d  or  hops[%d] old:%d != new:%d\n",
          i, vlist[i],vlistTmp[i],
          i, hops[i], hopsTmp[i]);
        good = false;
      }
    }
    
    if (!good) {
      fprintf(stdout, "%s:%d laterTime:%d interestingLaterTime:%d\n", __FILE__, breakLine, laterTime, interestingLaterTime);
      cheapAssert(false);
    }
  }

  if (use_both) { freeAndNULL(vlistTmp); freeAndNULL(hopsTmp); }
  
  return use_old ? result_old : result_new;
}

void MRIS_check_vertexNeighbours(MRIS* mris) {

  static bool laterTime, doTesting;
  if (!laterTime) {
    laterTime = true;
    doTesting = !!getenv("MRIS_check_vertexNeighbours");
    if (doTesting) fprintf(stdout, "%s:%d MRIS_check_vertexNeighbours\n", __FILE__, __LINE__);
  }
  if (!doTesting) return;  

  MRIS_VertexNeighbourInfo info0 = MRIS_VertexNeighbourInfo();
  MRIS_VertexNeighbourInfo info1 = MRIS_VertexNeighbourInfo();
  int vlist[MAX_NEIGHBORS], hops[MAX_NEIGHBORS];
  
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {

    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
    VERTEX          * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;
    
    MRIS_VertexNeighbourInfo_load_from_VERTEX(&info0, mris, vno);
    
    int size = MRISfindNeighborsAtVertex_new(mris, vno, vt->nsizeMax, MAX_NEIGHBORS, vlist, hops, true);
    
    MRIS_VertexNeighbourInfo_load_from_vlist(&info1, mris, vno, size, vlist, hops);
    
    MRIS_VertexNeighbourInfo_check(&info0, &info1);
  }
}

static int MRISfindNeighborsAtVertex_newWkr(
    MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops, bool noCache, bool debug = false)
{
/*
  Fills in v, vnum, v2num, v3num, etc. in the vertex.
  However nlinks may be much higher than these.

  There are two issues here
    1) adding or removing an edge should invalidate the cached vertex.v#num info
    2) changing a ripflag should invalidate the cached vertex.v#num info
    
  The old code did not detect either case and used the invalid cache entry for the answer.
  
  The new code usually asserts in when the edges have changed, and always asserts if a ripped vertex is encountered.
    In the future it may need to recompute if a ripped vertex is encountered.
*/

  // Get an empty set
  // Note: This code always clears the temp it uses after use, is quick since it knows which ones it set
  //
  unsigned char const Status_notInSet    = 0;   // assumed to be 0 below
  unsigned char const Status_inSet       = 1;

  typedef struct Temp {
    size_t          capacity;
    unsigned char * status;
  } Temp;
  
  static Temp tempForEachThread[_MAX_FS_THREADS];
  
  Temp* const temp = &tempForEachThread[omp_get_thread_num()];
  if (temp->capacity < (unsigned)mris->nvertices) {
    temp->capacity = mris->nvertices;
    temp->status   = (unsigned char*)realloc(temp->status, temp->capacity*sizeof(unsigned char));
    bzero(temp->status, temp->capacity*sizeof(unsigned char));
  }  
  
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
  VERTEX          * const v  = &mris->vertices         [vno];

  cheapAssert(!v->ripflag);

  // The following is easier with these in an array
  //
  cheapAssert(nlinks <= MRIS_MAX_NEIGHBORHOOD_LINKS);
  int vnums[MRIS_MAX_NEIGHBORHOOD_LINKS + 1]; 
  int nsize = 0;
  {
                           vnums[nsize++] = 0; 
    if (vt->nsizeMax >= 1) vnums[nsize++] = vt->vnum;
    if (!noCache) {
      if (vt->nsizeMax >= 2) vnums[nsize++] = vt->v2num; 
      if (vt->nsizeMax >= 3) vnums[nsize++] = vt->v3num;
    }
  }
  // vnums[nsize-1] is valid, vnums[nsize] is not
  
  // The center is assumed to be in the set, so it is not added again
  //
  temp->status[vno] = Status_inSet; // already in the set

  // Add the known rings
  //
  unsigned int neighborCount = 0;

  cheapAssert(vt->nsizeMax > 0);

  // add the known rings, 
  // during this loop ringLinks is the ring being added
  //
  int ringLinks;
  for (ringLinks = 1; ringLinks < nsize && ringLinks <= nlinks; ringLinks++) {
  
    int ringBegin = vnums[ringLinks-1];
    int ringEnd   = vnums[ringLinks  ];

    // Add all the elements known to be in the ring
    //
    int n;
    for (n = ringBegin; n < ringEnd; n++) {
      int vnoCandidate = vt->v[n];
      VERTEX* vCandidate = &mris->vertices[vnoCandidate];

      // TODO cope with added edges
      
      if (vCandidate->ripflag) { 
        nsize = ringLinks;      // cause this ring to get rewritten and no further rings to be used
        continue;               // other vCandidates in this ring are still okay
      }
      
      temp->status[vnoCandidate] = Status_inSet;

      if (neighborCount < listCapacity) {
        vlist[neighborCount] = vnoCandidate;
        hops [neighborCount] = ringLinks;
      }
      if (debug) fprintf(stdout, "  vnoCandidate:%d known to be in ring:%d\n", vnoCandidate, ringLinks-1);
      neighborCount++;
    }
  }
  
  // ringLinks is the first ring to be added, which will be 2 or more because 1
  // is the immediate neighbors.
  //
  cheapAssert(ringLinks >= 2);
  
  // compute the new rings
  //
  for (; ringLinks <= nlinks; ringLinks++) {
    //
    int knownRingBegin = vnums[ringLinks-2];
    int knownRingEnd   = vnums[ringLinks-1];
    
    // Scan all the vertexs in the current border ring
    // Add their immediate neighbors that are further away
    //
    int i;        
    for (i = knownRingBegin; i < knownRingEnd; i++) {
      int                     const vnoRing = vlist[i];
      VERTEX_TOPOLOGY const * const vtRing  = &mris->vertices_topology[vnoRing];
          
      int j;
      for (j = 0; j < vtRing->vnum; j++) {
        int      const vnoCandidate = vtRing->v[j];
        VERTEX * const vCandidate   = &mris->vertices[vnoCandidate];

        if (vCandidate->ripflag) continue;
        if (temp->status[vnoCandidate] != Status_notInSet) continue;

        temp->status[vnoCandidate] = Status_inSet;

        if (neighborCount < listCapacity) {
          vlist[neighborCount] = vnoCandidate;
          hops [neighborCount] = ringLinks;
        }
        if (debug) fprintf(stdout, "  vnoCandidate:%d is neighbour of vno:%d in ring:%d\n", vnoCandidate, vnoRing, i);
        neighborCount++;
      }
    }
    vnums[ringLinks] = MIN(listCapacity,neighborCount);
  }

  // Make nsize the highest current valid vt->nsizeCur
  //
  // nsize is the number of valid rings, which is 1 more than the last valid index
  // it might be 1 when there is a ripped immediate neighbour
  //
  cheapAssert(nsize > 0);
  nsize -= 1;
    
  // Update the cache
  //
  if (!noCache) {

    int const newPossibleNsizeMax = MIN(3, ringLinks-1);
    if (nsize < newPossibleNsizeMax) {

      int oldSize = vnums[nsize];
      int newSize = vnums[newPossibleNsizeMax];    
      resizeVertexV(mris,vno, newSize, oldSize);

      int i;
      for (i = oldSize; i < newSize; i++) vt->v[i] = vlist[i];

      int cachedRing;
      for (cachedRing = nsize+1; cachedRing <= newPossibleNsizeMax; cachedRing++) {
        switch (cachedRing) {
        case 1: modVnum(mris,vno,vnums[cachedRing],true);   break;   // happens when encounters ripped vertexs
        case 2: vt->v2num = vnums[cachedRing];              break;
        case 3: vt->v3num = vnums[cachedRing];              break;
        default: cheapAssert(false);
        }
        vt->nsizeMax = newPossibleNsizeMax; vt->nsizeMaxClock = mris->nsizeMaxClock; 
      }

      // cheapAssert(vt->nsizeCur >= 0); Always true for unsigned
      cheapAssert(vt->nsizeCur <= vt->nsizeMax);
      vt->vtotal = vnums[vt->nsizeCur];
    }
  }
    
  // Clear the temp for reuse later
  //
  temp->status[vno] = Status_notInSet;
  auto const usedNeighborCount = MIN(listCapacity,neighborCount);
  for (unsigned int i = 0; i < usedNeighborCount; i++) {
    temp->status[vlist[i]] = Status_notInSet;
  }
  
  // Done
  //
  return neighborCount;
}

static int MRISfindNeighborsAtVertex_new(MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops, bool noCache)
{
  auto neighborCount = MRISfindNeighborsAtVertex_newWkr(mris, vno, nlinks, listCapacity, vlist, hops, noCache);
  if (neighborCount > listCapacity) {
    fprintf(stdout, "MRISfindNeighborsAtVertex: vno:%d has %d neighbours, but listCapacity is %ld\n", vno, neighborCount, listCapacity);
    MRISfindNeighborsAtVertex_newWkr(mris, vno, nlinks, listCapacity, vlist, hops, noCache, true);
    mrisCheckVertexVertexTopologyWkr(__FILE__,__LINE__,mris, true);
    cheapAssert(false);
  } 
  return neighborCount;
}


static int MRISfindNeighborsAtVertex_old(MRIS *mris, int vno, int nlinks, size_t listCapacity, int* vlist, int* hops)
{
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];    
  VERTEX          * const v  = &mris->vertices         [vno];

  int m, n, vtotal = 0, link_dist, ring_total;
  if (v->ripflag) return (0);

  if (1) {
    static bool laterTime = false;
    if (!laterTime) { laterTime = true;
      fprintf(stdout, "%s:%d checking marks initially clear\n", __FILE__, __LINE__); 
    }
    for (n = 0; n < mris->nvertices; n++) {
      cheapAssert(mris->vertices[n].marked == 0);
    }
  }
  
  v->marked = -1;
  for (n = 0; n < vt->vtotal; n++) {
    vlist[n] = vt->v[n];
    mris->vertices[vt->v[n]].marked = n < vt->vnum ? 1 : (n < vt->v2num ? 2 : 3);
  }

  if (nlinks <= mris->nsize) {
    cheapAssert(nlinks <= vt->nsizeCur);
    switch (nlinks) {
      case 1:
        vtotal = vt->vnum;
        break;
      case 2:
        vtotal = vt->v2num;
        break;
      case 3:
        vt->vtotal = vt->v3num;
        break;
      default:
        vtotal = 0;
        ErrorExit(ERROR_BADPARM, "MRISfindNeighborsAtVertex: nlinks=%d invalid", nlinks);
        break;
    }
  }
  else  // bigger than biggest neighborhood held at each vertex
  {
    link_dist = vt->nsizeCur;
    vtotal    = vt->vtotal;
    // at each iteration mark one more ring with the ring distance
    while (link_dist < nlinks) {
      link_dist++;
      ring_total = 0;
      for (n = 0; n < vtotal; n++) {
        VERTEX_TOPOLOGY const * const vnt = &mris->vertices_topology[vlist[n]];    
        VERTEX          const * const vn  = &mris->vertices         [vlist[n]];
        if (vn->ripflag) continue;
        for (m = 0; m < vnt->vnum; m++)  // one more ring out
        {
          if (mris->vertices[vnt->v[m]].marked == 0) {
            vlist[vtotal + ring_total] = vnt->v[m];
            mris->vertices[vnt->v[m]].marked = link_dist;
            ring_total++;
          }
        }
      }
      vtotal += ring_total;
    }
  }
  
  for (n = 0; n < vtotal; n++) {
    hops[n] = mris->vertices[vlist[n]].marked;
    mris->vertices[vlist[n]].marked = 0;
  }
  mris->vertices[vno].marked = 0;
  
  return (vtotal);
}


static int mrisInitializeNeighborhood(MRI_SURFACE *mris, int vno)
{
  int vtmp[MAX_NEIGHBORS], vnum, i, j, n, neighbors, nsize;

  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
  VERTEX          * const v  = &mris->vertices         [vno];
  if (vno == Gdiag_no) {
    DiagBreak();
  }

  vt->nsizeCur = 1;
  vt->nsizeMax = 1;
  vt->vtotal   = vt->vnum;
  
  if (v->ripflag || !vt->vnum) {
    return (ERROR_BADPARM);
  }
  
  cheapAssert(vt->vnum < MAX_NEIGHBORS);
  memmove(vtmp, vt->v, vt->vnum * sizeof(int));

  /* mark center so not included */
  v->marked = 1;

  /* mark 1-neighbors so we don't count them twice */
  vnum = neighbors = vt->vnum;

  cheapAssert(mris->nsize > 0);

  for (nsize = 2; nsize <= mris->nsize; nsize++) {
    /* mark all current neighbors */
    vnum = neighbors; /* neighbors will be incremented during loop */
    for (i = 0; i < neighbors; i++) {
      mris->vertices[vtmp[i]].marked = 1;
    }
    /* look at the candidates for the next ring out */
    for (i = 0; neighbors < MAX_NEIGHBORS && i < vnum; i++) {
      n = vtmp[i];
      VERTEX_TOPOLOGY const * const vnbt = &mris->vertices_topology[n];
      VERTEX          const * const vnb  = &mris->vertices         [n];
      if (vnb->ripflag) {
        continue;
      }

      for (j = 0; j < vnbt->vnum; j++) {
        VERTEX * const vnb2 = &mris->vertices[vnbt->v[j]];
        if (vnb2->ripflag || vnb2->marked) {
          continue;
        }
        vtmp[neighbors] = vnbt->v[j];
        vnb2->marked = 1;
        if (++neighbors >= MAX_NEIGHBORS) {
          fprintf(stdout, "vertex %d has too many neighbors!\n", vno);
          break;
        }
      }
    }

    // Fill in the next layer's details
    switch (nsize) {
      case 2:
        vt->v2num = neighbors;
        break;
      case 3:
        vt->v3num = neighbors;
        break;
      default:
        cheapAssert(false);
        break;
    }
    vt->nsizeMax = nsize;
    vt->nsizeCur = nsize;
    vt->vtotal   = neighbors;
  }
  /*
    now reallocate the v->v structure and place the 2-connected neighbors
    sequentially after the 1-connected neighbors.
  */
  free(vt->v);
  vt->v = (int *)calloc(neighbors, sizeof(int));
  if (!vt->v)
    ErrorExit(ERROR_NO_MEMORY,
              "MRISsetNeighborhoodSize: could not allocate list of %d "
              "nbrs at v=%d",
              neighbors,
              vno);

  v->marked = 0;
  for (n = 0; n < neighbors; n++) {
    vt->v[n] = vtmp[n];
    mris->vertices[vtmp[n]].marked = 0;
  }

  for (n = 0; n < neighbors; n++)
    for (i = 0; i < neighbors; i++)
      if (i != n && vt->v[i] == vt->v[n])
        fprintf(stdout, "warning: vertex %d has duplicate neighbors %d and %d!\n", vno, i, n);

  if ((vno == Gdiag_no) && (Gdiag & DIAG_SHOW) && DIAG_VERBOSE_ON) {
    fprintf(stdout, "v %d: vnum=%d, v2num=%d, vtotal=%d\n", vno, vt->vnum, vt->v2num, vt->vtotal);
    for (n = 0; n < neighbors; n++) {
      fprintf(stdout, "v[%d] = %d\n", n, vt->v[n]);
    }
  }

  return (NO_ERROR);
}


// Another implementation is found in
// utils/mrisurf_vals.c:int MRISsampleDistances(MRI_SURFACE *mris, int *nbrs, int max_nbhd) {



//=============================================================================
// Faces
//
bool mrisCheckVertexFaceTopologyWkr(const char* file, int line, MRIS const * mris, bool always) {

  if (!always && !spendTimeCheckingForLogicProblem(file, line)) return true;

  int reported = 0;
  #define shouldReport (shouldReportWkr(__LINE__) && (reported=1))    // sets reported only if should report
  
  if (!mrisCheckVertexVertexTopologyWkr(file,line,mris,true)) reported = 1;
  
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE const * const f = &mris->faces[fno];

    int prevVno = f->v[VERTICES_PER_FACE-1];
    int n;
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      int const vno = f->v[n];
      if (f->ripflag) continue;
      if ((f->v[0]|f->v[1]|f->v[2]) == 0) continue; // HACK this is what an unattached face looks like
      
      VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

      int i = -1;

      // The vertex points to the face exactly once
      //
      if (v->num && !v->f) {
        if (shouldReport) { 
          fprintf(stdout, "nullptr in [vno:%d].f when num:%d\n", vno, v->num);
        }
      } else {
        int iTrial;
        for (iTrial = 0; iTrial < v->num; iTrial++) {
          if (v->f[iTrial] == fno) {
            if (i == v->num) {
              if (shouldReport) {
                fprintf(stdout, "fno:%d found twice in [vno:%d].f[i:%d && iTrial:%d]\n", fno, vno, i, iTrial);
              }
            }
            i = iTrial;
          }
        }
      }
      if (i < 0) {
        if (shouldReport) {
          fprintf(stdout, "fno:%d not found in [vno:%d].f[*]\n", fno, vno);
        }
      }

      if (!v->n) {
        if (shouldReport) {
          fprintf(stdout, "nullptr in [vno:%d].n\n", vno);
        }
      } else {
        if (v->n[i] != n) {
          if (shouldReport) {
            fprintf(stdout, "[fno:%d].v[n:%d] holds vno:%d but [vno:%d].n[i:%d]:%d != n:%d\n", 
              fno, n, vno, vno, i, v->n[i], n);
          }
        }
      }
            
      // The vertices are neighbours
      //
      if (!mrisVerticesAreNeighbors(mris, vno, prevVno)) {
        if (shouldReport) {
          fprintf(stdout, "[fno:%d] holds adjacent vno:%d and vno:%d but they are not neighbours\n", 
            fno, vno, prevVno);
        }
      }
      
      prevVno = vno;
    }
  }
  
  return reported == 0;
}


int MRISisSurfaceValid(MRIS *mris, int patch, int verbose)
{
  int n, m, p, q, vnop, nfound, mark;
  int euler, nedges, nvf, nffm, nffs, cf, nvc;
  FACE *fm;

  // check if every vertex has the same number of vertices and faces
  if (verbose == 2) {
    fprintf(stderr, "\nchecking for border faces...\n");
  }
  nvf = 0;
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vn = &mris->vertices_topology[n];
    if (vn->vnum != vn->num) {
      nvf++;
      if (verbose == 2 && patch == 0) {
        fprintf(stderr, "vertex %d : vnum = %d num=%d\n", n, vn->vnum, vn->num);
      }
    }
  }

  // check if some faces have more than 2 neighbors
  if (verbose == 2) {
    fprintf(stderr, "\nchecking for single or multiple faces...\n");
  }
  nffm = nffs = 0;
  nedges = 0;
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vn = &mris->vertices_topology[n];
    nedges += vn->vnum;
    for (p = 0; p < vn->vnum; p++) {
      vnop = vn->v[p];
      /* counts the # of faces that contain the edge n<-->vnop */
      cf = 0;
      for (m = 0; m < vn->num; m++) {
        fm = &mris->faces[vn->f[m]];
        // check if fm contains vnop
        for (q = 0; q < 3; q++) {
          if (fm->v[q] == vnop) {
            cf++;
            break;
          }
        }
      }
      if (cf < 2) {
        nffs++;
      }
      if (cf > 2) {
        nffm++;
      }
      if (verbose == 2) {
        if (cf < 2 && patch == 0) {
          fprintf(stderr, "edge %d <--> %d : single face\n", n, vnop);
        }
        if (cf > 2) {
          fprintf(stderr, "edge %d <--> %d : multiple face\n", n, vnop);
        }
      }
    }
  }
  nedges /= 2;
  euler = mris->nvertices - nedges + mris->nfaces;

  // finally check the # of connected neighbors of each vertex
  if (verbose == 2) {
    fprintf(stderr, "\nchecking for corner configurations...\n");
  }
  nvc = 0;
  for (n = 0; n < mris->nvertices; n++) {
    VERTEX_TOPOLOGY const * const vn = &mris->vertices_topology[n];
    for (p = 0; p < vn->vnum; p++) {
      vnop = vn->v[p];
      mris->vertices[vnop].marked = 0;
    }
    // mark first vertex
    vnop = vn->v[0];
    mris->vertices[vnop].marked = 1;
    // find connected neighbors
    nfound = 1;
    while (nfound) {
      nfound = 0;
      for (m = 0; m < vn->num; m++) {
        fm = &mris->faces[vn->f[m]];
        // check if fm contains a marked vertex
        mark = 0;
        for (q = 0; q < 3; q++) {
          if (fm->v[q] == n) {
            continue;
          }
          if (mris->vertices[fm->v[q]].marked) {
            mark = 1;
            break;
          }
        }
        if (mark) {
          for (q = 0; q < 3; q++) {
            if (fm->v[q] == n) {
              continue;
            }
            if (mris->vertices[fm->v[q]].marked) {
              continue;
            }
            mris->vertices[fm->v[q]].marked = 1;
            nfound = 1;
            break;
          }
        }
      }
    }
    // now check if there is a remaining vertex that is not marked
    for (p = 0; p < vn->vnum; p++) {
      vnop = vn->v[p];
      if (mris->vertices[vnop].marked == 0) {
        nvc++;
        if (verbose == 2) {
          fprintf(stderr, "vertex %d : corner configuration\n", n);
        }
        break;
      }
    }
  }
  if (verbose)
    fprintf(stderr,
            "Surface Diagnostics:       \n"
            "   eno=%d (nv=%d, nf=%d, ne=%d)\n"
            "   # of border vertices [ #v ~ #f ]     %d \n"
            "   # of edges with single face          %d \n"
            "   # of edges with more than 2 faces    %d \n"
            "   # of corner configurations           %d\n",
            euler,
            mris->nvertices,
            mris->nfaces,
            nedges,
            nvf,
            nffs,
            nffm,
            nvc);

  if (patch && (nffm || nvc)) {
    return 0;
  }
  else if (nvf || nffs || nffm || nvc) {
    return 0;
  }
  return 1;
}


int mrisValidFaces(MRIS *mris)
{
  int fno, nfaces, nvalid;

  nfaces = mris->nfaces;
  for (fno = nvalid = 0; fno < nfaces; fno++)
    if (!mris->faces[fno].ripflag) {
      nvalid++;
    }

  return (nvalid);
}


int mrisCountAttachedFaces(MRIS* mris, int vno0, int vno1) {
  int count = 0;
  
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno0];

  int n;
  for (n = 0; n < vt->num; n++) {                   // For each face attached to vno0
    int const fno = vt->f[n];
    FACE const * const f = &mris->faces[fno];
    int m;
    for (m = 0; m < VERTICES_PER_FACE; m++) {       // For each vertex of the face
      if (f->v[m] == vno1) count++;                 // Count the faces sharing this edge
    }
  }

  return count;
}


bool mrisCanAttachFaceToVertices(MRIS* mris, int vno0, int vno1, int vno2) {
  
  // Should not even be contemplating this!
  //
  cheapAssert(vno0 != vno1 && vno0 != vno2 && vno1 != vno2);
  
  // None of the three edges can have two faces already using them
  //
  return 
    mrisCountAttachedFaces(mris,vno0,vno1) < 2 &&
    mrisCountAttachedFaces(mris,vno0,vno2) < 2 &&
    mrisCountAttachedFaces(mris,vno1,vno2) < 2;
}


// Return the index of the triangular face vno0-->vno1-->vno2, or -1
//
int findFace(MRIS *mris, int vno0, int vno1, int vno2)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno0];

  int n;
  for (n = 0; n < vt->num; n++) {       // For each face attached to vno0
    int const fno = vt->f[n];
    
    FACE *f = &mris->faces[fno];        // For each vertex of the face
    int n1;
    for (n1 = 0; n1 < VERTICES_PER_FACE; n1++) {
      int const vno = f->v[n1];

      // If the vertex is not in the triangle, reject this face
      //
      if (vno != vno0 && vno != vno1 && vno != vno2) goto NotMatch;
    }
    
    return (fno);   // not rejected
    
NotMatch:;
  }

  return -1;        // no matching face
}


bool isFace(MRIS *mris, int vno0, int vno1, int vno2)
{
  return findFace(mris, vno0, vno1, vno2) >= 0;
}


void mrisSetVertexFaceIndex(MRIS *mris, int vno, int fno)
  // HACK - external usage of this should be eliminated!
{
  FACE const *      const f  = &mris->faces[fno];
  VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];

  int n;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    if (f->v[n] == vno) break;
  }
  cheapAssert(n < VERTICES_PER_FACE);

  int i;
  for (i = 0; i < vt->num; i++)
    if (vt->f[i] == fno) {
      vt->n[i] = n;
    }
}


int mrisVertexFaceIndex(MRIS *mris, int vno, int fno) {
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];
  int i;
  for (i = 0; i < v->num; i++) {
    if (v->f[i] == fno) return i;
  }
  return -1;
}

int mrisFaceVertexIndex(MRIS *mris, int fno, int vno) {
  FACE const * const face = &mris->faces[fno];
  int i;
  for (i = 0; i < VERTICES_PER_FACE; i++) {
    if (face->v[i] == vno) return i;
  }
  return -1;
}

int vertexInFace(MRIS *mris, int vno, int fno)
{
  VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[vno];

  int n;
  for (n = 0; n < v->num; n++)
    if (v->f[n] == fno) {
      return (1);
    }
  return (0);
}


int findOtherEdgeFace(MRIS const *mris, int fno, int vno, int vn1)
{
  int n, m;

  VERTEX_TOPOLOGY const * const v1 = &mris->vertices_topology[vno];
  VERTEX_TOPOLOGY const * const v2 = &mris->vertices_topology[vn1];
  for (n = 0; n < v1->num; n++) {
    if (v1->f[n] == fno) {
      continue;
    }
    for (m = 0; m < v2->num; m++)
      if (v1->f[n] == v2->f[m]) {
        return (v1->f[n]);
      }
  }

  // fprintf(WHICH_OUTPUT,"edge (%d<-->%d) does not have two faces\n",vno,vn1);

  return fno;
}


short FACE_vertexIndex_find(FACE *pFace, int avertex)
{
  //
  // PRECONDITIONS
  // o <avertex> denotes a vertex number to lookup in the <pFace>.
  //
  // POSTCONDITIONS
  // o The vertex index (0, 1, 2) containing the <avertex> is returned
  //   or -1 if not found.
  //
  int vertex = 0;
  int ret = -1;
  for (vertex = 0; vertex < VERTICES_PER_FACE; vertex++)
    if (pFace->v[vertex] == avertex) {
      ret = vertex;
    }
  return ret;
}


int FACE_vertexIndexAtMask_find(FACE *apFACE_I, VECTOR *apv_verticesCommon)
{
  //
  // PRECONDITIONS
  //  o Called after <apv_verticesCommon> has been processed by
  //    VERTICES_commonInFaces_find().
  //
  // POSTCONDITIONS
  //  o The vertex index in <apFACE_I> corresponding to the '-1'
  //    in <apv_verticesCommon> is returned. For two bordering faces
  //    that share an edge defined by <apv_verticesCommon>, this function
  //    determines the index of the vertex that is *not* on this shared
  //    edge.

  int face = 0;
  int vertex = 0;
  short b_inCommon = 0;
  int ret = -1;

  for (face = 0; face < 3; face++) {
    vertex = apFACE_I->v[face];
    b_inCommon = VECTOR_elementIndex_find(apv_verticesCommon, (float)vertex);
    if (!b_inCommon) {
      ret = vertex;
      break;
    }
  }

  return ret;
}

short VERTICES_commonInFaces_find(FACE *apFACE_I, FACE *apFACE_J, VECTOR *apv_verticesCommon)
{
  //
  // PRECONDITIONS
  //  o The <apFACE>s must be triangles with 3 vertices each.
  //  o It is assumed (but not mandatory) that the FACES share
  //    at least one common vertex - or more often a common
  //    edge, i.e. two common vertices.
  //  o The <apv_uncommon> VECTOR's memory must be managed by the
  //    caller, i.e. created and freed, and should be a 3x1 VECTOR.
  //
  // POSTCONDITIONS
  //  o The number of vertices that are in common between the two
  //    faces are returned in the function name. This can be either
  //    (0, 1, 2, 3).
  //  o The indices of the common vertices are returned in
  //    apv_verticesCommon. This is a three element vector, with each
  //    element corresponding to a common vertex. Vertex indices that
  //    are not common between the two faces have a -1.
  //

  int i = 0;
  int j = 0;
  int k = 0;
  const char *pch_function = "VERTICES_commonInFaces_find";
  short b_hit = 0;
  float f_val = -1.;

  if (apv_verticesCommon->rows != 3 || apv_verticesCommon->cols != 1) {
    ErrorExit(-1, "%s: Return VECTOR must be 3x1.\n", pch_function);
  }
  V3_LOAD(apv_verticesCommon, -1, -1, -1);

  for (i = 0; i < 3; i++) {
    b_hit = 0;
    f_val = -1.;
    for (j = 0; j < 3; j++) {
      if (apFACE_J->v[j] == apFACE_I->v[i]) {
        b_hit = 1;
      }
    }
    if (b_hit) {
      f_val = apFACE_I->v[i];
      VECTOR_ELT(apv_verticesCommon, ++k) = f_val;
    }
  }
  return k;
}


static void mrisAddFaceToVertex(MRIS *mris, int vno, int fno, int n) {
  costlyAssert(mrisVertexFaceIndex(mris, vno, fno) < 0);     // check not already added
  FACE *            const f = &mris->faces[fno];
  VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno];
  f->v[n] = vno;
  int const i = v->num++;
  v->f = (int  *)realloc(v->f, v->num*sizeof(int));
  v->n = (uchar*)realloc(v->n, v->num*sizeof(int));
  v->f[i] = fno;
  v->n[i] = n;
}

static void mrisAttachFaceWkr(MRIS* mris, int fno, int vno0, int vno1, int vno2, bool edgesMustExist) {
  cheapAssertValidFno(mris,fno);
  cheapAssertValidVno(mris,vno0);
  cheapAssertValidVno(mris,vno1);
  cheapAssertValidVno(mris,vno2);
 
  cheapAssert(mrisCanAttachFaceToVertices(mris, vno0, vno1, vno2));
    //
    // This assertion was seen when a triangular tesselation written as a quad file did so by creating quad with two identical vertices
    // which then might be read back as two triangles, one of which has identical vertices!  This would ruin euler
    // calculations as well as create huge numbers of zero-area badly-defined-norm FACE.

  int vno[4]; vno[0] = vno0; vno[1] = vno1;  vno[2] = vno2;  vno[3] = vno0;
  if (edgesMustExist) {
    int i; 
    for (i = 0; i < 3; i++) {
      costlyAssert(mrisVerticesAreNeighbors(mris, vno[i], vno[i+1]));
    }
  } else {
    int i; 
    for (i = 0; i < 3; i++) {
      if (!mrisVerticesAreNeighbors(mris, vno[i], vno[i+1]))
           mrisAddEdgeWkr          (mris, vno[i], vno[i+1]);
    }
  }
  
  FACE * const f = &mris->faces[fno];
  if ((f->v[0]|f->v[1]|f->v[2]) == 0) {  // not currently attached
    int i;
    for (i = 0; i < 3; i++)
      mrisAddFaceToVertex(mris, vno[i], fno, i);
  } else {
    // assume whomever attached it got it right
    // That will be checked soon by mrisCheckVertexFaceTopology
  }
}


void setFaceAttachmentDeferred(MRIS* mris, bool to) {
  // currently mostly NYI until seen to be a performance problem
  if (!to) cheapAssert(mrisCheckVertexFaceTopology(mris));
}

void mrisAttachFaceToEdges   (MRIS* mris, int fno, int vno1, int vno2, int vno3) {
  mrisAttachFaceWkr(mris, fno, vno1, vno2, vno3, true);
}

void mrisAttachFaceToVertices(MRIS* mris, int fno, int vno1, int vno2, int vno3) {
  //
  // This is the preferred way to build a surface.
  // Create the vertices, attach the faces, and let this code create the necessary edges
  //
  mrisAttachFaceWkr(mris, fno, vno1, vno2, vno3, false);
}

int mrisRemoveFace(MRIS *mris, int fno);
  
int mrisRemoveLink(MRIS *mris, int vno1, int vno2)
{
  FACE *face;
  int vno, fno, nvalid;

  VERTEX_TOPOLOGY * const v1t = &mris->vertices_topology[vno1];
  VERTEX_TOPOLOGY * const v2t = &mris->vertices_topology[vno2];
  VERTEX          * const v1  = &mris->vertices         [vno1];
  VERTEX          * const v2  = &mris->vertices         [vno2];

  mrisRemoveEdge(mris, vno1, vno2);
  mrisRemoveEdge(mris, vno2, vno1);

  /* now remove all the faces which contain both edges */
  for (fno = 0; fno < v1t->num; fno++) {
    face = &mris->faces[v1t->f[fno]];
    for (vno = 0; vno < VERTICES_PER_FACE; vno++)
      if (face->v[vno] == vno2) /* found a face with both vertices */
      {
        face->ripflag = 1;
      }
    if (face->ripflag  >= 2) {
      mrisRemoveFace(mris, v1t->f[fno]);
    }
  }

  /* make sure vno1 and vno2 are still part of at least 1 valid face */
  for (nvalid = fno = 0; fno < v1t->num; fno++)
    if (!mris->faces[v1t->f[fno]].ripflag) {
      nvalid++;
    }
  if (nvalid <= 0) {
    v1->ripflag = 1;
  }

  for (nvalid = fno = 0; fno < v2t->num; fno++)
    if (!mris->faces[v2t->f[fno]].ripflag) {
      nvalid++;
    }
  if (nvalid <= 0) {
    v2->ripflag = 1;
  }

  return (NO_ERROR);
}


/*
  Contructs an MRIS instance from a set of vertices and faces.

  \param vertices An (nvertices x 3) array of vertex positions.
  \param nvertices Number of vertices in mesh.
  \param faces An (nfaces x 3) array of face indices.
  \param nfaces Number of faces in mesh.
*/
MRIS* MRISfromVerticesAndFaces(const float *vertices, int nvertices, const int *faces, int nfaces)
{
  MRIS* mris = MRISalloc(nvertices, nfaces);
  mris->type = MRIS_TRIANGULAR_SURFACE;

  // vertex positions
  const float* v = vertices;
  for (int n = 0 ; n < mris->nvertices ; n++, v += 3) MRISsetXYZ(mris, n, v[0], v[1], v[2]);

  // face positions
  const int* f = faces;
  setFaceAttachmentDeferred(mris, true);
  for (int n = 0 ; n < mris->nfaces ; n++, f += 3) {

    int vno[4];
    vno[0] = f[0];
    vno[1] = f[1];
    vno[2] = f[2];
    vno[3] = f[0];

    for (int i = 0; i < 3; i++) {
      if (!mrisVerticesAreNeighbors(mris, vno[i], vno[i+1])) mrisAddEdgeWkr(mris, vno[i], vno[i+1]);
    }

    FACE * const f = &mris->faces[n];
    if ((f->v[0]|f->v[1]|f->v[2]) == 0) {
      for (int i = 0; i < 3; i++) mrisAddFaceToVertex(mris, vno[i], n, i);
    }
  }
  setFaceAttachmentDeferred(mris, false);

  mrisCompleteTopology(mris);
  MRIScomputeMetricProperties(mris);

  return mris;
}


#define MAX_VERTEX_NEIGHBORS 50
#define MAX_FACES 50

static void mrisDivideFace(MRIS *mris, int fno, int vno1, int vno2, int vnew_no);

int mrisDivideEdgeTopologically(MRIS * const mris, int const vno1, int const vno2)
{
  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "dividing edge %d --> %d\n", vno1, vno2);
  }

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || mris->nvertices == Gdiag_no) {
    printf("dividing edge %d --> %d, adding vertex number %d\n", vno1, vno2, mris->nvertices);
    DiagBreak();
  }
  VERTEX_TOPOLOGY const * const v1t = &mris->vertices_topology[vno1];
  VERTEX          const * const v1  = &mris->vertices         [vno1];
  VERTEX_TOPOLOGY const * const v2t = &mris->vertices_topology[vno2];
  VERTEX          const * const v2  = &mris->vertices         [vno2];

  if (mris->nfaces >= mris->max_faces-1)
    DiagBreak() ;
  cheapAssert(!(v1->ripflag || v2->ripflag || mris->nvertices >= mris->max_vertices || mris->nfaces >= (mris->max_faces - 1)));

  /* check to make sure these vertices or the faces they are part of
     have enough room to expand.
  */
  cheapAssert(
       !(   v1t->vnum >= MAX_VERTEX_NEIGHBORS || v2t->vnum >= MAX_VERTEX_NEIGHBORS 
         || v1t->num  >= MAX_FACES            || v2t->num  >= MAX_FACES
       ));

  /* add 1 new vertex, 2 new faces, and 2 new edges */
  int const vnew_no = mris->nvertices;
  MRISgrowNVertices(mris, mris->nvertices+1);

  VERTEX_TOPOLOGY * const vnewt = &mris->vertices_topology[vnew_no];
  
  modVnum(mris,vnew_no,2,true); /* at least connected to two bisected vertices */

  /* count the # of faces that both vertices are part of */
  int n, m, n1, n2;

  int flist[100];
  for (n = 0; n < v1t->num; n++) {
    int const fno = v1t->f[n];
    FACE const * const face = &mris->faces[fno];
    for (m = 0; m < VERTICES_PER_FACE; m++)
      if (face->v[m] == vno2) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
          fprintf(stdout, " face %d shared.\n", fno);
        }
        flist[vnewt->num++] = fno;
        if (vnewt->num == 100) {
          ErrorExit(ERROR_BADPARM, "Too many faces to divide edge");
        }
        addVnum(mris,vnew_no,1);
        vnewt->vtotal = vnewt->vnum;
      }
  }

  /* will be part of two new faces also */
  // total array size is going to be vnew->num*2!
  if (vnewt->num >= 50) {
    ErrorExit(ERROR_BADPARM, "Too many faces to divide edge");
  }
  for (n = 0; n < vnewt->num; n++) {
    flist[vnewt->num + n] = mris->nfaces + n;
  }
  vnewt->num *= 2;
  vnewt->f = (int *)calloc(vnewt->num, sizeof(int));
  if (!vnewt->f) {
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth face list.\n", vnew_no);
  }
  vnewt->n = (uchar *)calloc(vnewt->num, sizeof(uchar));
  if (!vnewt->n) {
    ErrorExit(ERROR_NOMEMORY, "could not allocate %dth face list.\n", vnew_no);
  }
  vnewt->v = (int *)calloc(vnewt->vnum, sizeof(int));
  if (!vnewt->v) ErrorExit(ERROR_NOMEMORY, "could not allocate %dth vertex list.\n", vnew_no);

  vnewt->num = 0; clearVnum(mris,vnew_no);

  /* divide every face that both vertices are part of in two */
  for (n = 0; n < v1t->num; n++) {
    int const fno = v1t->f[n];
    FACE const * const face = &mris->faces[fno];
    for (m = 0; m < VERTICES_PER_FACE; m++)
      if (face->v[m] == vno2) {
        if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
          fprintf(stdout, "dividing face %d along edge %d-->%d.\n", fno, vno1, vno2);
        if (face->v[m] == Gdiag_no || vno2 == Gdiag_no) {
          DiagBreak();
        }
        mrisDivideFace(mris, fno, vno1, vno2, vnew_no);
      }
  }

  /* build vnew->f and vnew->n lists by going through all faces v1 and
     v2 are part of */
  int fno;
  for (fno = 0; fno < vnewt->num; fno++) {
    vnewt->f[fno] = flist[fno];
    FACE const * const face = &mris->faces[flist[fno]];
    for (n = 0; n < VERTICES_PER_FACE; n++)
      if (face->v[n] == vnew_no) {
        vnewt->n[fno] = (uchar)n;
      }
  }

  /* remove vno1 from vno2 list and visa-versa */
  for (n = 0; n < v1t->vnum; n++)
    if (v1t->v[n] == vno2) {
      v1t->v[n] = vnew_no;
      break;
    }
  for (n = 0; n < v2t->vnum; n++)
    if (v2t->v[n] == vno1) {
      v2t->v[n] = vnew_no;
      break;
    }
  /* build vnew->v list by going through faces it is part of and
     rejecting duplicates
  */
  for (fno = 0; fno < vnewt->num; fno++) {
    FACE const * const face = &mris->faces[vnewt->f[fno]];
    n1 = vnewt->n[fno] == 0 ? VERTICES_PER_FACE - 1 : vnewt->n[fno] - 1;
    n2 = vnewt->n[fno] == VERTICES_PER_FACE - 1 ? 0 : vnewt->n[fno] + 1;
    int vnoA = face->v[n1];
    int vnoB = face->v[n2];

    /* go through this faces vertices and see if they should be added to v[] */
    for (n = 0; n < vnewt->vnum; n++) {
      if (vnewt->v[n] == vnoA) {
        vnoA = -1;
      }
      if (vnewt->v[n] == vnoB) {
        vnoB = -1;
      }
    }
    if (vnoA >= 0) {
      vnewt->v[vnumAdd(mris,vnew_no,1)] = vnoA;
    }
    if (vnoB >= 0) {
      vnewt->v[vnumAdd(mris,vnew_no,1)] = vnoB;
    }
    vnewt->vtotal = vnewt->vnum;
  }
  if (0 && Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "%d edges and %d faces.\n", vnewt->vnum, vnewt->num);
  }

  if (!vnewt->vnum || !v1t->vnum || !v2t->vnum) {
    fprintf(stderr, "empty vertex (%d <-- %d --> %d!\n", vno1, vnew_no, vno2);
    DiagBreak();
  }
  if (vnewt->vnum != 4 || vnewt->num != 4) {
    DiagBreak();
  }

  mrisInitializeNeighborhood(mris, vnew_no);
  
  return vnew_no;
}


static void mrisDivideFace(MRIS *mris, int fno, int vno1, int vno2, int vnew_no)
{
  int fnew_no, n, vno3, flist[5000], vlist[5000], nlist[5000];

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || vnew_no == Gdiag_no) {
    DiagBreak();
  }

  /* divide this face in two, reusing one of the face indices, and allocating
     one new one
  */
  cheapAssert(mris->nfaces < mris->max_faces);

  fnew_no = mris->nfaces;
  MRISgrowNFaces(mris, fnew_no+1);

  FACE * const f1 = &mris->faces[fno];
  FACE * const f2 = &mris->faces[fnew_no];
  VERTEX_TOPOLOGY const * const v2   = &mris->vertices_topology[vno2];
  VERTEX_TOPOLOGY       * const vnew = &mris->vertices_topology[vnew_no];
  memmove(f2->v, f1->v, VERTICES_PER_FACE * sizeof(int));

  /* set v3 to be other vertex in face being divided */

  /* 1st construct f1 by replacing vno2 with vnew_no */
  for (vno3 = -1, n = 0; n < VERTICES_PER_FACE; n++) {
    if (f1->v[n] == vno2) /* replace it with vnew */
    {
      f1->v[n] = vnew_no;
      vnew->f[vnew->num] = fno;
      vnew->n[vnew->num++] = (uchar)n;
    }
    else if (f1->v[n] != vno1) {
      vno3 = f1->v[n];
    }
  }
  
  VERTEX_TOPOLOGY * const v3 = &mris->vertices_topology[vno3];

  if (vno1 == Gdiag_no || vno2 == Gdiag_no || vno3 == Gdiag_no) {
    DiagBreak();
  }

  /* now construct f2 */

  /*  replace vno1 with vnew_no in f2 */
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    if (f2->v[n] == vno1) /* replace it with vnew */
    {
      f2->v[n] = vnew_no;
      vnew->f[vnew->num] = fnew_no;
      vnew->n[vnew->num++] = (uchar)n;
    }
  }

  /* now replace f1 in vno2 with f2 */
  for (n = 0; n < v2->num; n++)
    if (v2->f[n] == fno) {
      v2->f[n] = fnew_no;
    }

  /* add new face and edge connected to new vertex to v3 */
  memmove(flist, v3->f, v3->num  * sizeof(v3->f[0]));
  memmove(vlist, v3->v, v3->vnum * sizeof(v3->v[0]));
  memmove(nlist, v3->n, v3->num  * sizeof(v3->n[0]));
  free(v3->f);
  free(v3->v);
  free(v3->n);
  v3->v = (int *)calloc(v3->vnum + 1, sizeof(int));
  if (!v3->v) ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d vertices", v3->vnum);
  v3->f = (int *)calloc(v3->num + 1, sizeof(int));
  if (!v3->f) ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d faces", v3->num);
  v3->n = (uchar *)calloc(v3->num + 1, sizeof(uchar));
  if (!v3->n) ErrorExit(ERROR_NO_MEMORY, "mrisDivideFace: could not allocate %d nbrs", v3->n);
  memmove(v3->f, flist, v3->num * sizeof(v3->f[0]));
  memmove(v3->n, nlist, v3->num * sizeof(v3->n[0]));
  memmove(v3->v, vlist, v3->vnum * sizeof(v3->v[0]));
  v3->v[vnumAdd(mris,vno3,1)] = vnew_no;
  v3->vtotal = v3->vnum;
  v3->f[v3->num] = fnew_no;

  /*  find position of v3 in new face f2 */
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    if (f2->v[n] == vno3) {
      v3->n[v3->num] = n;
      break;
    }
  }
  if (n >= VERTICES_PER_FACE) fprintf(stderr, "could not find v3 (%d) in new face %d!\n", vno3, fnew_no);
  v3->num++;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "face %d: (%d, %d, %d)\n", fno, f1->v[0], f1->v[1], f1->v[2]);
    fprintf(stdout, "face %d: (%d, %d, %d)\n", fnew_no, f2->v[0], f2->v[1], f2->v[2]);
  }

  // MRISfindNeighborsAtVertex needs to be called on all the vertices within some extended neighborhood of the added vertex
  //
  mrisInitializeNeighborhood(mris, vno3);
}

static void mrisCompleteTopology_old(MRI_SURFACE *mris);
static void mrisCompleteTopology_new(MRI_SURFACE *mris);

void mrisCompleteTopology(MRI_SURFACE *mris) {
  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d mrisCompleteTopology ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
  }

  static bool laterTime, use_new;
  if (!laterTime) {
    laterTime = true;
    use_new = !!getenv("mrisCompleteTopology_new");
  }
  
  if (use_new) mrisCompleteTopology_new(mris);
  else         mrisCompleteTopology_old(mris);

  mrisCheckVertexFaceTopology(mris);

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d mrisCompleteTopology ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
    mrisDumpShape(stdout, mris);
  }
}

static void mrisCompleteTopology_new(MRI_SURFACE *mris)
{
  setFaceAttachmentDeferred(mris,true);
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE const * face = &mris->faces[fno];
    mrisAttachFaceToVertices(mris, fno, face->v[0], face->v[1], face->v[2]); 
  }
  setFaceAttachmentDeferred(mris,false);

  int ntotal = 0, vtotal = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    vtotal += vt->vtotal;
    ntotal++;
  }
  
  if (ntotal == 0) ntotal = 1;
  mris->avg_nbrs = (float)vtotal / (float)ntotal;

}

static void mrisCompleteTopology_old(MRI_SURFACE *mris) // was mrisFindNeighbors
{
  int n0, n1, i, k, m, n, vno, vtotal, ntotal, vtmp[MAX_NEIGHBORS];
  FACE *f;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "finding surface neighbors...");
  }

  for (k = 0; k < mris->nvertices; k++) {
    if (k == Gdiag_no) {
      DiagBreak();
    }
    VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[k];    
    clearVnum(mris,k);
    for (m = 0; m < vt->num; m++) {
      n = vt->n[m];               /* # of this vertex in the mth face that it is in */
      f = &mris->faces[vt->f[m]]; /* ptr to the mth face */
      /* index of vertex we are connected to */
      n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
      n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
      for (i = 0; i < vt->vnum && vtmp[i] != f->v[n0]; i++) {
        ;
      }
      if (i == vt->vnum) {
        vtmp[(int)vnumAdd(mris,k,1)] = f->v[n0];
      }
      for (i = 0; i < vt->vnum && vtmp[i] != f->v[n1]; i++) {
        ;
      }
      if (i == vt->vnum) {
        vtmp[(int)vnumAdd(mris,k,1)] = f->v[n1];
      }
    }
    if (mris->vertices_topology[k].v) {
      free(mris->vertices_topology[k].v);
    }
    mris->vertices_topology[k].v = (int *)calloc(mris->vertices_topology[k].vnum, sizeof(int));
    if (!mris->vertices_topology[k].v) ErrorExit(ERROR_NOMEMORY, "mrisFindNeighbors: could not allocate nbr array");

    vt->vtotal = vt->vnum;
    vt->nsizeMax = vt->nsizeCur = 1;
    for (i = 0; i < vt->vnum; i++) {
      vt->v[i] = vtmp[i];
    }

    /*
      if (vt->num != vt->vnum)
      printf("%d: num=%d vnum=%d\n",k,vt->num,vt->vnum);
    */
  }
  for (k = 0; k < mris->nfaces; k++) {
    f = &mris->faces[k];
    for (m = 0; m < VERTICES_PER_FACE; m++) {
      VERTEX_TOPOLOGY const * const v = &mris->vertices_topology[f->v[m]];
      for (i = 0; i < v->num && k != v->f[i]; i++) {
        ;
      }
      if (i == v->num) /* face has vertex, but vertex doesn't have face */
        ErrorExit(ERROR_BADPARM,
                  "mrisFindNeighbors: %s: face[%d].v[%d] = %d, "
                  "but face %d not in vertex %d "
                  "face list\n",
                  mris->fname,
                  k,
                  m,
                  f->v[m],
                  k,
                  f->v[m]);
    }
  }

  for (vno = ntotal = vtotal = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];    
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    vtotal += vt->vtotal;
    ntotal++;
  }

  mris->avg_nbrs = (float)vtotal / (float)ntotal;
}





void MRISsetRipInFacesWithRippedVertices(MRIS *mris)
{
  int n, k;

  for (k = 0; k < mris->nfaces; k++) {
    mris->faces[k].ripflag = FALSE;
  }
  
  for (k = 0; k < mris->nfaces; k++) {
    face_type *f = &mris->faces[k];
    for (n = 0; n < VERTICES_PER_FACE; n++)
      if (mris->vertices[f->v[n]].ripflag) {
        f->ripflag = TRUE;
      }
  }

  for (k = 0; k < mris->nvertices; k++) {
    mris->vertices[k].border = FALSE;
  }
  
  for (k = 0; k < mris->nfaces; k++)
    if (mris->faces[k].ripflag) {
      face_type *f = &mris->faces[k];
      for (n = 0; n < VERTICES_PER_FACE; n++) {
        mris->vertices[f->v[n]].border = TRUE;
      }
    }
}


/*-----------------------------------------------------
  remove this face, as well as any links that only exist through this face.
  ------------------------------------------------------*/
int mrisRemoveFace(MRIS *mris, int fno)
{
  FACE* face = &mris->faces[fno];
  face->ripflag = 1;
  
  int vno;
  for (vno = 0; vno < VERTICES_PER_FACE; vno++) {
    int vno1 = face->v[vno];
    int vno2 = face->v[vno < VERTICES_PER_FACE - 1 ? vno + 1 : 0];
    if (!mrisCountValidLinks(mris, vno1, vno2)) {
      mrisRemoveEdge(mris, vno1, vno2);
      mrisRemoveEdge(mris, vno2, vno1);
    }
  }

  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Remove ripped vertices and faces from the v->v and the
  v->f arrays respectively.
  ------------------------------------------------------*/
static void removeRippedFaces(MRI_SURFACE *mris)
{
  int    vno, n, fno, *out_faces, out_fno, nfaces;
  FACE   *face;

  out_faces = (int *)calloc(mris->nfaces, sizeof(int)) ;
  nfaces = mris->nfaces ;
  for (out_fno = fno = 0; fno < mris->nfaces; fno++) 
  {
    face = &mris->faces[fno];
    if (fno == Gdiag_no)
      DiagBreak() ;
    if (face->ripflag) 
    {
      out_faces[fno] = -1 ;
      nfaces-- ;
    }
    else
    {
      out_faces[fno] = out_fno ;
      if (out_fno == Gdiag_no)
	DiagBreak() ;
      if (fno != out_fno) { // at least one compressed out already
	mris->faces                  [out_fno] = mris->faces                  [fno];    // should free memory here
	mris->faceNormCacheEntries   [out_fno] = mris->faceNormCacheEntries   [fno];
      }
      out_fno++ ;
    }
  }

  cheapAssert(nfaces == out_fno);
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    int num ;
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno];    
    num = v->num ; v->num = 0 ;
    for (n = 0 ; n < num ; n++)
    {
      int fno = out_faces[v->f[n]] ;
      if (fno == Gdiag_no)
	DiagBreak() ;
      if (fno < 0 || mris->faces[fno].ripflag == 1)
	continue ;

      v->f[v->num++] = out_faces[v->f[n]] ;
    }
  }

  MRISremovedFaces(mris, nfaces);

  free(out_faces) ;
}

static void removeRippedVertices(MRI_SURFACE *mris)
{
  int    vno, n, fno, *out_vnos, out_vno, nvertices;
  FACE   *face;

  out_vnos = (int *)calloc(mris->nvertices, sizeof(int)) ;
  nvertices = mris->nvertices ;
  for (out_vno = vno = 0; vno < mris->nvertices; vno++) 
  {
    VERTEX const * const v = &mris->vertices[vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag) 
    {
      nvertices-- ;
      out_vnos[vno] = -1 ;  // mark it as ripped - fixes boundary condition when coming to end of array
    }
    else
    {
      if (out_vno == Gdiag_no)
	DiagBreak() ;
      out_vnos[vno] = out_vno ;
      if (vno != out_vno)  // at least one compressed out already
	memcpy(mris->vertices+out_vno, mris->vertices+vno, sizeof(VERTEX)) ;    // GROSS HACK
      out_vno++ ;
    }
  }
  
  for (vno = 0 ; vno < nvertices ; vno++)
  {
    VERTEX_TOPOLOGY * const v = &mris->vertices_topology[vno] ;

    int vnum, v2num, v3num ;

    vnum = v->vnum ; v2num = v->v2num ; v3num = v->v3num ;
    v->v3num = v->v2num = 0 ;  clearVnum(mris,vno);
    for (n = 0 ; n < v->vtotal ; n++)
    {
      int vno2 = v->v[n] ;
      if (vno2 < 0 || mris->vertices[vno2].ripflag == 1)
	continue ;

      v->v[v->v3num++] = out_vnos[v->v[n]] ;
      if (n < v2num)
	v->v2num++ ;
      if (n < vnum)
	addVnum(mris,vno,1);
    }
    v->nsizeMax = v->nsizeCur;  // since the above loop only went to v->vtotal
    cheapAssert(mris->nsize <= v->nsizeMax);
    switch (mris->nsize) {
    default:
    case 1:
      v->vtotal = v->vnum; v->nsizeCur = 1;
      break;
    case 2:
      v->vtotal = v->v2num; v->nsizeCur = 2;
      break;
    case 3:
      v->vtotal = v->v3num; v->nsizeCur = 3;
      break;
    }
  }

  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      face->v[n] = out_vnos[face->v[n]] ;
  }

  MRISremovedVertices(mris, nvertices);

  free(out_vnos) ;
}


void MRISrenumberRemovingRippedFacesAndVertices(MRIS* mris) {

  cheapAssert(!mris->dist_alloced_flags);
  
  mrisCheckVertexFaceTopology(mris);

  removeRippedFaces(mris);
  removeRippedVertices(mris);

  mrisCheckVertexFaceTopology(mris);
}


/*-----------------------------------------------------
  Remove ripped vertices and faces from the 
  v->v and the v->f arrays
  ------------------------------------------------------*/
void MRISremoveRipped(MRIS *mris)
{
  float* distCache          = (float*)calloc(mris->nvertices,sizeof(float));
  float* distOrigCache      = (float*)calloc(mris->nvertices,sizeof(float));
  int* affectedVnosPlus2    = (int  *)calloc(mris->nvertices,sizeof(int));
    // 0 means unused
    // 1 means end of list
    // 2 or more is an entry in the list other than at the end
    
  int  headAffectedVnosPlus2  = 1;
  int  affectedSize           = 0;
  
  int  largestNsizeMax = 0;
  
  // For the non-ripped vno
  //    Calculate the largest known neighbourhood   (used later)
  // For all ripped vno
  //    Rip all the faces they are attached to, adding to those already marked
  //    Add to the affected neighbours list
  //
  {
    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX          const * const v  = &mris->vertices         [vno];

      if (!v->ripflag) {
        if (largestNsizeMax < vt->nsizeMax) largestNsizeMax = vt->nsizeMax;
        continue;
      }

      cheapAssert(vt->vtotal == VERTEXvnum(vt, vt->nsizeCur));
          // since this is what the following code will set to later

      int n;
      for (n = 0; n < vt->num; n++) {
        int const fno = vt->f[n];
        FACE * const face = &mris->faces[fno];
        face->ripflag = 1;
      }

      affectedVnosPlus2[vno] = headAffectedVnosPlus2;             // add to list
      headAffectedVnosPlus2  = vno+2;
      affectedSize++;
    }
  }
  int const startOfRippedTail = headAffectedVnosPlus2;
  
  // Grow the affected neighbours set to include all that are within largestNsizeMax
  //
  int nsize;
  int endOfList = 1;
  for (nsize = 0; nsize <= largestNsizeMax; nsize++) {
    int startOfList = headAffectedVnosPlus2;
    int vno, vnoPlus2;
    for (vnoPlus2 = headAffectedVnosPlus2; vnoPlus2 != endOfList; vnoPlus2 = affectedVnosPlus2[vno]) {
      vno = vnoPlus2 - 2;
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      if (affectedVnosPlus2[vno]) continue;                     // already noted, includes ripped
      int n;
      for (n = 0; n < vt->vnum; n++) {
        int const vnoNbr = vt->v[n];
        if (affectedVnosPlus2[vnoNbr]) continue;                // already noted, includes ripped
        affectedVnosPlus2[vnoNbr] = headAffectedVnosPlus2;      // add to list
        headAffectedVnosPlus2     = vnoNbr + 2;
        affectedSize++;
      }
    }
    endOfList = startOfList;
  }
  
  // Recompute the neighbourhoods of all the non-ripped vertices on the list
  // Clear out the list for reuse in the next step
  //
  { 
    int vlist[MAX_NEIGHBORS], hops[MAX_NEIGHBORS];

    int vnoPlus2 = headAffectedVnosPlus2; 
    headAffectedVnosPlus2 = 1;

    while (vnoPlus2 != startOfRippedTail) {
      int const vno = vnoPlus2 - 2;

      VERTEX_TOPOLOGY * const vt = &mris->vertices_topology[vno];
      VERTEX          * const v  = &mris->vertices         [vno];

      int* const old_vt_v     = vt->v;
      int  const old_nsizeMax = vt->nsizeMax;
      int  const old_nsizeCur = vt->nsizeCur;
      int  const old_vtotal   = vt->vtotal;
      int  const old_cached   = VERTEXvnum(vt, old_nsizeMax);
      int  const old_vno_beyond_nsizeMax =   // used for a consistency check
        (old_vtotal > old_cached) ? vt->v[old_cached] : -1;
      
      // Must preserve the dist and dist_orig for all the vtotal entries
      // even though the v list will get reordered.  The ripped nodes may
      // change the #hops to a neighbour from 2 to 4 so it is not adequate
      // to simply strip the ripped nodes out of the v list, because that
      // might keep the neighbour at the wrong distance.
      //
      // Also must preserve the non-ripped nodes after old_nsizeMax up to vtotal because
      // these are the sampled ones, and the same sample must be used to keep their
      // dist and distOrig valid.  These can not have moved closer because removing 
      // nodes just increases hops.
      //
      int i;
      if (v->dist)      for (i = 0; i < vt->vtotal; i++) distCache    [vt->v[i]] = v->dist     [i];
      if (v->dist_orig) for (i = 0; i < vt->vtotal; i++) distOrigCache[vt->v[i]] = v->dist_orig[i];
      
      MRISfindNeighborsAtVertex_new(mris, vno, old_nsizeMax, MAX_NEIGHBORS, vlist, hops, false);    
        // Note: this code copes with ripped and removes them, even out of the vnum portion!
        // Fortunately it doesn't realloc v->v
        
      MRIS_setNsizeCur(mris, vno, old_nsizeCur);
      
      if (old_vtotal > old_cached) {
        cheapAssert(old_vt_v == vt->v);
        cheapAssert(old_vno_beyond_nsizeMax == -1 || old_vno_beyond_nsizeMax == vt->v[old_cached]);
        
        int new_vtotal = VERTEXvnum(vt, old_nsizeMax);
      
        for (i = old_cached; i < old_vtotal; i++) {
          VERTEX const * const vnbr = &mris->vertices[vt->v[i]];
          if (!vnbr->ripflag) vt->v[new_vtotal++] = vt->v[i];
        }
        
        cheapAssert(old_vtotal >= new_vtotal);
        vt->vtotal = new_vtotal;
      }
      
      if (v->dist)      for (i = 0; i < vt->vtotal; i++) v->dist     [i] = distCache    [vt->v[i]];
      if (v->dist_orig) for (i = 0; i < vt->vtotal; i++) v->dist_orig[i] = distOrigCache[vt->v[i]];
      
      // Clear this entry
      //
      vnoPlus2 = affectedVnosPlus2[vno];
      affectedVnosPlus2[vno] = 0;
      affectedSize--;
    }
    
    while (vnoPlus2 != 1) {
      int const vno = vnoPlus2 - 2;
      vnoPlus2 = affectedVnosPlus2[vno];
      affectedVnosPlus2[vno] = 0;
      affectedSize--;
    }
  }
  
  // Above should have emptied the list
  //
  cheapAssert(affectedSize == 0);
  cheapAssert(headAffectedVnosPlus2 == 1);

  // For all the ripped faces
  //    create a set of their vertices, including their ripped ones (so don't need to read the VERTEX in this loop)
  //    delete all their knowledge of their vertices
  //
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE * const face = &mris->faces[fno];
    if (!face->ripflag) continue;
    int n;
    for (n = 0; n < VERTICES_PER_FACE; n++) {
      int vno = face->v[n];
      cheapAssertValidVno(mris,vno);
      
      face->v[n] = 0;                                       // should be -1 but uninit'ed faces have 0's
      
      if (affectedVnosPlus2[vno]) continue;                 // already noted

      affectedVnosPlus2[vno] = headAffectedVnosPlus2;       // add to list
      headAffectedVnosPlus2 = vno + 2;
      affectedSize++;
    }
  }
  cheapAssert(affectedSize <= mris->nvertices);

  // Fix the face info of all the affected vno
  //
  int count = 0;
  int vnoPlus2;
  for (vnoPlus2 = headAffectedVnosPlus2; vnoPlus2 >= 2; vnoPlus2 = affectedVnosPlus2[vnoPlus2 - 2]) {
    
    count++;
    cheapAssert(count <= mris->nvertices);
    
    int const vno = vnoPlus2 - 2;
    VERTEX_TOPOLOGY       * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;                               // ripped
    int num = 0;
    int i;
    for (i = 0; i < vt->num; i++) {
      int fno = vt->f[i];
      FACE const * const f = &mris->faces[fno];
      if (f->ripflag) continue;                             // entry to be deleted
      if (num != i) {
        vt->f[num] = fno;
        vt->n[num] = vt->n[i];                              // the face will not have had any of its vertexs removed
      }                                                     // since it is not being ripped...
      num++;
    }
    vt->num = num;
  }

  mrisCheckVertexFaceTopology(mris);

  freeAndNULL(affectedVnosPlus2);
  freeAndNULL(distCache);
  freeAndNULL(distOrigCache);
}


//======================================================================
// Functions that deal with the way the faces are attached to the vertices
// 
int computeOrientation(MRIS *mris, int f, int v0, int v1)
{
  FACE *face;

  face = &mris->faces[f];
  if (face->v[0] == v0) {
    if (face->v[1] == v1) {
      return 1;
    }
    else {
      return -1;
    }
  };
  if (face->v[1] == v0) {
    if (face->v[2] == v1) {
      return 1;
    }
    else {
      return -1;
    }
  };
  if (face->v[0] == v1) {
    return 1;
  }
  else {
    return -1;
  }
}


/*!
  Reverse order of the vertices in each face. 
  This is needed when changing the sign of the x surface coord.
*/
void MRISreverseFaceOrder(MRIS *mris)
{
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *f = &mris->faces[fno];
    int vno0 = f->v[0];
    int vno1 = f->v[1];
    int vno2 = f->v[2];
    f->v[0] = vno2;
    f->v[1] = vno1;
    f->v[2] = vno0;
    mrisSetVertexFaceIndex(mris, vno0, fno);
    mrisSetVertexFaceIndex(mris, vno1, fno);
    mrisSetVertexFaceIndex(mris, vno2, fno);
  }
  mrisCheckVertexFaceTopology(mris);
}


int MRISevertSurface(MRIS *mris)
{
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *face = &mris->faces[fno];
    if (face->ripflag) continue;
    int vno0 = face->v[0];
    int vno1 = face->v[1];
    face->v[0] = vno1;
    face->v[1] = vno0;
    mrisSetVertexFaceIndex(mris, vno0, fno);
    mrisSetVertexFaceIndex(mris, vno1, fno);
  }
  mrisCheckVertexFaceTopology(mris);

  return (NO_ERROR);
}


static short FACES_aroundVertex_reorder(MRIS *apmris, int avertex, VECTOR *pv_geometricOrder);

int MRIS_facesAtVertices_reorder(MRIS *apmris)
{
  //
  // PRECONDITIONS
  //  o <apmris> is a valid surface.
  //
  // POSTCONDITIONS
  //  o The 'f' FACE array at each vertex has its indices reordered
  //    so that bordering face indices index (i) and index (i+1)
  //    correspond to the actual geometric order of the faces about
  //    each vertex.
  //  o Note that the 'f' FACE array is changed at each vertex by
  //    this function.
  //
  cheapAssert(MRIScountAllMarked(apmris) == 0);
  
  const char * const pch_function = "MRIS_facesAtVertices_reorder";
  DebugEnterFunction((pch_function));
  fprintf(stderr, "\n");

  int ret = 1;
  
  int vno = 0;
  for (vno = 0; vno < apmris->nvertices; vno++) {
    MRIS_vertexProgress_print(apmris, vno, "Determining geometric order for vno faces...");

    VERTEX_TOPOLOGY const * const vt = &apmris->vertices_topology[vno];
    VERTEX                * const v  = &apmris->vertices         [vno];
    int const nfaces = vt->num;
    
    VECTOR *pv_geometricOrderIndx = VectorAlloc(nfaces, MATRIX_REAL);
    VECTOR *pv_logicalOrderFace   = VectorAlloc(nfaces, MATRIX_REAL);

    ret = FACES_aroundVertex_reorder(apmris, vno, pv_geometricOrderIndx);

    if (ret < 0) {

      v->marked = 1;

    } else {

      int n;
      for (n = 0; n < nfaces; n++) {
        VECTOR_ELT(pv_logicalOrderFace, n + 1) = vt->f[n];
      }
      for (n = 0; n < nfaces; n++) {
        int const orderedIndex = VECTOR_ELT(pv_geometricOrderIndx, n + 1);
        int const orderedFno   = VECTOR_ELT(pv_logicalOrderFace, orderedIndex + 1);
        
        vt->f[n] = orderedFno;
        vt->n[n] = mrisFaceVertexIndex(apmris, orderedFno, vno);
      }
    }
        
    VectorFree(&pv_geometricOrderIndx);
    VectorFree(&pv_logicalOrderFace);
  }
  
  MRISdilateMarked(apmris, 1);  // neighbors of vertices we couldn't process are also suspect and should be skipped

  xDbg_PopStack();
  
  mrisCheckVertexFaceTopology(apmris);

  return ret;
}


static void dumpFacesAroundVertex(MRIS* mris, int vno)
{
  using namespace SurfaceFromMRIS::Topology;
  fprintf(stdout, "Dumping faces around vno:%d\n",vno);
  Vertex  vertex(mris,vno);
  auto const numFaces = vertex.num(); 
  for (size_t vi = 0; vi < numFaces; vi++) {
    auto face = vertex.f(vi);
    fprintf(stdout, "  vi:%d is face:%d with vertices",vno,face.fno());
    for (size_t vi = 0; vi < VERTICES_PER_FACE; vi++) {
      auto v = face.v(vi);
      fprintf(stdout, " %d",v.vno());
    }
    fprintf(stdout, "\n");
  }
}

    
static short FACES_aroundVertex_reorder(MRIS *apmris, int avertex, VECTOR *pv_geometricOrder)
{
  //
  // PRECONDITIONS
  //  o <avertex> is a valid vertex on the surface.
  //  o <pll_faces> should not be allocated.
  //
  // POSTCONDITIONS
  //  o The face indices about vertex <avertex>, starting with face 0
  //    are returned in geometric order in the <pv_geometricOrder> vector.
  //    By geometric order is implied that the "next" index denotes
  //    the next face that directly borders the current face.
  //  o The number of connected faces is returned in the function name, or
  //    0 is there is some error.

  const char *pch_function = "FACES_aroundVertex_reorder";
  int nfaces = 0;
  int *pFaceIndex = NULL;
  int packedCount = 1;
  int i = 0;
  int I = 0;
  int j = 0;
  int k = 0;
  FACE *pFACE_I;
  FACE *pFACE_J;
  VECTOR *pv_commonVertices = NULL;  // Vector housing vertices that
  // are common between two
  // neighbouring faces.
  int commonVertices = 0;  // number of vertices in common
  // between two faces
  short b_borderFound = 0;

  pv_commonVertices = VectorAlloc(3, MATRIX_REAL);
  VERTEX_TOPOLOGY const * const pVERTEXt = &apmris->vertices_topology[avertex];
  nfaces     = pVERTEXt->num;
  pFaceIndex = pVERTEXt->f;

  DebugEnterFunction((pch_function));

  if (!nfaces) {
    fprintf(stderr,
            "\n\nFATAL ERROR encountered in function ''%s'.\nMesh structural error. Vertex %d has no faces.\n\n",
            pch_function,
            avertex);
    exit(1);
  }

  for (i = 1; i <= nfaces; i++) {
    VECTOR_ELT(pv_geometricOrder, i) = -1;
    pFACE_I = &apmris->faces[pFaceIndex[i - 1]];
  }
  VECTOR_ELT(pv_geometricOrder, 1) = 0;
  for (i = 0; i < nfaces; i++) {
    if (packedCount == nfaces) {
      break;
    }
    I = VECTOR_ELT(pv_geometricOrder, i + 1);
    pFACE_I = &apmris->faces[pFaceIndex[I]];
    for (j = 0; j < nfaces; j++) {
      k = (i + j) % nfaces;
      pFACE_J = &apmris->faces[pFaceIndex[k]];
      commonVertices = VERTICES_commonInFaces_find(pFACE_I, pFACE_J, pv_commonVertices);
      if (commonVertices == 2) {
        if (!VECTOR_elementIndex_find(pv_geometricOrder, k)) {
          VECTOR_ELT(pv_geometricOrder, i + 2) = k;
          b_borderFound = 1;
          packedCount++;
          break;
        }
      }
    }
    if (j == nfaces) DiagBreak();
  }
  VectorFree(&pv_commonVertices);
  xDbg_PopStack();

  static size_t count, limit = 1;  
  if (count++ == limit && getenv("BEVIN_TESTING_dumpFacesAroundVertex")) {
    limit *= 2;
    dumpFacesAroundVertex(apmris,avertex);
  }
  
  if (packedCount != nfaces) {
    dumpFacesAroundVertex(apmris,avertex);
    ErrorReturn(-4,
                (-4,
                 "%s: packed / faces mismatch; vertex = %d, faces = %d, packed = %d",
                 pch_function,
                 avertex,
                 nfaces,
                 packedCount));
  }
  
  return 1;
}
