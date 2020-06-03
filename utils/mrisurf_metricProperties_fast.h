// Not a header file - is included into mrisurf_metricProperties.c
//
// By having mrisurf_metricProperties_{fast,slow}.h it is easier to keep the two versions in sync


// WEIRD DISCOVERY - v_area is written twice!
//      once  by MRISMP_computeNormals 
//      later by MRISMP_computeTriangleProperties
// and I don't think it is used in between
// so the first calculation is a waste of time, except that it also writes origarea which the later one seems not to...
//
// There is a lot of overlap between these two functions that needs to be rationalized


void MRISMP_ctr(MRIS_MP* mp) {
  bzero(mp, sizeof(*mp));
}

void MRISMP_dtr(MRIS_MP* mp) {
  // Faces
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) freeAndNULL(mp->f_##N);
  MRIS_MP__LIST_F_IN SEP MRIS_MP__LIST_F_OUT
#undef ELT
#undef ELTX
#undef SEP

  // Vertices
  //
  if (mp->v_dist) {
    int vno;
    for (vno = 0; vno < mp->nvertices; vno++) {
      freeAndNULL(mp->v_dist[vno]);
    }
  }
    
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) freeAndNULL(mp->v_##N);
  MRIS_MP__LIST_V_IN SEP MRIS_MP__LIST_V_IN_OUT SEP MRIS_MP__LIST_V_OUT
#undef ELT
#undef ELTX
#undef SEP

  if (mp->v_dist_buffer) freeAndNULL(mp->v_dist_buffer);

  // MRIS
  //
  {
    // hack because mrisurf doesn't have a FACE_TOPOLOGY yet
    // When fixing, fix the dtr also!
    //
    free((void*)mp->faces_topology);
    *(FACE_TOPOLOGY**)(&mp->faces_topology) = NULL;
  }

  bzero(mp, sizeof(*mp));
}

static void MRISMP_makeDist2(MRIS_MP* mp, int vno, int neededCapacity) {
  int const vSize = mp->v_VSize[vno];

  if (neededCapacity < vSize) neededCapacity = vSize;

  if (neededCapacity <= mp->v_dist_capacity[vno]) {
    mp->v_dist[vno] = mp->v_dist_buffer[vno];
    return;
  }
  
  if (!mp->v_dist_buffer) {
    mp->v_dist_buffer = (float**)calloc(mp->nvertices,sizeof(float*));
  }
  
  if (mp->v_dist_buffer[vno]) {
    mp->v_dist[vno] = mp->v_dist_buffer[vno] = (float*)realloc(mp->v_dist_buffer[vno], neededCapacity*sizeof(float));
  } else {
    mp->v_dist[vno] = mp->v_dist_buffer[vno] = mrisStealDistStore(mp->underlyingMRIS, vno, neededCapacity);
  }
  
  mp->v_dist_capacity[vno] = neededCapacity;
}

static void MRISMP_makeDist(MRIS_MP* mp, int vno) {
  MRISMP_makeDist2(mp, vno, mp->v_dist_capacity[vno]);
}

void MRISMP_copy(MRIS_MP* dst, MRIS_MP* src, 
  bool only_inputs,
  bool no_need_to_copy_xyz) {    // NYI
  
  dst->underlyingMRIS = src->underlyingMRIS;

  if (dst->in_src != src) {
    if (dst->in_src) {
      dst->in_src->in_ref_count--;              // no longer being referenced by dst
      dst->in_src = NULL;
    }
    src->in_ref_count++;                        // stop src from destructing
    dst->in_src = src;                          // remember so can decrement later
  }

  // MRIS
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) *(T*)&(dst->N) = src->N;                   // support initializing the const members
  MRIS_MP__LIST_MRIS_IN SEP MRIS_MP__LIST_MRIS_IN_OUT
  if (!only_inputs) MRIS_MP__LIST_MRIS_OUT
#undef ELT
#undef ELTX
#undef SEP

  // Vertices
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) dst->v_##N = src->v_##N;
  MRIS_MP__LIST_V_IN 
#undef ELT
#undef ELTX
#undef SEP

#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) T* v_##N = (T*)realloc((void*)dst->v_##N, dst->nvertices*sizeof(T)); dst->v_##N = v_##N;
  MRIS_MP__LIST_V_IN_OUT SEP MRIS_MP__LIST_V_OUT
#undef ELT
#undef ELTX
#undef SEP

  if (dst->status != MRIS_PLANE) { freeAndNULL(dst->v_neg); v_neg = NULL; }

  int vno;
  for (vno = 0; vno < src->nvertices; vno++) {
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
    v_dist         [vno] = NULL;                        // will make when needed
#define ELT(C,T,N) v_##N[vno] = src->v_##N[vno];
    if (!no_need_to_copy_xyz) { MRIS_MP__LIST_V_IN_OUT_XYZ }
    MRIS_MP__LIST_V_IN_OUT_NOXYZ
    if (!only_inputs) { 
      if (v_neg) v_neg[vno] = src->v_neg[vno];
      if (src->v_dist[vno]) {
        MRISMP_makeDist2(dst, vno, src->v_dist_capacity[vno]);
        memcpy(v_dist[vno], src->v_dist[vno], src->v_dist_capacity[vno]*sizeof(*v_dist[vno]));
      }
      MRIS_MP__LIST_V_OUT 
    }
#undef ELT
#undef ELTX
#undef SEP
  }

  // Faces
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) T* f_##N = (T*)realloc((void*)dst->f_##N, dst->nfaces*sizeof(T)); dst->f_##N = f_##N;
  MRIS_MP__LIST_F_OUT
#undef ELT
#undef ELTX
#undef SEP

#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) dst->f_##N = src->f_##N;
  MRIS_MP__LIST_F_IN 
#undef ELT
#undef ELTX
#undef SEP

#ifdef MRIS_MP__LIST_F_IN_OUT
  move test inside the loop when there are MRIS_MP__LIST_F_OUT to process
#endif
  if (!only_inputs) {
    int fno;
    for (fno = 0; fno < dst->nfaces; fno++) {
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
#define ELT(C,T,N) f_##N[fno] = src->f_##N[fno];
      f_normSet[fno] = src->f_normSet[fno];
      f_norm   [fno] = src->f_norm   [fno];
      copyAnglesPerTriangle(f_angle[fno],src->f_angle[fno]); // not so special
      MRIS_MP__LIST_F_OUT 
#undef ELT
#undef ELTX
#undef SEP
    }
  }
}


void MRISMP_load(MRIS_MP* mp, MRIS* mris,
  bool loadOutputs,
  float * dx_or_NULL, float * dy_or_NULL, float * dz_or_NULL) {

  cheapAssert(!dx_or_NULL == !dy_or_NULL);
  cheapAssert(!dx_or_NULL == !dz_or_NULL);
  
  MRISMP_dtr(mp);
  MRISMP_ctr(mp);

  mp->underlyingMRIS = mris;

  // MRIS
  //
#define SEP
#define ELTX(C,T,N)
  {
    // hack because mrisurf doesn't have a FACE_TOPOLOGY yet
    // When fixing, fix the dtr also!
    //
    FACE_TOPOLOGY* ft = (FACE_TOPOLOGY*)malloc(mris->nfaces * sizeof(FACE_TOPOLOGY));
    cheapAssert(sizeof(ft->v) == sizeof(mris->faces[0].v));
    int i;
    for (i = 0; i < mris->nfaces; i++) memcpy(ft[i].v, mris->faces[i].v, sizeof(ft->v));
    *(FACE_TOPOLOGY**)(&mp->faces_topology) = ft;
  }

#define ELT(C,T,N) *(T*)&(mp->N) = mris->N;                   // support initializing the const members
  MRIS_MP__LIST_MRIS_IN SEP MRIS_MP__LIST_MRIS_IN_OUT
  if (loadOutputs) {
    ELT(,double,avg_vertex_dist)
    MRIS_MP__LIST_MRIS_OUT
  }
#undef ELT
#undef ELTX
#undef SEP

  // Vertices
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) T* v_##N = (T*)malloc(mris->nvertices*sizeof(T)); mp->v_##N = v_##N;
  MRIS_MP__LIST_V_IN SEP 
#define ELTX(C,T,N) ELT(C,T,N)
  MRIS_MP__LIST_V_IN_OUT SEP MRIS_MP__LIST_V_OUT
#undef ELT
#undef ELTX
#undef SEP

  if (mris->status != MRIS_PLANE) { freeAndNULL(mp->v_neg); v_neg = NULL; }

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const v = &mris->vertices[vno];
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
    v_dist_orig[vno]     = v->dist_orig;
    v_VSize[vno] = mrisVertexVSize(mris, vno);
    v_dist [vno] = NULL;
    v_dist_capacity[vno] = 0;
#define ELT(C,T,N) v_##N[vno] = v->N;
    MRIS_MP__LIST_V_IN SEP MRIS_MP__LIST_V_IN_OUT
    if (loadOutputs) {
      if (v_neg) v_neg[vno] = v->neg;
      MRIS_MP__LIST_V_OUT
      MRISMP_makeDist2(mp, vno, v->dist_capacity);
      memcpy(v_dist [vno], v->dist,  v_VSize[vno]*sizeof(*v_dist[vno]));
    }
#undef ELT
#undef ELTX
#undef SEP
    if (!dx_or_NULL) continue;
    if (v->ripflag) {
      dx_or_NULL[vno] = 0.0f;
      dy_or_NULL[vno] = 0.0f;
      dz_or_NULL[vno] = 0.0f;
    } else {
      dx_or_NULL[vno] = v->dx;
      dy_or_NULL[vno] = v->dy;
      dz_or_NULL[vno] = v->dz;
    }
  }
    
  // Faces
  //
#define SEP
#define ELTX(C,T,N) ELT(C,T,N)
#define ELT(C,T,N) T* f_##N = (T*)malloc(mris->nfaces*sizeof(T)); mp->f_##N = f_##N;
  MRIS_MP__LIST_F_IN SEP MRIS_MP__LIST_F_OUT
#undef ELT
#undef ELTX
#undef SEP

  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE const * const f = &mris->faces[fno];
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
    f_norm_orig_area [fno] = fNorm->orig_area;
    copyAnglesPerTriangle(f_orig_angle[fno],f->orig_angle);
#define ELT(C,T,N) f_##N[fno] = f->N;
    MRIS_MP__LIST_F_IN
    f_normSet[fno] = false;
    if (loadOutputs) {
      MRIS_MP__LIST_F_OUT
      DMATRIX* d = f->norm;
      if (!d) {
        f_norm[fno].x = 0.0;
        f_norm[fno].y = 0.0;
        f_norm[fno].z = 0.0;
      } else {
        cheapAssert(d->rows == 3);
        cheapAssert(d->cols == 1);
        f_norm[fno].x = d->rptr[0][0];
        f_norm[fno].y = d->rptr[1][0];
        f_norm[fno].z = d->rptr[2][0];
      }
      copyAnglesPerTriangle(f_angle[fno],f->angle);
    }
#undef ELT
#undef ELTX
#undef SEP
  }
}


#define comparison(NO,LHS,RHS) { comparisonWkr((NO), (LHS), (RHS), #LHS, __LINE__, &errorCount); }
static void comparisonWkr(int vnoOrFno, double lhs, double rhs, const char* expr, int line, int* errorCount) {
  if (lhs == rhs) {
    // if (vnoOrFno < 1) fprintf(stdout, "%d %s matches\n", vnoOrFno,expr);
    return;
  }
  if ((*errorCount)++ > 10) return;
  
  fprintf(stdout, "no:%d  %s  %f != %f  during count_MRIScomputeMetricProperties_calls:%d\n", 
    vnoOrFno,expr,lhs,rhs,count_MRIScomputeMetricProperties_calls);
    
  static int count;
  count++;
  if (count == 1) fprintf(stdout, "b %s:%d  needed\n",__FILE__,__LINE__);
}

void MRISMP_unload(MRIS* mris, MRIS_MP* mp, bool check) {
  int errorCount = 0;
  
  // MRIS
  //
#define SEP
#define ELTX(C,T,N)
  if (check) {
    // hack because mrisurf doesn't have a FACE_TOPOLOGY yet
    // When fixing, fix the dtr also!
    //
    FACE_TOPOLOGY const * ft = mp->faces_topology;
    int fno;
    for (fno = 0; fno < mris->nfaces; fno++) {
      comparison(fno, 0, memcmp(ft[fno].v, mris->faces[fno].v, sizeof(ft->v)));
    }
    comparison(-1, mris->avg_vertex_dist, mp->avg_vertex_dist);
  } else {
    mrisSetAvgInterVertexDist(mris, mp->avg_vertex_dist);
  }
  
#define ELT(C,T,N) comparison(-1, 1, mris->N==mp->N);
  MRIS_MP__LIST_MRIS_IN 
#undef ELT
#define ELT(C,T,N) if (check) comparison(-1,mris->N,mp->N) else mris->N = mp->N;
  MRIS_MP__LIST_MRIS_IN_OUT SEP MRIS_MP__LIST_MRIS_OUT
#undef ELT
#undef ELTX
#undef SEP

  // Vertices
  //
  errorCount = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (mp->v_ripflag[vno] && v->ripflag) continue;
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
    if (mp->v_neg) { if (check) comparison(vno,v->neg,mp->v_neg[vno]) else v->neg = mp->v_neg[vno]; }
    if (mp->v_dist[vno]) {
      if (!check) {                                                         // this requires that the dist capacity match, not just the vtotal
        mrisSetDist(mris,vno,mp->v_dist[vno],mp->v_dist_capacity[vno]);     // why assign and free when you can just move?
        mp->v_dist[vno] = mp->v_dist_buffer[vno] = NULL;
        mp->v_dist_capacity[vno] = 0;
      } else {
      
        float* const mp_dist = mp->v_dist[vno];
        float* const mris_dist = v->dist;
        
        int const vtotal = vt->vtotal;
        cheapAssert(vtotal <= v->dist_capacity);

        int n;
        for (n = 0; n < vtotal; n++) {
          if (check) comparison(vno,mris_dist[n],mp_dist[n]) else mris_dist[n] = mp_dist[n];
        }
      }
    }
#define ELT(C,T,N) comparison(vno,v->N,mp->v_##N[vno]);
    MRIS_MP__LIST_V_IN
#undef ELT
#define ELT(C,T,N) if (check) comparison(vno,v->N,mp->v_##N[vno]) else v->N = mp->v_##N[vno];
    MRIS_MP__LIST_V_IN_OUT SEP MRIS_MP__LIST_V_OUT
#undef ELT
#undef ELTX
#undef SEP
  }
    
  // Faces
  //
  errorCount = 0;
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE * const f = &mris->faces[fno];
    if (mp->f_ripflag[fno] && mris->faces[fno].ripflag) continue;
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
#define SEP
#define ELTX(C,T,N) // these are the special cases dealt with here
    if (check) comparison(fno,fNorm->orig_area,mp->f_norm_orig_area[fno]) else setFaceOrigArea(mris, fno, mp->f_norm_orig_area[fno]);
    if (check) comparison(fno,cmpAnglesPerTriangle(f->orig_angle,mp->f_orig_angle[fno]),0) else copyAnglesPerTriangle(f->orig_angle,mp->f_orig_angle[fno]);
    if (mp->f_normSet[fno]) { 
      FloatXYZ* fn = &mp->f_norm[fno]; 
      if (!check) {
        setFaceNorm(mris, fno,  fn->x,fn->y,fn->z);
      } else {
        comparison(fno,fNorm->nx,fn->x);
        comparison(fno,fNorm->ny,fn->y);
        comparison(fno,fNorm->nz,fn->z);
      }
    }
#define ELT(C,T,N) if (check) comparison(fno,f->N,mp->f_##N[fno]);
    MRIS_MP__LIST_F_IN
#undef ELT
#define ELT(C,T,N) if (check) comparison(fno,f->N,mp->f_##N[fno]) else f->N = mp->f_##N[fno];
    MRIS_MP__LIST_F_OUT
#undef ELT
#undef ELTX
#undef SEP
  }
}

#undef comparison




//
//
// These should become invalid when XYZ are changed - NYI
//
static void MRISMP_computeSurfaceDimensions(MRIS_MP* mris)
{
  float xlo, ylo, zlo, xhi, yhi, zhi;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    float const x = mris->v_x[vno];
    float const y = mris->v_y[vno];
    float const z = mris->v_z[vno];
    if (x > xhi) xhi = x;
    if (x < xlo) xlo = x;
    if (y > yhi) yhi = y;
    if (y < ylo) ylo = y;
    if (z > zhi) zhi = z;
    if (z < zlo) zlo = z;
  }
  
  mris->xlo = xlo;
  mris->xhi = xhi;
  mris->ylo = ylo;
  mris->yhi = yhi;
  mris->zlo = zlo;
  mris->zhi = zhi;
  
  mris->xctr = 0.5f * (float)((double)xlo + (double)xhi);
  mris->yctr = 0.5f * (float)((double)ylo + (double)yhi);
  mris->zctr = 0.5f * (float)((double)zlo + (double)zhi);
}


static float mrismp_TriangleArea(MRIS_MP *mris, int fno, int n)
{
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];

  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;

  float v0[3], v1[3];

  int const vno0 = ft->v[n0];
  int const vno1 = ft->v[n1];
  int const vno2 = ft->v[n ];
  
  v0[0] = mris->v_x[vno2] - mris->v_x[vno0];
  v0[1] = mris->v_y[vno2] - mris->v_y[vno0];
  v0[2] = mris->v_z[vno2] - mris->v_z[vno0];
  v1[0] = mris->v_x[vno1] - mris->v_x[vno2];
  v1[1] = mris->v_y[vno1] - mris->v_y[vno2];
  v1[2] = mris->v_z[vno1] - mris->v_z[vno2];
  
  float d1 = -v1[1] * v0[2] + v0[1] * v1[2];
  float d2 =  v1[0] * v0[2] - v0[0] * v1[2];
  float d3 = -v1[0] * v0[1] + v0[0] * v1[1];
  return sqrt(d1 * d1 + d2 * d2 + d3 * d3) / 2;
}

static void mrismp_setFaceNorm(MRIS_MP* mris, int fno, float nx, float ny, float nz) {
  mris->f_normSet[fno] = true;
  mris->f_norm   [fno].x = nx;
  mris->f_norm   [fno].y = ny;
  mris->f_norm   [fno].z = nz;
}

static void mrismp_setFaceAreaNormal(MRIS_MP* mris, int fno)
{
  // The code seems to not care about the length of the face normal
  // so use the cross product of any two edges
  //
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
  int  const * const pv = ft->v;
  
  int  const vno0 = pv[0];
  int  const vno1 = pv[1];
  int  const vno2 = pv[2];
  
  float const x0 = mris->v_x[vno0];
  float const y0 = mris->v_y[vno0];
  float const z0 = mris->v_z[vno0];
  float const x1 = mris->v_x[vno1];
  float const y1 = mris->v_y[vno1];
  float const z1 = mris->v_z[vno1];
  float const x2 = mris->v_x[vno2];
  float const y2 = mris->v_y[vno2];
  float const z2 = mris->v_z[vno2];
  
  float v20[3], v12[3];

  // Any two sides of the triangle will give the same answer
  v20[0] = x2 - x0;
  v20[1] = y2 - y0;
  v20[2] = z2 - z0;
  v12[0] = x1 - x2;
  v12[1] = y1 - y2;
  v12[2] = z1 - z2;

  // compute cross product  
  float nx = -v12[1]*v20[2] + v20[1]*v12[2];
  float ny =  v12[0]*v20[2] - v20[0]*v12[2];
  float nz = -v12[0]*v20[1] + v20[0]*v12[1];

  mrismp_setFaceNorm(mris, fno, nx, ny, nz);
}


/*!
  \fn int mrismp_NormalFace(MRIS *mris, int fac, int n,float norm[])
  \brief Computes the normal to a triangle face. The normal will not
  have a unit length unless global variable UnitizeNormalFace=1. fac
  is the face index in mris->faces. n is the nth (0,1,2) vertex 
  
  The definition here results in the length being the sin of the angle at the vertex
  and is used to bias the sum of the face normal vectors used to compute a vertex normal vector
  
 */ 
static void mrismp_NormalFace(MRIS_MP* mris, int fno, int n, float norm[])
{
  FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
  int  const * const pv = ft->v;
  
  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  int const n2 = n;
  
  int  const vno0 = pv[n0];
  int  const vno1 = pv[n1];
  int  const vno2 = pv[n2];
  
  float const x0 = mris->v_x[vno0];
  float const y0 = mris->v_y[vno0];
  float const z0 = mris->v_z[vno0];
  float const x1 = mris->v_x[vno1];
  float const y1 = mris->v_y[vno1];
  float const z1 = mris->v_z[vno1];
  float const x2 = mris->v_x[vno2];
  float const y2 = mris->v_y[vno2];
  float const z2 = mris->v_z[vno2];

  float v0[3], v1[3];
  v0[0] = x2 - x0;
  v0[1] = y2 - y0;
  v0[2] = z2 - z0;
  v1[0] = x1 - x2;
  v1[1] = y1 - y2;
  v1[2] = z1 - z2;

  mrisNormalize(v0);
  mrisNormalize(v1);

  // compute cross product
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] =  v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];

  // Note: cross product is not a unit vector even if inputs
  // are. Inputs do not need to be unit.  Until Oct 2017, this
  // function always returned a non-unitized vector. UnitizeNormalFace is a
  // global variable defined in mrisurf.h that can turn unitization on
  // and off to test its effect. Note: skull stripping uses this
  // function.
  if (UnitizeNormalFace) mrisNormalize(norm);
}


static void MRISMP_computeNormals(MRIS_MP* mris, bool check)
{
  static const double RAN = 0.001; /* one thousandth of a millimeter */

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
    // mris_print_hash(stdout, mris, "mris ", "");
    fprintf(stdout, "\n");
  }

  int k;

  // For every face, 
  // if   it is ripped then mark its vertices as .border
  // else compute its norm so we don't have to compute it later
  //
  ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
  for (k = 0; k < mris->nfaces; k++) {
    ROMP_PFLB_begin

    FACE_TOPOLOGY const * const ft = &mris->faces_topology[k];
    
    if (mris->f_ripflag[k]) {
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++) {
#ifdef HAVE_OPENMP
        #pragma omp critical
#endif
        {
          mris->v_border[ft->v[n]] = TRUE;
        }
      }
    } else {
      
      // The old code only adjusts the face norm
      // if the face is adjacent to a non-ripped vertex.
      //
      // Mimic that behavior here, for no good reason other than compatibility
      //
      int vi;
      for (vi = 0; vi < VERTICES_PER_FACE; vi++) {
          if (!mris->v_ripflag[ft->v[vi]]) break;
      }
      
      if (vi < VERTICES_PER_FACE) {
        mrismp_setFaceAreaNormal(mris, k);
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  // Build the initial pending list
  //
  int  pendingCapacity = mris->nvertices;
  int* pending         = (int*)malloc(pendingCapacity * sizeof(int));
  int* nextPending     = (int*)malloc(pendingCapacity * sizeof(int));
  int  pendingSize     = 0;

  for (k = 0; k < mris->nvertices; k++) {
    if (mris->v_ripflag[k]) continue;
    pending[pendingSize++] = k;
  }
     
  // Process the pending vertices, keeping those that need more work on the nextPending list
  // Try only a few times, because the problem might be insoluable
  //
  int trial = 0;
  for (trial = 0; (trial < 5) && (pendingSize > 0); trial++) {

    int nextPendingSize = 0;
  
    if (debugNonDeterminism) {
      fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
      // mris_print_hash(stdout, mris, "mris ", "\n");
    }

    int p;
    ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
    for (p = 0; p < pendingSize; p++) {
      ROMP_PFLB_begin

      int const     k = pending[p];
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      //RTEX                * const v  = &mris->vertices         [k];

      // calculate the vertex area (sum of the face areas)
      // and       the average of the face normals
      //
      float snorm[3];
      snorm[0] = snorm[1] = snorm[2] = 0;

      float area = 0;

      int count = 0;

      int n;
      for (n = 0; n < vt->num; n++) {
        int const fno = vt->f[n];
        //FACE* face = &mris->faces[fno];
        if (mris->f_ripflag[fno]) continue;
        
        count++;
        
        float norm[3];
        mrismp_NormalFace(mris, vt->f[n], (int)vt->n[n], norm);
            // The normal is NOT unit length OR area length
            // Instead it's length is the sin of the angle of the vertex
            // The vertex normal is biased towards being perpendicular to 90degree contributors...

        snorm[0] += norm[0];
        snorm[1] += norm[1];
        snorm[2] += norm[2];

        area += mrismp_TriangleArea(mris, vt->f[n], (int)vt->n[n]);
      }
      
      if (!count || mrisNormalize(snorm) > 0.0) {       // Success?

        if (fix_vertex_area)
          area = area / 3.0;                            // Since each face is added to three vertices...
        else
          area = area / 2.0;

        mris->v_area[k] = area;
        
        if (mris->v_origarea[k] < 0)                   // has never been set
            mris->v_origarea[k] = area;

        mris->v_nx[k] = snorm[0];
        mris->v_ny[k] = snorm[1];
        mris->v_nz[k] = snorm[2];
        
        ROMP_PFLB_continue;
      }
      
        
#ifdef HAVE_OPENMP
      #pragma omp critical                              // Retry after various adjustments below
#endif
      {
        nextPending[nextPendingSize++] = k;
      }
      
      ROMP_PFLB_end
    }
    ROMP_PF_end

    // The test I was using ALWAYS took this path!
    //
    if (nextPendingSize == 0) {
        pendingSize = 0;
        break;
    }
    
    // THIS HAS BEEN OBSERVED
    
    // Sort the nextPending list because the above appends are not in order
    //
    qsort(nextPending, nextPendingSize, sizeof(int), int_compare);
    
    // Randomly move nextPending vertices and their neighbors
    //
    // This can not be done easily in parallel because they share faces
    // If this is a performance problem, can be fixed, but I doubt it will be
    //
    for (p = 0; p < nextPendingSize; p++) {
      int     const k = nextPending[p];
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      //RTEX                * const v  = &mris->vertices         [k];

      // Warn
      //      
      if ((mris->status == MRIS_SPHERICAL_PATCH      || 
           mris->status == MRIS_PARAMETERIZED_SPHERE ||
           mris->status == MRIS_SPHERE) &&
          DIAG_VERBOSE_ON) {
        fprintf(stderr, "vertex %d: degenerate normal\n", k);
      }
    
      // randomly move x and y
      // When written, assumed being done in parallel, hence fnv_hash giving reproducible movement
      //
      int random_counter = 0;

      mris->v_x[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      mris->v_y[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      
      if (mris->status == MRIS_PLANE || mris->status == MRIS_CUT) {
        // nothing
      } else {
        mris->v_z[k] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        
        // I don't know why this was not done for the MRIS_PLANE and _CUT
        //
        int n;
        for (n = 0; n < vt->vnum; n++) { /* if (!mris->faces[v->f[n]].ripflag) */
          int const vnon = vt->v[n];
          //RTEX* vn = &mris->vertices[vnon];
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
            printf("   k=%5d %d nbr = %5d / %d\n", k, n, vt->v[n], vt->vnum);
          mris->v_x[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          mris->v_y[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          mris->v_z[vnon] += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        }
        
        // Recompute the face norms for the affected faces
        // Note: this will recompute some, but I suspect the number is too small to be relevant
        //
        for (n = 0; n < vt->num; n++) {
          int const fno = vt->f[n];
          if (mris->f_ripflag[fno]) continue;
          mrismp_setFaceAreaNormal(mris, fno);
        }
      }
    }
    
    // Try the moved ones again - although not their moved neighbors, weirdly enough
    //
    int* tempPending = pending; pending = nextPending; nextPending = tempPending;
    pendingSize = nextPendingSize;

  } // trials

#if 0  
  if (pendingSize > 0) {
    fprintf(stderr, "%s:%d MRISMP_computeNormals could not do all vertices after %d attempts, %d remain\n",
      __FILE__, __LINE__, trial, pendingSize);
  }
#endif

  freeAndNULL(nextPending);
  freeAndNULL(pending);
}


#define COMPILING_MRIS_MP


#define FUNCTION_NAME MRISMP_computeVertexDistancesWkr
#define INPUT_STATUS status
#define INPUT_X v_x
#define INPUT_Y v_y
#define INPUT_Z v_z
#define OUTPUT_DIST dist
#define OUTPUT_MAKER MRISMP_makeDist
#include "mrisComputeVertexDistancesWkr_extracted.h"

static void MRISMP_computeVertexDistances(MRIS_MP *mris) {
  MRISMP_computeVertexDistancesWkr(mris, mris->nsize, false);
  mris->dist_nsize = mris->nsize;
}

#define FUNCTION_NAME MRISMP_computeTriangleProperties
#include "MRIScomputeTriangleProperties_extracted.h"


/*-------------------------------------------------------------
  MRIScomputeAvgInterVertexDist() - computes the average and stddev of
  the distance between neighboring vertices. If StdDev is NULL,
  it is ignored. Requires that mrisComputeVertexDistances()
  have been run in order to compute vertex->dist[n].
  -------------------------------------------------------------*/
static void MRISMP_computeAvgInterVertexDist(MRIS_MP *mris, double *StdDev)
{
  bool const showHashs = false || debugNonDeterminism;

  if (showHashs) {
    fprintf(stdout, "%s:%d MRISMP_computeAvgInterVertexDist starting ",__FILE__,__LINE__);
    //mris_print_hash(stdout, mris, "mris ", "");
    fprintf(stdout, "\n");
  }
  
  double Sum = 0, Sum2 = 0;

  double N = 0.0;

  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  Sum
  #define ROMP_SUMREDUCTION1  Sum2
  #define ROMP_SUMREDUCTION2  N
    
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define Sum  ROMP_PARTIALSUM(0)
    #define Sum2 ROMP_PARTIALSUM(1)
    #define N    ROMP_PARTIALSUM(2)

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];

    if (mris->v_ripflag[vno]) {
      continue;
    }
    int const vnum = vt->vnum;
    int m;
    for (m = 0; m < vnum; m++) {
      int const vno2 = vt->v[m];
      
      if (mris->v_ripflag[vno2]) {
        continue;
      }
      float* dist = mris->v_dist[vno];
      double d = dist[m];
      Sum  += d;
      Sum2 += (d * d);
      N    += 1;
    }

    #undef Sum
    #undef Sum2
    #undef N

  #include "romp_for_end.h"
  
  // NOTE - This is a poor algorithm for computing the std dev because of how the floating point errors accumulate
  // but it seems to work for us because the double has enough accuracy to sum the few hundred thousand small but not
  // too small floats that we have
  double Avg = Sum / N;
  if (StdDev != NULL) {
    *StdDev = sqrt(N * (Sum2 / N - Avg * Avg) / (N - 1));
  }

  // printf("\n\nN = %ld, Sum = %g, Sum2 = %g, Avg=%g, Std = %g\n\n",
  // N,Sum,Sum2,Avg,*StdDev);

  mris->avg_vertex_dist = Avg;
}


static void mrismp_OrientEllipsoid(MRIS_MP *mris)
{
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f;

  double total_area = 0.0, neg_area = 0.0, neg_orig_area = 0.0;

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces

  #define ROMP_SUMREDUCTION0  total_area
  #define ROMP_SUMREDUCTION1  neg_area
  #define ROMP_SUMREDUCTION2  neg_orig_area

  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible

#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define total_area      ROMP_PARTIALSUM(0)
    #define neg_area        ROMP_PARTIALSUM(1)
    #define neg_orig_area   ROMP_PARTIALSUM(2)
    
    FACE_TOPOLOGY const * const ft = &mris->faces_topology[fno];
    
    if (mris->f_ripflag[fno]) {
      ROMP_PFLB_continue;
    }

    FloatXYZ const * const fNorm = &mris->f_norm[fno];
    float const nx = fNorm->x;
    float const ny = fNorm->y;
    float const nz = fNorm->z;

    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
    */
    int const vno0 = ft->v[0];
    int const vno1 = ft->v[1];
    int const vno2 = ft->v[2];
    
    float x0 = mris->v_x[vno0];
    float y0 = mris->v_y[vno0];
    float z0 = mris->v_z[vno0];
    float x1 = mris->v_x[vno1];
    float y1 = mris->v_y[vno1];
    float z1 = mris->v_z[vno1];
    float x2 = mris->v_x[vno2];
    float y2 = mris->v_y[vno2];
    float z2 = mris->v_z[vno2];

    float   const xc = (x0 + x1 + x2) /* / 3 */;   // These divides by three are a waste of time
    float   const yc = (y0 + y1 + y2) /* / 3 */;   // since we only use the magnitude of the dot product
    float   const zc = (z0 + z1 + z2) /* / 3 */;

    float   const dot = xc * nx + yc * ny + zc * nz;
    
    if (dot < 0.0f) /* not in same direction, area < 0 and reverse n */
    {
      mris->f_area[fno]   *= -1.0f;

      mris->f_norm[fno].x *= -1.0f;
      mris->f_norm[fno].y *= -1.0f;
      mris->f_norm[fno].z *= -1.0f;

      angles_per_triangle_t* angle = &mris->f_angle[fno];
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        (*angle)[ano]      *= -1.0f;
      }
    }
    
    if (mris->f_area[fno] >= 0.0f) {
      total_area += mris->f_area[fno];
    } else {
      neg_area      += -mris->f_area[fno];
      neg_orig_area +=  mris->f_norm_orig_area[fno];
    }

    #undef total_area
    #undef neg_area
    #undef neg_orig_area

  #include "romp_for_end.h"

  mris->total_area    = total_area;
  mris->neg_area      = neg_area;
  mris->neg_orig_area = neg_orig_area;
}



static void MRISMP_updateEllipsoidSurface(MRIS_MP* mris)
{
  if (mris->status != MRIS_UNORIENTED_SPHERE) {
    mrismp_OrientEllipsoid(mris); /* orient the normals and angles */
  }
}



static void mrismp_OrientPlane(MRIS_MP* mris)
{
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f;

  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {

    if (mris->f_ripflag[fno]) continue;

    // give the area an orientation: if the unit normal is pointing
    // downwards in the plane then the area should be negative.
    //
    FloatXYZ* const fNorm = &mris->f_norm[fno];
    float const nx = fNorm->x;
    float const ny = fNorm->y;
    float const nz = fNorm->z;
    
    if (nz < 0.0f) {
      /* not in same direction, area < 0 and reverse n */
      mris->f_area[fno]   *= -1.0f;
      
      fNorm->x = -nx;
      fNorm->y = -ny;
      fNorm->z = -nz;

      angles_per_triangle_t* angle = &mris->f_angle[fno];
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        (*angle)[ano]      *= -1.0f;
      }
    }
    
    if (mris->f_area[fno] >= 0.0f) {
      mris->total_area += mris->f_area[fno];
    } else {
      mris->neg_area      += -mris->f_area[fno];
      mris->neg_orig_area +=  mris->f_norm_orig_area[fno];
    }
    
  }

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    
    if (mris->v_ripflag[vno]) continue;

    if (mris->v_nz[vno] < 0) {
      mris->v_nz [vno] *= -1.0f;
      mris->v_neg[vno] = 1;
    } else {
      mris->v_neg[vno] = 0;
    }
    float v_area = 0.0f;
    int fn;
    for (fn = 0; fn < vt->num; fn++) {
      int fno = vt->f[fn];
      v_area += mris->f_area[fno];
    }
    if (fix_vertex_area) {
      v_area /= 3.0;
    } else {
      v_area /= 2.0;
    }
    mris->v_area[vno] = v_area;
  }
}


int mrismp_OrientSurface(MRIS_MP *mris)
{
  switch (mris->status) {
    case MRIS_RIGID_BODY:
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
    case MRIS_ELLIPSOID:
    case MRIS_SPHERICAL_PATCH:
      MRISMP_updateEllipsoidSurface(mris);
      break;
    case MRIS_PLANE:
      mrismp_OrientPlane(mris);
      break;
    default:
      /*    MRISupdateSurface(mris) ;*/
      break;
  }
  return (NO_ERROR);
}


void MRIScomputeMetricProperties(MRIS_MP* mris) {
  // fprintf(stdout,"%s:%d %s\n",__FILE__,__LINE__,__MYFUNCTION__);

  MRISMP_computeNormals(mris, false);           // changes XYZ
  MRISMP_computeSurfaceDimensions(mris);
  MRISMP_computeVertexDistances(mris);
  MRISMP_computeTriangleProperties(mris);       // compute face areas and normals
  
  mris->avg_vertex_area = mris->total_area / mris->nvertices;

  MRISMP_computeAvgInterVertexDist(mris, &mris->std_vertex_dist);

  mrismp_OrientSurface(mris);
  // See also MRISrescaleMetricProperties()
  
  if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_RIGID_BODY || mris->status == MRIS_SPHERE) {
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0;
  }
}

#undef COMPILING_MRIS_MP
