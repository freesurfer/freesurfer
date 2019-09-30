// These should become invalid when XYZ are changed - NYI
//
int mrisComputeSurfaceDimensions(MRIS *mris)
{
  float xlo, ylo, zlo, xhi, yhi, zhi;

  xhi = yhi = zhi = -10000;
  xlo = ylo = zlo = 10000;

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX *vertex = &mris->vertices[vno];
    float const x = vertex->x;
    float const y = vertex->y;
    float const z = vertex->z;
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
  
  return (NO_ERROR);
}


static float mrisTriangleArea(MRIS* mris, int fno, int n)
{
  FACE const * const f = &mris->faces[fno];

  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;

  float v0[3], v1[3];

  int const vno0 = f->v[n0];
  int const vno1 = f->v[n1];
  int const vno2 = f->v[n ];
  
  v0[0] = mris->vertices[vno2].x - mris->vertices[vno0].x;
  v0[1] = mris->vertices[vno2].y - mris->vertices[vno0].y;
  v0[2] = mris->vertices[vno2].z - mris->vertices[vno0].z;
  v1[0] = mris->vertices[vno1].x - mris->vertices[vno2].x;
  v1[1] = mris->vertices[vno1].y - mris->vertices[vno2].y;
  v1[2] = mris->vertices[vno1].z - mris->vertices[vno2].z;
  
  float d1 = -v1[1] * v0[2] + v0[1] * v1[2];
  float d2 =  v1[0] * v0[2] - v0[0] * v1[2];
  float d3 = -v1[0] * v0[1] + v0[0] * v1[1];
  return sqrt(d1 * d1 + d2 * d2 + d3 * d3) / 2;
}


static void mrisFaceAreaNormal(MRIS* mris, int fno, float norm[])
{
  FACE const * const f  = &mris->faces[fno];
  int  const * const pv = f->v;
  
  int  const vno0 = pv[0];
  int  const vno1 = pv[1];
  int  const vno2 = pv[2];
  
  VERTEX const * const v0 = &mris->vertices[vno0];
  VERTEX const * const v1 = &mris->vertices[vno1];
  VERTEX const * const v2 = &mris->vertices[vno2];
  
  float const x0 = v0->x;
  float const y0 = v0->y;
  float const z0 = v0->z;
  float const x1 = v1->x;
  float const y1 = v1->y;
  float const z1 = v1->z;
  float const x2 = v2->x;
  float const y2 = v2->y;
  float const z2 = v2->z;
  
  float v20[3], v12[3];

  // Any two sides of the triangle will give the same answer
  v20[0] = x2 - x0;
  v20[1] = y2 - y0;
  v20[2] = z2 - z0;
  v12[0] = x1 - x2;
  v12[1] = y1 - y2;
  v12[2] = z1 - z2;

  // compute cross product
  norm[0] = -v12[1]*v20[2] + v20[1]*v12[2];
  norm[1] =  v12[0]*v20[2] - v20[0]*v12[2];
  norm[2] = -v12[0]*v20[1] + v20[0]*v12[1];
}


/*!
  \fn int mrisNormalFace(MRIS *mris, int fac, int n,float norm[])
  \brief Computes the normal to a triangle face. The normal will not
  have a unit length unless global variable UnitizeNormalFace=1. fac
  is the face index in mris->faces. n is the nth (0,1,2) vertex 
  
  The definition here results in the length being the sin of the angle at the vertex
  and is used to bias the sum of the face normal vectors used to compute a vertex normal vector
  
 */ 
int mrisNormalFace(MRIS* mris, int fno, int n, float norm[])
{
  FACE const * const f  = &mris->faces[fno];
  int  const * const pv = f->v;
  
  int const n0 = (n == 0) ? VERTICES_PER_FACE - 1 : n - 1;
  int const n1 = (n == VERTICES_PER_FACE - 1) ? 0 : n + 1;
  int const n2 = n;
  
  int  const vno0 = pv[n0];
  int  const vno1 = pv[n1];
  int  const vno2 = pv[n2];
  
  VERTEX const * const vtx0 = &mris->vertices[vno0];
  VERTEX const * const vtx1 = &mris->vertices[vno1];
  VERTEX const * const vtx2 = &mris->vertices[vno2];

  float const x0 = vtx0->x;
  float const y0 = vtx0->y;
  float const z0 = vtx0->z;
  float const x1 = vtx1->x;
  float const y1 = vtx1->y;
  float const z1 = vtx1->z;
  float const x2 = vtx2->x;
  float const y2 = vtx2->y;
  float const z2 = vtx2->z;

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

  return (NO_ERROR);
}


int MRIScomputeNormals(MRIS *mris)
{
  static const double RAN = 0.001; /* one thousandth of a millimeter */

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d stdout ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
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
    FACE* f = &mris->faces[k];
    
    if (f->ripflag) {
      int n;
      for (n = 0; n < VERTICES_PER_FACE; n++) {
#ifdef HAVE_OPENMP
        #pragma omp critical
#endif
        {
          mris->vertices[f->v[n]].border = TRUE;
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
          if (!mris->vertices[f->v[vi]].ripflag) break;
      }
      
      if (vi < VERTICES_PER_FACE) {
        float norm[3];

        mrisFaceAreaNormal(mris, k, norm);
      
        // The code seems to not care about the length of the face normal
        // so use the cross product of any two edges
        //
        setFaceNorm(mris, k, norm[0], norm[1], norm[2]);
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
    VERTEX* v = &mris->vertices[k];
    if (v->ripflag) continue;
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
      mris_print_hash(stdout, mris, "mris ", "\n");
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
      VERTEX                * const v  = &mris->vertices         [k];

      // calculate the vertex area (sum of the face areas)
      // and       the average of the face normals
      //
      float snorm[3];
      snorm[0] = snorm[1] = snorm[2] = 0;

      float area = 0;

      int count = 0;

      int n;
      for (n = 0; n < vt->num; n++) {
        FACE* face = &mris->faces[vt->f[n]];
        if (face->ripflag) continue;
        
        count++;
        
        float norm[3];
        mrisNormalFace(mris, vt->f[n], (int)vt->n[n], norm);
            // The normal is NOT unit length OR area length
            // Instead it's length is the sin of the angle of the vertex
            // The vertex normal is biased towards being perpendicular to 90degree contributors...

        snorm[0] += norm[0];
        snorm[1] += norm[1];
        snorm[2] += norm[2];

        area += mrisTriangleArea(mris, vt->f[n], (int)vt->n[n]);
      }
      
      if (!count || mrisNormalize(snorm) > 0.0) {       // Success?

        if (fix_vertex_area)
          v->area = area / 3.0;            // Since each face is added to three vertices...
        else
          v->area = area / 2.0;

        if (v->origarea < 0)                            // has never been set
          v->origarea = v->area;

        v->nx = snorm[0];
        v->ny = snorm[1];
        v->nz = snorm[2];
        
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
      VERTEX                * const v  = &mris->vertices         [k];

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

      v->x += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      v->y += fnv_hash(trial, k, &random_counter, -RAN, RAN);
      
      if (mris->status == MRIS_PLANE || mris->status == MRIS_CUT) {
        // nothing
      } else {
        v->z += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        
        // I don't know why this was not done for the MRIS_PLANE and _CUT
        //
        int n;
        for (n = 0; n < vt->vnum; n++) { /* if (!mris->faces[v->f[n]].ripflag) */
          VERTEX* vn = &mris->vertices[vt->v[n]];
          if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) 
            printf("   k=%5d %d nbr = %5d / %d\n", k, n, vt->v[n], vt->vnum);
          vn->x += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          vn->y += fnv_hash(trial, k, &random_counter, -RAN, RAN);
          vn->z += fnv_hash(trial, k, &random_counter, -RAN, RAN);
        }
        
        // Recompute the face norms for the affected faces
        // Note: this will recompute some, but I suspect the number is too small to be relevant
        //
        for (n = 0; n < vt->num; n++) {
          int const fno = vt->f[n];
          FACE* f = &mris->faces[fno];
          if (f->ripflag) continue;
          float norm[3];
          mrisFaceAreaNormal(mris, fno, norm);
          setFaceNorm(mris, fno, norm[0], norm[1], norm[2]);
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
    fprintf(stderr, "%s:%d MRIScomputeNormals_new could not do all vertices after %d attempts, %d remain\n",
      __FILE__, __LINE__, trial, pendingSize);
  }
#endif

  freeAndNULL(nextPending);
  freeAndNULL(pending);
  
  return (NO_ERROR);
}

#define FUNCTION_NAME MRIScomputeTriangleProperties
#include "MRIScomputeTriangleProperties_extracted.h"


/*-------------------------------------------------------------
  MRIScomputeAvgInterVertexDist() - computes the average and stddev of
  the distance between neighboring vertices. If StdDev is NULL,
  it is ignored. Requires that mrisComputeVertexDistances()
  have been run in order to compute vertex->dist[n].
  -------------------------------------------------------------*/
void MRIScomputeAvgInterVertexDist(MRIS *mris, double *StdDev)
{
  mrisCheckVertexFaceTopology(mris);
  
  bool const showHashs = false || debugNonDeterminism;

  if (showHashs) {
    fprintf(stdout, "%s:%d MRIScomputeAvgInterVertexDist starting ",__FILE__,__LINE__);
    mris_print_hash(stdout, mris, "mris ", "\n");
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
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    int const vnum = vt->vnum;
    int m;
    for (m = 0; m < vnum; m++) {
      int const vno2 = vt->v[m];
      
      VERTEX const * const v2 = &mris->vertices[vno2];
      
      if (v2->ripflag) {
        continue;
      }
      
      double d = v->dist[m];
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

  mrisSetAvgInterVertexDist(mris, Avg);
}


static int mrisOrientEllipsoid(MRI_SURFACE *mris)
{
  int fno;

  ROMP_PF_begin		// mris_fix_topology uses this
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
  for (fno = 0; fno < mris->nfaces; fno++) {
    ROMP_PFLB_begin
    
    FACE* const face = &mris->faces[fno];
    
    if (face->ripflag) {
      ROMP_PFLB_continue;
    }

    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

    /* now give the area an orientation: if the unit normal is pointing
       inwards on the ellipsoid then the area should be negative.
    */
    VERTEX const * const v0 = &mris->vertices[face->v[0]];
    VERTEX const * const v1 = &mris->vertices[face->v[1]];
    VERTEX const * const v2 = &mris->vertices[face->v[2]];

    float   const xc = (v0->x + v1->x + v2->x) /* / 3 */;   // These divides by three are a waste of time
    float   const yc = (v0->y + v1->y + v2->y) /* / 3 */;   // since we only use the magnitude of the dot product
    float   const zc = (v0->z + v1->z + v2->z) /* / 3 */;

    float   const dot = xc * fNorm->nx + yc * fNorm->ny + zc * fNorm->nz;
    
    if (dot < 0.0f) /* not in same direction, area < 0 and reverse n */
    {
      face->area *= -1.0f;

      setFaceNorm(mris, fno, -fNorm->nx, -fNorm->ny, -fNorm->nz);

      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        face->angle[ano] *= -1.0f;
      }
    }
    
    ROMP_PFLB_end
  }
  ROMP_PF_end

/* now recompute the total surface area, ignoring negative areas */
#if 0
  if ((mris->status != MRIS_PARAMETERIZED_SPHERE) || (!mris->total_area))
#endif
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

      FACE const * const face = &mris->faces[fno];
      if (face->ripflag) {
        ROMP_PF_continue;
      }
      if (face->area >= 0.0f) {
        total_area += face->area;
      }
      else {
        FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
        neg_area      += -face->area;
        neg_orig_area +=  fNorm->orig_area;
      }
      
      #undef total_area
      #undef neg_area
      #undef neg_orig_area

    #include "romp_for_end.h"

    mris->total_area    = total_area;
    mris->neg_orig_area = neg_orig_area;
    mris->neg_area      = neg_area;
  }

  return (NO_ERROR);
}



int MRISupdateEllipsoidSurface(MRI_SURFACE *mris)
{
  if (mris->status != MRIS_UNORIENTED_SPHERE) {
    mrisOrientEllipsoid(mris); /* orient the normals and angles */
  }
  return (NO_ERROR);
}



int mrisOrientPlane(MRI_SURFACE *mris)
{
  int fno;
  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *face = &mris->faces[fno];

    if (face->ripflag) continue;

    /* now give the area an orientation: if the unit normal is pointing
       downwards in the plane then the area should be negative.
    */
    FaceNormCacheEntry const * const fNorm =  getFaceNorm(mris, fno);
    if (fNorm->nz < 0.0f) {
      /* not in same direction, area < 0 and reverse n */
      face->area *= -1.0f;
      setFaceNorm(mris, fno, -fNorm->nx, -fNorm->ny, -fNorm->nz);
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        face->angle[ano] *= -1.0f;
      }
    }
  }

  /* now recompute the total surface area, ignoring negative areas */
  mris->total_area = mris->neg_orig_area = mris->neg_area = 0.0f;

  for (fno = 0; fno < mris->nfaces; fno++) {
    FACE *face = &mris->faces[fno];
    
    if (face->ripflag) continue;

    if (face->area >= 0.0f) {
      mris->total_area += face->area;
    } else {
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
      mris->neg_area += -face->area;
      mris->neg_orig_area += fNorm->orig_area;
    }
  }

  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    
    if (v->ripflag) continue;

    if (v->nz < 0) {
      v->nz *= -1;
      v->neg = 1;
    } else {
      v->neg = 0;
    }
    v->area = 0;
    for (fno = 0; fno < vt->num; fno++) {
      FACE const * face = &mris->faces[vt->f[fno]];
      v->area += face->area;
    }
    if (fix_vertex_area) {
      v->area /= 3.0;
    } else {
      v->area /= 2.0;
    }
  }

  return (NO_ERROR);
}


int mrisOrientSurface(MRI_SURFACE *mris)
{
  switch (mris->status) {
    case MRIS_RIGID_BODY:
    case MRIS_PARAMETERIZED_SPHERE:
    case MRIS_SPHERE:
    case MRIS_ELLIPSOID:
    case MRIS_SPHERICAL_PATCH:
      MRISupdateEllipsoidSurface(mris);
      break;
    case MRIS_PLANE:
      mrisOrientPlane(mris);
      break;
    default:
      /*    MRISupdateSurface(mris) ;*/
      break;
  }
  return (NO_ERROR);
}


static void MRIScomputeMetricPropertiesWkr(MRIS *mris)
{
  // fprintf(stdout,"%s:%d %s\n",__FILE__,__LINE__,__MYFUNCTION__);

  mrisCheckVertexFaceTopology(mris);

  MRIScomputeNormals(mris);             // changes XYZ

  mrisComputeSurfaceDimensions(mris);
  mrisComputeVertexDistances(mris);

  MRIScomputeTriangleProperties(mris);  // compute areas and normals
  
  mris->avg_vertex_area = mris->total_area / mris->nvertices;
  MRIScomputeAvgInterVertexDist(mris, &mris->std_vertex_dist);
  mrisOrientSurface(mris);
  // See also MRISrescaleMetricProperties()
  
  if (mris->status == MRIS_PARAMETERIZED_SPHERE || mris->status == MRIS_RIGID_BODY || mris->status == MRIS_SPHERE) {
    mris->total_area = M_PI * mris->radius * mris->radius * 4.0;
  }
}

float mrisSampleMinimizationEnergy(
    MRIS *mris, int vno, INTEGRATION_PARMS *parms, float cx, float cy, float cz)
{
  VERTEX const * const v = &mris->vertices[vno];
  
  float xw, yw, zw, dx, dy, dz, thick_sq, xp, yp, zp;

  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
  MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);

  dx = xp - xw;
  dy = yp - yw;
  dz = zp - zw;
  thick_sq = dx * dx + dy * dy + dz * dz;

  return (thick_sq);
}

float mrisSampleNormalEnergy(
    MRIS *mris, int vno, INTEGRATION_PARMS *parms, float cx, float cy, float cz)
{
  VERTEX const * const v = &mris->vertices[vno];
  
  float dx, dy, dz, len, xw, yw, zw, xp, yp, zp, pnx, pny, pnz;
  double sse;

  if (v - mris->vertices == Gdiag_no) DiagBreak();

  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
  MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);

  dx = xp - xw;
  dy = yp - yw;
  dz = zp - zw;
  len = sqrt(dx * dx + dy * dy + dz * dz);
  if (len < 0.01)  // can't reliably estimate normal. Probably not in cortex
  {
    return (0.0);
  }
  dx /= len;
  dy /= len;
  dz /= len;
  len = dx * v->wnx + dy * v->wny + dz * v->wnz;  // dot product
  len = 1 - len;
  sse = sqrt(len * len);
  //  sse = SQR(dx-v->wnx) + SQR(dy-v->wny) + SQR(dz-v->wnz) ;

  MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_NORMALS, &pnx, &pny, &pnz);
  len = sqrt(pnx * pnx + pny * pny + pnz * pnz);
  if (len < 0.01)  // can't reliably estimate normal. Probably not in cortex
    return (0.0);
  pnx /= len;
  pny /= len;
  pnz /= len;
  len = dx * pnx + dy * pny + dz * pnz;  // dot product
  len = 1 - len;
  sse += sqrt(len * len);
  //  sse += SQR(dx-pnx) + SQR(dy-pny) + SQR(dz-pnz) ;

  if (!devFinite(sse)) DiagBreak();

  return (sse);
}

float mrisSampleSpringEnergy(
    MRIS *mris, int const vno, float cx, float cy, float cz, INTEGRATION_PARMS *parms)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX          const * const v  = &mris->vertices         [vno];

  float xn, yn, zn, xp, yp, zp, xc, yc, zc;
  double sse, fdist, vdist = mris->avg_vertex_dist;
  int n, num, fno;

  FACE *face;
  MHT *mht = (MHT *)(parms->mht);

  if (vdist == 0.0) vdist = 1.0;

  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  if (v - mris->vertices == Gdiag_no) DiagBreak();

  MHTfindClosestFaceGeneric(mht, mris, cx, cy, cz, 4, 4, 1, &face, &fno, &fdist);
  if (fno < 0) MHTfindClosestFaceGeneric(mht, mris, cx, cy, cz, 1000, -1, -1, &face, &fno, &fdist);
  MRISsampleFaceCoords(mris, fno, cx, cy, cz, PIAL_VERTICES, CANONICAL_VERTICES, &xp, &yp, &zp);

  xc = yc = zc = 0;
  for (num = n = 0; n < vt->vnum; n++) {
    VERTEX const * const vn = &mris->vertices[vt->v[n]];
    if (vn->ripflag) continue;
    MRISsampleFaceCoords(mris, vn->fno, vn->x, vn->y, vn->z, PIAL_VERTICES, CANONICAL_VERTICES, &xn, &yn, &zn);
    xc += xn;
    yc += yn;
    zc += zn;
    num++;
  }
  if (num > 0) {
    xc /= num;
    yc /= num;
    zc /= num;
  }

  sse = (SQR(xc - xp) + SQR(yc - yp) + SQR(zc - zp)) / vdist;

  if (!devFinite(sse)) DiagBreak();

  return (sse);
}


float mrisSampleParallelEnergyAtVertex(MRIS *mris, int const vno, INTEGRATION_PARMS *parms)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];

  float xw, yw, zw, xp, yp, zp, dx, dy, dz, dxn, dyn, dzn, len, dxn_total, dyn_total, dzn_total;
  int n, num;
  double sse;

  MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
  MRISsampleFaceCoords(mris, v->fno, v->x, v->y, v->z, PIAL_VERTICES, CANONICAL_VERTICES, &xp, &yp, &zp);

  dx = xp - xw;
  dy = yp - yw;
  dz = zp - zw;
  len = sqrt(dx * dx + dy * dy + dz * dz);
  if (FZERO(len)) return (0.0);

  dx /= len;
  dy /= len;
  dz /= len;

  // compute average of neighboring vectors connecting white and pial surface, store it in d[xyz]n_total
  dxn_total = dyn_total = dzn_total = 0.0;
  for (num = 0, sse = 0.0, n = 0; n < vt->vnum; n++) {
    VERTEX * const vn = &mris->vertices[vt->v[n]];
    if (vn->ripflag) continue;

    MRISvertexCoord2XYZ_float(vn, WHITE_VERTICES, &xw, &yw, &zw);
    if (vn->fno >= 0) {
      MRISsampleFaceCoords(mris, vn->fno, vn->x, vn->y, vn->z, PIAL_VERTICES, CANONICAL_VERTICES, &xp, &yp, &zp);
    }
    else {
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, vn->x, vn->y, vn->z, PIAL_VERTICES, &xp, &yp, &zp);
    }

    dxn = xp - xw;
    dyn = yp - yw;
    dzn = zp - zw;
    len = sqrt(dxn * dxn + dyn * dyn + dzn * dzn);
    if (FZERO(len)) continue;

    dxn /= len;
    dyn /= len;
    dzn /= len;

    dxn_total += dxn;
    dyn_total += dyn;
    dzn_total += dzn;
#if 0
    sse += SQR(dxn-dx) + SQR(dyn-dy) + SQR(dzn-dz) ;
#else
    len = dx * dxn + dy * dyn + dz * dzn;  // dot product
    len = 1 - len;
    sse += sqrt(len * len);
#endif

    num++;
  }

#if 0
  len = sqrt(dxn_total*dxn_total + dyn_total*dyn_total + dzn_total*dzn_total) ;
  if (len > 0)
  {
    dxn_total /= len ; dyn_total /= len ; dzn_total /= len ;
    sse = SQR(dxn_total-dx) + SQR(dyn_total-dy) + SQR(dzn_total-dz) ;
    //    sse /= num ;
  }
#else
  if (num > 0) sse /= num;
#endif

  return (sse);
}

float mrisSampleParallelEnergy(
    MRIS *mris, int const vno, INTEGRATION_PARMS *parms, float cx, float cy, float cz)
{
  VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
  VERTEX                * const v  = &mris->vertices         [vno];

  int fno = 0, old_fno, num;
  double sse, fdist;
  FACE *face;
  float x, y, z;
  MHT *mht = (MHT *)(parms->mht);
  int n;

  project_point_onto_sphere(cx, cy, cz, mris->radius, &cx, &cy, &cz);
  x = v->x;
  y = v->y;
  z = v->z;          // store old coordinates
  old_fno = v->fno;  // store old face
  MHTfindClosestFaceGeneric(mht, mris, cx, cy, cz, 4, 4, 1, &face, &fno, &fdist);
  if (fno < 0) {
    MHTfindClosestFaceGeneric(mht, mris, cx, cy, cz, 1000, -1, -1, &face, &fno, &fdist);
  }
  v->fno = fno;

  v->x = cx;
  v->y = cy;
  v->z = cz;  // change coords to here and compute effects on sse
  sse = mrisSampleParallelEnergyAtVertex(mris, vno, parms);
#if 1
  for (num = 1, n = 0; n < vt->vnum; n++) {
    int const vnno = vt->v[n];
    VERTEX *vn;
    vn = &mris->vertices[vnno];
    if (vn->ripflag) continue;

    sse += mrisSampleParallelEnergyAtVertex(mris, vnno, parms);
    num++;
  }
#else
  num = 1;
#endif
  v->x = x;
  v->y = y;
  v->z = z;          // restore old coordinates
  v->fno = old_fno;  // restore old face
  sse /= (num);
  return (sse);
}

int MRISsampleFaceCoordsCanonical(
    MHT *mht, MRIS *mris, float x, float y, float z, int which, float *px, float *py, float *pz)
{
  float xv, yv, zv;
  double lambda[3], fdist, norm;
  int n, ret = NO_ERROR, fno;
  FACE *face;
  VERTEX *v;

  norm = sqrt(x * x + y * y + z * z);
  if (!FEQUAL(norm, mris->radius))  // project point onto sphere
  {
    DiagBreak();
    project_point_onto_sphere(x, y, z, mris->radius, &x, &y, &z);
  }

  xv = yv = zv = 0.0;  // to get rid of mac warnings
  MHTfindClosestFaceGeneric(mht, mris, x, y, z, 8, 8, 1, &face, &fno, &fdist);
  if (fno < 0) {
    DiagBreak();
    MHTfindClosestFaceGeneric(mht, mris, x, y, z, 1000, -1, -1, &face, &fno, &fdist);
    lambda[0] = lambda[1] = lambda[2] = 1.0 / 3.0;
  }
  else {
    ret = face_barycentric_coords(mris, fno, CANONICAL_VERTICES, x, y, z, &lambda[0], &lambda[1], &lambda[2]);
  }

  *px = *py = *pz = 0;
  for (n = 0; n < VERTICES_PER_FACE; n++) {
    v = &mris->vertices[face->v[n]];
    MRISvertexCoord2XYZ_float(v, which, &xv, &yv, &zv);
    *px += lambda[n] * xv;
    *py += lambda[n] * yv;
    *pz += lambda[n] * zv;
  }

  return (ret);
}

