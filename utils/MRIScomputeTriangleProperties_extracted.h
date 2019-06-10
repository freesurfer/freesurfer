

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  each face has 2 triangles defined by it:

  V0       d      V3
  o--------------o
  |              |
  | A0           |
  a |              | c
  |              |
  |           A1 |
  o--------------o
  V1      b        V2

  a = V1 - V0
  d = V3 - V0
  e = V3 - V1
  A0 = 0.5 (a x d) . n

  b = V1 - V2
  c = V3 - V2
  A1 = 0.5 (c x b) . n


  each face has 1 triangle defined by it:
  V0    b     V2
  o----------o
  |         /
  | A0    /
  a |     /
  |   /
  | /
  o
  V1

  a = V1 - V0
  b = V2 - V0
  A0 = 0.5 (a x b) . n


  ------------------------------------------------------*/

#ifdef COMPILING_MRIS_MP
static void FUNCTION_NAME(MRIS_MP* mris)
#else
int FUNCTION_NAME(MRIS *mris)
#endif
{
  int const max_threads = omp_get_max_threads();
  
  VECTOR *v_a[_MAX_FS_THREADS], *v_b[_MAX_FS_THREADS], *v_n[_MAX_FS_THREADS];
  int tno;
  for (tno = 0; tno < max_threads; tno++) {
    v_a[tno] = VectorAlloc(3, MATRIX_REAL);
    v_b[tno] = VectorAlloc(3, MATRIX_REAL);
    v_n[tno] = VectorAlloc(3, MATRIX_REAL); /* normal vector */
  }

  double reduction_total_area = 0.0f;

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d %s stdout ",__FILE__,__LINE__, __MYFUNCTION__);
#ifndef COMPILING_MRIS_MP
    mris_print_hash(stdout, mris, "mris ", "");
#endif
    fprintf(stdout, "\n");
  }

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  reduction_total_area
    
  #define ROMP_FOR_LEVEL      ROMP_level_shown_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define reduction_total_area ROMP_PARTIALSUM(0)
  
    // these are done in the loop to allow parallelization
#ifndef COMPILING_MRIS_MP
    FACE *face = &mris->faces[fno];
    if (face->ripflag) continue;
    FACE const * ft = face;
#else
    if (mris->f_ripflag[fno]) continue;
    FACE_TOPOLOGY const * ft = &mris->faces_topology[fno];
#endif

    int const vno0 = ft->v[0];
    int const vno1 = ft->v[1];
    int const vno2 = ft->v[2];

#ifndef COMPILING_MRIS_MP
    VERTEX const 
      *v0 = &mris->vertices[vno0], 
      *v1 = &mris->vertices[vno1], 
      *v2 = &mris->vertices[vno2];
#else
    FloatXYZ v0xyz; v0xyz.x = mris->v_x[vno0]; v0xyz.y = mris->v_y[vno0]; v0xyz.z = mris->v_z[vno0]; 
    FloatXYZ v1xyz; v1xyz.x = mris->v_x[vno1]; v1xyz.y = mris->v_y[vno1]; v1xyz.z = mris->v_z[vno1]; 
    FloatXYZ v2xyz; v2xyz.x = mris->v_x[vno2]; v2xyz.y = mris->v_y[vno2]; v2xyz.z = mris->v_z[vno2]; 
    FloatXYZ const
     *v0 = &v0xyz,
     *v1 = &v1xyz,
     *v2 = &v2xyz;
#endif


    int const tid = omp_get_thread_num();

    VERTEX_EDGE(v_a[tid], v0, v1);
    VERTEX_EDGE(v_b[tid], v0, v2);

    /* compute metric properties of first triangle */
    {
      V3_CROSS_PRODUCT(v_a[tid], v_b[tid], v_n[tid]);

      float const area = V3_LEN(v_n[tid]) * 0.5f;
      // float const dot  = V3_DOT(v_a[tid], v_b[tid]);

#ifndef COMPILING_MRIS_MP
      face->area = area;
#else
      mris->f_area[fno] = area;
#endif
      if (area < 0) DiagBreak();

      if (!devFinite(area)) DiagBreak();

      reduction_total_area += area;
    }
    
    V3_NORMALIZE(v_n[tid], v_n[tid]); /* make it a unit vector */

    float const nx = V3_X(v_n[tid]), ny = V3_Y(v_n[tid]), nz = V3_Z(v_n[tid]);

#ifndef COMPILING_MRIS_MP
    setFaceNorm(mris, fno, nx, ny, nz);
#else
    mrismp_setFaceNorm(mris, fno, nx, ny, nz);
#endif

    /* now compute angles */
    VECTOR_LOAD(v_n[tid], nx, ny, nz);

    float dz;
    if ((V3_X(v_n[tid]) < V3_Y(v_n[tid])) && (V3_X(v_n[tid]) < V3_Z(v_n[tid]))) {
      dz = fabs(V3_X(v_n[tid]));
    }
    else if (V3_Y(v_n[tid]) < V3_Z(v_n[tid])) {
      dz = fabs(V3_Y(v_n[tid]));
    }
    else {
      dz = fabs(V3_Z(v_n[tid]));
    }
    
    int ano;
    for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
      
#ifndef COMPILING_MRIS_MP
      VERTEX const
#else
      FloatXYZ const
#endif
        *va, 
        *vb, 
        *vo;

      switch (ano) /* vertices for triangle 1 */
      {
        default:
        case 0:
          vo = v0;
          va = v2;
          vb = v1;
          break;
        case 1:
          vo = v1;
          va = v0;
          vb = v2;
          break;
        case 2:
          vo = v2;
          va = v1;
          vb = v0;
          break;
      }

      VERTEX_EDGE(v_a[tid], vo, va);
      VERTEX_EDGE(v_b[tid], vo, vb);
      float cross = VectorTripleProduct(v_b[tid], v_a[tid], v_n[tid]);
      float dot   = V3_DOT(v_a[tid], v_b[tid]);
      float angle = fastApproxAtan2f(cross, dot);
      
      
#ifndef COMPILING_MRIS_MP
      face->angle[ano] = angle;
#else
      angles_per_triangle_t* f_angle = &mris->f_angle[fno];
      (*f_angle)[ano] = angle;
#endif
    }

     #undef reduction_total_area       
   #include "romp_for_end.h"

  for (tno = 0; tno < max_threads; tno++) {
    VectorFree(&v_a[tno]);
    VectorFree(&v_b[tno]);
    VectorFree(&v_n[tno]);
  }
  
  mris->total_area = (float)reduction_total_area;

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d %s stdout ",__FILE__,__LINE__,__MYFUNCTION__);
#ifndef COMPILING_MRIS_MP
    mris_print_hash(stdout, mris, "mris ", "");
#endif
    fprintf(stdout, "\n");
  }

/* calculate the "area" of the vertices */
  int vno;
  ROMP_PF_begin		// mris_fix_topology
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(shown_reproducible)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
#ifndef COMPILING_MRIS_MP
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;
#else
    if (mris->v_ripflag[vno]) continue;
#endif

    float area = 0.0;
    int fn;
    for (fn = 0; fn < vt->num; fn++) {
      int const fno = vt->f[fn];
#ifndef COMPILING_MRIS_MP
      FACE * const face = &mris->faces[fno];
      if (face->ripflag == 0) area += face->area;
#else
      if (mris->f_ripflag[fno] == 0) area += mris->f_area[fno];
#endif
    }
    if (fix_vertex_area)
      area /= 3.0;
    else
      area /= 2.0;
#ifndef COMPILING_MRIS_MP
    v->area = area;
#else
    mris->v_area[vno] = area;
#endif
    ROMP_PFLB_end
  }
  ROMP_PF_end

#ifdef COMPILING_MRIS_MP
#else
  return (NO_ERROR);
#endif
}

#undef FUNCTION_NAME 
