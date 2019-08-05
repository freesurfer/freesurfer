#pragma once

#include "mrisurf_aaa.h"
#include "mrisurf_MRIS_MPPropertiesInVectors.h"

// Just the information needed to compute the metric properties
// Some read-only information is obtained from the underlyingMRIS


  // MRIS
  //
  // In
  //
#define ELT(C,T,N) T C N;

  #define MRIS_MP__LIST_MRIS_IN \
    ELT(const,  MRIS_Status, status         ) SEP \
    ELT(const,  MRIS_Status, origxyz_status ) SEP \
    ELT(const,  int,         nvertices      ) SEP \
    ELT(const,  int,         nfaces         ) SEP \
    ELT(const,  char,        nsize          ) SEP \
    ELT(const,  double,      radius         ) SEP \
    ELT(const,  VERTEX_TOPOLOGY const *, vertices_topology) \
    ELTX(const, FACE_TOPOLOGY   const *, faces_topology)

  // In out
  #define MRIS_MP__LIST_MRIS_IN_OUT \
    ELT(,       char,    dist_nsize ) \

  // Out
  #define MRIS_MP__LIST_MRIS_OUT \
    ELT(,       float,   xlo        ) SEP   \
    ELT(,       float,   xhi        ) SEP   \
    ELT(,       float,   ylo        ) SEP   \
    ELT(,       float,   yhi        ) SEP   \
    ELT(,       float,   zlo        ) SEP   \
    ELT(,       float,   zhi        ) SEP   \
    ELT(,       float,   xctr       ) SEP   \
    ELT(,       float,   yctr       ) SEP   \
    ELT(,       float,   zctr       ) SEP   \
    ELT(,       float,   total_area ) SEP   \
    ELT(,       double,  avg_vertex_area ) SEP   \
    ELTX(,      double,  avg_vertex_dist ) SEP   \
    ELT(,       double,  std_vertex_dist ) SEP   \
    ELT(,       float,   neg_orig_area   ) SEP   \
    ELT(,       float,   neg_area        ) 

#undef ELT

  // Vertices
  //
#define ELT(C,T,N) C T * v_##N;

  // In
  #define MRIS_MP__LIST_V_IN                \
    ELT(const,  char,   ripflag     ) SEP   \
    ELTX(const, int,    VSize       ) SEP   \
    /* following needed for SSE */          \
    ELT(const,  float,  whitex      ) SEP   \
    ELT(const,  float,  whitey      ) SEP   \
    ELT(const,  float,  whitez      ) SEP   \
    ELT(const,  float,  pialx       ) SEP   \
    ELT(const,  float,  pialy       ) SEP   \
    ELT(const,  float,  pialz       ) SEP   \
    ELT(const,  float,  wnx         ) SEP   \
    ELT(const,  float,  wny         ) SEP   \
    ELT(const,  float,  wnz         ) SEP   \
    ELTX(const,  float*, dist_orig)      /* note: these keep pointing to the original ones in the MRIS - change if code wants to change these values */
    
  // In out
  #define MRIS_MP__LIST_V_IN_OUT_XYZ        \
    ELT(,       float,  x           ) SEP   \
    ELT(,       float,  y           ) SEP   \
    ELT(,       float,  z           )

  #define MRIS_MP__LIST_V_IN_OUT_NOXYZ      \
    ELTX(,      int,    dist_capacity) SEP  \
    ELT(,       char,   border      ) SEP   \
    ELT(,       float,  cx          ) SEP   \
    ELT(,       float,  cy          ) SEP   \
    ELT(,       float,  cz          ) SEP   \
    ELT(,       float,  curv        ) SEP   \
    ELT(,       float,  origarea    ) SEP   \
    ELT(,       int,    fno         ) 

  #define MRIS_MP__LIST_V_IN_OUT            \
    MRIS_MP__LIST_V_IN_OUT_XYZ SEP MRIS_MP__LIST_V_IN_OUT_NOXYZ
  
  // Out
  #define MRIS_MP__LIST_V_OUT               \
    ELT(,       float,  area        ) SEP   \
    ELT(,       float,  nx          ) SEP   \
    ELT(,       float,  ny          ) SEP   \
    ELT(,       float,  nz          ) SEP   \
    ELTX(,      char,   neg         ) SEP   \
    ELTX(,      float*, dist        )

#undef ELT

  // Faces
  //
#define ELT(C,T,N) C T * f_##N;
  #define MRIS_MP__LIST_F_IN                                    \
    ELT(const,  char,                   ripflag         ) SEP   \
    ELTX(const, float,                  norm_orig_area  ) SEP   \
    ELTX(const, angles_per_triangle_t,  orig_angle      )
    
  #define MRIS_MP__LIST_F_OUT                                   \
    ELT(,       float,                  area            ) SEP   \
    ELTX(,      char,                   normSet         ) SEP   \
    ELTX(,      FloatXYZ,               norm            ) SEP   \
    ELTX(,      angles_per_triangle_t,  angle)

#undef ELT

#undef ELTX
#undef SEP


// Should turn these into member functions now we are using C++
//
void MRISMP_ctr(MRIS_MP* mp);
void MRISMP_dtr(MRIS_MP* mp);
void MRISMP_copy(MRIS_MP* dst, MRIS_MP* src, 
    bool only_inputs,
    bool no_need_to_copy_xyz);

void MRISMP_load(
    MRIS_MP* mp,                        // output
    MRIS*    mris,                      // input
    bool     loadOutputs = false,       // used when doing stuff based on the outputs - see MRIScomputeSSE_asThoughGradientApplied
    float *  dx_or_NULL = nullptr,      // loaded if not nullptr. The dx,dy,dz for ripped vertices are set to zero
    float *  dy_or_NULL = nullptr, 
    float *  dz_or_NULL = nullptr); 

void MRISMP_unload(MRIS* mris, MRIS_MP* mp, bool check);
