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
#include "mrisurf_sseTerms.h"
#include "mrisurf_project.h"

#include "mrisurf_base.h"
#include "mrisurf_MRIS_MP.h"
#include "mrisurf_SurfaceFromMRIS_MP_generated.h"

#include "mrishash_SurfaceFromMRIS.h"


#define MAX_VOXELS          mrisurf_sse_MAX_VOXELS
#define MAX_DISPLACEMENT    mrisurf_sse_MAX_DISPLACEMENT 
#define DISPLACEMENT_DELTA  mrisurf_sse_DISPLACEMENT_DELTA
#define DEFAULT_STD         mrisurf_sse_DEFAULT_STD


// The SSE terms are either computed by iterating over the vertices or the faces
// Ideally there would be one pass over each, to 
//      a) minimize the scanning logic
//      b) minimize the cache traffic
// Furthermore it would be more efficient to do these using the dense MRIS_MP format rather than the MRIS format
// We are working towards this ideal...
//
// For now, only the terms used in recon-all bert have been moved here
// The rest will follow


static bool const trace = false;

struct SseTermsBase {
#ifndef METRIC_SCALE
    static double const area_scale = 1.0;
    static double const dist_scale = 1.0;
#endif
};


// This template is for all the terms that can be computed from any Surface 
//
template <class _Surface>
struct SseTerms_DistortedSurfaces : public SseTermsBase {
    typedef          _Surface        Surface;
    typedef typename Surface::Face   Face;
    typedef typename Surface::Vertex Vertex;
    
    Surface surface;
#if METRIC_SCALE
    double const area_scale;
    double const dist_scale;
#endif

    int const vnoBegin;
    int const vnoEnd;
    int const fnoBegin;
    int const fnoEnd;

    SseTerms_DistortedSurfaces(Surface surface, int selector) 
      : surface(surface), 
#if METRIC_SCALE                                                         
        area_scale((surface.patch())                               ? 1.0 :  // it is weird that this doesn't have all the alternatives the dist has
                                                                         (surface.orig_area() / (surface.total_area()                     ))
                  ),
        dist_scale((surface.patch())                               ? 1.0 : 
                   (surface.status() == MRIS_PARAMETERIZED_SPHERE) ? sqrt(surface.orig_area() / (surface.total_area()                     )) :
                   (surface.neg_area() < surface.total_area())     ? sqrt(surface.orig_area() / (surface.total_area() - surface.neg_area())) :
                                                                     sqrt(surface.orig_area() / (surface.total_area()                     ))
                  ),
#endif
        vnoBegin(selector >= 0 ? selector   : 0),
        vnoEnd  (selector >= 0 ? selector+1 : surface.nvertices()),
        fnoBegin(selector >= 0 ? selector   : 0),
        fnoEnd  (selector >= 0 ? selector+1 : surface.nfaces())
    {}


    //========================================
    // Energy terms that iterate over the vertices, but not their neighbours
    //

    //========================================
    // Energy terms that iterate over the vertices, and for each iterate over all their immediate (1-hop) neighbours
    //
    double RepulsiveRatioEnergy( double l_repulse)
    {
        if (FZERO(l_repulse))
            return (0.0);

        double sse_repulse = 0.0;
        for (int vno = vnoBegin; vno < vnoEnd; vno++) {
            auto const v = surface.vertices(vno);
            if (v.ripflag()) continue;

            double const x  = v.x(),  y  = v.y(),  z  = v.z();
            double const cx = v.cx(), cy = v.cy(), cz = v.cz();

            double v_sse = 0.0;
            for (int n = 0; n < v.vnum(); n++) {
                auto const vn = v.v(n);
                if (vn.ripflag()) continue;

                double const  dx =  x - vn.x(),  dy  =  y - vn.y(),  dz  = z  - vn.z();
                double const cdx = cx - vn.cx(), cdy = cy - vn.cy(), cdz = cz - vn.cz();

                double const dist       = sqrt(dx*dx   + dy*dy   + dz*dz);
                double const canon_dist = sqrt(cdx*cdx + cdy*cdy + cdz*cdz) + REPULSE_E;

                double const adjusted_dist = dist/canon_dist + REPULSE_E;
                v_sse += REPULSE_K / (adjusted_dist * adjusted_dist);
            }
            sse_repulse += v_sse;
        }

        return l_repulse*sse_repulse;
    }
    
    double SpringEnergy()
    {
        double sse_spring = 0.0;
        for (int vno = vnoBegin; vno < vnoEnd; vno++) {
            auto const v = surface.vertices(vno);
            if (v.ripflag()) continue;
            double v_sse = 0.0;
            for (int n = 0; n < v.vnum(); n++) {
                v_sse += square(v.dist(n));
            }
            sse_spring += area_scale * v_sse;
        }
        return sse_spring;
    }

    double LaplacianEnergy()
    {
        //  Note: this function assumes that the mris surface has the original
        //  (i.e. after global rotational alignment)
        //  spherical coordinates in the TMP2_VERTICES
        //
        double sse_lap = 0.0;
        for (int vno = vnoBegin; vno < vnoEnd; vno++) {
            auto const v = surface.vertices(vno);
            if (v.ripflag()) continue;

            double const vx = v.x() - v.t2x(), vy = v.y() - v.t2y(), vz = v.z() - v.t2z();
            double v_sse = 0.0; 
            for (int n = 0; n < v.vnum(); n++) {
                auto const vn = v.v(n);
                double const vnx = vn.x() - vn.t2x(), vny = vn.y() - vn.t2y(), vnz = vn.z() - vn.t2z();
                double const dx  = vnx - vx, dy = vny - vy, dz = vnz - vz;
                double const error = square(dx) + square(dy) + square(dz);
                v_sse += error;
            }
            sse_lap += area_scale * v_sse;
        }

        return sse_lap;
    }

    double TangentialSpringEnergy()
    {
        double sse_spring = 0.0;
        for (int vno = vnoBegin; vno < vnoEnd; vno++) {
            auto const v = surface.vertices(vno);
            if (v.ripflag()) continue;

            float const x = v.x(), v_nx = v.nx();
            float const y = v.y(), v_ny = v.ny();
            float const z = v.z(), v_nz = v.nz();

            double v_sse = 0.0;
            for (int n = 0; n < v.vnum(); n++) {
                auto const vn = v.v(n);
                
                float dx = vn.x() - x;
                float dy = vn.y() - y;
                float dz = vn.z() - z;
                
                float const nc = dx * v_nx + dy * v_ny + dz * v_nz;
                dx -= nc * v_nx;
                dy -= nc * v_ny;
                dz -= nc * v_nz;
                
                float const dist_sq = square(dx) + square(dy) + square(dz);
                
                v_sse += dist_sq;
            }
            sse_spring += area_scale * v_sse;
        }
        return sse_spring;
    }

    //========================================
    // Error terms that iterate over the vertices, but not their neighbours
    //

    //========================================
    // Error terms that iterate over the vertices, and for each iterate over all their neighbours
    //
    
    //========================================
    // misc
    //
    double NonlinearDistanceSSE()
    {
        double sse_dist = 0.0;

        for (int vno = vnoBegin; vno < vnoEnd; vno++) {
            auto const v = surface.vertices(vno);
            if (v.ripflag()) continue;

            double v_sse = 0.0;
            for (int n = 0; n < v.vtotal(); n++) {                              // why is this vtotal and others are vnum?
              if (FZERO(v.dist_orig(n))) continue;

              double const ratio = dist_scale * v.dist(n) / v.dist_orig(n);
              double const delta = log(1 + exp(ratio));
              v_sse += delta;
            }

            sse_dist += v_sse;
        }

        return sse_dist;
    }


    double NonlinearAreaSSE()
    {
        double sse = 0.0;
  
        #define ROMP_VARIABLE       fno
        #define ROMP_LO             fnoBegin
        #define ROMP_HI             fnoEnd

        #define ROMP_SUMREDUCTION0  sse

        #define ROMP_FOR_LEVEL      ROMP_level_fast     // it is ROMP_level_assume_reproducible but that might change results

#ifdef ROMP_SUPPORT_ENABLED
        const int romp_for_line = __LINE__;
#endif
        #include "romp_for_begin.h"
        ROMP_for_begin

            #define sse  ROMP_PARTIALSUM(0)

            auto const face = surface.faces(fno);
            if (face.ripflag()) ROMP_PF_continue;

            double const ratio = min(max(-MAX_NEG_RATIO, area_scale * face.area()), MAX_NEG_RATIO);
            double const error = (log(1.0 + exp(NEG_AREA_K * ratio)) / NEG_AREA_K) - ratio;

            sse += error;

            #undef sse
        #include "romp_for_end.h"

        return sse;
    }

};


// Choose the Surface et al to be used for MRIS and for MRIS_MP
//
typedef SurfaceFromMRIS::Distort::Surface    SSE_Surface_type_MRIS;
typedef SurfaceFromMRIS_MP::Distort::Surface SSE_Surface_types_MRIS_MP;



// This maps the MRIS and MRIS_MP to the set of functions that are implemented for them
// whether done using the SseTerms_DistortedSurfaces implementations or specific implementations or not at all
//
typedef SseTerms_DistortedSurfaces<SSE_Surface_type_MRIS> SseTerms_Template_for_SurfaceFromMRIS;

struct SseTerms_MRIS : public SseTerms_Template_for_SurfaceFromMRIS {
    MRIS* const mris;
    SseTerms_MRIS(MRIS* const mris, int selector) : SseTerms_Template_for_SurfaceFromMRIS(Surface(mris),selector), mris(mris) {}
    
    #define MRIS_PARAMETER          
    #define MRIS_PARAMETER_COMMA
    #define NOCOMMA_SELECTOR
    #define COMMA_SELECTOR
    #define SEP 
    #define ELT(NAME, SIGNATURE, CALL)    double NAME SIGNATURE;
    LIST_OF_SSETERMS
    #undef ELT
    #undef SEP
    #undef NOCOMMA_SELECTOR
    #undef COMMA_SELECTOR
    #undef MRIS_PARAMETER_COMMA
    #undef MRIS_PARAMETER
};


typedef SseTerms_DistortedSurfaces<SSE_Surface_types_MRIS_MP> SseTerms_Template_for_SurfaceFromMRIS_MP;

struct SseTerms_MRIS_MP : public SseTerms_Template_for_SurfaceFromMRIS_MP {
    MRIS_MP* const mris;
    SseTerms_MRIS_MP(MRIS_MP* const mris, int selector) : SseTerms_Template_for_SurfaceFromMRIS_MP(Surface(mris),selector), mris(mris) {}
    
    #define MRIS_PARAMETER          
    #define MRIS_PARAMETER_COMMA
    #define NOCOMMA_SELECTOR
    #define COMMA_SELECTOR
    #define SEP 
    // #define ELT(NAME, SIGNATURE, CALL) double NAME SIGNATURE;
    // andrew temporarily adding these definitions since they don't get defined anywhere else
    #define ELT(NAME, SIGNATURE, CALL) inline double NAME SIGNATURE { fs::fatal() << #NAME << "() has not been implemented for SseTerms_MRIS_MP"; return 0; };
    LIST_OF_SSETERMS
    #undef ELT
    #undef SEP
    #undef NOCOMMA_SELECTOR
    #undef COMMA_SELECTOR
    #undef MRIS_PARAMETER_COMMA
    #undef MRIS_PARAMETER
};



//=============
// Energy Terms
//
double SseTerms_MRIS::RepulsiveRatioEnergy(double l_repulse)
{
    return SseTerms_Template_for_SurfaceFromMRIS::RepulsiveRatioEnergy(l_repulse);
}


double SseTerms_MRIS::SpringEnergy()
{
    return SseTerms_Template_for_SurfaceFromMRIS::SpringEnergy();
}

double SseTerms_MRIS::LaplacianEnergy()
{
    return SseTerms_Template_for_SurfaceFromMRIS::LaplacianEnergy();
}

double SseTerms_MRIS::TangentialSpringEnergy()
{
    return SseTerms_Template_for_SurfaceFromMRIS::TangentialSpringEnergy();
}

// Error terms
//


// Misc
//
double SseTerms_MRIS::NonlinearDistanceSSE()
{
    return SseTerms_Template_for_SurfaceFromMRIS::TangentialSpringEnergy();
}

double SseTerms_MRIS::NonlinearAreaSSE()
{
    return SseTerms_Template_for_SurfaceFromMRIS::NonlinearAreaSSE();
}




//========================================
// Ones that have not been tidied up yet


/*-----------------------------------------------------
  Fit a 2-d quadratic to the surface locally and compute the SSE as
  the square of the constant term (the distance the quadratic fit surface
  is from going through the central vertex)
  ------------------------------------------------------*/
double SseTerms_MRIS::QuadraticCurvatureSSE( double l_curv)            // BEVIN mris_make_surfaces 3
{
  if (FZERO(l_curv)) {
    return (NO_ERROR);
  }

  mrisComputeTangentPlanes(mris);
  
  typedef struct Reused {
    VECTOR * v_n  ;
    VECTOR * v_P  ;
    VECTOR * v_e1 ;
    VECTOR * v_e2 ;
    VECTOR * v_nbr;
  } Reused;
  
#ifdef HAVE_OPENMP
  int const maxThreads = omp_get_max_threads();
#else
  int const maxThreads = 1;
#endif
  Reused* reusedByThread = (Reused*)calloc(maxThreads, sizeof(Reused));
  
  double sse = 0.0;

#ifdef BEVIN_MRISCOMPUTEQUADRATICCURVATURESSE_REPRODUCIBLE

  #define ROMP_VARIABLE       vno 
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse ROMP_PARTIALSUM(0)

#else
  
  int vno;
  
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+:sse)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif

#ifdef HAVE_OPENMP
    int const tid = omp_get_thread_num();
#else
    int const tid = 0;
#endif

    Reused* reused = reusedByThread + tid;
    #define REUSE(NAME,DIM) \
        VECTOR* NAME = reused->NAME; if (!NAME) NAME = reused->NAME = VectorAlloc(DIM, MATRIX_REAL);
    REUSE(v_n   ,3)
    REUSE(v_P   ,5)
    REUSE(v_e1  ,3)
    REUSE(v_e2  ,3)
    REUSE(v_nbr ,3)
    #undef REUSE

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    VECTOR* v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL);    /* heights above TpS */
    VECTOR_LOAD(v_n,  v->nx,  v->ny,  v->nz);
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z);
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z);

    MATRIX* m_X = MatrixAlloc(vt->vtotal, 5, MATRIX_REAL); /* 2-d quadratic fit */
    int n;
    for (n = 0; n < vt->vtotal; n++) /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      VERTEX_EDGE(v_nbr, v, vn);
      VECTOR_ELT(v_Y, n + 1) = V3_DOT(v_nbr, v_n);
      float ui = V3_DOT(v_e1, v_nbr);
      float vi = V3_DOT(v_e2, v_nbr);

      *MATRIX_RELT(m_X, n + 1, 1) = ui * ui;
      *MATRIX_RELT(m_X, n + 1, 2) = vi * vi;
      *MATRIX_RELT(m_X, n + 1, 3) = ui;
      *MATRIX_RELT(m_X, n + 1, 4) = vi;
      *MATRIX_RELT(m_X, n + 1, 5) = 1;
    }

    MATRIX* m_X_inv = MatrixPseudoInverse(m_X, NULL);
    if (!m_X_inv) {
      MatrixFree(&m_X);
      VectorFree(&v_Y);
      continue;
    }
    
    v_P = MatrixMultiply(m_X_inv, v_Y, v_P);
    //oat a = VECTOR_ELT(v_P, 1);
    float e = VECTOR_ELT(v_P, 5);

    sse += e * e;
    if (vno == Gdiag_no) printf("v %d: e=%2.2f, curvature sse %2.2f\n", vno, e, e * e);

    MatrixFree(&m_X);
    MatrixFree(&m_X_inv);
    VectorFree(&v_Y);
    
#ifdef BEVIN_MRISCOMPUTEQUADRATICCURVATURESSE_REPRODUCIBLE
    #undef sse
  #include "romp_for_end.h"
#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

  { int tid;
    for (tid = 0; tid < maxThreads; tid++) {
      Reused* reused = reusedByThread + tid;
      if (reused->v_n)   VectorFree(&reused->v_n);
      if (reused->v_e1)  VectorFree(&reused->v_e1);
      if (reused->v_e2)  VectorFree(&reused->v_e2);
      if (reused->v_nbr) VectorFree(&reused->v_nbr);
      if (reused->v_P)   VectorFree(&reused->v_P);
  } }
  
  return (NO_ERROR);
}




//====================================================================================
// VERTEX SSE TERMS
//
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double SseTerms_MRIS::ThicknessMinimizationEnergy( double l_thick_min, INTEGRATION_PARMS *parms)
{
  int vno;
  double sse_tmin;
  static int cno = 0;
  static double last_sse[MAXVERTICES];

  if (FZERO(l_thick_min)) {
    return (0.0);
  }

  if (cno == 0) {
    memset(last_sse, 0, sizeof(last_sse));
  }
  cno++;

  sse_tmin = 0.0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse_tmin)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    float thick_sq;
    VERTEX *v;
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    thick_sq = mrisSampleMinimizationEnergy(mris, vno, parms, v->x, v->y, v->z);

    if (vno < MAXVERTICES && thick_sq > last_sse[vno] && cno > 1 && vno == Gdiag_no) DiagBreak();

    if (vno < MAXVERTICES && (thick_sq > last_sse[vno] && cno > 1)) DiagBreak();

    if (vno < MAXVERTICES) last_sse[vno] = thick_sq;
    // diagnostics end

    v->curv = sqrt(thick_sq);
    sse_tmin += thick_sq;
    if (Gdiag_no == vno) {
      printf("E_thick_min:  v %d @ (%2.2f, %2.2f, %2.2f): thick = %2.5f\n", vno, v->x, v->y, v->z, v->curv);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  sse_tmin /= 2;
  return (sse_tmin);
}


double SseTerms_MRIS::ThicknessParallelEnergy( double l_thick_parallel, INTEGRATION_PARMS *parms)
{
  int vno, max_vno;
  double sse_tparallel, max_inc;
  static int cno = 0;
  static double last_sse[MAXVERTICES];

  if (FZERO(l_thick_parallel)) {
    return (0.0);
  }
  if (cno == 0) {
    memset(last_sse, 0, sizeof(last_sse));
  }
  cno++;

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time

  max_inc = 0;
  max_vno = 0;
  sse_tparallel = 0.0;
  ROMP_PF_begin
  // ifdef HAVE_OPENMP
  // pragma omp parallel for if_ROMP(experimental) reduction(+:sse_tparallel)
  // endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    
    double sse;

    VERTEX * const v = &mris->vertices[vno];
    if (v->ripflag) continue;
    if (vno == Gdiag_no) DiagBreak();

    sse = mrisSampleParallelEnergy(mris, vno, parms, v->x, v->y, v->z);
    if ((vno < MAXVERTICES) && (sse > last_sse[vno] && cno > 1 && vno == Gdiag_no)) DiagBreak();

    if ((vno < MAXVERTICES) && (sse > last_sse[vno] && cno > 1)) {
      if (sse - last_sse[vno] > max_inc) {
        max_inc = sse - last_sse[vno];
        max_vno = vno;
      }
      DiagBreak();
    }

    if (vno < MAXVERTICES) last_sse[vno] = sse;
    sse_tparallel += sse;
    if (vno == Gdiag_no) {
      printf("E_parallel: vno = %d, E = %f\n", vno, sse);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  sse_tparallel /= 2;
  return (sse_tparallel);
}


double SseTerms_MRIS::ThicknessSmoothnessEnergy( double l_tsmooth, INTEGRATION_PARMS *parms)
{
  int vno, n;
  double sse_tsmooth, v_sse, dn, dx, dy, dz, d0;
  float xp, yp, zp;

  if (FZERO(l_tsmooth)) {
    return (0.0);
  }

  for (sse_tsmooth = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);

    d0 = SQR(xp - v->whitex) + SQR(yp - v->whitey) + SQR(zp - v->whitez);
    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (!vn->ripflag) {
        MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, vn->x, vn->y, vn->z, PIAL_VERTICES, &xp, &yp, &zp);

        dx = xp - vn->whitex;
        dy = yp - vn->whitey;
        dz = zp - vn->whitez;
        dn = (dx * dx + dy * dy + dz * dz);
        v_sse += (dn - d0) * (dn - d0);
      }
    }
    sse_tsmooth += v_sse;
  }
  return (sse_tsmooth);
}


static double big_sse = 10.0;
double SseTerms_MRIS::ThicknessNormalEnergy( double l_thick_normal, INTEGRATION_PARMS *parms)
{
  int vno;
  double sse_tnormal;
  static int cno = 0;
  static double last_sse[MAXVERTICES];

  if (FZERO(l_thick_normal)) return (0.0);

  if (cno == 0) memset(last_sse, 0, sizeof(last_sse));

  cno++;

  sse_tnormal = 0.0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse_tnormal)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    double sse;
    VERTEX *v;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    sse = mrisSampleNormalEnergy(mris, vno, parms, v->x, v->y, v->z);
    if (sse > big_sse) DiagBreak();

    if (vno < MAXVERTICES && ((sse > last_sse[vno] && cno > 1 && vno == Gdiag_no) || (sse > last_sse[vno] && cno > 1)))
      DiagBreak();

    sse_tnormal += sse;
    if (vno < MAXVERTICES) last_sse[vno] = sse;
    if (Gdiag_no == vno) {
      float E;
      float dx, dy, dz, len, xw, yw, zw, xp, yp, zp, cx, cy, cz;

      cx = v->x;
      cy = v->y;
      cz = v->z;
      E = mrisSampleNormalEnergy(mris, vno, parms, v->x, v->y, v->z);
      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);
      dx = xp - xw;
      dy = yp - yw;
      dz = zp - zw;
      len = sqrt(dx * dx + dy * dy + dz * dz);
      if (FZERO(len) == 0) {
        dx /= len;
        dy /= len;
        dz /= len;
      }
      printf("E_thick_normal: vno %d, E=%f, N = (%2.2f, %2.2f, %2.2f), D = (%2.2f, %2.2f, %2.2f), dot= %f\n",
             vno,
             E,
             v->wnx,
             v->wny,
             v->wnz,
             dx,
             dy,
             dz,
             v->wnx * dx + v->wny * dy + v->wnz * dz);
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  sse_tnormal /= 2;
  return (sse_tnormal);
}


double SseTerms_MRIS::ThicknessSpringEnergy( double l_thick_spring, INTEGRATION_PARMS *parms)
{
  int vno;
  double sse_spring, sse;
  VERTEX *v;

  if (FZERO(l_thick_spring)) {
    return (0.0);
  }

  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) continue;

    sse = mrisSampleSpringEnergy(mris, vno, v->x, v->y, v->z, parms);

    sse_spring += sse;
    if (Gdiag_no == vno) {
      float E;
      float dx, dy, dz, len, xw, yw, zw, xp, yp, zp, cx, cy, cz;

      cx = v->x;
      cy = v->y;
      cz = v->z;
      E = mrisSampleSpringEnergy(mris, vno, v->x, v->y, v->z, parms);
      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, cx, cy, cz, PIAL_VERTICES, &xp, &yp, &zp);
      dx = xp - xw;
      dy = yp - yw;
      dz = zp - zw;
      len = sqrt(dx * dx + dy * dy + dz * dz);
      if (FZERO(len) == 0) {
        dx /= len;
        dy /= len;
        dz /= len;
      }
      printf("E_thick_spring: vno %d, E=%f, N = (%2.2f, %2.2f, %2.2f), D = (%2.2f, %2.2f, %2.2f), dot= %f\n",
             vno,
             E,
             v->wnx,
             v->wny,
             v->wnz,
             dx,
             dy,
             dz,
             v->wnx * dx + v->wny * dy + v->wnz * dz);
    }
  }
  sse_spring /= 2;
  return (sse_spring);
}


/*!
  \fn double SseTerms_MRIS::RepulsiveEnergy( double l_repulse, MHT *mht, MHT *mht_faces)
  \brief The repulsive term causes vertices to push away from each
  other based on the distance in 3D space (does not apply to nearest
  neighbors). This helps to prevent self-intersection. The force is
  inversely proportional to the distance to the 6th power (hidden
  parameter). Sets v->{dx,dy,dz}. 
  Hidden parameters:
    REPULSE_K - scaling term
    REPULSE_E - sets minimum distance
    4 - scaling term
*/
double SseTerms_MRIS::RepulsiveEnergy( double l_repulse, MHT *mht, MHT *mht_faces)
{
  int vno, num, min_vno, i, n;
  float dist, dx, dy, dz, x, y, z, min_d;
  double sse_repulse, v_sse;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }
  
  long hash_count = -32768-11000-8192-1101-105, hash_limit = 1, hash_limit_max = 10;  
  auto hash = fnv_init();

  min_d = 1000.0;
  min_vno = 0;
  sse_repulse = 0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) 
      continue;

    x = v->x;
    y = v->y;
    z = v->z;

    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket)
      continue;
    
    MHB *bin;
    for (v_sse = 0.0, bin = bucket->bins, num = i = 0; i < bucket->nused; i++, bin++) {

      bool const debugBin = debugNonDeterminism && (0 <= hash_count && hash_count < hash_limit_max);
      if (debugBin) {
        fprintf(stdout, "%s:%d sse_repulse bin->fno:%d\n",__FILE__,__LINE__,bin->fno);
      }

      /* don't be repelled by myself */
      if (bin->fno == vno)
        continue; 

      /* don't be repelled by a neighbor */
      for (n = 0; n < vt->vtotal; n++){
        if (vt->v[n] == bin->fno) {
          break;
        }
      }
      if (n < vt->vtotal) {
        if (debugBin) {
          fprintf(stdout, "%s:%d sse_repulse bin->fno:%d n:%d is nbr\n",__FILE__,__LINE__,bin->fno,n);
        }
        continue;
      }
      
      VERTEX const * const vn = &mris->vertices[bin->fno];
      if (!vn->ripflag) {
        dx = vn->x - x;
        dy = vn->y - y;
        dz = vn->z - z;
        dist = sqrt(dx * dx + dy * dy + dz * dz) + REPULSE_E;
        if (debugBin) {
          fprintf(stdout, "%s:%d sse_repulse bin->fno:%d n:%d dist:%f \n",__FILE__,__LINE__,bin->fno,n,dist);
        }

        if (vno == Gdiag_no) {
          if (dist - REPULSE_E < min_d) {
            min_vno = bin->fno;
            min_d = dist - REPULSE_E;
          }
        }

        dist = dist * dist * dist;
        dist *= dist; /* dist^6 */
	// dist = pow(dist,6.0); 
        v_sse += REPULSE_K / dist;
        if (debugBin) {
          fprintf(stdout, "%s:%d sse_repulse bin->fno:%d n:%d adding:%f\n",__FILE__,__LINE__,bin->fno,n,REPULSE_K/dist);
        }
      }
    } // loop over bucket

    sse_repulse += v_sse; // does not divide by the number of bins

    if (vno == Gdiag_no && !FZERO(v_sse)) {
      printf("v %d: repulse sse:    min_dist=%2.4f, v_sse %2.4f\n", vno, min_d, v_sse);
    }
    
    if (debugNonDeterminism && (hash_count < hash_limit_max)) {
      hash = fnv_add(hash, (unsigned char*)&sse_repulse, sizeof(sse_repulse));
      if (++hash_count >= hash_limit) {
        hash_limit += 1;
        fprintf(stdout, "%s:%d sse_repulse hash_count:%ld hash:%ld\n",__FILE__,__LINE__,hash_count,hash);
      }
    }
    
    MHTrelBucket(&bucket);
  }

  if (debugNonDeterminism) {
    fprintf(stdout, "%s:%d sse_repulse hash_count:%ld hash:%ld\n",__FILE__,__LINE__,hash_count,hash);
  }
  
  return (l_repulse * sse_repulse);
}


double SseTerms_MRIS::NonlinearSpringEnergy( INTEGRATION_PARMS *parms)
{
  int vno, n;
  double area_scale, sse_spring, E, F, f, rmin, rmax, ftotal;
  float dx, dy, dz, nc, r, lsq, mean_vdist;

  mean_vdist = MRIScomputeVertexSpacingStats(mris, NULL, NULL, NULL, NULL, NULL, CURRENT_VERTICES);
  lsq = mean_vdist * mean_vdist;
#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  rmin = parms->rmin;
  rmax = parms->rmax;
  if (FZERO(rmin) || FZERO(rmax))
    ErrorReturn(ERROR_BADPARM, (ERROR_BADPARM, "mrisComputeNonlinearSpringEnergy: rmin or rmax = 0!"));

  F = 6.0 / (1.0 / rmin - 1.0 / rmax);
  E = (1.0 / rmin + 1.0 / rmax) / 2;
  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    for (ftotal = r = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - v->x;
      dy = vn->y - v->y;
      dz = vn->z - v->z;
      //      lsq = dx*dx + dy*dy + dz*dz ;
      nc = dx * v->nx + dy * v->ny + dz * v->nz;
      dx = nc * v->nx;
      dy = nc * v->ny;
      dz = nc * v->nz;  // sn
      r = lsq / fabs(2.0 * nc);
      f = (1 + tanh(F * (1.0 / r - E)));
      ftotal += f * f;
    }
    if (vno == Gdiag_no) {
      printf("E_nlspring: f = %2.3f\n", ftotal / vt->vnum);
    }
    sse_spring += area_scale * ftotal / vt->vnum;
  }
  return (sse_spring);
}

double SseTerms_MRIS::SurfaceRepulsionEnergy( double l_repulse, MHT *mht)
{
  int vno, max_vno, i;
  float dx, dy, dz, x, y, z, sx, sy, sz, norm[3], dot;
  float max_scale, max_dot;
  double scale, sse;
  VERTEX *v, *vn;

  if (FZERO(l_repulse)) {
    return (NO_ERROR);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    x = v->x;
    y = v->y;
    z = v->z;
    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) {
      continue;
    }
    sx = sy = sz = 0.0;
    max_dot = max_scale = 0.0;
    max_vno = 0;
    MHB *bin = bucket->bins;
    for (i = 0; i < bucket->nused; i++, bin++) {
      vn = &mris->vertices[bin->fno];
      if (bin->fno == Gdiag_no) {
        DiagBreak();
      }
      if (vn->ripflag) {
        continue;
      }
      dx = x - vn->origx;
      dy = y - vn->origy;
      dz = z - vn->origz;
      mrisComputeOrigNormal(mris, bin->fno, norm);
      dot = dx * norm[0] + dy * norm[1] + dz * norm[2];
      if (dot > 1) {
        continue;
      }
      if (dot < 0 && vno == Gdiag_no) {
        DiagBreak();
      }
      if (dot > MAX_NEG_RATIO) {
        dot = MAX_NEG_RATIO;
      }
      else if (dot < -MAX_NEG_RATIO) {
        dot = -MAX_NEG_RATIO;
      }
#if 0
      scale = l_repulse / (1.0+exp(NEG_AREA_K*dot)) ;
#else
      scale = l_repulse * pow(1.0 - (double)dot, 4.0);
#endif
      if (scale > max_scale) {
        max_scale = scale;
        max_vno = bin->fno;
        max_dot = dot;
      }
      sx += (scale * v->nx);
      sy += (scale * v->ny);
      sz += (scale * v->nz);
    }

    sse += (sx * sx + sy * sy + sz * sz);
    if (vno == Gdiag_no) fprintf(stdout, "v %d inside repulse energy %2.3f\n", vno, (sx * sx + sy * sy + sz * sz));
    
    MHTrelBucket(&bucket);
  }
  return (sse);
}


double SseTerms_MRIS::AshburnerTriangleEnergy(
                                                 double l_ashburner_triangle,
                                                 INTEGRATION_PARMS *parms)
{
  int vno;
  double sse_ashburner;

  if (FZERO(l_ashburner_triangle)) return (0.0);

  mrisAssignFaces(mris, (MHT *)(parms->mht), CANONICAL_VERTICES);  // don't look it up every time
  sse_ashburner = 0.0;

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse_ashburner)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    double sse;
    VERTEX *v;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) ROMP_PF_continue;

    sse = mrisSampleAshburnerTriangleEnergy(mris, vno, parms, v->x, v->y, v->z);
    if (sse < 0) DiagBreak();

    sse_ashburner += sse;
    if (vno == Gdiag_no) printf("E_ash_triangle: vno = %d, E = %f\n", vno, sse);
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  sse_ashburner /= 2;
  return (sse_ashburner);
}

double SseTerms_MRIS::HistoNegativeLikelihood( INTEGRATION_PARMS *parms)
{
  double likelihood, entropy;
  int x, y, z, label;
  VECTOR *v_brain, *v_hires;
  MATRIX *m_brain_to_hires;
  double xv, yv, zv, val, pval;

  if (DZERO(parms->l_histo)) return (NO_ERROR);

  mrisCreateLikelihoodHistograms(mris, parms);

  v_brain = VectorAlloc(4, MATRIX_REAL);
  v_hires = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_brain, 4) = VECTOR_ELT(v_hires, 4) = 1.0;
  m_brain_to_hires = MRIgetVoxelToVoxelXform(parms->mri_brain, parms->mri_labels);
  for (likelihood = 0, x = 0; x < parms->mri_brain->width; x++) {
    V3_X(v_brain) = x;
    for (y = 0; y < parms->mri_brain->height; y++) {
      V3_Y(v_brain) = y;
      for (z = 0; z < parms->mri_brain->height; z++) {
        V3_Z(v_brain) = z;
        MatrixMultiply(m_brain_to_hires, v_brain, v_hires);
        xv = V3_X(v_hires);
        yv = V3_Y(v_hires);
        zv = V3_Z(v_hires);
        if (MRIindexNotInVolume(parms->mri_labels, xv, yv, zv)) {
          label = MRI_NONBRAIN;
        }
        else {
          label = MRIgetVoxVal(parms->mri_labels, nint(xv), nint(yv), nint(zv), 0);
        }

        val = MRIgetVoxVal(parms->mri_brain, x, y, z, 0);
        switch (label) {
          default:
          case MRI_NONBRAIN:
            pval = HISTOgetCount(parms->h_nonbrain, val);
            break;
          case MRI_WHITE_INTERIOR:
            pval = HISTOgetCount(parms->h_wm, val);
            break;
          case MRI_PIAL_INTERIOR:
            pval = HISTOgetCount(parms->h_gm, val);
            break;
        }
        likelihood += pval;
      }
    }
  }

  entropy = HISTOgetEntropy(parms->h_nonbrain) + HISTOgetEntropy(parms->h_gm);
  MatrixFree(&v_brain);
  MatrixFree(&v_hires);
  MatrixFree(&m_brain_to_hires);

  return (entropy);
}


double SseTerms_MRIS::NegativeLogPosterior( INTEGRATION_PARMS *parms, int *pnvox)
{
  MRI *mri = parms->mri_brain;
  double sse = 0.0, ll, wm_frac, gm_frac, out_frac, Ig, Ic, pval;
  float vmin, vmax, val;
  HISTOGRAM *hin, *hout;
  int x, y, z, nvox, label;
  MRI *mri_ll = NULL;
  static double last_sse = 0.0;

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    mri_ll = MRIcloneDifferentType(mri, MRI_FLOAT);
    int req = snprintf(fname, STRLEN, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(parms->mri_volume_fractions, fname);
  }

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }
  if (parms->hgm == NULL) {
    printf("creating intensity histograms\n");
    parms->hgm = hin = HISTOinit(parms->hgm, 256, (double)vmin, (double)vmax);
    parms->hout = hout = HISTOinit(parms->hout, 256, (double)vmin, (double)vmax);
    MRISclearMarks(mris);

    // build histogram estimates of PDFs of interior and exterior of ribbon
    for (x = 0; x < mri->width; x++)
      for (y = 0; y < mri->height; y++)
        for (z = 0; z < mri->depth; z++) {
          if (Gx == x && Gy == y && Gz == z) DiagBreak();
          if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
          val = MRIgetVoxVal(mri, x, y, z, 0);
          if (FZERO(val)) continue;

          wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
          gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1);
          out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2);
          if (parms->mri_aseg) {
            label = MRIgetVoxVal(parms->mri_aseg, x, y, z, 0);
            if (FZERO(gm_frac) && IS_CORTEX(label))  // aseg thinks it is but outside ribbon - ambiguous
              continue;
          }
          HISTOaddFractionalSample(hout, val, 0, 0, out_frac);
          HISTOaddFractionalSample(hin, val, 0, 0, gm_frac);
        }

    HISTOmakePDF(hin, hin);
    HISTOmakePDF(hout, hout);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      int req = snprintf(fname, STRLEN, "hin.%s.%3.3d.plt", 
			 parms->base_name, parms->t);  
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      HISTOplot(hin, fname);
      req = snprintf(fname, STRLEN, "hout.%s.%3.3d.plt", parms->base_name, parms->t);  
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      HISTOplot(hout, fname);
    }
  }
  else  // use previously computed ones
  {
    hin = parms->hgm;
    hout = parms->hout;
  }

  for (sse = 0.0, nvox = x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        if (Gx == x && Gy == y && Gz == z) DiagBreak();
        if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
        val = MRIgetVoxVal(mri, x, y, z, 0);
        if (FZERO(val)) continue;
        wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
        if (wm_frac > 0) continue;

        gm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 1);
        out_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 2);

        if (FZERO(out_frac))  // all gm
          pval = HISTOgetCount(hin, val);
        else if (FZERO(gm_frac))
          pval = HISTOgetCount(hout, val);
        else  // partial volume voxel
        {
          for (pval = 0.0, Ig = 0; Ig <= 256; Ig++) {
            Ic = (val - gm_frac * Ig) / out_frac;
            pval += HISTOgetCount(hout, Ic) * HISTOgetCount(hin, Ig);
          }
        }

        if (pval > 1 || pval < 0) DiagBreak();
        if (DZERO(pval)) pval = 1e-20;
        ll = -log10(pval);
        sse += ll;
        nvox++;
        if (mri_ll) MRIsetVoxVal(mri_ll, x, y, z, 0, ll);
        if (Gx == x && Gy == y && Gz == z)
          printf("voxel(%d, %d, %d) = %d, vfracs = (%2.1f, %2.1f, %2.1f), ll = %2.1f\n",
                 x,
                 y,
                 z,
                 nint(val),
                 wm_frac,
                 gm_frac,
                 out_frac,
                 ll);
      }

  if (!FZERO(last_sse) && sse > last_sse) DiagBreak();
  if (mri_ll) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.ll.%4.4d.mgz", parms->base_name, parms->t); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing log likelihood volume to %s\n", fname);
    MRIwrite(mri_ll, fname);
    MRIfree(&mri_ll);
  }
  //  HISTOfree(&hin) ; HISTOfree(&hout) ;
  if (pnvox) *pnvox = nvox;
  if (Gdiag_no >= 0) printf("E_logPosterior, %3.3d: %.1f (nvox=%d)\n", parms->t, sse, nvox);

  last_sse = sse;  // diagnostics
  return (sse);
}


double SseTerms_MRIS::NegativeLogPosterior2D( INTEGRATION_PARMS *parms, int *pnvox)
{
  MRI *mri = parms->mri_brain;
  MHT *mht;
  double sse = 0.0, ll, wm_frac, pval, dist, vdist, xs, ys, zs;
  float vmin, vmax, val;
  int x, y, z, nvox, label, xd, yd, zd, vno;
  VERTEX *v;
  MRI *mri_ll = NULL;
  static double last_sse = 0.0;
  HISTOGRAM2D *hs;
  MATRIX *m_vox2vox;
  VECTOR *v1, *v2;

  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 10);

  if (parms->mri_dtrans == NULL) {
    int nvox;
    MRI *mri_tmp;

    nvox = nint(
        MAX(MAX((mri->xsize / parms->resolution), (mri->ysize / parms->resolution)), (mri->zsize / parms->resolution)));
    mri_tmp = MRIcloneDifferentType(mri, MRI_FLOAT);
    parms->mri_dtrans = MRIupsampleN(mri_tmp, NULL, nvox);
    MRIfree(&mri_tmp);
    MRIScomputeDistanceToSurface(mris, parms->mri_dtrans, parms->resolution);
    if (parms->t > 0) DiagBreak();
    DiagBreak();
  }

  m_vox2vox = MRIgetVoxelToVoxelXform(mri, parms->mri_dtrans);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    mri_ll = MRIcloneDifferentType(mri, MRI_FLOAT);
    int req = snprintf(fname, STRLEN, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(parms->mri_volume_fractions, fname);
  }

  MRIvalRange(mri, &vmin, &vmax);
  if (mri->type == MRI_UCHAR) {
    vmin = 0;
    vmax = 255;
  }

  if (parms->h2d == NULL) {
    printf("creating 2D intensity histograms\n");
    parms->h2d = HISTO2Dinit(parms->h2d, 128, 101, (double)vmin, (double)vmax, -10, 10);
    MRISclearMarks(mris);

    // build histogram estimates of PDFs of interior and exterior of ribbon
    for (x = 0; x < mri->width; x++)
      for (y = 0; y < mri->height; y++)
        for (z = 0; z < mri->depth; z++) {
          if (Gx == x && Gy == y && Gz == z) DiagBreak();
          if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
          val = MRIgetVoxVal(mri, x, y, z, 0);
          wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
          if (wm_frac > 0) continue;

          V3_X(v1) = x;
          V3_Y(v1) = y;
          V3_Z(v1) = z;
          MatrixMultiply(m_vox2vox, v1, v2);
          xd = nint(V3_X(v2));
          yd = nint(V3_Y(v2));
          zd = nint(V3_Z(v2));
          if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd)) continue;
          dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
          if (FZERO(val) && (dist > 1 || dist < 0))  // don't allow 0s inside ribbon or too far out
            continue;
          if (FZERO(val)) continue;
          if (FZERO(val) && dist < 0) DiagBreak();

          // update distance with continuum measure
          MRIvoxelToSurfaceRAS(mri, x, y, z, &xs, &ys, &zs);
          MHTfindClosestVertexGeneric(mht, xs, ys, zs, 10, 4, &vno, &vdist);
          v = (vno < 0) ? nullptr : &mris->vertices[vno];
          if (v != NULL)  // compute distance from surface in normal direction
          {
            double dx, dy, dz;

            dx = xs - v->x;
            dy = ys - v->y;
            dz = zs - v->z;
            dist = dx * v->nx + dy * v->ny + dz * v->nz;
          }

          if (parms->mri_aseg) {
            label = MRIgetVoxVal(parms->mri_aseg, x, y, z, 0);
            if (dist > 1 && IS_CORTEX(label))  // aseg thinks it is but outside ribbon - ambiguous
              continue;
          }
          HISTO2DaddSample(parms->h2d, val, dist, 0, 0, 0, 0);
        }
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];

      int req = snprintf(fname, STRLEN, "h.%s.%3.3d.plt", parms->base_name, parms->t);
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(parms->h2d, fname);
    }

    //    hs = HISTO2DsoapBubbleZeros(parms->h2d, NULL, 100) ;
    hs = HISTO2DmakePDF(parms->h2d, NULL);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      int req = snprintf(fname, STRLEN, "h.%s.%3.3d.pdf.plt", parms->base_name, parms->t); 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(hs, fname);
    }
    //    HISTO2DsmoothBins2(hs, parms->h2d, 5) ;
    if (parms->h2d_out) HISTO2Dfree(&parms->h2d_out);
    parms->h2d_out = HISTO2DsmoothAnisotropic(hs, NULL, .25, 1);
    HISTO2Dsmooth(hs, parms->h2d, 1);
    HISTO2DmakePDF(parms->h2d_out, hs);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      int req = snprintf(fname, STRLEN, "h.%s.%3.3d.pdf.smooth.plt", parms->base_name, parms->t); 
      if( req >= STRLEN ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(hs, fname);
    }
    HISTO2Dfree(&parms->h2d);
    parms->h2d = hs;

    if (Gx > 0) {
      int b1, b2;
      float val, c;

      x = Gx;
      y = Gy;
      z = Gz;
      V3_X(v1) = x;
      V3_Y(v1) = y;
      V3_Z(v1) = z;
      MatrixMultiply(m_vox2vox, v1, v2);
      xd = nint(V3_X(v2));
      yd = nint(V3_Y(v2));
      zd = nint(V3_Z(v2));
      if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd))
        dist = 10;
      else
        dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
      val = MRIgetVoxVal(mri, x, y, z, 0);
      b1 = HISTO2DfindBin1(parms->h2d, val);
      b2 = HISTO2DfindBin2(parms->h2d, dist);
      c = HISTO2DgetCount(parms->h2d, val, dist);
      printf("voxel (%d, %d, %d) = %2.0f, dist=%2.2f, bin=%d,%d, count=%f\n", Gx, Gy, Gz, val, dist, b1, b2, c);

      DiagBreak();
    }
  }
  else  // use previously computed ones
  {
  }

  for (sse = 0.0, nvox = x = 0; x < mri->width; x++)
    for (y = 0; y < mri->height; y++)
      for (z = 0; z < mri->depth; z++) {
        if (Gx == x && Gy == y && Gz == z) DiagBreak();
        if (Gx2 == x && Gy2 == y && Gz2 == z) DiagBreak();
        val = MRIgetVoxVal(mri, x, y, z, 0);
        wm_frac = MRIgetVoxVal(parms->mri_volume_fractions, x, y, z, 0);
        if (wm_frac > 0) continue;

        V3_X(v1) = x;
        V3_Y(v1) = y;
        V3_Z(v1) = z;
        MatrixMultiply(m_vox2vox, v1, v2);
        xd = nint(V3_X(v2));
        yd = nint(V3_Y(v2));
        zd = nint(V3_Z(v2));
        if (MRIindexNotInVolume(parms->mri_dtrans, xd, yd, zd))
          dist = 10;
        else
          dist = MRIgetVoxVal(parms->mri_dtrans, xd, yd, zd, 0);
        if (FZERO(val) && dist > 1) continue;
        if (FZERO(val)) continue;
        if (FZERO(val)) DiagBreak();

        pval = HISTO2DgetCount(parms->h2d, val, dist);
        if (pval > 1 || pval < 0) DiagBreak();
        if (DZERO(pval)) pval = 1e-20;
        ll = -log10(pval);
        if (!std::isfinite(ll)) DiagBreak();
        sse += ll;
        nvox++;
        if (mri_ll) MRIsetVoxVal(mri_ll, x, y, z, 0, ll);
        if (Gx == x && Gy == y && Gz == z) {
          printf("voxel(%d, %d, %d) = %d, dist=%2.2f, ll = %2.1f, bins = %d, %d\n",
                 x,
                 y,
                 z,
                 nint(val),
                 dist,
                 ll,
                 HISTO2DfindBin1(parms->h2d, val),
                 HISTO2DfindBin2(parms->h2d, dist));
        }
      }

  if (!FZERO(last_sse) && sse > last_sse) DiagBreak();
  if (mri_ll) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.ll.%4.4d.mgz", parms->base_name, parms->t);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing log likelihood volume to %s\n", fname);
    MRIwrite(mri_ll, fname);
    MRIfree(&mri_ll);
  }
  //  HISTOfree(&hin) ; HISTOfree(&hout) ;
  if (pnvox) *pnvox = nvox;
  if (Gdiag_no >= 0) printf("E_logPosterior, %3.3d: %.1f (nvox=%d)\n", parms->t, sse, nvox);

  MHTfree(&mht);
  MatrixFree(&m_vox2vox);
  VectorFree(&v1);
  VectorFree(&v2);
  last_sse = sse;  // diagnostics
  return (sse);
}


//===================================================================================================
// Error functions
//
double SseTerms_MRIS::DistanceError( INTEGRATION_PARMS *parms)
{
  if (!(mris->dist_alloced_flags & 1)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceError","should have computed distances already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeVertexDistances(mris);
    }
  }
  if (!(mris->dist_alloced_flags & 2)) {
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceError","should have computed dist_origs already")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      mrisComputeOriginalVertexDistances(mris);
    }
  }

  if (false) {
    fprintf(stdout, "%s:%d calling mrisCheckDistOrig\n", __FILE__, __LINE__);
    if (!mrisCheckDistOrig(mris)) 
      fprintf(stdout, "  failed mrisCheckDistOrig\n");
  }
  
  int max_v, max_n, err_cnt, max_errs;
  volatile int count_dist_orig_zeros = 0;
  double dist_scale, sse_dist, max_del;

#if METRIC_SCALE
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else if (mris->status == MRIS_PARAMETERIZED_SPHERE) {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }
  else
    dist_scale = mris->neg_area < mris->total_area ? sqrt(mris->orig_area / (mris->total_area - mris->neg_area))
                                                   : sqrt(mris->orig_area / mris->total_area);
#else
  dist_scale = 1.0;
#endif
  max_del = -1.0;
  max_v = max_n = -1;

  err_cnt = 0;
  max_errs = 100;

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
  int trial;
  double sse_dist_trial0;
  for (trial = 0; trial < 2; trial++) {
#endif

  sse_dist = 0.0;
  
  int const acceptableNumberOfZeros = 
    ( mris->status == MRIS_PARAMETERIZED_SPHERE
    ||mris->status == MRIS_SPHERE)
    ? mris->nvertices * mris->avg_nbrs * 0.01       // only the direction from 000 has to match
    : mris->nvertices * mris->avg_nbrs * 0.001;     // xyz has to match

  const char* vertexRipflags = MRISexportVertexRipflags(mris);  
    // since we have to read them a lot, get them into 100KB in the L2 cache
    // rather than reading them in 6.4MB of cache lines
  
#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_REPRODUCIBLE

  #define ROMP_VARIABLE       vno 
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse_dist
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse_dist ROMP_PARTIALSUM(0)
    
#else
  int vno;
  
  ROMP_PF_begin         // mris_register

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
  #pragma omp parallel for if(trial==0) reduction(+ : sse_dist)
#else
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse_dist)
#endif
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif    

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) ROMP_PF_continue;

    if (vno == Gdiag_no) DiagBreak();

#if NO_NEG_DISTANCE_TERM
    if (v->neg) ROMP_PF_continue;
#endif

    double v_sse = 0.0;

    int n;
    for (n = 0; n < vt->vtotal; n++) {
      int const vn_vno = vt->v[n];
      if (vertexRipflags[vn_vno]) continue;
      
#if NO_NEG_DISTANCE_TERM
      if (mris->vertices[vn_vno].neg) continue;

#endif
      float const dist_orig_n = !v->dist_orig ? 0.0 : v->dist_orig[n];
      
      if (dist_orig_n >= UNFOUND_DIST) continue;

      if (DZERO(dist_orig_n) && (count_dist_orig_zeros++ > acceptableNumberOfZeros)) {
        fprintf(stderr, "v[%d]->dist_orig[%d] = %f!!!!, count_dist_orig_zeros:%d\n", vno, n, dist_orig_n, count_dist_orig_zeros);
        fflush(stderr);
        DiagBreak();
        if (++err_cnt > max_errs) {
          fprintf(stderr, ">>> dump of head zeroes \n");
          int dump_vno,dump_n,dump_count_dist_orig_zeros = 0;
          for (dump_vno = 0; dump_vno < mris->nvertices; dump_vno++) {
            VERTEX_TOPOLOGY const * const dump_vt = &mris->vertices_topology[dump_vno];
            VERTEX          const * const dump_v  = &mris->vertices         [dump_vno];
            if (dump_v->ripflag) continue;
            for (dump_n = 0; dump_n < dump_vt->vtotal; dump_n++) {
              int const dump_vn_vno = dump_vt->v[n];
              if (vertexRipflags[dump_vn_vno]) continue;
              float const dump_dist_orig_n = !dump_v->dist_orig ? 0.0 : dump_v->dist_orig[n];
              if (!DZERO(dump_dist_orig_n)) continue;
              dump_count_dist_orig_zeros++;
              fprintf(stderr, "v[%d]->dist_orig[%d] = %f!!!!, count_dist_orig_zeros:%d\n", dump_vno, dump_n, dump_dist_orig_n, dump_count_dist_orig_zeros);
              if (dump_count_dist_orig_zeros > 30) goto dump_done;
            }
          }
          dump_done:;
          ErrorExit(ERROR_BADLOOP, "mrisComputeDistanceError: Too many errors!\n");
        }
      }

      double delta = dist_scale * v->dist[n] - dist_orig_n;
      if (parms->vsmoothness)
        v_sse += (1.0 - parms->vsmoothness[vno]) * (delta * delta);
      else
        v_sse += delta * delta;

      if (!std::isfinite(delta) || !std::isfinite(v_sse)) DiagBreak();
    }
    if (v_sse > 10000) DiagBreak();

    if (parms->dist_error) parms->dist_error[vno] = v_sse;

    sse_dist += v_sse;
    if (!std::isfinite(sse_dist) || !std::isfinite(v_sse)) DiagBreak();
    
#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_REPRODUCIBLE
    #undef sse_dist 
  #include "romp_for_end.h"
#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

#ifdef BEVIN_MRISCOMPUTEDISTANCEERROR_CHECK
    if (trial == 0) {
        sse_dist_trial0 = sse_dist;
    } else { 
        if (sse_dist_trial0 != sse_dist) {
            fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
               sse_dist_trial0, sse_dist, sse_dist_trial0-sse_dist);
        }
    }
  } // trial
#endif

  freeAndNULL(vertexRipflags);  

  /*fprintf(stdout, "max_del = %f at v %d, n %d\n", max_del, max_v, max_n) ;*/
  return (sse_dist);
}


double MRIScomputeCorrelationError(MRI_SURFACE *mris, MRI_SP *mrisp_template, int fno)
{
  INTEGRATION_PARMS parms;
  float error;

  if (!mrisp_template) {
    return (0.0);
  }

  parms.mrisp_template = mrisp_template;
  parms.l_corr = 1.0f;
  parms.frame_no = fno;
  error = mrisComputeCorrelationError(mris, &parms, 1);
  return (sqrt(error / (double)MRISvalidVertices(mris)));
}

double SseTerms_MRIS::CorrelationError( INTEGRATION_PARMS *parms, int use_stds)
{
  float l_corr;

  l_corr = parms->l_corr + parms->l_pcorr; /* only one will be nonzero */
  if (FZERO(l_corr)) {
    return (0.0);
  }

  double sse = 0.0;
  
#ifdef BEVIN_MRISCOMPUTECORRELATIONERROR_REPRODUCIBLE

  #define ROMP_VARIABLE       vno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nvertices
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse  ROMP_PARTIALSUM(0)
    
#else
  int vno;

  ROMP_PF_begin         // Important during mris_register
 
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse)
#endif

  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin

#endif
    
    bool const vertexTrace = trace && (vno == 0);
    
    VERTEX *v = &mris->vertices[vno];
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    if (v->ripflag) {
      ROMP_PF_continue;
    }

    double src, target, delta, std;
    float x, y, z;

    x = v->x;
    y = v->y;
    z = v->z;
#if 0
    src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0, vertexTrace) ;
#else
    src = v->curv;
#endif
    target = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no, vertexTrace);

#define DISABLE_STDS 0
#if DISABLE_STDS
    std = 1.0f;
#else
    std = MRISPfunctionValTraceable(parms->mrisp_template, mris->radius, x, y, z, parms->frame_no + 1, vertexTrace);
    std = sqrt(std);
    if (FZERO(std)) {
      std = DEFAULT_STD /*FSMALL*/;
    }
    if (!use_stds) {
      std = 1.0f;
    }
#endif
    delta = (src - target) / std;
    if (!std::isfinite(target) || !std::isfinite(delta)) {
      DiagBreak();
    }
    if (parms->geometry_error) {
      parms->geometry_error[vno] = (delta * delta);
    }
    if (parms->abs_norm) {
      sse += fabs(delta);
    }
    else {
      sse += delta * delta;
    }
#ifdef BEVIN_MRISCOMPUTECORRELATIONERROR_REPRODUCIBLE

    #undef sse
  #include "romp_for_end.h"

#else
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif

  return (sse);
}


/*!
  \fn double SseTerms_MRIS::IntensityError( INTEGRATION_PARMS *parms)
  \brief Computes the sum of the squares of the value at a vertex minus the v->val.
   Ignores ripped vertices or any with v->val<0. Does not normalize by the number
   of vertices. Basically same computation as mrisRmsValError() but that func
   does normalize.
*/
double SseTerms_MRIS::IntensityError( INTEGRATION_PARMS *parms)
{
  int vno,nhits;
  VERTEX *v;
  double val0, xw, yw, zw;
  double sse, del0;

  if (FZERO(parms->l_intensity))
    return (0.0f);

  nhits = 0;
  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0)
      continue;
    nhits++;
    // Sample mri_brain at vertex
    MRISvertexToVoxel(mris, v, parms->mri_brain, &xw, &yw, &zw);
    MRIsampleVolume(parms->mri_brain, xw, yw, zw, &val0);
    del0 = v->val - val0;
    sse += (del0 * del0);
  }
  //printf("mrisComputeIntensityError() %f %d\n",sse,nhits);
  return (sse);
}

/*!
  \fn double MRISpointSetLocationError(MRIS *surf, double weight, json *pPointSet, int ComputeGradient)
  \brief Computes the SSError based on the distance between the set of
  points (as saved by freeview) and the surface. For each point, the
  closest vertex is found. This algorithm leaves a lot to be
  desired. The same vertex can be mapped to multiple points. A vertex
  may be very close to a point, but it is not included because it is
  slightly further away than another vertex. Could/should operate over
  a neighborhood.
 */
double MRISpointSetLocationError(MRIS *surf, double weight, json *pPointSet, int ComputeGradient)
{
  if(weight==0) return(0);
  json PointSet = *pPointSet;
  double max_spacing,sse=0;
  int max_vno;

  bool bTkReg = (PointSet["vox2ras"].get<std::string>() == std::string("tkreg"));
  json js_pts = PointSet["points"];
  printf("MRISpointSetLocationError(): TkReg=%d point set has %d points w=%g grad=%d\n",
    bTkReg,(int)js_pts.size(),weight,ComputeGradient);

  // The computation of the hash may be inefficient here as the same
  // hash may have been computed as some other point in the surface placement process
  MRIScomputeVertexSpacingStats(surf, NULL, NULL, &max_spacing, NULL, &max_vno, CURRENT_VERTICES);
  MRIS_HASH_TABLE *hash = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, max_spacing);

  // The conversion to TkReg could/should be done when it is read in,
  // but there are generally only a few points, so not too costly, and I
  // have not figured out how to change the vox2ras string
  MRI *mri = NULL;
  if(!bTkReg) mri = MRIallocFromVolGeom(&(surf->vg), MRI_UCHAR, 1,1);

  for(int i = 0; i < PointSet["points"].size(); i++) {
    double x,y,z,dx,dy,dz,mag;
    float distance;
    x = js_pts[i]["coordinates"]["x"].get<float>();
    y = js_pts[i]["coordinates"]["y"].get<float>();
    z = js_pts[i]["coordinates"]["z"].get<float>();
    if(!bTkReg){ // convert to TkReg
      double sx, sy, sz;
      MRIRASToSurfaceRAS(mri, x, y, z, &sx, &sy, &sz);
      x = sx; y = sy; z = sz;
    }
    // Find the closest vertex to this point
    int vno = MHTfindClosestVertexNoXYZ(hash, surf, x, y, z, &distance);
    if (vno < 0){
      //printf("Failed to find closest vertex in hash, using brute force\n");
      vno = MRISfindClosestVertex(surf, x, y, z, &distance, CURRENT_VERTICES);
    }
    VERTEX *v = &(surf->vertices[vno]);
    // Distance between the current location and the target
    dx = v->x - x;
    dy = v->y - y;
    dz = v->z - z;
    mag = dx * dx + dy * dy + dz * dz;
    sse += mag;
    if(ComputeGradient) {
      // cf mrisComputeTargetLocationTerm() (also note LOCATION_MOVE_LEN not use here)
      // negative is used below because dx here is -dx in mrisCompute...()
      double norm = sqrt(mag);
      v->dx -= weight * dx / norm;
      v->dy -= weight * dy / norm;
      v->dz -= weight * dz / norm;
    }
    //printf("%3d %6d (%6.3f %6.3f %6.3f) (%6.3f %6.3f %6.3f) %6.4lf %6.4lf\n",
    //	   i,vno,x,y,z,v->x,v->y,v->z,mag,sse);
  }
  MHTfree(&hash);
  if(mri) MRIfree(&mri);
  printf(" targetpointset sse = %g\n",sse);
  return(sse);
}

double SseTerms_MRIS::TargetPointSetError( INTEGRATION_PARMS *parms)
{
  if(FZERO(parms->l_targetpointset)) return (0.0f);
  //double sse = MRISpointSetLocationError(mris, parms->l_targetpointset, parms->TargetPointSet, 0);
  double sse = parms->TargetPointSet->CostAndGrad(parms->l_targetpointset,0);
  return(sse);
}

double SseTerms_MRIS::TargetLocationError( INTEGRATION_PARMS *parms)
{
  int vno, max_vno;
  VERTEX *v;
  double dx, dy, dz;
  double sse, mag, max_mag, last_mag;
  static double last_error[500000];

  if (FZERO(parms->l_location)) return (0.0f);

  last_mag = max_mag = 0;
  max_vno = -1;
  sse = 0.0;
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) continue;

    if (vno == Gdiag_no) DiagBreak();

    dx = v->x - v->targx;
    dy = v->y - v->targy;
    dz = v->z - v->targz;

    mag = dx * dx + dy * dy + dz * dz;
    if (mag > 50) DiagBreak();

    if (!devFinite(mag)) DiagBreak();

    // Not sure what this bit of code is for since it does not affect
    // the output of the function. Maybe just a diagnosis.  This is
    // not thread safe because of the static, but maybe that does not
    // matter. 
    if (mag > last_error[vno]) {
      if (mag > max_mag) {
        last_mag = last_error[vno];
        max_mag = mag;
        max_vno = vno;
      }
    }
    last_error[vno] = mag;
    sse += mag;
  }

  if (last_mag > 0) DiagBreak();
  return (sse);
}


double SseTerms_MRIS::RmsDistanceError()
{
  INTEGRATION_PARMS parms;
  double rms;

  parms.l_location = 1;
  rms = mrisComputeTargetLocationError(mris, &parms);
  return (sqrt(rms / MRISvalidVertices(mris)));
}


double SseTerms_MRIS::IntensityGradientError( INTEGRATION_PARMS *parms)
{
  int vno;
  VERTEX *v;
  float x, y, z;
  double mag0, xw, yw, zw, dx, dy, dz;
  double sse, del0;

  if (FZERO(parms->l_grad)) {
    return (0.0f);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

// MRIworldToVoxel(parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
#if 0
    MRISsurfaceRASToVoxelCached(mris, parms->mri_smooth, x, y, z, &xw, &yw, &zw) ;
#else
    MRISsurfaceRASToVoxelCached(mris, parms->mri_smooth, x, y, z, &xw, &yw, &zw);
#endif
    MRIsampleVolumeGradient(parms->mri_smooth, xw, yw, zw, &dx, &dy, &dz);
    mag0 = sqrt(dx * dx + dy * dy + dz * dz);

    del0 = v->mean - mag0;
    sse += (del0 * del0);
  }

  return (sse);
}



double SseTerms_MRIS::SphereError( double l_sphere, double r0)
{
  int vno;
  double sse, x0, y0, z0;

  if (FZERO(l_sphere)) {
    return (0.0f);
  }

  x0 = (mris->xlo + mris->xhi) / 2.0f;
  y0 = (mris->ylo + mris->yhi) / 2.0f;
  z0 = (mris->zlo + mris->zhi) / 2.0f;

  // This code appears to have, but does not have, a numeric stability problem
  // Typical inputs are
  //  x:7.50467 y:-43.4641 z:-9.50775 r:45.1203 del:154.88
  // But typical results are
  //  sse:4.33023e+09
  // I tried summing in groups of 1024 numbers then summing the sums, and got the same answer in %g format
  //        /Bevin
  //
  sse = 0.0;
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) reduction(+ : sse)
#endif
  for (vno = 0; vno < mris->nvertices; vno++) {
    ROMP_PFLB_begin
    
    VERTEX *v;
    double del, x, y, z, r;

    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = (double)v->x - x0;
    y = (double)v->y - y0;
    z = (double)v->z - z0;
    r = sqrt(x * x + y * y + z * z);

    del = r0 - r;
    sse += (del * del);
#if 0
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d sphere term: (%2.3f, %2.3f, %2.3f)\n",
              vno, v->dx, v->dy, v->dz) ;
#endif
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  return (sse);
}

double SseTerms_MRIS::DuraError( INTEGRATION_PARMS *parms)
{
  double dura_thresh = parms->dura_thresh, sse;
  MRI *mri_dura = parms->mri_dura;
  int vno;
  VERTEX *v;
  float x, y, z;
  double val0, xw, yw, zw, delV;

  if (FZERO(parms->l_dura)) {
    return (0.0);
  }

  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag || v->val < 0) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }

    x = v->x;
    y = v->y;
    z = v->z;

    MRISvertexToVoxel(mris, v, mri_dura, &xw, &yw, &zw);
    MRIsampleVolume(mri_dura, xw, yw, zw, &val0);
    if (val0 < dura_thresh) {
      continue;  // no effect
    }

    delV = dura_thresh - val0;
    sse += delV * delV;
  }

  return (sse);
}


double SseTerms_MRIS::VectorCorrelationError( INTEGRATION_PARMS *parms, int use_stds)
{
  double src, target, sse, delta, std;
  VERTEX *v;
  int n, vno, fno;
  float x, y, z, l_corr;

  double *vals, *corrs, *sses;
  int nframes, *frames, *ind;

  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    nframes++;
  }
  if (!nframes) {
    return 0.0;
  }

  corrs = (double *)malloc(nframes * sizeof(double));
  sses = (double *)malloc(nframes * sizeof(double));
  ind = (int *)malloc(nframes * sizeof(int));

  vals = (double *)malloc(2 * nframes * sizeof(double)); /* include the variances */
  frames = (int *)malloc(2 * nframes * sizeof(int));
  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;
    frames[2 * nframes] = fno;
    frames[2 * nframes + 1] = fno + 1;
    ind[nframes] = n;
    corrs[nframes] = l_corr;
    nframes++;
  }

  memset(sses, 0, nframes * sizeof(double)); /* set the vals to zero */
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    x = v->x;
    y = v->y;
    z = v->z;
    /* get the template values (variance included) */
    MRISPfunctionVectorVals(parms->mrisp_template, mris, x, y, z, frames, 2 * nframes, vals);
#if DISABLE_STDS
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = 1.0f;
      delta = (src - target) / std;
      sses[n] += delta * delta;
    }
#else
    for (n = 0; n < nframes; n++) {
      src = ((VALS_VP *)v->vp)->vals[ind[n]];
      target = vals[2 * n];
      std = sqrt(vals[2 * n + 1]);
      if (FZERO(std)) {
        std = DEFAULT_STD /*FSMALL*/; /* to be checked */
      }
      if (!use_stds) {
        std = 1.0f;
      }

      if (parms->fields[ind[n]].type) {
        std = MAX(0.01, std);
      }

      delta = (src - target) / std;

      // if(vno==0) fprintf(stderr,"[%d : %f , %f , %f ]",n, src,target,std);

      sses[n] += delta * delta;
    }
#endif
  }
  sse = 0.0f;
  for (n = 0; n < nframes; n++) {
    // fprintf(stderr,"(%d,%f,%f -> %f)\n",n,corrs[n],sses[n],corrs[n]*sses[n]);
    parms->fields[ind[n]].sse = sses[n];
    sse += corrs[n] * sses[n];
  }

  free(corrs);
  free(sses);
  free(ind);
  free(frames);
  free(vals);

  return (sse);

#if 0
  sse=0.0f; /* compiler warnings */
  for (n=0 ; n < parms->nfields ; n++)
  {

    l_corr = parms->fields[n].l_corr+parms->fields[n].l_pcorr ;

    if (FZERO(l_corr))
    {
      continue;
    }

    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
      {
        continue ;
      }

      x = v->x ;
      y = v->y ;
      z = v->z ;
#if 0
      src = MRISPfunctionVal(parms->mrisp, mris, x, y, z, 0) ;
#else
      src = ((VALS_VP*)v->vp)->vals[n] ;
#endif
      target = MRISPfunctionVal(parms->mrisp_template, mris, x, y, z, fno) ;

#define DISABLE_STDS 0
#if DISABLE_STDS
      std = 1.0f ;
#else
      std = MRISPfunctionVal(parms->mrisp_template,mris,x,y,z,fno+1);
      std = sqrt(std) ;
      if (FZERO(std))
      {
        std = DEFAULT_STD /*FSMALL*/ ;
      }
      if (!use_stds)
      {
        std = 1.0f ;
      }
#endif
      delta = (src - target) / std ;
      if (!finite(target) || !finite(delta))
      {
        DiagBreak() ;
      }
      sse += l_corr * delta * delta ;
    }
  }
  return(sse) ;
#endif
}


double SseTerms_MRIS::ExpandwrapError( MRI *mri_brain, double l_expandwrap, double target_radius)
{
  int vno;
  double xw, yw, zw, x, y, z, val, dx, dy, dz, sse, error, dist;
  VERTEX *v;
  float min_val, max_val, target_val, delta;

  if (FZERO(l_expandwrap)) {
    return (NO_ERROR);
  }

  mrisComputeSurfaceDimensions(mris);
  MRIvalRange(mri_brain, &min_val, &max_val);
  target_val = (min_val + max_val) / 2;
  for (sse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    target_val = v->val;
    x = v->x;
    y = v->y;
    z = v->z;
    MRISsurfaceRASToVoxelCached(mris, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &val);
    delta = (val - target_val);
    if (val < 0.25 * target_val) {
      dx = x - mris->xctr;
      dy = y - mris->yctr;
      dz = z - mris->zctr;
      dist = sqrt(dx * dx + dy * dy + dz * dz);
      error = (target_radius - dist);
      sse += error * error;
    }
    else {
      error = 0.0;
    }

    if (vno == Gdiag_no)
      fprintf(stdout,
              "v %d expandwrap error: %2.3f, target %2.1f, "
              "MRI %2.1f, del=%2.1f, \n",
              vno,
              error,
              target_val,
              val,
              delta);
  }
  return (sse);
}

double SseTerms_MRIS::ShrinkwrapError( MRI *mri_brain, double l_shrinkwrap)
{
#if 0
  static int iter = 100 ;
  int    vno ;
  double   xw, yw, zw, x, y, z, val ;
  VERTEX *v ;
  float  min_val, max_val, target_val, error ;
  double sse ;

  MRIvalRange(mri_brain, &min_val, &max_val) ;
  target_val = (min_val + max_val) / 2 ;
  sse = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    sse += iter ;
  }
  iter-- ;
  return(sse) ;
#else
  return (0.0);
#endif
}



double SseTerms_MRIS::Error(
                               INTEGRATION_PARMS *parms,
                               float *parea_rms,
                               float *pangle_rms,
                               float *pcurv_rms,
                               float *pdist_rms,
                               float *pcorr_rms)
{
  double rms, sse_area, sse_angle, sse_curv, delta, area_scale, sse_dist, sse_corr;
  int ano, fno, ntriangles, total_neighbors;
  FACE *face;
  float nv;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  sse_angle = sse_area = 0.0;
  for (ntriangles = fno = 0; fno < mris->nfaces; fno++) {
    face = &mris->faces[fno];
    if (face->ripflag) {
      continue;
    }
    FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);
    ntriangles++;
    delta = (double)(area_scale * face->area - fNorm->orig_area);
    sse_area += delta * delta;
    for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
      delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
      sse_angle += delta * delta;
    }
    if (!std::isfinite(sse_area) || !std::isfinite(sse_angle))
      ErrorExit(ERROR_BADPARM, "sse (%f, %f) is not finite at face %d!\n", sse_area, sse_angle, fno);
  }

  sse_corr = mrisComputeCorrelationError(mris, parms, 1);
  
  if (!DZERO(parms->l_dist)) {
    sse_dist = mrisComputeDistanceError(mris, parms);
  }
  else {
    sse_dist = 0;
  }

#if 0
  if (!FZERO(parms->l_spring))
  {
    sse_spring = mrisComputeSpringEnergy(mris) ;
  }
  if (!FZERO(parms->l_tspring))
  {
    sse_tspring = mrisComputeTangentialSpringEnergy(mris) ;
  }
#endif
  sse_curv = mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv);

  total_neighbors = mrisCountTotalNeighbors(mris);

  nv = (float)MRISvalidVertices(mris);
  if (mris->status != MRIS_PLANE) {
    *pcurv_rms = (float)sqrt(sse_curv / nv);
  }
  else {
    *pcurv_rms = 0.0f;
  }
  *pdist_rms = (float)sqrt(sse_dist / (double)total_neighbors);
  *parea_rms = (float)sqrt(sse_area / (double)ntriangles);
  *pangle_rms = (float)sqrt(sse_angle / (double)(ntriangles * ANGLES_PER_TRIANGLE));
  *pcorr_rms = (float)sqrt(sse_corr / (double)nv);

  rms = MRIScomputeSSE(mris, parms);

#if 0
  rms =
    *pdist_rms * parms->l_dist +
    *parea_rms * parms->l_area +
    *parea_rms * parms->l_parea +
    *pangle_rms * parms->l_angle +
    *pcorr_rms * (parms->l_corr+parms->l_pcorr) +
    sqrt(sse_spring/nv) * parms->l_spring ;
#endif
  return (rms);
}


int MRIScomputeDistanceErrors(MRI_SURFACE *mris, int nbhd_size, int max_nbrs)
{
  int vno, n, nvertices;
  double dist_scale, pct, dist, odist, mean, mean_error, smean, total_mean_error, total_mean;

  MRIScomputeMetricProperties(mris);
  if (mris->patch) {
    dist_scale = 1.0;
  }
  else {
    dist_scale = sqrt(mris->orig_area / mris->total_area);
  }

  total_mean = total_mean_error = mean = 0.0;
  for (pct = 0.0, nvertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      v->val = 0.0f;
      continue;
    }
    for (smean = mean_error = mean = 0.0, n = 0; n < vt->vtotal; n++) {
      nvertices++;
      dist = dist_scale * v->dist[n];
      odist = v->dist_orig[n];
#if 0
      mean += (dist - odist) * (dist - odist) ;
#else
      mean += odist;
#endif
      total_mean += odist;
      smean += dist - odist;
#define USE_FABS 1
#if USE_FABS
      mean_error += fabs(dist - odist);
      total_mean_error += fabs(dist - odist);
#else
      mean_error += (dist - odist) * (dist - odist);
      total_mean_error += (dist - odist) * (dist - odist);
#endif
      if (!FZERO(odist)) {
        pct += fabs(dist - odist) / odist;
      }
    }
    mean /= (double)vt->vtotal;
#if USE_FABS
    mean_error /= (double)vt->vtotal;
#else
    mean_error = sqrt(mean_error / (double)v->vtotal);
#endif
#if 0
    if (smean < 0.0f)
    {
      mean_error *= -1.0f ;
    }
#endif
    v->val = mean_error / mean;
  }

#if USE_FABS
  total_mean_error /= (double)nvertices;
#else
  total_mean_error = sqrt(total_mean_error / (double)nvertices);
#endif
  total_mean /= (double)nvertices;
  total_mean_error /= total_mean;
  fprintf(stdout, "mean dist = %2.3f, rms error = %2.2f%%\n", total_mean, 100.0 * total_mean_error);
  return (NO_ERROR);
}


//========================================
// misc
//
int mrisCreateLikelihoodHistograms(MRIS* mris, INTEGRATION_PARMS *parms)
{
  // Hard to convert to Surface because of saving an restoring vertex positions, cloning, etc.
  //
  int x, y, z, wlabel, plabel;
  VECTOR *v_brain, *v_hires;
  MATRIX *m_hires_to_brain;
  MRI *mri_pial;
  double xv, yv, zv, val, dist;

  if (parms->mri_white == NULL)  // white isn't moving, so only have to do it once
  {
    parms->mri_white = MRIupsampleN(parms->mri_brain, NULL, 3);
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, WHITE_VERTICES);
    MRISfillInterior(mris, parms->mri_white->xsize, parms->mri_white);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
  }
  mri_pial = MRIclone(parms->mri_white, NULL);
  MRISfillInterior(mris, mri_pial->xsize, mri_pial);
  if (parms->mri_labels == NULL) {
    parms->mri_labels = MRIclone(parms->mri_white, NULL);
  }
  parms->mri_dist = MRIdistanceTransform(mri_pial, parms->mri_dist, 1, 5 / mri_pial->xsize, DTRANS_MODE_SIGNED, NULL);

  parms->h_wm = HISTOinit(parms->h_wm, 256, 0, 255);
  parms->h_gm = HISTOinit(parms->h_gm, 256, 0, 255);
  parms->h_nonbrain = HISTOinit(parms->h_nonbrain, 256, 0, 255);
  m_hires_to_brain = MRIgetVoxelToVoxelXform(parms->mri_labels, parms->mri_brain);

  v_brain = VectorAlloc(4, MATRIX_REAL);
  v_hires = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v_brain, 4) = VECTOR_ELT(v_hires, 4) = 1.0;
  for (x = 0; x < mri_pial->width; x++) {
    V3_X(v_hires) = x;
    for (y = 0; y < mri_pial->height; y++) {
      V3_Y(v_hires) = y;
      for (z = 0; z < mri_pial->height; z++) {
        V3_Z(v_hires) = z;
        MatrixMultiply(m_hires_to_brain, v_hires, v_brain);
        xv = V3_X(v_brain);
        yv = V3_Y(v_brain);
        zv = V3_Z(v_brain);
        if (MRIindexNotInVolume(parms->mri_brain, xv, yv, zv)) {
          val = 0;
        }
        else {
          MRIsampleVolume(parms->mri_brain, xv, yv, zv, &val);
        }

        wlabel = MRIgetVoxVal(parms->mri_white, x, y, z, 0);
        plabel = MRIgetVoxVal(mri_pial, x, y, z, 0);
        dist = MRIgetVoxVal(parms->mri_dist, x, y, z, 0);
        if (dist > 3) {
          continue;  // don't consider the millions of voxels far from the surface
        }

        if (wlabel) {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_WHITE_INTERIOR);
          HISTOaddSample(parms->h_wm, val, 0, 255);
        }
        else if (plabel) {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_PIAL_INTERIOR);
          HISTOaddSample(parms->h_gm, val, 0, 255);
        }
        else {
          MRIsetVoxVal(parms->mri_labels, x, y, z, 0, MRI_NONBRAIN);
          HISTOaddSample(parms->h_nonbrain, val, 0, 255);
        }
      }
    }
  }
  HISTOmakePDF(parms->h_nonbrain, parms->h_nonbrain);
  HISTOmakePDF(parms->h_wm, parms->h_wm);
  HISTOmakePDF(parms->h_gm, parms->h_gm);
  MatrixFree(&m_hires_to_brain);
  MatrixFree(&v_brain);
  MatrixFree(&v_hires);
  MRIfree(&mri_pial);
  return (NO_ERROR);
}


double vlst_loglikelihood(MRIS *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM *hin, HISTOGRAM *hout)
{
  double ll = 0.0, dot, dx, dy, dz, pval, dist, Ig, Ic, gm_frac, out_frac;
  int i;
  float val;
  VERTEX *v;
  double xs, ys, zs;

  v = &mris->vertices[vno];
  xs = v->x + displacement * v->nx;
  ys = v->y + displacement * v->ny;
  zs = v->z + displacement * v->nz;
  for (i = 0; i < vl->nvox; i++) {
    dx = vl->xd[i] - xs;
    dy = vl->yd[i] - ys;
    dz = vl->zd[i] - zs;
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    dot = dx * v->nx + dy * v->ny + dz * v->nz;
    val = MRIgetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    if (dist < .5)  // distance to center<.5 --> distance to edge <1
    {
      if (dot > 0) {
        out_frac = dist + .5;
        gm_frac = 1 - out_frac;
      }
      else {
        gm_frac = dist + .5;
        out_frac = 1 - gm_frac;
      }
      for (pval = 0.0, Ig = 0; Ig <= 256; Ig++) {
        Ic = (val - gm_frac * Ig) / out_frac;
        pval += HISTOgetCount(hout, Ic) * HISTOgetCount(hin, Ig);
      }
    }
    else if (dot > 0)  // outside surface
      pval = HISTOgetCount(hout, val);
    else  // inside the surface
      pval = HISTOgetCount(hin, val);
    if (DZERO(pval)) pval = 1e-10;
    ll += -log(pval);
  }

  return (ll);
}


double vlst_loglikelihood2D(MRIS *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM2D *h, FILE *fp)
{
  double ll = 0.0, dot, dx, dy, dz, pval, dist;
  int i;
  float val;
  VERTEX *v;
  double xs, ys, zs;

  if (fp) fprintf(fp, "%f ", displacement);

  v = &mris->vertices[vno];
  xs = v->x + displacement * v->nx;
  ys = v->y + displacement * v->ny;
  zs = v->z + displacement * v->nz;
  for (i = 0; i < vl->nvox; i++) {
    dx = vl->xd[i] - xs;
    dy = vl->yd[i] - ys;
    dz = vl->zd[i] - zs;
    dist = sqrt(dx * dx + dy * dy + dz * dz);
    dot = dx * v->nx + dy * v->ny + dz * v->nz;
    val = MRIgetVoxVal(mri, vl->xi[i], vl->yi[i], vl->zi[i], 0);
    pval = HISTO2DgetCount(h, val, dot);
    if (DZERO(pval)) pval = 1e-10;
    if (fp) fprintf(fp, "%d %2.2f %2.2f ", (int)val, dot, -log(pval));
    ll += -log(pval);
  }

  if (fp) fprintf(fp, "\n");
  return (ll);
}


//==================================================================================
// MRIScomputeSSE and MRIScomputeSSEExternal are used for the numerical integration
//
// As such, they represent the exact error function being minimized, as opposed to computeError above.
//
// The performance critical ones are ...  
//
//      sse_area sse_neg_area       in this function                                here                iterates over faces computing face normals
//      sse_nl_area                 mrisComputeNonlinearAreaSSE(mris)               mrisurf_deform.c    iterates over faces doing a simple calculation
//      sse_dist                    mrisComputeNonlinearDistanceSSE(mris)           mrisurf_deform.c    iterates over vertices over their dist and dist_orig
//      sse_corr                    mrisComputeCorrelationError(mris, parms, 1)     mrisurf_deform.c    iterates over vertices using xyz calling mrisp.c MRISPfunctionValTraceable
//                                                                                                                          which does a lot more work than the others
//
// These are in the order the original code computed them, so that side effects are not reordered
// In older code the ashburner_triangle is computed but not used , here it is not computed at all

// The ELTS terms have a working overloading of mrisCompute### that can take a SurfaceFromMRIS_MP::XYZPositionConsequences::Surface as their first parameter
//      They are implemented below in template <class _Surface> struct SseTerms_DistortedSurfaces {...}
//
// The ELTM terms have a working overloading of mrisCompute### that can take a MRIS* as their first parameter
//      They also have an asserting overloading that can take a SurfaceFromMRIS_MP::XYZPositionConsequences::Surface as their first parameter
//      which will not be called because MRIScomputeSSE_canDo(MRIS_MP* usedOnlyForOverloadingResolution, INTEGRATION_PARMS *parms) returns false for these
//
#define SSE_TERMS \
      ELTM(sse_area                  , parms->l_parea,                            true,    computed_area                                                                   ) \
      ELTM(sse_neg_area              , parms->l_area,                             true,    computed_neg_area                                                               ) \
      ELTM(sse_repulse               , 1.0,                     (parms->l_repulse > 0),    mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current, mht_f_current)) \
      ELTM(sse_repulsive_ratio       , 1.0,                                       true,    mrisComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio)                   ) \
      ELTM(sse_tsmooth               , 1.0,                                       true,    mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms)             ) \
      ELTM(sse_thick_min             , parms->l_thick_min,                        true,    mrisComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms)         ) \
      ELTM(sse_ashburner_triangle    , parms->l_ashburner_triangle,               false,   mrisComputeAshburnerTriangleEnergy(mris, parms->l_ashburner_triangle, parms)    ) \
      ELTM(sse_thick_parallel        , parms->l_thick_parallel,                   true,    mrisComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms)        ) \
      ELTM(sse_thick_normal          , parms->l_thick_normal,                     true,    mrisComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms)            ) \
      ELTM(sse_thick_spring          , parms->l_thick_spring,                     true,    mrisComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms)            ) \
      ELTM(sse_nl_area               , parms->l_nlarea,        !FZERO(parms->l_nlarea),    mrisComputeNonlinearAreaSSE(mris)                                               ) \
      ELTM(sse_nl_dist               , parms->l_nldist,        !DZERO(parms->l_nldist),    mrisComputeNonlinearDistanceSSE(mris)                                           ) \
      ELTM(sse_dist                  , parms->l_dist,          !DZERO(parms->l_dist),      mrisComputeDistanceError(mris, parms)                                           ) \
      ELTM(sse_spring                , parms->l_spring,        !DZERO(parms->l_spring),    mrisComputeSpringEnergy(mris)                                                   ) \
      ELTM(sse_lap                   , parms->l_lap,           !DZERO(parms->l_lap),       mrisComputeLaplacianEnergy(mris)                                                ) \
      ELTM(sse_tspring               , parms->l_tspring,       !DZERO(parms->l_tspring),   mrisComputeTangentialSpringEnergy(mris)                                         ) \
      ELTM(sse_nlspring              , parms->l_nlspring,      !DZERO(parms->l_nlspring),  mrisComputeNonlinearSpringEnergy(mris, parms)                                   ) \
      ELTM(sse_curv                  , l_curv_scaled,          !DZERO(parms->l_curv),      mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv)                           ) \
      ELTM(sse_corr                  , l_corr,                 !DZERO(l_corr),             mrisComputeCorrelationError(mris, parms, 1)                                     ) \
      ELTM(sse_val                   , parms->l_intensity,     !DZERO(parms->l_intensity), mrisComputeIntensityError(mris, parms)                                          ) \
      ELTM(sse_loc                   , parms->l_location,      !DZERO(parms->l_location),  mrisComputeTargetLocationError(mris, parms)                                     ) \
      ELTM(sse_tps                   , parms->l_targetpointset,!DZERO(parms->l_targetpointset),  parms->TargetPointSet->CostAndGrad(parms->l_targetpointset, 0)            ) \
      ELTM(sse_dura                  , parms->l_dura,          !DZERO(parms->l_dura),      mrisComputeDuraError(mris, parms)                                               ) \
      ELTM(sse_histo                 , parms->l_histo,         !DZERO(parms->l_histo),     mrisComputeHistoNegativeLikelihood(mris, parms)                                 ) \
      ELTM(sse_map                   , parms->l_map,           !DZERO(parms->l_map),       mrisComputeNegativeLogPosterior(mris, parms, NULL)                              ) \
      ELTM(sse_map2d                 , parms->l_map2d,         !DZERO(parms->l_map2d),     mrisComputeNegativeLogPosterior2D(mris, parms, NULL)                            ) \
      ELTM(sse_grad                  , parms->l_grad,          !DZERO(parms->l_grad),      mrisComputeIntensityGradientError(mris, parms)                                  ) \
      ELTM(sse_sphere                , parms->l_sphere,        !DZERO(parms->l_sphere),    mrisComputeSphereError(mris, parms->l_sphere, parms->a)                         ) \
      ELTM(sse_shrinkwrap            , parms->l_shrinkwrap,    !DZERO(parms->l_shrinkwrap),mrisComputeShrinkwrapError(mris, parms->mri_brain, parms->l_shrinkwrap)         ) \
      ELTM(sse_expandwrap            , parms->l_expandwrap,    !DZERO(parms->l_expandwrap),mrisComputeExpandwrapError(mris, parms->mri_brain, parms->l_expandwrap, parms->target_radius)) \
      ELTM(sse_vectorCorrelationError, 1.0,                    use_multiframes,            mrisComputeVectorCorrelationError(mris, parms, 1)                               ) \
      // end of list

template <class Surface, class Some_MRIS>
double MRIScomputeSSE_template(Surface surface, Some_MRIS* mris, INTEGRATION_PARMS *parms)
{
  bool const debug = debugNonDeterminism;
  
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);
  double const l_curv_scaled    = (double)parms->l_curv * CURV_SCALE;
  double const area_scale =
#if METRIC_SCALE
    (surface.patch() || surface.noscale()) ? 1.0 : surface.orig_area() / surface.total_area();
#else
    1.0;
#endif

  double relevant_angle = 0, computed_neg_area = 0, computed_area = 0;

  if (!FZERO(parms->l_angle) || !FZERO(parms->l_area) || (!FZERO(parms->l_parea))) {

#ifdef BEVIN_MRISCOMPUTESSE_CHECK
    int trial; 
    double relevant_angle_trial0, computed_neg_area_trial0, computed_area_trial0;
    for (trial = 0; trial < 2; trial++) {

#endif

    relevant_angle = 0; computed_neg_area = 0; computed_area = 0;

    auto const nfaces = surface.nfaces();

#ifdef BEVIN_MRISCOMPUTESSE_REPRODUCIBLE

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             nfaces
    
  #define ROMP_SUMREDUCTION0  relevant_angle
  #define ROMP_SUMREDUCTION1  computed_neg_area
  #define ROMP_SUMREDUCTION2  computed_area
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define relevant_angle    ROMP_PARTIALSUM(0)
    #define computed_neg_area ROMP_PARTIALSUM(1)
    #define computed_area     ROMP_PARTIALSUM(2)

#else

    int fno;
    
    ROMP_PF_begin

#ifdef BEVIN_MRISCOMPUTESSE_CHECK
    #pragma omp parallel for if(trial==0) reduction(+ : relevant_angle, computed_neg_area, computed_area)
#else
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(fast) reduction(+ : relevant_angle, computed_neg_area, computed_area)
#endif
#endif
    for (fno = 0; fno < nfaces; fno++) {
      ROMP_PFLB_begin

#endif      
      auto face = surface.faces(fno);
      if (face.ripflag()) ROMP_PF_continue;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

      {
        auto const area = face.area();
        double const delta = (double)(area_scale * area - fNorm->orig_area);
#if ONLY_NEG_AREA_TERM
        if (area < 0.0f) computed_neg_area += delta * delta;
#endif
        computed_area += delta * delta;
      }
      
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        auto const angle = face.angle()[ano];
        double delta = deltaAngle(angle, face.orig_angle()[ano]);
#if ONLY_NEG_AREA_TERM
        if (angle >= 0.0f) delta = 0.0f;

#endif
        relevant_angle += delta * delta;
      }
      
      if (!std::isfinite(computed_area) || !std::isfinite(relevant_angle)) {
        ErrorExit(ERROR_BADPARM, "sse not finite at face %d!\n", fno);
      }
#ifdef BEVIN_MRISCOMPUTESSE_REPRODUCIBLE

    #undef relevant_angle
    #undef computed_neg_area
    #undef computed_area

  #include "romp_for_end.h"

#else
      ROMP_PFLB_end
    }
    ROMP_PF_end
#endif
    
#ifdef BEVIN_MRISCOMPUTESSE_CHECK

    if (trial == 0) {
       
      relevant_angle_trial0 = relevant_angle;
      computed_neg_area_trial0   = computed_neg_area;
      computed_area_trial0       = computed_area;
    } else { 
      if (relevant_angle_trial0 != relevant_angle) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           relevant_angle_trial0, relevant_angle, relevant_angle_trial0-relevant_angle);
      }
      if (computed_neg_area_trial0 != computed_neg_area) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           computed_neg_area_trial0, computed_neg_area, computed_neg_area_trial0-computed_neg_area);
      }
      if (computed_area_trial0 != computed_area) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           computed_area_trial0, computed_area, computed_area_trial0-computed_area);
      }
    }
    
    } // trial
#endif

  }

  MHT* mht_v_current = NULL;
  MHT* mht_f_current = NULL;
  if (!FZERO(parms->l_repulse)) {
    double vmean, vsigma;
    vmean = MRIScomputeTotalVertexSpacingStats     (mris, &vsigma, NULL, NULL, NULL, NULL);
    mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, vmean);
    mht_f_current = MHTcreateFaceTable_Resolution  (mris, CURRENT_VERTICES, vmean);
  }


#define ELTS(NAME, MULTIPLIER, COND, EXPR) double const NAME = (COND) ? (EXPR) : 0.0;
#define ELTM(NAME, MULTIPLIER, COND, EXPR) double const NAME = (COND) ? (EXPR) : 0.0;
    SSE_TERMS
#undef ELTM
#undef ELTS

  if (parms->l_thick_spring > 0 || parms->l_thick_min > 0 || parms->l_thick_parallel > 0 /* && DIAG_VERBOSE_ON*/)
    printf("min=%2.3f, parallel=%2.4f, normal=%2.4f, spring=%2.4f, ashburner=%2.3f, tsmooth=%2.3f\n",
           sse_thick_min            / (float)mris->nvertices,
           sse_thick_parallel       / (float)mris->nvertices,
           sse_thick_normal         / (float)mris->nvertices,
           sse_thick_spring         / (float)mris->nvertices,
           sse_ashburner_triangle   / (float)mris->nvertices,
           sse_tsmooth              / (float)mris->nvertices);
           
  double sse_init = 0;

  if (gMRISexternalSSE) {
    sse_init = (*gMRISexternalSSE)(mris, parms);
  }
  
  double sse = sse_init
#define ELTS(NAME, MULTIPLIER, COND, EXPR) + (MULTIPLIER) * (NAME)
#define ELTM(NAME, MULTIPLIER, COND, EXPR) + (MULTIPLIER) * (NAME)
    SSE_TERMS
#undef ELTM
#undef ELTS
    ;

  if (debug) {
    double sum = 0;
    #define ELTS(NAME, MULTIPLIER, COND, EXPR) fprintf(stdout, "new %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);
    #define ELTM(NAME, MULTIPLIER, COND, EXPR) fprintf(stdout, "new %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);
    ELTS(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sum);
    #undef ELTM
    #undef ELTS
  }
  
  // This code matches code Bevin added to the previous good code to compare old and new runs
  //
  static int logSSECount, logSSE;
  if (!logSSECount) { logSSE = !!getenv("FREESURFER_logSSE"); }
  logSSECount++;
  
  if (false || logSSE) {
    fprintf(stdout, "logSSE:%d \n", logSSECount);
    
    if (parms->l_dist) {
      bool dist_avail  = 
#ifdef COMPILING_MRIS_MP
        !!mris->v_dist[0];
#else
        !!(mris->dist_alloced_flags & 1);
#endif
      #define ELT(X) fprintf(stdout, " %s:%f\n", #X, (float)(X));
      ELT(dist_avail)
      if (dist_avail) {
        VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[0];
        VERTEX          const * const v  = &mris->vertices         [0];
        int n;
        for (n = 0; n < vt->vtotal; n++) {
          float const dist_n      = !v->dist      ? 0.0 : v->dist     [n];
          float const dist_orig_n = !v->dist_orig ? 0.0 : v->dist_orig[n];
          ELT(dist_n);
          ELT(dist_orig_n);
        }
      }
      ELT(mris->patch)
      ELT(mris->status)
      ELT(mris->orig_area)
      ELT(mris->total_area)
      ELT(mris->neg_area)
#undef ELT
    }

#define ELTS(NAME, MULTIPLIER, COND, EXPR) ELTM(NAME, MULTIPLIER, COND, EXPR)
#define ELTM(NAME, MULTIPLIER, COND, EXPR) \
    { double term = (MULTIPLIER) * (NAME); \
      if (term != 0.0) { fprintf(stdout, "new %s : %f \n", #NAME, term);  } \
    }
    ELTM(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sse);
#undef ELTM
#undef ELTS

  }
  //
  // end of Bevin added

  if (mht_v_current) MHTfree(&mht_v_current);
  if (mht_f_current) MHTfree(&mht_f_current);

  if (!devFinite(sse)) {
    DiagBreak();
  }

  return sse;
}

bool MRIScomputeSSE_canDo(MRIS* usedOnlyForOverloadingResolution, INTEGRATION_PARMS *parms)
{
  return true;
}

double MRIScomputeSSE(MRIS* mris, INTEGRATION_PARMS *parms)
{
  SurfaceFromMRIS::XYZPositionConsequences::Surface surface(mris);
  return MRIScomputeSSE_template(surface,mris,parms);    
}


bool MRIScomputeSSE_canDo(MRIS_MP* usedOnlyForOverloadingResolution, INTEGRATION_PARMS *parms)
{
#if 1
  return false;
#else
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);

  bool result = true;
#define ELTS(NAME, MULTIPLIER, COND, EXPR)
#define ELTM(NAME, MULTIPLIER, COND, EXPR) \
  if (COND) { static bool reported = false; \
    if (!reported) { reported = true; fprintf(stdout, "%s:%d can't do %s %s\n", __FILE__,__LINE__,#NAME,#EXPR); } \
    result = false; \
  }
  SSE_TERMS
  ELTM(sse_init,1.0,gMRISexternalSSE,)
#undef ELTM
#undef ELTS

  return result;
#endif
}


double MRIScomputeSSE(MRIS_MP* mris_mp, INTEGRATION_PARMS *parms)
{
#if 1
  return 0.0;
#else
  SurfaceFromMRIS_MP::XYZPositionConsequences::Surface surface(mris_mp);
  return MRIScomputeSSE_template(surface,mris_mp,parms);    
#endif
}

#undef SSE_TERMS

double MRIScomputeSSEExternal(MRIS* mris, INTEGRATION_PARMS *parms, double *ext_sse)
{
  double sse;

  if (gMRISexternalSSE) {
    sse = (*gMRISexternalSSE)(mris, parms);
  }
  else {
    sse = 0;
  }
  *ext_sse = sse;
  sse = MRIScomputeSSE(mris, parms); /* throw out ext_sse
                                        as it will be recomputed */

  return (sse);
}



// Generate all the jackets
//
#define MRIS_PARAMETER          MRIS* mris
#define MRIS_PARAMETER_COMMA    MRIS_PARAMETER ,
#define NOCOMMA_SELECTOR        int selector
#define COMMA_SELECTOR          , NOCOMMA_SELECTOR
#define SEP 
#define ELT(NAME, SIGNATURE, CALL)    double mrisCompute##NAME SIGNATURE { SseTerms_MRIS sseTerms(mris,selector); return sseTerms.NAME CALL; }
LIST_OF_SSETERMS
#undef ELT
#undef SEP
#undef COMMA_SELECTOR
#undef NOCOMMA_SELECTOR
#undef MRIS_PARAMETER_COMMA
#undef MRIS_PARAMETER

#define MRIS_PARAMETER          MRIS_MP* mris
#define MRIS_PARAMETER_COMMA    MRIS_PARAMETER ,
#define NOCOMMA_SELECTOR        int selector
#define COMMA_SELECTOR          , NOCOMMA_SELECTOR
#define SEP 
#define ELT(NAME, SIGNATURE, CALL)    double mrisCompute##NAME SIGNATURE { SseTerms_MRIS_MP sseTerms(mris,selector); return sseTerms.NAME CALL; }
LIST_OF_SSETERMS
#undef ELT
#undef SEP
#undef COMMA_SELECTOR
#undef NOCOMMA_SELECTOR
#undef MRIS_PARAMETER_COMMA
#undef MRIS_PARAMETER
