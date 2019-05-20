/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2019 The General Hospital Corporation (Boston, MA) "MGH"
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

// The SSE terms are either computed by iterating over the vertices or the faces
// Ideally there would be one pass over each, to 
//      a) minimize the scanning logic
//      b) minimize the cache traffic
// Furthermore it would be more efficient to do these using the dense MRIS_MP format rather than the MRIS format
// We are working towards this ideal...
//
// For now, only the terms used in recon-all bert have been moved here
// The rest will follow



// Utilities
//


//====================================================================================
// VERTEX SSE TERMS
//
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeNonlinearAreaSSE(MRI_SURFACE *mris)
{
  double area_scale;

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

  double sse;

#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
  int trial; 
  double sse_trial0;
  for (trial = 0; trial < 2; trial++) {
#endif

  sse = 0;
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_REPRODUCIBLE
  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  sse
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse  ROMP_PARTIALSUM(0)

#else
  int fno;
  
  ROMP_PF_begin     // mris_register
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
  #pragma omp parallel for if(trial==0) reduction(+ : sse)
#else
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(fast) reduction(+ : sse)
#endif
#endif
  for (fno = 0; fno < mris->nfaces; fno++) {
    ROMP_PFLB_begin

#endif
    
    double error, ratio;
    FACE *face;

    face = &mris->faces[fno];
    if (face->ripflag) {
      ROMP_PF_continue;
    }
#define SCALE_NONLINEAR_AREA 0
#if SCALE_NONLINEAR_AREA
    if (!FZERO(face->orig_area)) {
      ratio = area_scale * face->area / face->orig_area;
    }
    else {
      ratio = 0.0f;
    }
#else
    ratio = area_scale * face->area;
#endif
    if (ratio > MAX_NEG_RATIO) {
      ratio = MAX_NEG_RATIO;
    }
    else if (ratio < -MAX_NEG_RATIO) {
      ratio = -MAX_NEG_RATIO;
    }
#if 0
    error = (1.0 / NEG_AREA_K) * log(1.0+exp(-NEG_AREA_K*ratio)) ;
#else
    error = (log(1.0 + exp(NEG_AREA_K * ratio)) / NEG_AREA_K) - ratio;
#endif

    sse += error;
    if (!isfinite(sse) || !isfinite(error)) {
      ErrorExit(ERROR_BADPARM, "nlin area sse not finite at face %d!\n", fno);
    }
    
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_REPRODUCIBLE

    #undef sse
  #include "romp_for_end.h"
#else
  
    ROMP_PFLB_end
  }
  ROMP_PF_end
#endif
  
#ifdef BEVIN_MRISCOMPUTENONLINEARAREASSE_CHECK
    if (trial == 0) {
        sse_trial0 = sse;
    } else { 
        if (sse_trial0 != sse) {
            fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
               sse_trial0, sse, sse_trial0-sse);
        }
    }
  } // trial
#endif

  return (sse);
}


double mrisComputeDistanceError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
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
    switch (copeWithLogicProblem("FREESURFER_fix_mrisComputeDistanceTerm","should have computed dist_origs already")) {
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

      if (!isfinite(delta) || !isfinite(v_sse)) DiagBreak();
    }
    if (v_sse > 10000) DiagBreak();

    if (parms->dist_error) parms->dist_error[vno] = v_sse;

    sse_dist += v_sse;
    if (!isfinite(sse_dist) || !isfinite(v_sse)) DiagBreak();
    
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

  memset(&parms, 0, sizeof(parms));
  parms.mrisp_template = mrisp_template;
  parms.l_corr = 1.0f;
  parms.frame_no = fno;
  error = mrisComputeCorrelationError(mris, &parms, 1);
  return (sqrt(error / (double)MRISvalidVertices(mris)));
}

double mrisComputeCorrelationErrorTraceable(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace);

double mrisComputeCorrelationError(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds) {
    return mrisComputeCorrelationErrorTraceable(mris, parms, use_stds, false);
}

double mrisComputeCorrelationErrorTraceable(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace)
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
#define DEFAULT_STD 4.0f
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
    if (!isfinite(target) || !isfinite(delta)) {
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


double mrisComputeRepulsiveRatioEnergy(MRI_SURFACE *mris, double l_repulse)
{
  int vno, n;
  double sse_repulse, v_sse, dist, dx, dy, dz, x, y, z, canon_dist, cdx, cdy, cdz;

  if (FZERO(l_repulse))
    return (0.0);

  for (sse_repulse = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];

    if (v->ripflag) 
      continue;

    x = v->x;
    y = v->y;
    z = v->z;
    for(v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      if (!vn->ripflag) {
        dx = x - vn->x;
        dy = y - vn->y;
        dz = z - vn->z;
        dist = sqrt(dx * dx + dy * dy + dz * dz);
        cdx = vn->cx - v->cx;
        cdy = vn->cy - v->cy;
        cdz = vn->cz - v->cz;
        canon_dist = sqrt(cdx * cdx + cdy * cdy + cdz * cdz) + REPULSE_E;
        dist /= canon_dist;
        dist += REPULSE_E;
#if 0
        v_sse += REPULSE_K / (dist*dist*dist*dist) ;
#else
        v_sse += REPULSE_K / (dist * dist);
#endif
      }
    }
    sse_repulse += v_sse;
  }
  return (l_repulse * sse_repulse);
}


double mrisComputeSpringEnergy(MRI_SURFACE *mris)
{
  int vno, n;
  double area_scale, sse_spring, v_sse;

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

  for (sse_spring = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      v_sse += (v->dist[n] * v->dist[n]);
    }
    sse_spring += area_scale * v_sse;
  }
  return (sse_spring);
}


double mrisComputeThicknessMinimizationEnergy(MRIS *mris, double l_thick_min, INTEGRATION_PARMS *parms)
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


double mrisComputeThicknessParallelEnergy(MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms)
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


double mrisComputeThicknessSmoothnessEnergy(MRIS *mris, double l_tsmooth, INTEGRATION_PARMS *parms)
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


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static double big_sse = 10.0;
double mrisComputeThicknessNormalEnergy(MRIS *mris, double l_thick_normal, INTEGRATION_PARMS *parms)
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
/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double mrisComputeThicknessSpringEnergy(MRIS *mris, double l_thick_spring, INTEGRATION_PARMS *parms)
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


