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


#define MAX_VOXELS          mrisurf_sse_MAX_VOXELS
#define MAX_DISPLACEMENT    mrisurf_sse_MAX_DISPLACEMENT 
#define DISPLACEMENT_DELTA  mrisurf_sse_DISPLACEMENT_DELTA
#define DEFAULT_STD         mrisurf_sse_DEFAULT_STD


int mrisCreateLikelihoodHistograms(MRIS* mris, INTEGRATION_PARMS *parms)
{
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

struct SseTerms {
    MRIS* const mris;
    SseTerms(MRIS* const mris) : mris(mris) {}
    
    #define MRIS_PARAMETER          
    #define MRIS_PARAMETER_COMMA    
    #define SEP 
    #define ELT(NAME, SIGNATURE, CALL)    double NAME SIGNATURE;
    LIST_OF_SSETERMS
    #undef ELT
    #undef SEP
    #undef MRIS_PARAMETER_COMMA
    #undef MRIS_PARAMETER
};


/*-----------------------------------------------------
  Fit a 2-d quadratic to the surface locally and compute the SSE as
  the square of the constant term (the distance the quadratic fit surface
  is from going through the central vertex)
  ------------------------------------------------------*/
double SseTerms::QuadraticCurvatureSSE( double l_curv)            // BEVIN mris_make_surfaces 3
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
double SseTerms::NonlinearAreaSSE()
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


double SseTerms::RepulsiveRatioEnergy( double l_repulse)
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


double SseTerms::SpringEnergy()
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


double SseTerms::ThicknessMinimizationEnergy( double l_thick_min, INTEGRATION_PARMS *parms)
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


double SseTerms::ThicknessParallelEnergy( double l_thick_parallel, INTEGRATION_PARMS *parms)
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


double SseTerms::ThicknessSmoothnessEnergy( double l_tsmooth, INTEGRATION_PARMS *parms)
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
double SseTerms::ThicknessNormalEnergy( double l_thick_normal, INTEGRATION_PARMS *parms)
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
double SseTerms::ThicknessSpringEnergy( double l_thick_spring, INTEGRATION_PARMS *parms)
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
  \fn double SseTerms::RepulsiveEnergy( double l_repulse, MHT *mht, MHT *mht_faces)
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
double SseTerms::RepulsiveEnergy( double l_repulse, MHT *mht, MHT *mht_faces)
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

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description


  Note: this function assumes that the mris surface has the original
  (i.e. after global rotational alignment)
  spherical coordinates in the TMP2_VERTICES
  ------------------------------------------------------*/
double SseTerms::LaplacianEnergy()
{
  int vno, n;
  double area_scale, sse_lap, v_sse, dx, dy, dz, vx, vy, vz, vnx, vny, vnz, error;

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

  for (sse_lap = 0.0, vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }

    vx = v->x - v->t2x;
    vy = v->y - v->t2y;
    vz = v->z - v->t2z;
    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      vnx = vn->x - vn->t2x;
      vny = vn->y - vn->t2y;
      vnz = vn->z - vn->t2z;
      dx = vnx - vx;
      dy = vny - vy;
      dz = vnz - vz;
      error = dx * dx + dy * dy + dz * dz;
      v_sse += error;
    }
    sse_lap += area_scale * v_sse;
  }
  return (sse_lap);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
double SseTerms::TangentialSpringEnergy()
{
  int vno, n;
  double area_scale, sse_spring, v_sse;
  float dx, dy, dz, x, y, z, nc, dist_sq;

#if METRIC_SCALE
  if (mris->patch) {
    area_scale = 1.0;
  }
  else {
    area_scale = FZERO(mris->total_area) ? 1.0 : mris->orig_area / mris->total_area;
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

    x = v->x;
    y = v->y;
    z = v->z;

    for (v_sse = 0.0, n = 0; n < vt->vnum; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]];
      dx = vn->x - x;
      dy = vn->y - y;
      dz = vn->z - z;
      nc = dx * v->nx + dy * v->ny + dz * v->nz;
      dx -= nc * v->nx;
      dy -= nc * v->ny;
      dz -= nc * v->nz;
      dist_sq = dx * dx + dy * dy + dz * dz;
      v_sse += dist_sq;
    }
    sse_spring += area_scale * v_sse;
  }
  return (sse_spring);
}

double SseTerms::NonlinearSpringEnergy( INTEGRATION_PARMS *parms)
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

double SseTerms::SurfaceRepulsionEnergy( double l_repulse, MHT *mht)
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

double SseTerms::NonlinearDistanceSSE()
{
  int vno, n, nvertices, max_v, max_n;
  double dist_scale, sse_dist, delta, v_sse, max_del, ratio;

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
  for (sse_dist = 0.0, nvertices = vno = 0; vno < mris->nvertices; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag) {
      continue;
    }
    nvertices++;
    for (v_sse = 0.0, n = 0; n < vt->vtotal; n++) {
      if (FZERO(v->dist_orig[n])) {
        continue;
      }
      ratio = dist_scale * v->dist[n] / v->dist_orig[n];
      delta = log(1 + exp(ratio));
      v_sse += delta;
      if (!isfinite(delta) || !isfinite(v_sse)) {
        DiagBreak();
      }
    }
    sse_dist += v_sse;
    if (!isfinite(sse_dist) || !isfinite(v_sse)) {
      DiagBreak();
    }
  }

  return (sse_dist);
}


double SseTerms::AshburnerTriangleEnergy(
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

double SseTerms::HistoNegativeLikelihood( INTEGRATION_PARMS *parms)
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


double SseTerms::NegativeLogPosterior( INTEGRATION_PARMS *parms, int *pnvox)
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
    sprintf(fname, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t);
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
      sprintf(fname, "hin.%s.%3.3d.plt", parms->base_name, parms->t);
      HISTOplot(hin, fname);
      sprintf(fname, "hout.%s.%3.3d.plt", parms->base_name, parms->t);
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
    sprintf(fname, "%s.ll.%4.4d.mgz", parms->base_name, parms->t);
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


double SseTerms::NegativeLogPosterior2D( INTEGRATION_PARMS *parms, int *pnvox)
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
    sprintf(fname, "%s.vfrac.%4.4d.mgz", parms->base_name, parms->t);
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

      sprintf(fname, "h.%s.%3.3d.plt", parms->base_name, parms->t);
      printf("writing histogram %s\n", fname);
      HISTO2Dwrite(parms->h2d, fname);
    }

    //    hs = HISTO2DsoapBubbleZeros(parms->h2d, NULL, 100) ;
    hs = HISTO2DmakePDF(parms->h2d, NULL);
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      sprintf(fname, "h.%s.%3.3d.pdf.plt", parms->base_name, parms->t);
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
      sprintf(fname, "h.%s.%3.3d.pdf.smooth.plt", parms->base_name, parms->t);
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
        if (!isfinite(ll)) DiagBreak();
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
    sprintf(fname, "%s.ll.%4.4d.mgz", parms->base_name, parms->t);
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


/*-----------------------------------------------------
  Description
  MRIScomputeSSE and MRIScomputeSSEExternal
  are used for the numerical integration.
  As such, they represent the exact error function being minimized, as opposed to computeError above.
  ------------------------------------------------------*/
#include "mrisurf_deform_computeSSE.h"

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


//===================================================================================================
// Error functions
//
double SseTerms::DistanceError( INTEGRATION_PARMS *parms)
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

double SseTerms::CorrelationError( INTEGRATION_PARMS *parms, int use_stds)
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


/*!
  \fn double SseTerms::IntensityError( INTEGRATION_PARMS *parms)
  \brief Computes the sum of the squares of the value at a vertex minus the v->val.
   Ignores ripped vertices or any with v->val<0. Does not normalize by the number
   of vertices. Basically same computation as mrisRmsValError() but that func
   does normalize.
*/
double SseTerms::IntensityError( INTEGRATION_PARMS *parms)
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
/*! -----------------------------------------------------
  \fn static double SseTerms::TargetLocationError( INTEGRATION_PARMS *parms)
  \brief Computes the distance squared between v->{xyz} and v->targ{xyz} and sums up over
  all unripped vertices. See also mrisRmsDistanceError(mris).
  ------------------------------------------------------*/
double SseTerms::TargetLocationError( INTEGRATION_PARMS *parms)
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


double SseTerms::RmsDistanceError()
{
  INTEGRATION_PARMS parms;
  double rms;

  memset(&parms, 0, sizeof(parms));
  parms.l_location = 1;
  rms = mrisComputeTargetLocationError(mris, &parms);
  return (sqrt(rms / MRISvalidVertices(mris)));
}


double SseTerms::IntensityGradientError( INTEGRATION_PARMS *parms)
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



double SseTerms::SphereError( double l_sphere, double r0)
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

double SseTerms::DuraError( INTEGRATION_PARMS *parms)
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


double SseTerms::VectorCorrelationError( INTEGRATION_PARMS *parms, int use_stds)
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


double SseTerms::ExpandwrapError( MRI *mri_brain, double l_expandwrap, double target_radius)
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

double SseTerms::ShrinkwrapError( MRI *mri_brain, double l_shrinkwrap)
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


double SseTerms::Error(
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
    if (!isfinite(sse_area) || !isfinite(sse_angle))
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


// Generate all the jackets
//
#define MRIS_PARAMETER          MRIS* mris          
#define MRIS_PARAMETER_COMMA    MRIS_PARAMETER ,
#define SEP 
#define ELT(NAME, SIGNATURE, CALL)    double mrisCompute##NAME SIGNATURE { SseTerms sseTerms(mris); return sseTerms.NAME CALL; }
LIST_OF_SSETERMS
#undef ELT
#undef SEP
#undef MRIS_PARAMETER_COMMA
#undef MRIS_PARAMETER
