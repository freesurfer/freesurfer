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
#include "mrisurf_integrate.h"

#include "mrisurf_sseTerms.h"
#include "mrisurf_compute_dxyz.h"

#include "mrisurf_base.h"


static int mrisIntegrationEpoch     (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_avgs);
static double mrisLineMinimize      (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
static double mrisLineMinimizeSearch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

/*-----------------------------------------------------*/
int mrisLogIntegrationParms(FILE *fp, MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  char *cp, host_name[STRLEN];
  static int first = 1;

  if (!fp) {
    return (NO_ERROR);
  }

  cp = getenv("HOST");
  if (cp) {
    strcpy(host_name, cp);
  }
  else {
    strcpy(host_name, "unknown");
  }

  fprintf(fp,
          "tol=%2.1e, sigma=%2.1f, host=%5.5s, nav=%d, nbrs=%d",
          (float)parms->tol,
          parms->sigma,
          host_name,
          parms->n_averages,
          mris->nsize);
  if (!FZERO(parms->l_area)) {
    fprintf(fp, ", l_area=%2.3f", parms->l_area);
  }
  if (!FZERO(parms->l_map)) {
    fprintf(fp, ", l_map=%2.3f", parms->l_map);
  }
  if (!FZERO(parms->l_map2d)) {
    fprintf(fp, ", l_map2d=%2.3f", parms->l_map2d);
  }
  if (!FZERO(parms->l_dura)) {
    fprintf(fp, ", l_dura=%2.3f", parms->l_dura);
  }
  if (!FZERO(parms->l_expandwrap)) {
    fprintf(fp, ", l_expandwrap=%2.3f", parms->l_expandwrap);
  }
  if (!FZERO(parms->l_shrinkwrap)) {
    fprintf(fp, ", l_shrinkwrap=%2.3f", parms->l_shrinkwrap);
  }
  if (!FZERO(parms->l_external)) {
    fprintf(fp, ", l_extern=%2.3f", parms->l_external);
  }
  if (!FZERO(parms->l_parea)) {
    fprintf(fp, ", l_parea=%2.3f", parms->l_parea);
  }
  if (!FZERO(parms->l_nlarea)) {
    fprintf(fp, ", l_nlarea=%2.3f", parms->l_nlarea);
  }
  if (!FZERO(parms->l_nldist)) {
    fprintf(fp, ", l_nldist=%2.3f", parms->l_nldist);
  }
  if (!FZERO(parms->l_angle)) {
    fprintf(fp, ", l_angle=%2.3f", parms->l_angle);
  }
  if (!FZERO(parms->l_repulse)) {
    fprintf(fp, ", l_repulse=%2.3f", parms->l_repulse);
  }
  if (!FZERO(parms->l_repulse_ratio)) {
    fprintf(fp, ", l_repulse_ratio=%2.3f", parms->l_repulse_ratio);
  }
  if (!FZERO(parms->l_surf_repulse)) {
    fprintf(fp, ", l_surf_repulse=%2.3f", parms->l_surf_repulse);
  }
  if (!FZERO(parms->l_corr)) {
    fprintf(fp, ", l_corr=%2.3f", parms->l_corr);
  }
  if (!FZERO(parms->l_spring)) {
    fprintf(fp, ", l_spring=%2.3f", parms->l_spring);
  }
  if (!FZERO(parms->l_spring_norm)) {
    fprintf(fp, ", l_spring_norm=%2.3f", parms->l_spring_norm);
  }
  if (!FZERO(parms->l_tspring)) {
    fprintf(fp, ", l_tspring=%2.3f", parms->l_tspring);
  }
  if (!FZERO(parms->l_nltspring)) {
    fprintf(fp, ", l_nltspring=%2.3f", parms->l_nltspring);
  }
  if (!FZERO(parms->l_histo)) {
    fprintf(fp, ", l_histo=%2.3f", parms->l_histo);
  }
  if (!FZERO(parms->l_nlspring)) {
    fprintf(fp, ", l_nlspring=%2.3f", parms->l_nlspring);
  }
  if (!FZERO(parms->l_nspring)) {
    fprintf(fp, ", l_nspring=%2.3f", parms->l_nspring);
  }
  if (!FZERO(parms->l_dist)) {
    fprintf(fp, ", l_dist=%2.3f", parms->l_dist);
  }
  if (!FZERO(parms->l_intensity)) {
    fprintf(fp, ", l_intensity=%2.3f", parms->l_intensity);
  }
  if (!FZERO(parms->l_location)) {
    fprintf(fp, ", l_location=%2.3f", parms->l_location);
  }
  if (!FZERO(parms->l_grad)) {
    fprintf(fp, ", l_grad=%2.3f", parms->l_grad);
  }
  if (!FZERO(parms->l_sphere)) {
    fprintf(fp, ", l_sphere=%2.3f", parms->l_sphere);
  }
  if (!FZERO(parms->l_expand)) {
    fprintf(fp, ", l_expand=%2.3f", parms->l_expand);
  }
  if (!FZERO(parms->l_curv)) {
    fprintf(fp, ", l_curv=%2.3f", parms->l_curv);
  }
  if (!FZERO(parms->l_convex)) {
    fprintf(fp, ", l_convex=%2.3f", parms->l_convex);
  }
  if (!FZERO(parms->l_boundary)) {
    fprintf(fp, ", l_boundary=%2.3f", parms->l_boundary);
  }
  if (!FZERO(parms->l_neg)) {
    fprintf(fp, ", l_neg=%2.3f", parms->l_neg);
  }
  if (!FZERO(parms->l_tsmooth)) {
    fprintf(fp, ", l_tsmooth=%2.3f", parms->l_tsmooth);
  }
  if (!FZERO(parms->l_thick_min)) {
    fprintf(fp, ", l_thick_min=%2.3f", parms->l_thick_min);
  }
  if (!FZERO(parms->l_thick_parallel)) {
    fprintf(fp, ", l_thick_parallel=%2.3f", parms->l_thick_parallel);
  }
  if (!FZERO(parms->l_ashburner_triangle)) {
    fprintf(fp, ", l_ashburner_triangle=%2.3f", parms->l_ashburner_triangle);
  }
  if (!FZERO(parms->l_thick_normal)) {
    fprintf(fp, ", l_thick_normal=%2.3f", parms->l_thick_normal);
  }
  if (!FZERO(parms->l_thick_spring)) {
    fprintf(fp, ", l_thick_spring=%2.3f", parms->l_thick_spring);
  }
  fprintf(fp, "\n");
  switch (parms->integration_type) {
    case INTEGRATE_LM_SEARCH:
      fprintf(fp, "using binary search line minimization\n");
      break;
    case INTEGRATE_LINE_MINIMIZE:
      fprintf(fp, "using quadratic fit line minimization\n");
      break;
    case INTEGRATE_ADAPTIVE:
      fprintf(fp,
              "mom=%2.2f, dt=%2.2f, base_dt=%2.3f, dt_inc=%2.2f, "
              "dt_dec=%2.2f, err_rat=%2.2f\n",
              (float)parms->momentum,
              (float)parms->dt,
              (float)parms->base_dt,
              (float)parms->dt_increase,
              (float)parms->dt_decrease,
              (float)parms->error_ratio);
      break;
    default:
    case INTEGRATE_MOMENTUM:
      fprintf(fp, "mom=%2.2f, dt=%2.2f\n", (float)parms->momentum, (float)parms->dt);
      break;
  }
#if 0
  fprintf(fp, "nbhd_size=%d, max_nbrs=%d ", parms->nbhd_size,parms->max_nbrs);
#endif
  if (parms->desired_rms_height > 0.0) {
    fprintf(fp, "desired rms height=%2.3f", parms->desired_rms_height);
  }
  if (first) {
    fprintf(fp, "complete_dist_mat %d\n", parms->complete_dist_mat);
    fprintf(fp, "rms %g\n", parms->rms);
    fprintf(fp, "smooth_averages %d\n", parms->smooth_averages);
    fprintf(fp, "remove_neg %d\n", parms->remove_neg);
    fprintf(fp, "ico_order %d\n", parms->ico_order);
    fprintf(fp, "which_surface %d\n", parms->which_surface);
    fprintf(fp, "target_radius %f\n", parms->target_radius);
    fprintf(fp, "nfields %d\n", parms->nfields);
    fprintf(fp, "scale %lf\n", parms->scale);
    fprintf(fp, "desired_rms_height %lf\n", parms->desired_rms_height);
    fprintf(fp, "momentum %lf\n", parms->momentum);
    fprintf(fp, "nbhd_size %d\n", parms->nbhd_size);
    fprintf(fp, "max_nbrs %d\n", parms->max_nbrs);
    fprintf(fp, "niterations %d\n", parms->niterations);
    fprintf(fp, "nsurfaces %d\n", parms->nsurfaces);
    fprintf(fp, "SURFACES %d\n", (int)(SURFACES));
    fprintf(fp, "flags %d (%x)\n", parms->flags, parms->flags);
    fprintf(fp, "use curv %d\n", (parms->flags & IP_USE_CURVATURE));
    fprintf(fp, "no sulc %d\n", (parms->flags & IP_NO_SULC));
    fprintf(fp, "no rigid align %d\n", (parms->flags & IP_NO_RIGID_ALIGN));
    fprintf(fp, "mris->nsize %d\n", mris->nsize);
    fprintf(fp, "mris->hemisphere %d\n", mris->hemisphere);
    fprintf(fp, "randomSeed %ld\n", getRandomSeed());
    fprintf(fp, "\n");
    first = 0;
  }
  fflush(fp);
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  ------------------------------------------------------*/
static int mrisIntegrationEpoch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int base_averages)
{
  int total_steps, done, steps, n_averages, old_averages;
  const char *snum, *sdenom;
  float ratio, *pdenom, *pnum;

  if (!FZERO(parms->l_corr)) {
    sdenom = "corr";
    pdenom = &parms->l_corr;
  }
  else {
    sdenom = "dist";
    pdenom = &parms->l_dist;
  }

  if (!FZERO(parms->l_nlarea)) {
    snum = "nlarea";
    pnum = &parms->l_nlarea;
  }
  else if (!FZERO(parms->l_area)) {
    snum = "area";
    pnum = &parms->l_area;
  }
  else if (!FZERO(parms->l_parea)) {
    snum = "parea";
    pnum = &parms->l_parea;
  }
  else {
    snum = "spring";
    pnum = &parms->l_spring;
  }

  if (Gdiag & DIAG_SHOW) {
    mrisLogIntegrationParms(stderr, mris, parms);
  }
  if (Gdiag & DIAG_WRITE) {
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  if (!FZERO(*pdenom)) {
    ratio = *pnum / *pdenom;
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "%s/%s = %2.3f\n", snum, sdenom, ratio);
    }
    if (Gdiag & DIAG_WRITE) {
      char fname[STRLEN];
      if (!parms->fp) {
        int req = snprintf(fname, STRLEN, "%s.%s.out", 
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);   
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}

        if (!parms->start_t) {
          INTEGRATION_PARMS_openFp(parms, fname, "w");
        }
        else {
          INTEGRATION_PARMS_openFp(parms, fname, "a");
        }
        if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
      }
      fprintf(parms->fp, "%s/%s = %2.3f\n", snum, sdenom, ratio);
    }
  }

  old_averages = parms->n_averages;
  for (done = total_steps = 0, n_averages = base_averages; !done; n_averages /= 4) {
    parms->n_averages = n_averages;
    steps = MRISintegrate(mris, parms, n_averages);

    if (n_averages > 0 && parms->flags & IP_RETRY_INTEGRATION &&
        ((parms->integration_type == INTEGRATE_LINE_MINIMIZE) || (parms->integration_type == INTEGRATE_LM_SEARCH))) {
      int niter = parms->niterations;
      int integration_type = parms->integration_type;

      fprintf(stdout, "taking momentum steps...\n");
      parms->integration_type = INTEGRATE_MOMENTUM;
      parms->niterations = 10;
      parms->start_t += steps;
      total_steps += steps;
      steps = MRISintegrate(mris, parms, n_averages);
      parms->integration_type = integration_type;
      parms->niterations = niter;
      parms->start_t += steps;
      total_steps += steps;
      steps = MRISintegrate(mris, parms, n_averages);
    }
    parms->start_t += steps;
    total_steps += steps;
    done = n_averages == parms->min_averages;
    if (mris->status == MRIS_SPHERE) {
      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
        MRISprintTessellationStats(mris, stderr);
      }
      parms->scale *= parms->dt_decrease;
      if (parms->scale < 1.0f) {
        parms->scale = 1.0f;
      }
    }
  }
#if 0
  MRIScomputeNormals(mris) ;
  mrisComputeVertexDistances(mris) ;
  MRIScomputeTriangleProperties(mris) ;  /* compute areas and normals */
  mrisOrientSurface(mris) ;
#endif
  parms->n_averages = old_averages; /* hack, but no time to clean up now */
  return (total_steps);
}


static int mrisComputeVariableSmoothnessCoefficients(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  VERTEX *v;
  int vno;
  float vsmooth;

  if (parms->vsmoothness == NULL) {
    return (NO_ERROR);
  }

  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    v->val = parms->l_dist * parms->dist_error[vno] - parms->l_corr * parms->geometry_error[vno];
  }

#define VDELTA 0.005
  MRISaverageVals(mris, 64);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (v->ripflag) {
      continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    vsmooth = parms->vsmoothness[vno];
    if (v->val > 0) {
      vsmooth += VDELTA;
    }
    else {
      vsmooth -= VDELTA;
    }
    vsmooth = MAX(0, vsmooth);
    vsmooth = MIN(1, vsmooth);  // make it in [0 1]
    parms->vsmoothness[vno] = vsmooth;
  }

  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description

  ------------------------------------------------------*/
int MRISintegrate(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int n_averages)
{
  int t, write_iterations, niterations, nsmall, neg;
  double l_dist, l_area, l_spring, sse, old_sse, delta_t, total_small = 0.0;
  double sse_thresh, pct_neg, pct_neg_area, total_vertices, tol;
  /*, scale, last_neg_area */;
  MHT *mht_v_current = NULL;

  if (Gdiag & DIAG_WRITE && parms->fp == NULL) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.%s.out", 
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name); 
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w");
    }
    else {
      INTEGRATION_PARMS_openFp(parms, fname, "a");
    }
    if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  l_spring = parms->l_spring;
  l_dist = parms->l_dist;
  l_area = parms->l_area;
  write_iterations = parms->write_iterations;
  niterations = parms->niterations;
  if (((parms->flags & IPFLAG_QUICK) == 0) && ((parms->flags & IPFLAG_NOSCALE_TOL) == 0)) {
    tol = parms->tol * sqrt(((double)n_averages + 1.0) / 1024.0);
  }
  else {
    tol = parms->tol;
  }
  sse_thresh = tol;
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "integrating with navgs=%d and tol=%2.3e\n", n_averages, tol);
  }

  mrisProjectSurface(mris);

  MRIScomputeMetricProperties(mris);            // this can change XYZ slightly

#if AVERAGE_AREAS
  MRISreadTriangleProperties(mris, mris->fname);
  mrisAverageAreas(mris, n_averages, ORIG_AREAS);
#endif

  parms->starting_sse = sse = old_sse = MRIScomputeSSE(mris, parms);

  delta_t = 0.0;
  niterations += parms->start_t;
  parms->t = parms->start_t;

  if (!parms->start_t) {
    mrisLogStatus(mris, parms, stderr, 0.0f, -1);
    if (Gdiag & DIAG_WRITE) {
      mrisLogStatus(mris, parms, parms->fp, 0.0f, -1);
      if (write_iterations > 0) {
        mrisWriteSnapshot(mris, parms, 0);
      }
    }
  }

  if (Gdiag_no >= 0) {
    fprintf(stdout,
            "v %d curvature = %2.5f, position = (%2.3f,%2.3f,%2.3f)\n",
            Gdiag_no,
            mris->vertices[Gdiag_no].H,
            mris->vertices[Gdiag_no].x,
            mris->vertices[Gdiag_no].y,
            mris->vertices[Gdiag_no].z);
    if (parms->l_thick_normal > 0 || parms->l_thick_parallel > 0 || parms->l_thick_min > 0 || parms->l_thick_spring) {
      float dx, dy, dz, xp, yp, zp, xw, yw, zw;
      VERTEX *v;

      v = &mris->vertices[Gdiag_no];

      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);

      dx = xp - xw;
      dy = yp - yw;
      dz = zp - zw;
      printf("v %d: C=(%2.1f %2.1f %2.1f), W=(%2.1f %2.1f %2.1f), P=(%2.1f %2.1f %2.1f), DX=(%2.2f %2.2f %2.2f)\n",
             Gdiag_no,
             v->x,
             v->y,
             v->z,
             xw,
             yw,
             zw,
             xp,
             yp,
             zp,
             dx,
             dy,
             dz);
    }
  }
  total_small = 0.0;
  nsmall = 0;

  total_vertices = (double)MRISvalidVertices(mris);
  neg = MRIScountNegativeTriangles(mris);
  pct_neg = (double)neg / total_vertices;
  pct_neg_area = (float)mris->neg_area / (float)(mris->total_area + mris->neg_area);

  if (parms->l_thick_min > 0 && (Gdiag & DIAG_WRITE) && (parms->write_iterations > 0)) {
    char fname[STRLEN];
    int vno;
    MRI *mri_vector;

    mri_vector = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);

    for (vno = 0; vno < mris->nvertices; vno++) {
      float xw, yw, zw, xp, yp, zp;
      VERTEX *v;

      v = &mris->vertices[vno];
      if (vno == Gdiag_no) DiagBreak();

      if (v->ripflag) continue;
      MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
      MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);
      MRIsetVoxVal(mri_vector, vno, 0, 0, 0, xp - xw);
      MRIsetVoxVal(mri_vector, vno, 0, 0, 1, yp - yw);
      MRIsetVoxVal(mri_vector, vno, 0, 0, 2, zp - zw);
    }

    int req = snprintf(fname, STRLEN, "%s.%s.%3.3d.mgz",
		       mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms->base_name, parms->t);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    
    printf("writing vectors to %s\n", fname);
    MRIwrite(mri_vector, fname);
    MRIfree(&mri_vector);
  }

  /*  mrisClearMomentum(mris) ;*/
  for (parms->t = t = parms->start_t; t < niterations; t++) {
    if (!FZERO(parms->l_repulse_ratio)) {
      MHTfree(&mht_v_current);
      mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 3.0);
    }
    
    if (!FZERO(parms->l_curv)) {
      MRIScomputeSecondFundamentalForm(mris);
    }

    MRISclearGradient(mris); /* clear old deltas */
    mrisComputeVariableSmoothnessCoefficients(mris, parms);
    mrisComputeDistanceTerm(mris, parms);
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0 && DIAG_VERBOSE_ON) {
      MRI *mri;
      char fname[STRLEN];
      sprintf(fname, "dist_dx%04d.mgz", parms->t);
      mri = MRISwriteIntoVolume(mris, NULL, VERTEX_DX);
      MRIwrite(mri, fname);
      sprintf(fname, "dist_dy%04d.mgz", parms->t);
      MRISwriteIntoVolume(mris, mri, VERTEX_DY);
      MRIwrite(mri, fname);
      sprintf(fname, "dist_dz%04d.mgz", parms->t);
      MRISwriteIntoVolume(mris, mri, VERTEX_DZ);
      MRIwrite(mri, fname);
      MRIfree(&mri);
    }
    mrisComputeAngleAreaTerms(mris, parms);
    mrisComputeCorrelationTerm(mris, parms);
    mrisComputePolarCorrelationTerm(mris, parms);
    /*    mrisComputeSpringTerm(mris, parms->l_spring) ;*/
    /* vectorial registration */
    if (parms->flags & IP_USE_MULTIFRAMES) {
      mrisComputeVectorCorrelationTerm(mris, parms);
      mrisComputePolarVectorCorrelationTerm(mris, parms);
    }

    mrisComputeRepulsiveRatioTerm(mris, parms->l_repulse_ratio, mht_v_current);

    mrisComputeLaplacianTerm(mris, parms->l_lap);
    MRISaverageGradients(mris, n_averages);
    mrisComputeSpringTerm(mris, parms->l_spring);
    mrisComputeThicknessMinimizationTerm(mris, parms->l_thick_min, parms);
    mrisComputeThicknessParallelTerm(mris, parms->l_thick_parallel, parms);
    mrisComputeThicknessNormalTerm(mris, parms->l_thick_normal, parms);
    mrisComputeThicknessSpringTerm(mris, parms->l_thick_spring, parms);
    mrisComputeAshburnerTriangleTerm(mris, parms->l_ashburner_triangle, parms);

    mrisComputeTangentialSpringTerm(mris, parms->l_tspring);
    mrisComputeNonlinearTangentialSpringTerm(mris, parms->l_nltspring, parms->min_dist);
    mrisComputeNonlinearSpringTerm(mris, parms->l_nlspring, parms);
    mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv);
    if (Gdiag & DIAG_WRITE && parms->write_iterations > 0 && DIAG_VERBOSE_ON) {
      MRI *mri;
      char fname[STRLEN];
      sprintf(fname, "dx%04d.mgz", parms->t);
      mri = MRISwriteIntoVolume(mris, NULL, VERTEX_DX);
      MRIwrite(mri, fname);
      sprintf(fname, "dy%04d.mgz", parms->t);
      MRISwriteIntoVolume(mris, mri, VERTEX_DY);
      MRIwrite(mri, fname);
      sprintf(fname, "dz%04d.mgz", parms->t);
      MRISwriteIntoVolume(mris, mri, VERTEX_DZ);
      MRIwrite(mri, fname);
      MRIfree(&mri);
    }

    switch (parms->integration_type) {
      case INTEGRATE_LM_SEARCH:
        delta_t = mrisLineMinimizeSearch(mris, parms);
        break;
      default:
      case INTEGRATE_LINE_MINIMIZE:
        delta_t = mrisLineMinimize(mris, parms);
        break;
      case INTEGRATE_MOMENTUM:
        delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, tol, parms->n_averages);
        break;
      case INTEGRATE_ADAPTIVE:
        delta_t = mrisAdaptiveTimeStep(mris, parms);
        break;
    }

    if (!FZERO(parms->l_thick_min)) {
      int vno;
      VERTEX *v;
      double mn, sq;

      for (mn = 0.0, vno = 0; vno < mris->nvertices; vno++) {
        v = &mris->vertices[vno];
        if (v->ripflag) {
          continue;
        }
        sq = mrisSampleMinimizationEnergy(mris, vno, parms, v->x, v->y, v->z);
        mn += sqrt(sq);
      }
      mn /= (double)MRISvalidVertices(mris);
      printf("mean thickness = %2.4f\n", mn);
      if ((Gdiag & DIAG_WRITE) && (parms->write_iterations > 0) && !((t + 1) % parms->write_iterations)) {
        char fname[STRLEN];
        int vno;
        MRI *mri_vector;

        mri_vector = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, 3);

        for (vno = 0; vno < mris->nvertices; vno++) {
          float xw, yw, zw, xp, yp, zp;
          VERTEX *v;

          v = &mris->vertices[vno];
          if (vno == Gdiag_no) DiagBreak();

          if (v->ripflag) continue;

          MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
          MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);
          MRIsetVoxVal(mri_vector, vno, 0, 0, 0, xp - xw);
          MRIsetVoxVal(mri_vector, vno, 0, 0, 1, yp - yw);
          MRIsetVoxVal(mri_vector, vno, 0, 0, 2, zp - zw);
        }

        int req = snprintf(fname,
			   STRLEN,
			   "%s.%s.%3.3d.vec.mgz",
			   mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh",
			   parms->base_name,
			   parms->t + 1);   
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}

        printf("writing vectors to %s\n", fname);
        MRIwrite(mri_vector, fname);
        MRIfree(&mri_vector);
      }
    }

    if (!FZERO(parms->l_pcorr) && (Gdiag & DIAG_SHOW)) {
      float alpha, beta, gamma;

      alpha = DEGREES(delta_t * mris->alpha);
      beta = DEGREES(delta_t * mris->beta);
      gamma = DEGREES(delta_t * mris->gamma);
      fprintf(stdout, "rotating brain by (%2.1f, %2.1f, %2.1f)\n", alpha, beta, gamma);
    }

    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    if (parms->remove_neg && mris->neg_area > 0) {
      INTEGRATION_PARMS p;
      //      printf("removing overlap with smoothing\n") ;
      p.niterations = 150;
      MRISremoveOverlapWithSmoothing(mris, &p);
    }
    if (Gdiag_no >= 0 && DIAG_VERBOSE_ON)
      fprintf(stdout, "v %d curvature = %2.3f\n", Gdiag_no, mris->vertices[Gdiag_no].H);
    /* only print stuff out if we actually took a step */
    sse = MRIScomputeSSE(mris, parms);
    if (!FZERO(old_sse) && ((old_sse - sse) / (old_sse) < sse_thresh)) {
      if (++nsmall > MAX_SMALL) {
        break;
      }
      if (++total_small > TOTAL_SMALL) {
        break;
      }
    }
    else {
      if (total_small > 0.0) /* if error increases more
                                than 1/4 time quit */
      {
        total_small -= .25;
      }
      nsmall = 0;
    }

    parms->t++;
    mrisLogStatus(mris, parms, stderr, delta_t, old_sse);
    if (Gdiag & DIAG_WRITE) {
      mrisLogStatus(mris, parms, parms->fp, delta_t, -1);
    }

    if ((Gdiag & DIAG_SHOW) && !((t + 1) % 10) && DIAG_VERBOSE_ON) {
      MRISprintTessellationStats(mris, stderr);
    }

    if ((write_iterations > 0) && !((t + 1) % write_iterations) && (Gdiag & DIAG_WRITE) &&
        (FZERO(parms->l_thick_min))) {
      mrisWriteSnapshot(mris, parms, t + 1);
    }
    if (mris->status == MRIS_PLANE && mris->neg_area > 4 * mris->total_area) {
      fprintf(stdout, "flipping flattened patch...\n");
      mrisClearMomentum(mris);
      mrisFlipPatch(mris);
      MRIScomputeMetricProperties(mris);
    }

    if (FZERO(sse)) {
      break;
    }
    if ((parms->integration_type == INTEGRATE_LINE_MINIMIZE) || (parms->integration_type == INTEGRATE_LM_SEARCH) ||
        (parms->integration_type == INTEGRATE_MOMENTUM)) {
      if ((100 * (old_sse - sse) / sse) < tol) {
        break;
      }
    }
    old_sse = sse;
    if (FZERO(delta_t)) /* reached the minimum */
    {
      break;
    }
  }

  if (!FZERO(parms->l_repulse)) {
    MHTfree(&mht_v_current);
  }

  parms->ending_sse = MRIScomputeSSE(mris, parms);
  /*  mrisProjectSurface(mris) ;*/

  return (parms->t - parms->start_t); /* return actual # of steps taken */
}

/*
  Note that at the start of this function, the ORIGINAL_VERTICES must
  contain the surface that has the metric properties to be preserved (e.g.
  smoothwm, white), and the CANONICAL_VERTICES must contain the uniform
  spherical surface (e.g. ?h.sphere).
*/
int MRISregister(MRI_SURFACE *mris,
                 MRI_SP *mrisp_template,
                 INTEGRATION_PARMS *parms,
                 int max_passes,
                 float min_degrees,
                 float max_degrees,
                 int nangles)
{
  float sigma /*, target_sigma, dof*/;
  int i, start_t, sno, ino, msec, min_averages = 0, nsurfaces, using_big_averages = 0;
  MRI_SP *mrisp;
  char fname[STRLEN], base_name[STRLEN], path[STRLEN];
  double base_dt;
  int first = 1;
  INTEGRATION_PARMS saved_parms;

  printf("MRISregister() -------\n");
  printf("max_passes = %d \n", max_passes);
  printf("min_degrees = %f \n", min_degrees);
  printf("max_degrees = %f \n", max_degrees);
  printf("nangles = %d \n", nangles);
  mrisLogIntegrationParms(stdout, mris, parms);
  printf("--------------------\n");

  // ATH: this is an ugly hack to restore the old (and incorrect) computation of the distance
  // term, which had previously used an incorrect value for avg_nbrs
  MRISresetNeighborhoodSize(mris, 3);
  float incorrect_avg_nbrs = mris->avg_nbrs;
  MRISresetNeighborhoodSize(mris, 1);

  INTEGRATION_PARMS_copy(&saved_parms, parms);

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }
  Timer start;
  FileNamePath(mris->fname, path);
  int req = snprintf(base_name, STRLEN, "%s/%s.%s",
		     path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms->base_name);  
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  

  mrisComputeOriginalVertexDistances(mris);
  
  if (parms->nbhd_size > 3) {
    int nbrs[MAX_NBHD_SIZE];

    printf("sampling long-range distances\n");
    memset(nbrs, 0, MAX_NBHD_SIZE * sizeof(nbrs[0]));
    for (i = mris->nsize + 1; i <= parms->nbhd_size; i++) {
      nbrs[i] = parms->max_nbrs;
    }
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
    MRISsampleDistances(mris, nbrs, parms->nbhd_size);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
    MRIScomputeMetricProperties(mris);
  }
  
  base_dt = parms->dt;
  if (Gdiag & DIAG_WRITE) {
    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w");
    }
    else {
      INTEGRATION_PARMS_openFp(parms, fname, "a");
    }
    if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  //  if (Gdiag & DIAG_SHOW)
  mrisLogIntegrationParms(stderr, mris, parms);

  if (parms->flags & IP_NO_SULC) {
    fprintf(stderr, "will not use the sulcal depth map\n");
    fprintf(stderr, "will not rigidly align the surface\n");
    first = 0;
    sno = 2;
  }
  else if (parms->flags & IP_USE_INFLATED) {
    sno = 0;
  }
  else {
    sno = 1;
  }

  if (parms->nsurfaces > 0) {
    nsurfaces = parms->nsurfaces;
  }
  else {
    nsurfaces = SURFACES;
  }

  using_big_averages = ((parms->start_t == 0) && (parms->first_pass_averages > 0));
  for (; sno < nsurfaces; sno++) {
    if (!first && ((parms->flags & IP_USE_CURVATURE) == 0)) {
      break;
    }

    ino = parms->frame_no = sno * IMAGES_PER_SURFACE;
    if (curvature_names[sno]) /* read in precomputed curvature file */
    {
      int req = snprintf(fname, STRLEN, "%s.%s", 
			 mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", curvature_names[sno]); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      printf("%d Reading %s\n", sno, fname);
      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
        ErrorExit(Gerror, "%s: could not read curvature file '%s'\n", "MRISregister", fname);
      MRISnormalizeCurvature(mris, parms->which_norm);
      if (parms->nonmax)
	MRISnonmaxSuppress(mris) ;
      if (parms->trinarize_thresh > 0) {
        MRIStrinarizeCurvature(mris, parms->trinarize_thresh);
      }
    }
    else /* compute curvature of surface */
    {
      int req = snprintf(fname, STRLEN, "%s", mrisurf_surface_names[sno]);
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      MRISsaveVertexPositions(mris, TMP_VERTICES);
      printf("%d Reading %s\n", sno, fname);
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", "MRISregister", fname);

      MRISresetNeighborhoodSize(mris, -1); /* back to max */
      MRIScomputeMetricProperties(mris);
      MRIScomputeSecondFundamentalForm(mris);
      MRISuseMeanCurvature(mris);
      MRISnormalizeCurvature(mris, parms->which_norm);
      if (parms->nonmax) {
	MRISnonmaxSuppress(mris) ;
      }
      MRISresetNeighborhoodSize(mris, 1); /*only use nearest neighbor distances*/

      MRISrestoreVertexPositions(mris, TMP_VERTICES);
      MRIScomputeMetricProperties(mris);
    }
    MRISstoreMeanCurvature(mris);  // store current curv target in H

    // ATH: see above explanation for incorrect_avg_nbrs
    if (getenv("MRIS_REGISTER_NEW_BEHAVIOR") == nullptr) mris->avg_nbrs = incorrect_avg_nbrs;

    if (Gdiag & DIAG_SHOW) {
      if (curvature_names[sno]) {
        printf("reading precomputed curvature from %s\n", fname);
      }
      else {
        printf("calculating curvature of %s surface\n", fname);
      }
    }

    if (Gdiag & DIAG_WRITE) {
      if (curvature_names[sno]) {
        fprintf(parms->fp, "using precomputed curvature from %s\n", fname);
      }
      else {
        fprintf(parms->fp, "calculating curvature of %s surface\n", fname);
      }
    }

    if (sno == 2 && parms->flags & IP_USE_CURVATURE) {
      /* only small adjustments needed after 1st time around */
      parms->tol *= 2.0f;
      parms->l_corr /= 20.0f;  // should be more adaptive - used to be 20
#if 1
      parms->l_spring = .5;  // regularize mesh
#else
      MRISclearOrigDistances(mris);  // replicates old bug
#endif
      if (Gdiag & DIAG_WRITE) {
        mrisLogIntegrationParms(parms->fp, mris, parms);
      }
      //      if (Gdiag & DIAG_SHOW)
      mrisLogIntegrationParms(stderr, mris, parms);
    }
    else if (!first) /* don't do curvature alignment */
    {
      break; /* finished */
    }

    if (sno == 0)  // doing inflated - make it fairly rigid
    {
      /* only small adjustments needed after 1st time around */
      min_averages = parms->min_averages;
      parms->min_averages = 256;
      parms->l_corr /= 3.0f; /* should be more adaptive */
      if (Gdiag & DIAG_WRITE) {
        mrisLogIntegrationParms(parms->fp, mris, parms);
      }
      //      if (Gdiag & DIAG_SHOW)
      mrisLogIntegrationParms(stderr, mris, parms);
    }

    for (i = 0; i < nsigmas; i++) /* for each spatial scale (blurring) */
    {
      parms->sigma = sigma = sigmas[i];
      parms->dt = base_dt;
      if (Gdiag & DIAG_SHOW) fprintf(stdout, "\nblurring surfaces with sigma=%2.2f...\n", sigma);
      if (Gdiag & DIAG_WRITE) fprintf(parms->fp, "\ncorrelating surfaces with with sigma=%2.2f\n", sigma);
      if ((Gdiag & DIAG_WRITE) && !i && (!parms->start_t || sno <= 1)) {
        fprintf(parms->fp,
                "writing target curvature, i=%d, start_t=%d, "
                "sno=%d, v0=(%2.1f, %2.1f, %2.1f)\n",
                i,
                parms->start_t,
                sno,
                mris->vertices[0].cx,
                mris->vertices[0].cy,
                mris->vertices[0].cz);
        fflush(parms->fp);
        MRISsaveVertexPositions(mris, TMP_VERTICES);
        MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
        MRIScomputeMetricProperties(mris);
        MRISfromParameterization(mrisp_template, mris, ino);
        MRISnormalizeCurvature(mris, parms->which_norm);
        int req = snprintf(fname,
			   STRLEN,
			   "%s/%s.%s.target%d",
			   path,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->base_name,
			   sno);    
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
	
        MRISwriteCurvature(mris, fname);
        MRISrestoreVertexPositions(mris, TMP_VERTICES);
        MRIScomputeMetricProperties(mris);
      }
      MRISuseMeanCurvature(mris);  // restore current target
      mrisp = MRIStoParameterization(mris, NULL, 1, 0);
#if 1
      parms->mrisp = MRISPblur(mrisp, NULL, sigma, 0);
      parms->mrisp_template = MRISPblur(mrisp_template, NULL, sigma, ino);
      MRISPblur(mrisp_template, parms->mrisp_template, sigma, ino + 1); /* variances */
#else
      dof = *IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2);
      if (dof < 1) {
        dof = 1;
      }
      target_sigma = sigma / dof;
      target_sigma = sigma;
      printf("template has %2.0f dofs - setting target sigma to %2.2f\n", dof, target_sigma);
      parms->mrisp = MRISPiterative_blur(mris, mrisp, NULL, sigma, 0);
      parms->mrisp_template = MRISPiterative_blur(mris, mrisp_template, NULL, target_sigma, ino);
      MRISPiterative_blur(mris, parms->mrisp_template, parms->mrisp_template, target_sigma, ino + 1); /* variances */
#endif
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "done.\n");
      }
      /* normalize curvature intensities for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino);
      MRISnormalizeCurvature(mris, parms->which_norm);
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino);

#if 0
      /* normalize variances for both source and target */
      MRISfromParameterization(parms->mrisp_template, mris, ino+1);
      MRISnormalizeCurvature(mris, parms->which_norm) ;
      MRIStoParameterization(mris, parms->mrisp_template, 1, ino+1) ;
#endif

      if (Gdiag & DIAG_WRITE) {
#if 1
        MRISsaveVertexPositions(mris, TMP_VERTICES);
        MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);
        MRIScomputeMetricProperties(mris);
#endif
        MRISfromParameterization(mrisp_template, mris, ino);
        int req = snprintf(fname,
			   STRLEN,
			   "%s/%s.%s.sno%d_target_blur%2.2f",
			   path,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->base_name,
			   sno,
			   sigma);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
	
        MRISwriteCurvature(mris, fname);
#if 1
        MRISrestoreVertexPositions(mris, TMP_VERTICES);
        MRIScomputeMetricProperties(mris);
#endif
      }

      MRISfromParameterization(parms->mrisp, mris, 0);
      MRISnormalizeCurvature(mris, parms->which_norm);
      MRIStoParameterization(mris, parms->mrisp, 1, 0);
      MRISPfree(&mrisp);

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        MRISPwrite(parms->mrisp, "mrisp_blur.hipl");
        MRISPwrite(parms->mrisp_template, "mrisp_template_blur.hipl");
      }
      mris->vp = (void *)parms->mrisp; /* hack to get it
                                to projectSurface */

      if (Gdiag & DIAG_WRITE) {
        int req = snprintf(fname,
			   STRLEN,
			   "%s/%s.%s.sno%d_blur%2.2f",
			   path,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->base_name,
			   sno,
			   sigma);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}

        MRISwriteCurvature(mris, fname);
      }
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        int req = snprintf(fname,
			   STRLEN,
			   "%s/%s.%s.%4.4dblur%2.2f",
			   path,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->base_name,
			   parms->start_t,
			   sigma);
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
        if (Gdiag & DIAG_SHOW) {
          fprintf(stdout, "writing curvature file %s...", fname);
        }
        MRISwriteCurvature(mris, fname);
        if (Gdiag & DIAG_SHOW) {
          fprintf(stdout, "done.\n");
        }
        req = snprintf(fname, STRLEN, "target.%s.%4.4d.hipl",
		       parms->base_name, parms->start_t);   
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
       
        if (Gdiag & DIAG_SHOW) {
          fprintf(stdout, "writing parameterization file %s...", fname);
        }
        MRISPwrite(parms->mrisp_template, fname);
        if (Gdiag & DIAG_SHOW) {
          fprintf(stdout, "done.\n");
        }
      }
      if (first && i == 0) /* only do rigid alignment first time through */
      {
        if (sno >= 1)  // do global rigid for both inflated and sulc
        {
          first = 0;
        }
        if ((parms->flags & IP_NO_RIGID_ALIGN) == 0) {
          if (Gdiag & DIAG_SHOW) {
            fprintf(stdout, "finding optimal rigid alignment\n");
          }
          if (Gdiag & DIAG_WRITE) {
            fprintf(parms->fp, "finding optimal rigid alignment\n");
          }
          MRISrigidBodyAlignGlobal(mris, parms, min_degrees, max_degrees, nangles);
          /* MRISrigidBodyAlignGlobal(mris, parms, 0.5f, 32.0f, 8) ;*/
          if (Gdiag & DIAG_WRITE && parms->write_iterations != 0) {
            char fname[STRLEN], path[STRLEN];
            FileNamePath(mris->fname, path);
            int req = snprintf(fname,
			       STRLEN,
			       "%s/%s.%s.sno%d.rotated",
			       path,
			       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			       parms->base_name,
			       sno);
	    if( req >= STRLEN ) {
	      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	    }
	    
            printf("writing rigid aligned surface to %s\n", fname);
            MRISwrite(mris, fname);
          }
          MRISsaveVertexPositions(mris, TMP2_VERTICES);
          if (parms->niterations == 0)  // only rigid
          {
            break;
          }
        }
      }

      mrisClearMomentum(mris);

      if (using_big_averages) {
        float sigma = 4.0;
        MRISsetRegistrationSigmas(&sigma, 1);
        mrisIntegrationEpoch(mris, parms, parms->first_pass_averages);
        MRISsetRegistrationSigmas(NULL, 0);
        using_big_averages = 0;
      }
      mrisIntegrationEpoch(mris, parms, parms->n_averages);
    }
    if (parms->niterations == 0)  // only rigid
    {
      break;
    }
    if (sno == 0)  // doing inflated - was fairly rigid - restore orig values
    {
      /* only small adjustments needed after 1st time around */
      parms->min_averages = min_averages;
      parms->l_corr *= 3.0f; /* should be more adaptive */
      if (Gdiag & DIAG_WRITE) {
        mrisLogIntegrationParms(parms->fp, mris, parms);
      }
      if (Gdiag & DIAG_SHOW) {
        mrisLogIntegrationParms(stderr, mris, parms);
      }
    }
  }

#if 1
  if (mris->neg_area > 0) {
    parms->tol /= 10; /* remove everything possible pretty much */
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "\nRemoving remaining folds...\n");
    }
    if (Gdiag & DIAG_WRITE) {
      fprintf(parms->fp, "removing remaining folds...\n");
    }
    parms->l_nlarea *= 100;
    parms->n_averages = 64;  // don't let averaging effect too much of surface
    parms->l_parea /= 100;
    parms->l_spring /= 100;
    parms->l_corr /= 100;
    parms->l_dist /= 100;
    mrisIntegrationEpoch(mris, parms, parms->n_averages);
  }
#else
  parms->l_nlarea = 1;
  parms->l_corr /= 10.0;
  parms->l_area = parms->l_parea = parms->l_spring = 0.0;
  mrisRemoveNegativeArea(mris, parms, parms->n_averages, MAX_NEG_AREA_PCT, 3);
#endif
#if 0
  MRISPfree(&parms->mrisp) ;
  MRISPfree(&parms->mrisp_template) ;
#endif
  msec = start.milliseconds();
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "registration took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->fp, "registration took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
    INTEGRATION_PARMS_closeFp(parms);
  }

  start_t = parms->start_t;
  INTEGRATION_PARMS_copy(parms, &saved_parms);
  parms->start_t = start_t;
  printf("MRISregister() return, current seed %ld\n", getRandomSeed());
  fflush(stdout);
  return (NO_ERROR);
}


int MRISvectorRegister(MRI_SURFACE *mris,
                       MRI_SP *mrisp_template,
                       INTEGRATION_PARMS *parms,
                       int max_passes,
                       float min_degrees,
                       float max_degrees,
                       int nangles)
{
  float sigma;
  int i, /*steps,*/ msec;
  MRI_SP *mrisp;
  VERTEX *v;
  char fname[STRLEN], base_name[STRLEN], path[STRLEN];
  double base_dt;
  static int first = 1;
  int n, fno, ncorrs;
  int *frames, nframes, nf, *indices;
  float l_corr;
  VALS_VP *vp;

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }
  Timer start;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES);
  FileNamePath(mris->fname, path);
  int req = snprintf(base_name, STRLEN, "%s/%s.%s",
		     path, mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh", parms->base_name);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }


  base_dt = parms->dt;
  if (Gdiag & DIAG_WRITE) {
    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);  
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }

    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w");
      if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    }
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  if (Gdiag & DIAG_SHOW) {
    mrisLogIntegrationParms(stderr, mris, parms);
  }

  /*    for ( nframes = n = 0 ; n < parms->nfields ; n++){ */
  /*            parms->fields[n].l_corr=0.0f; */
  /*            parms->fields[n].l_pcorr=0.0f; */
  /*    } */
  /*    parms->l_corrs[5]=1.0f; // only one structure at a time */

  /* excluding frames with zero correlation coefficients */
  for (nframes = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    nframes++;
  }
  if (!nframes) {
    return NO_ERROR;
  }

  fprintf(stderr, "MRISvectorRegister will use %d fields\n", nframes);

  indices = (int *)malloc(nframes * sizeof(int));
  frames = (int *)malloc(2 * nframes * sizeof(int));
  for (nf = n = 0; n < parms->nfields; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue;
    }
    fno = parms->fields[n].frame * IMAGES_PER_SURFACE;
    frames[nf] = fno;               /* mean */
    frames[nf + nframes] = fno + 1; /* variance */
    indices[nf] = n;
    nf++;
  }

  ncorrs = parms->nfields;
  /* allocate the VALS_VP structure */
  for (n = 0; n < mris->nvertices; n++) {
    v = &mris->vertices[n];
    vp = (VALS_VP *)calloc(1, sizeof(VALS_VP));
    vp->nvals = ncorrs;
    vp->orig_vals = (float *)calloc(ncorrs, sizeof(float)); /* before blurring */
    vp->vals = (float *)malloc(ncorrs * sizeof(float));     /* values used by
                                                                       MRISintegrate */
    v->vp = (void *)vp;
  }

  /* load the fields into vertex->vp */
  for (n = 0; n < ncorrs; n++) {
    l_corr = parms->fields[n].l_corr + parms->fields[n].l_pcorr;
    if (FZERO(l_corr)) {
      continue; /* don't load useless fields */
    }

    fprintf(stderr,
            "  -loading field %d with correlation coefficients "
            "( %2.1f , %2.1f )...\n",
            n,
            parms->fields[n].l_corr,
            parms->fields[n].l_pcorr);
    if (parms->fields[n].name != NULL) {
      char path[STRLEN];
      FileNamePath(mris->fname, path);
      if (parms->overlay_dir == NULL) {
        int req = snprintf(fname,
			   STRLEN,
			   "%s/../label/%s.%s",
			   path,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->fields[n].name);  
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}

      }
      else {
        int req = snprintf(fname,
			   STRLEN,
			   "%s/../%s/%s.%s",
			   path,
			   parms->overlay_dir,
			   mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			   parms->fields[n].name);    
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
	
      }
      printf("reading overlay file %s...\n", fname);
      if (MRISreadValues(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_BADPARM, "%s: could not read overlay file %s", Progname, fname);
      if (mris->ct) {
        MRISripMedialWall(mris);
      }
      MRIScopyValuesToCurvature(mris);
    }
    else if (ReturnFieldName(parms->fields[n].field)) {
      /* read in precomputed curvature file */
      int req = snprintf(fname, STRLEN, "%s.%s", 
			 mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", ReturnFieldName(parms->fields[n].field));
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR) {
        fprintf(stderr, "%s: could not read curvature file '%s'\n", "MRISvectorRegister", fname);
        fprintf(stderr, "setting up correlation coefficient to zero\n");
        parms->fields[n].l_corr = parms->fields[n].l_pcorr = 0.0;
        continue;
      }
    }
    else {
      /* compute curvature of surface */
      int req = snprintf(fname, STRLEN, "%s", mrisurf_surface_names[parms->fields[n].field]);   
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }

      /*                        if(parms->fields[n].field==0) */
      /*                                sprintf(fname, "inflated") ; */
      /*                        else */
      /*                                sprintf(fname, "smoothwm") ; */
      MRISsaveVertexPositions(mris, TMP_VERTICES);
      if (MRISreadVertexPositions(mris, fname) != NO_ERROR) {
        ErrorPrintf(ERROR_NOFILE, "%s: could not read surface file %s", "MRISvectorRegister", fname);
        fprintf(stderr, "setting up correlation coefficient to zero\n");
        parms->fields[n].l_corr = parms->fields[n].l_pcorr = 0.0;
        continue;
      }
      MRISsetNeighborhoodSizeAndDist(mris, -1); /* back to max */
      MRIScomputeMetricProperties(mris);
      MRIScomputeSecondFundamentalForm(mris);
      MRISuseMeanCurvature(mris);
      MRISresetNeighborhoodSize(mris, 1); /*only use nearest neighbor distances*/
      MRISrestoreVertexPositions(mris, TMP_VERTICES);
    }
    MRISnormalizeField(mris, parms->fields[n].type, parms->fields[n].which_norm);
    MRISsetCurvaturesToOrigValues(mris, n);
  }
  MRISaverageCurvatures(mris, parms->fields[n].navgs);

  /* multiscale registration */
  parms->mrisp_template = MRISPclone(mrisp_template);
  parms->mrisp = MRISPclone(mrisp_template);
  mris->vp = (void *)parms->mrisp; /* hack to get it to projectSurface */

#define DEBUG_NO_BLURRING 0

  for (i = 0; i < nsigmas; i++) /* for each spatial scale (blurring) */
  {
#if (DEBUG_NO_BLURRING)
    i = nsigmas - 1;
    first = 0;
#endif

    if (Gdiag & DIAG_WRITE && !i && !parms->start_t) {
      MRISfromParameterizations(parms->mrisp_template, mris, frames, indices, nframes);
      //      MRISnormalizeCurvature(mris, parms->fields[0].which_norm) ;
      int req = snprintf(fname, STRLEN, "%s/%s.%s.target",
			 path, mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name); 
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "writing curvature file %s...\n", fname);
      }
      MRISwriteCurvature(mris, fname);
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "done.\n");
      }
    }
    if ((parms->flags & IP_NO_RIGID_ALIGN)) /* no rigid alignment */
    {
      first = 0;
    }

    parms->sigma = sigma = sigmas[i];
    parms->dt = base_dt;
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "\nblurring surfaces with sigma=%2.2f...\n", sigma);
    }
    if (Gdiag & DIAG_WRITE) fprintf(parms->fp, "\ncorrelating surfaces with with sigma=%2.2f\n", sigma);

    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
      int req = snprintf(fname,
			 STRLEN,
			 "%s/%s.%s.%4.4dtarget%2.2f",
			 path,
			 mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh",
			 parms->base_name,
			 parms->start_t,
			 sigma);
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "writing curvature file %s...", fname);
      }
      MRISwriteCurvature(mris, fname);
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "done.\n");
      }
    }

#if (!DEBUG_NO_BLURRING)
    /* blurring only the used frames */
    fprintf(stderr, "blurring %d target fields (means and variances)...\n", nframes);
    if (*IMAGEFseq_pix(mrisp_template->Ip, 0, 0, 2) <= 1.0)                           /* 1st time */
      MRISPblurFrames(mrisp_template, parms->mrisp_template, sigma, frames, nframes); /* means only */
    else
      MRISPblurFrames(mrisp_template, parms->mrisp_template, sigma, frames, 2 * nframes); /* means and variances */
#else
    parms->mrisp_template = MRISPclone(mrisp_template);
#endif

    /* normalize mean (only) intensities for target */
    MRISfromParameterizations(parms->mrisp_template, mris, frames, indices, nframes);
    if (Gdiag & DIAG_WRITE) {
      sprintf(fname, "%s.target_%d_%d", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", n, (int)(2 * sigma));
      fprintf(stderr, "writing target curvature for field %d\n", n);
      MRISwriteCurvature(mris, fname);
    }
    for (n = 0; n < nframes; n++) {
      fprintf(stderr,
              "normalized target field %d (frame = %d(%d) - field = %d)...\n",
              indices[n],
              parms->fields[indices[n]].frame,
              frames[n],
              parms->fields[indices[n]].field);
      MRISsetValuesToCurvatures(mris, indices[n]);
      MRISnormalizeField(mris, parms->fields[indices[n]].type, parms->fields[indices[n]].which_norm);
      MRISsetCurvaturesToValues(mris, indices[n]);

      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
        sprintf(fname, "%s.target_%d_%d", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", n, (int)(2 * sigma));
        fprintf(stderr, "writing target curvature for field %d\n", n);
        MRISwriteCurvature(mris, fname);
      }
    }
    MRIStoParameterizations(mris, parms->mrisp_template, 1, frames, indices, nframes);

    for (n = 0; n < nframes; n++) {
      MRISsetOrigValuesToValues(mris, indices[n]);
    }

    mrisp = MRISPclone(mrisp_template);
    MRIStoParameterizations(mris, mrisp, 1, frames, indices, nframes);
#if (!DEBUG_NO_BLURRING)
    /* blur source intensities for frame #n */
    fprintf(stderr, "blurring %d source field (means only)...\n", nframes);
    MRISPblurFrames(mrisp, parms->mrisp, sigma, frames, nframes); /* only mean fields */
#else
    parms->mrisp = MRISPclone(mrisp);
#endif
    MRISPfree(&mrisp);

    /* normalize mean intensities for source */
    MRISfromParameterizations(parms->mrisp, mris, frames, indices, nframes);
    for (n = 0; n < nframes; n++) {
      fprintf(stderr,
              "normalized source field %d (frame = %d(%d) - field = %d)...\n",
              indices[n],
              parms->fields[indices[n]].frame,
              frames[n],
              parms->fields[indices[n]].field);
      MRISsetValuesToCurvatures(mris, indices[n]);
      MRISnormalizeField(mris, parms->fields[indices[n]].type, parms->fields[indices[n]].which_norm);
      MRISsetCurvaturesToValues(mris, indices[n]);
      if (Gdiag & DIAG_WRITE) {
        sprintf(fname, "%s.source_%d_%d", mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", n, (int)(2 * sigma));
        fprintf(stderr, "writing source curvature for field %d\n", n);
        MRISwriteCurvature(mris, fname);
      }
    }
    /*MRIStoParameterizations(mris, parms->mrisp,
      1, frames,indices,nframes) ; */

    /* use the frame #indices[0] (sulc) to write out snapshots */
    MRISsetValuesToCurvatures(mris, indices[0]);
    if (mris->ct) {
      MRISripMedialWall(mris);
    }
    MRIStoParameterization(mris, parms->mrisp, 1, 0);

    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "done.\n");
    }

#if 0
    if (!first)
    {
      // TO BE CHECKED XXX
      if (pdone)
      {
        /* only once */
        /* only small adjustments needed after 1st time around */
        parms->tol *= 2.0f ;

        for ( n = 0 ; n < nframes ; n++)
        {
          parms->fields[indices[n]].l_corr /= 20.0f; /* should be more
                                                        adaptive */
          parms->fields[indices[n]].l_pcorr /= 20.0f; /* should be more
                                                         adaptive */
        }
      }
    }
#endif

    if (first) {
      /* do rigid registration the first time only */
      first = 0;
      if ((parms->flags & IP_NO_RIGID_ALIGN) == 0) {
        if (Gdiag & DIAG_SHOW) {
          fprintf(stdout, "finding optimal rigid alignment\n");
        }
        if (Gdiag & DIAG_WRITE) {
          fprintf(parms->fp, "finding optimal rigid alignment\n");
        }

        MRISrigidBodyAlignVectorGlobal(mris, parms, min_degrees, max_degrees, nangles);
        MRISsaveVertexPositions(mris, TMP2_VERTICES);
        if (Gdiag & DIAG_WRITE && parms->write_iterations != 0) {
          MRISwrite(mris, "rotated");
        }
      }
    }

    mrisClearMomentum(mris);

    mrisIntegrationEpoch(mris, parms, parms->n_averages);
  }

  parms->tol /= 10; /* remove everything possible pretty much */
  if (Gdiag & DIAG_SHOW) {
    fprintf(stdout, "\nRemoving remaining folds...\n");
  }
  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->fp, "\nRemoving remaining folds...\n");
  }
  parms->l_nlarea *= 5;
  mrisIntegrationEpoch(mris, parms, parms->n_averages);

#if 0
  /* free everything */
  MRISPfree(&parms->mrisp) ;
  MRISPfree(&parms->mrisp_template) ;

  /* free the VALS_VP structure */
  for ( n = 0; n < mris->nvertices ; n++)
  {
    v=&mris->vertices[n];
    vp=(VALS_VP*)v->vp;
    free(vp->orig_vals);
    free(vp->vals);
    free(vp);
    v->vp=NULL;
  }
#endif

  free(frames);
  free(indices);

  msec = start.milliseconds();
  if (Gdiag & DIAG_SHOW) fprintf(stdout, "registration took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
  if (Gdiag & DIAG_WRITE)
    fprintf(parms->fp, "registration took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Fit a quadratic form to the error function in the gradient
  direction and use it to predict the location of the minimum.
  Pick the dt which minimizes the error function among the
  sampled points, including the predicted one.
  ------------------------------------------------------*/
#define MAX_ENTRIES 100
static double mrisLineMinimize(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  FILE *fp = NULL;
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "%s%4.4d.dat", FileName(parms->base_name), parms->t + 1);
    fp = fopen(fname, "w");
  }

  /* compute the magnitude of the gradient, and the max delta */
  double max_delta = 0.0, sum_squaredDelta = 0.0, sum_delta = 0.0;
  int n = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }

    float dx = vertex->dx;
    float dy = vertex->dy;
    float dz = vertex->dz;
    double deltaSquared = dx * dx + dy * dy + dz * dz;
    double delta        = sqrt(deltaSquared);
    
    sum_squaredDelta += deltaSquared;
    sum_delta        += delta;
    n++;

    if (max_delta < delta) {
      max_delta = delta;
    }
  }
  
  double const grad = sqrt(sum_squaredDelta), mean_delta = sum_delta/((double)(n));

  if (FZERO(max_delta)) {
    return (0.0); /* at a local minimum */
  }


  /* limit the size of the largest time step */
  double max_dt;
  switch (parms->projection) {
    case PROJECT_SPHERE:
    case PROJECT_ELLIPSOID:
      max_dt = MAX_MM / mean_delta;
      break;
    case NO_PROJECTION:
      max_dt = 100.0 * MAX_MM / max_delta;
      break;
    default:
    case PROJECT_PLANE:
      max_dt = MAX_PLANE_MM / max_delta;
      break;
  }
  double const min_dt = MIN_MM / mean_delta;

  double dt_in[MAX_ENTRIES], sse_out[MAX_ENTRIES];
  int N = 0;
  {
    MRIScomputeSSE_asThoughGradientApplied_ctx sseCtx;

    double const starting_sse = MRIScomputeSSE(mris, parms);

    /* write out some data on supposed quadratic form */
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {

      char fname[STRLEN];
      sprintf(fname, "nn%s%4.4d.dat", FileName(parms->base_name), parms->t + 1);
      FILE *fp2 = fopen(fname, "w");

      double delta = max_dt / 100.0;
      double delta_t;
      for (delta_t = delta; delta_t <= max_dt; delta_t += delta) {
        double predicted_sse = starting_sse - grad * delta_t;
        
        double sse = MRIScomputeSSE_asThoughGradientApplied(mris, delta_t, parms, sseCtx);
        
        fprintf(fp2, "%f  %f  %f\n", delta_t, sse, predicted_sse);
        fflush(fp2);
      }

      fclose(fp2);
    }


    double min_sse   = starting_sse;
    double min_delta = 0.0f; /* to get rid of compiler warning */

    /* pick starting step size */
    double delta_t;
    for (delta_t = min_dt; delta_t < max_dt; delta_t *= 10.0) {
      double sse = MRIScomputeSSE_asThoughGradientApplied(mris, delta_t, parms, sseCtx);
      
      if (sse <= min_sse) /* new minimum found */
      {
        min_sse   = sse;
        min_delta = delta_t;
      }

    }

    if (FZERO(min_delta)) /* dt=0 is min starting point, look mag smaller */
    {
      delta_t = min_dt / 10.0; /* start at smallest step */

      double sse = MRIScomputeSSE_asThoughGradientApplied(mris, delta_t, parms, sseCtx);
    
      min_sse = sse;
      min_delta = delta_t;
    }

    delta_t = min_delta;

    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      fprintf(stdout,
              "grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
              "starting dt=%2.3f, min_dt=%2.3f\n",
              (float)grad,
              (float)max_delta,
              mean_delta,
              (float)max_dt,
              (float)delta_t,
              min_dt);

    /* fit a quadratic form to it, and predict location of minimum */
    /* bracket the minimum by sampling on either side */
    double const dt0 = min_delta - (min_delta / 2);
    double const dt2 = min_delta + (min_delta / 2);
  
    double sse0 = MRIScomputeSSE_asThoughGradientApplied(mris, dt0, parms, sseCtx);
    double sse2 = MRIScomputeSSE_asThoughGradientApplied(mris, dt2, parms, sseCtx);

    /* now fit a quadratic form to these values */

    sse_out[0] = sse0;          dt_in[0] = dt0;
    sse_out[1] = min_sse;       dt_in[1] = min_delta;
    sse_out[2] = sse2;          dt_in[2] = dt2;

    N = 3;    // 3 entries in sse_out and dt_in so far
    cheapAssert(N <= MAX_ENTRIES);
    
    MATRIX* mX = MatrixAlloc(N, 3, MATRIX_REAL);
    VECTOR* vY = VectorAlloc(N, MATRIX_REAL);

    for (int i = 1; i <= N; i++) {
      *MATRIX_RELT(mX, i, 1) = dt_in[i - 1] * dt_in[i - 1];
      *MATRIX_RELT(mX, i, 2) = 2 * dt_in[i - 1];
      *MATRIX_RELT(mX, i, 3) = 1.0f;

      VECTOR_ELT(vY, i) = sse_out[i - 1];
    }

    MATRIX* m_xT      = MatrixTranspose(mX,       NULL);
    MATRIX* m_xTx     = MatrixMultiply (m_xT, mX, NULL);
    MATRIX* m_xTx_inv = MatrixInverse  (m_xTx,    NULL);

    if (!m_xTx_inv) {
      fprintf(stderr, "singular matrix in quadratic form\n");
    } else {

      MATRIX* m_xTy = MatrixMultiply(m_xT, vY, NULL);
      MATRIX* mP    = MatrixMultiply(m_xTx_inv, m_xTy, NULL);

      double a = RVECTOR_ELT(mP, 1);
      double b = RVECTOR_ELT(mP, 2);
      double c = RVECTOR_ELT(mP, 3);

      if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
        fprintf(stdout, "(a,b,c) = (%2.3f, %2.3f, %2.3f), predicted min at %2.3f\n", a, b, c, -b / a);
      if (!std::isfinite(a)) {
        DiagBreak();
      }

      MatrixFree(&mP);
      MatrixFree(&m_xTy);

      dt_in  [N] = 0;
      sse_out[N] = starting_sse;
      N++;
      cheapAssert(N <= MAX_ENTRIES);

      if (std::isfinite(a) && !FZERO(a)) {
        float new_min_delta = -b / a;

        if (new_min_delta < 10.0f * min_delta && new_min_delta > min_delta / 10.0f) {
          double sse = MRIScomputeSSE_asThoughGradientApplied(mris, new_min_delta, parms, sseCtx);
	  
          dt_in  [N] = new_min_delta;
          sse_out[N] = sse;
          N++;
          cheapAssert(N <= MAX_ENTRIES);
        }
      }
    }

    if (m_xTx_inv) MatrixFree(&m_xTx_inv);
    MatrixFree(&m_xTx);
    MatrixFree(&m_xT);
    VectorFree(&vY);
    MatrixFree(&mX);
  }


  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "sses: %2.2f  ", sse_out[0]);
  }

  int    min_i = 0;
  double min_sse = sse_out[min_i];
  for (int i = 1; i < N; i++) {
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
      fprintf(stdout, "%2.2f  ", sse_out[i]);
    }
    if (min_sse > sse_out[i]) {
      min_sse = sse_out[i];
      min_i = i;
    }
  }

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON) {
    fprintf(stdout, "min %d (%2.3f)\n", min_i, dt_in[min_i]);
  }

  if (mris->status == MRIS_PLANE)  // remove global translation component
  {
    double dx = 0.0, dy = 0.0, dz = 0.0;
    int nv = 0;

    int vno;
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX const * const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      dx += v->dx;
      dy += v->dy;
      dz += v->dz;
      nv++;
    }
    
    if (nv == 0) nv = 1;
    
    dx /= nv;
    dy /= nv;
    dz /= nv;
    
    for (vno = 0; vno < mris->nvertices; vno++) {
      VERTEX* const v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      v->dx -= dx;
      v->dy -= dy;
      v->dz -= dz;
    }
  }

  if (mris->status == MRIS_SPHERICAL_PATCH && parms->flags & IPFLAG_PRESERVE_TOPOLOGY_CONVEXHULL) {
    mrisApplyTopologyPreservingGradient(mris, dt_in[min_i], 0);
  }
  else if (parms->flags & IPFLAG_PRESERVE_SPHERICAL_POSITIVE_AREA) {
    mrisApplyGradientPositiveAreaPreserving(mris, dt_in[min_i]);
  }
  else if (parms->flags & IPFLAG_MAXIMIZE_SPHERICAL_POSITIVE_AREA) {
    mrisApplyGradientPositiveAreaMaximizing(mris, dt_in[min_i]);
  }
  else {
    MRISapplyGradient(mris, dt_in[min_i]);
  }
  
  return (dt_in[min_i]);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  Use a binary search in the gradient direction to find the
  location of the minimum.
  ------------------------------------------------------*/
static double mrisLineMinimizeSearch(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  FILE *fp = NULL;
  if ((Gdiag & DIAG_WRITE) && DIAG_VERBOSE_ON) {
    char fname[STRLEN];
    sprintf(fname, "%s%4.4d.dat", FileName(parms->base_name), parms->t + 1);
    fp = fopen(fname, "w");
  }

  /* compute the magnitude of the gradient, and the max delta */
  double max_delta = 0.0, sum_deltaSquared = 0.0f, sum_delta = 0.0;
  int n = 0;
  for (int vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * const vertex = &mris->vertices[vno];
    if (vertex->ripflag) {
      continue;
    }
    
    float dx = vertex->dx;
    float dy = vertex->dy;
    float dz = vertex->dz;
    double deltaSquared = dx * dx + dy * dy + dz * dz;
    double delta        = sqrt(deltaSquared);
    
    sum_deltaSquared += deltaSquared;
    sum_delta        += delta;
    if (!FZERO(delta)) {
      n++;
    }
    
    if (max_delta < delta) {
      max_delta = delta;
    }
  }
  
  double const mean_delta = sum_delta/float(n);
  double const grad       = sqrt(sum_deltaSquared);

  if (FZERO(max_delta)) {
    return (0.0); /* at a local minimum */
  }

  /* limit the size of the largest time step */
  double max_dt;
  switch (parms->projection) {
    case PROJECT_SPHERE:
    case PROJECT_ELLIPSOID:
      max_dt = MAX_MM / mean_delta;
      break;
    case NO_PROJECTION:
      max_dt = 100.0 * MAX_MM / max_delta;
      break;
    default:
    case PROJECT_PLANE:
      max_dt = MAX_PLANE_MM / max_delta;
      break;
  }
  double const min_dt = MIN_MM / mean_delta;

  double const starting_sse = MRIScomputeSSE(mris, parms);
  double min_sse = starting_sse;

  /* pick starting step size */
  double min_delta = 0.0f; /* to get rid of compiler warning */
  for (double delta_t = min_dt; delta_t < max_dt; delta_t *= 10.0) {

    MRISapplyGradient(mris, delta_t);
    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    double sse = MRIScomputeSSE(mris, parms);
    MRISrestoreOldPositions(mris);

    if (sse <= min_sse) /* new minimum found */
    {
      min_sse = sse;
      min_delta = delta_t;
    }
  }

  if (FZERO(min_delta)) /* dt=0 is min starting point, look mag smaller */
  {
    min_delta = min_dt / 10.0; /* start at smallest step */

    MRISapplyGradient(mris, min_delta);
    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    double sse = MRIScomputeSSE(mris, parms);
    MRISrestoreOldPositions(mris);

    min_sse = sse; 
  }


  /* now search for minimum in gradient direction */
  double delta_t = min_delta;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    fprintf(stdout,
            "grad=%2.3f, max_del=%2.3f, mean=%2.3f, max_dt=%2.1f, "
            "starting dt=%2.3f, min_dt=%2.3f\n",
            (float)grad,
            (float)max_delta,
            mean_delta,
            (float)max_dt,
            (float)delta_t,
            min_dt);

  int increasing = 1;
  double total_delta = 0.0;
  
  min_sse = starting_sse;
  bool done = false;
  while (!done) {
  
    MRISapplyGradient(mris, delta_t);
    mrisProjectSurface(mris);
    MRIScomputeMetricProperties(mris);
    double sse = MRIScomputeSSE(mris, parms);

    if (sse <= min_sse) /* new minimum found */
    {
      if ((parms->projection == PROJECT_ELLIPSOID) && (total_delta + delta_t > max_dt)) {
        increasing = 0; /* limit size of largest time step */
      }
      min_sse = sse;
      total_delta += delta_t; /* keep track of total time step */
    }
    else /* error increased - undo it and decrease time step */
    {
      if (increasing) /* went too far - reverse delta_t change */
      {
        increasing = 0;
      }

      MRISrestoreOldPositions(mris);
      mrisProjectSurface(mris);
      MRIScomputeMetricProperties(mris);
    }
    if (total_delta + delta_t >= 10.0 * min_delta) {
      increasing = 0;
    }
    if (increasing) {
      delta_t *= 2.0; /* increase time step and search further out */
    }
    else /* decreasing - reduce time step */
    {
      delta_t *= 0.5;
    }
    done = delta_t < min_dt;
  }

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON && fp) {
    fclose(fp);
  }

  return (total_delta);
}


/*
  Function to find pial surface vertex locations that minimize the
  specified energy functional. Note that spherical coordinates, white and
  pial socoordinates all need to be loaded, and that pial coordinates
  will be replaced with the new ones computed by this function.
*/

int MRISminimizeThicknessFunctional(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, float max_thick)
{
  int vno, navgs = 5;
  VERTEX *v;
  double ending_sse, vsize;

  //  parms->integration_type = INTEGRATE_MOMENTUM ;
  //  parms->integration_type = INTEGRATE_LM_SEARCH ;
  parms->flags |= IPFLAG_NOSCALE_TOL;  // don't scale tol down with decreasing # of averages
  parms->niterations = 10000;
  parms->projection = PROJECT_ELLIPSOID;

  MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, navgs);  // will be used in normal term
  MRIScomputeSurfaceNormals(mris, PIAL_VERTICES, navgs);  // will be used in normal term

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON) {
    //    int  n ;
    char fname[STRLEN], base_name[STRLEN];
    FileNameRemoveExtension(mris->fname, base_name);
#if 0
    for (n = 0 ; n < 100 ; n+=5)
    {
      MRIScomputeSurfaceNormals(mris, WHITE_VERTICES, n) ;  // will be used in normal term
      sprintf(fname, "%s.wnormals.a%3.3d.%s.mgz", base_name, n, parms->base_name) ;
      printf("writing normals to %s\n", fname) ;
      MRISwriteWhiteNormals(mris, fname) ;
    }
#endif
    int req = snprintf(fname, STRLEN, "%s.wnormals.a%3.3d.%s.mgz",
		       base_name, navgs, parms->base_name);   
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    printf("writing smoothed white matter surface normals to %s\n", fname);
    MRISwriteWhiteNormals(mris, fname);
  }

  MRISstoreMetricProperties(mris);
  //  mris->status = MRIS_SPHERE ;
  vsize = (mris->vg.xsize + mris->vg.ysize + mris->vg.zsize) / 3;
  mris->mht = parms->mht =
      (void *)MHTcreateFaceTable_Resolution(mris, CANONICAL_VERTICES, vsize);  // to lookup closest face

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);  
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->fp) {
      if (!parms->start_t) {
        INTEGRATION_PARMS_openFp(parms, fname, "w");
      }
      else {
        INTEGRATION_PARMS_openFp(parms, fname, "a");
      }

      if (!parms->fp) ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n", fname);
    }
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }

/*
  the first time through if we don't care about negative vertices, start things
  out with min dist (which will cause negative triangles).

*/
#if 1                                                    // NJS: else code causes mris_expand test failure Aug2012
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);  // start with pial coord of this vtx
#else
#if 1
  if (parms->start_t == 0)  // && parms->remove_neg == 0
  {
    MRISfindClosestPialVerticesCanonicalCoords(mris, mris->nsize);
  }
  else
#endif
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES);  // start with pial coord of this vtx
  mris->status = MRIS_SPHERE;
#endif
  MRIScomputeSecondFundamentalForm(mris);
  MRISrestoreRipFlags(mris);
  mrisClearMomentum(mris);
  if (parms->l_thick_min > 0 || parms->l_thick_normal > 0 || parms->l_thick_spring > 0 || parms->l_tsmooth > 0 ||
      parms->l_thick_parallel > 0)
    mrisIntegrationEpoch(mris, parms, 0);
  parms->niterations = 150;
  if (parms->remove_neg) {
    MRISremoveOverlapWithSmoothing(mris, parms);
  }
  ending_sse = MRIScomputeSSE(mris, parms);
  printf("ending sse = %f\n", ending_sse);

  // compute thickness and put it in v->curv field
  for (vno = 0; vno < mris->nvertices; vno++) {
    float xw, yw, zw, xp, yp, zp, thick;

    v = &mris->vertices[vno];
    if (vno == Gdiag_no) DiagBreak();

    if (v->ripflag) {
      v->tx = v->whitex;
      v->ty = v->whitey;
      v->tz = v->whitez;
      continue;
    }
    MRISvertexCoord2XYZ_float(v, WHITE_VERTICES, &xw, &yw, &zw);
    MRISsampleFaceCoordsCanonical((MHT *)(parms->mht), mris, v->x, v->y, v->z, PIAL_VERTICES, &xp, &yp, &zp);
    thick = sqrt(SQR(xp - xw) + SQR(yp - yw) + SQR(zp - zw));
    v->curv = thick;
    v->tx = xp;
    v->ty = yp;
    v->tz = zp;
    v->pialx = xp;
    v->pialy = yp;
    v->pialz = zp;
  }

  {
    MHT *mht = ((MHT *)(parms->mht));
    MHTfree(&mht);
    parms->mht = NULL;
  }
  return (NO_ERROR);
}


/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
static float area_coefs[] = {1.0f, 1.0f, 0.1f};
static float dist_coefs[] = {0.1f, 1.0f, 1.0f};

#define NCOEFS sizeof(area_coefs) / sizeof(area_coefs[0])

#define MAX_NBHD_SIZE 200
#define NBR_COEF (M_PI * 1.0f)

MRI_SURFACE *MRISunfold(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int base_averages, i, nbrs[MAX_NBHD_SIZE], niter, passno, msec, use_nl_area;
  double starting_sse, ending_sse, l_area, pct_error;

  printf("MRISunfold() max_passes = %d -------\n", max_passes);
  mrisLogIntegrationParms(stdout, mris, parms);
  printf("--------------------\n");

  use_nl_area = (!FZERO(parms->l_nlarea));

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }
  Timer start;
  starting_sse = ending_sse = 0.0f; /* compiler warning */
  memset(nbrs, 0, MAX_NBHD_SIZE * sizeof(nbrs[0]));
#if 0
  if (mris->nsize < 2)
  {
    nbrs[2] = nint(NBR_COEF*2.0) ;
  }
  for (i = 4 ; i <= parms->nbhd_size ; i*= 2)
  {
    nbrs[i] = nint(NBR_COEF*(float)i) ;
  }
#else
  for (i = mris->nsize + 1; i <= parms->nbhd_size; i++) {
    nbrs[i] = parms->max_nbrs;
  }
#endif

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name); 
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->fp) {
      if (!parms->start_t) {
        INTEGRATION_PARMS_openFp(parms, fname, "w");
      }
      else {
        INTEGRATION_PARMS_openFp(parms, fname, "a");
      }

      if (!parms->fp) ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n", fname);
    }
    mrisLogIntegrationParms(parms->fp, mris, parms);
    int valid_vertices = MRISvalidVertices(mris);
    if (parms->complete_dist_mat)
      fprintf(parms->fp, "using complete distance matrix (%d x %d)\n", valid_vertices, valid_vertices);
    else {
      for (i = mris->nsize + 1; i <= parms->nbhd_size; i++)
        if (nbrs[i]) {
          fprintf(parms->fp, "%d: %d | ", i, nbrs[i]);
        }
      fprintf(parms->fp, "\n");
    }
  }
  if (Gdiag & DIAG_SHOW) {
    int valid_vertices = MRISvalidVertices(mris);
    if (parms->complete_dist_mat)
      printf("using complete distance matrix (%d x %d)\n", valid_vertices, valid_vertices);
    else {
      for (i = mris->nsize + 1; i <= parms->nbhd_size; i++)
        if (nbrs[i]) {
          fprintf(stdout, "%d: %d | ", i, nbrs[i]);
        }
      fprintf(stdout, "\n");
    }
    mrisLogIntegrationParms(stderr, mris, parms);
  }

  /*  parms->start_t = 0 ;*/
  /*
    integrate until no improvement can be made at ANY scale, or until
    the error is effectively zero.
  */
  base_averages = parms->n_averages;
  l_area = parms->l_area;
  niter = parms->niterations;
  passno = 0;
  for (passno = 0; passno < max_passes; passno++) {
#if 0
    if (mris->nsize < parms->nbhd_size)  /* resample distances on surface */
#endif
    {
      if (Gdiag & DIAG_SHOW) {
        fprintf(stdout, "resampling long-range distances");
      }
      MRISsaveVertexPositions(mris, TMP_VERTICES);
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
      MRIScomputeMetricProperties(mris);
      if (parms->complete_dist_mat) {
        MRIScomputeAllDistances(mris);
      }
      else {
        MRISsampleDistances(mris, nbrs, parms->nbhd_size);
      }
      MRISrestoreVertexPositions(mris, TMP_VERTICES);
      MRIScomputeMetricProperties(mris);
      mrisClearMomentum(mris);
    }

    {
      char *cp;
      int vno, n;
      float d;
      FILE *fp;

      cp = getenv("FS_MEASURE_DISTANCES");
      if (cp) {
        fprintf(stdout, "outputting distance errors to distance.log...\n");
        fp = fopen("distance.log", "w");
        for (vno = 0; vno < mris->nvertices; vno++) {
          VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
          VERTEX          const * const v  = &mris->vertices         [vno];
          if (v->ripflag) {
            continue;
          }
          for (n = 0; n < vt->vtotal; n++) {
            VERTEX const * const vn = &mris->vertices[vt->v[n]];
            if (vn->ripflag) {
              continue;
            }
            d = sqrt(SQR(vn->origx - v->origx) + SQR(vn->origy - v->origy) + SQR(vn->origz - v->origz));
            fprintf(fp, "%2.4f  %2.4f\n", d, v->dist_orig[n]);
          }
        }
        fclose(fp);
        exit(1);
      }
    }
    {
      char *cp;
      double max_pct;

      cp = getenv("FS_DISTURB_DISTANCES");
      if (cp) {
        max_pct = atof(cp);
        fprintf(stdout, "disturbing distances by %%%2.1f\n", (float)max_pct);
        if (Gdiag & DIAG_WRITE) fprintf(parms->fp, "disturbing distances by %%%2.1f\n", (float)max_pct);
        MRISdisturbOriginalDistances(mris, max_pct);
      }
    }

    {
      char *cp;

      cp = getenv("FS_SPHERE");
      if (cp) {
        MRISstoreAnalyticDistances(mris, MRIS_SPHERE);
      }
      cp = getenv("FS_PLANE");
      if (cp) {
        MRISstoreAnalyticDistances(mris, MRIS_PLANE);
      }
    }

    //    if (!passno && ((parms->flags & IPFLAG_QUICK) == 0))
    if (((parms->flags & IPFLAG_QUICK) == 0)) {
      double tol = parms->tol;
      parms->tol = 0.5;
      if (niter > 30) {
        parms->niterations = 30;
      }
      printf("  mrisRemoveNegativeArea()\n");
      mrisRemoveNegativeArea(mris, parms, base_averages, MAX_NEG_AREA_PCT, 2);
      parms->niterations = niter;
      parms->tol = tol;
    }

    for (i = 0; (unsigned)i < NCOEFS; i++) {
      if (mris->status == MRIS_SPHERE && i == NCOEFS - 1) {
        continue;
      }

      pct_error = MRISpercentDistanceError(mris);
      if (Gdiag & DIAG_WRITE)
        fprintf(parms->fp,
                "pass %d: epoch %d of %d starting distance error %%%2.2f\n",
                passno,
                i + 1,
                (int)(NCOEFS),
                (float)pct_error);
      fprintf(stdout,
              "pass %d: epoch %d of %d starting distance error %%%2.2f\n",
              passno + 1,
              i + 1,
              (int)(NCOEFS),
              (float)pct_error);

      parms->l_dist = dist_coefs[i];
      if (use_nl_area) {
        parms->l_nlarea = area_coefs[i];
      }
      else {
        parms->l_area = area_coefs[i];
      }
      parms->l_angle = ANGLE_AREA_SCALE * parms->l_area;
      if (i == NCOEFS - 1) /* see if distance alone
                                can make things better */
      {
        starting_sse = MRIScomputeSSE(mris, parms);
      }
      mrisIntegrationEpoch(mris, parms, base_averages);
    }

    if (use_nl_area) {
      parms->l_nlarea = area_coefs[NCOEFS - 1];
    }
    else {
      parms->l_area = area_coefs[NCOEFS - 1];
    }
    parms->l_dist = dist_coefs[NCOEFS - 1];
    ending_sse = MRIScomputeSSE(mris, parms);
    if (Gdiag & DIAG_SHOW) {
#if 0
      fprintf(stdout, "pass %d: start=%2.1f, end=%2.1f, ratio=%2.3f\n",
              passno+1, starting_sse, ending_sse,
              (starting_sse-ending_sse)/starting_sse) ;
      if (Gdiag & DIAG_WRITE)
        fprintf
        (parms->fp,
         "pass %d: start=%2.4f, end=%2.4f, ratio=%2.4f\n",
         passno+1, starting_sse, ending_sse,
         (starting_sse-ending_sse)/starting_sse) ;
#endif
    }
  }
#if 0
  while (
    !FZERO(ending_sse) &&
    (((starting_sse-ending_sse)/starting_sse) > parms->tol) &&
    (++passno < max_passes)
  )
  {
    ;
  }
#endif

  fprintf(stdout, "unfolding complete - removing small folds...\n");
  pct_error = MRISpercentDistanceError(mris);
  if (Gdiag & DIAG_WRITE) fprintf(parms->fp, "starting distance error %%%2.2f\n", (float)pct_error);
  fprintf(stdout, "starting distance error %%%2.2f\n", (float)pct_error);

  /* finally, remove all the small holes */
  parms->l_nlarea = 1.0f;
  parms->l_area = 0.0;
  parms->l_dist = 0.1f; /* was 0.001 */
  parms->l_angle = ANGLE_AREA_SCALE * parms->l_nlarea;
  parms->niterations = niter;
  parms->tol = 1e-2; // try and remove as much negative stuff as possible

  // ATH: this is a pretty ugly hack that's caused much despair, it should really be cleaned up
  mrisStoreVtotalInV3num(mris); // hack to speed up neg. area removal
  fprintf(stdout, "removing remaining folds...\n");
  float incorrect_avg_nbrs = mris->avg_nbrs;
  MRISresetNeighborhoodSize(mris, 1);
  if (getenv("MRIS_SPHERE_NEW_BEHAVIOR") == nullptr) mris->avg_nbrs = incorrect_avg_nbrs;
  mrisRemoveNegativeArea(mris, parms, base_averages > 32 ? 32 : base_averages, MAX_NEG_AREA_PCT, 2);
  MRISresetNeighborhoodSize(mris, 3);

  if (mris->status == MRIS_PLANE && parms->complete_dist_mat == 0) /* smooth out remaining folds */
  {
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "smoothing final surface...\n");
    }
    if (Gdiag & DIAG_WRITE) {
      fprintf(parms->fp, "smoothing final surface...\n");
    }
    parms->l_spring = 1.0f;
    parms->l_area = parms->l_nlarea = 0.0f;
    parms->niterations = 5;
    parms->integration_type = INTEGRATE_MOMENTUM;
    parms->dt = 0.5f;
    parms->momentum = 0.0f;
    parms->n_averages = 0;
    MRISintegrate(mris, parms, 0);
    /* mrisRemoveNegativeArea(mris, parms, 0, MAX_NEG_AREA_PCT, 1);*/
  }

  pct_error = MRISpercentDistanceError(mris);
  fprintf(stdout, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris);
  msec = start.milliseconds();
  if (Gdiag & DIAG_SHOW) {
    mrisLogStatus(mris, parms, stderr, 0, -1);
    fprintf(stdout, "optimization complete.\n");
    fprintf(stdout, "unfolding took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
  }
  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->fp, "unfolding took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
    mrisLogStatus(mris, parms, parms->fp, 0, -1);
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
    INTEGRATION_PARMS_closeFp(parms);
  }
  printf("MRISunfold() return, current seed %ld\n", getRandomSeed());
  fflush(stdout);

  return (mris);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
MRI_SURFACE *MRISquickSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int max_passes)
{
  int niter, passno, msec, nbrs[MAX_NBHD_SIZE], i, use_dists, base_averages;
  double pct_error, orig_k, last_sse, sse, pct_change;

  orig_k = NEG_AREA_K;

  Timer start;

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }

  use_dists = (!DZERO(parms->l_dist) || !DZERO(parms->l_nldist)) && (parms->nbhd_size > mris->nsize);

  memset(nbrs, 0, MAX_NBHD_SIZE * sizeof(nbrs[0]));
  for (i = mris->nsize + 1; i <= parms->nbhd_size; i++) {
    nbrs[i] = parms->max_nbrs;
  }

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->fp) {
      if (!parms->start_t) {
        INTEGRATION_PARMS_openFp(parms, fname, "w");
      }
      else {
        INTEGRATION_PARMS_openFp(parms, fname, "a");
      }
      if (!parms->fp) ErrorExit(ERROR_NOFILE, "MRISquickSphere: could not open log file %s\n", fname);
    }
    mrisLogIntegrationParms(parms->fp, mris, parms);
    if (use_dists) {
      for (i = mris->nsize + 1; i <= parms->nbhd_size; i++)
        if (nbrs[i]) {
          fprintf(parms->fp, "%d: %d | ", i, nbrs[i]);
        }
      fprintf(parms->fp, "\n");
    }
  }
  if (Gdiag & DIAG_SHOW) {
    if (use_dists) {
      for (i = mris->nsize + 1; i <= parms->nbhd_size; i++)
        if (nbrs[i]) {
          fprintf(stdout, "%d: %d | ", i, nbrs[i]);
        }
      fprintf(stdout, "\n");
    }
    mrisLogIntegrationParms(stderr, mris, parms);
  }

  /* resample distances on surface */
  if (use_dists && mris->nsize < parms->nbhd_size) {
    if (Gdiag & DIAG_SHOW) {
      fprintf(stdout, "resampling long-range distances");
    }
    MRISsaveVertexPositions(mris, TMP_VERTICES);
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES);
    MRISsampleDistances(mris, nbrs, parms->nbhd_size);
    MRISrestoreVertexPositions(mris, TMP_VERTICES);
    mrisClearMomentum(mris);
  }

  /*
    integrate until no improvement can be made at ANY scale, or until
    the error has asymptoted.
  */
  base_averages = parms->n_averages;
  parms->flags |= IP_RETRY_INTEGRATION;
  niter = parms->niterations;
  passno = 0;

#if 1
  if ((parms->flags & IPFLAG_QUICK) == 0) {
    parms->tol = parms->tol * 1024 / (sqrt((double)base_averages + 1));
  }
#endif

  for (i = 0, NEG_AREA_K = orig_k; i < 4; NEG_AREA_K *= 4, i++) {
    passno = 0;
    do {
      last_sse = MRIScomputeSSE(mris, parms);
      printf("epoch %d (K=%2.1f), pass %d, starting sse = %2.2f\n", i + 1, NEG_AREA_K, passno + 1, last_sse);
      niter = mrisIntegrationEpoch(mris, parms, base_averages);
      sse = MRIScomputeSSE(mris, parms);
      pct_change = (last_sse - sse) / (last_sse * niter); /* per time step */
      passno++;
      printf("pass %d complete, delta sse/iter = %2.2f/%d = %2.5f\n",
             passno,
             (last_sse - sse) / last_sse,
             niter,
             pct_change);
    } while (pct_change > parms->tol);
#if 0
    if (passno == 1)   /* couldn't make any progress at all */
    {
      break ;
    }
#endif
  }

  NEG_AREA_K = orig_k;
  pct_error = MRISpercentDistanceError(mris);
  fprintf(stdout, "final distance error %%%2.2f\n", (float)pct_error);
  mrisProjectSurface(mris);
  msec = start.milliseconds();
  if (Gdiag & DIAG_SHOW) {
    mrisLogStatus(mris, parms, stderr, 0, -1);
    fprintf(stdout, "optimization complete.\n");
    fprintf(stdout, "unfolding took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
  }
  if (Gdiag & DIAG_WRITE) {
    fprintf(parms->fp, "unfolding took %2.2f hours\n", (float)msec / (1000.0f * 60.0f * 60.0f));
    mrisLogStatus(mris, parms, parms->fp, 0, -1);
    fprintf(parms->fp, "final distance error %%%2.2f\n", pct_error);
#if 0
    fclose(parms->fp) ;
#endif
  }

  return (mris);
}



int MRISinflateBrain(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int write_iterations = parms->write_iterations;

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }

#if 0
  {
    FACE *face ;
    int eno, nvertices, nfaces, nedges, v0, v1, v2 ;

    
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;
    face = &mris->faces[0] ;

    v0 = face->v[0] ;
    v1 = face->v[1] ;
    v2 = face->v[2] ;
    mrisRemoveLink(mris, v0, v1) ;
    MRISremoveRipped(mris) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "after removing %d --> %d, euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            v0, v1, nvertices, nedges, nfaces, eno, 2-eno) ;
    mrisRemoveLink(mris, v0, v2) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "after removing %d --> %d, euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            v0, v2, nvertices, nedges, nfaces, eno, 2-eno) ;
    MRISremoveRipped(mris) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "after removing %d --> %d, euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            v0, v2, nvertices, nedges, nfaces, eno, 2-eno) ;

    mrisRemoveLink(mris, v1, v2) ;
    MRISremoveRipped(mris) ;
    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "after removing %d --> %d, euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            v1, v2, nvertices, nedges, nfaces, eno, 2-eno) ;
    DiagBreak() ;
  }
#endif

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];

    int req = snprintf(fname, STRLEN, "%s.out", parms->base_name);
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->start_t) {
      INTEGRATION_PARMS_openFp(parms, fname, "w");
    }
    else {
      INTEGRATION_PARMS_openFp(parms, fname, "a");
    }
    if (!parms->fp) ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname, fname);
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  if (Gdiag & DIAG_SHOW) {
    mrisLogIntegrationParms(stderr, mris, parms);
  }

  MRIScomputeMetricProperties(mris);    // changes XYZ
  
  int    const niterations        = parms->niterations;
  double const desired_rms_height = parms->desired_rms_height;
  // fprintf(stdout, "inflating to desired rms height = %2.3f\n",
  // desired_rms_height);

  /* write out initial surface */
  if (!parms->start_t && (parms->write_iterations > 0) && (Gdiag & DIAG_WRITE)) {
    mrisWriteSnapshot(mris, parms, 0);
  }

  bool useOldBehaviour = false;
  if (useOldBehaviour) {
    switch (copeWithLogicProblem("FREESURFER_fix_inflateBrain","should set origx et al here")) {
    case LogicProblemResponse_old: 
      break;
    case LogicProblemResponse_fix:
      useOldBehaviour = false;
    }
  }
  if (!useOldBehaviour) {
    MRISsetOriginalXYZfromXYZ(mris);
    mrisComputeOriginalVertexDistances(mris);
  }
  
  MRIScomputeSSE(mris, parms);        // WHAT DOES THIS ACHIEVE?
  double rms_height = MRISrmsTPHeight(mris);
  
  if (!parms->start_t) {
    if (Gdiag & DIAG_SHOW)
      fprintf(stdout, "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d\n", 0, 0.0f, desired_rms_height, parms->n_averages);
    else
      fprintf(stdout, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", 0, rms_height, desired_rms_height);
    if (Gdiag & DIAG_WRITE) {
      fprintf(
          parms->fp, "%3.3d: dt: %2.4f, rms height=%2.3f, avgs=%d", 0, 0.0f, rms_height, parms->n_averages);
      fflush(parms->fp);
    }

    MRISclearCurvature(mris); /* curvature will be used to calculate sulc */
  }

  double const l_dist = parms->l_dist;
  
  int n_averages;
  for (n_averages = parms->n_averages; n_averages >= 0; n_averages /= 2) {
    parms->l_dist = l_dist * sqrt(n_averages);

    int n;
    for (n = parms->start_t; n < parms->start_t + niterations; n++) {
      mrisTearStressedRegions(mris, parms)  ;
      
      if (parms->explode_flag)
      {
        MRISsetOrigArea(mris);  // used to happen inside MRISrenumberRemovingRippedFacesAndVertices
	MRISrenumberRemovingRippedFacesAndVertices(mris);
      }
  
      MRISclearGradient(mris);
      mrisComputeDistanceTerm(mris, parms);
      mrisComputeSphereTerm(mris, parms->l_sphere, parms->a, parms->explode_flag);
      mrisComputeExpansionTerm(mris, parms->l_expand);

      MRISaverageGradients(mris, n_averages);
      mrisComputeNormalSpringTerm(mris, parms->l_nspring);
      mrisComputeNonlinearSpringTerm(mris, parms->l_nlspring, parms);
      mrisComputeTangentialSpringTerm(mris, parms->l_tspring);
      mrisComputeNonlinearTangentialSpringTerm(mris, parms->l_nltspring, parms->min_dist);
      mrisComputeQuadraticCurvatureTerm(mris, parms->l_curv);
      mrisComputeSpringTerm(mris, parms->l_spring);
      mrisComputeLaplacianTerm(mris, parms->l_lap);
      mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm);
      
      double delta_t;
      switch (parms->integration_type) {
        case INTEGRATE_LM_SEARCH:
          delta_t = mrisLineMinimizeSearch(mris, parms);
          break;
        default:
        case INTEGRATE_LINE_MINIMIZE:
          delta_t = mrisLineMinimize(mris, parms);
          break;
        case INTEGRATE_MOMENTUM:
          delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, parms->tol, 0 /*parms->n_averages*/);
          break;
        case INTEGRATE_ADAPTIVE:
          delta_t = mrisAdaptiveTimeStep(mris, parms);
          break;
      }
      
      mrisTrackTotalDistanceNew(mris); /* update sulc */
      MRIScomputeMetricProperties(mris);
      
      MRIScomputeSSE(mris, parms);  // WHAT DOES THIS ACHIEVE?
      
      if (0) mris_print_hash(stdout, mris, "\nBefore calling MRISrmsTPHeight", "\n");
      rms_height = MRISrmsTPHeight(mris);
      if (0) mris_print_hash(stdout, mris, "\nAfter calling MRISrmsTPHeight",  "\n");
      
      if (!((n + 1) % 5)) /* print some diagnostics */
      {
        if (Gdiag & DIAG_SHOW)
          fprintf(stdout,
                  "%3.3d: dt: %2.4f, "
                  "rms height=%2.3f, avgs=%d, l_dist=%2.2f\n",
                  n + 1,
                  (float)delta_t,
                  (float)rms_height,
                  n_averages,
                  parms->l_dist);
        else {
          fprintf(stdout, "\rstep %3.3d: RMS=%2.3f (target=%2.3f)   ", n + 1, rms_height, desired_rms_height);
          fflush(stdout);
        }
        if (Gdiag & DIAG_WRITE) {
          fprintf(parms->fp,
                  "%3.3d: dt: %2.4f, rms height=%2.3f, "
                  "avgs=%d, l_dist=%2.2f\n",
                  n + 1,
                  (float)delta_t,
                  (float)rms_height,
                  n_averages,
                  parms->l_dist);
          fflush(parms->fp);
        }
      }

      if (parms->scale > 0) {
        MRIScomputeMetricProperties(mris);
        printf(
            "rescaling brain to retain original "
            "surface area %2.0f (%2.2f), current %2.1f\n",
            mris->orig_area,
            sqrt(mris->orig_area / (mris->total_area + mris->neg_area)), mris->total_area);
        MRISscaleBrainArea(mris);
        MRIScomputeMetricProperties(mris);
//	printf("after rescaling, surface area %2.1f\n", mris->total_area);
	MRISprintTessellationStats(mris, stderr);
      }
      if ((parms->write_iterations > 0) && !((n + 1) % write_iterations) && (Gdiag & DIAG_WRITE)) {
        mrisWriteSnapshot(mris, parms, n + 1);
      }
      if (rms_height < desired_rms_height) {
        break;
      }
    }

    parms->start_t = n;
    if (!n_averages || rms_height < desired_rms_height) {
      break;
    }
  }

  fprintf(stdout, "\ninflation complete.\n");
  if (Gdiag & DIAG_WRITE) {
    INTEGRATION_PARMS_closeFp(parms);
  }

  return (NO_ERROR);
}

/*-----------------------------------------------------
  Parameters:

  Returns value:

  Description
  ------------------------------------------------------*/
int MRISinflateToSphere(MRI_SURFACE *mris, INTEGRATION_PARMS *parms)
{
  int n_averages, n, write_iterations, niterations, base_averages;
  double delta_t = 0.0, rms_radial_error, sse, base_dt;
  MHT *mht_v_current = NULL;

  printf("Entering MRISinflateToSphere()\n");

  if (IS_QUADRANGULAR(mris)) {
    MRISremoveTriangleLinks(mris);
  }
  write_iterations = parms->write_iterations;
  n_averages = parms->n_averages;

  if (Gdiag & DIAG_WRITE) {
    char fname[STRLEN];
    int req = snprintf(fname, STRLEN, "%s.%s.out",
		       mris->hemisphere == RIGHT_HEMISPHERE ? "rh" : "lh", parms->base_name);  
    if( req >= STRLEN ) {   
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!parms->fp) {
      if (!parms->start_t) 
        INTEGRATION_PARMS_openFp(parms, fname, "w");
      else
        INTEGRATION_PARMS_openFp(parms, fname, "a");
      if(!parms->fp) ErrorExit(ERROR_NOFILE, "MRISunfold: could not open log file %s\n", fname);
    }
    mrisLogIntegrationParms(parms->fp, mris, parms);
  }
  if(Gdiag & DIAG_SHOW) 
    mrisLogIntegrationParms(stderr, mris, parms);

  MRIScomputeMetricProperties(mris);

  /*  parms->start_t = 0 ;*/
  niterations = parms->niterations;
  MRISstoreMetricProperties(mris);
  rms_radial_error = sqrt(mrisComputeSphereError(mris, 1.0, parms->a) / mris->nvertices);
  fprintf(stdout, "inflating to sphere (rms error < %2.2f)\n", parms->tol * parms->a / 100.0f);

  /* write out initial surface */
  if (!parms->start_t && (parms->write_iterations > 0) && (Gdiag & DIAG_WRITE)) {
    mrisWriteSnapshot(mris, parms, 0);
  }

  sse = MRIScomputeSSE(mris, parms);
  if (!parms->start_t) {
    fprintf(
        stdout, "%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 0, 0.0f, (float)rms_radial_error, n_averages);
    if (Gdiag & DIAG_WRITE) {
      fprintf(parms->fp, "%3.3d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n", 0, 0.0f, rms_radial_error, n_averages);
      fflush(parms->fp);
    }

    MRISclearCurvature(mris); /* curvature will be used to
                       calculate sulc */
  }

  base_averages = parms->n_averages;
  base_dt = parms->dt;
  for (n_averages = base_averages; n_averages >= 0; n_averages /= 4) {
    parms->n_averages = n_averages;
    parms->dt = (sqrt((float)n_averages) + 1) * base_dt;
    for (n = parms->start_t; n < parms->start_t + niterations; n++) {
      if (!FZERO(parms->l_repulse_ratio)) {
        MHTfree(&mht_v_current);
        mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 3.0f);
      }
      MRISclearGradient(mris);
      mrisComputeDistanceTerm(mris, parms);
      mrisComputeSphereTerm(mris, parms->l_sphere, parms->a, parms->explode_flag);
      mrisComputeExpansionTerm(mris, parms->l_expand);
      mrisComputeRepulsiveRatioTerm(mris, parms->l_repulse_ratio, mht_v_current);
      mrisComputeConvexityTerm(mris, parms->l_convex);

      mrisComputeLaplacianTerm(mris, parms->l_lap);
      mrisComputeNonlinearSpringTerm(mris, parms->l_nlspring, parms);
      mrisComputeTangentialSpringTerm(mris, parms->l_tspring);
      mrisComputeNonlinearTangentialSpringTerm(mris, parms->l_nltspring, parms->min_dist);
      MRISaverageGradients(mris, n_averages);
      mrisComputeSpringTerm(mris, parms->l_spring);
      mrisComputeNormalizedSpringTerm(mris, parms->l_spring_norm);
      switch (parms->integration_type) {
        case INTEGRATE_LM_SEARCH:
          delta_t = mrisLineMinimizeSearch(mris, parms);
          break;
        default:
        case INTEGRATE_LINE_MINIMIZE:
          delta_t = mrisLineMinimize(mris, parms);
          break;
        case INTEGRATE_MOMENTUM:
          delta_t = MRISmomentumTimeStep(mris, parms->momentum, parms->dt, parms->tol, 0 /*parms->n_averages*/);
          break;
        case INTEGRATE_ADAPTIVE:
          delta_t = mrisAdaptiveTimeStep(mris, parms);
          break;
      }
      mrisTrackTotalDistance(mris); /* update sulc */
      MRIScomputeMetricProperties(mris);
      sse = MRIScomputeSSE(mris, parms);
      rms_radial_error = sqrt(mrisComputeSphereError(mris, 1.0, parms->a) / mris->nvertices);
      if (!((n + 1) % 5)) /* print some diagnostics */
      {
        fprintf(stdout,
                "%3.3d/%d: dt: %2.4f, rms radial error=%2.3f, avgs=%d\n",
                n + 1,
                parms->start_t + niterations,
                (float)delta_t,
                (float)rms_radial_error,
                n_averages);
        if (Gdiag & DIAG_WRITE) {
          fprintf(parms->fp,
                  "%3.3d: dt: %2.4f, rms radial error=%2.3f, "
                  "avgs=%d\n",
                  n + 1,
                  (float)delta_t,
                  (float)rms_radial_error,
                  n_averages);
          fflush(parms->fp);
        }
      }

      if ((parms->write_iterations > 0) && !((n + 1) % write_iterations) && (Gdiag & DIAG_WRITE)) {
        mrisWriteSnapshot(mris, parms, n + 1);
      }
      if (100.0 * rms_radial_error / parms->a < parms->tol) {
        break;
      }
    }
    parms->start_t = n;
    if (!n_averages || (100.0 * rms_radial_error / parms->a < parms->tol)) {
      break;
    }
    printf("reducing # of averages to %d\n", n_averages / 4);
  }

  fprintf(stdout, "\nspherical inflation complete.\n");
  if (Gdiag & DIAG_WRITE) {
    INTEGRATION_PARMS_closeFp(parms);
  }
  if (!FZERO(parms->l_repulse)) {
    MHTfree(&mht_v_current);
  }

  return (NO_ERROR);
}
