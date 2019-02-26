// Do not include, intended only to be included into mrisurf_deform.c

static double MRIScomputeSSE_new(MRIS *mris, INTEGRATION_PARMS *parms, bool debug);
static double MRIScomputeSSE_old(MRIS *mris, INTEGRATION_PARMS *parms, bool debug);

double MRIScomputeSSE(MRIS *mris, INTEGRATION_PARMS *parms) {
    static const bool useOld = false;
    static int count = 0;
    double old_sse = useOld ? MRIScomputeSSE_old(mris, parms, false) : 0.0;
    double new_sse =          MRIScomputeSSE_new(mris, parms, false);
    if (useOld && !closeEnough((float)old_sse,(float)new_sse)) {
        fprintf(stdout, "%s:%d count:%d old_sse:%f new_sse:%f\n", __FILE__, __LINE__, count, old_sse, new_sse);
        MRIScomputeSSE_old(mris, parms, true);
        MRIScomputeSSE_new(mris, parms, true);
        exit(1);
    }
    return new_sse;
}

// These are in the order the original code computed them, so that side effects are not reordered
// In older code the ashburner_triangle is computed but not used , here it is not computed at all
//
// Important ones are ...  
//      sse_area sse_neg_area       in this function                                here                iterates over faces computing face normals
//      sse_nl_area                 mrisComputeNonlinearAreaSSE(mris)               mrisurf_deform.c    iterates over faces doing a simple calculation
//      sse_dist                    mrisComputeNonlinearDistanceSSE(mris)           mrisurf_deform.c    iterates over vertices over their dist and dist_orig
//      sse_corr                    mrisComputeCorrelationError(mris, parms, 1)     mrisurf_deform.c    iterates over vertices using xyz calling mrisp.c MRISPfunctionValTraceable
//                                                                                                                          which does a lot more work than the others
//
#define COMPUTE_DISTANCE_ERROR mrisComputeDistanceError(mris, parms)

#define SSE_TERMS \
      ELT(sse_area                  , parms->l_parea,                            true,    computed_area                                                                   ) SEP \
      ELT(sse_neg_area              , parms->l_area,                             true,    computed_neg_area                                                               ) SEP \
      ELT(sse_repulse               , 1.0,                     (parms->l_repulse > 0),    mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current, mht_f_current)) SEP \
      ELT(sse_repulsive_ratio       , 1.0,                                       true,    mrisComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio)                   ) SEP \
      ELT(sse_tsmooth               , 1.0,                                       true,    mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms)             ) SEP \
      ELT(sse_thick_min             , parms->l_thick_min,                        true,    mrisComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms)         ) SEP \
      ELT(sse_ashburner_triangle    , parms->l_ashburner_triangle,               false,   mrisComputeAshburnerTriangleEnergy(mris, parms->l_ashburner_triangle, parms)    ) SEP \
      ELT(sse_thick_parallel        , parms->l_thick_parallel,                   true,    mrisComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms)        ) SEP \
      ELT(sse_thick_normal          , parms->l_thick_normal,                     true,    mrisComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms)            ) SEP \
      ELT(sse_thick_spring          , parms->l_thick_spring,                     true,    mrisComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms)            ) SEP \
      ELT(sse_nl_area               , parms->l_nlarea,        !FZERO(parms->l_nlarea),    mrisComputeNonlinearAreaSSE(mris)                                               ) SEP \
      ELT(sse_nl_dist               , parms->l_nldist,        !DZERO(parms->l_nldist),    mrisComputeNonlinearDistanceSSE(mris)                                           ) SEP \
      ELT(sse_dist                  , parms->l_dist,          !DZERO(parms->l_dist),      COMPUTE_DISTANCE_ERROR                                                          ) SEP \
      ELT(sse_spring                , parms->l_spring,        !DZERO(parms->l_spring),    mrisComputeSpringEnergy(mris)                                                   ) SEP \
      ELT(sse_lap                   , parms->l_lap,           !DZERO(parms->l_lap),       mrisComputeLaplacianEnergy(mris)                                                ) SEP \
      ELT(sse_tspring               , parms->l_tspring,       !DZERO(parms->l_tspring),   mrisComputeTangentialSpringEnergy(mris)                                         ) SEP \
      ELT(sse_nlspring              , parms->l_nlspring,      !DZERO(parms->l_nlspring),  mrisComputeNonlinearSpringEnergy(mris, parms)                                   ) SEP \
      ELT(sse_curv                  , l_curv_scaled,          !DZERO(parms->l_curv),      mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv)                           ) SEP \
      ELT(sse_corr                  , l_corr,                 !DZERO(l_corr),             mrisComputeCorrelationError(mris, parms, 1)                                     ) SEP \
      ELT(sse_val                   , parms->l_intensity,     !DZERO(parms->l_intensity), mrisComputeIntensityError(mris, parms)                                          ) SEP \
      ELT(sse_loc                   , parms->l_location,      !DZERO(parms->l_location),  mrisComputeTargetLocationError(mris, parms)                                     ) SEP \
      ELT(sse_dura                  , parms->l_dura,          !DZERO(parms->l_dura),      mrisComputeDuraError(mris, parms)                                               ) SEP \
      ELT(sse_histo                 , parms->l_histo,         !DZERO(parms->l_histo),     mrisComputeHistoNegativeLikelihood(mris, parms)                                 ) SEP \
      ELT(sse_map                   , parms->l_map,           !DZERO(parms->l_map),       mrisComputeNegativeLogPosterior(mris, parms, NULL)                              ) SEP \
      ELT(sse_map2d                 , parms->l_map2d,         !DZERO(parms->l_map2d),     mrisComputeNegativeLogPosterior2D(mris, parms, NULL)                            ) SEP \
      ELT(sse_grad                  , parms->l_grad,          !DZERO(parms->l_grad),      mrisComputeIntensityGradientError(mris, parms)                                  ) SEP \
      ELT(sse_sphere                , parms->l_sphere,        !DZERO(parms->l_sphere),    mrisComputeSphereError(mris, parms->l_sphere, parms->a)                         ) SEP \
      ELT(sse_shrinkwrap            , parms->l_shrinkwrap,    !DZERO(parms->l_shrinkwrap),mrisComputeShrinkwrapError(mris, parms->mri_brain, parms->l_shrinkwrap)         ) SEP \
      ELT(sse_expandwrap            , parms->l_expandwrap,    !DZERO(parms->l_expandwrap),mrisComputeExpandwrapError(mris, parms->mri_brain, parms->l_expandwrap, parms->target_radius)) SEP \
      ELT(sse_vectorCorrelationError, 1.0,                    use_multiframes,            mrisComputeVectorCorrelationError(mris, parms, 1)                               )     \
      // end of list

double MRIScomputeSSE_new(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, bool debug)
{
  bool   const use_multiframes  = !!(parms->flags & IP_USE_MULTIFRAMES);
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);
  double const l_curv_scaled    = (double)parms->l_curv * CURV_SCALE;
  double const area_scale =
#if METRIC_SCALE
    (mris->patch || mris->noscale) ? 1.0 : mris->orig_area / mris->total_area;
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

#ifdef BEVIN_MRISCOMPUTESSE_REPRODUCIBLE

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
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
    
    ROMP_PF_begin       // mris_register

#ifdef BEVIN_MRISCOMPUTESSE_CHECK
    #pragma omp parallel for if(trial==0) reduction(+ : relevant_angle, computed_neg_area, computed_area)
#else
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(fast) reduction(+ : relevant_angle, computed_neg_area, computed_area)
#endif
#endif
    for (fno = 0; fno < mris->nfaces; fno++) {
      ROMP_PFLB_begin

#endif      
      FACE const * const face = &mris->faces[fno];
      if (face->ripflag) ROMP_PF_continue;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

      {
        double const delta = (double)(area_scale * face->area - fNorm->orig_area);
#if ONLY_NEG_AREA_TERM
        if (face->area < 0.0f) computed_neg_area += delta * delta;
#endif
        computed_area += delta * delta;
      }
      
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        double delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f) delta = 0.0f;

#endif
        relevant_angle += delta * delta;
      }
      
      if (!isfinite(computed_area) || !isfinite(relevant_angle)) {
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
    vmean = MRIScomputeTotalVertexSpacingStats(mris, &vsigma, NULL, NULL, NULL, NULL);
    mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, vmean);
    mht_f_current = MHTcreateFaceTable_Resolution  (mris, CURRENT_VERTICES, vmean);
  }


#define SEP
#define ELT(NAME, MULTIPLIER, COND, EXPR) double const NAME = (COND) ? (EXPR) : 0.0;
    SSE_TERMS
#undef ELT
#undef SEP

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
  
  double sse = sse_init +
#define SEP +
#define ELT(NAME, MULTIPLIER, COND, EXPR) (MULTIPLIER) * (NAME)
    SSE_TERMS ;
#undef ELT
#undef SEP

  if (debug) {
    #define SEP
    double sum = 0;
    #define ELT(NAME, MULTIPLIER, COND, EXPR) fprintf(stdout, "new %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);
    ELT(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "new sum = %f \n", sum);
    #undef ELT
    #undef SEP
  }
  

  if (mht_v_current) MHTfree(&mht_v_current);
  if (mht_f_current) MHTfree(&mht_f_current);

  if (!devFinite(sse)) {
    DiagBreak();
  }

#undef COMPUTE_DISTANCE_ERROR

  return sse;
}


// OLD CODE

double MRIScomputeSSE_old(MRI_SURFACE *mris, INTEGRATION_PARMS *parms, bool debug)
{
  double const l_corr           = (double)(parms->l_corr + parms->l_pcorr);
  double const l_curv_scaled    = (double)parms->l_curv * CURV_SCALE;

  double 
      sse_curv, sse_spring, sse_dist, area_scale, sse_corr, sse_val,
      sse_sphere, sse_thick_min, sse_thick_parallel, sse_ashburner_triangle, sse_grad, sse_nl_area, sse_nl_dist,
      sse_tspring, sse_repulse, sse_tsmooth, sse_loc, sse_thick_spring, sse_repulsive_ratio, sse_shrinkwrap,
      sse_expandwrap, sse_lap, sse_dura, sse_nlspring, sse_thick_normal, sse_histo, sse_map, sse_map2d;

#if METRIC_SCALE
  if (mris->patch || mris->noscale) {
    area_scale = 1.0;
  }
  else {
    area_scale = mris->orig_area / mris->total_area;
  }
#else
  area_scale = 1.0;
#endif

  sse_repulse = sse_repulsive_ratio = sse_tsmooth = 
  sse_thick_parallel = sse_thick_normal = sse_thick_spring =
      sse_nl_area = sse_nl_dist = sse_corr = sse_ashburner_triangle = sse_val = sse_sphere =
          sse_shrinkwrap = sse_expandwrap = sse_dura = sse_histo = sse_lap = sse_spring = sse_curv =
              sse_dist = sse_tspring = sse_loc = sse_nlspring = sse_thick_min =
	      sse_map = sse_map2d = sse_grad = 0.0;

  MHT *mht_v_current = NULL;
  MHT *mht_f_current = NULL;
  if (!FZERO(parms->l_repulse)) {
    double vmean, vsigma;
    vmean = MRIScomputeTotalVertexSpacingStats(mris, &vsigma, NULL, NULL, NULL, NULL);
    mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, vmean);
    mht_f_current = MHTcreateFaceTable_Resolution  (mris, CURRENT_VERTICES, vmean);
  }

  double sse_angle = 0, sse_neg_area = 0, sse_area = 0;

  if (!FZERO(parms->l_angle) || !FZERO(parms->l_area) || (!FZERO(parms->l_parea))) {

#ifdef BEVIN_MRISCOMPUTESSE_CHECK
    int trial; 
    double sse_angle_trial0, sse_neg_area_trial0, sse_area_trial0;
    for (trial = 0; trial < 2; trial++) {

#endif

    sse_angle = 0; sse_neg_area = 0; sse_area = 0;

#ifdef BEVIN_MRISCOMPUTESSE_REPRODUCIBLE

  #define ROMP_VARIABLE       fno
  #define ROMP_LO             0
  #define ROMP_HI             mris->nfaces
    
  #define ROMP_SUMREDUCTION0  sse_angle
  #define ROMP_SUMREDUCTION1  sse_neg_area
  #define ROMP_SUMREDUCTION2  sse_area
    
  #define ROMP_FOR_LEVEL      ROMP_level_assume_reproducible
    
#ifdef ROMP_SUPPORT_ENABLED
  const int romp_for_line = __LINE__;
#endif
  #include "romp_for_begin.h"
  ROMP_for_begin
    
    #define sse_angle    ROMP_PARTIALSUM(0)
    #define sse_neg_area ROMP_PARTIALSUM(1)
    #define sse_area     ROMP_PARTIALSUM(2)

#else

    int fno;
    
    ROMP_PF_begin       // mris_register

#ifdef BEVIN_MRISCOMPUTESSE_CHECK
    #pragma omp parallel for if(trial==0) reduction(+ : sse_angle, sse_neg_area, sse_area)
#else
#ifdef HAVE_OPENMP
    #pragma omp parallel for if_ROMP(fast) reduction(+ : sse_angle, sse_neg_area, sse_area)
#endif
#endif
    for (fno = 0; fno < mris->nfaces; fno++) {
      ROMP_PFLB_begin

#endif      
      FACE *face;
      double delta;

      face = &mris->faces[fno];
      if (face->ripflag) ROMP_PF_continue;
      FaceNormCacheEntry const * const fNorm = getFaceNorm(mris, fno);

      delta = (double)(area_scale * face->area - fNorm->orig_area);
#if ONLY_NEG_AREA_TERM
      if (face->area < 0.0f) sse_neg_area += delta * delta;

#endif
      sse_area += delta * delta;
      int ano;
      for (ano = 0; ano < ANGLES_PER_TRIANGLE; ano++) {
        delta = deltaAngle(face->angle[ano], face->orig_angle[ano]);
#if ONLY_NEG_AREA_TERM
        if (face->angle[ano] >= 0.0f) delta = 0.0f;

#endif
        sse_angle += delta * delta;
      }
      if (!isfinite(sse_area) || !isfinite(sse_angle)) {
        ErrorExit(ERROR_BADPARM, "sse not finite at face %d!\n", fno);
      }
#ifdef BEVIN_MRISCOMPUTESSE_REPRODUCIBLE

    #undef sse_angle
    #undef sse_neg_area
    #undef sse_area

  #include "romp_for_end.h"

#else
      ROMP_PFLB_end
    }
    ROMP_PF_end
#endif
    
#ifdef BEVIN_MRISCOMPUTESSE_CHECK

    if (trial == 0) {
       
      sse_angle_trial0    = sse_angle;
      sse_neg_area_trial0 = sse_neg_area;
      sse_area_trial0     = sse_area;
    } else { 
      if (sse_angle_trial0 != sse_angle) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           sse_angle_trial0, sse_angle, sse_angle_trial0-sse_angle);
      }
      if (sse_neg_area_trial0 != sse_neg_area) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           sse_neg_area_trial0, sse_neg_area, sse_neg_area_trial0-sse_neg_area);
      }
      if (sse_area_trial0 != sse_area) {
        fprintf(stderr, "%s:%d diff thread count, diff result %g %g %g\n",__FILE__,__LINE__,
           sse_area_trial0, sse_area, sse_area_trial0-sse_area);
      }
    }
    
    } // trial
#endif

  }
  
  if (parms->l_repulse > 0)
    sse_repulse = mrisComputeRepulsiveEnergy(mris, parms->l_repulse, mht_v_current, mht_f_current);
  sse_repulsive_ratio = mrisComputeRepulsiveRatioEnergy(mris, parms->l_repulse_ratio);
  sse_tsmooth = mrisComputeThicknessSmoothnessEnergy(mris, parms->l_tsmooth, parms);
  sse_thick_min = mrisComputeThicknessMinimizationEnergy(mris, parms->l_thick_min, parms);
  sse_ashburner_triangle = mrisComputeAshburnerTriangleEnergy(mris, parms->l_ashburner_triangle, parms);
  sse_thick_parallel = mrisComputeThicknessParallelEnergy(mris, parms->l_thick_parallel, parms);
  sse_thick_normal = mrisComputeThicknessNormalEnergy(mris, parms->l_thick_normal, parms);
  sse_thick_spring = mrisComputeThicknessSpringEnergy(mris, parms->l_thick_spring, parms);
  if (parms->l_thick_spring > 0 || parms->l_thick_min > 0 || parms->l_thick_parallel > 0 /* && DIAG_VERBOSE_ON*/)
    printf("min=%2.3f, parallel=%2.4f, normal=%2.4f, spring=%2.4f, ashburner=%2.3f, tsmooth=%2.3f\n",
           sse_thick_min / (float)mris->nvertices,
           sse_thick_parallel / (float)mris->nvertices,
           sse_thick_normal / (float)mris->nvertices,
           sse_thick_spring / (float)mris->nvertices,
           sse_ashburner_triangle / (float)mris->nvertices,
           sse_tsmooth / (float)mris->nvertices);
  if (!FZERO(parms->l_nlarea)) {
    sse_nl_area = mrisComputeNonlinearAreaSSE(mris);
  }
  if (!DZERO(parms->l_nldist)) {
    sse_nl_dist = mrisComputeNonlinearDistanceSSE(mris);
  }
  if (!DZERO(parms->l_dist)) {
    sse_dist = mrisComputeDistanceError(mris, parms);
  }
  if (!DZERO(parms->l_spring)) {
    sse_spring = mrisComputeSpringEnergy(mris);
  }
  if (!DZERO(parms->l_lap)) {
    sse_lap = mrisComputeLaplacianEnergy(mris);
  }
  if (!DZERO(parms->l_tspring)) {
    sse_tspring = mrisComputeTangentialSpringEnergy(mris);
  }
  if (!DZERO(parms->l_nlspring)) {
    sse_nlspring = mrisComputeNonlinearSpringEnergy(mris, parms);
  }
  if (!DZERO(parms->l_curv)) {
    sse_curv = mrisComputeQuadraticCurvatureSSE(mris, parms->l_curv);
  }
  if (!DZERO(l_corr)) {
    sse_corr = mrisComputeCorrelationError(mris, parms, 1);
  }
  if (!DZERO(parms->l_intensity)) {
    sse_val = mrisComputeIntensityError(mris, parms);
  }
  if (!DZERO(parms->l_location)) {
    sse_loc = mrisComputeTargetLocationError(mris, parms);
  }
  if (!DZERO(parms->l_dura)) {
    sse_dura = mrisComputeDuraError(mris, parms);
  }
  if (!DZERO(parms->l_histo)) {
    sse_histo = mrisComputeHistoNegativeLikelihood(mris, parms);
  }
  if (!DZERO(parms->l_map)) {
    sse_map = mrisComputeNegativeLogPosterior(mris, parms, NULL);
  }
  if (!DZERO(parms->l_map2d)) {
    sse_map2d = mrisComputeNegativeLogPosterior2D(mris, parms, NULL);
  }
  if (!DZERO(parms->l_grad)) {
    sse_grad = mrisComputeIntensityGradientError(mris, parms);
  }
  if (!DZERO(parms->l_sphere)) {
    sse_sphere = mrisComputeSphereError(mris, parms->l_sphere, parms->a);
  }
  if (!DZERO(parms->l_shrinkwrap))
    sse_shrinkwrap = mrisComputeShrinkwrapError(mris, parms->mri_brain, parms->l_shrinkwrap);
  if (!DZERO(parms->l_expandwrap))
    sse_expandwrap = mrisComputeExpandwrapError(mris, parms->mri_brain, parms->l_expandwrap, parms->target_radius);

  double sse_init = 0;

  if (gMRISexternalSSE) {
    sse_init = (*gMRISexternalSSE)(mris, parms);
  }
  
  double sse = sse_init +
         (double)parms->l_area * sse_neg_area + 
         sse_repulse + 
         (double)parms->l_tsmooth * sse_tsmooth +
         (double)parms->l_thick_min * sse_thick_min + 
	 (double)parms->l_thick_parallel * sse_thick_parallel +
         (double)parms->l_thick_spring * sse_thick_spring + 
	 (double)parms->l_thick_normal * sse_thick_normal +
         (double)parms->l_sphere * sse_sphere + sse_repulsive_ratio + 
	 (double)parms->l_intensity * sse_val +
         (double)parms->l_location * sse_loc + 
	 (double)parms->l_lap * sse_lap +
         //    (double)parms->l_ashburner_triangle * sse_ashburner_triangle +
         (double)parms->l_shrinkwrap * sse_shrinkwrap + 
	 (double)parms->l_expandwrap * sse_expandwrap +
         (double)parms->l_map * sse_map + 
	 (double)parms->l_map2d * sse_map2d + 
	 (double)parms->l_grad * sse_grad +
         (double)parms->l_parea * sse_area + 
	 (double)parms->l_nldist * sse_nl_dist +
         (double)parms->l_nlarea * sse_nl_area + 
	 (double)parms->l_angle * sse_angle + 
	 (double)parms->l_dist * sse_dist +
         (double)parms->l_nlspring * sse_nlspring + 
	 (double)parms->l_spring * sse_spring +
	 (double)parms->l_histo * sse_histo +
         (double)parms->l_dura * sse_dura + 
	 (double)parms->l_tspring * sse_tspring + 
	 (double)l_corr * sse_corr + 
	 (double)parms->l_curv * CURV_SCALE * sse_curv;

  double sse_vectorCorrelationError = 0.0;
  if (parms->flags & IP_USE_MULTIFRAMES) {
    sse += sse_vectorCorrelationError = (double)mrisComputeVectorCorrelationError(mris, parms, 1);
  }
  if (debug) {
    fprintf(stdout, "old terms\n");
    fprintf(stdout, " %f\n", sse_init );
    fprintf(stdout, " %f\n", (double)parms->l_area * sse_neg_area );
    fprintf(stdout, " %f\n", sse_repulse );
    fprintf(stdout, " %f\n", (double)parms->l_tsmooth * sse_tsmooth );
    fprintf(stdout, " %f\n", (double)parms->l_thick_min * sse_thick_min );
    fprintf(stdout, " %f\n", (double)parms->l_thick_parallel * sse_thick_parallel );
    fprintf(stdout, " %f\n", (double)parms->l_thick_spring * sse_thick_spring );
    fprintf(stdout, " %f\n", (double)parms->l_thick_normal * sse_thick_normal );
    fprintf(stdout, " %f\n", (double)parms->l_sphere * sse_sphere );
    fprintf(stdout, " %f\n", sse_repulsive_ratio );
    fprintf(stdout, " %f\n", (double)parms->l_intensity * sse_val );
    fprintf(stdout, " %f\n", (double)parms->l_location * sse_loc );
    fprintf(stdout, " %f\n", (double)parms->l_lap * sse_lap );
    fprintf(stdout, " %f\n", (double)parms->l_ashburner_triangle * sse_ashburner_triangle );
    fprintf(stdout, " %f\n", (double)parms->l_shrinkwrap * sse_shrinkwrap );
    fprintf(stdout, " %f\n", (double)parms->l_expandwrap * sse_expandwrap );
    fprintf(stdout, " %f\n", (double)parms->l_map * sse_map );
    fprintf(stdout, " %f\n", (double)parms->l_map2d * sse_map2d );
    fprintf(stdout, " %f\n", (double)parms->l_grad * sse_grad );
    fprintf(stdout, " %f\n", (double)parms->l_parea * sse_area );
    fprintf(stdout, " %f\n", (double)parms->l_nldist * sse_nl_dist );
    fprintf(stdout, " %f\n", (double)parms->l_nlarea * sse_nl_area );
    fprintf(stdout, " %f\n", (double)parms->l_angle * sse_angle );
    fprintf(stdout, " %f\n", (double)parms->l_dist * sse_dist );
    fprintf(stdout, " %f\n", (double)parms->l_nlspring * sse_nlspring );
    fprintf(stdout, " %f\n", (double)parms->l_spring * sse_spring );
    fprintf(stdout, " %f\n", (double)parms->l_histo * sse_histo );
    fprintf(stdout, " %f\n", (double)parms->l_dura * sse_dura );
    fprintf(stdout, " %f\n", (double)parms->l_tspring * sse_tspring );
    fprintf(stdout, " %f\n", (double)l_corr * sse_corr );
    fprintf(stdout, " %f\n", (double)parms->l_curv * CURV_SCALE * sse_curv);

    #define SEP
    double sum = 0;
    #define ELT(NAME, MULTIPLIER, COND, EXPR) fprintf(stdout, "old %s : %f \n", #NAME, (MULTIPLIER) * (NAME));  sum += (MULTIPLIER) * (NAME);
    ELT(sse_init, 1, true, sse_init)
    SSE_TERMS
    fprintf(stdout, "old sum = %f,  sse = %f\n", sum, sse);
    #undef ELT
    #undef SEP
  }
  
  if (mht_v_current) {
    MHTfree(&mht_v_current);
  }
  if (mht_f_current) {
    MHTfree(&mht_f_current);
  }

  if (!devFinite(sse)) {
    DiagBreak();
  }
  return (sse);
}

#undef COMPUTE_DISTANCE_ERROR
#undef SSE_TERMS
