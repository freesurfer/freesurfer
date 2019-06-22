#pragma once
/*
 * @file The various ways of calculating vertex dx dy dz
 *
 */
/*
 * Author: Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett 
 *
 * $ © copyright-2019 The General Hospital Corporation (Boston, MA) "MGH"
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
 
//  Organize the ComputeSSE terms in a way that enables calculation of several of them in a single pass over either the vertices or the faces
//  and using either the MRIS or another representation of the Surface Face Vertex information
//
#include "mrisurf_metricProperties.h"


int    mrisComputeAngleAreaTerms                (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeAshburnerTriangleTerm         (MRI_SURFACE *mris, double l_ashburner_triangle, INTEGRATION_PARMS *parms);
int    mrisComputeBorderTerm                    (MRI_SURFACE *mris, double l_border);
int    mrisComputeConvexityTerm                 (MRI_SURFACE *mris, double l_convex);
int    mrisComputeCorrelationTerm               (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeDistanceTerm                  (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeExpandwrapTerm                (MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap);
int    mrisComputeExpansionTerm                 (MRI_SURFACE *mris, double l_expand);
int    mrisComputeIntensityGradientTerm         (MRI_SURFACE *mris, double l_grad,      MRI *mri_brain, MRI *mri_smooth);
int    mrisComputeIntensityTerm                 (MRI_SURFACE *mris, double l_intensity, MRI *mri_brain, MRI *mri_smooth, double sigma, INTEGRATION_PARMS *parms);
int    mrisComputeLinkTerm                      (MRI_SURFACE *mris, double l_spring, int    pial);
int    mrisComputeMaxSpringTerm                 (MRI_SURFACE *mris, double l_spring);
int    mrisComputeNonlinearAreaTerm             (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeNonlinearDistanceTerm         (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeNonlinearSpringTerm           (MRI_SURFACE *mris, double l_nlspring,           INTEGRATION_PARMS *parms);
int    mrisComputeNonlinearTangentialSpringTerm (MRI_SURFACE *mris, double l_spring, double min_dist);
int    mrisComputeNormalSpringTerm              (MRI_SURFACE *mris, double l_spring);
int    mrisComputePolarCorrelationTerm          (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputePolarVectorCorrelationTerm    (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeQuadraticCurvatureTerm        (MRI_SURFACE *mris, double l_curv);
int    mrisComputeRepulsiveRatioTerm            (MRI_SURFACE *mris, double l_repulse, MHT *mht_v);
int    mrisComputeRepulsiveTerm                 (MRI_SURFACE *mris, double l_repulse, MHT *mht_v, MHT *mht_f);
int    mrisComputeShrinkwrapTerm                (MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap);
int    mrisComputeSphereTerm                    (MRI_SURFACE *mris, double l_sphere, float radius, int    explode_flag);
int    mrisComputeSurfaceNormalIntersectionTerm (MRI_SURFACE *mris, MHT *mht, double l_norm, double max_dist);
int    mrisComputeSurfaceRepulsionTerm          (MRI_SURFACE *mris, double l_repulse, MHT *mht);
int    mrisComputeTangentialSpringTerm          (MRI_SURFACE *mris, double l_spring);
int    mrisComputeTargetLocationTerm            (MRI_SURFACE *mris, double l_location,           INTEGRATION_PARMS *parms);
int    mrisComputeThicknessMinimizationTerm     (MRI_SURFACE *mris, double l_thick_min,          INTEGRATION_PARMS *parms);
int    mrisComputeThicknessNormalTerm           (MRI_SURFACE *mris, double l_thick_normal,       INTEGRATION_PARMS *parms);
int    mrisComputeThicknessParallelTerm         (MRI_SURFACE *mris, double l_thick_parallel,     INTEGRATION_PARMS *parms);
int    mrisComputeThicknessSmoothnessTerm       (MRI_SURFACE *mris, double l_tsmooth,            INTEGRATION_PARMS *parms);
int    mrisComputeThicknessSpringTerm           (MRI_SURFACE *mris, double l_thick_spring,       INTEGRATION_PARMS *parms);
int    mrisComputeVectorCorrelationTerm         (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
int    mrisComputeWhichSurfaceRepulsionTerm     (MRI_SURFACE *mris, double l_repulse, MHT *mht, int which, float max_dot);
int    mrisComputePlaneTerm                     (MRI_SURFACE *mris, double l_plane, double l_spacing);
int    mrisComputeBorderTerm                    (MRI_SURFACE *mris, double l_border);
int    mrisComputeConvexityTerm                 (MRI_SURFACE *mris, double l_convex);
int    mrisComputeMaxSpringTerm                 (MRI_SURFACE *mris, double l_max_spring);
int    mrisComputeTargetLocationTerm            (MRI_SURFACE *mris, double l_location, INTEGRATION_PARMS *parms);
int mrisComputeShrinkwrapTerm                   (MRI_SURFACE *mris, MRI *mri_brain, double l_shrinkwrap);
int mrisComputeExpandwrapTerm                   (MRI_SURFACE *mris, MRI *mri_brain, double l_expandwrap);
int mrisComputeIntensityGradientTerm            (MRI_SURFACE *mris, double l_grad, MRI *mri_brain, MRI *mri_smooth);
int mrisComputeSurfaceRepulsionTerm             (MRI_SURFACE *mris, double l_repulse, MHT *mht);
int mrisComputeHistoTerm                        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
int mrisComputePosteriorTerm                    (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
int mrisComputePosterior2DTerm                  (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
int mrisComputeWhichSurfaceRepulsionTerm        (MRI_SURFACE *mris, double l_repulse, MHT *mht, int which, float dot_thresh);
int mrisComputeLaplacianTerm                    (MRI_SURFACE *mris, double l_lap);
int mrisComputeNormalizedSpringTerm             (MRI_SURFACE *mris, double const l_spring);
int mrisComputeRepulsiveTerm                    (MRI_SURFACE *mris, double l_repulse, MHT *mht, MHT *mht_faces);
int mrisComputeThicknessSmoothnessTerm          (MRI_SURFACE *mris, double l_tsmooth, INTEGRATION_PARMS *parms);
int mrisComputeThicknessMinimizationTerm        (MRI_SURFACE *mris, double l_thick_min, INTEGRATION_PARMS *parms);
int mrisComputeThicknessParallelTerm            (MRI_SURFACE *mris, double l_thick_parallel, INTEGRATION_PARMS *parms);
int mrisComputeNormalSpringTerm                 (MRI_SURFACE *mris, double l_spring);
int mrisComputeQuadraticCurvatureTerm           (MRI_SURFACE *mris, double l_curv);
int mrisComputeNonlinearSpringTerm              (MRI_SURFACE *mris, double l_nlspring, INTEGRATION_PARMS *parms);
int mrisComputeTangentialSpringTerm             (MRI_SURFACE *mris, double l_spring);
int mrisComputeNonlinearTangentialSpringTerm    (MRI_SURFACE *mris, double l_spring, double min_dist);
int mrisComputeLinkTerm                         (MRI_SURFACE *mris, double l_link, int pial);
int mrisComputeSurfaceNormalIntersectionTerm    (MRI_SURFACE *mris, MHT *mht, double l_norm, double max_dist);
