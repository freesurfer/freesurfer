#pragma once
/*
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

void mrisDxyzSetLocationMoveLen(double newval);
int mrisComputeAngleAreaTerms(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeAshburnerTriangleTerm(MRIS *mris, double l_ashburner_triangle, INTEGRATION_PARMS *parms);
int mrisComputeBorderTerm(MRIS *mris, double l_border);
int mrisComputeConvexityTerm(MRIS *mris, double l_convex);
int mrisComputeCorrelationTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeDistanceTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeExpandwrapTerm(MRIS *mris, MRI *mri_brain, double l_expandwrap);
int mrisComputeExpansionTerm(MRIS *mris, double l_expand);
int mrisComputeHistoTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeIntensityGradientTerm(MRIS *mris, double l_grad, MRI *mri_brain, MRI *mri_smooth);
int mrisComputeIntensityTerm(MRIS *mris, double l_intensity, MRI *mri_brain, MRI *mri_smooth, double sigma, INTEGRATION_PARMS *parms);
int mrisComputeLaplacianTerm(MRIS *mris, double l_lap);
int mrisComputeLinkTerm(MRIS *mris, double l_link, int pial);
int mrisComputeMaxSpringTerm(MRIS *mris, double l_max_spring);
int mrisComputeNonlinearAreaTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeNonlinearDistanceTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeNonlinearSpringTerm(MRIS *mris, double l_nlspring, INTEGRATION_PARMS *parms);
int mrisComputeNonlinearTangentialSpringTerm(MRIS *mris, double l_spring, double min_dist);
int mrisComputeNormalizedSpringTerm(MRIS *mris, double const l_spring);
int mrisComputePlaneTerm(MRIS *mris, double l_plane, double l_spacing);
int mrisComputePolarCorrelationTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputePolarVectorCorrelationTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputePosterior2DTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputePosteriorTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeQuadraticCurvatureTerm(MRIS *mris, double l_curv);
int mrisComputeRepulsiveRatioTerm(MRIS *mris, double l_repulse, MHT *mht_v);
int mrisComputeRepulsiveTerm(MRIS *mris, double l_repulse, MHT *mht_v, MHT *mht_f);
int mrisComputeShrinkwrapTerm(MRIS *mris, MRI *mri_brain, double l_shrinkwrap);
int mrisComputeSphereTerm(MRIS *mris, double l_sphere, float radius, int    explode_flag);
int mrisComputeSurfaceNormalIntersectionTerm(MRIS *mris, MHT *mht, double l_norm, double max_dist);
int mrisComputeSurfaceRepulsionTerm(MRIS *mris, double l_repulse, MHT *mht);
int mrisComputeTargetLocationTerm(MRIS *mris, double l_location, INTEGRATION_PARMS *parms);
int mrisComputeThicknessMinimizationTerm(MRIS *mris, double l_thick_min, INTEGRATION_PARMS *parms);
int mrisComputeThicknessNormalTerm(MRIS *mris, double l_thick_normal, INTEGRATION_PARMS *parms);
int mrisComputeThicknessParallelTerm(MRIS *mris, double l_thick_parallel, INTEGRATION_PARMS *parms);
int mrisComputeThicknessSmoothnessTerm(MRIS *mris, double l_tsmooth, INTEGRATION_PARMS *parms);
int mrisComputeThicknessSpringTerm(MRIS *mris, double l_thick_spring, INTEGRATION_PARMS *parms);
int mrisComputeVectorCorrelationTerm(MRIS *mris, INTEGRATION_PARMS *parms);
int mrisComputeWhichSurfaceRepulsionTerm(MRIS *mris, double l_repulse, MHT *mht, int which, float max_dot);

// normal spring
void mrisComputeNormalSpringTerm(MRIS *mris, double l_spring);
void vertexComputeNormalSpringTerm(MRIS* mris, int vno, float* dx, float* dy, float* dz, double l_spring = 1.0);

// tangential spring
void mrisComputeTangentialSpringTerm(MRIS *mris, double l_spring);
void vertexComputeTangentialSpringTerm(MRIS* mris, int vno, float* dx, float* dy, float* dz, double l_spring = 1.0);
