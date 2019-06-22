#pragma once
/*
 * @file The SSE terms
 *
 */
/*
 * Author: Bevin R Brett
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


// misc
//
static const unsigned int   mrisurf_sse_MAX_VOXELS = 1500;
static const unsigned int   mrisurf_sse_MAX_DISPLACEMENT = 5;
static const double         mrisurf_sse_DISPLACEMENT_DELTA = 0.1;
static const float          mrisurf_sse_DEFAULT_STD = 4.0f;

double vlst_loglikelihood(
    MRI_SURFACE *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM *hin, HISTOGRAM *hout);

double vlst_loglikelihood2D(
    MRI_SURFACE *mris, MRI *mri, int vno, double displacement, VOXEL_LIST *vl, HISTOGRAM2D *h, FILE *fp);

int mrisCreateLikelihoodHistograms(MRI_SURFACE *mris, INTEGRATION_PARMS *parms);

void mrisComputeFaceRelevantAngleAndArea   (MRIS*    mris, INTEGRATION_PARMS *parms, double* relevant_angle, double* computed_neg_area, double* computed_area);
void mrismp_ComputeFaceRelevantAngleAndArea(MRIS_MP* mris, INTEGRATION_PARMS *parms, double* relevant_angle, double* computed_neg_area, double* computed_area);
int mrisAverageSignedGradients             (MRIS*    mris, int num_avgs);
int mrisComputePositioningGradients        (MRIS*    mris, INTEGRATION_PARMS *parms);


// error measurements
//
double mrisComputeCorrelationError              (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds);
double mrisComputeCorrelationErrorTraceable     (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds, bool trace);
double mrisComputeDistanceError                 (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeDuraError                     (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeExpandwrapError               (MRI_SURFACE *mris,                             MRI *mri_brain, double l_expandwrap, double target_radius);
double mrisComputeIntensityError                (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeIntensityGradientError        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeShrinkwrapError               (MRI_SURFACE *mris,                             MRI *mri_brain, double l_shrinkwrap);
double mrisComputeSphereError                   (MRI_SURFACE *mris, double l_sphere, double a);
double mrisComputeTargetLocationError           (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeVectorCorrelationError        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int use_stds);

double mrisComputeError                         (MRI_SURFACE *mris, INTEGRATION_PARMS *parms,
                                                                                            float *parea_rms,
                                                                                            float *pangle_rms,
                                                                                            float *pcurv_rms,
                                                                                            float *pdist_rms,
                                                                                            float *pcorr_rms);
double mrisRmsDistanceError                     (MRI_SURFACE *mris);

double mrisComputeHistoNegativeLikelihood       (MRI_SURFACE *mris, INTEGRATION_PARMS *parms);
double mrisComputeNegativeLogPosterior          (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox);
double mrisComputeNegativeLogPosterior2D        (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, int *pnvox);


// MEF support
//
double mrisRmsValError_mef                      (MRI_SURFACE *mris, MRI *mri_30, MRI *mri_5, float weight30, float weight5);
double mrisComputeSSE_MEF                       (MRI_SURFACE *mris, INTEGRATION_PARMS *parms, MRI *mri30, MRI *mri5, double weight30, double weight5, MHT *mht);

int    mrisComputeIntensityTerm_mef             (MRI_SURFACE *mris,
                                                                    double l_intensity,
                                                                    MRI *mri_30,
                                                                    MRI *mri_5,
                                                                    double sigma_global,
                                                                    float weight30,
                                                                    float weight5,
                                                                    INTEGRATION_PARMS *parms);
int mrisComputeIntensityTerm_mef                (MRI_SURFACE *mris,
                                                                    double l_intensity,
                                                                    MRI *mri_30,
                                                                    MRI *mri_5,
                                                                    double sigma_global,
                                                                    float weight30,
                                                                    float weight5,
                                                                    INTEGRATION_PARMS *parms);

// energy measurements
//
double mrisComputeAshburnerTriangleEnergy       (MRI_SURFACE *mris, double l_ashburner_triangle, INTEGRATION_PARMS *parms);
double mrisComputeLaplacianEnergy               (MRI_SURFACE *mris);
double mrisComputeNonlinearSpringEnergy         (MRI_SURFACE *mris,                              INTEGRATION_PARMS *parms);
double mrisComputeRepulsiveEnergy               (MRI_SURFACE *mris, double l_repulse, MHT *mht_v_current, MHT *mht_f_current);
double mrisComputeRepulsiveRatioEnergy          (MRI_SURFACE *mris, double l_repulse);
double mrisComputeSpringEnergy                  (MRI_SURFACE *mris);
double mrisComputeSurfaceRepulsionEnergy        (MRI_SURFACE *mris, double l_repulse, MHT *mht);
double mrisComputeTangentialSpringEnergy        (MRI_SURFACE *mris);
double mrisComputeThicknessMinimizationEnergy   (MRI_SURFACE *mris, double l_thick_min,          INTEGRATION_PARMS *parms);
double mrisComputeThicknessNormalEnergy         (MRI_SURFACE *mris, double l_thick_normal,       INTEGRATION_PARMS *parms);
double mrisComputeThicknessParallelEnergy       (MRI_SURFACE *mris, double l_thick_parallel,     INTEGRATION_PARMS *parms);
double mrisComputeThicknessSmoothnessEnergy     (MRI_SURFACE *mris, double l_repulse,            INTEGRATION_PARMS *parms);
double mrisComputeThicknessSpringEnergy         (MRI_SURFACE *mris, double l_thick_spring,       INTEGRATION_PARMS *parms);

#define LIST_OF_PER_VERTEX_SSETERMS \
    ELT(CorrelationError              , (MRIS *mris,                            INTEGRATION_PARMS *parms, int use_stds              )) SEP \
    ELT(RepulsiveRatioEnergy          , (MRIS *mris, double l_repulse                                                               )) SEP \
    ELT(SpringEnergy                  , (MRIS *mris                                                                                 )) SEP \
    ELT(ThicknessMinimizationEnergy   , (MRIS *mris, double l_thick_min,        INTEGRATION_PARMS *parms                            )) SEP \
    ELT(ThicknessNormalEnergy         , (MRIS *mris, double l_thick_normal,     INTEGRATION_PARMS *parms                            )) SEP \
    ELT(ThicknessSpringEnergy         , (MRIS *mris, double l_thick_spring,     INTEGRATION_PARMS *parms                            )) SEP \
    ELT(ThicknessParallelEnergy       , (MRIS *mris, double l_thick_parallel,   INTEGRATION_PARMS *parms                            )) SEP \
    ELT(ThicknessSmoothnessEnergy     , (MRIS *mris, double l_tsmooth,          INTEGRATION_PARMS *parms                            )) SEP \
    ELT(DistanceError                 , (MRIS *mris,                            INTEGRATION_PARMS *parms                            ))     \
    // end of macro
    
#define LIST_OF_PER_FACE_SSETERMS \
    ELT(NonlinearAreaSSE              , (MRIS *mris                                                                                 )) \
    // end of macro

#define LIST_OF_SSETERMS LIST_OF_PER_VERTEX_SSETERMS SEP LIST_OF_PER_FACE_SSETERMS


#define SEP 
#define ELT(NAME, SIGNATURE)    double mrisCompute##NAME SIGNATURE;
LIST_OF_SSETERMS
#undef ELT
#undef SEP

#define SEP 
#define MRIS MRIS_MP
#define ELT(NAME, SIGNATURE)    double mrismp_Compute##NAME SIGNATURE;
LIST_OF_SSETERMS
#undef ELT
#undef MRIS
#undef SEP


#define SEP
#define ELT(NAME, SIGNATURE)    MRIS_sseTermFlag_##NAME##_bitPosition,
enum {
    LIST_OF_SSETERMS
    MRIS_sseTermFlag_bitPosition_LAST };
#undef ELT
#undef SEP
#define SEP ,
#define ELT(NAME, SIGNATURE)    MRIS_sseTermFlag_##NAME = (1L << (int)MRIS_sseTermFlag_##NAME##_bitPosition)
typedef enum {
    LIST_OF_SSETERMS } MRIS_sseTermFlags;
#undef ELT
#undef SEP


typedef struct MRIS_sseTerms {
    double term[MRIS_sseTermFlag_bitPosition_LAST];
} MRIS_sseTerms;
