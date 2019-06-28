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
#define LIST_OF_PER_VERTEX_SSETERMS_NYI \
    ELT(AshburnerTriangleEnergy       , (MRIS_PARAMETER_COMMA double l_ashburner_triangle, INTEGRATION_PARMS *parms                         ), (l_ashburner_triangle,parms)) SEP \
    ELT(LaplacianEnergy               , (MRIS_PARAMETER                                                                                     ), ()) SEP \
    ELT(NonlinearSpringEnergy         , (MRIS_PARAMETER_COMMA                              INTEGRATION_PARMS *parms                         ), (parms)) SEP \
    ELT(RepulsiveEnergy               , (MRIS_PARAMETER_COMMA double l_repulse, MHT *mht_v_current, MHT *mht_f_current                      ), (l_repulse, mht_v_current, mht_f_current)) SEP \
    ELT(SurfaceRepulsionEnergy        , (MRIS_PARAMETER_COMMA double l_repulse, MHT *mht                                                    ), (l_repulse, mht)) SEP \
    ELT(TangentialSpringEnergy        , (MRIS_PARAMETER                                                                                     ), ())     \
    // end of macro

#define LIST_OF_PER_VERTEX_SSETERMS_Implemented \
    ELT(RepulsiveRatioEnergy          , (MRIS_PARAMETER_COMMA double l_repulse                                                              ), (l_repulse)) SEP \
    ELT(SpringEnergy                  , (MRIS_PARAMETER                                                                                     ), ()) SEP \
    ELT(ThicknessMinimizationEnergy   , (MRIS_PARAMETER_COMMA double l_thick_min,        INTEGRATION_PARMS *parms                           ), (l_thick_min,parms)) SEP \
    ELT(ThicknessNormalEnergy         , (MRIS_PARAMETER_COMMA double l_thick_normal,     INTEGRATION_PARMS *parms                           ), (l_thick_normal,parms)) SEP \
    ELT(ThicknessSpringEnergy         , (MRIS_PARAMETER_COMMA double l_thick_spring,     INTEGRATION_PARMS *parms                           ), (l_thick_spring,parms)) SEP \
    ELT(ThicknessParallelEnergy       , (MRIS_PARAMETER_COMMA double l_thick_parallel,   INTEGRATION_PARMS *parms                           ), (l_thick_parallel,parms)) SEP \
    ELT(ThicknessSmoothnessEnergy     , (MRIS_PARAMETER_COMMA double l_tsmooth,          INTEGRATION_PARMS *parms                           ), (l_tsmooth,parms)) SEP \
    /**/ \
    ELT(NonlinearDistanceSSE          , (MRIS_PARAMETER                                                                                     ), ()) SEP \
    ELT(QuadraticCurvatureSSE         , (MRIS_PARAMETER_COMMA double l_curv                                                                 ), (l_curv)) SEP \
    ELT(HistoNegativeLikelihood       , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(NegativeLogPosterior          , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms, int *pnvox               ), (parms,pnvox)) SEP \
    ELT(NegativeLogPosterior2D        , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms, int *pnvox               ), (parms,pnvox)) SEP \
    /**/ \
    ELT(CorrelationError              , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms, int use_stds             ), (parms,use_stds)) SEP \
    ELT(DistanceError                 , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(DuraError                     , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(IntensityError                , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(TargetLocationError           , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(IntensityGradientError        , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms                           ), (parms)) SEP \
    ELT(VectorCorrelationError        , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms, int use_stds             ), (parms,use_stds)) SEP \
    ELT(ExpandwrapError               , (MRIS_PARAMETER_COMMA                            MRI *mri_brain, double l_expandwrap, double target_radius), (mri_brain,l_expandwrap,target_radius)) SEP \
    ELT(ShrinkwrapError               , (MRIS_PARAMETER_COMMA                            MRI *mri_brain, double l_shrinkwrap                ), (mri_brain,l_shrinkwrap)) SEP \
    ELT(SphereError                   , (MRIS_PARAMETER_COMMA double l_sphere, double r0                                                    ), (l_sphere,r0)) SEP \
    ELT(RmsDistanceError              , (MRIS_PARAMETER                                                                                     ), ()) SEP \
    ELT(Error                         , (MRIS_PARAMETER_COMMA                            INTEGRATION_PARMS *parms, \
                                                                                            float *parea_rms, \
                                                                                            float *pangle_rms,\
                                                                                            float *pcurv_rms, \
                                                                                            float *pdist_rms, \
                                                                                            float *pcorr_rms                                ), (parms,parea_rms,pangle_rms,pcurv_rms,pdist_rms,pcorr_rms))     \
    // end of macro

#define LIST_OF_PER_VERTEX_SSETERMS LIST_OF_PER_VERTEX_SSETERMS_NYI SEP LIST_OF_PER_VERTEX_SSETERMS_Implemented

#define LIST_OF_PER_FACE_SSETERMS \
    ELT(NonlinearAreaSSE              , (MRIS_PARAMETER                                                                                     ), ())     \
    // end of macro


#define LIST_OF_SSETERMS LIST_OF_PER_VERTEX_SSETERMS SEP LIST_OF_PER_FACE_SSETERMS


#define MRIS_PARAMETER          MRIS* mris
#define MRIS_PARAMETER_COMMA    MRIS_PARAMETER,
#define SEP 
#define ELT(NAME, SIGNATURE, CALL)    double mrisCompute##NAME SIGNATURE;
LIST_OF_SSETERMS
#undef ELT
#undef SEP
#undef MRIS_PARAMETER_COMMA
#undef MRIS_PARAMETER

#define MRIS_PARAMETER          MRIS* mris
#define MRIS_PARAMETER_COMMA    MRIS_PARAMETER,
#define SEP 
#define MRIS MRIS_MP
#define ELT(NAME, SIGNATURE, CALL)    double mrismp_Compute##NAME SIGNATURE;
LIST_OF_PER_VERTEX_SSETERMS_Implemented SEP LIST_OF_PER_FACE_SSETERMS
#undef ELT
#undef MRIS
#undef SEP
#undef MRIS_PARAMETER_COMMA
#undef MRIS_PARAMETER


#define SEP
#define ELT(NAME, SIGNATURE,CALL)    MRIS_sseTermFlag_##NAME##_bitPosition,
enum {
    LIST_OF_SSETERMS
    MRIS_sseTermFlag_bitPosition_LAST };
#undef ELT
#undef SEP
#define SEP ,
#define ELT(NAME, SIGNATURE,CALL)    MRIS_sseTermFlag_##NAME = (1L << (int)MRIS_sseTermFlag_##NAME##_bitPosition)
typedef enum {
    LIST_OF_SSETERMS } MRIS_sseTermFlags;
#undef ELT
#undef SEP


typedef struct MRIS_sseTerms {
    double term[MRIS_sseTermFlag_bitPosition_LAST];
} MRIS_sseTerms;
