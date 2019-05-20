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
//

#include "mrisurf_metricProperties.h"

void mrisComputeFaceRelevantAngleAndArea   (MRIS    *mris, INTEGRATION_PARMS *parms, double* relevant_angle, double* computed_neg_area, double* computed_area);
void mrismp_ComputeFaceRelevantAngleAndArea(MRIS_MP *mris, INTEGRATION_PARMS *parms, double* relevant_angle, double* computed_neg_area, double* computed_area);

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
