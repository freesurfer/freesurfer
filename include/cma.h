/**
 * @brief constants for neuroanatomical structures.
 *
 * constants and macros for neuroanatomical and some vascular structures.
 * Originally it was just cma labels, but it has been generalized.
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef CMA_H
#define CMA_H

#include "stats.h"
#include "mri2.h"

/*
 * colors for these labels defined in distribution/FreeSurferColorLUT.txt
 */

#define Unknown                        0
#define Left_Cerebral_Exterior         1
#define Left_Cerebral_White_Matter     2
#define Left_Cerebral_Cortex           3
#define Left_Lateral_Ventricle         4
#define Left_Inf_Lat_Vent              5
#define Left_Cerebellum_Exterior       6
#define Left_Cerebellum_White_Matter   7
#define Left_Cerebellum_Cortex         8
// 9 used to be thalamus and 10 thalamus_proper
#define Left_Thalamus_Proper          10
#define Left_Thalamus                 10
#define Left_Caudate                  11
#define Left_Putamen                  12
#define Left_Pallidum                 13
#define Third_Ventricle               14
#define Fourth_Ventricle              15
#define Brain_Stem                    16
#define Left_Hippocampus              17
#define Left_Amygdala                 18
#define Left_Insula                   19
#define Left_Operculum                20
#define Line_1                        21
#define Line_2                        22
#define Line_3                        23
#define CSF                           24
#define Left_Lesion                   25
#define Left_Accumbens_area           26
#define Left_Substancia_Nigra         27
#define Left_VentralDC                28
#define Left_undetermined             29
#define Left_vessel                   30
#define Left_choroid_plexus           31
#define Left_F3orb                    32
#define Left_lOg                      33
#define Left_aOg                      34
#define Left_mOg                      35
#define Left_pOg                      36
#define Left_Stellate                 37
#define Left_Porg                     38
#define Left_Aorg                     39
#define Right_Cerebral_Exterior       40
#define Right_Cerebral_White_Matter   41
#define Right_Cerebral_Cortex         42
#define Right_Lateral_Ventricle       43
#define Right_Inf_Lat_Vent            44
#define Right_Cerebellum_Exterior     45
#define Right_Cerebellum_White_Matter 46
#define Right_Cerebellum_Cortex       47
// 48 used to be thalamus and 49 thalamus_proper
#define Right_Thalamus_Proper         49
#define Right_Thalamus                49
#define Right_Caudate                 50
#define Right_Putamen                 51
#define Right_Pallidum                52
#define Right_Hippocampus             53
#define Right_Amygdala                54
#define Right_Insula                  55
#define Right_Operculum               56
#define Right_Lesion                  57
#define Right_Accumbens_area          58
#define Right_Substancia_Nigra        59
#define Right_VentralDC               60
#define Right_undetermined            61
#define Right_vessel                  62
#define Right_choroid_plexus          63
#define Right_F3orb                   64
#define Right_lOg                     65
#define Right_aOg                     66
#define Right_mOg                     67
#define Right_pOg                     68
#define Right_Stellate                69
#define Right_Porg                    70
#define Right_Aorg                    71
#define Fifth_Ventricle               72
#define Left_Interior                 73
#define Right_Interior                74
//#define Left_Lateral_Ventricles       75
//#define Right_Lateral_Ventricles      76
#define WM_hypointensities            77
#define Left_WM_hypointensities       78
#define Right_WM_hypointensities      79
#define non_WM_hypointensities        80
#define Left_non_WM_hypointensities   81
#define Right_non_WM_hypointensities  82
#define Left_F1                       83
#define Right_F1                      84
#define Optic_Chiasm                  85
#define Left_future_WMSA              87
#define Right_future_WMSA             88
#define future_WMSA                   89
#define non_WM_hypointensities        80
#define Left_Amygdala_Anterior        96
#define Right_Amygdala_Anterior       97

#define Left_Claustrum                138  //                          65  135 20  0
#define Right_Claustrum               139 //                          65  135 20  0

#define IS_FUTURE_WMSA(l) (((l) == Left_future_WMSA) || ((l) == Right_future_WMSA) || ((l) == future_WMSA))

/*
 * no brain labels after this please unless you fix the IS_BRAIN macro
 */

#define Dura              98
#define Epidermis        118
#define Conn_Tissue      119
#define SC_FAT_MUSCLE    120
#define Cranium          121
#define CSF_SA           122
#define Muscle           123
#define Ear              124
#define Fatty_Tissue     125
#define Spinal_Cord      126
#define Soft_Tissue      127
#define Nerve            128
#define Bone             129
#define Air              130
#define Orbit            131
#define Tongue           132
#define Nasal_Structures 133
#define Globe            134
#define Teeth            135

#define Skull            165
#define Right_Temporal_Cerebral_White_Matter 186
#define Left_Temporal_Cerebral_White_Matter  187
#define Fat              189
#define Bright_Unknown   190
#define Dark_Unknown     191
#define Corpus_Callosum  192

#define alveus                    201
#define perforant_pathway         202
#define parasubiculum             203
#define presubiculum              204
#define subiculum                 205
#define CA1                       206
#define CA2                       207
#define CA3                       208
#define CA4                       209
#define GC_DG                     210
#define HATA                      211
#define fimbria                   212
#define lateral_ventricle         213
#define molecular_layer_HP        214
#define hippocampal_fissure       215
#define entorhinal_cortex         216
#define molecular_layer_subiculum 217
#define Amygdala                  218
#define Cerebral_White_Matter     219
#define Cerebral_Cortex           220
#define Inf_Lat_Vent              221

#if 1
#define Left_hippocampal_fissure  193
#define Left_CADG_head            194
#define Left_subiculum            195
#define Left_fimbria              196
#define Right_hippocampal_fissure 197
#define Right_CADG_head           198
#define Right_subiculum           199
#define Right_fimbria             200
#endif

// reserve 230-249 for more CC divisions

#define MIN_CC_EXTRA      230
#define MAX_CC_EXTRA      249

#define Fornix            250
#define CC_Posterior      251
#define CC_Mid_Posterior  252
#define CC_Central        253
#define CC_Mid_Anterior   254
#define CC_Anterior       255

#define IS_CC(l) ((l >= CC_Posterior && l <= CC_Anterior) || (l >= MIN_CC_EXTRA && l <= MAX_CC_EXTRA))

#define VOXEL_UNCHANGED   256
#define CSF_ExtraCerebral 257
#define Head_ExtraCerebral 258

// vascular and lymph labels (from Alex G)
#define Aorta                    331
#define Left_Common_IliacA       332
#define Right_Common_IliacA      333
#define Left_External_IliacA     334
#define Right_External_IliacA    335
#define Left_Internal_IliacA     336
#define Right_Internal_IliacA    337
#define Left_Lateral_SacralA     338
#define Right_Lateral_SacralA    339
#define Left_ObturatorA          340
#define Right_ObturatorA         341
#define Left_Internal_PudendalA  342
#define Right_Internal_PudendalA 343
#define Left_UmbilicalA          344
#define Right_UmbilicalA         345
#define Left_Inf_RectalA         346
#define Right_Inf_RectalA        347
#define Left_Common_IliacV       348
#define Right_Common_IliacV      349
#define Left_External_IliacV     350
#define Right_External_IliacV    351
#define Left_Internal_IliacV     352
#define Right_Internal_IliacV    353
#define Left_ObturatorV          354
#define Right_ObturatorV         355
#define Left_Internal_PudendalV  356
#define Right_Internal_PudendalV 357
#define Pos_Lymph                358
#define Neg_Lymph                359

// Brodmann Areas
#define BA17 400  //                                    206 62  78  0
#define BA18 401  //                                    121 18  134 0
#define BA44 402  //                                    199 58  250 0
#define BA45 403  //                                    1   148 0   0
#define BA4a 404  //                                    221 248 164 0
#define BA4p 405  //                                    231 148 34  0
#define BA6 406  //                                     1   118 14  0
#define BA2 407  //                                     120 118 14  0
#define BAun1 408  //                                   123 186 220 0
#define BAun2 409  //                                   238 13  176 0

// HiRes Hippocampus labeling
#define right_CA2_3 500  //                             17  85  136 0
#define right_alveus 501  //                            119 187 102 0
#define right_CA1 502  //                               204 68  34  0
#define right_fimbria 503  //                           204 0   255 0
#define right_presubiculum 504  //                      221 187 17  0
#define right_hippocampal_fissure 505  //               153 221 238 0
#define right_CA4_DG 506  //                            51  17  17  0
#define right_subiculum 507  //                         0   119 85  0
#define right_fornix 508  //                            20  100 200 0

#define left_CA2_3 550  //                              17  85  137 0
#define left_alveus 551  //                             119 187 103 0
#define left_CA1 552  //                                204 68  35  0
#define left_fimbria 553  //                            204 0   254 0
#define left_presubiculum 554  //                       221 187 16  0
#define left_hippocampal_fissure 555  //                153 221 239 0
#define left_CA4_DG 556  //                             51  17  18  0
#define left_subiculum 557  //                          0   119 86  0
#define left_fornix 558  //                             20  100 201 0

#define Tumor       600 //                              253 253 253 0

// Cerebellar parcellation labels from SUIT (matches FreeSurferColorLUT.txt)
#define Cbm_Left_I_IV     601    // 70  130 180 0
#define Cbm_Right_I_IV    602    // 245 245 245 0
#define Cbm_Left_V        603    // 205 62  78  0
#define Cbm_Right_V       604    // 120 18  134 0
#define Cbm_Left_VI       605    // 196 58  250 0
#define Cbm_Vermis_VI     606    // 0   148 0   0
#define Cbm_Right_VI      607    // 220 248 164 0
#define Cbm_Left_CrusI    608    // 230 148 34  0
#define Cbm_Vermis_CrusI  609    // 0   118 14  0
#define Cbm_Right_CrusI   610    // 0   118 14  0
#define Cbm_Left_CrusII   611    // 122 186 220 0
#define Cbm_Vermis_CrusII 612    // 236 13  176 0
#define Cbm_Right_CrusII  613    // 12  48  255 0
#define Cbm_Left_VIIb     614    // 204 182 142 0
#define Cbm_Vermis_VIIb   615    // 42  204 164 0
#define Cbm_Right_VIIb    616    // 119 159 176 0
#define Cbm_Left_VIIIa    617    // 220 216 20  0
#define Cbm_Vermis_VIIIa  618    // 103 255 255 0
#define Cbm_Right_VIIIa   619    // 80  196 98  0
#define Cbm_Left_VIIIb    620    // 60  58  210 0
#define Cbm_Vermis_VIIIb  621    // 60  58  210 0
#define Cbm_Right_VIIIb   622    // 60  58  210 0
#define Cbm_Left_IX       623    // 60  58  210 0
#define Cbm_Vermis_IX     624    // 60  60  60  0
#define Cbm_Right_IX      625    // 255 165 0   0
#define Cbm_Left_X        626    // 255 165 0   0
#define Cbm_Vermis_X      627    // 0   255 127 0
#define Cbm_Right_X       628    // 165 42  42  0


#define SUSPICIOUS 999 //                               255 100 100 0


// Tracula labeling
#define fmajor 5100 //                                  204 102 102 0
#define fminor 5101 //                                  204 102 102 0
#define lh_atr 5102 //                                  255 255 102 0
#define lh_cab 5103 //                                  153 204 0   0
#define lh_ccg 5104 //                                  0   153 153 0
#define lh_cst 5105 //                                  204 153 255 0
#define lh_ilf 5106 //                                  255 153 51  0
#define lh_slfp 5107 //                                 204 204 204 0
#define lh_slft 5108 //                                 153 255 255 0
#define lh_unc 5109 //                                  102 153 255 0
#define rh_atr 5110 //                                  255 255 102 0
#define rh_cab 5111 //                                  153 204 0   0
#define rh_ccg 5112 //                                  0   153 153 0
#define rh_cst 5113 //                                  204 153 255 0
#define rh_ilf 5114 //                                  255 153 51  0
#define rh_slfp 5115 //                                 204 204 204 0
#define rh_slft 5116 //                                 153 255 255 0
#define rh_unc 5117 //                                  102 153 255 0
#define lh_ifof 5118 //                                 153 255 255 0
#define rh_ifof 5119 //                                 153 255 255 0
#define lh_fornix 5120 //                               204 102 153 0
#define rh_fornix 5121 //                               204 102 153 0


/*
# Below is the color table for the cortical labels of the seg volume
# created by mri_aparc2aseg in which the aseg cortex label is replaced
# by the labels in the aparc. It also supports wm labels that will
# eventually be created by mri_aparc2aseg. Otherwise, the aseg labels
# do not change from above. The cortical lables are the same as in
# colortable_desikan_killiany.txt, except that left hemisphere has
# 1000 added to the index and the right has 2000 added.  The label
# names are also prepended with ctx-lh or ctx-rh. The white matter
# labels are the same as in colortable_desikan_killiany.txt, except
# that left hemisphere has 3000 added to the index and the right has
# 4000 added. The label names are also prepended with wm-lh or wm-rh.
# Centrum semiovale is also labled with 5001 (left) and 5002 (right).
# Even further below are the color tables for aparc.a2005s and aparc.a2009s.

*/

#define wm_lh_unknown        3000
#define wm_rh_unknown        4000
#define Left_Unsegmented_WM  5001
#define Right_Unsegmented_WM 5002
#define MIN_CORTICAL_PARCELLATION   1000

#define    ctx_lh_unknown  1000 //                      25  5   25  0
#define    ctx_lh_bankssts  1001 //                     25  100 40  0
#define    ctx_lh_caudalanteriorcingulate  1002 //      125 100 160 0
#define    ctx_lh_caudalmiddlefrontal  1003 //          100 25  0   0
#define    ctx_lh_corpuscallosum  1004 //               120 70  50  0
#define    ctx_lh_cuneus  1005 //                       220 20  100 0
#define    ctx_lh_entorhinal  1006 //                   220 20  10  0
#define    ctx_lh_fusiform  1007 //                     180 220 140 0
#define    ctx_lh_inferiorparietal  1008 //             220 60  220 0
#define    ctx_lh_inferiortemporal  1009 //             180 40  120 0
#define    ctx_lh_isthmuscingulate  1010 //             140 20  140 0
#define    ctx_lh_lateraloccipital  1011 //             20  30  140 0
#define    ctx_lh_lateralorbitofrontal  1012 //         35  75  50  0
#define    ctx_lh_lingual  1013 //                      225 140 140 0
#define    ctx_lh_medialorbitofrontal  1014 //          200 35  75  0
#define    ctx_lh_middletemporal  1015 //               160 100 50  0
#define    ctx_lh_parahippocampal  1016 //              20  220 60  0
#define    ctx_lh_paracentral  1017 //                  60  220 60  0
#define    ctx_lh_parsopercularis  1018 //              220 180 140 0
#define    ctx_lh_parsorbitalis  1019 //                20  100 50  0
#define    ctx_lh_parstriangularis  1020 //             220 60  20  0
#define    ctx_lh_pericalcarine  1021 //                120 100 60  0
#define    ctx_lh_postcentral  1022 //                  220 20  20  0
#define    ctx_lh_posteriorcingulate  1023 //           220 180 220 0
#define    ctx_lh_precentral  1024 //                   60  20  220 0
#define    ctx_lh_precuneus  1025 //                    160 140 180 0
#define    ctx_lh_rostralanteriorcingulate  1026 //     80  20  140 0
#define    ctx_lh_rostralmiddlefrontal  1027 //         75  50  125 0
#define    ctx_lh_superiorfrontal  1028 //              20  220 160 0
#define    ctx_lh_superiorparietal  1029 //             20  180 140 0
#define    ctx_lh_superiortemporal  1030 //             140 220 220 0
#define    ctx_lh_supramarginal  1031 //                80  160 20  0
#define    ctx_lh_frontalpole  1032 //                  100 0   100 0
#define    ctx_lh_temporalpole  1033 //                 70  70  70  0
#define    ctx_lh_transversetemporal  1034 //           150 150 200 0
#define    ctx_lh_insula  1035 //                       255 192 32  0

#define    ctx_rh_unknown  2000 //                      25  5   25  0
#define    ctx_rh_bankssts  2001 //                     25  100 40  0
#define    ctx_rh_caudalanteriorcingulate  2002 //      125 100 160 0
#define    ctx_rh_caudalmiddlefrontal  2003 //          100 25  0   0
#define    ctx_rh_corpuscallosum  2004 //               120 70  50  0
#define    ctx_rh_cuneus  2005 //                       220 20  100 0
#define    ctx_rh_entorhinal  2006 //                   220 20  10  0
#define    ctx_rh_fusiform  2007 //                     180 220 140 0
#define    ctx_rh_inferiorparietal  2008 //             220 60  220 0
#define    ctx_rh_inferiortemporal  2009 //             180 40  120 0
#define    ctx_rh_isthmuscingulate  2010 //             140 20  140 0
#define    ctx_rh_lateraloccipital  2011 //             20  30  140 0
#define    ctx_rh_lateralorbitofrontal  2012 //         35  75  50  0
#define    ctx_rh_lingual  2013 //                      225 140 140 0
#define    ctx_rh_medialorbitofrontal  2014 //          200 35  75  0
#define    ctx_rh_middletemporal  2015 //               160 100 50  0
#define    ctx_rh_parahippocampal  2016 //              20  220 60  0
#define    ctx_rh_paracentral  2017 //                  60  220 60  0
#define    ctx_rh_parsopercularis  2018 //              220 180 140 0
#define    ctx_rh_parsorbitalis  2019 //                20  100 50  0
#define    ctx_rh_parstriangularis  2020 //             220 60  20  0
#define    ctx_rh_pericalcarine  2021 //                120 100 60  0
#define    ctx_rh_postcentral  2022 //                  220 20  20  0
#define    ctx_rh_posteriorcingulate  2023 //           220 180 220 0
#define    ctx_rh_precentral  2024 //                   60  20  220 0
#define    ctx_rh_precuneus  2025 //                    160 140 180 0
#define    ctx_rh_rostralanteriorcingulate  2026 //     80  20  140 0
#define    ctx_rh_rostralmiddlefrontal  2027 //         75  50  125 0
#define    ctx_rh_superiorfrontal  2028 //              20  220 160 0
#define    ctx_rh_superiorparietal  2029 //             20  180 140 0
#define    ctx_rh_superiortemporal  2030 //             140 220 220 0
#define    ctx_rh_supramarginal  2031 //                80  160 20  0
#define    ctx_rh_frontalpole  2032 //                  100 0   100 0
#define    ctx_rh_temporalpole  2033 //                 70  70  70  0
#define    ctx_rh_transversetemporal  2034 //           150 150 200 0
#define    ctx_rh_insula  2035 //                       255 192 32  0

// be sure to update MAX_LABEL if additional labels are added!

#define MAX_LABEL rh_unc
#define MAX_CMA_LABEL (MAX_LABEL)
#define MAX_CMA_LABELS (MAX_CMA_LABEL+1)


#define IS_UNKNOWN(label)  (((label) == Unknown) || (label == 255) || (label == Bright_Unknown) || (label == Dark_Unknown))

#define IS_BRAIN(label)  ((!IS_UNKNOWN(label) && label < Dura) || IS_CC(label))

// (label == Brain_Stem) ||
#define IS_WM(label) (((label) == Left_Cerebral_White_Matter) || ((label) == Right_Cerebral_White_Matter) || ((label) == Left_Temporal_Cerebral_White_Matter) || ((label) == Right_Temporal_Cerebral_White_Matter) || IS_CC(label) || (label == Left_VentralDC) || (label == Right_VentralDC) || IS_CEREBELLAR_WM(label))
#define IS_HYPO(label) (((label) == WM_hypointensities)  || ((label) == Left_WM_hypointensities)  || ((label) == Right_WM_hypointensities) || IS_FUTURE_WMSA(label))
#define IS_WMSA(label) IS_HYPO(label)
#define IS_WMH(label) (IS_WM(label) || IS_HYPO(label))
#define IS_THALAMUS(label)  (((label) == Left_Thalamus) || ((label) == Left_Thalamus) || ((label) == Right_Thalamus) || ((label) == Right_Thalamus))
#define IS_GM(label) (((label) == Left_Cerebral_Cortex) || ((label) == Right_Cerebral_Cortex))
#define IS_VENTRAL_DC(l)  (((l) == Left_VentralDC) || ((l) == Right_VentralDC))



#define IS_CEREBELLAR_WM(label) (((label) == Left_Cerebellum_White_Matter) || ((label) == Right_Cerebellum_White_Matter))
#define IS_CEREBELLAR_GM(label) (((label) == Left_Cerebellum_Cortex) || ((label) == Right_Cerebellum_Cortex))
#define IS_CEREBELLUM(label) (IS_CEREBELLAR_WM(label) || IS_CEREBELLAR_GM(label))

#define IS_HIPPO(l) (((l) == Left_Hippocampus) || ((l) == Right_Hippocampus))
#define IS_AMYGDALA(l) (((l) == Left_Amygdala) || ((l) == Right_Amygdala))
#define IS_MTL(l)   (IS_HIPPO(l) || IS_AMYGDALA(l))
#define IS_CORTEX(l) (((l) == Left_Cerebral_Cortex) || \
                      ((l) == Right_Cerebral_Cortex))
#define IS_LAT_VENT(l) (((l) == Left_Lateral_Ventricle) || \
                        ((l) == Right_Lateral_Ventricle) || \
                        ((l) == Right_Inf_Lat_Vent) || \
                        ((l) == Left_Inf_Lat_Vent))
#define IS_CSF(l) (IS_LAT_VENT(l) || ((l) == CSF) || ((l) == CSF_SA) || ((l) == Third_Ventricle) || ((l) == Fourth_Ventricle))

#define IS_INF_LAT_VENT(l)  (((l) == Left_Inf_Lat_Vent) || ((l) == Right_Inf_Lat_Vent))
#define IS_VENTRICLE(l)  (IS_LAT_VENT(l) || IS_INF_LAT_VENT(l) || ((l) == Third_Ventricle) || ((l) == Fourth_Ventricle))
#define IS_CAUDATE(l) (((l) == Left_Caudate) || ((l) == Right_Caudate))
#define IS_PUTAMEN(l) (((l) == Left_Putamen) || ((l) == Right_Putamen))
#define IS_CLAUSTRUM(l) (((l) == Left_Claustrum) || ((l) == Right_Claustrum))
#define IS_ACCUMBENS(l) (((l) == Left_Accumbens_area) || ((l) == Right_Accumbens_area))
#define IS_PALLIDUM(l) (((l) == Left_Pallidum) || ((l) == Right_Pallidum))

#define LABEL_WITH_NO_TOPOLOGY_CONSTRAINT(l) (\
   ((l) == Right_non_WM_hypointensities) || \
   ((l) == Left_non_WM_hypointensities) || \
   ((l) == Left_WM_hypointensities) || \
   ((l) == Right_WM_hypointensities) || \
   ((l) == Left_choroid_plexus) || \
   ((l) == Right_choroid_plexus) || \
   ((l) == WM_hypointensities) || \
   ((l) == Dura) || \
   ((l) == Bone) || \
   ((l) == CSF_SA) || \
   ((l) == Epidermis) || \
   ((l) == SC_FAT_MUSCLE) || \
   ((l) == non_WM_hypointensities))

/* --- below: see cma.c --- */

#define MAX_OUTLINE_CLAIMS 10

#define CMA_FILL_NONE     0
#define CMA_FILL_OUTLINE  1
#define CMA_FILL_INTERIOR 2

typedef struct
{
  int n_claims;
  int interior_claim_flag;
  short claim_labels[MAX_OUTLINE_CLAIMS];
  float claim_values[MAX_OUTLINE_CLAIMS];
  float no_label_claim;
}
CMAoutlineClaim;

typedef struct
{
  int width, height;
  CMAoutlineClaim **claim_field;
  //unsigned char **fill_field;
  short **fill_field;
  unsigned char **outline_points_field;
}
CMAoutlineField;

CMAoutlineField *CMAoutlineFieldAlloc(int width, int height);
int CMAfreeOutlineField(CMAoutlineField **of);

int CMAclearFillField(CMAoutlineField *field);
int CMAfill(CMAoutlineField *field, short seed_x, short seed_y);

int CMAclaimPoints(CMAoutlineField *field,
                   short label, short *points,
                   int n_points, short seed_x, short seed_y);
int CMAassignLabels(CMAoutlineField *field);

int CMAvalueClaims(CMAoutlineClaim *claim);
int CMAvalueAllClaims(CMAoutlineField *field);

/* the weights to give claims of nearby points -- see cma.c */
#define OTL_CLAIM_WEIGHT_CENTER   1.000
#define OTL_CLAIM_WEIGHT_SIDE     0.500
#define OTL_CLAIM_WEIGHT_DIAGONAL 0.707

short CMAtotalClaims(CMAoutlineField *field, int x, int y);
int CMAaddWeightedTotals(CMAoutlineClaim *claim,
                         float weight, float *claim_totals);

int CMAzeroOutlines(CMAoutlineField *field);
const char *cma_label_to_name(int label) ;
int IsSubCorticalGray(int SegId);
#include "mri.h"
double SupraTentorialVolCorrection(MRI *aseg, MRI *ribbon);
double CorticalGMVolCorrection(MRI *aseg, MRI *ribbon, int hemi);
MRI *MRIfixAsegWithRibbon(MRI *aseg, MRI *ribbon, MRI *asegfixed);

std::vector<double> ComputeBrainVolumeStats(const std::string& subject, const std::string& subjdir);
std::vector<double> ComputeBrainVolumeStats2(const std::string& subject, const std::string& subjdir, const int KeepCSF);
void CacheBrainVolumeStats(const std::vector<double>& stats, const std::string& subject, const std::string& subjdir);
std::vector<double> ReadCachedBrainVolumeStats(const std::string& subject, const std::string& subjdir);

#define IS_FIMBRIA(l) ((l) == left_fimbria || (l) == right_fimbria || (l) == fimbria)
#define CSF_CLASS        0
#define GM_CLASS         1
#define WM_CLASS         2
#define NTISSUE_CLASSES  3
#define LH_CLASS         4
#define RH_CLASS         5
#define FLUID_CLASS      6
#define OTHER_CLASS      7

#define IS_GRAY_MATTER(l) (         \
   ((l) == Left_Cerebral_Cortex) || \
   ((l) == Left_Hippocampus) || \
   ((l) == Left_Amygdala) || \
   ((l) == Left_Putamen) || \
   ((l) == Left_Caudate) || \
   ((l) == Left_Cerebellum_Cortex) || \
   ((l) == Left_Accumbens_area) || \
   ((l) == Right_Cerebral_Cortex) || \
   ((l) == Right_Hippocampus) || \
   ((l) == Right_Amygdala) || \
   ((l) == Right_Putamen) || \
   ((l) == Right_Caudate) || \
   ((l) == Right_Cerebellum_Cortex) || \
   ((l) == Right_Accumbens_area))

#define IS_WHITE_MATTER(l) (         \
   ((l) == Left_Cerebellum_White_Matter) || \
   ((l) == Left_Cerebral_White_Matter) || \
   ((l) == Right_Cerebellum_White_Matter) || \
   ((l) == Right_Cerebral_White_Matter))

#define IS_FLUID(l) (         \
   ((l) == Left_Lateral_Ventricle) || \
   ((l) == Left_Inf_Lat_Vent) || \
   ((l) == Right_Lateral_Ventricle) || \
   ((l) == Right_Inf_Lat_Vent) || \
   ((l) == Unknown))


#define IS_LH_CLASS(l) (\
   ((l) == Left_Cerebral_Cortex) || \
   ((l) == Left_Cerebral_White_Matter) || \
   ((l) == Left_Hippocampus) || \
   ((l) == Left_Amygdala) || \
   ((l) == Left_Putamen) || \
   ((l) == Left_Pallidum) || \
   ((l) == Left_Caudate) || \
   ((l) == Left_Lateral_Ventricle) || \
   ((l) == Left_Inf_Lat_Vent) || \
   ((l) == Left_Cerebellum_Cortex) || \
   ((l) == Left_Cerebellum_White_Matter) || \
   ((l) == Left_Thalamus) || \
   ((l) == Left_vessel) || \
   ((l) == Left_choroid_plexus) || \
   ((l) == Left_VentralDC) || \
   ((l) == Left_Accumbens_area))

#define IS_RH_CLASS(l) (\
   ((l) == Right_Cerebral_Cortex) || \
   ((l) == Right_Cerebral_White_Matter) || \
   ((l) == Right_Hippocampus) || \
   ((l) == Right_Amygdala) || \
   ((l) == Right_Putamen) || \
   ((l) == Right_Pallidum) || \
   ((l) == Right_Caudate) || \
   ((l) == Right_Lateral_Ventricle) || \
   ((l) == Right_Inf_Lat_Vent) || \
   ((l) == Right_Cerebellum_Cortex) || \
   ((l) == Right_Cerebellum_White_Matter) || \
   ((l) == Right_Thalamus) || \
   ((l) == Right_vessel) || \
   ((l) == Right_choroid_plexus) || \
   ((l) == Right_VentralDC) || \
   ((l) == Right_Accumbens_area))


// only structures that are 'pure' of one tissue type
#define IS_GRAY_CLASS(l) (IS_GM(l) || ((l) == Left_Caudate) || ((l) == Right_Caudate))

#define IS_WHITE_CLASS(l) (((l) == Left_Cerebral_White_Matter) || ((l) == Right_Cerebral_White_Matter))

#define IS_CSF_CLASS(l) (((l) == Left_Lateral_Ventricle) || ((l) == Right_Lateral_Ventricle) || ((l) == CSF) || ((l) == CSF_SA) || ((l) == Third_Ventricle) || ((l) == Fourth_Ventricle) || ((l) == Fifth_Ventricle) || ((l) == left_hippocampal_fissure) || ((l) == right_hippocampal_fissure) || ((l) == hippocampal_fissure))

#define IS_CLASS(l,c) (c == CSF_CLASS ? IS_CSF_CLASS(l) : c == GM_CLASS ? IS_GRAY_CLASS(l) : IS_WHITE_CLASS(l))

#include "mrisurf.h"
int insert_ribbon_into_aseg(MRI *mri_src_aseg, MRI *mri_aseg,
                            MRI_SURFACE *mris_white, MRI_SURFACE *mris_pial,
                            int hemi) ;
int MRIasegContraLatLabel(int id);
MRI *MRIlrswapAseg(MRI *aseg);
MRI *MRIseg2TissueType(MRI *seg, COLOR_TABLE *ct, MRI *tt);
int CheckSegTissueType(MRI *seg, COLOR_TABLE *ct);
MRI *MRIextractTissueTypeSeg(MRI *seg, COLOR_TABLE *ct, int tt, MRI *ttseg);
MRI **MRIdilateSegWithinTT(MRI *seg, int nDils, COLOR_TABLE *ct, MRI **r);
SEGSTAT *Seg2NbrNonBrain(MRI *seg, COLOR_TABLE *ctab, double threshmm);
int Seg2NbrNonBrainWrapper(char *subject, char *segname, COLOR_TABLE *ctab, char *statname, double threshmm);

#endif
