#ifndef CMA_H
#define CMA_H

#define    Unknown   0 /*                                 0    0    0     0*/
#define    Left_Cerebral_Exterior  1/*              205   62   78     0 */
#define    Left_Cerebral_White_Matter  2/*          225  225  225     0 */
#define    Left_Cerebral_Cortex  3/*                205   62   78     0 */
#define    Left_Lateral_Ventricle  4/*              120   18  134     0 */
#define    Left_Inf_Lat_Vent  5/*                    196   58  250     0 */
#define    Left_Cerebellum_Exterior  6/*               0  148    0     0 */
#define    Left_Cerebellum_White_Matter  7/*         220  248  164     0 */
#define    Left_Cerebellum_Cortex  8/*               230  148   34     0 */
#define    Left_Thalamus  9/*                          0  118   14     0 */
#define   Left_Thalamus_Proper  10/*                   0  118   14     0 */
#define   Left_Caudate  11/*                         122  186  220    0  */
#define   Left_Putamen  12/*                         236   13  176    0  */
#define   Left_Pallidum  13/*                         12   48  255    0  */
#define   Third_Ventricle  14/*                        204  182  142    0  */
#define   Fourth_Ventricle  15/*                         42  204  164    0  */
#define   Brain_Stem  16/*                           119  159  176    0  */
#define   Left_Hippocampus  17/*                     220  216   20     0 */
#define   Left_Amygdala  18/*                        103  255  255    0  */
#define   Left_Insula  19/*                           80  196   98     0 */
#define   Left_Operculum  20/*                        60   58  210    0  */
#define   Line_1  21/*                                60   58  210    0  */
#define   Line_2  22/*                                60   58  210    0  */
#define   Line_3  23/*                                60   58  210    0  */
#define   CSF  24/*                                   60   60   60     0 */
#define   Left_Lesion  25/*                          255  165    0    0  */
#define   Left_Accumbens_area  26/*                  255  165    0    0  */
#define   Left_Substancia_Nigra  27/*                  0  255  127    0  */
#define   Left_VentralDC  28/*                       165   42   42    0  */
#define   Left_undetermined  29/*                    135  206  235    0  */
#define   Left_vessel  30/*                          160   32  240    0  */
#define   Left_choroid_plexus  31/*                    0  255  255    0  */
#define   Left_F3orb  32/*                           100   50  100    0  */
#define   Left_lOg  33/*                             135   50   74     0 */
#define   Left_aOg  34/*                             122  135   50     0 */
#define   Left_mOg  35/*                              51   50  135    0  */
#define   Left_pOg  36/*                              74  155   60     0 */
#define   Left_Stellate  37/*                        120   62   43    0  */
#define   Left_Porg  38/*                             74  155   60     0 */
#define   Left_Aorg  39/*                            122  135   50      0 */
#define   Right_Cerebral_Exterior  40/*              205   62   78      0 */
#define   Right_Cerebral_White_Matter  41/*            0  225    0     0  */
#define   Right_Cerebral_Cortex  42/*                205   62   78      0  */
#define   Right_Lateral_Ventricle  43/*              120   18  134     0  */
#define   Right_Inf_Lat_Vent  44/*                   196   58  250     0   */
#define   Right_Cerebellum_Exterior  45/*              0  148    0     0 */
#define   Right_Cerebellum_White_Matter  46/*        220  248  164     0  */
#define   Right_Cerebellum_Cortex  47/*              230  148   34      0  */
#define   Right_Thalamus  48/*                         0  118   14      0  */
#define   Right_Thalamus_Proper  49/*                  0  118   14      0  */
#define   Right_Caudate  50/*                        122  186  220     0  */
#define   Right_Putamen  51/*                        236   13  176     0  */
#define   Right_Pallidum  52/*                       255   48  255     0   */
#define   Right_Hippocampus  53/*                    220  216   20      0  */
#define   Right_Amygdala  54/*                       103  255  255     0   */
#define   Right_Insula  55/*                          80  196   98      0  */
#define   Right_Operculum  56/*                       60   58  210     0   */
#define   Right_Lesion  57/*                         255  165    0     0   */
#define   Right_Accumbens_area  58/*                 255  165    0     0   */
#define   Right_Substancia_Nigra  59/*                 0  255  127     0   */
#define   Right_VentralDC  60/*                      165   42   42     0   */
#define   Right_undetermined  61/*                   135  206  235     0   */
#define   Right_vessel  62/*                         160   32  240     0   */
#define   Right_choroid_plexus  63/*                   0  255  255     0   */
#define   Right_F3orb  64/*                          100   50  100     0   */
#define   Right_lOg  65/*                            135   50   74     0   */
#define   Right_aOg  66/*                            122  135   50     0   */
#define   Right_mOg  67/*                             51   50  135     0   */
#define   Right_pOg  68/*                             74  155   60     0   */
#define   Right_Stellate  69/*                       120   62   43     0   */
#define   Right_Porg  70/*                            74  155   60     0   */
#define   Right_Aorg  71/*                          122  135   50     0*/
#define   Fifth_Ventricle 72/*                        156   25  250    0 */
#define   Left_Interior 73  /*                       122  135   50    0*/
#define   Right_Interior 74  /*                      122  135   50    0*/
#define   Left_Lateral_Ventricles 75  /*             120   18  134    0*/
#define   Right_Lateral_Ventricles 76  /*            120   18  134    0*/
#define   WM_hypointensities 77  /*                  124  140  178    0*/
#define   Left_WM_hypointensities 78  /*             124  140  178    0*/
#define   Right_WM_hypointensities 79  /*            124  140  178    0*/
#define   non_WM_hypointensities 80  /*              164  108  226    0*/
#define   Left_non_WM_hypointensities 81  /*         164  108  226    0*/
#define   Right_non_WM_hypointensities 82  /*        164  108  226    0*/
#define   Left_F1 83  /*                             255  218  185    0*/
#define   Right_F1 84  /*                            255  218  185    0*/
#define   Optic_Chiasm   85  /*                        234  169   30    0*/
#define   Left_Amygdala_Anterior  96 /*              205   10  125    0*/
#define   Right_Amygdala_Anterior 97 /*              205   10  125    0*/

/* no brain labels after this please unless you fix the IS_BRAIN macro */
#define   Dura           98
#define   Epidermis      118
#define   Conn_Tissue    119
#define   SC_FAT_MUSCLE  120
#define   Cranium        121
#define   CSF_SA         122
#define   Muscle         123
#define   Ear            124
#define   Fatty_Tissue   125
#define   Spinal_Cord    126
#define   Soft_Tissue    127
#define   Nerve          128
#define   Bone           129                            
#define   Air            130
#define   Orbit          131
#define   Tongue         132
#define   Nasal_Structures 133
#define   Globe          134
#define   Teeth          135

#define   Right_Temporal_Cerebral_White_Matter  186      /*  240  240  240     0 */
#define   Left_Temporal_Cerebral_White_Matter  187       /* 240  240  240     0 */
#define   Fat         189 /*  255 255 255 0 */
#define   Bright_Unknown  190   /* 100 240 240 0 */
#define   Dark_Unknown    191   /* 20  100  100     0  */
#define   Corpus_Callosum 192   /* 255 255 255 */

#define  alveus                   201 /*                    254	254	254	0*/
#define  perforant_pathway	      202 /*	 255	128	128	0*/
#define  parasubiculum	          203 /*	 255	255	0	0*/
#define  presubiculum		          204 /*	 64	0	64	0*/
#define  subiculum		            205 /*	 0	0	255	0*/
#define  CA1			                206 /*	 255	0	0	0*/
#define  CA2			                207 /*	 128	128	255	0*/
#define  CA3			                208 /*	 0	128	0	0*/
#define  CA4			                209 /*	 196	160	128	0*/
#define GC_DG		                  210 /*	 255	128	255	0*/
#define HATA			                211 /*	 128	255	128	0*/
#define fimbria                   212 /*	 0	0	0	0*/
#define lateral_ventricle	        213 /*	 255	0	255	0*/
#define molecular_layer_HP        214 /*  128    0       0       0*/
#define hippocampal_fissure       215 /*  64	255	196	0*/
#define entorhinal_cortex         216 /*  255	204	102	0*/
#define molecular_layer_subiculum 217 /*  0 0 0 0 */
#define Amygdala                  218 /*  103  255  255    0    */
#define Cerebral_White_Matter     219 /*  0  225    0     0   */
#define Cerebral_Cortex           220 /*	 205   62   78     0 */
#define Inf_Lat_Vent              221 /*  196   58  250     0 */

#define Left_hippocampal_fissure 193 //		0	196	255	0
#define Left_CADG_head 194 //			255	164	164	0
#define Left_subiculum 195 //			196	196	0	0
#define Left_fimbria 196 //			0	100	255	0
#define Right_hippocampal_fissure 197 //		128	196	164	0
#define Right_CADG_head 198 //			0	196	0	0
#define Right_subiculum 199 //			128	96	64	0
#define Right_fimbria 200 //			0	50	128	0

// vascular and lymph labels (from Alex G)
#define Aorta 331 // 255 0 0 0
#define Left_Common_IliacA 332 // 255 80 0 0
#define Right_Common_IliacA 333 // 255 160 0 0
#define Left_External_IliacA 334 // 255 255 0 0
#define Right_External_IliacA 335 // 0 255 0 0
#define Left_Internal_IliacA 336 // 255 0 160 0
#define Right_Internal_IliacA 337 // 255 0 255 0
#define Left_Lateral_SacralA 338 // 255 50 80 0
#define Right_Lateral_SacralA 339 // 80 255 50 0
#define Left_ObturatorA 340 // 160 255 50 0
#define Right_ObturatorA 341 // 160 200 255 0
#define Left_Internal_PudendalA 342 // 0 255 160 0
#define Right_Internal_PudendalA 343 // 0 0 255 0
#define Left_UmbilicalA 344 // 80 50 255 0
#define Right_UmbilicalA 345 // 160 0 255 0
#define Left_Inf_RectalA 346 // 255 210 0 0
#define Right_Inf_RectalA 347 // 0 160 255 0
#define Left_Common_IliacV 348 // 255 200 80 0
#define Right_Common_IliacV 349 // 255 200 160 0
#define Left_External_IliacV 350 // 255 80 200 0
#define Right_External_IliacV 351 // 255 160 200 0
#define Left_Internal_IliacV 352 // 30 255 80 0
#define Right_Internal_IliacV 353 // 80 200 255 0
#define Left_ObturatorV 354 // 80 255 200 0
#define Right_ObturatorV 355 // 195 255 200 0
#define Left_Internal_PudendalV 356 // 120 200 20 0
#define Right_Internal_PudendalV 357 // 170 10 200 0
#define Pos_Lymph 358 // 20 130 180 0
#define Neg_Lymph 359 // 20 180 130 0

#define MAX_LABEL Neg_Lymph
#define MAX_CMA_LABEL (MAX_LABEL)
#define MAX_CMA_LABELS (MAX_CMA_LABEL+1)

#define   Not_Set         255


#define IS_UNKNOWN(label)  (((label) == Unknown) || (label == 255) || (label == Bright_Unknown) || (label == Dark_Unknown))

#define IS_BRAIN(label)  (!IS_UNKNOWN(label) && label < Dura)

#define IS_WM(label) (((label) == Left_Cerebral_White_Matter) || ((label) == Right_Cerebral_White_Matter) || ((label) == Left_Temporal_Cerebral_White_Matter) || ((label) == Right_Temporal_Cerebral_White_Matter))
#define IS_HYPO(label) (((label) == WM_hypointensities)  || ((label) == Left_WM_hypointensities)  || ((label) == Left_WM_hypointensities))
#define IS_WMH(label) (IS_WM(label) || IS_HYPO(label))
#define IS_THALAMUS(label)  (((label) == Left_Thalamus) || ((label) == Left_Thalamus_Proper) || ((label) == Right_Thalamus) || ((label) == Right_Thalamus_Proper))
#define IS_GM(label) (((label) == Left_Cerebral_Cortex) || ((label) == Right_Cerebral_Cortex))


#define IS_CEREBELLAR_WM(label) (((label) == Left_Cerebellum_White_Matter) || ((label) == Right_Cerebellum_White_Matter))
#define IS_CEREBELLAR_GM(label) (((label) == Left_Cerebellum_Cortex) || ((label) == Right_Cerebellum_Cortex))

#define IS_HIPPO(l) (((l) == Left_Hippocampus) || ((l) == Right_Hippocampus))
#define IS_AMYGDALA(l) (((l) == Left_Amygdala) || ((l) == Right_Amygdala))
#define IS_CORTEX(l) (((l) == Left_Cerebral_Cortex) || \
                      ((l) == Right_Cerebral_Cortex))
#define IS_LAT_VENT(l) (((l) == Left_Lateral_Ventricle) || \
                        ((l) == Right_Lateral_Ventricle) || \
                        ((l) == Right_Inf_Lat_Vent) || \
                        ((l) == Left_Inf_Lat_Vent))
#define IS_CSF(l) (IS_LAT_VENT(l) || ((l) == CSF) || ((l) == CSF_SA) || ((l) == Third_Ventricle) || ((l) == Fourth_Ventricle))

#define IS_INF_LAT_VENT(l)  (((l) == Left_Inf_Lat_Vent) || ((l) == Right_Inf_Lat_Vent))
#define IS_CAUDATE(l) (((l) == Left_Caudate) || ((l) == Right_Caudate))
#define IS_PUTAMEN(l) (((l) == Left_Putamen) || ((l) == Right_Putamen))
#define IS_PALLIDUM(l) (((l) == Left_Pallidum) || ((l) == Right_Pallidum))

#define LABEL_WITH_NO_TOPOLOGY_CONSTRAINT(l) (\
   ((l) == Right_non_WM_hypointensities) || \
   ((l) == Left_non_WM_hypointensities) || \
   ((l) == Left_WM_hypointensities) || \
   ((l) == Right_WM_hypointensities) || \
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

typedef struct {
  int n_claims;
  int interior_claim_flag;
  short claim_labels[MAX_OUTLINE_CLAIMS];
  float claim_values[MAX_OUTLINE_CLAIMS];
  float no_label_claim;
} CMAoutlineClaim;

typedef struct {
  int width, height;
  CMAoutlineClaim **claim_field;
  unsigned char **fill_field;
  unsigned char **outline_points_field;
} CMAoutlineField;

CMAoutlineField *CMAoutlineFieldAlloc(int width, int height);
int CMAfreeOutlineField(CMAoutlineField **of);

int CMAclearFillField(CMAoutlineField *field);
int CMAfill(CMAoutlineField *field, short seed_x, short seed_y);

int CMAclaimPoints(CMAoutlineField *field, short label, short *points, int n_points, short seed_x, short seed_y);
int CMAassignLabels(CMAoutlineField *field);

int CMAvalueClaims(CMAoutlineClaim *claim);
int CMAvalueAllClaims(CMAoutlineField *field);

/* the weights to give claims of nearby points -- see cma.c */
#define OTL_CLAIM_WEIGHT_CENTER   1.000
#define OTL_CLAIM_WEIGHT_SIDE     0.500
#define OTL_CLAIM_WEIGHT_DIAGONAL 0.707

short CMAtotalClaims(CMAoutlineField *field, int x, int y);
int CMAaddWeightedTotals(CMAoutlineClaim *claim, float weight, float *claim_totals);

int CMAzeroOutlines(CMAoutlineField *field);
char *cma_label_to_name(int label) ;

#define IS_FIMBRIA(l) ((l) == Left_fimbria || (l) == Right_fimbria || (l) == fimbria)
#define CSF_CLASS        0
#define GM_CLASS         1
#define WM_CLASS         2
#define NTISSUE_CLASSES  3


// only structures that are 'pure' of one tissue type
#define IS_GRAY_CLASS(l) (IS_GM(l) || ((l) == Left_Caudate) || ((l) == Right_Caudate))

#define IS_WHITE_CLASS(l) (((l) == Left_Cerebral_White_Matter) || ((l) == Right_Cerebral_White_Matter))

#define IS_CSF_CLASS(l) (((l) == Left_Lateral_Ventricle) || ((l) == Right_Lateral_Ventricle) || ((l) == CSF) || ((l) == CSF_SA) || ((l) == Third_Ventricle) || ((l) == Fourth_Ventricle) || ((l) == Fifth_Ventricle) || ((l) == Left_hippocampal_fissure) || ((l) == Right_hippocampal_fissure) || ((l) == hippocampal_fissure))

#define IS_CLASS(l,c) (c == CSF_CLASS ? IS_CSF_CLASS(l) : c == GM_CLASS ? IS_GRAY_CLASS(l) : IS_WHITE_CLASS(l))
#endif
