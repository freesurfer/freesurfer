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
#define   Bone        80 /*  100 100 100 0 */                         
#define   Fat         81 /*  255 255 255 0 */
#define   Bright_Unknown  82   /* 100 240 240 0 */
#define   Dark_Unknown    83   /* 20  100  100     0  */



#define IS_UNKNOWN(label)  (((label) == Unknown) || (label < 0) || (label == Bright_Unknown) || (label == Dark_Unknown))

#define IS_BRAIN(label)  (!IS_UNKNOWN(label) && (label != Fat) && (label != Bone))

#define MAX_CMA_LABEL 83

#endif
