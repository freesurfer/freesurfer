/**
 * @file  cvector.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:08:59 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#ifndef CVECTOR_H
#define CVECTOR_H

#include "mrisurf.h"

#define STAT_T           0
#define STAT_F           1
#define STAT_MEAN        2

int   cvector_scalar_mul(float *v, float m, int num) ;
double   cvector_compute_pvalues(float *c1_mean, float *c1_var,
                                 float *c2_mean, float *c2_var,
                                 int num_class1, int num_class2,
                                 float *pvals, int num, int stat_type,
                                 int *pvno) ;
double cvector_average_in_label(float *v, LABEL *area, int num) ;
int   cvector_extract_best_avg(float *vbest_avgs, float *vsrc,
                               float *vdst, int navgs, int num) ;
int   cvector_normalize(float *v, float norm, int num) ;
int   cvector_accumulate(float *v, float *vtotal, int num) ;
int   cvector_accumulate_square(float *v, float *vtotal, int num) ;
int   cvector_compute_variance(float *var,float *mean,int norm,int num);
double   cvector_compute_t_test(float *c1_mean, float *c1_var,
                                float *c2_mean, float *c2_var,
                                int num_class1, int num_class2,
                                float *pvals, int num, int *pvno) ;
double   cvector_compute_mean_diff(float *c1_mean, float *c2_mean,
                                   float *pvals, int num, int *pvno) ;
int    cvector_copy(float *v1, float *v2, int num) ;
double cvector_compute_dist_free_snr(float **c1_thickness,
                                     int num_class1,
                                     float **c2_thickness,
                                     int num_class2,
                                     float *c1_mean, float *c2_mean,
                                     float *vsnr, int num, int *pi);
double cvector_compute_snr(float *c1_mean, float *c2_mean,
                           float *verror, float *snr, int num, int *pi,
                           float bonferroni, int stat_type);
double cvector_compute_snr_F(float *c1_mean, float *c2_mean,
                             float *verror, float *snr, int num, int *pi,
                             float bonferroni);


float *cvector_alloc(int num) ;
int   cvector_clear(float *v, int num) ;
int   cvector_set(float *v, float val, int num) ;
int   cvector_add_variances(float *c1_var, float *c2_var,
                            int num_class1, int num_class2,
                            float *total_var, int nvertices) ;
int   cvector_multiply_variances(float *c1_var, float *c2_var,
                                 int num_class1, int num_class2,
                                 float *total_var, int nvertices) ;
int   cvector_track_best_snr(float *vsnr, float *vbest_snr,
                             float *vbest_avgs,
                             float *c1_mean, float *c2_mean,
                             float *c1_best_mean, float *c2_best_mean,
                             float *c1_var, float *c2_var,
                             float *c1_best_var, float *c2_best_var,
                             float **c1_avg_thickness,
                             float **c2_avg_thickness,
                             float **c1_best_thicknesses, int nc1,
                             float **c2_best_thicknesses, int nc2,
                             int avgs, int num,
                             float fthresh, int *pnum_found) ;

int   cvector_track_best_stats(float *vsnr, float *vbest_snr,
                               float *vbest_avgs,
                               float *c1_mean, float *c2_mean,
                               float *c1_best_mean, float *c2_best_mean,
                               float *c1_var, float *c2_var,
                               float *c1_best_var, float *c2_best_var,
                               float **c1_avg_thickness,
                               float **c2_avg_thickness,
                               float **c1_best_thicknesses, int nc1,
                               float **c2_best_thicknesses, int nc2,
                               int avgs, int num,
                               float fthresh, int *pnum_found) ;
#endif
