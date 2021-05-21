/**
 * @brief linear registration to a gca atlas
 *
 * Various utilities
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

#ifndef EM_REGISTER_UTILS_H
#define EM_REGISTER_UTILS_H

#include "mri.h"
#include "gca.h"
#include "matrix.h"

extern int exvivo;
extern int robust;
extern float G_wm_mean, G_gm_mean, G_fluid_mean;

double local_GCAcomputeLogSampleProbability( GCA *gca,
                                             GCA_SAMPLE *gcas,
                                             MRI *mri,
                                             MATRIX *m_L,
                                             int nsamples,
                                             int exvivo, double clamp );

int compute_tissue_modes( MRI *mri_inputs,
                          GCA *gca,
                          GCA_SAMPLE *gcas,
                          TRANSFORM *transform,
                          int nsamples,
                          double *pwm, double *pgm, double *pfluid );

#endif
