/**
 * @brief linear registration to a gca atlas
 *
 * Header file for findtranslation.cpp
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

#ifndef EM_REGISTER_FIND_TRANSLATION_H
#define EM_REGISTER_FIND_TRANSLATION_H

#include "mri.h"
#include "gca.h"
#include "matrix.h"

double find_optimal_translation( GCA *gca,
                              GCA_SAMPLE *gcas,
                              MRI *mri,
                              int nsamples,
                              MATRIX *m_L,
                              float min_trans,
                              float max_trans,
                              float trans_steps,
                              int nreductions,
                              double clamp);

#endif
