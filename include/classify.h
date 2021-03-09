/*
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


#ifndef CLASSIFY_H
#define CLASSIFY_H


#define CLASSIFIER_GAUSSIAN        0
#define CLASSIFIER_BACKPROP        1
#define CLASSIFIER_ARTMAP          2
#define CLASSIFIER_GAUSSIAN_ARRAY  3
#define CLASSIFIER_RBF             4   /* radial basis function */
#define CLASSIFIER_SVM             5
#define CLASSIFIER_RFOREST         6

// different types of deep learning networks
#define NET_DBN  0
#define NET_RBM  1
#define NET_CDBN 2
#define NET_SAE  3

#endif
