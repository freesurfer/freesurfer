/**
 * @file RobustGaussian.h
 * @brief A class to esimate a robust Gaussian (using median and mad)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.8 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#ifndef RobustGaussian_H
#define RobustGaussian_H

#include <utility>


template <class T>
class RobustGaussian
{
public:
  static T median(T a[],int n);
  static T kth_smallest(T a[], int n, int k);
  static T quick_select(T a[], int n, int k);
  static std::pair <T, T> medianI(T a[],int n);
  static std::pair <T, int> kth_smallestI(T a[], int n, int k);
  static std::pair <T, int> quick_selectI(T a[], int n, int k);
  static T mad(T a[],int n, T d =  1.4826);
};

#include "RobustGaussian.cpp"

#endif
