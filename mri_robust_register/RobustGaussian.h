/**
 * @file RobustGaussian.h
 * @brief A class to esimate a robust Gaussian (using median and mad)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/09/21 23:05:16 $
 *    $Revision: 1.10 $
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

/** \class RobustGaussian
 * \brief A static class with routines for robust compuations (median, MAD, ...)
 */
template<class T>
class RobustGaussian
{
public:
  //! Find median of array a with length n
  static T median(T a[], int n);
  //! Find k-th smallest of array a with length n
  static T kth_smallest(T a[], int n, int k);
  //! Find k-th smallest of array a with length n (using quick select)
  static T quick_select(T a[], int n, int k);
  //! Find median of array a with length n and the index
  static std::pair<T, T> medianI(T a[], int n);
  //! Find k-th smallest of array a with length n and the index
  static std::pair<T, int> kth_smallestI(T a[], int n, int k);
  //! Find k-th smallest of array a with length n and the index (using quick select)
  static std::pair<T, int> quick_selectI(T a[], int n, int k);
  //! Find median absolute deviation
  static T mad(T a[], int n, T d = 1.4826);
};

#include "RobustGaussian.cpp"

#endif
