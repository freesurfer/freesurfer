/**
 * @file RobustGaussian.h
 * @brief A class to esimate a robust Gaussian (using median and mad)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/07/17 02:35:08 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2008-2009
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
