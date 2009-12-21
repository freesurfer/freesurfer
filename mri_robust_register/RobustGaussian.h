/**
 * @file RobustGaussian.h
 * @brief A class to esimate a robust Gaussian (using median and mad)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/12/21 19:46:30 $
 *    $Revision: 1.6 $
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


class RobustGaussian
{
public:
  static double median(double a[],int n);
  static double kth_smallest(double a[], int n, int k);
  static double quick_select(double a[], int n, int k);
  static std::pair <double, double> medianI(double a[],int n);
  static std::pair <double, int> kth_smallestI(double a[], int n, int k);
  static std::pair <double, int> quick_selectI(double a[], int n, int k);
  static double mad(double a[],int n, double d =  1.4826);
};


#endif
