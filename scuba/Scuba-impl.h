/**
 * @file  Scuba-impl.h
 * @brief Implementations for templated classes
 *
 * This file is necessary for creating the instantiations for template
 * classes. Note that we actually include the .cpp files here.  Then
 * one line for each template class or function using a specific
 * type. See
 * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:36 $
 *    $Revision: 1.15 $
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

#include "Array2.cpp"
#include "Path.cpp"
#include "Point2.h"
#include "Point3.h"
#include "ShortestPathFinder.h"
#include "Volume3.cpp"
#include "VolumeCollection.h"

using namespace std;

template class Volume3<bool>;
template class Array2<float>;
template class Array2<int>;
template class Array2<bool>;
template class Path<float>;
template class Volume3<Point3<int> >;
template class Array2<Point3<float> >;
template class Array2<VolumeLocation*>;
template class Volume3<VolumeLocation*>;
