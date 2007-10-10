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
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/10 18:43:47 $
 *    $Revision: 1.10 $
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

#include "Volume3.cpp"
#include "Point2.cpp"
#include "Point3.cpp"
#include "Array2.cpp"
#include "Path.cpp"
#include "ShortestPathFinder.h"
#include "VolumeCollection.h"

using namespace std;

template class Volume3<bool>;
template class Point3<int>;
template ostream& operator << ( ostream&, Point3<int> const& );
template class Point3<float>;
template ostream& operator << ( ostream&, Point3<float> const& );
template class Point2<int>;
template ostream& operator << ( ostream&, Point2<int> const& );
template class Point2<float>;
template ostream& operator << ( ostream&, Point2<float> const& );
template class Point2<Point3<float> >;
template class Array2<float>;
template class Array2<int>;
template class Array2<bool>;
template class Array2<ShortestPathFinder::listElement>;
template class Path<float>;
template class Volume3<Point3<int> >;
template class Array2<Point3<float> >;
template class Array2<VolumeLocation*>;
template class Volume3<VolumeLocation*>;
