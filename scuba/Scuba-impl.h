/**
 * @file  Scuba-impl.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:14 $
 *    $Revision: 1.9 $
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


//
// Scuba-impl.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2006/12/29 02:09:14 $
// Revision       : $Revision: 1.9 $

// This file is necessary for creating the instantiations for template
// classes. Note that we actually include the .cpp files here.  Then
// one line for each template class or function using a specific
// type. See
// http://www.parashift.com/c++-faq-lite/containers-and-templates.html#faq-34.14



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
template ostream& operator << ( ostream&, Point3<int> );
template class Point3<float>;
template ostream& operator << ( ostream&, Point3<float> );
template class Point2<int>;
template ostream& operator << ( ostream&, Point2<int> );
template class Point2<float>;
template ostream& operator << ( ostream&, Point2<float> );
template class Array2<float>;
template class Array2<int>;
template class Array2<bool>;
template class Array2<listElement>;
template class Path<float>;
DeclareIDTracker(Path<float>); //ugh
template class Volume3<Point3<int> >;
template class Array2<Point3<float> >;
template class Array2<VolumeLocation*>;
template class Volume3<VolumeLocation*>;
