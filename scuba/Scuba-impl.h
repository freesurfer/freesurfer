//
// Scuba-impl.h
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: kteich $
// Revision Date  : $Date: 2004/04/12 19:21:21 $
// Revision       : $Revision: 1.2 $

// This file is necessary for creating the instantiations for template
// classes. Note that we actually include the .cpp files here.  Then
// one line for each template class or function using a specific
// type. See
// http://www.parashift.com/c++-faq-lite/containers-and-templates.html#faq-34.14



#include "Volume3.cpp"
#include "Point2.cpp"
#include "Point3.cpp"
#include "Array2.cpp"
#include "ShortestPathFinder.h"

using namespace std;

template Volume3<bool>;
template Point3<int>;
template ostream& operator << ( ostream&, Point3<int> );
template Point3<float>;
template ostream& operator << ( ostream&, Point3<float> );
template Point2<int>;
template ostream& operator << ( ostream&, Point2<int> );
template Point2<float>;
template ostream& operator << ( ostream&, Point2<float> );
template Array2<float>;
template Array2<int>;
template Array2<bool>;
template Array2<listElement>;
