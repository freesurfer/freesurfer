#include "Volume3.cpp"
#include "Point2.cpp"
#include "Point3.cpp"
#include "Utilities.h"

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
