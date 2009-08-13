
#include "Quaternion.h"


// ---------------------------------------------------------- global functions


// Multiplication of scalar and Quaternion from left.
Quaternion operator*(const double& scalar, const Quaternion& vect)
{
  return vect*scalar;
}

std::ostream& operator<<(std::ostream& os, const Quaternion& q)
{
  q.write(os);
  return os;
}
