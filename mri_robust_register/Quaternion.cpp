/**
 * @file Quaternion.cpp
 * @brief A class representing Quaterions
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:24 $
 *    $Revision: 1.3 $
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
