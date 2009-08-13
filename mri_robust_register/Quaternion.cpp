/**
 * @file Quaternion.cpp
 * @brief A class representing Quaterions
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2009/08/13 02:51:19 $
 *    $Revision: 1.2 $
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
