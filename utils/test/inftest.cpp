/**
 * @file  inftest.cpp
 * @brief infinity tests
 *
 */
/*
 * Original Author: Y. Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/10 19:47:42 $
 *    $Revision: 1.5.2.1 $
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

#include <iostream>
#include <math.h>

extern "C"
{
#include "mghendian.h"
#include "utils.h"
const char *Progname = "inftest";
}

using namespace std;

int main(int argc, char *argv[])
{
  int fails=0;
  float v = 0.f;
  int res = devIsnan(0.f/v);
  if (res != 1)
  {
    cerr << res << " should be nan when 0.f/0.f" << endl;
    fails++;
  }
  res= devIsnan(1.f);
  if (res != 0)
  {
    cerr << res << " should not be nan when 0.f/0.f" << endl;
    fails++;
  }
  res= devIsinf(1.f);
  if (res != 0)
  {
    cerr << res << " should not be inf when 1.f" << endl;
    fails++;
  }
#ifndef SunOS
  res = devIsinf(HUGE_VALF);
  if (res != 1)
  {
    cerr << res << " should be +inf for HUGE_VALF" << endl;
    fails++;
  }
  res = devIsinf(-HUGE_VALF);
  if (res != -1)
  {
    cerr << res << " should be -inf for HUGE_VALF" << endl;
    fails++;
  }
#endif
  if (fails) return 77; // don't indicate out-right failure
}
