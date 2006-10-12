// 
// inftest.cpp
//

#include <iostream>
#include <math.h>

extern "C" {
#include "mghendian.h"
#include "utils.h"
  char *Progname = "inftest";
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
  if (fails) return 77; // don't indicate out-right failure
}
