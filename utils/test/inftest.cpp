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
  float v = 0.f;
  int res = devIsnan(0.f/v);
  if (res != 1)
  {
    cerr << "should be nan when 0.f/0.f" << endl;
    return -1;
  }
  res= devIsnan(1.f);
  if (res != 0)
  {
    cerr << "should not be nan when 0.f/0.f" << endl;
    return -1;
  }
  res= devIsinf(1.f);
  if (res != 0)
  {
    cerr << "should not be inf when 1.f" << endl;
    return -1;
  }
  res = devIsinf(HUGE_VALF);
  if (res != 1)
  {
    cerr << "should be +inf for HUGE_VALF" << endl;
    return -1;
  }
  res = devIsinf(-HUGE_VALF);
  if (res != -1)
  {
    cerr << "should be -inf for HUGE_VALF" << endl;
    return -1;
  }
}
