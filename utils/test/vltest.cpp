//
// vltest.cpp
//

#include <iostream>
#include <fstream>
#include <mcheck.h>

extern "C" {
#include "label.h"
#include "vlabels.h"
#include "transform.h"
  char *Progname = "vltest";
}

using namespace std;

int main(int, char **)
{
  VLI *pVli = 0;
  VL ***voxel_labels;
  VOXEL_LABELS *vl;

  mtrace();
  for (size_t i=0; i < 10; ++i)
  {
    pVli = VLread(const_cast<char *>("vltest.dat"));
    voxel_labels = pVli->vl;
    vl = &voxel_labels[0][0][0];
    VLfree(&pVli);
  }
}
