//
// ctest.c
//
// purpose: usability of cpp library exported C function
//
#include <stdio.h>
#include <stdlib.h>
#include "mri.h"
#include "fastmarching.h"

char *Progname = "ctest";

int main(int argc, char *argv[])
{
  MRI *src;
  MRI *dst;
  int label=2;
  float max_distance=10.;
  int mode=1;
  
  if (argc < 1)
  {
    fprintf(stderr, "Usage: ctest <volume>\n");
    return -1;
  }
  src = MRIread(argv[1]);
  dst = MRIextractDistanceMap(src, NULL, label, max_distance, mode);
  return 0;
}

