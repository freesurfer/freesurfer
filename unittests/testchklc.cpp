#include <stdio.h>
#include <stdlib.h>

#include "diag.h"
#include "chklc.h"

// usage: testchklc [license-file] [0|1]
int main(int argc, char *argv[])
{
  Gdiag_no = 1;

  if (argc > 1)
  {
    // testchklc <license-file>
    setenv("FS_LICENSE", argv[1], 1);

    const char *fslic = getenv("FS_LICENSE");
    if (fslic != NULL)
      printf("FS_LICENSE = %s\n", fslic);
  }

  int mode = 0;   // 1=freeview, 0=other freesurfer tools
  if (argc > 2)   // testchklc <license-file> [0|1]
    mode = atoi(argv[2]);

  if (mode)
  {
    printf("!!!CHECK LICENSE IN FREEVIEW MODE!!!\n");
    char license_msg[2000];
    if (!chklc(license_msg))
      printf("License Error: %s", license_msg);
  }
  else
    chklc();

  exit(0);
}
