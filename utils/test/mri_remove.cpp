//
// mri_remove
//

#include <iostream>

extern "C" {
#include "mri.h"
#include "version.h"
  char *Progname = "mri_remove";
}

using namespace std;

int get_option(int argc, char *argv[], int *yval)
{
  int nargs = 0;
  char *option;
  option=argv[1] + 1;
  if (!strcmp(option, "-help"))
  {
    cout << "Usage: mri_remove [--ypos <pos>] <vol> <volchopped>" << endl;
    cout << "where <volchopped> has voxel val set to zero higher that <pos>" << endl;
    cout << "      if --ypos is not given, use pos = 170" << endl;
    exit (-1);
  }
  else if (!strcmp(option, "-ypos"))
  {
    *yval = atoi(argv[2]);
    nargs = 1;
  }
  return nargs;
}

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: mri_remove [--ypos <pos>] <vol> <volchopped>" << endl;
    return -1;
  }
  int nargs;
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_remove.cpp,v 1.1 2004/10/12 20:37:48 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  int yval = 170;
  // argument handling
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv, &yval) ;
    argc -= nargs ;
    argv += nargs ;
  }

  MRI *mri=MRIread(argv[1]);
  for (int z =0; z < mri->depth; z++)
    for (int y=0; y < mri->height; y++)
      for (int x=0; x < mri->width; x++)
      {
	if (y > yval)
	  MRIvox(mri, x, y, z) = 0;
      }

  MRIwrite(mri, argv[2]);
}
