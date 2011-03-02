/**
 * @file  mri_remove.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
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


//
// mri_remove
//

#include <iostream>

extern "C"
{
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
  nargs = handle_version_option (argc, argv, "$Id: mri_remove.cpp,v 1.3 2011/03/02 00:04:55 nicks Exp $", "$Name:  $");
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
