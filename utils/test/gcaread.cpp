/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
// gcaread
//

#include <iostream>

extern "C"
{
#include "mri.h"
#include "gca.h"

  char *Progname = "gcaread";
}

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: gcaread <gcafile>" << endl;
    return -1;
  }
  GCA *gca = GCAread(argv[1]);
  if (gca == 0)
  {
    cout << "could not open file " << argv[1] << endl;
  }
  else
    GCAfree(&gca);
}
