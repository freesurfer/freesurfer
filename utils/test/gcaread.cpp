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

#include "gca.h"
#include "mri.h"

const char *Progname = "gcaread";

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Usage: gcaread <gcafile>" << std::endl;
    return -1;
  }
  GCA *gca = GCAread(argv[1]);
  if (gca == 0) {
    std::cout << "could not open file " << argv[1] << std::endl;
  } else
    GCAfree(&gca);
}
