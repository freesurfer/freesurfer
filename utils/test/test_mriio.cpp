/**
 * @brief mriio test routines
 *
 */
/*
 * Original Author: Y. Tosa
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

#include <iostream>
#include <sstream>
#include <stdexcept>

#include "NrrdIO.h"
#include "mri.h"

const char *Progname = NULL;

using namespace std;

#define Assert(x, s)                                                           \
  if (!(x)) {                                                                  \
    stringstream ss;                                                           \
    ss << "Line " << __LINE__ << ": " << s;                                    \
    throw std::runtime_error(ss.str());                                        \
  }

class MriioTester {
public:
  void TestNrrdIO();
};

void MriioTester::TestNrrdIO() {
  std::cerr << "Check that mriNrrdRead reads foolc.nrrd...";
  Assert(MRIreadType("test_mriio_data/foolc.nrrd", NRRD_FILE) != NULL,
         "failed.");
  std::cerr << "passed." << std::endl;
}

int main(int argc, char **argv) {
  // Progname = argv[0];

  std::cerr << "Beginning tests..." << std::endl;

  try {

    MriioTester tester;
    tester.TestNrrdIO();

  } catch (runtime_error &e) {
    std::cerr << "failed " << std::endl
              << "exception: " << e.what() << std::endl;
    exit(1);
  } catch (...) {
    std::cerr << "failed" << std::endl;
    exit(1);
  }

  std::cerr << "Success" << std::endl;

  exit(0);
}
