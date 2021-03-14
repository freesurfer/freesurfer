/**
 * @brief test routines
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

#include "mri.h"
#include "mri_identify.h"

const char *Progname = NULL;

using namespace std;

#define Assert(x, s)                                                           \
  if (!(x)) {                                                                  \
    stringstream ss;                                                           \
    ss << "Line " << __LINE__ << ": " << s;                                    \
    throw std::runtime_error(ss.str());                                        \
  }

class mri_identifyTester {
public:
  void Test();
};

void mri_identifyTester::Test() {
  std::cerr << "Check that mri_identify() returns NRRD_FILE for a .nrrd...";
  Assert(mri_identify((char *)"test_mri_identify_data/nrrd_basic.nrrd") ==
             NRRD_FILE,
         "Failed to recognize NRRD file w/ correct extension.");
  std::cerr << "passed." << std::endl;

  std::cerr << "Check that a .nrrd is recognized by is_nrrd()...";
  Assert(is_nrrd((char *)"test_mri_identify_data/nrrd_basic.nrrd"),
         "Failed to recognize NRRD file w/ correct extension.");
  std::cerr << "passed." << std::endl;

  std::cerr
      << "Check recognition of NRRD by magic number when extension is bad...";
  Assert(is_nrrd((char *)"test_mri_identify_data/nrrd_bad_extension.something"),
         "Failed to identify NRRD w/ bad extension but correct magic.");
  std::cerr << "passed." << std::endl;

  std::cerr
      << "Check recognition of not NRRD by magic number when extension is "
         "bad...";
  Assert(!is_nrrd((char *)"test_mri_identify_data/nrrd_bad_magic.something"),
         "Failed to identify bad magic NRRD as not a NRRD.");
  std::cerr << "passed." << std::endl;
}

int main(int argc, char **argv) {
  Progname = argv[0];

  std::cerr << "Beginning tests..." << std::endl;

  try {

    mri_identifyTester tester;
    tester.Test();

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
