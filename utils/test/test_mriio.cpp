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


#include <stdexcept>
#include <sstream>
#include <iostream>

extern "C"
{
#include "mri.h"
#include "NrrdIO.h"
}

using namespace std;

const char *Progname = NULL;

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  throw runtime_error( ss.str() ); \
  }

class MriioTester
{
public:
  void TestNrrdIO ();
};


void
MriioTester::TestNrrdIO ()
{
  cerr << "Check that mriNrrdRead reads foolc.nrrd...";
  Assert(MRIreadType("test_mriio_data/foolc.nrrd", NRRD_FILE) != NULL,
         "failed.");
  cerr << "passed." << endl;
}


int main ( int argc, char** argv )
{
  //Progname = argv[0];

  cerr << "Beginning tests..." << endl;

  try
  {

    MriioTester tester;
    tester.TestNrrdIO();

  }
  catch ( runtime_error& e )
  {
    cerr << "failed " << endl << "exception: " << e.what() << endl;
    exit( 1 );
  }
  catch (...)
  {
    cerr << "failed" << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;

  exit( 0 );
}

