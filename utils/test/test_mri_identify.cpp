/**
 * @file  test_mri_identify.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:46 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdexcept>
#include <sstream>
#include <iostream>

extern "C"
{
#include "mri_identify.h"
#include "mri.h"
}

using namespace std;

char *Progname = NULL; // to satisfy mriio.c's extern Progname declaration

#define Assert(x,s)   \
  if(!(x)) { \
  stringstream ss; \
  ss << "Line " << __LINE__ << ": " << s; \
  throw runtime_error( ss.str() ); \
  }

class mri_identifyTester
{
public:
  void Test ();
};


void
mri_identifyTester::Test ()
{
  cerr << "Check that mri_identify() returns NRRD_FILE for a .nrrd...";
  Assert(mri_identify("test_mri_identify_data/nrrd_basic.nrrd") == NRRD_FILE,
         "Failed to recognize NRRD file w/ correct extension.");
  cerr << "passed." << endl;

  cerr << "Check that a .nrrd is recognized by is_nrrd()...";
  Assert(is_nrrd("test_mri_identify_data/nrrd_basic.nrrd"),
         "Failed to recognize NRRD file w/ correct extension.");
  cerr << "passed." << endl;

  cerr << "Check recognition of NRRD by magic number when extension is bad...";
  Assert(is_nrrd("test_mri_identify_data/nrrd_bad_extension.something"),
         "Failed to identify NRRD w/ bad extension but correct magic.");
  cerr << "passed." << endl;

  cerr << "Check recognition of not NRRD by magic number when extension is bad...";
  Assert(!is_nrrd("test_mri_identify_data/nrrd_bad_magic.something"),
         "Failed to identify bad magic NRRD as not a NRRD.");
  cerr << "passed." << endl;
}


int main ( int argc, char** argv )
{
  Progname = argv[0];

  cerr << "Beginning tests..." << endl;

  try
  {

    mri_identifyTester tester;
    tester.Test();

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

