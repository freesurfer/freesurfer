/**
 * @file  test_mriio.cpp
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
#include "mri.h"
#include "NrrdIO.h"
}

using namespace std;

char *Progname = NULL; // to satisfy mriio.c's extern Progname declaration

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

