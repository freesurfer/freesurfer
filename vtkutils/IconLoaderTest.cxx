/**
 * @file  IconLoaderTest.cxx
 * @brief Tests IconLoader class
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/06/17 20:43:19 $
 *    $Revision: 1.1.2.1 $
 *
 * Copyright (C) 2008,
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

#include <iostream>
#include <stdexcept>

#include "IconLoader.h"

using namespace std;

class IconLoaderTest
{
public:

  static int Test () {
    // Init the icon loader with the app and load our icons.
    try {
      IconLoader::Initialize( NULL );
      
      try {
        IconLoader::LoadIconsFromFile( "./IconLoaderTestIcons.txt" );
      }
      catch(...) {
        char* pfnFreesurferDir = getenv( "FREESURFER_HOME" );
        if( NULL != pfnFreesurferDir ) {
          string fnIcons =
            string(pfnFreesurferDir) + "/lib/resource/QdecIcons.txt";
          IconLoader::LoadIconsFromFile( fnIcons.c_str() );
        }
      }
    }
    catch( exception& e ) {
      cerr << "Error loading icons: " << e.what() << endl;
      return 1;
    }
    return 0;
  };
};

int main ( int argc, char** argv ) {

  return IconLoaderTest::Test();
}
