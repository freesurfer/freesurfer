/**
 * @file  gcaread.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:44 $
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
