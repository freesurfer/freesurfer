/**
 * @file  mris_label_calc.cpp
 * @brief A programm to calculate stuff on surface labels
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/02/12 01:48:12 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2008-2009
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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

// all other software are all in "C"
#ifdef __cplusplus
extern "C"
{
#endif
#include "error.h"
#include "macros.h"
#include "mri.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"
#include "label.h"

#ifdef __cplusplus
}
#endif

using namespace std;

//static char vcid[] = "$Id: mris_label_calc.cpp,v 1.1 2010/02/12 01:48:12 mreuter Exp $";
char *Progname = NULL;


int main(int argc, char *argv[])
{

//   // Default initialization
//   int nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
//   if (nargs && argc - nargs == 1) exit (0);
//   argc -= nargs;
  Progname = argv[0] ;
//   argc --;
//   argv++;
//  if (vcid)
//  {};

  if (argc < 3)
  {
    cout << endl;
    cout << argv[0] << " command input1 input2 output" << endl;
    cout << endl;
    cout << "    command: " << endl;
    cout << "       union       union of both input labels" << endl;
    cout << "       intersect   intersection of both input labels" << endl;
//    cout << "       invert label surface   compute inverse of label (need surface)" << endl;
    cout << endl;
    exit(1);
  }
  string comm = argv[1];
  string if1   = argv[2];
  string if2   = argv[3];
	string of    = argv[4];

	
	if (comm == "union")
	{
    LABEL *l1 = LabelRead(NULL,if1.c_str());
	  LABEL *l2 = LabelRead(NULL,if2.c_str());
	  assert (l1 != NULL);
	  assert (l2 != NULL);
	  LABEL * ret = LabelCombine(l1,l2);
		LabelRemoveDuplicates(ret);
		LabelWrite(ret,of.c_str());
	}
	else if (comm == "intersect")
	{
    LABEL *l1 = LabelRead(NULL,if1.c_str());
	  LABEL *l2 = LabelRead(NULL,if2.c_str());
	  assert (l1 != NULL);
	  assert (l2 != NULL);
		LabelIntersect(l1,l2);
		LabelWrite(l1,of.c_str());
	
	}
	else
	{
	  cerr << " Command: " << comm << " unknown !" << endl;
		exit(1);
  }

}
