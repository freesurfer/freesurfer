/**
 * @file  mri_gradient_info.cpp
 * @brief A programm to compute gradient information
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2010/01/14 22:55:17 $
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
#include "histo.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

#ifdef __cplusplus
}
#endif

using namespace std;

//static char vcid[] = "$Id: mri_gradient_info.cpp,v 1.1 2010/01/14 22:55:17 mreuter Exp $";
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

  if (argc < 2)
  {
    cout << endl;
    cout << argv[0] << " image.mgz" << endl;
    cout << endl;
//    cout << "    norm-div  (=1)  divide final distance by this (e.g. step adjustment)" << endl;
//    cout << "    dist-type " << endl;
//    cout << "       1  (default) Rigid Transform Distance (||log(R)|| + ||T||)" << endl;
//    cout << "       2            Affine Transform Distance (RMS) " << endl;
//    cout << "       3            8-corners mean distance after transform " << endl;
//		cout << "    invert1         1 true, 0 false (default)" << endl;
//    cout << endl;
    exit(1);
  }
  string mrif = argv[1];

  MRI* mri_in = MRIread(mrif.c_str());
	
	MRI* mri_mag  = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT);
	MRI* mri_grad = MRIsobel(mri_in, NULL, mri_mag);
	MRIfree(&mri_grad);
	
	int n = 5;
  HISTOGRAM *h = MRIhistogram(mri_mag,n);
	
	cout << h->counts[n-1] << endl;	

}
