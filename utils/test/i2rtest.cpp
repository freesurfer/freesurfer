/**
 * @file  i2rtest.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 01:49:44 $
 *    $Revision: 1.6 $
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
// i2rtest.cpp
//
// Purpose:
// Testing routine to verify extract_i_to_r
//
// Requires:
//    rot0.mgh
//

#include <iostream>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

#include <mri.h>
#include <matrix.h>

#ifdef __cplusplus
}
#endif

// libutils needs Progname defined
char *Progname = "i2rtest";

using namespace std;

int main(int argc, char *argv[])
{
  string command;

  int width, height, depth;
  Real c_r, c_a, c_s;
  Real xsize, ysize, zsize;

  string filename;
  if (argc < 2)
    filename = "./rot0.mgh";
  else
    filename = argv[1];

  MRI *mri;
  mri = MRIread(const_cast<char *> ( filename.c_str() ));

  width = mri->width;
  height= mri->height;
  depth = mri->depth;

  xsize = mri->xsize;
  ysize = mri->ysize;
  zsize = mri->zsize;

  c_r = mri->c_r;
  c_a = mri->c_a;
  c_s = mri->c_s;

  cout << "input file = " << filename.c_str() << endl;
  cout << "Input data ------------------------------------------------" << endl;
  cout << "width= " << width << "  height= " << height << "  depth= " << depth << endl;
  cout << "xsize= " << xsize << "  ysize = " << ysize  << "  zsize= " << zsize << endl;
  cout << "c_r  = " << c_r   << "  c_a   = " << c_a    << "  c_s  = " << c_s << endl;
  cout << "ras_good_flag = " << mri->ras_good_flag << endl;
  cout << endl;

  MATRIX *m = extract_i_to_r(mri);
  VECTOR *c = VectorAlloc(4, MATRIX_REAL);
  c->rptr[1][1] = ((Real) width)/2.;
  c->rptr[2][1] = ((Real) height)/2.;
  c->rptr[3][1] = ((Real) depth)/2.;
  c->rptr[4][1] = 1.;

  // check the definition
  cout << "Check   C = M * c with c = (width/2, height/2, depth/2) " << endl;
  cout << "--------------------------------------------------------" << endl;
  VECTOR *C = MatrixMultiply(m, c, NULL);

  Real C_r = C->rptr[1][1];
  Real C_a = C->rptr[2][1];
  Real C_s = C->rptr[3][1];
  cout << "Calculated values are" << endl;
  cout << "C_r  = " << C_r   << "  C_a   = " << C_a    << "  C_s  = " << C_s << endl;

  MatrixFree(&m);
  VectorFree(&c);
  VectorFree(&C);
  MRIfree(&mri);

  double tolerance = 0.000001;
  double crdiff = C_r - c_r;
  double cadiff = C_a - c_a;
  double csdiff = C_s - c_s;
  if ( crdiff > tolerance || cadiff > tolerance || csdiff > tolerance )
  {
    cout << "***********************Error in extract_i_to_r()*********************** " << endl;
    cout << "cr diff = " << crdiff << "  ca diff = " << cadiff << "   cs diff = " << csdiff << endl;
    return -1;
  }
  cout << "No problem found." << endl;
  return 0;
}
