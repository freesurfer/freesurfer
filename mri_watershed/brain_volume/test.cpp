/**
 * @file  test.cpp
 * @brief tests TVector
 *
 */
/*
 * Original Author: Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/02/01 18:16:43 $
 *    $Revision: 1.1 $
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
 *
 */


#include "TVector.h"
#include <iostream>

using namespace std;

void callfunc(int i, void (*func)(const TVector &a)) {
  TVector a(1,2,3);
  func(a);
}

void calc(const TVector &a) {
  cout << a << endl;
  double c = 5.0;
  cout << "\t c*a = " << c*a << endl;
  cout << "\t a*c = " << a*c << endl;
}

int main() {
  TVector VA(1,1,1);
  TVector VB(2,1,1);
  TVector VC(3,1,1);

  TVector v = VA + VB;
  cout << "Testing the sum of two vectors" << endl;
  cout << "\t VA " << VA << endl;
  cout << "\t VB " << VB << endl;
  cout << "\t VA+VB " << v << endl;

  cout << "Multiplying double and vector" << endl;
  double c = 5.0;
  double d = 6.0;

  cout << "\t VA = " << VA << "  c= " << c << endl;
  cout << "\t VA*c " << VA*c << endl;
  cout << "\t VA = " << VA << "  c= " << c << endl;
  cout << "\t c*VA " << c*VA << endl;
  cout << "Testing == operator" << endl;
  cout << "\tusing VA " << VA << ", VB " << VB << endl;
  cout << "\t check VA==VB (true=1, false = 0) : " << (VA==VB) << endl;

  cout << "Inner product and outer product of two vectors " << endl;
  cout << "\t VA*VB " <<  VA*VB << " (VA^VB) " << (VA^VB) << endl;

  cout << "Testing the combination" << endl;
  cout << "\tVA = " << VA << " c = " << c << "  d = " << d << endl;
  cout << "\t\t\tc*d*VA = " << c*d*VA << endl;
  cout << "\tVA = " << VA << " c = " << c << "  d = " << d << endl;
  cout << "\t\t\tc*VA*d = " << c*VA*d << endl;
  cout << "\tVA = " << VA << " c = " << c << "  d = " << d << endl;
  cout << "\t\t\tVA*c*d = " << VA*c*d << endl;
  cout << "\tVA = " << VA << " c = " << c << "  d = " << d << endl;

  cout << "Testing interesting combination" << endl;
  cout << "\tVA = " << VA << "  VB = " << VB << "  VC = " << VC << endl;
  cout << "\t\t\tVA*VB = " << VA*VB << endl;
  cout << "\t\t\tVB*VC = " << VB*VC << endl;
  cout << "\t\t\tVA*VB*VC = " << VA*VB*VC << endl;
  cout << "\t\t\tVA*(VB*VC) = " << VA*(VB*VC) << endl;
  cout << "\tVA = " << VA << " VB = " << VB << endl;
  cout << "\t\t\t(VA*VB)*VC = " << (VA*VB)*VC << endl;

  cout << "Testing the composite operator" << endl;
  cout << "\tVA = " << VA << "  c = " << c << endl;
  cout << "\t\t\t(VA*=c)" << (VA*=c) << endl;
  cout << "\tVA = " << VA << "  c = " << c << endl;
  cout << "\t\t\t(VA/=c) " << (VA/=c) << endl;

  cout << "\tVA = " << VA << "  VB = " << VB << endl;
  cout << "\t\t\t(VA+=VB) " << (VA+=VB) << endl;
  cout << "\tVA = " << VA << "  VB = " << VB << endl;
  cout << "\t\t\t(VA-=VB) " << (VA-=VB) << endl;
  cout << "\tVA = " << VA << "  VB = " << VB << endl;

  cout << "\nTesting the functional with TVector argument" << endl;
  callfunc(1, calc);

  return 0;
}
