/**
 * @brief tests TVector
 *
 */
/*
 * Original Author: Tosa
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

#include "TVector.h"

void callfunc(int i, void (*func)(const TVector &a)) {
  TVector a(1, 2, 3);
  func(a);
}

void calc(const TVector &a) {
  std::cout << a << std::endl;
  double c = 5.0;
  std::cout << "\t c*a = " << c * a << std::endl;
  std::cout << "\t a*c = " << a * c << std::endl;
}

int main() {
  TVector VA(1, 1, 1);
  TVector VB(2, 1, 1);
  TVector VC(3, 1, 1);

  TVector v = VA + VB;
  std::cout << "Testing the sum of two vectors" << std::endl;
  std::cout << "\t VA " << VA << std::endl;
  std::cout << "\t VB " << VB << std::endl;
  std::cout << "\t VA+VB " << v << std::endl;

  std::cout << "Multiplying double and vector" << std::endl;
  double c = 5.0;
  double d = 6.0;

  std::cout << "\t VA = " << VA << "  c= " << c << std::endl;
  std::cout << "\t VA*c " << VA * c << std::endl;
  std::cout << "\t VA = " << VA << "  c= " << c << std::endl;
  std::cout << "\t c*VA " << c * VA << std::endl;
  std::cout << "Testing == operator" << std::endl;
  std::cout << "\tusing VA " << VA << ", VB " << VB << std::endl;
  std::cout << "\t check VA==VB (true=1, false = 0) : " << (VA == VB)
            << std::endl;

  std::cout << "Inner product and outer product of two vectors " << std::endl;
  std::cout << "\t VA*VB " << VA * VB << " (VA^VB) " << (VA ^ VB) << std::endl;

  std::cout << "Testing the combination" << std::endl;
  std::cout << "\tVA = " << VA << " c = " << c << "  d = " << d << std::endl;
  std::cout << "\t\t\tc*d*VA = " << c * d * VA << std::endl;
  std::cout << "\tVA = " << VA << " c = " << c << "  d = " << d << std::endl;
  std::cout << "\t\t\tc*VA*d = " << c * VA * d << std::endl;
  std::cout << "\tVA = " << VA << " c = " << c << "  d = " << d << std::endl;
  std::cout << "\t\t\tVA*c*d = " << VA * c * d << std::endl;
  std::cout << "\tVA = " << VA << " c = " << c << "  d = " << d << std::endl;

  std::cout << "Testing interesting combination" << std::endl;
  std::cout << "\tVA = " << VA << "  VB = " << VB << "  VC = " << VC
            << std::endl;
  std::cout << "\t\t\tVA*VB = " << VA * VB << std::endl;
  std::cout << "\t\t\tVB*VC = " << VB * VC << std::endl;
  std::cout << "\t\t\tVA*VB*VC = " << VA * VB * VC << std::endl;
  std::cout << "\t\t\tVA*(VB*VC) = " << VA * (VB * VC) << std::endl;
  std::cout << "\tVA = " << VA << " VB = " << VB << std::endl;
  std::cout << "\t\t\t(VA*VB)*VC = " << (VA * VB) * VC << std::endl;

  std::cout << "Testing the composite operator" << std::endl;
  std::cout << "\tVA = " << VA << "  c = " << c << std::endl;
  std::cout << "\t\t\t(VA*=c)" << (VA *= c) << std::endl;
  std::cout << "\tVA = " << VA << "  c = " << c << std::endl;
  std::cout << "\t\t\t(VA/=c) " << (VA /= c) << std::endl;
  std::cout << "\tVA = " << VA << "  VB = " << VB << std::endl;
  std::cout << "\t\t\t(VA+=VB) " << (VA += VB) << std::endl;
  std::cout << "\tVA = " << VA << "  VB = " << VB << std::endl;
  std::cout << "\t\t\t(VA-=VB) " << (VA -= VB) << std::endl;
  std::cout << "\tVA = " << VA << "  VB = " << VB << std::endl;

  std::cout << "\nTesting the functional with TVector argument" << std::endl;
  callfunc(1, calc);

  return 0;
}
