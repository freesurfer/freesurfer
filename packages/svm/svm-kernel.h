/*
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


////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: KernelParam (abstract), LinearKernelParam, PolyKernelParam,
//       RbfKernelParam, Kernel
//
// This file defines various kernel types and functionality. Different
// KernelParam derived classes implement various types of kernels.
// Kernel is the class that should be used in the user code. It can
// initialize itself from a parameter string or a file, compute its
// value given two inputs and differentiate itself with respect to the
// first parameter.
//
// I/O format: type followed by the parameters. Right now, we have defined:
//
// "1" - linear
// "2 gamma" - rbf of width (2*signma^2) gamma
// "3 d" - polynomial of degree d
//
//  As new types are defines, Kernel:parse function should be expanded
//  to handle a new type.
//
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////

#ifndef __SVM_KERNEL_H__
#define __SVM_KERNEL_H__

#include <stdio.h>
#include <iostream>
#include <math.h>

#include "svm-kernel-param.h"


class Kernel {
  KernelParam* _p;

  // This function checks if the pointer has been set and prints an error
  // message if it's still null.
  bool isDefined() const {
    if ( _p == NULL ) {
      std::cerr << "Kerel error: attempting to use a kernel that has not "
      << "been defined yet.\n";
      return false;
    }
    return true;
  }


public:

  // Constructors and destructors
  Kernel(): _p(NULL) {};

  Kernel(const char *paramString):  _p(NULL) {
    parse(paramString);
  }

  Kernel(const Kernel& kernel): _p(NULL) {
    parse(kernel.getString());
  }

  Kernel& operator = (const Kernel& kernel) {
    parse(kernel.getString());
    return *this;
  }

  ~Kernel() {
    delete _p;
  };




  // Kernel type
  int type() const {
    if ( !isDefined() )
      return NO_KERNEL;
    return _p->type();
  }

  // Compute value
  double operator() (const SvmReal* v1, const SvmReal* v2, int n) const {
    if ( !isDefined() )
      return 0;

    return _p->operator()(v1,v2,n);
  }

  double operator() (SvmReal dist) const {
    if ( !isDefined() )
      return 0;
    return _p->operator()(dist);
  }


  // Compute derivatives

  SvmReal d10(int index, const SvmReal* x, const SvmReal* y, int n) const {
    if ( !isDefined() )
      return 0;
    return _p->d10(index,x,y,n);
  }

  bool d10(SvmReal* res, const SvmReal* x, const SvmReal* y, int n) const {
    if ( !isDefined() )
      return false;
    _p->d10(res,x,y,n);
    return true;
  }



  // I/O

  bool parse (const char* paramString);
  const char* getString() const {
    if ( !isDefined() )
      return "No kernel";
    return _p->getString();
  }

  bool read (FILE *f, bool binary = false);
  bool write (FILE *f, bool binary = false) const {
    if ( !isDefined() )
      return false;
    _p->write(f);
    return true;
  }
};




#endif // __SVM_KERNEL_H__







