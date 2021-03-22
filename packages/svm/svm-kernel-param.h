/**
 * @brief defines various kernel types and functionality.
 *
 * Name: KernelParam (abstract), LinearKernelParam, PolyKernelParam,
 *       RbfKernelParam, Kernel
 *
 * This file defines various kernel types and functionality. Different
 * KernelParam derived classes implement various types of kernels.
 * Kernel is the class that should be used in the user code. It can
 * initialize itself from a parameter string or a file, compute its
 * value given two inputs and differentiate itself with respect to the
 * first parameter.
 *
 * I/O format: type followed by the parameters. Right now, we have defined:
 *
 * "1" - linear
 * "2 d" - polynomial of degree d
 * "3 gamma" - rbf of width gamma (gamma = (2*signma^2))
 *
 * As new types are defines, Kernel:parse function should be exanded
 * to handle a new type.
 */
/*
 * Original Author: Polina Golland polina@ai.mit.edu
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

#ifndef __SVM_KERNEL_PARAM_H__
#define __SVM_KERNEL_PARAM_H__

#include <stdio.h>
#include <math.h>

#include "svm-element-type.h"
#include "svm-similarity.h"


const int NO_KERNEL     = 0;
const int LINEAR_KERNEL = 1;
const int POLY_KERNEL   = 2;
const int RBF_KERNEL    = 3;

class KernelParam {
protected:
  int _type;
  char _str[100];

public:
  KernelParam() : _type(NO_KERNEL) {};
  virtual ~KernelParam() {};

  // Value
  virtual double operator() (const SvmReal* v1, const SvmReal* v2, int n) const {
    return 1;
  };

  virtual double operator() (SvmReal dist) const {
    return 1;
  };


  // First order derivative functions: return either a particular component of
  // the derivative or a vector of all components.

  virtual SvmReal d10(int index, const SvmReal* x, const SvmReal* y, int n) const = 0;
  virtual void d10(SvmReal* res, const SvmReal* x, const SvmReal* y, int n) const = 0;


  // Second order derivatives, mixed: u_i v_j
  virtual SvmReal d11(int index1, int index2, const SvmReal* x, const SvmReal* y, int n) const = 0;

  // Accessors

  int type() const {
    return _type;
  }

  const char* getString() const {
    return _str;
  }

  // I/O

  void write (FILE *f) const {
    fprintf(f,"%s\n", getString());
  }
};



class LinearKernelParam: public KernelParam {
protected:

public:

  LinearKernelParam() {
    _type = LINEAR_KERNEL;
    sprintf(_str, "1 ");
  }

  // Value

  virtual double operator() (const SvmReal* v1, const SvmReal* v2, int n) const {
    return 1.0+Product(v1,v2,n);
  }

  virtual double operator() (SvmReal dist) const {
    return 1.0+dist;
  }


  // First order derivatives

  virtual SvmReal d10(int index, const SvmReal *x, const SvmReal* y, int n) const {
    return y[index];
  };

  virtual void d10(SvmReal* res, const SvmReal *x, const SvmReal* y, int n) const {
    for ( int i = 0; i < n; i++ )
      res[i] = y[i];
  }


  // Second order derivatives, mixed: u_i v_j
  virtual SvmReal d11(int index1, int index2, const SvmReal* x, const SvmReal* y, int n) const {
    return (index1==index2);
  }

};




class PolyKernelParam: public KernelParam {
  int _d;
public:

  PolyKernelParam(int d) {
    _type = POLY_KERNEL;
    _d = d;
    sprintf(_str,"%d %d ", _type, _d);
  }


  // Value

  virtual double operator() (const SvmReal* v1, const SvmReal* v2, int n) const {
    return pow(1.0+Product(v1,v2,n),_d);
  }

  virtual double operator() (SvmReal dist) const {
    return pow(1.0+dist,_d);
  };


  // First order derivatives

  virtual SvmReal d10(int index, const SvmReal *x, const SvmReal* y, int n) const {
    return y[index]*_d*pow(1+Product(x,y,n),_d-1);
  };

  virtual void d10(SvmReal* res, const SvmReal *x, const SvmReal* y, int n) const {
    double coef = _d*(1+Product(x,y,n),_d-1);
    for ( int i = 0; i < n; i++ )
      res[i] = (SvmReal)(y[i]*coef);
  };

  // Second order derivatives, mixed: u_i v_j
  virtual SvmReal d11(int index1, int index2, const SvmReal* x, const SvmReal* y, int n) const {
    double dataProduct = Product(x,y,n);
    return (SvmReal)(_d*(_d-1)*pow(1+dataProduct,_d-2)*x[index2]*y[index1]+
                     (index1==index2)*_d*pow(1+dataProduct,_d-1));
  }

};


class RbfKernelParam: public KernelParam {
  double _gamma;

public:

  RbfKernelParam(SvmReal gamma) {
    _type = RBF_KERNEL;
    _gamma = gamma;
    sprintf(_str,"%d %.12g ", _type, _gamma);
  }


  // Value

  virtual double operator() (const SvmReal* v1, const SvmReal* v2, int n) const {
    return exp(-Distance(v1,v2,n)/_gamma);
  }

  virtual double operator() (SvmReal dist) const {
    return exp(-dist/_gamma);
  }


  // First order derivatives

  virtual SvmReal d10(int index, const SvmReal *x, const SvmReal* y, int n) const {
    return (SvmReal)(-2/_gamma*operator()(x,y,n)*(x[index]-y[index]));
  }

  virtual void d10(SvmReal* res, const SvmReal *x, const SvmReal* y, int n) const {
    double val = -2*operator()(x,y,n)/_gamma;
    for ( int i = 0; i < n; i++ )
      res[i] = (SvmReal)((x[i]-y[i])*val);
  }

  // Second order derivatives, mixed: u_i v_j
  virtual SvmReal d11(int index1, int index2, const SvmReal* x, const SvmReal* y, int n) const {
    return (SvmReal)((index1==index2)*2/_gamma);
  }

};



#endif // __SVM_KERNEL_PARAM_H__







