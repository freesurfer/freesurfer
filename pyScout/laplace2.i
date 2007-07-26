%module Laplace
%{

#include "laplace2.h"

class Locator
{
public:
  Locator();
  void SetInputData(WrapperPolyData*);
  int FindClosestPoint(double,double,double);
};

%}

%include "std_vector.i"

namespace std {
  %template(DoubleVector) vector<double>;
  %template(WrapperLine)  vector< vector<double> >;
  %template(WrapperLineSet) vector< vector< vector<double> > >;
}

%ignore tracer.h;
%include "laplace2.h"

%include "WrapperPolyData.h"

class Locator
{
public:
  Locator();
  void SetInputData(WrapperPolyData*);
  int FindClosestPoint(double,double,double);
};
