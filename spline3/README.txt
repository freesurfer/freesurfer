
Martin Reuter
created Nov. 4 2014

Spline3

The class Spline3 allows for fast cubic spline interpolation and evaluation.
It can pre-compute and cache anything related to the x vector as well as
the xnew vector (of evaluation locations). Whenever a y vector comes in,
it can then quickly fit the spline and evaluate the spline at the xnew positions.


Spline3_test is a program to test the functionality. Both the class
and the test program are well documented.

Spline3 is independent of the rest of FreeSurfer and kept here for future 
updates. It is used in OCT projects.

