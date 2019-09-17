/*
 * Original Author: Martin Reuter
 * Nov 4, 2014
 *
 */
#include "Spline3.h"
#include <math.h>


int main () 
{
   // Set 10 x, xnew and y values
   int N = 10;
   Spline3 <double > S;
   std::vector < double > x(N);
   std::vector < double > xnew(N);
   std::vector < double > ynew;
   double snm1 = sqrt(N-1);
   for (int i = 0; i<N; i++)
   {
     x[i] = sqrt(i); // sqrt spacing
     xnew[i] = snm1 *i/(N-1); // equal spacing
   }   
   static const double arr [] = {2,3,4,9,4,2,4,7,8,13}; // y data
   // Spline3 uses std::vector . Here is how to set if from an array:
   std::vector<double> y (arr, arr + sizeof(arr) / sizeof(arr[0]) );  
   // Attention, data is copied ! To avoid loss of speed, it is best if the
   // main program IO directly works with std::vector to avoid copying data.
   // You can then (if you need) get an array from the std::vector via, e.g.
   // double* px = &x[0];
   // without copying data. This array will be valid until the vector
   // gets deleted or re-allocates (when appending too-much data to it).
   
   // print data:
   std::cout << " x    : "; for (int i = 0;i<N;i++) std::cout << " " << x[i] ; std::cout << std::endl;
   std::cout << " xnew : "; for (int i = 0;i<N;i++) std::cout << " " << xnew[i] ;std::cout << std::endl;
   std::cout << " y    : "; for (int i = 0;i<N;i++) std::cout << " " << y[i] ;std::cout << std::endl;
 
   // Cache everything related to x spacing
   std::cout << "run preCacheX" << std::endl;
   S.preCacheX(x);

   // Cache everything related to xnew spacing
   std::cout << "run preCacheXnew" << std::endl;
   S.preCacheXnew(xnew);
   
   // In a loop with lots of different y data
   // First interpolate spline (x,y):
   std::cout << "run interp" << std::endl;
   S.interp(y); // using cached x

   // Then evaluate at xnew location
   // note, y needs to be passed again, as it is not
   // cached in interp above.   
   std::cout << "run eval" << std::endl;
   ynew = S.eval(y);  // using cached xnew
   
   // print interpolated results:
   std::cout << " Ynew : "; for (int i = 0;i<N;i++) std::cout << " " << ynew[i] ;std::cout << std::endl;
   
   
   // NOW test array interface and use float:
   float* yarr = (float*)malloc( N*sizeof(float));
   for (int i = 0;i<N;i++) yarr[i] = (float)arr[i];
   float* xarr = (float*)malloc( N*sizeof(float));
   for (int i = 0;i<N;i++) xarr[i] = (float)x[i];
   float* xnewarr = (float*)malloc( N*sizeof(float));
   for (int i = 0;i<N;i++) xnewarr[i] = (float)xnew[i];
   float *ynewarr = (float*)malloc(N*sizeof(float));

   Spline3 <float> S2;
   S2.preCacheX(N,xarr);
   S2.preCacheXnew(N,xnewarr); // could be different length M
   S2.interp(N,yarr); // using cached x
   // Then evaluate at xnew location
   // note, y needs to be passed again, as it is not
   // cached in interp above.   
   ynewarr = S2.eval(N,yarr,ynewarr);  // using cached xnew, ynewarr needs to be same length as xnewarr
   
   // print interpolated results:
   std::cout << " Ynewarr : "; for (int i = 0;i<N;i++) std::cout << " " << ynewarr[i] ;std::cout << std::endl;
   
   // free everything, but here we just quit
     
   return 0;
}
