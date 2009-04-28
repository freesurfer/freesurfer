#include "RobustGaussian.h"

#include <cmath>
#include <iostream>

using namespace std;

/*
  * This Quickselect routine is based on the algorithm described in
  * "Numerical recipes in C", Second Edition,
  * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  */
#define ELEM_SWAP(a,b) { double t=a;a=b;b=t; }
double RobustGaussian::quick_select(double arr[], int n, int k)
{
     int low, high ;
     int median;
     int middle, ll, hh;
     low = 0 ; high = n-1 ; median = k-1;
     for (;;)
     {
         if (high <= low) /* One element only */
             return arr[median] ;
         if (high == low + 1)
	 { /* Two elements only */
             if (arr[low] > arr[high])
                 ELEM_SWAP(arr[low], arr[high]) ;
             return arr[median] ;
         }
     /* Find median of low, middle and high items; swap into position low */
     middle = (low + high) / 2;
     if (arr[middle] > arr[high])     ELEM_SWAP(arr[middle], arr[high]) ;
     if (arr[low] > arr[high])        ELEM_SWAP(arr[low], arr[high]) ;
     if (arr[middle] > arr[low])      ELEM_SWAP(arr[middle], arr[low]) ;
     /* Swap low item (now in position middle) into position (low+1) */
     ELEM_SWAP(arr[middle], arr[low+1]) ;
     /* Nibble from each end towards middle, swapping items when stuck */
     ll = low + 1;
     hh = high;
     for (;;)
     {
         do ll++; while (arr[low] > arr[ll]) ;
         do hh--; while (arr[hh] > arr[low]) ;
         if (hh < ll)
         break;
         ELEM_SWAP(arr[ll], arr[hh]) ;
     }
     /* Swap middle item (in position low) back into correct position */
     ELEM_SWAP(arr[low], arr[hh]) ;
     /* Re-set active partition */
     if (hh <= median)
         low = ll;
     if (hh >= median)
        high = hh - 1;
    }
}

/*---------------------------------------------------------------------------
   Function :   kth_smallest()
   In       :   array of elements, # of elements in the array, rank k
   Out      :   one element
   Job      :   find the kth smallest element in the array

                Reference:

                  Author: Wirth, Niklaus 
                   Title: Algorithms + data structures = programs 
               Publisher: Englewood Cliffs: Prentice-Hall, 1976 
    Physical description: 366 p. 
                  Series: Prentice-Hall Series in Automatic Computation 

 ---------------------------------------------------------------------------*/
double RobustGaussian::kth_smallest(double a[], int n, int k)
{
    int i,j,l,m ;
    int kk = k-1;
    double x ;

    l=0 ; m=n-1 ;
    while (l<m) {
        x=a[kk] ;
        i=l ;
        j=m ;
        do {
            while (a[i]<x) i++ ;
            while (x<a[j]) j-- ;
            if (i<=j) {
                ELEM_SWAP(a[i],a[j]) ;
                i++ ; j-- ;
            }
        } while (i<=j) ;
        if (j<kk) l=i ;
        if (kk<i) m=j ;
    }
    return a[kk] ;
}

#undef ELEM_SWAP

double RobustGaussian::median(double t[],int n)
// compute median
// in situ, t will be reordered
{
    
  double qs;
  if (n%2 == 1) //odd
  {
    qs = kth_smallest(t,n,(n+1)/2);
    //cout << " n: " << n << "   " << qs << endl;
    //free(t);
    return qs;
   }

  //  else even:

//  qs = kth_smallest(t,n,n/2);
  qs = quick_select(t,n,n/2);
 // double qs2 = kth_smallest(t,n,n/2 + 1);
  double qs2 = quick_select(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << qs << "   " << qs2 << endl;
  qs =  0.5 * ( qs + qs2);
  return qs;
 
}

double RobustGaussian::mad(double a[], int n, double d)
// robust estimate for sigma (using median absolute deviation)
// 1.4826 med_i(|r_i - med_j(r_j)|)
// array a will be reordered!
{
  double medi = median(a,n);
  //cout << " median: " << medi << endl;
  double* t = (double *)calloc(n, sizeof(double)); 
  for (int i=0;i<n;i++)
  {
     t[i] = fabs(a[i] -medi);
     //cout  << t[i] << " " << flush;
  }
  //cout <<endl;
  double mm = median(t,n);
  //cout << " mmedian: " << mm << endl;
  free(t);
  return d * mm;
}
