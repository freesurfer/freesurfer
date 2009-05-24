#include "RobustGaussian.h"

#include <cmath>
#include <iostream>

using namespace std;

/*
  * This Quickselect routine is based on the algorithm described in
  * "Numerical recipes in C", Second Edition,
  * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  */
#define ELEM_SWAPD(a,b) { double t=a;a=b;b=t; }
#define ELEM_SWAPI(a,b) { int t=a;a=b;b=t; }
pair < double , int > RobustGaussian::quick_selectI(double arr[], int n, int k)
{
  int low, high ;
  int median;
  int middle, ll, hh;

  int* pos = (int *)calloc(n, sizeof(int));
  for (int i = 0;i<n;i++) pos[i] = i;

  low = 0 ;
  high = n-1 ;
  median = k-1;
  for (;;)
  {
    if (high <= low) /* One element only */
    {
      int pp = pos[median];
      free(pos);
      return pair < double , int > (arr[median],pp) ;
    }
    if (high == low + 1)
    { /* Two elements only */
      if (arr[low] > arr[high])
      {
        ELEM_SWAPD(arr[low], arr[high]) ;
	ELEM_SWAPI(pos[low], pos[high]) ;
      }
      int pp = pos[median];
      free(pos);
      return pair < double , int > (arr[median],pp)  ;
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])
    {
      ELEM_SWAPD(arr[middle], arr[high]) ;
      ELEM_SWAPI(pos[middle], pos[high]) ;
    }  
    if (arr[low] > arr[high])
    {
      ELEM_SWAPD(arr[low], arr[high]) ;
      ELEM_SWAPI(pos[low], pos[high]) ;
    }
    if (arr[middle] > arr[low])
    {
      ELEM_SWAPD(arr[middle], arr[low]) ;
      ELEM_SWAPI(pos[middle], pos[low]) ;
    }
      
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAPD(arr[middle], arr[low+1]) ;
    ELEM_SWAPI(pos[middle], pos[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;)
    {
      do ll++;
      while (arr[low] > arr[ll]) ;
      do hh--;
      while (arr[hh] > arr[low]) ;
      if (hh < ll)
        break;
      ELEM_SWAPD(arr[ll], arr[hh]) ;
      ELEM_SWAPI(pos[ll], pos[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAPD(arr[low], arr[hh]) ;
    ELEM_SWAPI(pos[low], pos[hh]) ;
    /* Re-set active partition */
    if (hh <= median)
      low = ll;
    if (hh >= median)
      high = hh - 1;
  }
}

double RobustGaussian::quick_select(double arr[], int n, int k)
{
  int low, high ;
  int median;
  int middle, ll, hh;


  low = 0 ;
  high = n-1 ;
  median = k-1;
  for (;;)
  {
    if (high <= low) /* One element only */
    {
      return arr[median] ;
    }
    if (high == low + 1)
    { /* Two elements only */
      if (arr[low] > arr[high])
      {
        ELEM_SWAPD(arr[low], arr[high]) ;
      }
      return arr[median] ;
    }
    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])
    {
      ELEM_SWAPD(arr[middle], arr[high]) ;
    }  
    if (arr[low] > arr[high])
    {
      ELEM_SWAPD(arr[low], arr[high]) ;
    }
    if (arr[middle] > arr[low])
    {
      ELEM_SWAPD(arr[middle], arr[low]) ;
    }
      
    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAPD(arr[middle], arr[low+1]) ;
    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;)
    {
      do ll++;
      while (arr[low] > arr[ll]) ;
      do hh--;
      while (arr[hh] > arr[low]) ;
      if (hh < ll)
        break;
      ELEM_SWAPD(arr[ll], arr[hh]) ;
    }
    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAPD(arr[low], arr[hh]) ;
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
std::pair < double , int> RobustGaussian::kth_smallestI(double a[], int n, int k)
// reorders a[] (of length n)
// returns k_th smallest and original position in array
{
  int i,j,l,m ;
  int kk = k-1;
  double x ;
  
  int* pos = (int *)calloc(n, sizeof(int));
  for (int i = 0;i<n;i++) pos[i] = i;
  
  l=0 ;
  m=n-1 ;
  while (l<m)
  {
    x=a[kk] ;
    i=l ;
    j=m ;
    do
    {
      while (a[i]<x) i++ ;
      while (x<a[j]) j-- ;
      if (i<=j)
      {
        ELEM_SWAPD(a[i],a[j]) ;
	ELEM_SWAPI(pos[i],pos[j]);
        i++ ;
        j-- ;
      }
    }
    while (i<=j) ;
    if (j<kk) l=i ;
    if (kk<i) m=j ;
  }
  int pp = pos[kk];
  free(pos);
  return pair < double, int > (a[kk],pp) ;
}

double RobustGaussian::kth_smallest(double a[], int n, int k)
// reorders a[] (of length n)
// returns k_th smallest and original position in array
{
  int i,j,l,m ;
  int kk = k-1;
  double x ;
    
  l=0 ;
  m=n-1 ;
  while (l<m)
  {
    x=a[kk] ;
    i=l ;
    j=m ;
    do
    {
      while (a[i]<x) i++ ;
      while (x<a[j]) j-- ;
      if (i<=j)
      {
        ELEM_SWAPD(a[i],a[j]) ;
        i++ ;
        j-- ;
      }
    }
    while (i<=j) ;
    if (j<kk) l=i ;
    if (kk<i) m=j ;
  }
  return a[kk];
}


#undef ELEM_SWAPD
#undef ELEM_SWAPI

pair < double, double > RobustGaussian::medianI(double t[],int n)
// compute median
// in situ, t will be reordered
// index (return.second ) is double because with even no of elements, we lie between two indices
{

  pair < double, int > qs;
  if (n%2 == 1) //odd
  {
    qs = kth_smallestI(t,n,(n+1)/2);
    //cout << " n: " << n << "   " << qs << endl;
    //free(t);
    return pair < double , double > (qs.first,(double)qs.second);
  }

  //  else even:

//  qs = kth_smallest(t,n,n/2);
  qs = quick_selectI(t,n,n/2);
// double qs2 = kth_smallest(t,n,n/2 + 1);
  pair < double , int > qs2 = quick_selectI(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << qs << "   " << qs2 << endl;

  return pair < double, double > (0.5 * (qs.first + qs2.first), 0.5 * (qs.second + qs2.second));

}

double RobustGaussian::median(double t[],int n)
// compute median
// in situ, t will be reordered
{

  double q;
  if (n%2 == 1) //odd
  {
    q = kth_smallest(t,n,(n+1)/2);
    return q;
  }

  //  else even:

//  q = kth_smallest(t,n,n/2);
  q = quick_select(t,n,n/2);
// double q2 = kth_smallest(t,n,n/2 + 1);
  double q2 = quick_select(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << q << "   " << q2 << endl;

  return 0.5 * (q + q2);

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
