/**
 * @file RobustGaussian.cpp
 * @brief A class to esimate a robust Gaussian (using median and mad)
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2011/08/23 18:53:41 $
 *    $Revision: 1.12.2.1 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "RobustGaussian.h"

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#ifdef __cplusplus
extern "C"
{
#endif
   #include "error.h"
#ifdef __cplusplus
}
#endif

using namespace std;

/*
  * This Quickselect routine is based on the algorithm described in
  * "Numerical recipes in C", Second Edition,
  * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  * This code is based on code by Nicolas Devillard - 1998. Public domain.  
  * see also http://ndevilla.free.fr/median/median.pdf
  * or http://ndevilla.free.fr/median/median/index.html
  * modifications: - instead of only selecting the median, select the k-th smallest
  *                - additionally keep track of the index in the original array
  */
#define ELEM_SWAPD(a,b) { T t=a;a=b;b=t; }
#define ELEM_SWAPI(a,b) { int t=a;a=b;b=t; }
template <class T>
pair < T , int > RobustGaussian<T>::quick_selectI(T arr[], int n, int k)
{
  int low, high ;
  int median;
  int middle, ll, hh;

  int* pos = (int *)calloc(n, sizeof(int));
  if (pos == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"RobustGaussian<T>::quick_selectI could not allocate memory for pos") ;
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
      return pair < T , int > (arr[median],pp) ;
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
      return pair < T , int > (arr[median],pp)  ;
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

template <class T>
T RobustGaussian<T>::quick_select(T arr[], int n, int k)
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

    Implementation based on code by N.Devillard. Public Domain.

    kth_smallestI is modified: we keep track of the position inside the original 
    array of the k_th smallest element

 ---------------------------------------------------------------------------*/
template <class T>
std::pair < T , int> RobustGaussian<T>::kth_smallestI(T a[], int n, int k)
// reorders a[] (of length n)
// returns k_th smallest and original position in array
{
  int i,j,l,m ;
  int kk = k-1;
  T x ;
  
  int* pos = (int *)calloc(n, sizeof(int));
  if (pos == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"RobustGaussian<T>::kth_smallestI could not allocate memory for pos") ;
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
  return pair < T, int > (a[kk],pp) ;
}

template <class T>
T RobustGaussian<T>::kth_smallest(T a[], int n, int k)
// reorders a[] (of length n)
// returns k_th smallest and original position in array
{
  int i,j,l,m ;
  int kk = k-1;
  T x ;
    
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

template <class T>
pair < T, T > RobustGaussian<T>::medianI(T t[],int n)
// compute median
// in situ, t will be reordered
// index (return.second ) is T (double or float) because with even no of elements, we lie between two indices
{

  pair < T, int > qs;
  if (n%2 == 1) //odd
  {
    qs = kth_smallestI(t,n,(n+1)/2);
    //cout << " n: " << n << "   " << qs << endl;
    //free(t);
    return pair < T , T > (qs.first,( T )qs.second);
  }

  //  else even:

//  qs = kth_smallest(t,n,n/2);
  qs = quick_selectI(t,n,n/2);
// double qs2 = kth_smallest(t,n,n/2 + 1);
  pair < T , int > qs2 = quick_selectI(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << qs << "   " << qs2 << endl;

  return pair < T, T > (0.5 * (qs.first + qs2.first), 0.5 * (qs.second + qs2.second));

}

template <class T>
T RobustGaussian<T>::median(T t[],int n)
// compute median
// in situ, t will be reordered
{

  T q;
  if (n%2 == 1) //odd
  {
    q = kth_smallest(t,n,(n+1)/2);
    return q;
  }

  //  else even:

//  q = kth_smallest(t,n,n/2);
  q = quick_select(t,n,n/2);
// double q2 = kth_smallest(t,n,n/2 + 1);
  T q2 = quick_select(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << q << "   " << q2 << endl;

  return 0.5 * (q + q2);

}

template <class T>
void mmm (T a[], int n)
{

	 if (n <= 0) return;
   T min = a[0];
	 T max = a[0];
	 double mean = 0.0;
   for (int i = 1;i<n;i++)
	 {
	    if (a[i] < min) min = a[i];
			if (a[i] > max) max = a[i];
			mean += a[i];
	 }
	 mean /= n;
	 cout << " min: " << min << "  max: " << max << "  mean: " << mean << endl;
}

template <class T>
T RobustGaussian<T>::mad(T a[], int n, T d)
// robust estimate for sigma (using median absolute deviation)
// 1.4826 med_i(|r_i - med_j(r_j)|)
// array a will be reordered!
{
  T medi = median(a,n);
//	mmm(a,n);
//	cout << " median: " << medi << endl;
	
  T* t = (T *)calloc(n, sizeof(T));
  if (t == NULL) 
     ErrorExit(ERROR_NO_MEMORY,"RobustGaussian<T>::mad could not allocate memory for t") ;
     
//	double min = fabs(a[0] -medi);
//	double max = fabs(a[0] -medi);
//	double mean = 0.0;
  for (int i=0;i<n;i++)
  {
    t[i] = fabs(a[i] -medi);
//		if (t[i] < min) min = t[i];
//		if (t[i] > max) max = t[i];
//		mean += t[i];
//    //cout  << t[i] << " " << flush;
  }
//	mean /= n;
//	cout << " min: " << min << "  max: " << max << "  mean: " << mean << endl;

  T mm = median(t,n);
  //cout << endl <<" RobustGaussian<T>::mad median: " << mm << endl;
  free(t);
  return d * mm;
}
