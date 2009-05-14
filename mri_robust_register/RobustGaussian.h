#ifndef RegPowell_H
#define RegPowell_H

#include <utility>


class RobustGaussian
{
public:
  static double median(double a[],int n){return medianI(a,n).first;};
  static double kth_smallest(double a[], int n, int k){return kth_smallestI(a,n,k).first;};
  static double quick_select(double a[], int n, int k){return quick_selectI(a,n,k).first;};
  static std::pair <double, double> medianI(double a[],int n);
  static std::pair <double, int> kth_smallestI(double a[], int n, int k);
  static std::pair <double, int> quick_selectI(double a[], int n, int k);
  static double mad(double a[],int n, double d =  1.4826);
};


#endif
