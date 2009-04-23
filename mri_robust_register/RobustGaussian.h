#ifndef RegPowell_H
#define RegPowell_H


class RobustGaussian
{
  public:
     static double median(double a[],int n);
     static double mad(double a[],int n, double d =  1.4826);
     static double kth_smallest(double a[], int n, int k);
     static double quick_select(double a[], int n, int k);
};


#endif
