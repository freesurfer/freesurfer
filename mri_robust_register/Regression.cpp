
#include "Regression.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

MATRIX * Regression::getRobustEst(double sat, double sig)
{
   pair < MATRIX* , MATRIX* > pw = getRobustEstW(sat,sig);
   MatrixFree(&pw.second);
   return pw.first;
}

// void printusage ()
// {
//  char buf[30];
//         snprintf(buf, 30, "/proc/%u/statm", (unsigned)getpid());
//         FILE* pf = fopen(buf, "r");
//         if (pf) {
//             unsigned size; //       total program size
//             //unsigned resident;//   resident set size
//             //unsigned share;//      shared pages
//             //unsigned text;//       text (code)
//             //unsigned lib;//        library
//             //unsigned data;//       data/stack
//             //unsigned dt;//         dirty pages (unused in Linux 2.6)
//             //fscanf(pf, "%u" /* %u %u %u %u %u"*/, &size/*, &resident, &share, &text, &lib, &data*/);
//             fscanf(pf, "%u" , &size);
//             //DOMSGCAT(MSTATS, std::setprecision(4) << size / (1024.0) << "MB mem used");
// 	    cout <<  size / (1024.0) << " MB mem used" << endl;
//         }
//         fclose(pf);
// while (1)
// {
//     if ('n' == getchar())
//        break;
// }	
// }

pair < MATRIX *, MATRIX *> Regression::getRobustEstW(double sat, double sig)
// solves overconstrained system A p = b using
// M estimators (Beaton and Tukey's biweigt function)
// returns pair < p, w  > of vectors
// where p is the robust solution and w the used weights
// Method: iterative reweighted least squares
{
//   cout << "  Regression::getRobustEstW( "<<sat<<" , "<<sig<<" ) " << endl;

  // constants
  int MAXIT = 100;
  //double EPS = 2e-16;
  double EPS = 2e-12;
  
  // variables
  vector < double > err(MAXIT);
  err[0] = numeric_limits<double>::infinity();
  err[1] = 1e20;
  int count = 1, rr, cc;
  double sigma;
  
  MATRIX * w = MatrixConstVal(1.0,A->rows,1,NULL);
  MATRIX * r = MatrixConstVal(1.0,A->rows,1,NULL);
  MATRIX * p = MatrixZero(A->cols, 1, NULL);
  
  MATRIX * lastp = MatrixCopy(p,NULL);
  MATRIX * lastw = MatrixCopy(w,NULL);
  MATRIX * lastr = MatrixCopy(B,NULL);  // error lastr = b-A*p = b  , since p = 0 initially

  MATRIX * wAi = NULL, *v = NULL;
  MATRIX * wA = MatrixAlloc(A->rows, A->cols, MATRIX_REAL);
  MATRIX * wr = MatrixAlloc(A->rows, 1, MATRIX_REAL);

  // compute error (lastr = b-A*p)
  //v     = MatrixMultiply(A,p,NULL);
  //lastr = MatrixSubtract(B,v,lastr);
  //MatrixFree(&v); v = NULL;
  
  // iteration until we increase the error, we reach maxit or we have no error
  while ( err[count-1] > err[count] && count < MAXIT && err[count] > EPS )
  {
    // save previous values   (not necessary for first step)  
    MatrixCopy(p, lastp);
    MatrixCopy(w, lastw);
    MatrixCopy(r, lastr);

    count++;
    
    //cout << " count: "<< count << endl;
        
    // recompute weights
    if (count > 2)
    {
      sigma = getSigmaMAD(lastr);
      //cout << " sigma: " << sigma << endl;
      MatrixScalarMul(lastr,(float)1.0/sigma,r); // use r temporarily here
      getSqrtTukeyDiaWeights(r, sat, w); // here we get sqrt of weights into w
    }
    
    
    // solve WLS with the new weights 
    // compute wA  where w = sqrt(W);
    for (rr = 1;rr<=wA->rows;rr++)
    for (cc = 1;cc<=wA->cols;cc++)
  	wA->rptr[rr][cc] = A->rptr[rr][cc] * w->rptr[rr][1];

    // compute wr  = ((b-A*p).*w);
    for (rr = 1;rr<=wr->rows;rr++)
  	wr->rptr[rr][1] = lastr->rptr[rr][1] * w->rptr[rr][1];

    // compute new p = lastp + pinv(wA)*wr
    //wAi = MatrixPseudoInverse(wA,wAi);
    if (wAi != NULL) MatrixFree(&wAi);
    wAi = MatrixSVDPseudoInverse(wA,NULL);
    
   if (wAi == NULL)
  {
     cerr << "    Regression::getRobustEstW  could not compute pseudo inverse!" << endl;
     //cerr << " wa: " << endl ;
      //   MatrixPrintFmt(stdout,"% 2.8f",wAi);
	//cerr << endl;
     assert(wAi != NULL);
  }
   //cout << " got inverse" << endl;
   // MatrixPrintFmt(stdout,"% 2.8f",wAi);
   // cout << endl;
   // MatrixPrintFmt(stdout,"% 2.8f",wr);
   // cout << endl;

    v   = MatrixMultiply(wAi,wr,v);
   // cout << " asdf" << endl;
    MatrixAdd(lastp,v,p);
    
        
    // compute new residual and total errors (using new p)    
    // r = b - (A*p)
    MatrixMultiply(A,p,r);
    MatrixSubtract(B,r,r);
    // err = sum (w r^2) / sum (w)
    float swr = 0;
    float sw  = 0;
    int rr;
    for (rr = 1;rr<=r->rows;rr++)
    {
      float t1 = w->rptr[rr][1];
      float t2 = r->rptr[rr][1];
      t1 *= t1; // remember w is the sqrt of the weights
      t2 *= t2;
      sw  += t1;
      swr += t1*t2;
    }
    err[count] = swr/sw;
    
//    cout << " Step: " << count-1 << " ERR: "<< err[count]<< endl;
//    cout << " p  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",p);
//    cout << endl;
       
   // do not save values here, to be able to use lastp when exiting earlier due to increasing error
   // MatrixCopy(p, lastp);
   // MatrixCopy(w, lastw);
   // MatrixCopy(r, lastr);
    
  }
  
  if (err[count]>err[count-1])
  {
   // take previous values (since actual values made the error to increase)
   // cout << " last step was no improvement, taking values : "<<  count-1 << endl;
    MatrixFree(&p);
    MatrixFree(&w);
    p = lastp;
    w = lastw;
    cout << "     Step: " << count-2 << " ERR: "<< err[count-1]<< endl;
    lasterror = err[count-1];
  }
  else
  {
    MatrixFree(&lastp);
    MatrixFree(&lastw);
    cout << "     Step: " << count-1 << " ERR: "<< err[count]<< endl;
    lasterror = err[count];
  } 
  
  MatrixFree(&lastr);
  MatrixFree(&r);
  MatrixFree(&wA);
  MatrixFree(&wr);
  MatrixFree(&wAi);
  MatrixFree(&v);

//    cout << " p-final  : "<< endl;
//    MatrixPrintFmt(stdout,"% 2.8f",p);

   double d=0.0;
   double dd=0.0;
   double ddcount = 0;
   int zcount =  0;
   double val;
   for (int i = 1;i<=w->rows;i++)
   {
      val = (double) *MATRIX_RELT(w, i, 1);
      d +=val ;
      if (fabs(*MATRIX_RELT(B, i, 1)) > 0.00001)
      {
         dd+= val ;
	 ddcount++;
         if (val < 0.1) zcount++;
      }
   }
   d /= w->rows;
   dd /= ddcount;
   cout << "          weights average: " << dd << "  zero: " << (double)zcount/ddcount << flush;
   //"  on significant b vals ( " << ddcount << " ): " << dd <<endl;
   lastweight = dd;
   lastzero   = (double)zcount/ddcount;
  return pair < MATRIX*, MATRIX*> (p,w);

}

MATRIX* Regression::getLSEst(MATRIX* outX)
{
  //cout << " Regression::getLSEst " << endl;
  lastweight = -1;
  lastzero   = -1;
  MATRIX* Ai = MatrixSVDPseudoInverse(A,NULL);
  if (Ai == NULL)
  {
     cerr << "    Regression::getLSEst   could not compute pseudo inverse!" << endl;
     //cerr << "    A: " << endl ;
     //MatrixPrintFmt(stdout,"% 2.8f",A);
     //cerr << endl;
     assert(Ai != NULL);
  }
  outX   = MatrixMultiply(Ai,B,outX);
  MatrixFree(&Ai);
  
  // compute error:
  MATRIX* R = MatrixMultiply(A,outX,NULL);
  R = MatrixSubtract(R,B,R);
  double serror = 0; 
  int cc, rr;
    for (cc = 1;cc<=R->cols;cc++)
    for (rr = 1;rr<=R->rows;rr++)
    {
      serror += R->rptr[rr][cc] * R->rptr[rr][cc];
    }
  MatrixFree(&R);
//  cout << "     squared error: " << serror << endl;
  lasterror = serror;
  
  return outX;
}

MATRIX* getTukeyBiweight(MATRIX* r, double sat, MATRIX* w)
{
	int n = r->rows;
	if (w == NULL) w = MatrixAlloc(n,1,  MATRIX_REAL); 
	else (assert (n == w->rows));
	assert(r->cols ==1 );

	double a , b;
	for (int i = 0;i<n;i++)
	{
		if (fabs(r->rptr[i][1]) > sat) w->rptr[i][1] = sat*sat/2.0;
		else
		{
	
			a = r->rptr[i][1]/sat;
			b = 1.0 - a * a;
			w->rptr[i][1] = (sat*sat/2.0)*(1.0-b*b*b);
		}
	}
	return w;
}


MATRIX * Regression::getSqrtTukeyDiaWeights(MATRIX * r, double sat, MATRIX * w)
// computes sqrt(weights) for a reweighted least squares approach 
// returns elements of diag matrix W (as a column vector)
{
	//cout << " getTukeyDiaWeights  r size: " << r->rows << " , " << r->cols << endl;

	assert(r->cols == 1);	
	int n = r->rows;
	if (w == NULL) w = MatrixAlloc(n,1,  MATRIX_REAL); 
	else (assert (r->rows == w->rows));
	assert(r->cols ==1 );
	
	double t1, t2;
	int rr,cc;
	int count = 1;
	int ocount = 0;
	for (rr = 1;rr<=r->rows;rr++)
	for (cc = 1;cc<=r->cols;cc++)
	{
	   // cout << " fabs: " << fabs(r->rptr[rr][cc]) << " sat: " << sat << endl;
		if (fabs(r->rptr[rr][cc]) > sat)
		{
		   w->rptr[count][1] = 0;
		   ocount++;
		}
		else
		{
			t1 = r->rptr[rr][cc]/sat;
			t2 = 1.0 - t1 *t1;
			//w->rptr[count][1] = t2 * t2;
			w->rptr[count][1] = t2 ;
		}
		count++;
	}
	//cout << " over threshold: " << ocount << " times ! " << endl;
	return w;
}


/*
  * This Quickselect routine is based on the algorithm described in
  * "Numerical recipes in C", Second Edition,
  * Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
  */
#define ELEM_SWAP(a,b) { double t=a;a=b;b=t; }
double Regression::quick_select(double arr[], int n, int k)
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
double Regression::kth_smallest(double a[], int n, int k)
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


double Regression::VectorMedian(MATRIX *v)
{
  int n = v->rows * v->cols;
  double* t = (double *)calloc(n, sizeof(double)); 
  
  int r,c,cc;
  // copy array to t
  cc = 0;
  for (r=1; r <= v->rows; r++)
    for (c=1; c <= v->cols; c++)
    {
      t[cc] = v->rptr[r][c];
      cc++;
    }
    //for (int i = 0;i<n;i++) cout << " " << t[i];  cout << endl;
    
  double qs;
  if (n%2 == 1) //odd
  {
    qs = kth_smallest(t,n,(n+1)/2);
    //cout << " n: " << n << "   " << qs << endl;
    free(t);
    return qs;
   }

  //  else even:

//  qs = kth_smallest(t,n,n/2);
  qs = quick_select(t,n,n/2);
 // double qs2 = kth_smallest(t,n,n/2 + 1);
  double qs2 = quick_select(t,n,n/2 + 1);
  //cout << " n: " << n << "   " << qs << "   " << qs2 << endl;
  qs =  0.5 * ( qs + qs2);
  free(t);
  return qs;
 
}

double Regression::getSigmaMAD(MATRIX *r, double d)
// robust estimate for sigma (using median absolute deviation)
// 1.4826 med_i(|r_i - med_j(r_j)|)
{
  double medi = VectorMedian(r);
  MATRIX* s = MatrixAlloc(r->rows, r->cols, MATRIX_REAL);
  int rr,cc;
  for (rr=1; rr <= s->rows; rr++)
  for (cc=1; cc <= s->cols; cc++)
  {
    s->rptr[rr][cc] = fabs(r->rptr[rr][cc] - medi);
  }
  
  //MatrixPrint(stdout, s);
  double mm = VectorMedian(s);
  MatrixFree(&s);
  return d * mm;

}
