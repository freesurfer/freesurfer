#include <stdio.h>
#include <math.h>

double legendre_Pmm(int m, double x)
{
   if(m == 0)
   {
     return 1.0;
   }
   else
   {
     double p_mm = 1.0;
     double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
     double fact_coeff = 1.0;
     int i;
     for(i=1; i<=m; i++)
     {
       p_mm *= -fact_coeff * root_factor;
       fact_coeff += 2.0;
     }
     return p_mm;
   }
}

double gsl_sf_legendre_Plm_e(const int l, const int m, const double x)
{
   /* If l is large and m is large, then we have to worry
    * about overflow. Calculate an approximate exponent which
    * measures the normalization of this thing.
    */
   //const double dif = l-m;
   //const double sum = l+m;
   //const double t_d = ( dif == 0.0 ? 0.0 : 0.5 * dif * (log(dif)-1.0) );
   //const double t_s = ( dif == 0.0 ? 0.0 : 0.5 * sum * (log(sum)-1.0) );

   /* CHECK_POINTER(result) */

   if(m < 0 || l < m || x < -1.0 || x > 1.0) {
     printf("Domain error\n");
     return(0);
   }
   /* Account for the error due to the
    * representation of 1-x.
    */

   /* P_m^m(x) and P_{m+1}^m(x) */
   double p_mm   = legendre_Pmm(m, x);
   double p_mmp1 = x * (2*m + 1) * p_mm;

   if(l == m){
     return p_mm;
   }
   else if(l == m + 1) {
     return p_mmp1;
   }
   /* upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
    * start at P(m,m), P(m+1,m)
    */

   double p_ellm2 = p_mm;
   double p_ellm1 = p_mmp1;
   double p_ell = 0.0;
   int ell;

   for(ell=m+2; ell <= l; ell++){
     p_ell = (x*(2*ell-1)*p_ellm1 - (ell+m-1)*p_ellm2) / (ell-m);
     p_ellm2 = p_ellm1;
     p_ellm1 = p_ell;
   }

   return(p_ell);
}

double factorial(int n)
{
   int m;
   double f=1;
   for(m=2; m<=n; m++) f *= m;
   return(f);
}
