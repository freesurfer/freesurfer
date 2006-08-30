#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "volcluster.h"

#if USE_SC_GSL_REPLACEMENT
  #include <gsl_wrapper.h>
#else
  #include "gsl/gsl_randist.h"
#endif

char *Progname = "mri_cht2p";
int n;
CHT *cht, *cht2;


/*----------------------------------------*/
int main(int argc, char **argv)
{
  double p;

#if USE_SC_GSL_REPLACEMENT
  p = sc_ran_binomial_pdf (5, .7, 20) ;
#else
  p = gsl_ran_binomial_pdf (5, .7, 20) ;
#endif
  printf("p = %lf\n",p);
  exit(1);

  cht = CHTalloc(5,2,7,5,50,60);
  cht->nsim = 100000;
  cht->nvox = 100;
  cht->nsmooth = 200;
  cht->fwhm = 5.7;
  cht->totsize = 456.7;
  n = CHTsetSignString(cht, "pos");
  if(n) return(1);

  CHTprint(stdout, cht);
  CHTwrite("tmp.cht", cht);

  cht2 = CHTread("tmp.cht");
  CHTwrite("tmp2.cht", cht2);

  if(CHTcompare(cht,cht2)) printf("Different\n");
  else                     printf("Same\n");

  CHTfree(&cht);
  CHTfree(&cht2);

  return(0);
}
/*----------------------------------------*/
/*----------------------------------------*/

/*----------------------------------------*/


/*----------------------------------------*/
//int binomialconf(int n, double theta, double confint)
//{
//}
