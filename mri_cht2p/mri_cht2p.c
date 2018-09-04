/**
 * @file  mri_cht2p.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.5 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "volcluster.h"
#include "numerics.h"

const char *Progname = "mri_cht2p";
int n;
CHT *cht, *cht2;


/*----------------------------------------*/
int main(int argc, char **argv) {
  double p;

  p = sc_ran_binomial_pdf (5, .7, 20) ;
  printf("p = %lf\n",p);
  exit(1);

  cht = CHTalloc(5,2,7,5,50,60);
  cht->nsim = 100000;
  cht->nvox = 100;
  cht->nsmooth = 200;
  cht->fwhm = 5.7;
  cht->totsize = 456.7;
  n = CHTsetSignString(cht, "pos");
  if (n) return(1);

  CHTprint(stdout, cht);
  CHTwrite("tmp.cht", cht);

  cht2 = CHTread("tmp.cht");
  CHTwrite("tmp2.cht", cht2);

  if (CHTcompare(cht,cht2)) printf("Different\n");
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
