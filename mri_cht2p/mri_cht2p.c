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
 *    $Date: 2006/12/29 02:09:06 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "volcluster.h"
#include "numerics.h"

char *Progname = "mri_cht2p";
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
