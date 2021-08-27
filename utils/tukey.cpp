/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <stdlib.h>
#include <math.h>

#include "tukey.h"


double tukey_biweight(double residual, double C)
{
  if (abs(residual) > C) {
    return (C * C / 2);
  } else {
    double p = residual / C;
    p *= p;
    p = 1 - p;
    p = p * p * p;
    return (((C * C) / 2) * (1 - p));
  }
}
