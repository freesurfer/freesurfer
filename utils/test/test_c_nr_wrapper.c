/**
 * @brief test some of the numerical recipes replacements
 *
 */
/*
 * Original Author: Nick Schmansky
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
#include <time.h>

#include "numerics.h"

const char* Progname = "test_c_nr_wrapper";

int main(int argc, char *argv[])
{
  // just trying to make sure that we can call the functions from c
  long seed = -1L * (long)( abs( (int)time(NULL) ) );

  printf("random number: %f\n", OpenRan1( &seed ));

  return 0;
}
