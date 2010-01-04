/**
 * @file  test_c_nr_wrapper.c
 * @brief test some of the numerical recipes replacements
 *
 */
/*
 * Original Author: Nick Schmansky
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2010/01/04 15:40:57 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
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
