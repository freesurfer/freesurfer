/**
 * @file  mris_extract_main_component.c
 * @brief extract the main connected component from a surface
 *
 */
/*
 * Original Author: Florent Segonne
 * CVS Revision Info:
 *    $Author: segonne $
 *    $Date: 2007/02/13 17:17:30 $
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "MRIio_old.h"
#include "mri.h"
#include "mrisurf.h"
#include "gca.h"
#include "MC.h"

char *Progname;

int main(int argc, char *argv[]) {
  MRIS *mris_in,*mris_out;

  Progname=argv[0];

  if (argc < 3) {
    fprintf
      (stderr,
       "\n\nUSAGE: mris_extract_main_component "
       "input_surface output_surface\n\n");
    exit(-1);
  }

  mris_in=MRISread(argv[1]);

  mris_out=MRISextractMainComponent(mris_in,0,1,0);

  MRISwrite(mris_out,argv[2]);

  MRISfree(&mris_out);
  MRISfree(&mris_in);
  fprintf(stderr,"\ndone\n\n");
  return 0;
}


