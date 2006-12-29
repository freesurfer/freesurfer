/**
 * @file  mris_extract_main_component.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:07 $
 *    $Revision: 1.2 $
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


/**
 * @file   dummy.c
 * @author Yasunari Tosa
 * @date   Wed Oct 13 11:47:18 2004
 *
 * @brief  sample dummy program
 *
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
    fprintf(stderr,"\n\nUSAGE: mris_extract_main_component input_surface output_surface\n\n");
    exit(-1);
  }

  mris_in=MRISread(argv[1]);

  mris_out=MRISextractMainComponent(mris_in,0);

  MRISwrite(mris_out,argv[2]);

  MRISfree(&mris_out);
  MRISfree(&mris_in);
  fprintf(stderr,"\ndone\n\n");
  return 0;
}


