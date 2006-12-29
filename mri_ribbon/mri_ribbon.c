/**
 * @file  mri_ribbon.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:08 $
 *    $Revision: 1.7 $
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "version.h"

#define IMGSIZE 256

char *Progname;

int main(int argc, char *argv[]) {
  char *inner_mris_fname,*outer_mris_fname,*input_mri_pref,*output_mri_pref;
  MRI *mri=0,*mri_src=0;
  MRI_SURFACE *inner_mris=0,*outer_mris=0;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_ribbon.c,v 1.7 2006/12/29 02:09:08 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  /* Set command-line parameters */
  if (argc!=5) {
    printf("Usage: mri_ribbon inner_surface_fname outer_surface_fname input_volume_pref output_volume_pref\n");
    exit(1);
  }
  Progname=argv[0];
  inner_mris_fname=argv[1];
  outer_mris_fname=argv[2];
  input_mri_pref=argv[3];
  output_mri_pref=argv[4];

  /* Read surface information from inner surface file */
  printf("Reading surface file %s.\n",inner_mris_fname);
  inner_mris=MRISread(inner_mris_fname);
  if (!inner_mris) {
    printf("Could not read surface file %s.\n",inner_mris_fname);
    exit(1);
  }

  /* Read surface information from outer surface file */
  printf("Reading surface file %s.\n",outer_mris_fname);
  outer_mris=MRISread(outer_mris_fname);
  if (!outer_mris) {
    printf("Could not read surface file %s.\n",outer_mris_fname);
    exit(1);
  }

  /* Read example volume from file */
  printf("Reading MRI volume directory %s.\n",input_mri_pref);
  mri_src=MRIread(input_mri_pref);
  if (!mri_src) {
    printf("Could not read MRI volume directory %s.\n",input_mri_pref);
    exit(1);
  }

  /* Extract ribbon */
  printf("Extracting ribbon.\n");
  mri=MRISribbon(inner_mris,outer_mris,mri_src,NULL);

  /* Save MRI volume to directory */
  printf("Writing volume file %s.\n",output_mri_pref);
  MRIwrite(mri,output_mri_pref);

  MRIfree(&mri);
  MRIfree(&mri_src);
  MRISfree(&inner_mris);
  MRISfree(&outer_mris);

  return 0;
}

