#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"

#define IMGSIZE 256

char *Progname;

int main(int argc, char *argv[])
{
  char *inner_mris_fname,*outer_mris_fname,*input_mri_pref,*output_mri_pref,*mask_mri_pref;
  MRI *mri,*mri_src,*mri_mask;
  MRI_SURFACE *inner_mris,*outer_mris;
  int hemisphere;

  /* Set command-line parameters */
  if ((argc!=5)&&(argc!=7)) {
    printf("Usage: mri_partial_ribbon inner_surface_fname outer_surface_fname input_volume_pref output_volume_pref [mask_fname left/right]\n");
    exit(1);
  }
  Progname=argv[0];

  inner_mris_fname=argv[1];
  outer_mris_fname=argv[2];
  input_mri_pref=argv[3];
  output_mri_pref=argv[4];
  if (argc==7) {
    mask_mri_pref=argv[5];
    if (strcasecmp(argv[6],"left")==0)
      hemisphere=MRI_LEFT_HEMISPHERE;
    else
      if (strcasecmp(argv[6],"right")==0)
        hemisphere=MRI_RIGHT_HEMISPHERE;
      else {
        fprintf(stderr,"Hemisphere must be specified as left or right.\n");
        exit(1);
      }
  }
  else {
    mask_mri_pref=NULL;
    hemisphere=0;
  }

  /* Read surface information from inner surface file */
  printf("Reading surface file %s.\n",inner_mris_fname);
  inner_mris=MRISread(inner_mris_fname);
  if (!inner_mris) {
    fprintf(stderr,"Could not read surface file %s.\n",inner_mris_fname);
    exit(1);
  }

  /* Read surface information from outer surface file */
  printf("Reading surface file %s.\n",outer_mris_fname);
  outer_mris=MRISread(outer_mris_fname);
  if (!outer_mris) {
    fprintf(stderr,"Could not read surface file %s.\n",outer_mris_fname);
    exit(1);
  }

  /* Read example volume from file */
  printf("Reading MRI volume directory %s.\n",input_mri_pref);
  mri_src=MRIread(input_mri_pref);
  if (!mri_src) {
    fprintf(stderr,"Could not read MRI volume directory %s.\n",input_mri_pref);
    exit(1);
  }

  /* Read example volume from file */
  if (argc==7) {
    printf("Reading MRI volume directory %s.\n",mask_mri_pref);
    mri_mask=MRIread(mask_mri_pref);
    if (!mri_mask) {
      fprintf(stderr,"Could not read MRI volume directory %s.\n",mask_mri_pref);
      exit(1);
    }
  }
  else
    mri_mask=NULL;

  /* Extract ribbon */
  printf("Extracting ribbon.\n");
  mri=MRISpartialribbon(inner_mris,outer_mris,mri_src,NULL,mri_mask,hemisphere);

  /* Save MRI volume to directory */
  printf("Writing volume file %s.\n",output_mri_pref);
  MRIwrite(mri,output_mri_pref);

  MRIfree(&mri);
  MRIfree(&mri_src);
  MRIfree(&mri_mask);
  MRISfree(&inner_mris);
  MRISfree(&outer_mris);

  return 0;
}

