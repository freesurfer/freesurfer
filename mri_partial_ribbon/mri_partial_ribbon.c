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
  char *inner_mris_fname_lh,*outer_mris_fname_lh,*inner_mris_fname_rh,*outer_mris_fname_rh,*input_mri_pref,*output_mri_pref,*mask_mri_pref;
  MRI *mri,*mri_src,*mri_mask;
  MRI_SURFACE *inner_mris_lh,*outer_mris_lh,*inner_mris_rh,*outer_mris_rh;

  /* Set command-line parameters */
  if ((argc!=7)&&(argc!=8)) {
    printf("Usage: mri_partial_ribbon inner_surface_fname_lh outer_surface_fname_lh inner_surface_fname_rh outer_surface_fname_rh input_volume_pref output_volume_pref [mask_fname]\n");
    exit(1);
  }
  Progname=argv[0];

  inner_mris_fname_lh=argv[1];
  outer_mris_fname_lh=argv[2];
  inner_mris_fname_rh=argv[3];
  outer_mris_fname_rh=argv[4];
  input_mri_pref=argv[5];
  output_mri_pref=argv[6];
  if (argc==8)
    mask_mri_pref=argv[7];
  else
    mask_mri_pref=NULL;

  /* Read surface information from left inner surface file */
  printf("Reading left inner surface file %s.\n",inner_mris_fname_lh);
  inner_mris_lh=MRISread(inner_mris_fname_lh);
  if (!inner_mris_lh) {
    fprintf(stderr,"Could not read surface file %s.\n",inner_mris_fname_lh);
    exit(1);
  }

  /* Read surface information from left outer surface file */
  printf("Reading left outer surface file %s.\n",outer_mris_fname_lh);
  outer_mris_lh=MRISread(outer_mris_fname_lh);
  if (!outer_mris_lh) {
    fprintf(stderr,"Could not read surface file %s.\n",outer_mris_fname_lh);
    exit(1);
  }

  /* Read surface information from right inner surface file */
  printf("Reading right inner surface file %s.\n",inner_mris_fname_rh);
  inner_mris_rh=MRISread(inner_mris_fname_rh);
  if (!inner_mris_rh) {
    fprintf(stderr,"Could not read surface file %s.\n",inner_mris_fname_rh);
    exit(1);
  }

  /* Read surface information from right outer surface file */
  printf("Reading right outer surface file %s.\n",outer_mris_fname_rh);
  outer_mris_rh=MRISread(outer_mris_fname_rh);
  if (!outer_mris_rh) {
    fprintf(stderr,"Could not read surface file %s.\n",outer_mris_fname_rh);
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
  if (argc==8) {
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
  mri=MRISpartialribbon(inner_mris_lh,outer_mris_lh,inner_mris_rh,outer_mris_rh,mri_src,NULL,mri_mask);

  /* Save MRI volume to directory */
  printf("Writing volume file %s.\n",output_mri_pref);
  MRIwrite(mri,output_mri_pref);

  MRIfree(&mri);
  MRIfree(&mri_src);
  MRIfree(&mri_mask);
  MRISfree(&inner_mris_lh);
  MRISfree(&outer_mris_lh);
  MRISfree(&inner_mris_rh);
  MRISfree(&outer_mris_rh);

  return 0;
}

