/**
 * @brief I/O and analysis algorithms for FCDs (focal cortical dysplasias)
 *
 */
/*
 * Original Author: Bruce Fischl
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


#ifndef FCD_H
#define FCD_H

#include "mri.h"
#include "mrisurf.h"
#include "label.h"
#include "utils.h"
#include "label.h"
#define MAX_FCD_LABELS 100

typedef struct
{
  MRI_SURFACE *mris_lh ;
  MRI_SURFACE *mris_rh ;
  MRI_SURFACE *mris_lh_pial;
  MRI_SURFACE *mris_rh_pial;
  MRI_SURFACE *mris_lh_sphere_d1;
  MRI_SURFACE *mris_rh_sphere_d1;
  MRI         *mri_aseg ;
  MRI         *mri_aparc ;
  MRI         *mri_norm ;
  MRI         *mri_flair ;
  MRI         *mri_t2 ;
  MRI         *mri_thickness_increase ;
  MRI         *mri_thickness_decrease ;
  MRI         *mri_thickness_difference;
  MRI         *lh_thickness_on_lh ;  // thickness of lh mapped to lh
  MRI         *lh_thickness_on_rh ;  // thickness of lh mapped to rh
  MRI         *rh_thickness_on_lh ;  // thickness of rh mapped to lh
  MRI         *rh_thickness_on_rh ;  // thickness of rh mapped to rh
  double      thickness_threshold ;
  int         thickness_smooth_steps ;
  int         nlabels ;
  LABEL       *labels[MAX_FCD_LABELS] ;
  char        label_names[MAX_FCD_LABELS][STRLEN] ;
} FCD_DATA ;


FCD_DATA   *FCDloadData(char *sdir, char *subject, char* suffix);
int         FCDcomputeThicknessLabels(FCD_DATA *fcd,
                                      double thickness_thresh,
                                      double sigma,
                                      int size_thresh) ;
int         FCDwriteLabels(FCD_DATA *fcd, char* dir) ;
int         FCDfree(FCD_DATA **pfcd) ;

#endif
