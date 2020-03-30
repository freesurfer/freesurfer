/*
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef JOINT_HISTO_H
#define JOINT_HISTO_H

#include <stdio.h>

#ifndef UCHAR
#define UCHAR  unsigned char
#endif

#define MAX_BINS  10000
typedef struct
{
  int     nbins_1 ;
  int     nbins_2 ;
  //int     total_bins ; 
  int     sample_count ; 
  float   *counts ; /* # of voxels which map to this bin */

  //float   bin_size ;
  float   min ;
  float   max ;     // min and max vals in the histo
}
JOINT_HISTOGRAM, JHISTO ;

int       JHISTOfree(JOINT_HISTOGRAM **pjhisto) ;
int       JHISTOdump(JOINT_HISTOGRAM *jhisto, FILE *fp) ;
int       JHISTOwriteInto(JOINT_HISTOGRAM *jhisto, FILE *fp) ;
void      JHISTOfill(MRI* mri1, MRI* mri2, JOINT_HISTOGRAM* jhisto);
double    JHISTOgetEntropy(JOINT_HISTOGRAM* jhisto);
JOINT_HISTOGRAM* JHISTOreadFrom(FILE *fp) ;
JOINT_HISTOGRAM* JHISTOalloc(int nbins_1, int nbins_2) ;
JOINT_HISTOGRAM* JHISTOrealloc(JOINT_HISTOGRAM *jhisto, int nbins_1, int nbins_2) ;
double MRIcomputeMi(MRI *mri1, MRI *mri2, int bins1, int bins2);

#endif
