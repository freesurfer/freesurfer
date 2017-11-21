/**
 * @file  mri_exclude_minc_tests.c
 * @brief compares some of the exclude_minc algorithms with the minc algorithms
 *
 */
/*
 * Original Author: Bevin Brett
 *
 * Copyright © 2017 The General Hospital Corporation (Boston, MA) "MGH"
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


#include <stdio.h>
#include "minc_volume_io.h"


int main(int argc, char *argv[]) {

    printf("mri_exclude_minc_tests\n");
#if defined(BEVIN_EXCLUDE_MINC)
    volatile void* deleter = delete_general_transform;
    printf("BEVIN_EXCLUDE_MINC is defined, hence can't compare test_convert_transform_to_starts_and_steps with minc, but deleter=%p\n",deleter);
#else
    General_transform general_transform;

    VIO_Status status = input_transform_file("../talairach_afd/good_tal.xfm",  &general_transform );
    if (status != OK) { printf("input_transform_file failed\n"); return 1; }
 
    
    double step_signs                [VIO_N_DIMENSIONS] = { 1.0, 1.5, 5.1} ;
    int    spatial_axes              [VIO_N_DIMENSIONS] = { 3,   2,   1  } ;
    double starts		     [VIO_N_DIMENSIONS] = { 0.4, 0.5, 0.7} ;
    double steps                     [VIO_N_DIMENSIONS];
    double dir_cosines[3 /* X Y Z */][VIO_N_DIMENSIONS];
   
    int trial;
    for (trial = 0; trial < 2; trial++) {
        printf(trial?"test\n":"good\n");
        (trial ? test_convert_transform_to_starts_and_steps 
	       :      convert_transform_to_starts_and_steps) (
    	&general_transform,
    	VIO_N_DIMENSIONS,	// inp
    	step_signs,		// inp
    	spatial_axes,		// inp
    	starts,			// inp
    	steps,			// out
    	dir_cosines);		// out
    
        printf("steps = ");
	int i,j;
	for (i = 0; i < VIO_N_DIMENSIONS; i++) printf("%d:%lf ", i, steps[i]); 
        printf("\n dir_cosines[X] = ");
	for (i = 0; i < 3; i++) {
	    for (j = 0; j < 3; j++) {
	        printf("%d,%d:%lf ", i,j, dir_cosines[i][j]);
	    }
	    printf("\n");
	} 
    }
#endif  
    return 0;
}
