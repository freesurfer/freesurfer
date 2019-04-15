/**
 * @file  cpputils.h
 * @brief include files for cpp utils
 *
 */
/*
 * Original Author: Krish Subramaniam
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.3 $
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

// include guard
#ifndef cpputils_h
#define cpputils_h

MRI* MRISsignedFixedDistanceTransform(MRI_SURFACE *mris, MRI *mri_dist, double distance);
MRI* MRISfillInterior2(MRI_SURFACE *mris, MRI* mri_interior);

#endif 
