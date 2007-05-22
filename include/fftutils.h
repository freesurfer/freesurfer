/**
 * @file  fftutils.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nommert $
 *    $Date: 2007/05/22 14:12:50 $
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



#ifndef FFTUTILS_INC
#define FFTUTILS_INC


void CFFTforward(float* re, float* im, int length);
void CFFTbackward(float* re, float* im, int length);

void RFFTforward(float* data,int length, float* re, float* im );


void RFFT( float* data, int data_length, int length, int direction );
void FFTdebugAssert(int b, char *string);
int FFTisPowerOf2( int x ); 
int FFTpow2( int exponent );
float FFTdist(int x,int y,int z,float xsize,float ysize,float zsize,int len);
int FFTlog2( int x );
void FFTswitch_with_z (float *** vect, int dimension, int is_y);
void FFTmodarg_to_reim(float *** re_mod, float *** im_arg, int l);
float ***FFTinv_quarter(float *** vect, int dimension);
void FFTreim_to_modarg (float *** re_mod, float *** im_arg, int l);

#endif
