/**
 * @file  em_register_cuda.h
 * @brief Header file for em_register CUDA routines
 *
 * Contains CUDA routines for em_register
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:15 $
 *    $Revision: 1.2 $
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


#ifndef EM_REGISTER_CUDA_H
#define EM_REGISTER_CUDA_H

#if defined(__cplusplus)
extern "C" {
#endif
  
  //! Routine to set up the MRI and GCAsample on the GPU
  void CUDA_em_register_Prepare( GCA *gca,
				 GCA_SAMPLE *gcas,
				 const MRI *mri,
				 const int nSamples );

  //! Routine to release resources on the GPU
  void CUDA_em_register_Release( void );

  //! Routine to compute the probability of the given transform matrix
  float CUDA_ComputeLogSampleProbability( const MATRIX *m_L );

  //! Routine to find the best translation
  void CUDA_FindOptimalTranslation( const MATRIX *baseTransform,
				    const float minTrans,
				    const float maxTrans,
				    const unsigned int nTrans,
				    float *maxLogP,
				    float *dx,
				    float *dy,
				    float *dz );

  //! Routine to find the best transformation
  void CUDA_FindOptimalTransform( const MATRIX *baseTransform,
				  const MATRIX *originTranslation,
				  const float minTrans,
				  const float maxTrans,
				  const unsigned int nTrans,
				  const float minScale,
				  const float maxScale,
				  const unsigned nScale,
				  const float minRot,
				  const float maxRot,
				  const unsigned int nRot,
				  double *maxLogP,
				  double *dx,
				  double *dy,
				  double *dz,
				  double *sx,
				  double *sy,
				  double *sz,
				  double *rx,
				  double *ry,
				  double *rz );
#if defined(__cplusplus)
};
#endif


#endif
