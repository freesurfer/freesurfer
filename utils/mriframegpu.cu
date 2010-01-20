/**
 * @file  mriframegpu.cu
 * @brief Holds MRI frame template for the GPU
 *
 * Holds an MRI frame template type for the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/01/20 18:52:53 $
 *    $Revision: 1.5 $
 *
 * Copyright (C) 2002-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */


#include "mriframegpu.hpp"

// ====================================================

template<> int MRIframeGPU<unsigned char>::MRItype( void ) const {
  return( MRI_UCHAR );
}

template<> int MRIframeGPU<short>::MRItype( void ) const {
  return( MRI_SHORT );
}

template<> int MRIframeGPU<float>::MRItype( void ) const {
  return( MRI_FLOAT );
}


// ====================================================

template<>
void MRIframeGPU<unsigned char>::ExhumeRow( const MRI* src,
					    unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame ) const {
  // Do the copy
  memcpy( h_slab,
	  &MRIseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(unsigned char) );
}


template<>
void MRIframeGPU<short>::ExhumeRow( const MRI* src,
				    short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) const {
  // Do the copy
  memcpy( h_slab,
	  &MRISseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(short) );
}



template<>
void MRIframeGPU<float>::ExhumeRow( const MRI* src,
				    float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) const {
  
  // Do the copy
  memcpy( h_slab,
	  &MRIFseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(float) );
}




// ====================================================


template<>
void MRIframeGPU<unsigned char>::InhumeRow( MRI* dst,
					    const unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame ) const {
  // Do the copy
  memcpy( &MRIseq_vox( dst, 0, iy, iz, iFrame ),
	  h_slab,
	  dst->width*sizeof(unsigned char) );
}


template<>
void MRIframeGPU<short>::InhumeRow( MRI* dst,
				    const short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) const {
  // Do the copy
  memcpy( &MRISseq_vox( dst, 0, iy, iz, iFrame ),
	  h_slab,
	  dst->width*sizeof(short) );
}



template<>
void MRIframeGPU<float>::InhumeRow( MRI* dst,
				    const float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) const {
  // Do the copy
  memcpy( &MRIFseq_vox( dst, 0, iy, iz, iFrame ),
	  h_slab,
	  dst->width*sizeof(float) );
}




// ======================================================


static void NullFunction( void ) {
  MRIframeGPU<unsigned char> tstFrame;

  tstFrame.cpuDims.x = 0;
}
