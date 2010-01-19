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
 *    $Date: 2010/01/19 20:03:03 $
 *    $Revision: 1.3 $
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


template<> int GetAsMRItype<unsigned char>( const unsigned char tmp ) {
  return( MRI_UCHAR );
}


template<> int GetAsMRItype<short>( const short tmp ) {
  return( MRI_SHORT );
}

template<> int GetAsMRItype<float>( const float tmp ) {
  return( MRI_FLOAT );
}



// ====================================================

template<>
void CopyMRIrowToContiguous<unsigned char>( const MRI* src, unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame ) {
  // Do the copy
  memcpy( h_slab,
	  &MRIseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(unsigned char) );
}


template<>
void CopyMRIrowToContiguous<short>( const MRI* src, short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) {
  // Do the copy
  memcpy( h_slab,
	  &MRISseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(short) );
}



template<>
void CopyMRIrowToContiguous<float>( const MRI* src, float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) {
  
  // Do the copy
  memcpy( h_slab,
	  &MRIFseq_vox( src, 0, iy, iz, iFrame ),
	  src->width*sizeof(float) );
}




// ====================================================


template<>
void CopyMRIcontiguousToRow<unsigned char>( MRI* dst, const unsigned char* h_slab,
					    const unsigned int iy,
					    const unsigned int iz,
					    const unsigned int iFrame ) {
  // Do the copy
  memcpy( &MRIseq_vox( dst, 0, iy, iz, iFrame ),
	  h_slab,
	  dst->width*sizeof(unsigned char) );
}


template<>
void CopyMRIcontiguousToRow<short>( MRI* dst, const short* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) {
  // Do the copy
  memcpy( &MRISseq_vox( dst, 0, iy, iz, iFrame ),
	  h_slab,
	  dst->width*sizeof(short) );
}



template<>
void CopyMRIcontiguousToRow<float>( MRI* dst, const float* h_slab,
				    const unsigned int iy,
				    const unsigned int iz,
				    const unsigned int iFrame ) {
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
