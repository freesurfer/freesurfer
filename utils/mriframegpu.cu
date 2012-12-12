/**
 * @file  mriframegpu.cu
 * @brief Holds MRI frame template for the GPU
 *
 * Holds template specialisations for MRI frame on the GPU
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.9 $
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


#include "mriframegpu.hpp"

namespace GPU
{

namespace Classes
{

// ====================================================

template<> int MRIframeGPU<unsigned char>::MRItype( void ) const
{
  return( MRI_UCHAR );
}

template<> int MRIframeGPU<short>::MRItype( void ) const
{
  return( MRI_SHORT );
}

template<> int MRIframeGPU<float>::MRItype( void ) const
{
  return( MRI_FLOAT );
}

template<> int MRIframeGPU<int>::MRItype( void ) const
{
  return( MRI_INT );
}


// ====================================================

template<>
void MRIframeGPU<unsigned char>::ExhumeRow( const MRI* src,
    unsigned char* h_slab,
    const unsigned int iy,
    const unsigned int iz,
    const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( h_slab,
          &MRIseq_vox( src, 0, iy, iz, iFrame ),
          src->width*sizeof(unsigned char) );
}

// -----

template<>
void MRIframeGPU<short>::ExhumeRow( const MRI* src,
                                    short* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( h_slab,
          &MRISseq_vox( src, 0, iy, iz, iFrame ),
          src->width*sizeof(short) );
}

// -----

template<>
void MRIframeGPU<float>::ExhumeRow( const MRI* src,
                                    float* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const
{

  // Do the copy
  memcpy( h_slab,
          &MRIFseq_vox( src, 0, iy, iz, iFrame ),
          src->width*sizeof(float) );
}

// -----

template<>
void MRIframeGPU<int>::ExhumeRow( const MRI* src,
                                  int* h_slab,
                                  const unsigned int iy,
                                  const unsigned int iz,
                                  const unsigned int iFrame ) const
{

  // Do the copy
  memcpy( h_slab,
          &MRIIseq_vox( src, 0, iy, iz, iFrame ),
          src->width*sizeof(int) );
}



// ====================================================


template<>
void MRIframeGPU<unsigned char>::InhumeRow( MRI* dst,
    const unsigned char* h_slab,
    const unsigned int iy,
    const unsigned int iz,
    const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( &MRIseq_vox( dst, 0, iy, iz, iFrame ),
          h_slab,
          dst->width*sizeof(unsigned char) );
}

// -----

template<>
void MRIframeGPU<short>::InhumeRow( MRI* dst,
                                    const short* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( &MRISseq_vox( dst, 0, iy, iz, iFrame ),
          h_slab,
          dst->width*sizeof(short) );
}

// -----

template<>
void MRIframeGPU<float>::InhumeRow( MRI* dst,
                                    const float* h_slab,
                                    const unsigned int iy,
                                    const unsigned int iz,
                                    const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( &MRIFseq_vox( dst, 0, iy, iz, iFrame ),
          h_slab,
          dst->width*sizeof(float) );
}

// -----

template<>
void MRIframeGPU<int>::InhumeRow( MRI* dst,
                                  const int* h_slab,
                                  const unsigned int iy,
                                  const unsigned int iz,
                                  const unsigned int iFrame ) const
{
  // Do the copy
  memcpy( &MRIIseq_vox( dst, 0, iy, iz, iFrame ),
          h_slab,
          dst->width*sizeof(int) );
}


}
}
