/**
 * @file  cafactory.hpp
 * @brief Holds 'factory' object for textures
 *
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef CUDA_TEXTURE_FACTORY_HPP
#define CUDA_TEXTURE_FACTORY_HPP

#include "cudacheck.h"
#include "cudatypeutils.hpp"

#include "volumegpu.hpp"

namespace GPU
{

namespace Classes
{

//! Class to encapsulate a cudaArray and texture binding
/*!
  This is a helper class, which copies a given VolumeGPU
  into a cudaArray, and binds that array to the designated
  texture.
  The cudaArray is held internally to the class, and deleted
  by the destructor.
*/
class CTfactory
{
public:

  //! Construct from VolumeGPU
  template<typename T,typename U>
  CTfactory( const VolumeGPU<T>& src,
             U& texRef,
             const cudaTextureFilterMode fm = cudaFilterModePoint,
             const cudaTextureAddressMode am = cudaAddressModeClamp,
             const int norm = false ) : dca_data(NULL)
  {

    // Check for valid input
    if( src.d_data.ptr == NULL )
    {
      std::cerr << __FUNCTION__
                << ": Source has no data"
                << std::endl;
      abort();
    }

    // Allocate memory
    cudaChannelFormatDesc cd = cudaCreateChannelDesc<T>();
    cudaExtent tmpExtent = ExtentFromDims( src.dims );

    CUDA_SAFE_CALL( cudaMalloc3DArray( &(this->dca_data),
                                       &cd,
                                       tmpExtent ) );

    // Do the copy
    cudaMemcpy3DParms cp = {0};

    cp.srcPtr = src.d_data;
    cp.dstArray = this->dca_data;
    cp.extent = tmpExtent;
    cp.kind = cudaMemcpyDeviceToDevice;

    CUDA_SAFE_CALL( cudaMemcpy3D( &cp ) );


    // Bind the texture
    texRef.normalized = norm;
    texRef.addressMode[0] = am;
    texRef.addressMode[1] = am;
    texRef.addressMode[2] = am;
    texRef.filterMode = fm;

    CUDA_SAFE_CALL( cudaBindTextureToArray( texRef, this->dca_data ) );
  }

  //! Destructor releases array
  ~CTfactory( void )
  {
    if( this->dca_data != NULL )
    {
      // Actually, shouldn't exist without an array...
      cudaFreeArray( this->dca_data );
    }
  }

private:
  //! The actual CUDA array
  cudaArray *dca_data;

  //! Inhibit default construction
  CTfactory( void ) : dca_data(NULL)
  {
    std::cerr << __FUNCTION__
              << ": Default construction inhibited"
              << std::endl;
    abort();
  }

  //! Inhibit copy construction
  CTfactory( const CTfactory& src ) : dca_data(NULL)
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy"
              << std::endl;
    abort();
  }

  //! Inhibit assigment
  void operator=( const CTfactory& src )
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy"
              << std::endl;
    abort();
  }

};

}
}


#endif
