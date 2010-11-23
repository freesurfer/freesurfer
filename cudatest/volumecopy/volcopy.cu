#include <iostream>

#include "mri.h"

#include "mriframegpu.hpp"
#include "cudacheck.h"


#include "volcopy.hpp"





template<typename T>
void VolumeCopy( const MRI* src, MRI *dst ) {

  GPU::Classes::MRIframeGPU<T> srcGPU, dstGPU;


  srcGPU.Allocate( src );
  srcGPU.Send( src, 0 );

  dstGPU.Copy( srcGPU );

  dstGPU.Recv( dst, 0 );
}







void VolCopyTest( const MRI* src, MRI *dst ) {

  switch( src->type ) {
  case MRI_UCHAR:
    VolumeCopy<unsigned char>( src, dst );
    break;

  case MRI_FLOAT:
    VolumeCopy<float>( src, dst );
    break;

  case MRI_SHORT:
    VolumeCopy<short>( src, dst );
    break;
    
  default:
    std::cerr << __FUNCTION__
	      << "Unrecognised MRI type " << src->type
	      << std::endl;
    exit( EXIT_FAILURE );
  }
}
