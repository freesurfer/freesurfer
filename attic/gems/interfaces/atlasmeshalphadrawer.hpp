#pragma once

#include "itkImage.h"
#include "kvlAtlasMesh.h"

namespace kvl {
  namespace interfaces {
    class AtlasMeshAlphaDrawer {
    public: 
      typedef itk::Image<float,3> ImageType;
      virtual void SetRegions( const ImageType::RegionType&  region ) = 0;
      virtual void Interpolate( const kvl::AtlasMesh* mesh ) = 0;
      virtual const ImageType* GetImage() const = 0;
      virtual void SetClassNumber( const int classNumber ) = 0;
    };
  }
}
