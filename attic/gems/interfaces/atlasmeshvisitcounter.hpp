#pragma once

#include "itkImage.h"
#include "kvlAtlasMesh.h"

namespace kvl {
  namespace interfaces {
    class AtlasMeshVisitCounter {
    public:
      typedef itk::Image<int,3> ImageType;
      virtual void SetRegions( const ImageType::RegionType&  region ) = 0;
      virtual void VisitCount( const kvl::AtlasMesh* mesh ) = 0;
      virtual const ImageType* GetImage() const = 0;
    };
  }
}
