#pragma once

#include "atlasmeshalphadrawer.hpp"
#include "cudatetrahedralmesh.hpp"

#include "stopwatch.hpp"

namespace kvl {
  namespace cuda {
    class AtlasMeshAlphaDrawerCUDA : public kvl::interfaces::AtlasMeshAlphaDrawer {
    public:
      AtlasMeshAlphaDrawerCUDA() {}

      virtual void SetRegions( const ImageType::RegionType&  region ) override;
      virtual void Interpolate( const kvl::AtlasMesh* mesh ) override;
      virtual const AtlasMeshAlphaDrawerCUDA::ImageType* GetImage() const override;
      virtual void SetClassNumber( const int classNumber ) override;

      mutable Stopwatch tSetRegions;

      mutable Stopwatch tGetImage, tGetImageTransfer, tGetImageUnpack;
      
      mutable Stopwatch tInterpolate, tSendMesh, tKernel;
    private:
      int classNumber;

      CudaImage<float,3,unsigned short> d_Output;
      AtlasMeshAlphaDrawerCUDA::ImageType::Pointer image;
    };
  }
}
