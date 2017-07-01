#pragma once

#include "kvlAtlasMeshAlphaDrawer.h"
#include "atlasmeshalphadrawer.hpp"

#include "stopwatch.hpp"

namespace kvl {
  class AtlasMeshAlphaDrawerCPUWrapper : public interfaces::AtlasMeshAlphaDrawer {
  public:
    kvl::AtlasMeshAlphaDrawer::Pointer impl;

    AtlasMeshAlphaDrawerCPUWrapper() : impl(kvl::AtlasMeshAlphaDrawer::New()) {}

    virtual void SetRegions( const ImageType::RegionType&  region ) override;
    virtual void Interpolate( const kvl::AtlasMesh* mesh ) override;
    virtual const AtlasMeshAlphaDrawerCPUWrapper::ImageType* GetImage() const override;
    virtual void SetClassNumber( const int classNumber ) override;

    Stopwatch tSetRegions, tInterpolate;
  };
}
