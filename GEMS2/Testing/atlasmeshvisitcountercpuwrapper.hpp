#pragma once

#include "kvlAtlasMeshVisitCounterCPU.h"
#include "atlasmeshvisitcounter.hpp"

namespace kvl {
  class AtlasMeshVisitCounterCPUWrapper : public interfaces::AtlasMeshVisitCounter {
  public:
    AtlasMeshVisitCounterCPU::Pointer impl;

    AtlasMeshVisitCounterCPUWrapper() : impl(AtlasMeshVisitCounterCPU::New()) {}

    virtual void SetRegions( const ImageType::RegionType&  region ) override;
    virtual void VisitCount( const AtlasMesh* mesh ) override;
    virtual const AtlasMeshVisitCounterCPUWrapper::ImageType*  GetImage() const override;
  };
}
