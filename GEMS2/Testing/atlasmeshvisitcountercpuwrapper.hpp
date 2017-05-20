#pragma once

#include "kvlAtlasMeshVisitCounter.h"
#include "atlasmeshvisitcounter.hpp"

#include "stopwatch.hpp"

namespace kvl {
  class AtlasMeshVisitCounterCPUWrapper : public interfaces::AtlasMeshVisitCounter {
  public:
    kvl::AtlasMeshVisitCounter::Pointer impl;

    AtlasMeshVisitCounterCPUWrapper() : impl(kvl::AtlasMeshVisitCounter::New()) {}

    virtual void SetRegions( const ImageType::RegionType&  region ) override;
    virtual void VisitCount( const AtlasMesh* mesh ) override;
    virtual const AtlasMeshVisitCounterCPUWrapper::ImageType*  GetImage() const override;
    
    Stopwatch tSetRegions, tVisitCount;
  };
}
