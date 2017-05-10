#pragma once

#include "atlasmeshvisitcounter.hpp"
#include "stopwatch.hpp"

namespace kvl {
  namespace cuda {
    class AtlasMeshVisitCounterCUDA : public kvl::interfaces::AtlasMeshVisitCounter {
    public:
      AtlasMeshVisitCounterCUDA() {}

      virtual void SetRegions( const kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType& region ) override;
      virtual void VisitCount( const kvl::AtlasMesh* mesh ) override;
      virtual const AtlasMeshVisitCounterCUDA::ImageType*  GetImage() const override;
    };
  }
}
