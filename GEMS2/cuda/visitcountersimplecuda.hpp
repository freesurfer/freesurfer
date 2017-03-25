#pragma once

#include "atlasmeshvisitcounter.hpp"

namespace kvl {
  namespace cuda {
    template<typename T>
    class VisitCounterSimple : public kvl::interfaces::AtlasMeshVisitCounter {
      virtual void SetRegions( const kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType& region ) override {}

      virtual void VisitCount( const kvl::AtlasMesh* mesh ) override {

      };

      virtual const VisitCounterSimple<T>::ImageType* GetImage() const override {
	throw std::runtime_error("Not implemented");
      }

    };
  }
}
