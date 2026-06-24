#pragma once

namespace kvl {
  namespace cuda {
    void InitialiseCUDA(const int deviceID);
    void FinalizeCUDA();
  }
}
