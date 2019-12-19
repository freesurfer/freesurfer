#pragma once

#ifdef CUDA_FOUND
class CUDAGlobalFixture {
public:
  CUDAGlobalFixture();
};
#endif
