#include <cuda_runtime.h>
#include <iostream>

#include "cudacheck.hpp"
#include "cudacontroller.hpp"

#include "cudaglobalfixture.hpp"

CUDAGlobalFixture::CUDAGlobalFixture() {
  const int cudaDevice = 0;
  kvl::cuda::InitialiseCUDA(cudaDevice);
  
  cudaDeviceProp properties;
  CUDA_SAFE_CALL( cudaGetDeviceProperties(&properties, cudaDevice) );
  
  // It seems that we can't use BOOST_TEST_MESSAGE in the global fixture
  std::cout <<  "CUDA Device : " << properties.name << std::endl;
}
