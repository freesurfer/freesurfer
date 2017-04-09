#include <vector>

#include <boost/test/unit_test.hpp>

// This will only be compiled if CUDA_FOUND is defined
#include "dimensioncuda.hpp"
#include "cudaimage.hpp"

#include "testiosupport.hpp"

BOOST_AUTO_TEST_SUITE( CudaImage )

BOOST_AUTO_TEST_SUITE( SendSetReceive )

BOOST_AUTO_TEST_CASE( SendSetReceive1d )
{
  kvl::cuda::Dimension<1,unsigned long> srcDims, resultDims;

  srcDims[0] = 6347;

  std::vector<int> src, result;
  src.resize(srcDims.ElementCount());
  for( auto it=src.begin(); it!=src.end(); it++ ) {
    *it = 1;
  }

  kvl::cuda::CudaImage<int,1,size_t> d_image;

  d_image.Send( src, srcDims );
  d_image.SetMemory(0);
  d_image.Recv( result, resultDims );

  BOOST_CHECK_EQUAL( srcDims, resultDims );
  for( auto it=result.begin(); it!=result.end(); it++ ) {
    BOOST_CHECK_EQUAL( *it, 0 );
  }
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
