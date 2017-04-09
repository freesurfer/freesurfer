#include <vector>

#include <boost/test/unit_test.hpp>

// This will only be compiled if CUDA_FOUND is defined
#include "dimensioncuda.hpp"
#include "cudaimage.hpp"

#include "randomsupply.hpp"
#include "testiosupport.hpp"

RandomSupply gPRNG;

BOOST_AUTO_TEST_SUITE( CudaImage )

BOOST_AUTO_TEST_SUITE( SendSetReceive )

BOOST_AUTO_TEST_CASE( SendSetReceive1d )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<1,unsigned long> srcDims, resultDims;

    srcDims[0] = gPRNG.GetInteger<unsigned long>(1024*1024*4);
    
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
}

BOOST_AUTO_TEST_CASE( SendSetReceive2d )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<2,unsigned long> srcDims, resultDims;

    srcDims[0] = gPRNG.GetInteger<unsigned long>(4096);
    srcDims[1] = gPRNG.GetInteger<unsigned long>(4096);
    
    std::vector<int> src, result;
    src.resize(srcDims.ElementCount());
    for( auto it=src.begin(); it!=src.end(); it++ ) {
      *it = 1;
    }
    
    kvl::cuda::CudaImage<int,2,size_t> d_image;
    
    d_image.Send( src, srcDims );
    d_image.SetMemory(0);
    d_image.Recv( result, resultDims );
    
    BOOST_CHECK_EQUAL( srcDims, resultDims );
    for( auto it=result.begin(); it!=result.end(); it++ ) {
      BOOST_CHECK_EQUAL( *it, 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( SendSetReceive3d )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<3,unsigned long> srcDims, resultDims;

    srcDims[0] = gPRNG.GetInteger<unsigned long>(256);
    srcDims[1] = gPRNG.GetInteger<unsigned long>(256);
    srcDims[2] = gPRNG.GetInteger<unsigned long>(256);

    
    std::vector<int> src, result;
    src.resize(srcDims.ElementCount());
    for( auto it=src.begin(); it!=src.end(); it++ ) {
      *it = 1;
    }
    
    kvl::cuda::CudaImage<int,3,size_t> d_image;
    
    d_image.Send( src, srcDims );
    d_image.SetMemory(0);
    d_image.Recv( result, resultDims );
    
    BOOST_CHECK_EQUAL( srcDims, resultDims );
    for( auto it=result.begin(); it!=result.end(); it++ ) {
      BOOST_CHECK_EQUAL( *it, 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( SendSetReceive4d )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<4,unsigned long> srcDims, resultDims;

    srcDims[0] = gPRNG.GetInteger<unsigned long>(64);
    srcDims[1] = gPRNG.GetInteger<unsigned long>(64);
    srcDims[2] = gPRNG.GetInteger<unsigned long>(32);
    srcDims[3] = gPRNG.GetInteger<unsigned long>(32);
    
    std::vector<int> src, result;
    src.resize(srcDims.ElementCount());
    for( auto it=src.begin(); it!=src.end(); it++ ) {
      *it = 1;
    }
    
    kvl::cuda::CudaImage<int,4,size_t> d_image;
    
    d_image.Send( src, srcDims );
    d_image.SetMemory(0);
    d_image.Recv( result, resultDims );
    
    BOOST_CHECK_EQUAL( srcDims, resultDims );
    for( auto it=result.begin(); it!=result.end(); it++ ) {
      BOOST_CHECK_EQUAL( *it, 0 );
    }
  }
}

BOOST_AUTO_TEST_CASE( SendSetReceive5d )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<5,unsigned long> srcDims, resultDims;

    srcDims[0] = gPRNG.GetInteger<unsigned long>(64);
    srcDims[1] = gPRNG.GetInteger<unsigned long>(64);
    srcDims[2] = gPRNG.GetInteger<unsigned long>(32);
    srcDims[3] = gPRNG.GetInteger<unsigned long>(32);
    srcDims[4] = gPRNG.GetInteger<unsigned long>(32);
    
    std::vector<int> src, result;
    src.resize(srcDims.ElementCount());
    for( auto it=src.begin(); it!=src.end(); it++ ) {
      *it = 1;
    }
    
    kvl::cuda::CudaImage<int,5,size_t> d_image;
    
    d_image.Send( src, srcDims );
    d_image.SetMemory(0);
    d_image.Recv( result, resultDims );
    
    BOOST_CHECK_EQUAL( srcDims, resultDims );
    for( auto it=result.begin(); it!=result.end(); it++ ) {
      BOOST_CHECK_EQUAL( *it, 0 );
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
