#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/mpl/list.hpp>

// This will only be compiled if CUDA_FOUND is defined
#include "dimensioncuda.hpp"
#include "cudaimage.hpp"

#include "cudaimagetests.hpp"

#include "randomsupply.hpp"
#include "testiosupport.hpp"

typedef boost::mpl::list<unsigned char,char, unsigned short, short,unsigned int,int,unsigned long,long> IntegerTestTypes;

typedef boost::mpl::list<unsigned char, int, float, double> PlusTestTypes;

RandomSupply gPRNG;

template<typename T,int nDims>
void SendSetReceiveDriver() {
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<nDims,size_t> srcDims, resultDims;

    srcDims[0] = 128 + gPRNG.GetInteger<size_t>( 32 );
    for( int i=1; i<nDims; i++ ) {
      srcDims[i] = ( 2 << (9-nDims) ) + gPRNG.GetInteger<size_t>( 256 >> nDims );
    }

    std::vector<T> src, result;
    src.resize(srcDims.ElementCount());
    for( auto it=src.begin(); it!=src.end(); it++ ) {
      *it = 1;
    }
    BOOST_TEST_CHECKPOINT("Created src : " << srcDims);

    kvl::cuda::CudaImage<T,nDims,size_t> d_image;
    
    d_image.Send( src, srcDims );
    BOOST_TEST_CHECKPOINT("Send data : " << srcDims);
    d_image.SetMemory(0);
    BOOST_TEST_CHECKPOINT("Memory set");
    d_image.Recv( result, resultDims );
    BOOST_TEST_CHECKPOINT("Recv data : " << resultDims);
    
    BOOST_REQUIRE_EQUAL( srcDims, resultDims );
    BOOST_REQUIRE_EQUAL( resultDims.ElementCount(), result.size() );
    BOOST_REQUIRE_EQUAL( src.size(), result.size() );
    bool passed = true;
    for( auto it=result.begin(); it!=result.end(); it++ ) {
      passed = passed && ( *it == 0 );
    }
    BOOST_REQUIRE_EQUAL( passed, true );
  }
}

template<typename T,int nDims>
void PlusKernelDriver() {
   BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
     kvl::cuda::Dimension<nDims,size_t> srcDims, resultDims;
   
     const T value = gPRNG.GetInteger( 120 );

     srcDims[0] = 256 + gPRNG.GetInteger<size_t>( 512 );
     for( int i=1; i<nDims; i++ ) {
       srcDims[i] = ( 2 << (9-nDims) ) + gPRNG.GetInteger<size_t>( 512 >> nDims );
     }
     
     std::vector<T> src, result;
     src.resize(srcDims.ElementCount());
     for( auto it=src.begin(); it!=src.end(); it++ ) {
       *it = gPRNG.GetInteger<char>(120);
     }
     BOOST_TEST_CHECKPOINT("Created src : " << srcDims);
     
     kvl::cuda::CudaImage<T,nDims,size_t> d_src, d_dst;
     
     d_src.Send( src, srcDims );
     BOOST_TEST_CHECKPOINT("Send data : " << srcDims);
     runPlusTest( d_dst, d_src, value );
     BOOST_TEST_CHECKPOINT("Kernel complete");
     
     d_dst.Recv( result, resultDims );
     BOOST_TEST_CHECKPOINT("Recv data : " << resultDims);
    
     BOOST_REQUIRE_EQUAL( srcDims, resultDims );
     BOOST_REQUIRE_EQUAL( resultDims.ElementCount(), result.size() );
     BOOST_REQUIRE_EQUAL( src.size(), result.size() );
     bool passed = true;
     for( size_t idx = 0; idx < src.size(); idx++ ) {
       passed = passed && ( result.at(idx) == (value + src.at(idx)) );
     }
     BOOST_REQUIRE_EQUAL( passed, true );
   }
}



BOOST_AUTO_TEST_SUITE( CudaImage )

BOOST_AUTO_TEST_SUITE( SendSetReceive )

BOOST_AUTO_TEST_CASE_TEMPLATE( SendSetReceive1d, ElementType, IntegerTestTypes )
{
  SendSetReceiveDriver<ElementType,1>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SendSetReceive2d, ElementType, IntegerTestTypes )
{
  SendSetReceiveDriver<ElementType,2>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SendSetReceive3d, ElementType, IntegerTestTypes )
{
  SendSetReceiveDriver<ElementType,3>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SendSetReceive4d, ElementType, IntegerTestTypes )
{
  SendSetReceiveDriver<ElementType,4>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE( SendSetReceive5d, ElementType, IntegerTestTypes )
{
  SendSetReceiveDriver<ElementType,5>();
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE(IndexOperators)

BOOST_AUTO_TEST_CASE_TEMPLATE(PlusKernel1D, ElementType, PlusTestTypes )
{
  PlusKernelDriver<ElementType,1>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PlusKernel2D, ElementType, PlusTestTypes )
{
  PlusKernelDriver<ElementType,2>();
}

BOOST_AUTO_TEST_CASE_TEMPLATE(PlusKernel3D, ElementType, PlusTestTypes )
{
  PlusKernelDriver<ElementType,3>();
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
