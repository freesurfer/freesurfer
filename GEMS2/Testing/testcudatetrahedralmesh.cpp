#include <boost/test/unit_test.hpp>

#include "cudatetrahedralmesh.hpp"

#include "testiosupport.hpp"
#include "testfileloader.hpp"

// ========================

BOOST_AUTO_TEST_SUITE( CudaTetrahedralMesh )

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE( SendToGPU )
{
  kvl::cuda::CudaTetrahedralMesh<double,unsigned long> ctm;

  ctm.Send( mesh );
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
