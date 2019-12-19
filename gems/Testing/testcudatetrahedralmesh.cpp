#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/mpl/list.hpp>

#include "cudatetrahedralmesh.hpp"

#include "testiosupport.hpp"
#include "testfileloader.hpp"

// ========================

typedef boost::mpl::list<float,double> CoordTypes;

// ========================

BOOST_AUTO_TEST_SUITE( CudaTetrahedralMesh )

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE_TEMPLATE( SendToGPU, T, CoordTypes )
{
  kvl::cuda::CudaTetrahedralMesh<T,unsigned long,float> ctm;

  ctm.Send( mesh );
}

BOOST_AUTO_TEST_CASE( MeshIndexTypeTooSmall )
{
  kvl::cuda::CudaTetrahedralMesh<double,unsigned char,float> ctm;

  BOOST_CHECK_THROW( ctm.Send(mesh), std::out_of_range );
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
