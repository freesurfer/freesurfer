#include <boost/test/unit_test.hpp>

#include "itkImageRegionConstIteratorWithIndex.h"

#include "kvlAtlasMesh.h"
#include "atlasmeshvisitcounter.hpp"
#include "atlasmeshvisitcountercpuwrapper.hpp"
#ifdef CUDA_FOUND
#include "atlasmeshvisitcountercuda.hpp"
#endif

#include "testfileloader.hpp"

// --------------------

void CheckVisitCounter( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
			ImageType::ConstPointer targetImage,
			kvl::AtlasMesh::ConstPointer targetMesh ) {
  kvl::AtlasMeshVisitCounterCPU::Pointer  originalVisitCounter = kvl::AtlasMeshVisitCounterCPU::New();

  originalVisitCounter->SetRegions( targetImage->GetLargestPossibleRegion() );
  originalVisitCounter->Rasterize( targetMesh );
  BOOST_TEST_MESSAGE("Original counter complete");
  
  visitCounter->SetRegions( targetImage->GetLargestPossibleRegion() );
  visitCounter->VisitCount( targetMesh );
  BOOST_TEST_MESSAGE("AtlasMeshVisitCounterCPUWrapper complete");
  
  itk::ImageRegionConstIteratorWithIndex<kvl::interfaces::AtlasMeshVisitCounter::ImageType>  
    it( visitCounter->GetImage(), visitCounter->GetImage()->GetBufferedRegion() );
  itk::ImageRegionConstIteratorWithIndex<kvl::AtlasMeshVisitCounterCPU::ImageType>  
    itOrig( originalVisitCounter->GetImage(), originalVisitCounter->GetImage()->GetBufferedRegion() );
  
  for( ; !it.IsAtEnd(); ++it, ++itOrig ) {
    BOOST_TEST_CONTEXT( "Voxel Index: " << it.GetIndex() ) {
      BOOST_CHECK_EQUAL( it.Value(), itOrig.Value() );
    }
  }
}

// -------------------

BOOST_AUTO_TEST_SUITE( AtlasMeshVisitCounter )

BOOST_AUTO_TEST_SUITE( UnitCubeSingleTetrahedron )

BOOST_AUTO_TEST_CASE( OriginOnly )
{
  // Create unit cube
  kvl::interfaces::AtlasMeshVisitCounter::ImageType::RegionType region;
  kvl::interfaces::AtlasMeshVisitCounter::ImageType::IndexType start;
  kvl::interfaces::AtlasMeshVisitCounter::ImageType::SizeType size;
  start[0] = start[1] = start[2] = 0;
  size[0] = size[1] = size[2] = 1;

  region.SetSize(size);
  region.SetIndex(start);

  kvl::interfaces::AtlasMeshVisitCounter::ImageType::Pointer image;
  image->SetRegions(region);
  image->Allocate();

  for( int k=0; k<2; k++ ) {
    for( int j=0; j<2; j++ ) {
      for( int i=0; i<2; i++ ) {
	kvl::interfaces::AtlasMeshVisitCounter::ImageType::IndexType idx;
	idx[0] = k;
	idx[1] = j;
	idx[2] = i;

	image->SetPixel(idx, 0);
      }
    }
  }

  // Create mesh with single triangle
  kvl::AtlasMesh::Pointer mesh = kvl::AtlasMesh::New();
}

BOOST_AUTO_TEST_SUITE_END();

// --

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE( ReferenceImpl )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
 
  // Note that image and mesh are supplied by TestFileLoader
  CheckVisitCounter( &visitCounter, image, mesh );
}

#ifdef CUDA_FOUND
BOOST_AUTO_TEST_CASE( CUDAImpl )
{
  kvl::cuda::AtlasMeshVisitCounterCUDA visitCounter;
 
  // Note that image and mesh are supplied by TestFileLoader
  CheckVisitCounter( &visitCounter, image, mesh );
}
#endif

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
