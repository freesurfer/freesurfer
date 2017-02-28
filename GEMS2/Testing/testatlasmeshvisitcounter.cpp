#include <boost/test/unit_test.hpp>

#include "itkImageRegionConstIteratorWithIndex.h"

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

BOOST_FIXTURE_TEST_SUITE( AtlasMeshVisitCounter, TestFileLoader )

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
