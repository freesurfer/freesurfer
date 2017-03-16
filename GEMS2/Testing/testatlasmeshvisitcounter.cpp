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
  const int nx = 2;
  const int ny = 2;
  const int nz = 2;
  // Create unit cube

  typedef kvl::interfaces::AtlasMeshVisitCounter::ImageType ImageType;
  typedef itk::AutomaticTopologyMeshSource<kvl::AtlasMesh> MeshSource;
  typedef kvl::AtlasMesh Mesh;

  ImageType::RegionType region;
  ImageType::IndexType start;
  ImageType::SizeType size;
  start[0] = start[1] = start[2] = 0;
  size[0] = nx-1;
  size[1] = ny-1;
  size[2] = nz-1;

  region.SetSize(size);
  region.SetIndex(start);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();

  BOOST_TEST_CHECKPOINT("Cube allocated");

  for( int k=0; k<nz; k++ ) {
    for( int j=0; j<ny; j++ ) {
      for( int i=0; i<nx; i++ ) {
	ImageType::IndexType idx;
	idx[0] = i;
	idx[1] = j;
	idx[2] = k;

	image->SetPixel(idx, 0);
      }
    }
  }

  BOOST_TEST_CHECKPOINT("Cube set");

  // Create mesh
  MeshSource::Pointer meshSource = MeshSource::New();

  Mesh::Pointer mesh = Mesh::New();
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
