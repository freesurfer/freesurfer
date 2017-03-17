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
  typedef MeshSource::IdentifierType  IdentifierType;
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

  // Create mesh with single tetrahedron
  MeshSource::Pointer meshSource = MeshSource::New();

  // Define tetrahedron enclosing origin and (0,0,1), (0,1,0) and (1,0,0)
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  const float  p0[] = { -delta, -delta, -delta };
  const float  p1[] = { 1+(2*delta),- delta, -delta };
  const float  p2[] = { -delta, 1+(2*delta), -delta };
  const float  p3[] = { -delta, -delta, 1+(2*delta) };
  const IdentifierType  id0 = meshSource->AddPoint( p0 );
  const IdentifierType  id1 = meshSource->AddPoint( p1 );
  const IdentifierType  id2 = meshSource->AddPoint( p2 );
  const IdentifierType  id3 = meshSource->AddPoint( p3 );
  meshSource->AddTetrahedron( id0, id1, id2, id3 );

  Mesh::Pointer mesh = meshSource->GetOutput();

  BOOST_TEST_CHECKPOINT("Mesh created");

  kvl::AtlasMeshVisitCounterCPU::Pointer visitCounter = kvl::AtlasMeshVisitCounterCPU::New();

  visitCounter->SetRegions( image->GetLargestPossibleRegion() );
  visitCounter->Rasterize( mesh );

  BOOST_TEST_CHECKPOINT("VisitCounter Complete");

  // Check points in tetrahedron
  const ImageType* result = visitCounter->GetImage();
  for( int k=0; k<nz; k++ ) {
    for( int j=0; j<ny; j++ ) {
      for( int i=0; i<nx; i++ ) {
	ImageType::IndexType idx;
	idx[0] = i;
	idx[1] = j;
	idx[2] = k;

	// Recall that the tetrahedron is constructed to be the lower one
	// enclosing points with at most one index equal to 1
	int expected = 0;
	if( i+j+k <= 1 ) {
	  expected = 1;
	}

	BOOST_TEST_INFO( "(" << i << "," << j << "," << k << ")" );
	BOOST_CHECK_EQUAL( result->GetPixel(idx), expected );
      }
    }
  }
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
