#include <functional>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/mpl/list.hpp>

#include "itkImageRegionConstIteratorWithIndex.h"

#include "kvlAtlasMesh.h"
#include "atlasmeshvisitcounter.hpp"
#include "atlasmeshvisitcountercpuwrapper.hpp"
#ifdef CUDA_FOUND
#include "cudaimage.hpp"
#include "atlasmeshvisitcountercuda.hpp"
#include "visitcountersimplecuda.hpp"
#include "visitcountertetrahedralmeshcuda.hpp"
#endif

#include "imageutils.hpp"
#include "testfileloader.hpp"
#include "testiosupport.hpp"

#ifdef CUDA_FOUND
#ifdef GPU_ALL_PRECISIONS
typedef boost::mpl::list<
  kvl::cuda::VisitCounterSimple<float,float>,
  kvl::cuda::VisitCounterSimple<double,double>,
  kvl::cuda::VisitCounterSimple<float,double>,
  kvl::cuda::VisitCounterTetrahedralMesh
  > CUDAImplTypes;
#else
typedef boost::mpl::list<
  kvl::cuda::VisitCounterSimple<double,double>,
  kvl::cuda::VisitCounterTetrahedralMesh
  > CUDAImplTypes;
#endif
#endif

// --------------------

const int nDims = 3;
const int nVertices = 4;
const int nAlphas = 1;

// --------------------

static std::ostream& operator<<( std::ostream& os,
				 const float v[nVertices][nDims] ) {
  os << "[";

  for( unsigned int j=0; j<nVertices; j++ ) {
    os << "(" << v[j][0];
    for( unsigned int i=1; i<nDims; i++ ) {
      os << "," << v[j][i];
    }
    os << ")";
    if( j!=(nVertices-1) ) {
      os << ",";
    }
  }
  os << "]";
  
  return os;
}

static std::string TetrahedronToString( const float v[nVertices][nDims] ) {
  std::stringstream res;

  res << v;
  return res.str();
}

// --------------------

void CheckVisitCounter( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
			TestFileLoader::ImageType::ConstPointer targetImage,
			kvl::AtlasMesh::ConstPointer targetMesh ) {
  kvl::AtlasMeshVisitCounter::Pointer  originalVisitCounter = kvl::AtlasMeshVisitCounter::New();

  originalVisitCounter->SetRegions( targetImage->GetLargestPossibleRegion() );
  originalVisitCounter->Rasterize( targetMesh );
  BOOST_TEST_MESSAGE("Original counter complete");
  
  visitCounter->SetRegions( targetImage->GetLargestPossibleRegion() );
  visitCounter->VisitCount( targetMesh );
  BOOST_TEST_MESSAGE("Target Visit Counter complete");

  auto img = visitCounter->GetImage();
  itk::ImageRegionConstIteratorWithIndex<kvl::interfaces::AtlasMeshVisitCounter::ImageType>  
    it( img, img->GetBufferedRegion() );
  itk::ImageRegionConstIteratorWithIndex<kvl::AtlasMeshVisitCounter::ImageType>  
    itOrig( originalVisitCounter->GetImage(), originalVisitCounter->GetImage()->GetBufferedRegion() );
  
  for( ; !it.IsAtEnd(); ++it, ++itOrig ) {
    BOOST_TEST_CONTEXT( "Voxel Index: " << it.GetIndex() ) {
      BOOST_CHECK_EQUAL( it.Value(), itOrig.Value() );
    }
  }
}

// -------------------

typedef kvl::interfaces::AtlasMeshVisitCounter::ImageType ImageType;
typedef itk::AutomaticTopologyMeshSource<kvl::AtlasMesh> MeshSource;
typedef MeshSource::IdentifierType  IdentifierType;
typedef kvl::AtlasMesh Mesh;


// ----------------------------

void SingleTetrahedronUnitMesh( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
				float vertices[nVertices][nDims],
				std::function<int(int,int,int)> expectedCount ) {
  const int imageSize = 2;
  const int nx = imageSize;
  const int ny = imageSize;
  const int nz = imageSize;

  ImageType::Pointer image = kvl::Testing::CreateImageCube<ImageType>( imageSize, 0 );
  BOOST_TEST_CHECKPOINT("Image created");

  Mesh::Pointer mesh = kvl::Testing::CreateSingleTetrahedronMesh( vertices, nAlphas );
  BOOST_TEST_CHECKPOINT("Mesh created");

  visitCounter->SetRegions( image->GetLargestPossibleRegion() );
  visitCounter->VisitCount( mesh );

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

	BOOST_TEST_INFO( "(" << i << "," << j << "," << k << ")" );
	BOOST_CHECK_EQUAL( result->GetPixel(idx), expectedCount(i,j,k) );
      }
    }
  }
}

// -------------------

void LowerCorner( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron enclosing origin and (0,0,1), (0,1,0) and (1,0,0)
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { -delta, -delta, -delta },
    { 1+(2*delta), -delta, -delta },
    { -delta, 1+(2*delta), -delta },
    { -delta, -delta, 1+(2*delta) }
  };

  auto expectedCount = [](int i, int j, int k) {
    // Recall that the tetrahedron is constructed to be the lower one
    // enclosing points with at most one index equal to 1
    if( i+j+k <= 1 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void OriginOnly( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
   // Define tetrahedron enclosing origin only
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { -delta, -delta, -delta },
    { delta, 0, 0 },
    { 0, delta, 0 },
    { 0, 0, delta }
  };

  auto expectedCount = [](int i, int j, int k) {
    if( i+j+k == 0 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void XAxisOnly( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron enclosing origin and (1,0,0)
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { -delta, -delta, -delta },
    { 1+(2*delta), 0, 0 },
    { -delta, 2*delta, -delta },
    { -delta, -delta, 2*delta }
  };

  auto expectedCount = [](int i, int j, int k) {
    if( j+k == 0 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void FarCornerOnly( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron enclosing (1,1,1)
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { 1-delta, 1-delta, 1-delta },
    { 1+delta, 1, 1 },
    { 1, 1+delta, 1 },
    { 1, 1, 1+delta }
  };
  
  auto expectedCount = [](int i, int j, int k) {
    if( i+j+k == 3 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void UpperCornerOnly( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron enclosing (1,1,1), (0,1,1), (1,0,1), (1,1,0)
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { 1+delta, 1+delta, 1+delta },
    { -delta, 1, 1 },
    { 1, -delta, 1 },
    { 1, 1, -delta }
  };
  
  auto expectedCount = [](int i, int j, int k) {
    if( i+j+k >= 2 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void NoVertices( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron enclosing no vertices
  const float delta = 0.1f;
  BOOST_REQUIRE( delta < 0.5f );
  float verts[nVertices][nDims] = {
    { 0.5f+delta, 0.5f+delta, 0.5f+delta },
    { 0.5f-delta, 0.5f, 0.5f },
    { 0.5f, 0.5f-delta, 0.5f },
    { 0.5f, 0.5f, -delta }
  };

  auto expectedCount = [](int i, int j, int k) { return 0; };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void LowerCornerExact( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron in on exactly (0,0,0), (1,0,0), (0,1,0) and (0,0,1)
  // This is to go after some of the edge cases
  float verts[nVertices][nDims] = {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 0, 1, 0 },
    { 0, 0, 1 }
  };

  auto expectedCount = [](int i, int j, int k) {
    if( i+j+k == 0 ) {
      return 1;
    } else {
      return 0;
    }
  };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

void UpperCornerExact( kvl::interfaces::AtlasMeshVisitCounter* visitCounter ) {
  // Define tetrahedron in on exactly (1,1,1), (0,1,1), (1,0,1) and (1,1,0)
  // This is to go after some of the edge cases
  float verts[nVertices][nDims] = {
    { 1, 1, 1 },
    { 0, 1, 1 },
    { 1, 0, 1 },
    { 1, 1, 0 }
  };

  auto expectedCount = [](int i, int j, int k) { return 0; };

  SingleTetrahedronUnitMesh( visitCounter, verts, expectedCount );
}

// -------------------------------

void GenerateSpecificCornerTetrahedron( float verts[nVertices][nDims],
					const unsigned char corner,
					const float scale ) {
  // Generate a 'corner' tetrahedron, where the apex is specified by the bits of the 'corner' argument

  // Separate out the bits specifying the corner
  unsigned char mask = 4;
  for( unsigned int i=0; i<nDims; i++ ) {
    if( mask & corner ) {
      verts[0][i] = 1;
    } else {
      verts[0][i] = 0;
    }
    mask = mask >> 1;
  }
  
  // Compute the rest of the vertices, each has one index changed from the apex
  for( unsigned int j=1; j<nVertices; j++ ) {
    for( unsigned int i=0; i<nDims; i++ ) {
      if( i==(j-1) ) {
	verts[j][i] = 1 - verts[0][i];
      } else {
	verts[j][i] = verts[0][i];
      }
    }
  }

  // Scale everything
  for( unsigned int j=0; j<nVertices; j++ ) {
    for( unsigned int i=0; i<nDims; i++ ) {
      verts[j][i] *= scale;
    }
  }
}

ImageType::ConstPointer ApplyVisitCounterToMesh( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
						 const ImageType* targetImage,
						 Mesh::Pointer targetMesh ) {
  visitCounter->SetRegions( targetImage->GetLargestPossibleRegion() );
  visitCounter->VisitCount( targetMesh );
  BOOST_TEST_CHECKPOINT("VisitCount run");

  return visitCounter->GetImage();
}

void CompareVisitImages( const ImageType* standard,
			 const ImageType* compare ) {
  BOOST_CHECK( standard != NULL );
  BOOST_CHECK( compare != NULL );
  itk::ImageRegionConstIteratorWithIndex<kvl::interfaces::AtlasMeshVisitCounter::ImageType>  
    it( compare, compare->GetBufferedRegion() );
  itk::ImageRegionConstIteratorWithIndex<kvl::AtlasMeshVisitCounter::ImageType>  
    itOrig( standard, standard->GetBufferedRegion() );
  
  for( ; !it.IsAtEnd(); ++it, ++itOrig ) {
    BOOST_TEST_CONTEXT( "Voxel Index: " << it.GetIndex() ) {
      BOOST_CHECK_EQUAL( it.Value(), itOrig.Value() );
    }
  }
}

void CheckVisitCounterWithPermutations( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
					const ImageType* targetImage,
					const float tetrahedron[nVertices][nDims] ) {
  // Start by getting the 'standard' answers
  BOOST_TEST_CHECKPOINT("Starting CheckVisitCounterWithPermutations");
  ImageType::ConstPointer standardVisit = NULL;
  {
    kvl::AtlasMeshVisitCounterCPUWrapper origVisitCounter;
    BOOST_TEST_CHECKPOINT("Created reference VisitCounter");
    Mesh::Pointer baseMesh = kvl::Testing::CreateSingleTetrahedronMesh( tetrahedron, nAlphas );
    BOOST_TEST_CHECKPOINT("baseMesh Created");
    standardVisit = ApplyVisitCounterToMesh( &origVisitCounter, targetImage, baseMesh );
  }
  BOOST_TEST_CHECKPOINT("Created standardVisit image");

  // Now work on permuting the vertices of the tetrahedron

  // Create the permutation array
  std::vector<unsigned int> perm;
  for( unsigned int j=0; j<nVertices; j++ ) {
    perm.push_back(j);
  }

  // Iterate over the permutations
  do {
    // Create permuted tetrahedron
    float permTet[nVertices][nDims];
    for( unsigned int j=0; j<nVertices; j++ ) {
      for( unsigned int i=0; i<nDims; i++ ) {
	permTet[perm[j]][i] = tetrahedron[j][i];
      }
    }

    BOOST_TEST_CONTEXT( "Tetrahedron : " << TetrahedronToString(permTet) ) {
      // Generate the result image
      Mesh::Pointer mesh = kvl::Testing::CreateSingleTetrahedronMesh( permTet, nAlphas );
      BOOST_TEST_CHECKPOINT("Permuted mesh created");

      ImageType::ConstPointer currVisit = NULL;
      currVisit = ApplyVisitCounterToMesh( visitCounter, targetImage, mesh );
      BOOST_TEST_CHECKPOINT("Created currVisit image");

      CompareVisitImages( standardVisit, currVisit );
    }
  } while( std::next_permutation( perm.begin(), perm.end() ) );
} 

void AutoCorners( kvl::interfaces::AtlasMeshVisitCounter* visitCounter,
		  const float scaleTetrahedron,
		  const int imageSize ) {
  float baseVertices[nVertices][nDims];

  const unsigned char nCorners = 8;

  for( unsigned char corner=0; corner<nCorners; corner++ ) {
    BOOST_TEST_CONTEXT( "Corner : " << static_cast<unsigned int>(corner) ) {
      // Get the tetrahedron we want to test
      GenerateSpecificCornerTetrahedron( baseVertices, corner, scaleTetrahedron );
      BOOST_TEST_CHECKPOINT("Generated corner tetrahedron");

      // Add one to imageSize since we're specifying the number of points
      // on each edge
      ImageType::Pointer targetImage = kvl::Testing::CreateImageCube<ImageType>(imageSize+1,0);
      BOOST_TEST_CHECKPOINT("Generated target image");

      CheckVisitCounterWithPermutations( visitCounter, targetImage, baseVertices );
    }
  }
}

// ===========================================================

BOOST_AUTO_TEST_SUITE( AtlasMeshVisitCounter )

BOOST_AUTO_TEST_SUITE( SingleTetrahedron )

BOOST_AUTO_TEST_CASE( LowerCornerCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
 
  LowerCorner( &visitCounter );
}

BOOST_AUTO_TEST_CASE( OriginOnlyCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  OriginOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE( XAxisOnlyCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  XAxisOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE( FarCornerOnlyCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  FarCornerOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE( UpperCornerCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  UpperCornerOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE( NoVerticesCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  NoVertices( &visitCounter );
}

BOOST_AUTO_TEST_CASE( LowerCornerExactCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  LowerCornerExact( &visitCounter );
}

BOOST_AUTO_TEST_CASE( UpperCornerExactCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
  
  UpperCornerExact( &visitCounter );
}

BOOST_AUTO_TEST_CASE( AutoCornersBasicCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;

  // Generate unit tetrahedron and unit cube
  AutoCorners( &visitCounter, 1, 1  );
}

BOOST_AUTO_TEST_CASE( AutoCornersLargeImageLargeTetrahedronCPU )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;

  // Generate large tetrahedron on large cube
  AutoCorners( &visitCounter, 4, 4  );
}

#ifdef CUDA_FOUND
BOOST_AUTO_TEST_CASE_TEMPLATE( LowerCornerGPUSimple, ImplType, CUDAImplTypes  )
{
  ImplType visitCounter;
 
  LowerCorner( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( OriginOnlyGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  OriginOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( XAxisOnlyGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  XAxisOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( FarCornerOnlyGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  FarCornerOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( UpperCornerGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  UpperCornerOnly( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( NoVerticesGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  NoVertices( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( LowerCornerExactGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  LowerCornerExact( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( UpperCornerExactGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  UpperCornerExact( &visitCounter );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( AutoCornersGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  AutoCorners( &visitCounter, 1, 1 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( AutoCornersLargeImageLargeTetrahedronGPUSimple, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;

  // Generate large tetrahedron on large cube
  AutoCorners( &visitCounter, 4, 4  );
}
#endif

BOOST_DATA_TEST_CASE( ConsistencyCheck, boost::unit_test::data::xrange(1,2), scale )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;

  float unitVertices[nVertices][nDims];
  float scaleVertices[nVertices][nDims];

  GenerateSpecificCornerTetrahedron( unitVertices, 0, 1 );
  GenerateSpecificCornerTetrahedron( scaleVertices, 0, scale );
  BOOST_TEST_CHECKPOINT("Created tetrahedra");

  Mesh::Pointer unitMesh = kvl::Testing::CreateSingleTetrahedronMesh( unitVertices, nAlphas );
  Mesh::Pointer scaleMesh = kvl::Testing::CreateSingleTetrahedronMesh( scaleVertices, nAlphas );
  BOOST_TEST_CHECKPOINT("Created meshes");

  ImageType::Pointer targetImage = kvl::Testing::CreateImageCube<ImageType>(scale+1,0);
  BOOST_TEST_CHECKPOINT("Created target image");

  ImageType::ConstPointer unitVisit = ApplyVisitCounterToMesh( &visitCounter, targetImage, unitMesh );
  const ImageType* unitResult = visitCounter.GetImage();
  ImageType::ConstPointer scaleVisit = ApplyVisitCounterToMesh( &visitCounter, targetImage, scaleMesh );
  const ImageType* scaleResult = visitCounter.GetImage();
  BOOST_TEST_CHECKPOINT("VisitCounters complete");

  // Check the vertices of the image cube
  for( int k=0; k<2; k++ ) {
    for( int j=0; j<2; j++ ) {
      for( int i=0; i<2; i++ ) {
	ImageType::IndexType idx, idxScale;
	idx[0] = i;
	idx[1] = j;
	idx[2] = k;

	for( unsigned iDim=0; iDim<nDims; iDim++ ) {
	  idxScale[iDim] = idx[iDim]*scale;
	}

	BOOST_TEST_INFO( "(" << i << "," << j << "," << k << ")" );
	BOOST_TEST_INFO( "Unit Tetrahedron : " << TetrahedronToString(unitVertices) );
	BOOST_TEST_INFO( "Scale Tetrahedron : " << TetrahedronToString(scaleVertices) );
	BOOST_CHECK_EQUAL( unitResult->GetPixel(idx), scaleResult->GetPixel(idxScale) );
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();

// ---------------------------------------

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE( ReferenceImpl )
{
  kvl::AtlasMeshVisitCounterCPUWrapper visitCounter;
 
  // Note that image and mesh are supplied by TestFileLoader
  CheckVisitCounter( &visitCounter, image, mesh );
  
  BOOST_TEST_MESSAGE( "SetRegions Time  : " << visitCounter.tSetRegions );
  BOOST_TEST_MESSAGE( "VisitCounter Time: " << visitCounter.tVisitCount );
}

#ifdef CUDA_FOUND
BOOST_AUTO_TEST_CASE_TEMPLATE( SimpleCUDAImpl, ImplType, CUDAImplTypes )
{
  ImplType visitCounter;
 
  // Note that image and mesh are supplied by TestFileLoader
  CheckVisitCounter( &visitCounter, image, mesh );

  
  BOOST_TEST_MESSAGE( "SetRegions Time  : " << visitCounter.tSetRegions );
  BOOST_TEST_MESSAGE( "VisitCounter Time: " << visitCounter.tVisitCount );
  BOOST_TEST_MESSAGE( "       Pack : " << visitCounter.tVisitCountPack );
  BOOST_TEST_MESSAGE( "   Transfer : " << visitCounter.tVisitCountTransfer );
  BOOST_TEST_MESSAGE( "     Kernel : " << visitCounter.tVisitCountKernel );
  BOOST_TEST_MESSAGE( "GetImage Time    : " << visitCounter.tGetImage );
  BOOST_TEST_MESSAGE( "   Transfer : " << visitCounter.tGetImageTransfer );
  BOOST_TEST_MESSAGE( "     Unpack : " << visitCounter.tGetImageUnpack );
}
#endif

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
