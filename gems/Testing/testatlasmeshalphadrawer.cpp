#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "itkImageRegionConstIteratorWithIndex.h"

#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshAlphaDrawer.h"
#include "atlasmeshalphadrawer.hpp"
#include "atlasmeshalphadrawercpuwrapper.hpp"

#ifdef CUDA_FOUND
#include "cudaimage.hpp"
#include "atlasmeshalphadrawercuda.hpp"
#endif

#include "imageutils.hpp"
#include "testfileloader.hpp"
#include "testiosupport.hpp"

// ----------------------------------------------

const int nDims = 3;
const int nVertices = 4;

// -----------------------------------------

typedef kvl::interfaces::AtlasMeshAlphaDrawer::ImageType ImageType;
typedef kvl::AtlasMesh Mesh;

// ----------------------------------------------

void CheckAlphaDrawer( kvl::interfaces::AtlasMeshAlphaDrawer* ad,
		       TestFileLoader::ImageType::ConstPointer targetImage,
		       kvl::AtlasMesh::ConstPointer targetMesh,
		       const int classNumber,
		       const float percentTolerance ) {
  kvl::AtlasMeshAlphaDrawer::Pointer originalAD = kvl::AtlasMeshAlphaDrawer::New();

  originalAD->SetRegions( targetImage->GetLargestPossibleRegion() );
  originalAD->SetClassNumber( classNumber );
  originalAD->Rasterize( targetMesh );
  BOOST_TEST_MESSAGE("Original Alpha Drawer Complete");

  ad->SetRegions( targetImage->GetLargestPossibleRegion() );
  ad->SetClassNumber( classNumber );
  ad->Interpolate( targetMesh );
  BOOST_TEST_MESSAGE("Target Alpha Drawer Complete");

  auto img = ad->GetImage();
  itk::ImageRegionConstIteratorWithIndex<kvl::interfaces::AtlasMeshAlphaDrawer::ImageType>  
    it( img, img->GetBufferedRegion() );
  itk::ImageRegionConstIteratorWithIndex<kvl::AtlasMeshAlphaDrawer::ImageType>  
    itOrig( originalAD->GetImage(), originalAD->GetImage()->GetBufferedRegion() );
  
  for( ; !it.IsAtEnd(); ++it, ++itOrig ) {
    BOOST_TEST_CONTEXT( "Voxel Index: " << it.GetIndex() ) {
      // Crude test for small values
      if( fabs(itOrig.Value()) > percentTolerance ) {
	// If 'large' do a regular 'percentage' check
	BOOST_CHECK_CLOSE( it.Value(), itOrig.Value(), percentTolerance );
      } else {
	// Otherwise do absolute size, re-using the percentTolerance
	BOOST_CHECK_SMALL( it.Value(), percentTolerance );
      }
    }
  }
}


void SingleConstantTetrahedronContainedCube( kvl::interfaces::AtlasMeshAlphaDrawer* ad,
					     const int classNumber,
					     const int nAlphas,
					     const int imageSize ) {
  BOOST_REQUIRE( nAlphas > 0 );
  BOOST_REQUIRE( classNumber < nAlphas );
  BOOST_REQUIRE( imageSize > 1 );
  const float d = imageSize;

  // Set floating point tolerance as a percentage
  // A value of 1.0 means 1%
  const float percentTolerance = 0.0001;
  
  ImageType::Pointer image = kvl::Testing::CreateImageCube<ImageType>( imageSize, 0 );
  BOOST_TEST_CHECKPOINT("Image created");

  float verts[nVertices][nDims] = {
    { -1 , -1 , -1  },
    { 4*d, -1 , -1  },
    { -1 , 4*d, -1  },
    { -1 , -1 , 4*d }
  };

  Mesh::Pointer mesh = kvl::Testing::CreateSingleTetrahedronMesh( verts, nAlphas );
  BOOST_TEST_CHECKPOINT("Mesh created");

  ad->SetRegions( image->GetLargestPossibleRegion() );
  ad->SetClassNumber( classNumber );
  ad->Interpolate( mesh );
  BOOST_TEST_CHECKPOINT("AlphaDrawer complete");

  auto img = ad->GetImage();
  for( unsigned int k=0; k<imageSize; k++ ) {
    for( unsigned int j=0; j<imageSize; j++ ) {
      for( unsigned int i=0; i<imageSize; i++ ) {
	ImageType::IndexType idx;
	idx[0] = i;
	idx[1] = j;
	idx[2] = k;

	BOOST_TEST_INFO( "(" << i << "," << j << "," << k << ")" );
	float pxlValue = img->GetPixel(idx);
	BOOST_CHECK_CLOSE( img->GetPixel(idx), static_cast<float>(classNumber), percentTolerance );
      }
    }
  }
}


// ==========================================

BOOST_AUTO_TEST_SUITE( AtlasMeshAlphaDrawer )

BOOST_AUTO_TEST_SUITE( SingleTetrahedron )

const int nAlphas = 11;

BOOST_DATA_TEST_CASE( ContainedUnitCube,  boost::unit_test::data::xrange(nAlphas), classNumber )
{
  kvl::AtlasMeshAlphaDrawerCPUWrapper ad;

  SingleConstantTetrahedronContainedCube( &ad, classNumber, nAlphas, 2 );
}

BOOST_DATA_TEST_CASE( ContainedLargeCube,  boost::unit_test::data::xrange(nAlphas), classNumber )
{
  kvl::AtlasMeshAlphaDrawerCPUWrapper ad;

  SingleConstantTetrahedronContainedCube( &ad, classNumber, nAlphas, 23 );
}

#ifdef CUDA_FOUND
BOOST_DATA_TEST_CASE( ContainedUnitCubeGPU,  boost::unit_test::data::xrange(nAlphas), classNumber )
{
  kvl::cuda::AtlasMeshAlphaDrawerCUDA ad;

  SingleConstantTetrahedronContainedCube( &ad, classNumber, nAlphas, 2 );
}

BOOST_DATA_TEST_CASE( ContainedLargeCubeGPU,  boost::unit_test::data::xrange(nAlphas), classNumber )
{
  kvl::cuda::AtlasMeshAlphaDrawerCUDA ad;

  SingleConstantTetrahedronContainedCube( &ad, classNumber, nAlphas, 23 );
}
#endif

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE( ReferenceImpl )
{
  kvl::AtlasMeshAlphaDrawerCPUWrapper ad;
  const int classNumber = 1;

  // Set floating point tolerance as a percentage
  // A value of 1.0 means 1%
  const float percentTolerance = 0;

  // Note that image and mesh are supplied by TestFileLoader
  CheckAlphaDrawer( &ad, image, mesh, classNumber, percentTolerance );
  
  BOOST_TEST_MESSAGE( "SetRegions Time           : " << ad.tSetRegions );
  BOOST_TEST_MESSAGE( "Interpolate Time          : " << ad.tInterpolate );
  ad.tInterpolate.Reset();

  CheckAlphaDrawer( &ad, image, mesh, classNumber, percentTolerance );
  BOOST_TEST_MESSAGE( "Interpolate Time (repeat) : " << ad.tInterpolate );
}

#ifdef CUDA_FOUND
BOOST_AUTO_TEST_CASE( CudaImpl )
{
  kvl::cuda::AtlasMeshAlphaDrawerCUDA ad;
  const int classNumber = 1;

  // Set floating point tolerance as a percentage
  // A value of 1.0 means 1%
  const float percentTolerance = 0.0001;

  // Note that image and mesh are supplied by TestFileLoader
  CheckAlphaDrawer( &ad, image, mesh, classNumber, percentTolerance );
  
  BOOST_TEST_MESSAGE( "SetRegions Time           : " << ad.tSetRegions );
  BOOST_TEST_MESSAGE( "Interpolate Time          : " << ad.tInterpolate );
  BOOST_TEST_MESSAGE( "   Send Mesh Time         : " << ad.tSendMesh );
  BOOST_TEST_MESSAGE( "   Kernel Time            : " << ad.tKernel );
}
#endif

BOOST_AUTO_TEST_CASE( MeshInformation )
{
  size_t nOther = 0;
  std::vector<kvl::AtlasMesh::CellIdentifier> tetrahedronIds;

  for( auto cellIt = mesh->GetCells()->Begin();
       cellIt != mesh->GetCells()->End();
       ++cellIt ) {
    if( cellIt.Value()->GetType() == kvl::AtlasMesh::CellType::TETRAHEDRON_CELL ) {
      tetrahedronIds.push_back( cellIt.Index() );
    } else {
      nOther++;
    }
  }
  BOOST_TEST_MESSAGE("nTetrahedra : " << tetrahedronIds.size());
  BOOST_TEST_MESSAGE("nOther      : " << nOther);

  size_t minAlphas = std::numeric_limits<size_t>::max();
  size_t maxAlphas = 0;

  auto pointData = mesh->GetPointData();
  for( auto pointIt = pointData->Begin();
       pointIt != pointData->End();
       ++pointIt ) {
    size_t nAlphas = pointIt->Value().m_Alphas.size();
    if( nAlphas < minAlphas ) { minAlphas = nAlphas; }
    if( nAlphas > maxAlphas ) { maxAlphas = nAlphas; }
  }
  BOOST_TEST_MESSAGE("Total Points : " << pointData->size());
  BOOST_TEST_MESSAGE("minAlphas : " << minAlphas);
  BOOST_TEST_MESSAGE("maxAlphas : " << maxAlphas);
  
  std::set<size_t> pointIndices;
  for( int iTet=0; iTet<tetrahedronIds.size(); iTet++ ) {
    kvl::AtlasMesh::CellAutoPointer cell;
    mesh->GetCell( tetrahedronIds.at(iTet), cell );

    size_t pointCount = 0;
    for( auto pit = cell->PointIdsBegin();
	 pit != cell->PointIdsEnd();
	 ++pit ) {
      pointIndices.insert(*pit);
      pointCount++;
    }

    BOOST_CHECK_EQUAL(pointCount, 4);
  }

  BOOST_TEST_MESSAGE("Unique points referenced: " << pointIndices.size());
  BOOST_TEST_MESSAGE("Total Vertices          : " << tetrahedronIds.size() * 4);
}

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
