#include <boost/test/unit_test.hpp>

#include "itkImageRegionConstIteratorWithIndex.h"

#include "kvlAtlasMeshAlphaDrawer.h"
#include "atlasmeshalphadrawer.hpp"
#include "atlasmeshalphadrawercpuwrapper.hpp"

#include "testfileloader.hpp"
#include "testiosupport.hpp"

// ----------------------------------------------

void CheckAlphaDrawer( kvl::interfaces::AtlasMeshAlphaDrawer* ad,
		       TestFileLoader::ImageType::ConstPointer targetImage,
		       kvl::AtlasMesh::ConstPointer targetMesh,
		       const int classNumber ) {
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
      BOOST_CHECK_EQUAL( it.Value(), itOrig.Value() );
    }
  }
}


// ==========================================

BOOST_AUTO_TEST_SUITE( AtlasMeshAlphaDrawer )

BOOST_FIXTURE_TEST_SUITE( ActualImage, TestFileLoader )

BOOST_AUTO_TEST_CASE( ReferenceImpl )
{
  kvl::AtlasMeshAlphaDrawerCPUWrapper ad;
  const int classNumber = 1;

  // Note that image and mesh are supplied by TestFileLoader
  CheckAlphaDrawer( &ad, image, mesh, classNumber );
  
  BOOST_TEST_MESSAGE( "SetRegions Time           : " << ad.tSetRegions );
  BOOST_TEST_MESSAGE( "Interpolate Time          : " << ad.tInterpolate );
  ad.tInterpolate.Reset();

  CheckAlphaDrawer( &ad, image, mesh, classNumber );
  BOOST_TEST_MESSAGE( "Interpolate Time (repeat) : " << ad.tInterpolate );
}

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
