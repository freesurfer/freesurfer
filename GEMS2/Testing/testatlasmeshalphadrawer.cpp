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

BOOST_AUTO_TEST_SUITE_END();

BOOST_AUTO_TEST_SUITE_END();
