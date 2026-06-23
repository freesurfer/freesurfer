#pragma once

#include <boost/test/unit_test.hpp>

#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"
#include "itkImageFileReader.h"


class TestFileLoader {
public:
  typedef kvl::AtlasMeshToIntensityImageCostAndGradientCalculator::ImageType ImageType;

  ImageType::ConstPointer image;
  kvl::AtlasMeshCollection::Pointer meshCollection;
  kvl::AtlasMesh::ConstPointer mesh;

  //! Constructor reads in the two files
  TestFileLoader() {
    // Read a test image
    typedef itk::ImageFileReader<ImageType>  ReaderType;
    ReaderType::Pointer  reader = ReaderType::New();
    reader->SetFileName( "test.nii" );
    reader->Update();
    this->image = reader->GetOutput();
    BOOST_TEST_CHECKPOINT("Test image read");
    
    // Read a test atlas mesh
    this->meshCollection = kvl::AtlasMeshCollection::New();
    BOOST_CHECK( meshCollection->Read( "test.txt" ) );
    this->mesh = this->meshCollection->GetReferenceMesh();
    BOOST_TEST_CHECKPOINT("Test atlas mesh read");
  }
};
