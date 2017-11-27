#include <pybind11/pybind11.h>
#include "kvlCroppedImageReader.h"
#include "itkCastImageFilter.h"

namespace py = pybind11;


void kvlReadImage( const std::string imageFileName ) {
    // Read the image
    kvl::CroppedImageReader::Pointer reader = kvl::CroppedImageReader::New();
    reader->Read( imageFileName.c_str() );

    // Convert the image to float
    typedef itk::Image< float, 3 >  ImageType;
    typedef itk::CastImageFilter< kvl::CroppedImageReader::ImageType, ImageType >  CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( reader->GetImage() );
    caster->Update();


    // Store the image and transform in persistent memory
    const auto imageHandle = caster->GetOutput();
    typedef kvl::CroppedImageReader::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    reader->GetWorldToImageTransform()->GetInverse( transform );
    std::cout << "Read image: " << imageFileName << std::endl;
}


PYBIND11_MODULE(GEMS2Python, m) {
    m.def("kvlReadImage", &kvlReadImage);
}
