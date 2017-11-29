#include "pyKvlImage.h"
#include "itkCastImageFilter.h"
#include <pybind11/numpy.h>

KvlImage::KvlImage(const std::string &imageFileName) {
    // Read the image
    kvl::CroppedImageReader::Pointer reader = kvl::CroppedImageReader::New();
    reader->Read( imageFileName.c_str() );

    // Convert the image to float
    typedef itk::CastImageFilter< kvl::CroppedImageReader::ImageType, ImageType > CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( reader->GetImage() );
    caster->Update();

    // Store the image and transform in persistent memory
    imageHandle = caster->GetOutput();
    transform = TransformType::New();
    reader->GetWorldToImageTransform()->GetInverse( transform );
    std::cout << "Read image: " << imageFileName << std::endl;
}

py::array_t<double> KvlImage::getTransformMatrix() {
    double data[16];
    auto parameters = transform->GetParameters();
    for ( unsigned int row = 0; row < 3; row++ )
    {
        for ( unsigned int col = 0; col < 3; col++ )
        {
            data[ col * 4 + row ] = parameters[ row * 3 + col ];
        }
        data[ 12 + row ] = parameters[ 9 + row ];
    }
    for ( unsigned int col = 0; col < 3; col++ )
    {
        data[ col * 4 + 3 ] = 0.0f;
    }
    data[ 15 ] = 1.0f;
    auto result = py::array_t<double>(16, data);
    result.resize({4, 4});
    return result;
}


void KvlImage::greet() {
    std::cout << "hello from KvlImage" << std::endl;
}

void useImage(KvlImage* image) {
    image->greet();
}
