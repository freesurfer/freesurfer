#include "pyKvlImage.h"
#include "itkCastImageFilter.h"

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

void KvlImage::greet() {
    std::cout << "hello from KvlImage" << std::endl;
}

void useImage(KvlImage* image) {
    image->greet();
}
