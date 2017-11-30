#include "pyKvlImage.h"
#include "itkCastImageFilter.h"
#include "pyKvlNumpy.h"
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
    return transformToNumpy(transform);
}

py::array_t<float> KvlImage::getImageBuffer() {
    auto shape = imageHandle->GetBufferedRegion().GetSize();
    imageHandle->GetPixelContainer()->SetContainerManageMemory(false);
    return createNumpyArray({shape[2], shape[1], shape[0]}, imageHandle->GetPixelContainer()->GetImportPointer());
}
