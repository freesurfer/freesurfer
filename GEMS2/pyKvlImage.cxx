#include "pyKvlImage.h"
#include "itkCastImageFilter.h"
#include "pyKvlNumpy.h"
#include "pyKvlTransform.h"
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

std::unique_ptr<KvlTransform> KvlImage::GetTransformMatrix() {
    return std::unique_ptr<KvlTransform>(new KvlTransform(transform));
}

py::array_t<float> KvlImage::GetImageBuffer() {
    auto region = imageHandle->GetLargestPossibleRegion();
    auto shape = region.GetSize();
    auto* const buffer = new float[region.GetNumberOfPixels()];

    itk::ImageRegionConstIterator< ImageType > it( imageHandle, region );
    float* buffer_p = buffer;
    for ( ;!it.IsAtEnd(); ++it, ++buffer_p ){
        *buffer_p = it.Value();
    }
    return createNumpyArrayFStyle({shape[0], shape[1], shape[2]}, buffer);
}
