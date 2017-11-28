#include <pybind11/pybind11.h>
#include "kvlCroppedImageReader.h"
#include "itkCastImageFilter.h"

namespace py = pybind11;

typedef itk::Image< float, 3 >  ImageType;
typedef ImageType::Pointer ImagePointer;
typedef kvl::CroppedImageReader::TransformType  TransformType;
typedef TransformType::Pointer TransformPointer;

class KvlImage {
    ImagePointer imageHandle;
    TransformPointer transform;

    public:
        KvlImage(const std::string &imageFileName) {
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
        void greet() {
            std::cout << "hello from KvlImage" << std::endl;
        }
};

void useImage(KvlImage* image) {
    image->greet();
}

PYBIND11_MODULE(GEMS2Python, m) {
    py::class_<KvlImage>(m, "KvlImage")
    .def(py::init<const std::string &>())
    .def("greet", &KvlImage::greet);

    m.def("useImage", &useImage);
}
