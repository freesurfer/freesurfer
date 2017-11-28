#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
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

        py::array_t<double> getTransformMatrix() {
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

};

void useImage(KvlImage* image) {
    image->greet();
}



PYBIND11_MODULE(GEMS2Python, m) {
    py::class_<KvlImage>(m, "KvlImage")
    .def(py::init<const std::string &>())
    .def("greet", &KvlImage::greet)
    .def("getTransformMatrix", &KvlImage::getTransformMatrix);

    m.def("useImage", &useImage);
}
