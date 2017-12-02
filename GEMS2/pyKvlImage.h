#ifndef GEMS_KVLIMAGE_H
#define GEMS_KVLIMAGE_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "itkImage.h"
#include "pyKvlTransform.h"

namespace py = pybind11;

typedef itk::Image< float, 3 >  ImageType;
typedef ImageType::Pointer ImagePointer;

class KvlImage {
public:
    ImagePointer m_image;
    TransformPointer m_transform;
    KvlImage(const std::string &imageFileName);
    KvlImage(const py::array_t<float> &buffer);
    std::unique_ptr<KvlTransform> GetTransform();
    py::array_t<float> GetImageBuffer();
    void Write(std::string, KvlTransform &);
};



#endif //GEMS_KVLIMAGE_H
