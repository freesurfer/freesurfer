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
    std::vector<int> m_non_cropped_image_size;
    py::tuple m_crop_slices;

    // Implements kvlReadImage
    KvlImage(const std::string &imageFileName);
    // Implements kvlReadCroppedImage
    KvlImage(const std::string &imageFileName, const std::string &boundingFileName);
    KvlImage(const py::array_t<float> &buffer);
    std::unique_ptr<KvlTransform> GetTransform();
    py::array_t<float> GetImageBuffer();
    std::vector<int> GetNonCroppedImageSize();
    py::tuple GetCropSlices();
    void Write(std::string, KvlTransform &);
    void WriteImage(std::string);
    static py::array_t<float> smoothImageBuffer(const py::array_t<float>& imageBuffer, std::vector<double> sigmas);

private:
    static ImagePointer numpy_to_image(const py::array_t<float> &buffer);
    static py::array_t<float> image_to_numpy(ImagePointer image);
};



#endif //GEMS_KVLIMAGE_H
