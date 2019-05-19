#ifndef VOLUME_H
#define VOLUME_H

#include "numpy.h"
#include "mri.h"
#include "log.h"

void bindVolume(py::module &m);

typedef py::array_t<float, py::array::f_style | py::array::forcecast> affinematrix;

/**
  MRI subclass to allow interaction between c++ and python.
*/
class PyVolume : public MRI
{
public:
  PyVolume(py::array& array);
  PyVolume(py::array& array, affinematrix& affine);
  PyVolume(const std::string& filename) : MRI(filename) {};
  ~PyVolume();

  static py::array copyImage(MRI *mri);
  py::array copyImage() { return copyImage(this); };

  void setImage(const py::array& array);
  py::array getImage();

  void setAffine(const affinematrix& array);
  affinematrix getAffine();

  template <class T>
  void setBufferData(py::array_t<T, py::array::f_style | py::array::forcecast> array) {
    MRI::Shape inshape = MRI::Shape(array.request().shape);
    if (inshape != shape) logFatal(1) << "array " << shapeString(inshape) << " does not match volume shape " << shapeString(shape);
    T *dst = (T *)chunk;
    const T *src = array.data(0);
    for (unsigned int i = 0; i < vox_total ; i++, dst++, src++) *dst = *src;
  }

private:
  py::array buffer_array;
};

#endif
