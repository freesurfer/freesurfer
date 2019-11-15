#ifndef VOLUME_H
#define VOLUME_H

#include "pybind11/stl.h"

#include "numpy.h"
#include "mri.h"
#include "mri2.h"
#include "log.h"


typedef py::array_t<float, py::array::f_style | py::array::forcecast> affinematrix;


/**
  MRI subclass to allow interaction between c++ and python.
*/
class PyVolume
{
public:
  PyVolume(py::array& array);
  PyVolume(py::array& array, affinematrix& affine);
  PyVolume(const std::string& filename) { m_mri = new MRI(filename); }
  ~PyVolume();

  PyVolume(MRI* mri) { m_mri = mri; }
  operator MRI*() { return m_mri; }

  void write(const std::string& filename) { m_mri->write(filename); }

  static py::array copyImage(MRI* mri);
  py::array copyImage() { return copyImage(m_mri); };

  void setImage(const py::array& array);
  py::array getImage();

  void setAffine(const affinematrix& array);
  affinematrix getAffine();
  affinematrix computeVox2Surf();

  template <class T>
  void setBufferData(py::array_t<T, py::array::f_style | py::array::forcecast> array) {
    MRI::Shape inshape = MRI::Shape(array.request().shape);
    if (inshape != m_mri->shape) logFatal(1) << "array " << shapeString(inshape) << " does not match volume shape " << shapeString(m_mri->shape);
    T *dst = (T *)m_mri->chunk;
    const T *src = array.data(0);
    for (unsigned int i = 0; i < m_mri->vox_total ; i++, dst++, src++) *dst = *src;
  }

  void setEchoTime(float te) { m_mri->te = te; }
  float getEchoTime() { return m_mri->te; }

  void setRecoveryTime(float tr)  { m_mri->tr = tr; }
  float getRecoveryTime() { return m_mri->tr; }

  void setInversionTime(float ti)  { m_mri->ti = ti; }
  float getInversionTime() { return m_mri->ti; }

  void setFlipAngle(double flip) { m_mri->flip_angle = flip; }
  double getFlipAngle() { return m_mri->flip_angle; }

private:
  py::array buffer_array;
  MRI* m_mri;
};

inline void bindVolume(py::module &m)
{
  // PyVolume class
  py::class_<PyVolume>(m, "Volume")
    .def(py::init<py::array&>())
    .def(py::init<const std::string&>())
    .def("write", &PyVolume::write)
    .def("_compute_vox2surf", &PyVolume::computeVox2Surf)
    .def_property("image", &PyVolume::getImage, &PyVolume::setImage)
    .def_property("affine", &PyVolume::getAffine, &PyVolume::setAffine)
    .def_property("te", &PyVolume::getEchoTime, &PyVolume::setEchoTime)
    .def_property("tr", &PyVolume::getRecoveryTime, &PyVolume::setRecoveryTime)
    .def_property("ti",  &PyVolume::getInversionTime, &PyVolume::setInversionTime)
    .def_property("flip_angle", &PyVolume::getFlipAngle, &PyVolume::setFlipAngle)
  ;
}

#endif
