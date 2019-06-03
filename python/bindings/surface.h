#ifndef SURFACE_H
#define SURFACE_H

#include "numpy.h"
#include "mrisurf.h"


class PySurface : public MRIS
{
public:
  PySurface(const std::string &filename);

  // filesystem IO
  void write(const std::string &filename);

  // vertices getter/setter
  py::array_t<float> getVertices();
  void setVertices(py::array_t<float, py::array::c_style | py::array::forcecast> array);

  // faces getter/setter
  py::array_t<int> getFaces();
  void setFaces(py::array_t<int, py::array::c_style | py::array::forcecast> array);

  // wrapped utilities
  bool isSelfIntersecting();
  pybind11::array fillInterior();
};

py::array_t<int> readAnnotation(const std::string& filename);

inline void bindSurface(py::module &m)
{
  // PySurface class
  py::class_<PySurface>(m, "Surface")
    .def(py::init<const std::string &>())
    .def("write", &PySurface::write)
    .def("isSelfIntersecting", &PySurface::isSelfIntersecting)
    .def("fillInterior", &PySurface::fillInterior)
    .def_property("vertices", &PySurface::getVertices, &PySurface::setVertices)
    .def_property("faces", &PySurface::getFaces, &PySurface::setFaces)
    ;

  m.def("read_annotation", &readAnnotation);
}

#endif