#ifndef SURFACE_H
#define SURFACE_H

#include "numpy.h"
#include "mrisurf.h"
#include "volume.h"


typedef py::array_t<float, py::array::f_style | py::array::forcecast> affinematrix;


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

  // geometry
  void copyGeometry(PySurface *surf) { copyVolGeom(&surf->vg, &this->vg); }
  affinematrix computeSurf2Vox(PyVolume& vol);

  // wrapped utilities
  bool isSelfIntersecting();
  py::array fillInterior();

  // parameterization
  py::array_t<float> parameterizeBarycentric(const py::array_t<float, py::array::forcecast>& array, float scale);
  py::array_t<float> computeParameterizationMapBarycentric(const std::vector<int>& shape);
};

py::array_t<int> readAnnotation(const std::string &filename);

inline void bindSurface(py::module &m)
{
  // PySurface class
  py::class_<PySurface>(m, "Surface")
    .def(py::init<const std::string &>())
    .def("write", &PySurface::write)
    .def("_compute_surf2vox", &PySurface::computeSurf2Vox)
    .def("isSelfIntersecting", &PySurface::isSelfIntersecting)
    .def("fillInterior", &PySurface::fillInterior)
    .def("copy_geometry", (void (PySurface::*)(PySurface*)) &PySurface::copyGeometry)
    .def("parameterize", &PySurface::parameterizeBarycentric)
    .def("precompute_parameterization_map", &PySurface::computeParameterizationMapBarycentric)
    .def_property("vertices", &PySurface::getVertices, &PySurface::setVertices)
    .def_property("faces", &PySurface::getFaces, &PySurface::setFaces)
    ;

  m.def("read_annotation", &readAnnotation);
}

#endif