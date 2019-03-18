#ifndef SURFACE_H
#define SURFACE_H

#include "numpy.h"

// utils
#include "mrisurf.h"


void bindSurface(py::module &m);


class CoreSurface {
public:
  CoreSurface(const std::string &filename);
  ~CoreSurface();

  // mris pointer getter/setter
  void setMRIS(MRIS *mris);
  MRIS* getMRIS() { return m_mris; };

  // filesystem IO
  void read(const std::string &filename);
  void write(const std::string &filename);

  // vertices getter/setter
  py::array_t<float> getVertices();
  void setVertices(py::array_t<float, py::array::c_style | py::array::forcecast>);

  // wrapped utilities
  bool isSelfIntersecting();
  pybind11::array fillInterior();

private:
  MRIS *m_mris = nullptr;
};

#endif