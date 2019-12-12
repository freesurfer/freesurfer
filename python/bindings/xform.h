#pragma once

#include "numpy.h"
#include "transform.h"

py::object volGeomToPython(VOL_GEOM* vg);
py::object readLTA(const std::string& filename);

inline void bindTransform(py::module &m)
{
  m.def("read_lta", &readLTA);
}
