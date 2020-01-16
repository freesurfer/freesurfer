#pragma once

#include "numpy.h"
#include "transform.h"
#include "log.h"

namespace transform {


LTA* pythonToLTA(py::object transform);

void pythonToVolGeom(py::object geometry, VOL_GEOM* vg);
py::object volGeomToPython(VOL_GEOM* vg);

py::object readLTA(const std::string& filename);
void writeLTA(py::object transform, const std::string& filename);


// transform submodule binding
inline void bind(py::module &m)
{
  m.def("read_lta", &readLTA);
  m.def("write_lta", &writeLTA);
}


}  // end namespace transform
