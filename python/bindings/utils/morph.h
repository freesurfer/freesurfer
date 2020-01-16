#pragma once

#include "numpy.h"

namespace morph {


void writeGCAM(py::object warp, const std::string& filename);

// transform submodule binding
inline void bind(py::module &m)
{
  m.def("write_gcam", &writeGCAM);
}


}  // end namespace morph
