#include <pybind11/pybind11.h>

#include "volume.h"
#include "surface.h"


PYBIND11_MODULE(bindings, m) {
  throwExceptions(true);

  py::module mvol = m.def_submodule("vol");
  vol::bind(mvol);
}
