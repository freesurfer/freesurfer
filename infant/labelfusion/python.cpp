#include <pybind11/pybind11.h>

#include "labelfusion.h"


PYBIND11_MODULE(labelfusion, m) {
  m.def("maxflow", &CMF3D_ML);
  m.def("performFrontPropagation3D", &performFrontPropagation3D);
}
