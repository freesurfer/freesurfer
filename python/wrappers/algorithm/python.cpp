#include "python.h"
#include "numpy.h"


PYBIND11_MODULE(algorithm_python, m)
{
  m.def("maxflow", &CMF3D_ML);
  m.def("performFrontPropagation3D", &performFrontPropagation3D);
}
