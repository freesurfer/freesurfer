#include <pybind11/pybind11.h>

#include "bindings_mri.h"
#include "bindings_surface.h"
#include "bindings_transform.h"


PYBIND11_MODULE(fsbindings, m) {
  throwExceptions(true);

  // mri function bindings
  m.def("read_mri", &readMRI);
  m.def("write_mri", &writeMRI);

  // surface function bindings
  m.def("read_surf", &readSurface);
  m.def("write_surf", &writeSurface);
  m.def("compute_tangents", &computeTangents);
  m.def("compute_euler", &computeEulerNumber);
  m.def("count_intersections", &countIntersections);
  m.def("smooth_overlay", &smoothOverlay);
  m.def("quick_spherical_inflate", &quickSphericalInflate);

  // transform function bindings
  m.def("read_lta", &readLTA);
  m.def("write_lta", &writeLTA);
  m.def("write_gca_morph", &writeGCAMORPHfromPython);
}
