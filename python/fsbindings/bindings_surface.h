#pragma once

#include "mrisurf.h"
#include "log.h"

#include "bindings_numpy.h"
#include "bindings_mri.h"

using namespace pybind11::literals;


// conversion between MRIS cxx objects and surfa Mesh python objects
py::object MRIStoSurfaMesh(MRIS *mris, bool release);
MRIS* MRISfromSurfaMesh(py::object surface);

// wrapped functions
py::object readSurface(const std::string& filename);
void writeSurface(py::object surf, const std::string& filename);
py::object computeTangents(py::object surf);
int computeEulerNumber(py::object surf);
int countIntersections(py::object surf);
py::object smoothOverlay(py::object surf, py::object overlay, int steps);
py::object quickSphericalInflate(py::object surf, int max_passes, int n_averages, long seed);
