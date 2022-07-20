#pragma once

#include "transform.h"
#include "log.h"

#include "bindings_numpy.h"


void VOLGEOMfromSurfaImageGeometry(py::object geometry, VOL_GEOM* vg);
py::object VOLGEOMtoSurfaImageGeometry(VOL_GEOM* vg);
LTA* LTAfromSurfaAffine(py::object affine);
py::object LTAtoSurfaAffine(const LTA* lta);
py::object readLTA(const std::string& filename);
void writeLTA(py::object transform, const std::string& filename);
void writeGCAMORPHfromPython(py::object warp, py::object affine, py::object source,
                             py::object target, const std::string& filename);
