#pragma once

#include "mri.h"
#include "log.h"

#include "bindings_numpy.h"

// conversion between MRI cxx objects and surfa FramedArray python objects
py::object MRItoSurfaArray(MRI* mri, bool release);
MRI* MRIfromSurfaArray(py::object arr);

// wrapped functions
py::object readMRI(const std::string& filename);
void writeMRI(py::object vol, const std::string& filename);
