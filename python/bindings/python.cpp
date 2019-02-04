#include <pybind11/pybind11.h>

#include "volume.h"
#include "surface.h"


PYBIND11_MODULE(bindings, m) {
  bindVolume(m);
  bindSurface(m);
}
