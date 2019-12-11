#include <pybind11/pybind11.h>

#include "volume.h"
#include "surface.h"
#include "xform.h"


PYBIND11_MODULE(bindings, m) {
  bindVolume(m);
  bindSurface(m);
  bindTransform(m);
}
