#pragma once

#include "mrisurf.h"
#include "log.h"

#include "bindings_numpy.h"
#include "bindings_mri.h"

using namespace pybind11::literals;


// Bridge allows for easy conversion between MRIS objects
// and python Surface objects when writing c-bindings.
class MRISBridge
{
public:
  // MRIS instance conversions
  MRISBridge(MRIS* mris) { setMRIS(mris); }
  operator MRIS*() { return toMRIS(); }
  MRIS* toMRIS();

  // python object Surface conversions
  MRISBridge(py::object mesh) { source = mesh.attr("convert")("space"_a="surface"); }
  operator py::object() { return toPython(); }
  py::object source = py::none();

private:
  void setMRIS(MRIS* m) { p_mris = std::shared_ptr<MRIS>(m, [](MRIS* ptr) { MRISfree(&ptr); }); }
  py::object toPython();

  std::shared_ptr<MRIS> p_mris;
};

py::object readSurface(const std::string& filename);
void writeSurface(MRISBridge surf, const std::string& filename);
py::object computeTangents(MRISBridge surf);
int computeEulerNumber(MRISBridge surf);
int countIntersections(MRISBridge surf);
py::object surfaceDistance(MRISBridge surf1, MRISBridge surf2);
py::object parameterize(MRISBridge surf, py::object overlay, int scale, std::string interp);
py::object sampleParameterization(MRISBridge surf, py::object image, std::string interp);
py::object smoothOverlay(MRISBridge surf, MRIBridge overlay, int steps);
py::object quickSphericalInflate(Bridge inSurf, int max_passes, int n_averages, long seed);
