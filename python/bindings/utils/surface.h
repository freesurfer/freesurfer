#pragma once

#include "numpy.h"
#include "mrisurf.h"
#include "log.h"
#include "volume.h"

namespace surf {


// Bridge allows for easy conversion between MRIS objects
// and python Surface objects when writing c-bindings.
class Bridge
{
public:
  // MRIS instance conversions
  Bridge(MRIS* mris) { setmris(mris); }
  operator MRIS*() { return mris(); }
  MRIS* mris();

  // py::object Surface conversions
  Bridge(py::object src) : source(src) {}
  operator py::object() { return python(); }
  py::object source = py::none();
  void updateSource();

private:
  void setmris(MRIS* m) { p_mris = std::shared_ptr<MRIS>(m, [](MRIS* ptr) { MRISfree(&ptr); }); }
  py::object python();
  void transferParameters(py::object& pyobj);

  std::shared_ptr<MRIS> p_mris;
};

// IO
py::object read(const std::string& filename);
// py::object read_directly(const std::string& filename);
void write(Bridge surf, const std::string& filename);

// metrics
void computeNormals(Bridge surf);
void computeTangents(Bridge surf);
int computeEulerNumber(Bridge surf);
int countIntersections(Bridge surf);

py::object surfaceDistance(Bridge surf1, Bridge surf2);

// overlay utils
py::array parameterize(Bridge surf, const arrayf<float>& overlay, int scale, std::string interp);
py::array sampleParameterization(Bridge surf, const arrayf<float>& image, std::string interp);
py::object smoothOverlay(Bridge surf, vol::Bridge overlay, int steps);

// surface submodule binding
inline void bind(py::module &m)
{
  py::class_<Bridge>(m, "Bridge").def(py::init<py::object>());
  py::implicitly_convertible<py::object, Bridge>();

  m.def("read", &read);
  m.def("write", &write);
  m.def("compute_normals", &computeNormals);
  m.def("compute_tangents", &computeTangents);
  m.def("compute_euler", &computeEulerNumber);
  m.def("count_intersections", &countIntersections);
  m.def("parameterize", &parameterize);
  m.def("sample_parameterization", &sampleParameterization);
  m.def("smooth_overlay", &smoothOverlay);
  m.def("distance", &surfaceDistance);
  // m.def("read_directly", &read_directly);
}


}  // end namespace surf
