#pragma once

#include "numpy.h"
#include "mrisurf.h"
#include "log.h"


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

private:
  void setmri(MRIS* m) { p_mris = std::shared_ptr<MRIS>(m, [](MRIS* ptr) { MRIfree(&ptr); }); }
  py::object python();

  std::shared_ptr<MRIS> p_mris;
  py::object source = py::none();
};

// IO
py::object read(const std::string& filename);
py::object read_directly(const std::string& filename);
void write(Bridge surf, const std::string& filename);

// metrics
void computeNormals(py::object surf);


// surface submodule binding
inline void bind(py::module &m)
{
  py::class_<Bridge>(m, "Bridge").def(py::init<py::object>());
  py::implicitly_convertible<py::object, Bridge>();

  m.def("read", &read);
  m.def("write", &write);
  m.def("compute_normals", &computeNormals);
  m.def("read_directly", &read_directly);
}


}  // end namespace surf
