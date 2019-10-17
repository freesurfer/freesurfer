#pragma once

#include "numpy.h"
#include "mri.h"
#include "log.h"


namespace vol {


// Bridge allows for easy conversion between MRI objects
// and python objects when writing c-bindings.
class Bridge
{
public:
  // MRI instance conversions
  Bridge(MRI* mri) { setmri(mri); }
  operator MRI*() { return mri(); }
  MRI* mri();

  // py::object MultidimArray conversions
  Bridge(py::object src) : source(src) {}
  operator py::object() { return python(); }

private:
  void setmri(MRI* m) { p_mri = std::shared_ptr<MRI>(m, [](MRI* ptr) { MRIfree(&ptr); }); }
  py::object python();

  std::shared_ptr<MRI> p_mri;
  py::array mri_buffer;
  py::object source = py::none();
};


// MultidimArray IO
py::object read(const std::string& filename);
void write(Bridge vol, const std::string& filename);


// volume submodule binding
inline void bind(py::module &m)
{
  py::class_<Bridge>(m, "Bridge").def(py::init<py::object>());
  py::implicitly_convertible<py::object, Bridge>();

  m.def("read", &read);
  m.def("write", &write);
}


}  // end namespace vol
