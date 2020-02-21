#pragma once

#include "numpy.h"
#include "mri.h"
#include "log.h"

namespace vol {


// Bridge allows for easy conversion between MRI objects
// and python ArrayContainerTemplate objects when writing c-bindings.
class Bridge
{
public:
  // MRI instance conversions
  Bridge(MRI* mri) { setmri(mri); }
  operator MRI*() { return mri(); }
  MRI* mri();

  // py::object conversions
  Bridge(py::object src) : source(src) {}
  operator py::object() { return python(); }
  void updateSource();

private:
  void setmri(MRI* m) { p_mri = std::shared_ptr<MRI>(m, [](MRI* ptr) { MRIfree(&ptr); }); }
  py::object python();
  void transferParameters(py::object& pyobj);

  std::shared_ptr<MRI> p_mri;
  py::array mri_buffer;
  py::object source = py::none();
};


// IO
py::object read(const std::string& filename);
void write(Bridge vol, const std::string& filename);


// TODOC
using sample_array = py::array_t<float, py::array::f_style>;
void sampleIntoVolume(sample_array volume, sample_array weights, arrayc<float> coords, arrayc<float> values);
py::array_t<float> resampleVolume(py::array_t<float> source_vol, py::object target_shape, py::array_t<float> target2source);


// volume submodule binding
inline void bind(py::module &m)
{
  py::class_<Bridge>(m, "Bridge").def(py::init<py::object>());
  py::implicitly_convertible<py::object, Bridge>();

  m.def("read", &read);
  m.def("write", &write);
  m.def("sample_into_volume", &sampleIntoVolume);
  m.def("resample_volume", &resampleVolume);
}


}  // end namespace vol
