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


// resampling
void sampleIntoVolume(py::array volume, py::array weights, arrayc<double> coords, arrayc<double> values);
py::array resampleVolumeLinear(py::array source_vol, py::object target_shape, py::array_t<double> target2source, double fill = 0);
py::array resampleVolumeNearest(py::array source_vol, py::object target_shape, py::array_t<double> target2source, double fill = 0);


// volume submodule binding
inline void bind(py::module &m)
{
  py::class_<Bridge>(m, "Bridge").def(py::init<py::object>());
  py::implicitly_convertible<py::object, Bridge>();

  m.def("read", &read);
  m.def("write", &write);
  m.def("sample_into_volume", &sampleIntoVolume);
  m.def("resample_volume_linear", &resampleVolumeLinear);
  m.def("resample_volume_nearest", &resampleVolumeNearest);
}


}  // end namespace vol
