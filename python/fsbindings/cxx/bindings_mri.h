#pragma once

#include "mri.h"
#include "log.h"

#include "bindings_numpy.h"


// Bridge allows for easy conversion between MRI objects
// and python ArrayContainerTemplate objects when writing c-bindings.
class MRIBridge
{
public:
  // FS MRI instance conversions
  MRIBridge(MRI* mri) { setMRI(mri); }
  operator MRI*() { return toMRI(); }
  MRI* toMRI();

  // surfa FramedArray conversions
  MRIBridge(py::object src) : source(src) {}
  operator py::object() { return toPython(); }
  void updateSource();

private:
  void setMRI(MRI* m) { p_mri = std::shared_ptr<MRI>(m, [](MRI* ptr) { MRIfree(&ptr); }); }
  py::object toPython();
  void transferParameters(py::object& pyobj);

  std::shared_ptr<MRI> p_mri;
  py::array mri_buffer;
  py::object source = py::none();
};

py::object readMRI(const std::string& filename);
void writeMRI(MRIBridge vol, const std::string& filename);
