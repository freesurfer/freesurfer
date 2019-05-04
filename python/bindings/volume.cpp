#include <stdexcept>

#include "volume.h"


/**
  Binds the PyVolume class (renamed to Volume) in the python module.
*/
void bindVolume(py::module &m)
{
  py::class_<PyVolume>(m, "Volume")
    .def(py::init<py::array&>())
    .def(py::init<const std::string&>())
    .def("write", &PyVolume::write)
    .def_property("image", &PyVolume::getImage, &PyVolume::setImage)
  ;
}


/**
  Returns the associated MRI type of a numpy array.
*/
static int voltype(const py::array& array)
{
  if (py::isinstance<py::array_t<unsigned char>>(array)) {
    return MRI_UCHAR;
  } else if ((py::isinstance<py::array_t<int>>(array)) || (py::isinstance<py::array_t<bool>>(array))) {
    return MRI_INT;
  } else if (py::isinstance<py::array_t<long>>(array)) {
    return MRI_LONG;
  } else if ((py::isinstance<py::array_t<float>>(array)) || (py::isinstance<py::array_t<double>>(array))) {
    return MRI_FLOAT;
  } else if (py::isinstance<py::array_t<short>>(array)) {
    return MRI_SHORT;
  } else {
    throw py::value_error("unknown numpy array dtype");
  }
}


/**
  Returns the associated numpy dtype of an MRI.
*/
static py::dtype pytype(const MRI* mri)
{
  switch (mri->type) {
  case MRI_UCHAR:
    return py::dtype::of<unsigned char>(); break;
  case MRI_SHORT:
    return py::dtype::of<short>(); break;
  case MRI_INT:
    return py::dtype::of<int>(); break;
  case MRI_LONG:
    return py::dtype::of<long>(); break;
  case MRI_FLOAT:
    return py::dtype::of<float>(); break;
  default:
    throw py::value_error("unknown MRI data type");
  }
}


/**
  Constructs an MRI instance from a 3D or 4D numpy array.
*/
PyVolume::PyVolume(py::array& array) : MRI(array.request().shape, voltype(array))
{
  setImage(array);
}


PyVolume::~PyVolume()
{
  // if python is handling the buffer memory, null the chunk pointer so that
  // the MRI destructor doesn't actually delete the image data
  if (buffer_array.size() != 0) chunk = nullptr;
}


/**
  Copies the MRI image buffer into a new numpy array.
*/
py::array PyVolume::copyImage(MRI *mri)
{
  if (!mri->ischunked) logFatal(1) << "image is too large to fit into contiguous memory and cannot be supported by the python bindings";
  return py::array(pytype(mri), std::vector<ssize_t>(mri->shape), fstrides(mri->shape, mri->bytes_per_vox), mri->chunk);
}


/**
  Shares the underlying MRI image buffer as a numpy array.
*/
py::array PyVolume::getImage()
{
  // create the python buffer array if it hasn't been initialized already
  if (buffer_array.size() == 0) {
    if (!ischunked) logFatal(1) << "image is too large to fit into contiguous memory and cannot be supported by the python bindings";
    // python will now manage the buffer memory deletion
    py::capsule capsule(chunk, [](void *d) { if (d) free(d); } );
    buffer_array = py::array(pytype(this), std::vector<ssize_t>(shape), fstrides(shape, bytes_per_vox), chunk, capsule);
  }
  return buffer_array;
}


/**
  Sets the MRI image buffer from a numpy array. The input array must match the shape of the
  underlying data, but if the MRI has only 1 frame, then a 3D input is also allowed.
*/
void PyVolume::setImage(const py::array& array) {
  switch (type) {
    case MRI_UCHAR:
      setBufferData<unsigned char>(array); break;
    case MRI_INT:
      setBufferData<int>(array); break;
    case MRI_SHORT:
      setBufferData<short>(array); break;
    case MRI_LONG:
      setBufferData<long>(array); break;
    case MRI_FLOAT:
      setBufferData<float>(array); break;
    default:
      throw py::value_error("unknown MRI data type");
  }
}
