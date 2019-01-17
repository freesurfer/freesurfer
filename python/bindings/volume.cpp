#include <stdexcept>

#include "volume.h"


/*
  Volume-specific pybind module configuration.
*/
void bindVolume(py::module &m) {
  // CoreVolume class
  py::class_<CoreVolume>(m, "CoreVolume")
    .def(py::init<const std::string &>())
    .def(py::init<py::array>())
    .def("read", &CoreVolume::read)
    .def("write", &CoreVolume::write)
    .def_property("image", &CoreVolume::getImage, &CoreVolume::setImage)
  ;
}


/*
  Creates a numpy array handle from an MRI structure. If `copybuffer` is `false`, the
  returned handle will share the underlying buffer data of the MRI - otherwise, a new
  array buffer is allocated and the data is copied.
*/
py::array makeArray(MRI *mri, bool copybuffer) {
  // make sure mri buffer is chunked
  if (MRIchunk(&mri) != 0) throw py::value_error("could not chunk buffer data");

  // determine the appropriate numpy dtype from the mri type
  py::dtype dtype;
  switch (mri->type) {
    case MRI_UCHAR:
      dtype = py::dtype::of<unsigned char>(); break;
    case MRI_SHORT:
      dtype = py::dtype::of<short>(); break;
    case MRI_INT:
      dtype = py::dtype::of<int>(); break;
    case MRI_LONG:
      dtype = py::dtype::of<long>(); break;
    case MRI_FLOAT:
      dtype = py::dtype::of<float>(); break;
  }

  // compute strides (mri is always column-major)
  std::vector<ssize_t> shape = {mri->width, mri->height, mri->depth, mri->nframes};
  std::vector<ssize_t> strides = fstrides(shape, mri->bytes_per_vox);

  // make the numpy array
  py::array array;
  if (copybuffer) {
    // if no handle is provided, pybind will copy the buffer data
    array = py::array(dtype, shape, strides, mri->chunk);
  } else {
    // if an existing handle (simple capsule) is passed, the buffer data will be shared with the mri instance
    array =py::array(dtype, shape, strides, mri->chunk, py::capsule(mri->chunk));
  }

  return array;
}


/*
  Constructs a volume from an already allocated MRI structure. This does not copy
  the MRI instance, it only applies ownership of the pointer. Therefore, be careful
  not to free the MRI separately.
*/
CoreVolume::CoreVolume(MRI *mri) {
  setMRI(mri);
}


/*
  Constructs a volume from a volume file.
*/
CoreVolume::CoreVolume(const std::string &filename) {
  read(filename);
}


/*
  Constructs a volume from a 3D or 4D numpy array. The array dtype is used to
  determine the type of the underlying MRI structure.
*/
CoreVolume::CoreVolume(py::array array) {
  // determine input shape
  py::buffer_info info = array.request();
  int nframes = 1;
  if (info.ndim == 4) {
    nframes = info.shape[3];
  } else if (info.ndim != 3) {
    throw py::value_error("volume must be constructed with a 3D or 4D array");
  }
  int width = info.shape[0];
  int height = info.shape[1];
  int depth = info.shape[2];

  // determine mri type from numpy dtype
  int type;
  if (py::isinstance<py::array_t<unsigned char>>(array)) {
    type = MRI_UCHAR;
  } else if ((py::isinstance<py::array_t<int>>(array)) || (py::isinstance<py::array_t<bool>>(array))) {
    type = MRI_INT;
  } else if (py::isinstance<py::array_t<long>>(array)) {
    type = MRI_LONG;
  } else if ((py::isinstance<py::array_t<float>>(array)) || (py::isinstance<py::array_t<double>>(array))) {
    type = MRI_FLOAT;
  } else if (py::isinstance<py::array_t<short>>(array)) {
    type = MRI_SHORT;
  } else {
    throw py::value_error("unknown array dtype " + array.dtype().kind());
  }

  // allocate and set the mri
  setMRI(MRIallocChunk(width, height, depth, type, nframes));

  // set the underlying mri buffer data to match the input
  setImage(array);
}


/*
  Frees the underlying MRI before destructing.
*/
CoreVolume::~CoreVolume() {
  if (m_mri) MRIfree(&m_mri);
}


/*
  Sets the underlying MRI instance.
*/
void CoreVolume::setMRI(MRI *mri) {
  // make sure mri is chunked so we can access contiguous buffer data
  if (MRIchunk(&mri) != 0) throw py::value_error("could not chunk buffer data");
  // update the underlying mri pointer
  if (m_mri) MRIfree(&m_mri);
  m_mri = mri;
  // update the numpy array buffer
  imagebuffer = makeArray(m_mri, false);
}


/*
  Loads the underlying MRI structure from a volume file.
*/
void CoreVolume::read(const std::string &filename) {
  MRI *mri = MRIread(filename.c_str());
  if (!mri) throw py::value_error("could not load volume " + filename);
  setMRI(mri);
}


/*
  Writes the underlying MRI structure to a file.
*/
void CoreVolume::write(const std::string &filename) {
  MRIwrite(m_mri, filename.c_str());
}


/*
  Sets the MRI image buffer from a numpy array. The input array must match the shape of the
  underlying MRI, but if the MRI has only 1 frame, then a 3D input is also allowed.
*/
void CoreVolume::setImage(py::array_t<float, py::array::f_style | py::array::forcecast> array) {
  // get buffer info
  py::buffer_info info = array.request();
  
  // sanity check on input dimensions
  auto mrishape = imagebuffer.request().shape;
  py::value_error shapeError("array must be of shape " + shapeString(mrishape));
  
  // allow 3D arrays as long as the mri only has 1 frame 
  if (info.ndim == 3) {
    if (m_mri->nframes > 1) throw shapeError;
    mrishape.pop_back();
  }
  if (info.shape != mrishape) throw shapeError;

  // set type-specific image buffer data
  switch (m_mri->type) {
    case MRI_UCHAR:
      setBufferData<unsigned char>(array.data(0)); break;
    case MRI_SHORT:
      setBufferData<short>(array.data(0)); break;
    case MRI_INT:
      setBufferData<int>(array.data(0)); break;
    case MRI_LONG:
      setBufferData<long>(array.data(0)); break;
    case MRI_FLOAT:
      setBufferData<float>(array.data(0)); break;
  }
}
