#include <stdexcept>

#include "volume.h"


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
  Constructs an MRI instance from a 3/4D numpy array.
*/
PyVolume::PyVolume(py::array& array)
{
  m_mri = new MRI(array.request().shape, voltype(array));
  setImage(array);
};


/**
  Constructs an MRI instance from a 3/4D numpy array and an affine matrix.
*/
PyVolume::PyVolume(py::array& array, affinematrix& affine) : PyVolume(array)
{
  setAffine(affine);
};


PyVolume::~PyVolume()
{
  // if python is handling the buffer memory, null the chunk pointer so that
  // the MRI destructor doesn't actually delete the image data
  if (buffer_array.size() != 0) m_mri->chunk = nullptr;
  delete m_mri;
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
    if (!m_mri->ischunked) logFatal(1) << "image is too large to fit into contiguous memory and cannot be supported by the python bindings";
    // python will now manage the buffer memory deletion
    py::capsule capsule(m_mri->chunk, [](void *d) { if (d) free(d); } );
    buffer_array = py::array(pytype(m_mri), std::vector<ssize_t>(m_mri->shape), fstrides(m_mri->shape, m_mri->bytes_per_vox), m_mri->chunk, capsule);
  }
  return buffer_array;
}


/**
  Sets the MRI image buffer from a numpy array. The input array must match the shape of the
  underlying data, but if the MRI has only 1 frame, then a 3D input is also allowed.
*/
void PyVolume::setImage(const py::array& array)
{
  switch (m_mri->type) {
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


/**
  Gets the volume's affine vox->ras matrix.
*/
affinematrix PyVolume::getAffine()
{
  return affinematrix({4, 4}, fstrides({4, 4}, sizeof(float)), m_mri->i_to_r__->mat);
}


/**
  Sets the volume's affine vox->ras matrix.
*/
void PyVolume::setAffine(const affinematrix& array)
{
  const float *ptr = array.data(0);
  double xr = ptr[0];
  double xa = ptr[1];
  double xs = ptr[2];

  double yr = ptr[4];
  double ya = ptr[5];
  double ys = ptr[6];

  double zr = ptr[8];
  double za = ptr[9];
  double zs = ptr[10];

  double pr = ptr[12];
  double pa = ptr[13];
  double ps = ptr[14];

  double sizex = std::sqrt(xr * xr + xa * xa + xs * xs);
  double sizey = std::sqrt(yr * yr + ya * ya + ys * ys);
  double sizez = std::sqrt(zr * zr + za * za + zs * zs);

  m_mri->x_r = xr / sizex;
  m_mri->x_a = xa / sizex;
  m_mri->x_s = xs / sizex;

  m_mri->y_r = yr / sizey;
  m_mri->y_a = ya / sizey;
  m_mri->y_s = ys / sizey;

  m_mri->z_r = zr / sizez;
  m_mri->z_a = za / sizez;
  m_mri->z_s = zs / sizez;

  MRIp0ToCRAS(m_mri, pr, pa, ps);
}


/**
  Computes the volume's vox->surf matrix.
*/
affinematrix PyVolume::computeVox2Surf()
{
  MATRIX *vox2surf = MRIxfmCRS2XYZtkreg(m_mri);
  affinematrix pymat = affinematrix({4, 4}, cstrides({4, 4}, sizeof(float)), vox2surf->data);
  MatrixFree(&vox2surf);
  return pymat;
}
