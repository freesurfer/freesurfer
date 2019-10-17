#include "volume.h"
// #include <pybind11/stl.h>  // TODO MOVE THIS ELSEWHERE??

namespace vol {


typedef py::array_t<double, py::array::forcecast | py::array::c_style> affinematrix;


/*
  Creates a python MultidimArray object from an MRI instance. This really shouldn't
  be called more than once since that would create multiple objects sharing the same data array.
*/
py::object Bridge::python()
{
  // sanity checks
  if (!p_mri) throw std::runtime_error("MRI is null");
  if (!p_mri->ischunked) throw std::runtime_error("image is too large to fit into contiguous memory and cannot be supported by the python bindings");

  // if it exists, the mri_buffer numpy array must always match the MRI instance data
  if (mri_buffer.size() == 0) {
    // determine numpy dtype
    py::dtype dtype;
    switch (p_mri->type) {
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
    default:
      throw py::value_error("unknown MRI data type");
    }

    // configuration (make sure shape gets squeezed here)
    std::vector<ssize_t> shape = {p_mri->width};
    if (p_mri->height  > 1) shape.push_back(p_mri->height);
    if (p_mri->depth   > 1) shape.push_back(p_mri->depth);
    if (p_mri->nframes > 1) shape.push_back(p_mri->nframes);
    std::vector<ssize_t> strides = fstrides(shape, p_mri->bytes_per_vox);

    // create numpy array from chunked mri buffer
    py::capsule capsule(p_mri->chunk, [](void *p) { free(p); } );
    mri_buffer = py::array(dtype, shape, strides, p_mri->chunk, capsule);

    // remove data ownership from MRI so that MRIfree doesn't delete the array upon destruction
    p_mri->owndata = false;
  }

  // extract base dimensions
  int ndims = (p_mri->nframes > 1) ? mri_buffer.ndim() - 1 : mri_buffer.ndim();

  // build python instance
  py::module fsmodule = py::module::import("freesurfer");
  py::object vol;
  switch (ndims) {
  case 1: vol = fsmodule.attr("Overlay")(mri_buffer); break;
  case 2: vol = fsmodule.attr("Image")(mri_buffer); break;
  case 3: vol = fsmodule.attr("Volume")(mri_buffer); break;
  }
  
  // set filename
  vol.attr("filename") = py::str(p_mri->fname);

  // volume-specific parameters
  if (ndims == 3) {
    MATRIX *matrix = extract_i_to_r(p_mri.get());
    affinematrix affine = makeArray({4, 4}, MemoryOrder::Fortran, matrix->data, false);
    MatrixFree(&matrix);
    vol.attr("affine") = affine;
    vol.attr("voxsize") = py::make_tuple(p_mri->xsize, p_mri->ysize, p_mri->zsize);
  }

  return vol;
}


/*
  Returns the associated MRI instance. If one does not exist, it will be created
  from the cached python input MultidimArray.
*/
MRI* Bridge::mri()
{
  // return if the MRI instance has already been created
  if (p_mri) return p_mri.get();

  // make sure the source has been provided
  if (source.is(py::none())) throw py::value_error("cannot generate MRI instance without source object");

  // ensure buffer is in fortran order
  py::module np = py::module::import("numpy");
  mri_buffer = np.attr("asfortranarray")(source.attr("array"));

  // convert unsupported types
  if (py::isinstance<py::array_t<bool>>(mri_buffer)) mri_buffer = py::array_t<uchar>(mri_buffer);
  if (py::isinstance<py::array_t<double>>(mri_buffer)) mri_buffer = py::array_t<float>(mri_buffer);

  // extract MRI type from numpy datatype
  int dtype;
  if      (py::isinstance<py::array_t<uchar>>(mri_buffer)) { dtype = MRI_UCHAR; }
  else if (py::isinstance<py::array_t<int>>  (mri_buffer)) { dtype = MRI_INT; }
  else if (py::isinstance<py::array_t<short>>(mri_buffer)) { dtype = MRI_SHORT; }
  else if (py::isinstance<py::array_t<long>> (mri_buffer)) { dtype = MRI_LONG; }
  else if (py::isinstance<py::array_t<float>>(mri_buffer)) { dtype = MRI_FLOAT; }
  else {
    throw py::value_error("unsupported array dtype " + py::str(mri_buffer.attr("dtype")).cast<std::string>());
  }

  // initialize a header-only MRI structure with the known shape
  std::vector<int> expanded, shape = mri_buffer.attr("shape").cast<std::vector<int>>();
  int nframes = mri_buffer.attr("nframes").cast<int>();
  int ndims = mri_buffer.attr("basedims").cast<int>();
  switch (ndims) {
  case 1: expanded = {shape[0], 1,        1,        nframes}; break;
  case 2: expanded = {shape[0], shape[1], 1,        nframes}; break;
  case 3: expanded = {shape[0], shape[1], shape[2], nframes}; break;
  }
  MRI *mri = new MRI(expanded, dtype, false);

  // point MRI chunk at the numpy array
  mri->chunk = mri_buffer.mutable_data();
  mri->ischunked = true;
  mri->owndata = false;
  mri->initSlices();

  // filename metadata
  std::strcpy(mri->fname, source.attr("filename").cast<std::string>().c_str()); 

  // volume-specific parameters
  if (ndims == 3) {
    // voxel size
    std::vector<float> voxsize = source.attr("voxsize").cast<std::vector<float>>();
    mri->xsize = voxsize[0];
    mri->ysize = voxsize[1];
    mri->zsize = voxsize[2];

    // set affine matrix
    mri->ras_good_flag = 1;
    affinematrix pyaffine = source.attr("affine").cast<affinematrix>();
    const double *affine = pyaffine.data();
    double xr = affine[0]; double yr = affine[1]; double zr = affine[2];  double pr = affine[3];
    double xa = affine[4]; double ya = affine[5]; double za = affine[6];  double pa = affine[7];
    double xs = affine[8]; double ys = affine[9]; double zs = affine[10]; double ps = affine[11];
    double sx = std::sqrt(xr * xr + xa * xa + xs * xs);
    double sy = std::sqrt(yr * yr + ya * ya + ys * ys);
    double sz = std::sqrt(zr * zr + za * za + zs * zs);
    mri->x_r = xr / sx;
    mri->x_a = xa / sx;
    mri->x_s = xs / sx;
    mri->y_r = yr / sy;
    mri->y_a = ya / sy;
    mri->y_s = ys / sy;
    mri->z_r = zr / sz;
    mri->z_a = za / sz;
    mri->z_s = zs / sz;
    MRIp0ToCRAS(mri, pr, pa, ps);
  }

  // make sure to register the new MRI instance in the bridge
  setmri(mri);
  return mri;
}


/*
  Reads a python MultidimArray from file via the MRI bridge.
*/
py::object read(const std::string& filename)
{
  return Bridge(MRIread(filename.c_str()));
}


/*
  Writes a python MultidimArray to file via the MRI bridge.
*/
void write(Bridge vol, const std::string& filename)
{
  MRIwrite(vol, filename.c_str());
}


}  // end namespace vol

