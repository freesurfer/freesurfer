#include "volume.h"

namespace vol {


/*
  Builds a python Overlay, Image, or Volume object (whichever is appropriate given the dimensionality)
  from the internal MRI instance. Voxel data ownership is immediately transfered from the MRI to
  the python object once this is called.
*/
py::object Bridge::python()
{
  // sanity check on the MRI instance
  if (!p_mri) throw std::runtime_error("cannot bridge to python as MRI instance is null");
  if (!p_mri->ischunked) throw std::runtime_error("image is too large to fit into contiguous memory");

  // if the internal MRI instance was generated from a python object, it's chunked data
  // is already managed by the original buffer array; however, if the MRI was generated
  // elsewhere, we'll need to create this array now
  if (mri_buffer.size() == 0) {

    // determine numpy dtype from MRI type
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
      throw py::value_error("unknown MRI data type ID: " + std::to_string(p_mri->type));
    }

    // squeeze the 4D MRI to determine the actual represented shape
    std::vector<ssize_t> shape = {p_mri->width};
    if (p_mri->height  > 1) shape.push_back(p_mri->height);
    if (p_mri->depth   > 1) shape.push_back(p_mri->depth);
    if (p_mri->nframes > 1) shape.push_back(p_mri->nframes);
    std::vector<ssize_t> strides = fstrides(shape, p_mri->bytes_per_vox);

    // wrap a numpy array around the chunked MRI data
    py::capsule capsule(p_mri->chunk, [](void *p) { free(p); } );
    mri_buffer = py::array(dtype, shape, strides, p_mri->chunk, capsule);

    // the above array is now fully responsible for managing the chunk data,
    // so we must remove ownership from the MRI to avoid unwanted deletion
    p_mri->owndata = false;
  }

  // extract base dimensions (ignore frames) to determine whether the MRI
  // represents an overlay, image, or volume
  int ndims = (p_mri->nframes > 1) ? mri_buffer.ndim() - 1 : mri_buffer.ndim();

  // construct the appropriate python object knowing ndims
  py::module fsmodule = py::module::import("freesurfer");
  py::object pyobj;
  switch (ndims) {
  case 1: pyobj = fsmodule.attr("Overlay")(mri_buffer); break;
  case 2: pyobj = fsmodule.attr("Image")(mri_buffer); break;
  case 3: pyobj = fsmodule.attr("Volume")(mri_buffer); break;
  }

  // transfer the rest of the parameters
  transferParameters(pyobj);
  return pyobj;
}


/*
  Transfers parameters (except data array) from the internal MRI to the provided python instance.
*/
void Bridge::transferParameters(py::object& pyobj)
{
  // grab number of basedims
  const int ndims = pyobj.attr("basedims").cast<int>();

  // extract the affine transform if it's an image or volume
  if ((ndims == 2) || (ndims == 3)) {
    pyobj.attr("voxsize") = py::make_tuple(p_mri->xsize, p_mri->ysize, p_mri->zsize);
    if (p_mri->ras_good_flag == 1) {
      MATRIX *matrix = extract_i_to_r(p_mri.get());
      py::array affine = copyArray({4, 4}, MemoryOrder::C, matrix->data);
      MatrixFree(&matrix);
      pyobj.attr("affine") = affine;
    } else {
      pyobj.attr("affine") = py::none();
    }
  }

  // transfer lookup table if it exists
  pyobj.attr("lut") = py::none();
  if (p_mri->ct) {
    py::module fsmodule = py::module::import("freesurfer");
    py::object lut = fsmodule.attr("LookupTable")();
    py::object addfunc = lut.attr("add");
    for (int i = 0; i < p_mri->ct->nentries; i++) {
      if (!p_mri->ct->entries[i]) continue;
      std::string name = std::string(p_mri->ct->entries[i]->name);
      int r = p_mri->ct->entries[i]->ri;
      int g = p_mri->ct->entries[i]->gi;
      int b = p_mri->ct->entries[i]->bi;
      int a = p_mri->ct->entries[i]->ai;
      addfunc(i, name, py::make_tuple(r, g, b, a));
    }
    pyobj.attr("lut") = lut;
  }

  // transfer scan parameters if volume
  if (ndims == 3) {
    if (p_mri->te == 0) { pyobj.attr("te") = py::none(); } else { pyobj.attr("te") = p_mri->te; }
    if (p_mri->tr == 0) { pyobj.attr("tr") = py::none(); } else { pyobj.attr("tr") = p_mri->tr; }
    if (p_mri->ti == 0) { pyobj.attr("ti") = py::none(); } else { pyobj.attr("ti") = p_mri->ti; }
    if (p_mri->flip_angle == 0) { pyobj.attr("flip_angle") = py::none(); } else { pyobj.attr("flip_angle") = p_mri->flip_angle; }
  }
}


/*
  Updates the cached python object with the internal MRI instance.
*/
void Bridge::updateSource()
{
  // sanity check on the MRI instance and python source
  if (!p_mri) throw std::runtime_error("cannot bridge to python as MRI instance is null");
  if (source.is(py::none())) throw py::value_error("cannot update source if it does not exist");

  // update the data array and let transferParameters() do the rest
  source.attr("data") = mri_buffer;
  transferParameters(source);
}


/*
  Returns the internal MRI instance. If one does not exist, it will be created
  from the cached python source object.
*/
MRI* Bridge::mri()
{
  // return if the MRI instance has already been set or created
  if (p_mri) return p_mri.get();

  // make sure the source python object has been provided
  if (source.is(py::none())) throw py::value_error("cannot generate MRI instance without source object");

  // make sure buffer is in fortran order and cache the array in case we're converting back to python later
  py::module np = py::module::import("numpy");
  mri_buffer = np.attr("asfortranarray")(source.attr("data"));

  // convert unsupported data types
  if (py::isinstance<py::array_t<bool>>(mri_buffer)) mri_buffer = py::array_t<uchar>(mri_buffer);
  if (py::isinstance<py::array_t<long>>(mri_buffer)) mri_buffer = py::array_t<int>(mri_buffer);
  if (py::isinstance<py::array_t<double>>(mri_buffer)) mri_buffer = py::array_t<float>(mri_buffer);

  // determine valid MRI type from numpy datatype
  int dtype;
  if      (py::isinstance<py::array_t<uchar>>(mri_buffer)) { dtype = MRI_UCHAR; }
  else if (py::isinstance<py::array_t<int>>  (mri_buffer)) { dtype = MRI_INT; }
  else if (py::isinstance<py::array_t<short>>(mri_buffer)) { dtype = MRI_SHORT; }
  else if (py::isinstance<py::array_t<long>> (mri_buffer)) { dtype = MRI_LONG; }
  else if (py::isinstance<py::array_t<float>>(mri_buffer)) { dtype = MRI_FLOAT; }
  else {
    throw py::value_error("unsupported array dtype " + py::str(mri_buffer.attr("dtype")).cast<std::string>());
  }

  // initialize a header-only MRI structure with the known shape (expanded to 4D)
  std::vector<int> expanded, shape = mri_buffer.attr("shape").cast<std::vector<int>>();
  int nframes = source.attr("nframes").cast<int>();
  int ndims = source.attr("basedims").cast<int>();
  switch (ndims) {
  case 1: expanded = {shape[0], 1,        1,        nframes}; break;
  case 2: expanded = {shape[0], shape[1], 1,        nframes}; break;
  case 3: expanded = {shape[0], shape[1], shape[2], nframes}; break;
  }
  MRI *mri = new MRI(expanded, dtype, false);

  // point the MRI chunk to the numpy array data and finish initializing
  mri->chunk = mri_buffer.mutable_data();
  mri->ischunked = true;
  mri->owndata = false;
  mri->initSlices();
  mri->initIndices();

  if ((ndims == 2) || (ndims == 3)) {
    // voxel size
    std::vector<float> voxsize = source.attr("voxsize").cast<std::vector<float>>();
    mri->xsize = voxsize[0];
    mri->ysize = voxsize[1];
    mri->zsize = voxsize[2];

    // set the affine transform (must come after setting voxel size)
    py::object pyaffine = source.attr("affine");
    if (pyaffine.is(py::none())) {
      mri->ras_good_flag = 0;
    } else {
      mri->ras_good_flag = 1;
      // ensure it's a c-order double array
      arrayc<double> casted = pyaffine.cast<arrayc<double>>();
      const double *affine = casted.data();
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
  }

  // transfer lookup table if it exists
  // pretty crazy that it requires this many lines of code...
  py::object lookup = source.attr("lut");
  if (!lookup.is(py::none())) {
    py::dict lut = lookup;  // cast to dict
    COLOR_TABLE *ctab = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
    int maxidx = 0;
    for (auto item : lut) {
      int curridx = item.first.cast<int>();
      if (curridx > maxidx) maxidx = curridx;
    }
    ctab->nentries = maxidx + 1;
    ctab->entries = (COLOR_TABLE_ENTRY **)calloc(ctab->nentries, sizeof(COLOR_TABLE_ENTRY *));
    ctab->version = 2;
    for (auto item : lut) {
      int idx = item.first.cast<int>();
      std::string name = item.second.attr("name").cast<std::string>();
      std::vector<int> color = item.second.attr("color").cast<std::vector<int>>();
      if (color.size() != 4) throw std::runtime_error("lookup table colors must be RGBA");
      ctab->entries[idx] = (CTE *)malloc(sizeof(CTE));
      strncpy(ctab->entries[idx]->name, name.c_str(), sizeof(ctab->entries[idx]->name));
      ctab->entries[idx]->ri = color[0];
      ctab->entries[idx]->gi = color[1];
      ctab->entries[idx]->bi = color[2];
      ctab->entries[idx]->ai = color[3];
      ctab->entries[idx]->rf = (float)ctab->entries[idx]->ri / 255.0;
      ctab->entries[idx]->gf = (float)ctab->entries[idx]->gi / 255.0;
      ctab->entries[idx]->bf = (float)ctab->entries[idx]->bi / 255.0;
      ctab->entries[idx]->af = (float)ctab->entries[idx]->ai / 255.0;
      ctab->entries[idx]->TissueType = 0;
      ctab->entries[idx]->count = 0;
    }
    mri->ct = ctab;
  }

  // transfer scan parameters if volume
  if (ndims == 3) {
    if (!py::object(source.attr("te")).is(py::none())) mri->te = source.attr("te").cast<float>();
    if (!py::object(source.attr("tr")).is(py::none())) mri->tr = source.attr("tr").cast<float>();
    if (!py::object(source.attr("ti")).is(py::none())) mri->ti = source.attr("ti").cast<float>();
    if (!py::object(source.attr("flip_angle")).is(py::none())) mri->flip_angle = source.attr("flip_angle").cast<double>();
  }

  // make sure to register the new MRI instance in the bridge
  setmri(mri);
  return mri;
}


/*
  Reads a python object from file via the MRI bridge.
*/
py::object read(const std::string& filename)
{
  return Bridge(MRIread(filename.c_str()));
}


/*
  Writes a python object to file via the MRI bridge.
*/
void write(Bridge vol, const std::string& filename)
{
  MRIwrite(vol, filename.c_str());
}


}  // end namespace vol
