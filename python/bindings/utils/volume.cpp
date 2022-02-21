#include "volume.h"
#include "annotation.h"


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

  // extract resolution
  if (ndims == 2) pyobj.attr("pixsize") = py::make_tuple(p_mri->xsize, p_mri->ysize);
  if (ndims == 3) pyobj.attr("voxsize") = py::make_tuple(p_mri->xsize, p_mri->ysize, p_mri->zsize);
  
  // extract the affine transform if it's an image or volume
  if ((ndims == 2) || (ndims == 3)) {
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
  if (py::isinstance<py::array_t<double>>(mri_buffer)) mri_buffer = py::array_t<float>(mri_buffer);
  if (py::isinstance<py::array_t<bool>>(mri_buffer)) mri_buffer = py::array_t<uchar>(mri_buffer);
  if (py::isinstance<py::array_t<char>>(mri_buffer)) mri_buffer = py::array_t<int>(mri_buffer);
  if (py::isinstance<py::array_t<long>>(mri_buffer)) mri_buffer = py::array_t<int>(mri_buffer);

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

  // voxel size
  if (ndims == 2) {
    std::vector<float> pixsize = source.attr("pixsize").cast<std::vector<float>>();
    mri->xsize = pixsize[0];
    mri->ysize = pixsize[1];
    mri->zsize = 1.0;
  } else if (ndims == 3) {
    std::vector<float> voxsize = source.attr("voxsize").cast<std::vector<float>>();
    mri->xsize = voxsize[0];
    mri->ysize = voxsize[1];
    mri->zsize = voxsize[2];
  }

  // set the affine transform (must come after setting voxel size)
  if ((ndims == 2) || (ndims == 3)) {
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
  if (stringEndsWith(filename, ".annot")) {
    return Bridge(readAnnotationIntoSeg(filename));
  } else {
    return Bridge(MRIread(filename.c_str()));
  }
}


/*
  Writes a python object to file via the MRI bridge.
*/
void write(Bridge vol, const std::string& filename)
{
  if (stringEndsWith(filename, ".annot")) {
    writeAnnotationFromSeg(vol, filename);
  } else {
    MRIwrite(vol, filename.c_str());
  }
}


/*
  Templated function to interpolate and add values at given coordinates into a typed volume.
  See sampleIntoVolume() for full documentation.

  TODO right now 3D and 4D have separate implementations - which is not super ideal
*/
template <typename T>
static void sampleIntoTypedVolume(py::array_t<T> volume, py::array_t<T> weights, arrayc<double> coords, arrayc<double> values)
{
  const double * coords_ptr = coords.data();
  const double * vals_ptr = values.data();

  int xlim = volume.shape(0);
  int ylim = volume.shape(1);
  int zlim = volume.shape(2);

  int nvalues = values.shape(0); 

  if (volume.ndim() == 3) {

    // 3D sampling
    auto volume_m = volume.template mutable_unchecked<3>();
    auto weights_m = weights.template mutable_unchecked<3>();

    for (int n = 0; n < nvalues; n++) {
      double x = *coords_ptr++;
      double y = *coords_ptr++;
      double z = *coords_ptr++;

      // get low and high coordinates
      int x_low = int(floor(x));
      int y_low = int(floor(y));
      int z_low = int(floor(z));
      int x_high = x_low + 1;
      int y_high = y_low + 1;
      int z_high = z_low + 1;
      if (x_high == xlim) x_high = x_low;
      if (y_high == ylim) y_high = y_low;
      if (z_high == zlim) z_high = z_low;

      // get coordinate diff
      x -= x_low;
      y -= y_low;
      z -= z_low;
      double dx = 1.0f - x;
      double dy = 1.0f - y;
      double dz = 1.0f - z;

      // compute interpolation weights
      T w0 = dx * dy * dz;
      T w1 = x  * dy * dz;
      T w2 = dx * y  * dz;
      T w3 = dx * dy * z;
      T w4 = x  * dy * z;
      T w5 = dx * y  * z;
      T w6 = x  * y  * dz;
      T w7 = x  * y  * z;

      // add interpolated values
      T val = *vals_ptr++;
      volume_m(x_low , y_low , z_low ) += w0 * val;
      volume_m(x_high, y_low , z_low ) += w1 * val;
      volume_m(x_low , y_high, z_low ) += w2 * val;
      volume_m(x_low , y_low , z_high) += w3 * val;
      volume_m(x_high, y_low , z_high) += w4 * val;
      volume_m(x_low , y_high, z_high) += w5 * val;
      volume_m(x_high, y_high, z_low ) += w6 * val;
      volume_m(x_high, y_high, z_high) += w7 * val;

      // add weights
      weights_m(x_low , y_low , z_low ) += w0;
      weights_m(x_high, y_low , z_low ) += w1;
      weights_m(x_low , y_high, z_low ) += w2;
      weights_m(x_low , y_low , z_high) += w3;
      weights_m(x_high, y_low , z_high) += w4;
      weights_m(x_low , y_high, z_high) += w5;
      weights_m(x_high, y_high, z_low ) += w6;
      weights_m(x_high, y_high, z_high) += w7;
    }
  } 
  else if (volume.ndim() == 4) {

    // 4D sampling
    int nframes = volume.shape(3);
    auto volume_m = volume.template mutable_unchecked<4>();
    auto weights_m = weights.template mutable_unchecked<3>();

    for (int n = 0; n < nvalues; n++) {
      double x = *coords_ptr++;
      double y = *coords_ptr++;
      double z = *coords_ptr++;

      // get low and high coordinates
      int x_low = int(floor(x));
      int y_low = int(floor(y));
      int z_low = int(floor(z));
      int x_high = x_low + 1;
      int y_high = y_low + 1;
      int z_high = z_low + 1;
      if (x_high == xlim) x_high = x_low;
      if (y_high == ylim) y_high = y_low;
      if (z_high == zlim) z_high = z_low;

      // get coordinate diff
      x -= x_low;
      y -= y_low;
      z -= z_low;
      double dx = 1.0f - x;
      double dy = 1.0f - y;
      double dz = 1.0f - z;

      // compute interpolation weights
      T w0 = dx * dy * dz;
      T w1 = x  * dy * dz;
      T w2 = dx * y  * dz;
      T w3 = dx * dy * z;
      T w4 = x  * dy * z;
      T w5 = dx * y  * z;
      T w6 = x  * y  * dz;
      T w7 = x  * y  * z;

      // add interpolated values
      for (int fno = 0; fno < nframes; fno++) {
        T val = *vals_ptr++;
        volume_m(x_low , y_low , z_low , fno) += w0 * val;
        volume_m(x_high, y_low , z_low , fno) += w1 * val;
        volume_m(x_low , y_high, z_low , fno) += w2 * val;
        volume_m(x_low , y_low , z_high, fno) += w3 * val;
        volume_m(x_high, y_low , z_high, fno) += w4 * val;
        volume_m(x_low , y_high, z_high, fno) += w5 * val;
        volume_m(x_high, y_high, z_low , fno) += w6 * val;
        volume_m(x_high, y_high, z_high, fno) += w7 * val;
      }

      // add weights
      weights_m(x_low , y_low , z_low ) += w0;
      weights_m(x_high, y_low , z_low ) += w1;
      weights_m(x_low , y_high, z_low ) += w2;
      weights_m(x_low , y_low , z_high) += w3;
      weights_m(x_high, y_low , z_high) += w4;
      weights_m(x_low , y_high, z_high) += w5;
      weights_m(x_high, y_high, z_low ) += w6;
      weights_m(x_high, y_high, z_high) += w7;
    }
  }
  else {
    fs::fatal() << "target volume must be 3D or 4D";
  }
}


/*
  Interpolates and adds values at given coordinates into a volume and corresponding
  weights array. Since this function performs in-place modification of the volume and
  weights arrays, they MUST be double or float precision to avoid being copied.

  \param volume   Float or double array to add sampled values into.
  \param weights  Float or double array to add sampling weights into.
  \param coords   List of volume sampling coordinates. 
  \param values   List of values at each coordinate.
*/
void sampleIntoVolume(py::array volume, py::array weights, arrayc<double> coords, arrayc<double> values)
{
  if (py::isinstance<py::array_t<double>>(volume)) {
    // input is double array
    if (!py::isinstance<py::array_t<double>>(weights)) {
      fs::fatal() << "sampling weights volume must match volume type (double)";
    }
    sampleIntoTypedVolume<double>(volume, weights, coords, values);
  }
  else if (py::isinstance<py::array_t<float>>(volume)) {
    // input is float array
    if (!py::isinstance<py::array_t<float>>(weights)) {
      fs::fatal() << "sampling weights volume must match volume type (float)";
    }
    sampleIntoTypedVolume<float>(volume, weights, coords, values);
  }
  else {
    fs::fatal() << "sampling target volume must be of type float or double";
  }
}


/*
  Templated trilinear resampling function to support source inputs with either float or double precision.
  See resampleVolume() for full documentation.
*/
template <typename T>
static py::array_t<T> resampleTypedVolumeLinear(py::array_t<T> source_vol, py::object target_shape, py::array_t<double> target2source, T fill)
{
  // make target volume
  py::module numpy = py::module::import("numpy");
  py::array_t<T> target_vol = numpy.attr("zeros")(target_shape);  // TODO construct with dtype instead
  int nx = target_vol.shape(0);
  int ny = target_vol.shape(1);
  int nz = target_vol.shape(2);
  int nf = target_vol.shape(3);

  // get indexing functions
  auto source_m = source_vol.template unchecked<4>();
  auto target_m = target_vol.template mutable_unchecked<4>();

  // cache affine values
  auto t2s_m = target2source.unchecked<2>();
  double a11 = t2s_m(0, 0); double a12 = t2s_m(0, 1); double a13 = t2s_m(0, 2); double a14 = t2s_m(0, 3);
  double a21 = t2s_m(1, 0); double a22 = t2s_m(1, 1); double a23 = t2s_m(1, 2); double a24 = t2s_m(1, 3);
  double a31 = t2s_m(2, 0); double a32 = t2s_m(2, 1); double a33 = t2s_m(2, 2); double a34 = t2s_m(2, 3);

  // get source limits
  int sxlim = source_vol.shape(0);
  int sylim = source_vol.shape(1);
  int szlim = source_vol.shape(2);

  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {

        // compute source coordinate
        T sx = (a11 * x) + (a12 * y) + (a13 * z) + a14;
        T sy = (a21 * x) + (a22 * y) + (a23 * z) + a24;
        T sz = (a31 * x) + (a32 * y) + (a33 * z) + a34;

        // get low and high coords
        int sx_low = int(floor(sx));
        int sy_low = int(floor(sy));
        int sz_low = int(floor(sz));
        int sx_high = sx_low + 1;
        int sy_high = sy_low + 1;
        int sz_high = sz_low + 1;

        // get coordinate diff
        sx -= sx_low;
        sy -= sy_low;
        sz -= sz_low;
        T dsx = 1.0f - sx;
        T dsy = 1.0f - sy;
        T dsz = 1.0f - sz;

        // make sure voxels are within the source volume
        bool valid_xl = (sx_low >= 0 && sx_low < sxlim);
        bool valid_yl = (sy_low >= 0 && sy_low < sylim);
        bool valid_zl = (sz_low >= 0 && sz_low < szlim);
        bool valid_xh = (sx_high >= 0 && sx_high < sxlim);
        bool valid_yh = (sy_high >= 0 && sy_high < sylim);
        bool valid_zh = (sz_high >= 0 && sz_high < szlim);

        T w0 = dsx * dsy * dsz;
        T w1 = sx  * dsy * dsz;
        T w2 = dsx * sy  * dsz;
        T w3 = dsx * dsy * sz;
        T w4 = sx  * dsy * sz;
        T w5 = dsx * sy  * sz;
        T w6 = sx  * sy  * dsz;
        T w7 = sx  * sy  * sz;

        for (int f = 0; f < nf; f++) {

          // extract (valid) voxel values
          T v0 = (valid_xl && valid_yl && valid_zl) ? source_m(sx_low , sy_low , sz_low , f) : fill;
          T v1 = (valid_xh && valid_yl && valid_zl) ? source_m(sx_high, sy_low , sz_low , f) : fill;
          T v2 = (valid_xl && valid_yh && valid_zl) ? source_m(sx_low , sy_high, sz_low , f) : fill;
          T v3 = (valid_xl && valid_yl && valid_zh) ? source_m(sx_low , sy_low , sz_high, f) : fill;
          T v4 = (valid_xh && valid_yl && valid_zh) ? source_m(sx_high, sy_low , sz_high, f) : fill;
          T v5 = (valid_xl && valid_yh && valid_zh) ? source_m(sx_low , sy_high, sz_high, f) : fill;
          T v6 = (valid_xh && valid_yh && valid_zl) ? source_m(sx_high, sy_high, sz_low , f) : fill;
          T v7 = (valid_xh && valid_yh && valid_zh) ? source_m(sx_high, sy_high, sz_high, f) : fill;

          // interpolate
          T val = w0 * v0 +
                  w1 * v1 +
                  w2 * v2 +
                  w3 * v3 +
                  w4 * v4 +
                  w5 * v5 +
                  w6 * v6 +
                  w7 * v7;

          target_m(x, y, z, f) = val;
        }
      }
    }
  }

  return target_vol;
}


/*
  Templated nearest neighbor resampling function to support source inputs with either float or double precision.
  See resampleVolume() for full documentation.
*/
template <typename T>
static py::array_t<T> resampleTypedVolumeNearest(py::array_t<T> source_vol, py::object target_shape, py::array_t<double> target2source, T fill)
{
  // make target volume
  py::module numpy = py::module::import("numpy");
  py::array_t<T> target_vol = numpy.attr("zeros")(target_shape);  // TODO construct with dtype instead
  if (fill != 0) target_vol.attr("fill")(fill);
  int nx = target_vol.shape(0);
  int ny = target_vol.shape(1);
  int nz = target_vol.shape(2);
  int nf = target_vol.shape(3);

  // get indexing functions
  auto source_m = source_vol.template unchecked<4>();
  auto target_m = target_vol.template mutable_unchecked<4>();

  // cache affine values
  auto t2s_m = target2source.unchecked<2>();
  double a11 = t2s_m(0, 0); double a12 = t2s_m(0, 1); double a13 = t2s_m(0, 2); double a14 = t2s_m(0, 3);
  double a21 = t2s_m(1, 0); double a22 = t2s_m(1, 1); double a23 = t2s_m(1, 2); double a24 = t2s_m(1, 3);
  double a31 = t2s_m(2, 0); double a32 = t2s_m(2, 1); double a33 = t2s_m(2, 2); double a34 = t2s_m(2, 3);

  // get source limits
  int sxlim = source_vol.shape(0);
  int sylim = source_vol.shape(1);
  int szlim = source_vol.shape(2);

  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {

        // compute source coordinate
        T sx = (a11 * x) + (a12 * y) + (a13 * z) + a14;
        T sy = (a21 * x) + (a22 * y) + (a23 * z) + a24;
        T sz = (a31 * x) + (a32 * y) + (a33 * z) + a34;

        // get rounded coord and make sure voxels are within the source volume
        int rsx = std::round(sx);
        if (rsx < 0 || rsx >= sxlim) continue;

        int rsy = std::round(sy);
        if (rsy < 0 || rsy >= sylim) continue;

        int rsz = std::round(sz);
        if (rsz < 0 || rsz >= szlim) continue;

        for (int f = 0; f < nf; f++) target_m(x, y, z, f) = source_m(rsx, rsy, rsz, f);
      }
    }
  }

  return target_vol;
}


/*
  Resamples (trilinear) a volume array from one space to another given a target-to-source
  transformation matrix. The source vol can be of any type, but only float and
  double precision arrays are used directly. Arrays of any other dtype will be
  converted to float arrays.

  \param source_vol     Source array to sample from. Must be 4D.
  \param target_shape   Shape of the returned target array.
  \param target2source  4x4 affine matrix that transforms target coords to source coords.
*/
py::array resampleVolumeLinear(py::array source_vol, py::object target_shape, py::array_t<double> target2source, double fill)
{
  if (py::isinstance<py::array_t<double>>(source_vol)) {
    return resampleTypedVolumeLinear<double>(source_vol, target_shape, target2source, fill);
  } else {
    return resampleTypedVolumeLinear<float>(source_vol, target_shape, target2source, fill);
  }
}


/*
  Resamples (nearest) a volume array from one space to another given a target-to-source
  transformation matrix. The source vol can be of any type, but only float and
  double precision arrays are used directly. Arrays of any other dtype will be
  converted to float arrays.

  \param source_vol     Source array to sample from. Must be 4D.
  \param target_shape   Shape of the returned target array.
  \param target2source  4x4 affine matrix that transforms target coords to source coords.
*/
py::array resampleVolumeNearest(py::array source_vol, py::object target_shape, py::array_t<double> target2source, double fill)
{
  if (py::isinstance<py::array_t<double>>(source_vol)) {
    return resampleTypedVolumeNearest<double>(source_vol, target_shape, target2source, fill);
  } else {
    return resampleTypedVolumeNearest<float>(source_vol, target_shape, target2source, fill);
  }
}


}  // end namespace vol
