#include "annotation.h"
#include "mri_circulars.h"

#include "bindings_transform.h"
#include "bindings_mri.h"


/*
  Build a surfa Overlay, Slice, or Volume object (whichever is appropriate given the dimensionality)
  from an MRI instance. All data is copied between python and cxx, so any allocated MRI pointers
  will need to be freed, even after convert to python. However, if `release` is `true`, this function
  will free the MRI after converting.
*/
py::object MRItoSurfaArray(MRI* mri, bool release)
{
  // sanity check on the MRI instance
  if (!mri) throw std::runtime_error("MRItoSurfaArray: cannot convert to surfa - MRI input is null");
  if (!mri->ischunked) throw std::runtime_error("MRItoSurfaArray: image is too large to fit into contiguous memory");

  // determine numpy dtype from MRI type
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
  case MRI_USHRT:
    dtype = py::dtype::of<unsigned short>(); break;
  default:
    throw py::value_error("MRItoSurfaArray: unknown MRI data type ID: " + std::to_string(mri->type));
  }

  // squeeze the 4D MRI to determine the actual represented shape
  std::vector<ssize_t> shape = {mri->width};
  if (mri->height  > 1) shape.push_back(mri->height);
  if (mri->depth   > 1) shape.push_back(mri->depth);
  if (mri->nframes > 1) shape.push_back(mri->nframes);
  std::vector<ssize_t> strides = fstrides(shape, mri->bytes_per_vox);

  // wrap a numpy array around the chunked MRI data, then copy
  py::capsule capsule(mri->chunk);
  py::array buffer = py::array(dtype, shape, strides, mri->chunk, capsule).attr("copy")();

  // extract base dimensions (ignore frames) to determine whether the MRI
  // represents an overlay, image, or volume
  int ndims = (mri->nframes > 1) ? buffer.ndim() - 1 : buffer.ndim();

  // construct the appropriate FramedArray knowing ndims
  py::module fsmodule = py::module::import("surfa");
  py::object arr;
  switch (ndims) {
  case 1: arr = fsmodule.attr("Overlay")(buffer); break;
  case 2: arr = fsmodule.attr("Slice")(buffer);   break;
  case 3: arr = fsmodule.attr("Volume")(buffer);  break;
  }
  
  // extract the affine transform if it's an image or volume
  if (ndims > 1) {

    // extract resolution
    arr.attr("voxsize") = py::make_tuple(mri->xsize, mri->ysize, mri->zsize);
    
    // geometry
    if (mri->ras_good_flag == 1) {
      VOL_GEOM vg;
      MRIcopyVolGeomFromMRI(mri, &vg);
      arr.attr("geom") = VOLGEOMtoSurfaImageGeometry(&vg);
    }

    // extract scan parameters
    py::dict metadata = arr.attr("metadata");
    if (mri->te != 0) { metadata["te"] = mri->te; }
    if (mri->tr != 0) { metadata["tr"] = mri->tr; }
    if (mri->ti != 0) { metadata["ti"] = mri->ti; }
    if (mri->flip_angle != 0) { metadata["flip_angle"] = mri->flip_angle; }
  }

  // transfer lookup table if it exists
  if (mri->ct) {
    py::dict lookup = py::module::import("surfa").attr("LabelLookup")();
    for (int i = 0; i < mri->ct->nentries; i++) {
      if (!mri->ct->entries[i]) continue;
      std::string name = std::string(mri->ct->entries[i]->name);
      int r = mri->ct->entries[i]->ri;
      int g = mri->ct->entries[i]->gi;
      int b = mri->ct->entries[i]->bi;
      float a = mri->ct->entries[i]->ai / 255;
      lookup[py::int_{i}] = py::make_tuple(name, py::make_tuple(r, g, b, a));
    }
    arr.attr("labels") = lookup;
  }

  if (release) MRIfree(&mri);
}


/*
  Convert a surfa FramedArray to an MRI structure of appropriate dimensionality. The returned
  MRI pointer will need to be freed manually once it's done with.
*/
MRI* MRIfromSurfaArray(py::object arr)
{
  // type checking
  py::object arrclass = py::module::import("surfa").attr("core").attr("FramedArray");
  if (!py::isinstance(arr, arrclass)) throw py::value_error("MRIfromSurfaArray: cannot convert to MRI - input is not a surfa FramedArray");

  // make sure buffer is in fortran order and cache the array in case we're converting back to surfa later
  py::module np = py::module::import("numpy");
  py::array mri_buffer = np.attr("asfortranarray")(arr.attr("data"));

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
  else if (py::isinstance<py::array_t<unsigned short>>(mri_buffer)) { dtype = MRI_USHRT; }
  else {
    throw py::value_error("MRIfromSurfaArray: unsupported array dtype " + py::str(mri_buffer.attr("dtype")).cast<std::string>());
  }

  // initialize a header-only MRI structure with the known shape (expanded to 4D)
  std::vector<int> expanded, shape = mri_buffer.attr("shape").cast<std::vector<int>>();
  int nframes = arr.attr("nframes").cast<int>();

  // get dimensionality
  int ndims = arr.attr("basedim").cast<int>();
  switch (ndims) {
  case 1: expanded = {shape[0], 1,        1,        nframes}; break;
  case 2: expanded = {shape[0], shape[1], 1,        nframes}; break;
  case 3: expanded = {shape[0], shape[1], shape[2], nframes}; break;
  }
  MRI *mri = new MRI(expanded, dtype, false);

  // copy buffer data into MRI chunk
  memcpy(mri->chunk, mri_buffer.mutable_data(), mri->bytes_total);
  mri->ischunked = true;
  mri->initSlices();
  mri->initIndices();

  if (ndims > 1) {
    // set the image geometry
    mri->ras_good_flag = 1;
    VOL_GEOM vg;
    VOLGEOMfromSurfaImageGeometry(arr.attr("geom"), &vg);
    MRIcopyVolGeomToMRI(mri, &vg);

    // transfer scan parameters if volume
    py::dict metadata = arr.attr("metadata");
    if (!py::object(metadata.attr("get")("te")).is(py::none())) mri->te = metadata["te"].cast<float>();
    if (!py::object(metadata.attr("get")("tr")).is(py::none())) mri->tr = metadata["tr"].cast<float>();
    if (!py::object(metadata.attr("get")("ti")).is(py::none())) mri->ti = metadata["ti"].cast<float>();
    if (!py::object(metadata.attr("get")("flip_angle")).is(py::none())) mri->flip_angle = metadata["flip_angle"].cast<double>();
  }

  // transfer lookup table if it exists... pretty crazy that it requires this much code
  py::dict labels = arr.attr("labels");
  if (!labels.is(py::none())) {
    COLOR_TABLE *ctab = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
    int maxidx = 0;
    for (auto item : labels) {
      int curridx = item.first.cast<int>();
      if (curridx > maxidx) maxidx = curridx;
    }
    ctab->nentries = maxidx + 1;
    ctab->entries = (COLOR_TABLE_ENTRY **)calloc(ctab->nentries, sizeof(COLOR_TABLE_ENTRY *));
    ctab->version = 2;
    for (auto item : labels) {
      int idx = item.first.cast<int>();
      std::string name = item.second.attr("name").cast<std::string>();
      std::vector<float> color = item.second.attr("color").cast<std::vector<float>>();
      if (color.size() != 4) throw std::runtime_error("lookup table colors must be RGBA");
      ctab->entries[idx] = (CTE *)malloc(sizeof(CTE));
      strncpy(ctab->entries[idx]->name, name.c_str(), sizeof(ctab->entries[idx]->name));
      ctab->entries[idx]->ri = int(color[0]);
      ctab->entries[idx]->gi = int(color[1]);
      ctab->entries[idx]->bi = int(color[2]);
      ctab->entries[idx]->ai = int(color[3] * 255);
      ctab->entries[idx]->rf = (float)ctab->entries[idx]->ri / 255.0;
      ctab->entries[idx]->gf = (float)ctab->entries[idx]->gi / 255.0;
      ctab->entries[idx]->bf = (float)ctab->entries[idx]->bi / 255.0;
      ctab->entries[idx]->af = (float)ctab->entries[idx]->ai / 255.0;
      ctab->entries[idx]->TissueType = 0;
      ctab->entries[idx]->count = 0;
    }
    mri->ct = ctab;
  }

  return mri;
}


/*
  Read a surfa object from file via the MRI bridge. Surfa already covers array IO, so this
  is probably unnecessary, but might be useful at some point.
*/
py::object readMRI(const std::string& filename)
{
  if (stringEndsWith(filename, ".annot")) {
    return MRItoSurfaArray(readAnnotationIntoSeg(filename), true);
  } else {
    return MRItoSurfaArray(MRIread(filename.c_str()), true);
  }
}


/*
  Write a surfa object to file via the MRI bridge. Surfa already covers array IO, so this
  is probably unnecessary, but might be useful at some point.
*/
void writeMRI(py::object arr, const std::string& filename)
{
  MRI* mri = MRIfromSurfaArray(arr);
  if (stringEndsWith(filename, ".annot")) {
    writeAnnotationFromSeg(mri, filename);
  } else {
    MRIwrite(mri, filename.c_str());
  }
  MRIfree(&mri);
}
