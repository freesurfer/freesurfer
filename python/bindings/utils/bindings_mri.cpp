#include "annotation.h"
#include "mri_circulars.h"

#include "bindings_transform.h"
#include "bindings_mri.h"


/*
  Build a surfa Overlay, Slice, or Volume object (whichever is appropriate given the dimensionality)
  from the stored MRI instance. Voxel data ownership is immediately transfered from the MRI to
  the python object once this is called.
*/
py::object MRIBridge::toPython()
{
  // sanity check on the MRI instance
  if (!p_mri) throw std::runtime_error("cannot bridge to python as MRI instance is null");
  if (!p_mri->ischunked) throw std::runtime_error("image is too large to fit into contiguous memory");

  // if the stored MRI instance was generated from a python object, it's chunked data
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
    case MRI_USHRT:
      dtype = py::dtype::of<unsigned short>(); break;
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

  // construct the appropriate FramedArray knowing ndims
  py::module fsmodule = py::module::import("surfa");
  py::object framed;
  switch (ndims) {
  case 1: framed = fsmodule.attr("Overlay")(mri_buffer); break;
  case 2: framed = fsmodule.attr("Slice")(mri_buffer);   break;
  case 3: framed = fsmodule.attr("Volume")(mri_buffer);  break;
  }

  // transfer the rest of the parameters
  transferParameters(framed);
  return framed;
}


/*
  Transfer parameters (except data array) from the stored MRI to the provided surfa instance.
*/
void MRIBridge::transferParameters(py::object& framed)
{
  // grab number of basedims
  const int ndims = framed.attr("basedim").cast<int>();
  
  // extract the affine transform if it's an image or volume
  if (ndims > 1) {

    // extract resolution
    framed.attr("voxsize") = py::make_tuple(p_mri->xsize, p_mri->ysize, p_mri->zsize);
    
    // geometry
    if (p_mri->ras_good_flag == 1) {
      VOL_GEOM vg;
      MRIcopyVolGeomFromMRI(p_mri.get(), &vg);
      framed.attr("geom") = VOLGEOMtoSurfaImageGeometry(&vg);
    }

    // extract scan parameters
    py::dict metadata = framed.attr("metadata");
    if (p_mri->te != 0) { metadata["te"] = p_mri->te; }
    if (p_mri->tr != 0) { metadata["tr"] = p_mri->tr; }
    if (p_mri->ti != 0) { metadata["ti"] = p_mri->ti; }
    if (p_mri->flip_angle != 0) { metadata["flip_angle"] = p_mri->flip_angle; }
  }

  // transfer lookup table if it exists
  if (p_mri->ct) {
    py::dict lookup = py::module::import("surfa").attr("LabelLookup")();
    for (int i = 0; i < p_mri->ct->nentries; i++) {
      if (!p_mri->ct->entries[i]) continue;
      std::string name = std::string(p_mri->ct->entries[i]->name);
      int r = p_mri->ct->entries[i]->ri;
      int g = p_mri->ct->entries[i]->gi;
      int b = p_mri->ct->entries[i]->bi;
      float a = p_mri->ct->entries[i]->ai / 255;
      lookup[py::int_{i}] = py::make_tuple(name, py::make_tuple(r, g, b, a));
    }
    framed.attr("labels") = lookup;
  }
}


/*
  Updates the cached surfa object with the stored MRI instance.
*/
void MRIBridge::updateSource()
{
  // sanity check on the MRI instance and surfa source
  if (!p_mri) throw std::runtime_error("cannot bridge to python as MRI instance is null");
  if (source.is(py::none())) throw py::value_error("cannot update source if it does not exist");

  // update the data array and let transferParameters() do the rest
  source.attr("data") = mri_buffer;
  transferParameters(source);
}


/*
  Return the stored MRI instance. If one does not exist, it will be created
  from the cached surfa source object.
*/
MRI* MRIBridge::toMRI()
{
  // return if the MRI instance has already been set or created
  if (p_mri) return p_mri.get();

  // make sure the source surfa object has been provided
  if (source.is(py::none())) throw py::value_error("cannot generate MRI instance without source object");

  // make sure buffer is in fortran order and cache the array in case we're converting back to surfa later
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
  else if (py::isinstance<py::array_t<unsigned short>>(mri_buffer)) { dtype = MRI_USHRT; }
  else {
    throw py::value_error("unsupported array dtype " + py::str(mri_buffer.attr("dtype")).cast<std::string>());
  }

  // initialize a header-only MRI structure with the known shape (expanded to 4D)
  std::vector<int> expanded, shape = mri_buffer.attr("shape").cast<std::vector<int>>();
  int nframes = source.attr("nframes").cast<int>();

  // get dimensionality
  int ndims = source.attr("basedim").cast<int>();
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

  if (ndims > 1) {
    // set the image geometry
    mri->ras_good_flag = 1;
    VOL_GEOM vg;
    VOLGEOMfromSurfaImageGeometry(source.attr("geom"), &vg);
    MRIcopyVolGeomToMRI(mri, &vg);

    // transfer scan parameters if volume
    py::dict metadata = source.attr("metadata");
    if (!py::object(metadata.attr("get")("te")).is(py::none())) mri->te = source["te"].cast<float>();
    if (!py::object(metadata.attr("get")("tr")).is(py::none())) mri->tr = source["tr"].cast<float>();
    if (!py::object(metadata.attr("get")("ti")).is(py::none())) mri->ti = source["ti"].cast<float>();
    if (!py::object(metadata.attr("get")("flip_angle")).is(py::none())) mri->flip_angle = source["flip_angle"].cast<double>();
  }

  // transfer lookup table if it exists... pretty crazy that it requires this much code
  py::dict labels = source.attr("labels");
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

  // make sure to register the new MRI instance in the bridge
  setMRI(mri);
  return mri;
}


/*
  Read a surfa object from file via the MRI bridge. Surfa already covers array IO, so this
  is probably unnecessary, but might be useful at some point.
*/
py::object readMRI(const std::string& filename)
{
  if (stringEndsWith(filename, ".annot")) {
    return MRIBridge(readAnnotationIntoSeg(filename));
  } else {
    return MRIBridge(MRIread(filename.c_str()));
  }
}


/*
  Write a surfa object to file via the MRI bridge. Surfa already covers array IO, so this
  is probably unnecessary, but might be useful at some point.
*/
void writeMRI(MRIBridge vol, const std::string& filename)
{
  if (stringEndsWith(filename, ".annot")) {
    writeAnnotationFromSeg(vol, filename);
  } else {
    MRIwrite(vol, filename.c_str());
  }
}
