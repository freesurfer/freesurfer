#include "mrisurf_base.h"
#include "gcamorph.h"

#include "bindings_transform.h"

using namespace pybind11::literals;


/*
  Convert a surfa ImageGeometry to an FS VOL_GEOM. This function assumes vg
  is already allocated.
*/ 
void VOLGEOMfromSurfaImageGeometry(py::object geometry, VOL_GEOM* vg)
{
  if (geometry.is(py::none())) {
    // if the python objects is None, just make it an invalid geometry
    vg->valid = 0;
  
  } else if (!py::isinstance(geometry, py::module::import("surfa").attr("ImageGeometry"))) {
    // type checking
    throw py::value_error("VOLGEOMfromSurfaImageGeometry: cannot convert to VOL_GEOM - input is not a surfa ImageGeometry");
  
  } else {
    // extract shape
    std::vector<int> shape = geometry.attr("shape").cast<std::vector<int>>();
    MRI *tmp = new MRI(shape, MRI_UCHAR, false);

    // extract voxel size
    std::vector<float> vs = geometry.attr("voxsize").cast<std::vector<float>>();
    tmp->xsize = vs[0]; tmp->ysize = vs[1]; tmp->zsize = vs[2];

    // extract RAS rotation - make sure it's a c-ordered array
    py::object pyaffine = geometry.attr("rotation");
    arrayc<float> casted = pyaffine.cast<arrayc<float>>();
    const float *rot = casted.data();
    tmp->x_r = rot[0]; tmp->y_r = rot[1]; tmp->z_r = rot[2];
    tmp->x_a = rot[3]; tmp->y_a = rot[4]; tmp->z_a = rot[5];
    tmp->x_s = rot[6]; tmp->y_s = rot[7]; tmp->z_s = rot[8];

    // extract RAS center
    std::vector<float> center = geometry.attr("center").cast<std::vector<float>>();
    tmp->c_r = center[0]; tmp->c_a = center[1]; tmp->c_s = center[2];

    // copy data from temporary MRI and release
    MRIcopyVolGeomFromMRI(tmp, vg);
    MRIfree(&tmp);
    vg->valid = 1;
  }
}


/*
  Convert an FS VOL_GEOM to a surfa ImageGeometry. Returns None if geometry
  is not valid.
*/ 
py::object VOLGEOMtoSurfaImageGeometry(VOL_GEOM* vg)
{
  if (vg == nullptr) throw py::value_error("VOLGEOMtoSurfaImageGeometry: cannot convert to surfa ImageGeometry - VOL_GEOM input is null");

  // return None object is geom is invalid
  if (!vg->valid) return py::none();

  // extract basic parameters
  py::object shape = py::make_tuple(vg->width, vg->height, vg->depth);
  py::object voxsize = py::make_tuple(vg->xsize, vg->ysize, vg->zsize);

  // use temporary MRI structure so we can do easy rotation and center extraction
  MRI *tmp = new MRI({vg->width, vg->height, vg->depth}, MRI_UCHAR, false);
  MRIcopyVolGeomToMRI(tmp, vg);

  // extract rotation
  float* buffer = new float[9];
  buffer[0] = tmp->x_r; buffer[1] = tmp->y_r; buffer[2] = tmp->z_r;
  buffer[3] = tmp->x_a; buffer[4] = tmp->y_a; buffer[5] = tmp->z_a;
  buffer[6] = tmp->x_s; buffer[7] = tmp->y_s; buffer[8] = tmp->z_s;
  py::object rotation = makeArray({3, 3}, MemoryOrder::C, buffer);

  // extract RAS center
  py::object center = py::make_tuple(tmp->c_r, tmp->c_a, tmp->c_s);
  
  MRIfree(&tmp);

  // construct and return ImageGeometry
  py::object geometry = py::module::import("surfa").attr("ImageGeometry")(
    "shape"_a=shape,
    "voxsize"_a=voxsize,
    "rotation"_a=rotation,
    "center"_a=center);
  return geometry;
}


/*
  Convert a surfa Affine to an FS LTA.
*/ 
LTA* LTAfromSurfaAffine(py::object affine)
{
  // type checking
  py::object affineclass = py::module::import("surfa").attr("Affine");
  if (!py::isinstance(affine, affineclass)) throw py::value_error("LTAfromSurfaAffine: cannot convert to LTA - input is not a surfa Affine");

  // allocate an empty LTA
  LTA* lta = LTAalloc(1, nullptr);
  LINEAR_TRANSFORM *lt = &lta->xforms[0];

  // grab the affine matrix data and ensure c-order
  arrayc<double> casted = affine.attr("matrix").cast<arrayc<double>>();
  const double *matrix = casted.data();
  for (int i = 0; i < 16; i++) lt->m_L->data[i] = matrix[i];

  // extract "type" which is really the space of the affine transform
  py::object space = affine.attr("space");
  if (!space.is(py::none())) {
    std::string spacename = space.attr("name").cast<std::string>();
    if (spacename == "voxel") lta->type = LINEAR_VOX_TO_VOX;
    if (spacename == "world") lta->type = LINEAR_RAS_TO_RAS;
  }

  // finally extract the source and target information
  VOLGEOMfromSurfaImageGeometry(affine.attr("source"), &lt->src);
  VOLGEOMfromSurfaImageGeometry(affine.attr("target"), &lt->dst);

  return lta;
}


/*
  Convert an FS LTA to a surfa Affine.
*/ 
py::object LTAtoSurfaAffine(const LTA* lta)
{
  if (lta == nullptr) throw py::value_error("LTAtoSurfaAffine: cannot convert to surfa Affine - LTA input is null");

  // copy matrix data
  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  py::array matrix = copyArray({4, 4}, MemoryOrder::C, lt->m_L->data);

  // construct the affine then add other parameters
  py::object affine = py::module::import("surfa").attr("Affine")(matrix);

  // convert "type" to coordinate space
  if (lta->type == LINEAR_VOX_TO_VOX) {
    affine.attr("space") = "voxel";
  } else if (lta->type == LINEAR_RAS_TO_RAS) {
    affine.attr("space") = "word";
  }

  // grab the source and target information
  affine.attr("source") = VOLGEOMtoSurfaImageGeometry(&lt->src);
  affine.attr("target") = VOLGEOMtoSurfaImageGeometry(&lt->dst);

  return affine;
}


/*
  Read an affine with FS LTA code. Surfa already covers affine IO, so this
  is probably unnecessary, but might be useful at some point.
*/ 
py::object readLTA(const std::string& filename)
{
  LTA *lta = LTAread(filename.c_str());
  py::object transform = LTAtoSurfaAffine(lta);
  LTAfree(&lta);
  return transform;
}


/*
  Write an affine with FS LTA code. Surfa already covers affine IO, so this
  is probably unnecessary, but might be useful at some point.
*/ 
void writeLTA(py::object transform, const std::string& filename)
{
  LTA *lta = LTAfromSurfaAffine(transform);
  LTAwrite(lta, filename.c_str());
  LTAfree(&lta);
}


/*
  Write a warp as an FS GCAMORPH file. Inputs should include a
  3D warp field, an affine transformation, and source and target geometries. 
*/ 
void writeGCAMORPHfromPython(py::object warp,
                             py::object affine,
                             py::object source,
                             py::object target,
                             const std::string& filename)
{
  // import surfa
  py::module surfa = py::module::import("surfa");

  // get warp map and allocate GCAM
  if (py::isinstance(warp, surfa.attr("Volume"))) warp = warp.attr("data");
  arrayf<double> warpdata = warp.cast<arrayf<double>>();
  GCAM *gcam = GCAMalloc(warpdata.shape(0), warpdata.shape(1), warpdata.shape(2));

  // copy geometries
  VOLGEOMfromSurfaImageGeometry(source, &gcam->image);
  VOLGEOMfromSurfaImageGeometry(target, &gcam->atlas);

  // make LTA from matrix
  LTA* lta = LTAalloc(1, nullptr);
  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  lta->type = LINEAR_VOX_TO_VOX;
  if (py::isinstance(affine, surfa.attr("Affine"))) affine = affine.attr("matrix");
  arrayc<double> casted = affine.cast<arrayc<double>>();
  const double *matrix = casted.data();
  for (int i = 0; i < 16; i++) lt->m_L->data[i] = matrix[i];
  VOLGEOMfromSurfaImageGeometry(source, &lt->src);
  VOLGEOMfromSurfaImageGeometry(target, &lt->dst);

  // warp LTA in transform
  TRANSFORM transform;
  transform.type = LINEAR_VOX_TO_VOX;
  transform.xform = (void *)lta;

  // make dataless mri from source geometry and init morph
  MRI *image_header = MRIallocFromVolGeom(&gcam->image, MRI_UCHAR, 1, 1);
  GCAMinit(gcam, image_header, nullptr, &transform, 0);

  // transfer warp field
  double const *dptr = warpdata.data();
  int framesize = gcam->depth * gcam->height * gcam->width;
  for (int zp = 0; zp < gcam->depth; zp++) {
    for (int yp = 0; yp < gcam->height; yp++) {
      for (int xp = 0; xp < gcam->width; xp++) {
        if (*dptr > 0) {
          // assume that any coordinate less than zero should be skipped
          gcam->nodes[xp][yp][zp].x = *(dptr);
          gcam->nodes[xp][yp][zp].y = *(dptr + framesize);
          gcam->nodes[xp][yp][zp].z = *(dptr + framesize * 2);
        }
        dptr++;
      }
    }
  }

  // write file
  GCAMwrite(gcam, filename.c_str());

  // free
  LTAfree(&lta);
  MRIfree(&image_header);
  GCAMfree(&gcam);
}
