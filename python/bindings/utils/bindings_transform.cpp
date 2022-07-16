#include "mrisurf_base.h"
#include "xform.h"

namespace transform {


/*
  Copies a python Geometry object into a VOL_GEOM instance.
*/
void pythonToVolGeom(py::object geometry, VOL_GEOM* vg)
{
  if (geometry.is(py::none())) {
    vg->valid = 0;
  } else {
    std::vector<int> shape = geometry.attr("shape").cast<std::vector<int>>();
    MRI *tmp = new MRI(shape, MRI_UCHAR, false);

    // voxel size
    std::vector<float> voxsize = geometry.attr("voxsize").cast<std::vector<float>>();
    tmp->xsize = voxsize[0];
    tmp->ysize = voxsize[1];
    tmp->zsize = voxsize[2];

    // ensure affine is a c-order double array
    py::object pyaffine = geometry.attr("affine");
    arrayc<double> casted = pyaffine.cast<arrayc<double>>();
    const double *affine = casted.data();
    double xr = affine[0]; double yr = affine[1]; double zr = affine[2];  double pr = affine[3];
    double xa = affine[4]; double ya = affine[5]; double za = affine[6];  double pa = affine[7];
    double xs = affine[8]; double ys = affine[9]; double zs = affine[10]; double ps = affine[11];
    double sx = std::sqrt(xr * xr + xa * xa + xs * xs);
    double sy = std::sqrt(yr * yr + ya * ya + ys * ys);
    double sz = std::sqrt(zr * zr + za * za + zs * zs);
    tmp->x_r = xr / sx;
    tmp->x_a = xa / sx;
    tmp->x_s = xs / sx;
    tmp->y_r = yr / sy;
    tmp->y_a = ya / sy;
    tmp->y_s = ys / sy;
    tmp->z_r = zr / sz;
    tmp->z_a = za / sz;
    tmp->z_s = zs / sz;
    MRIp0ToCRAS(tmp, pr, pa, ps);

    MRIcopyVolGeomFromMRI(tmp, vg);
    MRIfree(&tmp);

    vg->valid = 1;
  }
}


/*
  Creates a python Geometry object from a VOL_GEOM instance.
*/
py::object volGeomToPython(VOL_GEOM* vg)
{
  if (!vg->valid) return py::none();

  py::object shape = py::make_tuple(vg->width, vg->height, vg->depth);
  py::object size = py::make_tuple(vg->xsize, vg->ysize, vg->zsize);

  MRI *tmp = new MRI({vg->width, vg->height, vg->depth}, MRI_UCHAR, false);
  MRIcopyVolGeomToMRI(tmp, vg);

  MATRIX *matrix = extract_i_to_r(tmp);
  py::array affine = copyArray({4, 4}, MemoryOrder::C, matrix->data);
  MatrixFree(&matrix);
  MRIfree(&tmp);

  py::module fsmodule = py::module::import("freesurfer");
  return fsmodule.attr("Geometry")(shape, size, affine);
}


/*
  Converts a python linear transform to an LTA.
*/
LTA* pythonToLTA(py::object transform)
{
  LTA* lta = LTAalloc(1, nullptr);
  LINEAR_TRANSFORM *lt = &lta->xforms[0];

  arrayc<double> casted = transform.attr("matrix").cast<arrayc<double>>();
  const double *affine = casted.data();
  for (int i = 0; i < 16; i++) lt->m_L->data[i] = affine[i];

  py::object type = transform.attr("type");
  if (!type.is(py::none())) lta->type = type.cast<int>();

  pythonToVolGeom(transform.attr("source"), &lt->src);
  pythonToVolGeom(transform.attr("target"), &lt->dst);

  return lta;
}


/*
  Reads a python transform from file via the LTA bridge.
*/
py::object readLTA(const std::string& filename)
{
  LTA *lta = LTAread(filename.c_str());

  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  py::array matrix = copyArray({4, 4}, MemoryOrder::C, lt->m_L->data);

  py::module fsmodule = py::module::import("freesurfer");
  py::object transform = fsmodule.attr("LinearTransform")(matrix);

  transform.attr("type") = lta->type;
  transform.attr("source") = volGeomToPython(&lt->src);
  transform.attr("target") = volGeomToPython(&lt->dst);

  LTAfree(&lta);
  return transform;
}


/*
  Writes a python transform to file via the LTA bridge.
*/
void writeLTA(py::object transform, const std::string& filename)
{
  LTA *lta = pythonToLTA(transform);
  LTAwrite(lta, filename.c_str());
  LTAfree(&lta);
}


}  // end transform namespace
