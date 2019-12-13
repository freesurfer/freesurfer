#include "mrisurf_base.h"
#include "xform.h"


typedef py::array_t<double, py::array::forcecast | py::array::c_style> DoubleArrayC;


py::object volGeomToPython(VOL_GEOM* vg)
{
  if (!vg->valid) return py::none();

  // TODOC
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
  Reads a python transform from file via the LTA bridge.
*/
py::object readLTA(const std::string& filename)
{
  LTA *lta = LTAread(filename.c_str());

  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  py::array matrix = copyArray({4, 4}, MemoryOrder::C, lt->m_L->data);

  py::module fsmodule = py::module::import("freesurfer");
  py::object transform = fsmodule.attr("LinearTransform")(matrix);

  transform.attr("type") = &lt->type;
  transform.attr("source") = volGeomToPython(&lt->src);
  transform.attr("target") = volGeomToPython(&lt->dst);

  LTAfree(&lta);
  return transform;
}
