#include "mrisurf_base.h"
#include "mrisurf_metricProperties.h"
#include "mrisutils.h"
#include "mrisp.h"

#include "bindings_surface.h"
#include "bindings_transform.h"



/*
  Convert an MRIS structure to a surfa Mesh. All data is copied between python and cxx, so any
  allocated MRIS pointers will need to be freed, even after convert to python. However, if
  `release` is `true`, this function will free the MRIS after converting.
*/
py::object MRIStoSurfaMesh(MRIS *mris, bool release)
{
  if (mris == nullptr) throw py::value_error("MRIStoSurfaMesh: cannot convert to surfa Mesh - MRIS input is null");

  // extract the vertex and face arrays in order to construct the surface object
  py::array vertices = makeArray({mris->nvertices, 3}, MemoryOrder::C, MRISgetVertexArray(mris));
  py::array faces = makeArray({mris->nfaces, 3}, MemoryOrder::C, MRISgetFaceArray(mris));
  py::object surface = py::module::import("surfa").attr("Mesh")(vertices, faces);

  // transfer source volume geometry
  surface.attr("geom") = VOLGEOMtoSurfaImageGeometry(&mris->vg);
  surface.attr("space") = "surface";

  if (release) MRISfree(&mris);

  return surface;
}


/*
  Convert a surfa Mesh to an MRIS structure. All data is copied between python and cxx, so the
  returned MRIS pointer will need to be freed manually once it's done with.
*/
MRIS* MRISfromSurfaMesh(py::object mesh)
{
  // type checking
  py::object meshclass = py::module::import("surfa").attr("Mesh");
  if (!py::isinstance(mesh, meshclass)) throw py::value_error("MRISfromSurfaMesh: cannot convert to MRIS - input is not a surfa Mesh");

  // construct the MRIS from vertex and face arrays
  arrayc<float> vertices = mesh.attr("vertices").cast<arrayc<float>>();
  arrayc<int> faces = mesh.attr("faces").cast<arrayc<int>>();
  MRIS *mris = MRISfromVerticesAndFaces(vertices.data(), vertices.shape(0), faces.data(), faces.shape(0));

  // transfer vol geometry info
  VOLGEOMfromSurfaImageGeometry(mesh.attr("geom"), &mris->vg);

  return mris;
}


/*
  Read a surfa Mesh using the FS code. Surfa already covers surface IO, so this
  is probably unnecessary, but might be useful at some point.
*/
py::object readSurface(const std::string& filename)
{
  return MRIStoSurfaMesh(MRISread(filename.c_str()), true);
}


/*
  Write a surfa Mesh using the FS code. Surfa already covers surface IO, so this
  is probably unnecessary, but might be useful at some point.
*/
void writeSurface(py::object surf, const std::string& filename)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  MRISwrite(mris, filename.c_str());
  MRISfree(&mris);
}


/*
  Run MRIScomputeSecondFundamentalForm() to compute vertex tangents along the primary
  curvature directions. Returns a tuple of two arrays, for both tangent directions.
*/
py::object computeTangents(py::object surf)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  MRIScomputeSecondFundamentalForm(mris);

  // pass along tangent vectors
  float * const e1_buffer = new float[mris->nvertices * 3];
  float * const e2_buffer = new float[mris->nvertices * 3];
  float * e1_ptr = e1_buffer;
  float * e2_ptr = e2_buffer;
  for (int n = 0 ; n < mris->nvertices ; n++) {
    VERTEX *v = &mris->vertices[n];
    *e1_ptr++ = v->e1x;
    *e1_ptr++ = v->e1y;
    *e1_ptr++ = v->e1z;
    *e2_ptr++ = v->e2x;
    *e2_ptr++ = v->e2y;
    *e2_ptr++ = v->e2z;
  }
  py::object tangent1 = makeArray({mris->nvertices, 3}, MemoryOrder::C, e1_buffer);
  py::object tangent2 = makeArray({mris->nvertices, 3}, MemoryOrder::C, e2_buffer);
  MRISfree(&mris);
  return py::make_tuple(tangent1, tangent2);
}


/*
  Compute the surfa Mesh euler number.
*/
int computeEulerNumber(py::object surf)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  int unused;
  int e = MRIScomputeEulerNumber(mris, &unused, &unused, &unused);
  MRISfree(&mris);
  return e;
}


/*
  Count number of surfa Mesh face intersections.
*/
int countIntersections(py::object surf)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  mrisMarkIntersections(mris,0);
  int nintersections = 0;
  for(int n = 0; n < mris->nvertices; n++) {
    if (mris->vertices[n].marked) nintersections++;
  }
  MRISfree(&mris);
  return nintersections;
}


/*
  Smooth a surfa Overlay along a Mesh topology.
*/
py::object smoothOverlay(py::object surf, py::object overlay, int steps)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  MRI *mri_overlay = MRIfromSurfaArray(overlay);
  MRI *mri_smoothed = MRISsmoothMRIFast(mris, mri_overlay, steps, nullptr, nullptr);
  MRISfree(&mris);
  MRIfree(&mri_overlay);
  return MRItoSurfaArray(mri_smoothed, true);
}


/* 
  Inflate surface with spherical unfolding. Returns a surfa Overlay.
*/
py::object quickSphericalInflate(py::object surf, int max_passes, int n_averages, long seed)
{
  MRIS *mris = MRISfromSurfaMesh(surf);
  MRISQuickSphericalInflate(max_passes, n_averages, seed, mris, NULL);
  return MRIStoSurfaMesh(mris, true);
}
