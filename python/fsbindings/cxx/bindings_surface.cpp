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
  mrisMarkIntersections(mris);
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
  Compute surface distance between two surfa Mesh objects. Returns a surfa Overlay.
*/
py::object surfaceDistance(py::object surf1, py::object surf2)
{
  MRIS *mris1 = MRISfromSurfaMesh(surf1);
  MRIS *mris2 = MRISfromSurfaMesh(surf2);
  MRISdistanceBetweenSurfacesExact(mris2, mris1);
  MRISfree(&mris1);
  MRISfree(&mris2);
  return MRItoSurfaArray(MRIcopyMRIS(NULL, mris2, 0, "curv"), true);
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


/*
  Parameterize an Overlay to an image Slice. Assumes input surface is a sphere.
  Interp methods can be 'nearest' or 'barycentric'.
*/
py::object parameterize(py::object surf, py::object overlay, int scale, std::string interp)
{
  // get frames and allocate mrisp
  int nframes = overlay.attr("nframes").cast<int>();
  MRI_SP *mrisp = MRISPalloc(scale, nframes);

  // configure projector
  MRIS *mris = MRISfromSurfaMesh(surf);
  SphericalProjector projector = SphericalProjector(mris, mrisp);
  SphericalProjector::InterpMethod interpmethod;
  if (interp == "nearest") {
    interpmethod = SphericalProjector::Nearest;
  } else if (interp == "barycentric") {
    interpmethod = SphericalProjector::Barycentric;
  } else {
    throw py::value_error("unknown parameterization interpolation method");
  }

  // parameterize the overlay
  const arrayf<float>& arr = overlay.attr("data").cast<arrayf<float>>();
  for (int frame = 0; frame < nframes ; frame++) {
    const float * buff = arr.data() + frame * arr.shape(0);
    projector.parameterizeOverlay(buff, frame, interpmethod);
  }

  // convert MRISP to numpy array
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  float *const buffer = new float[udim * vdim * nframes];
  float *ptr = buffer;
  for (int f = 0; f < nframes ; f++) {
    for (int v = 0; v < vdim; v++) {
      for (int u = 0; u < udim; u++) {
        *ptr++ = *IMAGEFseq_pix(mrisp->Ip, u, v, f);
      }
    }
  }
  MRISPfree(&mrisp);
  MRISfree(&mris);

  py::array result = makeArray({udim, vdim, nframes}, MemoryOrder::Fortran, buffer);
  return py::module::import("surfa").attr("Slice")(result);
}


/*
  Sample a parameterization into an Overlay. Assumes input surface is a sphere.
  Interp methods can be 'nearest' or 'barycentric'.
*/
py::object sampleParameterization(py::object surf, py::object image, std::string interp)
{
  // extract number of frames
  int nframes = image.attr("nframes").cast<int>();
  const arrayf<float>& arr = image.attr("data").cast<arrayf<float>>();

  // init MRISP
  int scale = int(arr.shape(0) / 256);
  MRI_SP *mrisp = MRISPalloc(scale, nframes);
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  if ((arr.shape(0) != udim) || (arr.shape(1) != vdim)) throw py::value_error("parameterization image does not match scale");

  // convert numpy array image to MRISP
  float const *iptr = arr.data();
  for (int f = 0; f < nframes ; f++) {
    for (int v = 0; v < vdim; v++) {
      for (int u = 0; u < udim; u++) {
        *IMAGEFseq_pix(mrisp->Ip, u, v, f) = *iptr++;
      }
    }
  }

  // init spherical projector
  MRIS *mris = MRISfromSurfaMesh(surf);
  SphericalProjector projector = SphericalProjector(mris, mrisp);
  SphericalProjector::InterpMethod interpmethod;
  if (interp == "nearest") {
    interpmethod = SphericalProjector::Nearest;
  } else if (interp == "barycentric") {
    interpmethod = SphericalProjector::Barycentric;
  } else {
    throw py::value_error("unknown parameterization interpolation method");
  }

  // sample MRISP
  float *const buffer = new float[mris->nvertices * nframes];
  float *vptr = buffer;
  for (int frame = 0; frame < nframes ; frame++) {
    projector.sampleParameterization(vptr, frame, interpmethod);
    vptr += mris->nvertices;
  }

  MRISfree(&mris);

  py::array result = makeArray({mris->nvertices, nframes}, MemoryOrder::Fortran, buffer);
  return py::module::import("surfa").attr("Overlay")(result);
}
