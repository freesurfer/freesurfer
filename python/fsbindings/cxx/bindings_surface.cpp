#include "mrisurf_base.h"
#include "mrisurf_metricProperties.h"
#include "mrisutils.h"
#include "mrisp.h"

#include "bindings_surface.h"
#include "bindings_transform.h"


/*
  Build a surfa Mesh object from the stored MRIS instance.
*/
py::object MRISBridge::toPython()
{
  // sanity check on the MRIS instance
  if (!p_mris) throw std::runtime_error("cannot bridge to python as MRIS instance is null");

  // extract the vertex and face arrays in order to construct the surface object
  py::array vertices = makeArray({p_mris->nvertices, 3}, MemoryOrder::C, MRISgetVertexArray(p_mris.get()));
  py::array faces = makeArray({p_mris->nfaces, 3}, MemoryOrder::C, MRISgetFaceArray(p_mris.get()));
  py::object surf = py::module::import("surfa").attr("Mesh")(vertices, faces);

  // transfer source volume geometry
  surf.attr("geom") = VOLGEOMtoSurfaImageGeometry(&p_mris->vg);
  surf.attr("space") = "surface";
  return surf;
}


/*
  Return the stored MRIS instance. If one does not exist, it will be created
  from the cached surfa Mesh (assuming it exists).
*/
MRIS* MRISBridge::toMRIS()
{
  // return if the MRI instance has already been set or created
  if (p_mris) return p_mris.get();

  // make sure the source python object has been provided
  if (source.is(py::none())) throw py::value_error("cannot generate MRIS instance without source object");

  // construct the MRIS from vertex and face arrays
  arrayc<float> vertices = source.attr("vertices").cast<arrayc<float>>();
  arrayc<int> faces = source.attr("faces").cast<arrayc<int>>();
  MRIS *mris = MRISfromVerticesAndFaces(vertices.data(), vertices.shape(0), faces.data(), faces.shape(0));

  // transfer vol geometry info
  VOLGEOMfromSurfaImageGeometry(source.attr("geom"), &mris->vg);

  // make sure to register the new MRIS instance in the bridge
  setMRIS(mris);
  return mris;
}


/*
  Read a surfa Mesh using the FS code. Surfa already covers surface IO, so this
  is probably unnecessary, but might be useful at some point.
*/
py::object readSurface(const std::string& filename)
{
  return MRISBridge(MRISread(filename.c_str()));
}


/*
  Write a surfa Mesh using the FS code. Surfa already covers surface IO, so this
  is probably unnecessary, but might be useful at some point.
*/
void writeSurface(MRISBridge surf, const std::string& filename)
{
  MRISwrite(surf, filename.c_str());
}


/*
  Run MRIScomputeSecondFundamentalForm() to compute vertex tangents along the primary
  curvature directions. Returns a tuple of two arrays, for both tangent directions.
*/
py::object computeTangents(MRISBridge surf)
{
  MRIScomputeSecondFundamentalForm(surf);

  // pass along tangent vectors
  MRIS *mris = surf;
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
  return py::make_tuple(tangent1, tangent2);
}


/*
  Compute the surfa Mesh euler number.
*/
int computeEulerNumber(MRISBridge surf)
{
  int unused;
  return MRIScomputeEulerNumber(surf, &unused, &unused, &unused);
}


/*
  Count number of surfa Mesh face intersections.
*/
int countIntersections(MRISBridge surf)
{
  MRIS *mris = surf;
  mrisMarkIntersections(mris);
  int nintersections = 0;
  for(int n = 0; n < mris->nvertices; n++) {
    if (mris->vertices[n].marked) nintersections++;
  }
  return nintersections;
}


/*
  Smooth a surfa Overlay along a Mesh topology.
*/
py::object smoothOverlay(MRISBridge surf, MRIBridge overlay, int steps)
{
  return MRIBridge(MRISsmoothMRIFast(surf, overlay, steps, nullptr, nullptr));
}


/*
  Compute surface distance between two surfa Mesh objects. Returns a surfa Overlay.
*/
py::object surfaceDistance(MRISBridge surf1, MRISBridge surf2)
{
  MRISdistanceBetweenSurfacesExact(surf2, surf1);
  return MRIBridge(MRIcopyMRIS(NULL, surf2, 0, "curv"));
}


/* 
  Inflate surface with spherical unfolding. Returns a surfa Overlay.
*/
py::object quickSphericalInflate(Bridge inSurf, int max_passes, int n_averages, long seed)
{
  // output is updated in the input inSurf
  MRISQuickSphericalInflate(max_passes, n_averages, seed, inSurf, NULL);
  inSurf.updateSource();
  return inSurf;
}


/*
  Parameterize an Overlay to an image Slice. Assumes input surface is a sphere.
  Interp methods can be 'nearest' or 'barycentric'.
*/
py::object parameterize(MRISBridge surf, py::object overlay, int scale, std::string interp)
{
  // get frames and allocate mrisp
  int nframes = overlay.attr("nframes").cast<int>();
  MRI_SP *mrisp = MRISPalloc(scale, nframes);

  // configure projector
  MRIS *mris = surf;
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

  py::array result = makeArray({udim, vdim, nframes}, MemoryOrder::Fortran, buffer);
  return py::module::import("surfa").attr("Slice")(result);
}


/*
  Sample a parameterization into an Overlay. Assumes input surface is a sphere.
  Interp methods can be 'nearest' or 'barycentric'.
*/
py::object sampleParameterization(MRISBridge surf, py::object image, std::string interp)
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
  MRIS *mris = surf;
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

  py::array result = makeArray({mris->nvertices, nframes}, MemoryOrder::Fortran, buffer);
  return py::module::import("surfa").attr("Overlay")(result);
}
