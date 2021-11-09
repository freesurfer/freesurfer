#include "mrisurf_base.h"
#include "surface.h"
#include "xform.h"
#include "mrisurf_metricProperties.h"


namespace surf {


/*
  Builds a python Surface object from the internal MRIS instance.
*/
py::object Bridge::python()
{
  // sanity check on the MRIS instance
  if (!p_mris) throw std::runtime_error("cannot bridge to python as MRIS instance is null");

  // extract the vertex and face arrays in order to construct the surface object
  py::array vertices = makeArray({p_mris->nvertices, 3}, MemoryOrder::C, MRISgetVertexArray(p_mris.get()));
  py::array faces = makeArray({p_mris->nfaces, 3}, MemoryOrder::C, MRISgetFaceArray(p_mris.get()));
  py::module fsmodule = py::module::import("freesurfer");
  py::object surf = fsmodule.attr("Surface")(vertices, faces);

  // transfer the rest of the parameters
  transferParameters(surf);
  return surf;
}


/*
  Transfers parameters (not vertices and faces) from the internal MRIS
  to the provided surface instance.
*/
void Bridge::transferParameters(py::object& surf)
{
  // pass along face and vertex normals
  surf.attr("vertex_normals") = makeArray({p_mris->nvertices, 3}, MemoryOrder::C, MRISgetVertexNormalArray(p_mris.get()));
  surf.attr("face_normals") = makeArray({p_mris->nfaces, 3}, MemoryOrder::C, MRISgetFaceNormalArray(p_mris.get()));

  // pass along hemi info
  switch (p_mris->hemisphere) {
    case LEFT_HEMISPHERE:  surf.attr("hemi") = "left";  break;
    case RIGHT_HEMISPHERE: surf.attr("hemi") = "right"; break;
    case BOTH_HEMISPHERES: surf.attr("hemi") = "both";  break;
    default:               surf.attr("hemi") = py::none(); break;
  }

  // transfer source volume geometry
  surf.attr("geom") = transform::volGeomToPython(&p_mris->vg);
}


/*
  Updates the cached python surface object with the internal MRIS instance.
*/
void Bridge::updateSource()
{
  // sanity check on the MRIS instance and python source
  if (!p_mris) throw std::runtime_error("cannot bridge to python as MRIS instance is null");
  if (source.is(py::none())) throw py::value_error("cannot update source if it does not exist");

  // update the extracted vertex and face arrays and let transferParameters() do the rest
  source.attr("vertices") = makeArray({p_mris->nvertices, 3}, MemoryOrder::C, MRISgetVertexArray(p_mris.get()));
  source.attr("faces") = makeArray({p_mris->nfaces, 3}, MemoryOrder::C, MRISgetFaceArray(p_mris.get()));
  transferParameters(source);
}


/*
  Returns the internal MRIS instance. If one does not exist, it will be created
  from the cached python source surface.
*/
MRIS* Bridge::mris()
{
  // return if the MRI instance has already been set or created
  if (p_mris) return p_mris.get();

  // make sure the source python object has been provided
  if (source.is(py::none())) throw py::value_error("cannot generate MRIS instance without source object");

  // construct the MRIS from vertex and face arrays
  arrayc<float> vertices = source.attr("vertices").cast<arrayc<float>>();
  arrayc<int> faces = source.attr("faces").cast<arrayc<int>>();
  MRIS *mris = MRISfromVerticesAndFaces(vertices.data(), vertices.shape(0), faces.data(), faces.shape(0));

  // transfer hemisphere info
  mris->hemisphere = NO_HEMISPHERE;
  py::object pyhemi = source.attr("hemi");
  if (!pyhemi.is(py::none())) {
    std::string hemi = pyhemi.cast<std::string>();
    if      (hemi == "left")  { mris->hemisphere = LEFT_HEMISPHERE; }
    else if (hemi == "right") { mris->hemisphere = RIGHT_HEMISPHERE; }
    else if (hemi == "both")  { mris->hemisphere = BOTH_HEMISPHERES; }
  }

  // transfer vol geometry info
  transform::pythonToVolGeom(source.attr("geom"), &mris->vg);

  // make sure to register the new MRIS instance in the bridge
  setmris(mris);
  return mris;
}


/*
  Reads a python surface from file via the MRIS bridge.
*/
py::object read(const std::string& filename)
{
  return Bridge(MRISread(filename.c_str()));
}


/*
  Writes a python surface to file via the MRIS bridge.
*/
void write(Bridge surf, const std::string& filename)
{
  MRISwrite(surf, filename.c_str());
}


/*
  Runs MRIScomputeMetricProperties() to compute face and vertex normals.
*/
void computeNormals(Bridge surf)
{
  MRIScomputeMetricProperties(surf);
  surf.updateSource();
}


/*
  Runs MRIScomputeSecondFundamentalForm() to compute vertex tangents along the
  primary curvature directions. No need to call surf.updateSource() at the
  end of this function since we're only updating tangents, which are parameters
  that might not be defined in the surface (so they shouldn't be handled in
  the transferParameters() function anyway).
*/
void computeTangents(Bridge surf)
{
  MRIScomputeSecondFundamentalForm(surf);

  // pass along tangent vectors
  MRIS *mris = surf.mris();
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
  surf.source.attr("vertex_tangents_1") = makeArray({mris->nvertices, 3}, MemoryOrder::C, e1_buffer);
  surf.source.attr("vertex_tangents_2") = makeArray({mris->nvertices, 3}, MemoryOrder::C, e2_buffer);
}


/*
  Computes the mesh euler number.
*/
int computeEulerNumber(Bridge surf)
{
  int unused;
  return MRIScomputeEulerNumber(surf, &unused, &unused, &unused);
}


/*
  Count number of face intersections.
*/
int countIntersections(Bridge surf)
{
  MRIS *mris = surf.mris();
  mrisMarkIntersections(mris);
  int nintersections = 0;
  for(int n = 0; n < mris->nvertices; n++) {
    if (mris->vertices[n].marked) nintersections++;
  }
  return nintersections;
}


/*
  Smoothes an overlay along mesh vertices.
*/
py::object smoothOverlay(Bridge surf, vol::Bridge overlay, int steps)
{
  return vol::Bridge(MRISsmoothMRIFast(surf, overlay, steps, nullptr, nullptr));
}


/*
  Smoothes an overlay along mesh vertices.
*/
py::object surfaceDistance(Bridge surf1, Bridge surf2)
{
  MRISdistanceBetweenSurfacesExact(surf2, surf1);
  return vol::Bridge(MRIcopyMRIS(NULL, surf2, 0, "curv"));
}


/*

namespace vec {

  static void make(const float* a, const float* b, float* ab) {
    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];
  }

  static void unit(float* x) {
    float norm = std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
    x[0] /= norm;
    x[1] /= norm;
    x[2] /= norm;
  }

  static void cross(const float* a, const float* b, float* x) {
    x[0] = - b[1] * a[2] + a[1] * b[2];
    x[1] =   b[0] * a[2] - a[0] * b[2];
    x[2] = - b[0] * a[1] + a[0] * b[1];
  }

  static void addto(const float* a, float* b) {
    b[0] += a[0];
    b[1] += a[1];
    b[2] += a[2];
  }

  static void copyto(const float* a, float* b) {
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
  }

  static void print(const float* x) {
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl << std::flush;
  }
}


void computeNormals(py::object surf)
{
  const float* vertices = surf.attr("vertices").cast<py::array_t<float>>().data();
  int nvertices = surf.attr("nvertices").cast<int>();
  int vsize = nvertices * 3;

  const int* faces = surf.attr("faces").cast<py::array_t<int>>().data();
  int nfaces = surf.attr("nfaces").cast<int>();
  int fsize = nfaces * 3;

  float* vertexnorms = new float[vsize]();
  float* facenorms = new float[fsize];

  for (int f = 0 ; f < fsize ; f += 3) {
    int va = 3 * faces[f];
    int vb = 3 * faces[f + 1];
    int vc = 3 * faces[f + 2];

    const float* a_pos = &vertices[va];
    const float* b_pos = &vertices[vb];
    const float* c_pos = &vertices[vc];

    float ab[3], ac[3], bc[3], cross[3];

    vec::make(a_pos, b_pos, ab);
    vec::make(a_pos, c_pos, ac);
    vec::make(b_pos, c_pos, bc);

    vec::unit(ab);
    vec::unit(ac);
    vec::unit(bc);

    // vertex a
    vec::cross(ab, ac, cross);
    vec::addto(cross, &vertexnorms[va]);

    // face normal
    vec::unit(cross);
    vec::copyto(cross, &facenorms[f]);

    // vertex b
    vec::cross(ab, bc, cross);
    vec::addto(cross, &vertexnorms[vb]);

    // vertex c
    vec::cross(ac, bc, cross);
    vec::addto(cross, &vertexnorms[vc]);
  }

  for (int n = 0 ; n < vsize ; n += 3) vec::unit(&vertexnorms[n]);

  surf.attr("vertex_normals") = makeArray({nvertices, 3}, MemoryOrder::C, vertexnorms);
  surf.attr("face_normals") = makeArray({nfaces, 3}, MemoryOrder::C, facenorms);
}


py::object read_directly(const std::string& filename)
{
  FILE* fp = fopen(filename.c_str(), "rb");
  if (!fp) ErrorExit(ERROR_BADFILE, "could not open file");

  int magic;
  fread3(&magic, fp);

  char line[STRLEN];
  fgets(line, 200, fp);
  fscanf(fp, "\n");

  int const nvertices = freadInt(fp);
  int const nfaces = freadInt(fp);

  if (nvertices < 0) ErrorExit(ERROR_BADFILE, "freadInt returned negative num vertices: %d!", nvertices);
  if (nfaces < 0) ErrorExit(ERROR_BADFILE, "freadInt returned negative num faces: %d!", nfaces);

  int vsize = nvertices * 3;
  float* vertices = new float[vsize];
  for (int n = 0; n < vsize; n++) vertices[n] = freadFloat(fp);

  int fsize = nfaces * 3;
  int* faces = new int[fsize];
  for (int n = 0; n < fsize; n++) faces[n] = freadInt(fp);

  fclose(fp);

  py::array varray = makeArray({nvertices, 3}, MemoryOrder::C, vertices);
  py::array farray = makeArray({nfaces, 3}, MemoryOrder::C, faces);
  py::module fsmodule = py::module::import("freesurfer");
  py::object surf = fsmodule.attr("Surface")(varray, farray);

  return surf;
}

*/

}  // end namespace surf
