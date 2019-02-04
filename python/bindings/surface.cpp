#include "surface.h"
#include "volume.h"


/*
  Surface-specific pybind module configuration.
*/
void bindSurface(py::module &m) {
  // CoreSurface class
  py::class_<CoreSurface>(m, "CoreSurface")
    .def(py::init<const std::string &>())
    .def("read", &CoreSurface::read)
    .def("write", &CoreSurface::write)
    .def("isSelfIntersecting", &CoreSurface::isSelfIntersecting)
    .def("fillInterior", &CoreSurface::fillInterior)
    .def_property("vertices", &CoreSurface::getVertices, &CoreSurface::setVertices)
    ;
}


/*
  Constructs a surface from a surface file.
*/
CoreSurface::CoreSurface(const std::string &filename) {
  read(filename);
}


/*
  Frees the underlying MRIS before destructing.
*/
CoreSurface::~CoreSurface() {
  if (m_mris) MRISfree(&m_mris);
}


/*
  Sets the underlying MRIS instance.
*/
void CoreSurface::setMRIS(MRIS *mris) {
  // update the underlying mri pointer
  if (m_mris) MRISfree(&m_mris);
  m_mris = mris;
}


/*
  Loads the underlying MRIS structure from a surface file.
*/
void CoreSurface::read(const std::string &filename) {
  MRIS *mris = MRISread(filename.c_str());
  if (!mris) throw py::value_error("could not load surface " + filename);
  setMRIS(mris);
}


/*
  Writes the underlying MRIS structure to a file.
*/
void CoreSurface::write(const std::string &filename) {
  MRISwrite(m_mris, filename.c_str());
}


/*
  Returns the underlying vertex data as an N x 3 numpy array.
*/
py::array_t<float> CoreSurface::getVertices() {
  double* const buffer = new double[m_mris->nvertices * 3];
  double* ptr = buffer;
  for (int v = 0 ; v < m_mris->nvertices ; v++) {
    *ptr++ = m_mris->vertices[v].x;
    *ptr++ = m_mris->vertices[v].y;
    *ptr++ = m_mris->vertices[v].z;
  }
  return makeArray({m_mris->nvertices, 3}, MemoryOrder::C, buffer);
}


/*
  Sets the underlying MRIS vertex data from an N x 3 numpy array, where N must match
  the surface's number of vertices.
*/
void CoreSurface::setVertices(py::array_t<float, py::array::c_style | py::array::forcecast> array) {
  // get buffer info
  py::buffer_info info = array.request();

  // sanity check on input dimensions
  if ((info.ndim != 2) ||(info.shape[0] != m_mris->nvertices) || (info.shape[1] != 3)) {
    throw py::value_error("vertex array must of shape " + shapeString({m_mris->nvertices, 3}));
  }

  // set vertex data
  const float *ptr = array.data(0);
  for (int v = 0 ; v < m_mris->nvertices ; v++) {
    VERTEX *vertex = &m_mris->vertices[v];
    vertex->x = *ptr++;
    vertex->y = *ptr++;
    vertex->z = *ptr++;
  }
}


/*
  Wrapper for IsMRISselfIntersecting. Returns true if surface is self-intersecting.
*/
bool CoreSurface::isSelfIntersecting() {
  return IsMRISselfIntersecting(m_mris);
}


/*
  Wrapper for MRISfillInterior. Returns a binary numpy array.
*/
py::array CoreSurface::fillInterior() {
  // allocate the destination volume
  MRI *mri = MRIallocChunk(256, 256, 256, MRI_INT, 1);
  MRISfillInterior(m_mris, 1, mri);
  // copy the mri buffer into a numpy array
  py::array array = makeArray(mri);
  // remember to free the mri before returning
  MRIfree(&mri);
  return array;
}
