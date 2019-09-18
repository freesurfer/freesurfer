#ifndef SURFACE_H
#define SURFACE_H

#include "numpy.h"
#include "mrisurf.h"
#include "volume.h"


typedef py::array_t<float, py::array::f_style | py::array::forcecast> affinematrix;

typedef py::array_t<int,   py::array::forcecast> IntArray;
typedef py::array_t<float, py::array::forcecast> FloatArray;
typedef py::array_t<float, py::array::forcecast | py::array::c_style> FloatArrayC;
typedef py::array_t<float, py::array::forcecast | py::array::f_style> FloatArrayF;


/*

*/
class PySurface
{
public:
  ~PySurface() { MRISfree(&m_ptr); }

  // MRIS conversions
  PySurface(MRIS *surf) { m_ptr = surf; }
  operator MRIS*() { return m_ptr; }
  MRIS* mris() { return m_ptr; }

  // filesystem IO
  PySurface(const std::string& filename)  { m_ptr = MRISread(filename.c_str()); }
  void write(const std::string& filename) { MRISwrite(m_ptr, filename.c_str()); }
  PySurface* copy() { return new PySurface(MRISclone(m_ptr)); }
  void status() { std::cout << m_ptr->status << std::endl; }

  // vertex information
  FloatArray getVertexPositions();
  void setVertexPositions(FloatArrayC& positions);
  FloatArray getVertexNormals();
  std::vector<IntArray> getVertexFaces();

  // face information
  IntArray getFaceVertices();

  // geometry
  void copyGeometry(PySurface& surf) { copyVolGeom(&surf.mris()->vg, &m_ptr->vg); }
  void copyGeometry(PyVolume& vol) { MRIScopyVolGeomFromMRI(m_ptr, vol); }
  FloatArrayF computeSurf2Vox(PyVolume& vol);

  // parameterization
  FloatArrayC parameterizeBarycentric(const FloatArrayF& array);
  FloatArrayC sampleParameterization(const FloatArrayC& param);
  FloatArrayC computeParameterizationMapBarycentric();

  // misc utilities
  bool isSelfIntersecting();
  py::array fillInterior();

  // deprecations - TO BE REMOVED 6/20
  FloatArray getVerticesDeprecated() {
    logWarning << "The Surface vertices property has been replaced by get_vertex_positions() and set_vertex_positions(). It will be completely removed soon."; return getVertexPositions(); }
  void setVerticesDeprecated(FloatArrayC &array) {
    logWarning << "The Surface vertices property has been replaced by get_vertex_positions() and set_vertex_positions(). It will be completely removed soon."; setVertexPositions(array); }
  FloatArray getFacesDeprecated() {
    logWarning << "The Surface faces property has been replaced by get_face_vertices(). It will be completely removed soon."; return getFaceVertices(); }

private:
  MRIS* m_ptr;
};


IntArray readAnnotation(const std::string &filename);


inline void bindSurface(py::module &m)
{
  // PySurface class
  py::class_<PySurface> surf(m, "Surface");
  // IO
  surf.def(py::init<const std::string &>());
  surf.def("write", &PySurface::write);
  // header info
  surf.def("copy", &PySurface::copy);
  surf.def("copy_geometry", (void (PySurface::*)(PySurface&)) &PySurface::copyGeometry);
  surf.def("copy_geometry", (void (PySurface::*)(PyVolume&)) &PySurface::copyGeometry);
  surf.def("_compute_surf2vox", &PySurface::computeSurf2Vox);
  // vertices
  surf.def("get_vertex_positions", &PySurface::getVertexPositions);
  surf.def("set_vertex_positions", &PySurface::setVertexPositions);
  surf.def("get_vertex_normals", &PySurface::getVertexNormals);
  surf.def("get_vertex_faces", &PySurface::getVertexFaces);
  surf.def_property("vertices", &PySurface::getVerticesDeprecated, &PySurface::setVerticesDeprecated);
  // faces
  surf.def("get_face_vertices", &PySurface::getFaceVertices);
  surf.def_property_readonly("faces", &PySurface::getFacesDeprecated);
  // parameterization
  surf.def("parameterize", &PySurface::parameterizeBarycentric);
  surf.def("sample_parameterization", &PySurface::sampleParameterization);
  surf.def("parameterization_map", &PySurface::computeParameterizationMapBarycentric);
  // misc utilities
  surf.def("isSelfIntersecting", &PySurface::isSelfIntersecting);
  surf.def("fillInterior", &PySurface::fillInterior);

  m.def("read_annotation", &readAnnotation);
}

#endif
