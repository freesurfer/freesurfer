#include "surface.h"
#include "volume.h"
#include "fio.h"
#include "mri_circulars.h"
#include "diag.h"


/*
  Returns vertex positions as an (nvertices, 3) array.
*/
FloatArray PySurface::getVertexPositions()
{
  double* const buffer = new double[m_ptr->nvertices * 3];
  double* ptr = buffer;
  for (int v = 0 ; v < m_ptr->nvertices ; v++) {
    *ptr++ = m_ptr->vertices[v].x;
    *ptr++ = m_ptr->vertices[v].y;
    *ptr++ = m_ptr->vertices[v].z;
  }
  return makeArray({m_ptr->nvertices, 3}, MemoryOrder::C, buffer);
}


/*
  Sets vertex positions from an (nvertices, 3) array.
*/
void PySurface::setVertexPositions(FloatArrayC &positions)
{
  if (positions.request().shape != std::vector<ssize_t>({m_ptr->nvertices, 3})) logFatal(1) << "vertex array shape must be (" << m_ptr->nvertices << ", 3)";
  const float *src = positions.data(0);
  for (int v = 0 ; v < m_ptr->nvertices ; v++) {
    m_ptr->vertices[v].x = *src++;
    m_ptr->vertices[v].y = *src++;
    m_ptr->vertices[v].z = *src++;
  }
}


/*
  Returns vertex normals as an (nvertices, 3) array.
*/
FloatArray PySurface::getVertexNormals()
{
  double* const buffer = new double[m_ptr->nvertices * 3];
  double* ptr = buffer;
  for (int v = 0 ; v < m_ptr->nvertices ; v++) {
    *ptr++ = m_ptr->vertices[v].nx;
    *ptr++ = m_ptr->vertices[v].ny;
    *ptr++ = m_ptr->vertices[v].nz;
  }
  return makeArray({m_ptr->nvertices, 3}, MemoryOrder::C, buffer);
}


/*
  Returns a list of faces per vertex.
*/
std::vector<IntArray> PySurface::getVertexFaces()
{
  std::vector<IntArray> vfaces;
  for (int vno = 0 ; vno < m_ptr->nvertices ; vno++) {
    VERTEX_TOPOLOGY *vt = &m_ptr->vertices_topology[vno];
    int *buffer = new int[vt->num];
    std::copy(vt->f, vt->f + vt->num, buffer);
    vfaces.push_back(makeArray({vt->num}, MemoryOrder::C, buffer));
  }
  return vfaces;
}


/*
  Returns face vertices as an (nfaces, 3) array.
*/
IntArray PySurface::getFaceVertices()
{
  int* const buffer = new int[m_ptr->nfaces * 3];
  int* ptr = buffer;
  for (int f = 0 ; f < m_ptr->nfaces ; f++) {
    *ptr++ = m_ptr->faces[f].v[0];
    *ptr++ = m_ptr->faces[f].v[1];
    *ptr++ = m_ptr->faces[f].v[2];
  }
  return makeArray({m_ptr->nfaces, 3}, MemoryOrder::C, buffer);
}


/*
  Computes the surface's surf->vox matrix
*/
FloatArrayF PySurface::computeSurf2Vox(PyVolume &vol)
{
  MATRIX *sras2ras;
  if (m_ptr->vg.valid) {
    MRI *tmp = MRIallocHeader(m_ptr->vg.width, m_ptr->vg.height, m_ptr->vg.depth, MRI_UCHAR, 1);
    MRIcopyVolGeomToMRI(tmp, &m_ptr->vg);
    sras2ras = RASFromSurfaceRAS_(tmp);
    MRIfree(&tmp);
  } else {
    // no valid geom - assume it came from the provided volume
    sras2ras = RASFromSurfaceRAS_(vol);
  }

  MATRIX *ras2vox = MRIgetRasToVoxelXform(vol);
  MATRIX *sras2vox = MatrixMultiply(ras2vox, sras2ras, NULL);

  FloatArrayF pymat = FloatArrayF({4, 4}, cstrides({4, 4}, sizeof(float)), sras2vox->data);
  MatrixFree(&ras2vox);
  MatrixFree(&sras2vox);
  MatrixFree(&sras2ras);
  return pymat;
}


/*
  Wrapper for IsMRISselfIntersecting. Returns true if surface is self-intersecting.
*/
bool PySurface::isSelfIntersecting()
{
  return IsMRISselfIntersecting(m_ptr);
}


/*
  Wrapper for MRISfillInterior. Returns a binary numpy array.
*/
py::array PySurface::fillInterior()
{
  MRI mri = MRI({256, 256, 256}, MRI_INT);
  MRISfillInterior(m_ptr, 1, &mri);
  return PyVolume::copyImage(&mri);
}


/*
  Reads an annotation file into a numpy array.
*/
IntArray readAnnotation(const std::string& filename)
{
  FILE *fp = fopen(filename.c_str(), "r");
  if (fp == nullptr) logFatal(1) << "could not read annotation file " << filename;

  // first int is the number of elements
  int num = freadInt(fp);

  std::vector<int> vnos, vals;

  // for each one, read in a vno and an int for the annotation value
  for (int i = 0; i < num; i++) {
    vnos.push_back(freadInt(fp));
    vals.push_back(freadInt(fp));
  }

  // get max vno
  ssize_t nvertices = *std::max_element(vnos.begin(), vnos.end()) + 1;
  int *labels = new int[nvertices];

  for (int i = 0; i < num; i++) {
    labels[vnos[i]] = vals[i];
  }

  fclose(fp);
  return py::array_t<int>(nvertices, labels);
}


py::array PySurface::parameterizeBarycentric(const FloatArrayF& overlay)
{
  // ensure overlay matches surfaces
  if (m_ptr->nvertices != overlay.shape(0)) throw py::value_error("overlay does not match surface vertices");

  // parameterization hangs without doing this
  MRISsaveVertexPositions(m_ptr, CANONICAL_VERTICES);

  // get frames and allocate mrisp
  int nframes = (overlay.ndim() == 2) ? overlay.shape(1) : 1;
  MRI_SP *mrisp = MRISPalloc(1, nframes);

  // parameterize
  for (int frame = 0; frame < nframes ; frame++) {
    const float* array = overlay.data() + frame * overlay.shape(0);
    MRIStoParameterizationBarycentric(m_ptr, mrisp, array, 1, frame);
  }

  // convert to numpy array
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
  return makeArray({udim, vdim, nframes}, MemoryOrder::Fortran, buffer);
}


py::array PySurface::sampleParameterization(const FloatArrayF& image)
{
  // get number of frames
  int nframes = (image.ndim() == 3) ? image.shape(2) : 1;

  // init MRISP
  MRI_SP *mrisp = MRISPalloc(1, nframes);
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  if ((image.shape(0) != udim) || (image.shape(1) != vdim)) throw py::value_error("parameterization image must be 256 x 512");

  // copy pixel values from image into MRISP
  float const *iptr = image.data();
  for (int f = 0; f < nframes ; f++) {
    for (int v = 0; v < vdim; v++) {
      for (int u = 0; u < udim; u++) {
        *IMAGEFseq_pix(mrisp->Ip, u, v, f) = *iptr++;
      }
    }
  }

  // sample and save results into numpy array
  float *const buffer = new float[m_ptr->nvertices * nframes];
  float *vptr = buffer;
  for (int frame = 0; frame < nframes ; frame++) {
    MRISfromParameterizationBarycentric(mrisp, m_ptr, frame);
    for (int v = 0 ; v < m_ptr->nvertices ; v++) *vptr++ = m_ptr->vertices[v].curv;
  }
  return makeArray({m_ptr->nvertices, nframes}, MemoryOrder::Fortran, buffer);
}


FloatArrayC PySurface::computeParameterizationMapBarycentric()
{
  float *buffer = new float[m_ptr->nvertices * 2];
  float *ptr = buffer;

  int udim = 256;
  int vdim = 512;

  float radius = MRISaverageRadius(m_ptr);
  for (int vno = 0; vno < m_ptr->nvertices; vno++) {

    VERTEX *vertex = &m_ptr->vertices[vno];
    float x = vertex->x;
    float y = vertex->y;
    float z = vertex->z;

    float theta = atan2(y/radius, x/radius);
    if (theta < 0.0f) theta = 2 * M_PI + theta;  // make it 0 -> 2PI

    float d = radius * radius - z * z;
    if (d < 0.0) d = 0.0;

    float phi = atan2(sqrt(d), z);
    if (phi < RADIANS(1)) DiagBreak();
    if (phi > M_PI) DiagBreak();

    float uf = udim * phi / PHI_MAX;
    float vf = vdim * theta / THETA_MAX;

    if (uf < 0) uf = -uf;
    if (uf >= udim) uf = udim - (uf - udim + 1);
    if (vf < 0) vf += vdim;
    if (vf >= vdim) vf -= vdim;

    *ptr++ = uf;
    *ptr++ = vf;
  }

  return makeArray({m_ptr->nvertices, 2}, MemoryOrder::C, buffer);
}
