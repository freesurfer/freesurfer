#include "surface.h"
#include "volume.h"
#include "fio.h"
#include "mri_circulars.h"
#include "diag.h"

// temporary utility for reading via constructor so that we don't
// have to make a base-class MRIS constructor (for now)
static void readIntoMRISPointer(MRIS *mris, const std::string &filename)
{
  *mris = *MRISread(filename.c_str());
}


/*
  Constructs a surface from a surface file.
*/
PySurface::PySurface(const std::string &filename)
{
  readIntoMRISPointer(this, filename);
}


/*
  Writes the underlying MRIS structure to a file.
*/
void PySurface::write(const std::string &filename)
{
  MRISwrite(this, filename.c_str());
}


/*
  Returns the underlying vertex data as an N x 3 numpy array.
*/
py::array_t<float> PySurface::getVertices()
{
  double* const buffer = new double[nvertices * 3];
  double* ptr = buffer;
  for (int v = 0 ; v < nvertices ; v++) {
    *ptr++ = vertices[v].x;
    *ptr++ = vertices[v].y;
    *ptr++ = vertices[v].z;
  }
  return makeArray({nvertices, 3}, MemoryOrder::C, buffer);
}


/*
  Sets the underlying vertex data from an N x 3 numpy array.
*/
void PySurface::setVertices(py::array_t<float, py::array::c_style | py::array::forcecast> array)
{
  if (array.request().shape != std::vector<ssize_t>({nvertices, 3})) logFatal(1) << "vertex array shape must be (" << nvertices << ", 3)";
  const float *src = array.data(0);
  for (int v = 0 ; v < nvertices ; v++) {
    vertices[v].x = *src++;
    vertices[v].y = *src++;
    vertices[v].z = *src++;
  }
}


/*
  Returns the underlying face data as an N x 3 numpy array.
*/
py::array_t<int> PySurface::getFaces()
{
  int* const buffer = new int[nfaces * 3];
  int* ptr = buffer;
  for (int f = 0 ; f < nfaces ; f++) {
    *ptr++ = faces[f].v[0];
    *ptr++ = faces[f].v[1];
    *ptr++ = faces[f].v[2];
  }
  return makeArray({nfaces, 3}, MemoryOrder::C, buffer);
}


/*
  Sets the underlying face data from an N x 3 numpy array.
*/
void PySurface::setFaces(py::array_t<int, py::array::c_style | py::array::forcecast> array)
{
  if (array.request().shape != std::vector<ssize_t>({nfaces, 3})) logFatal(1) << "vertex array shape must be (" << nfaces << ", 3)";
  const int *src = array.data(0);
  for (int f = 0 ; f < nfaces ; f++) {
    faces[f].v[0] = *src++;
    faces[f].v[1] = *src++;
    faces[f].v[2] = *src++;
  }
}


/*
  Computes the surface's surf->vox matrix
*/
affinematrix PySurface::computeSurf2Vox(PyVolume& vol)
{
  MATRIX *sras2ras;
  if (vg.valid) {
    MRI *tmp = MRIallocHeader(vg.width, vg.height, vg.depth, MRI_UCHAR, 1);
    MRIcopyVolGeomToMRI(tmp, &vg);
    sras2ras = RASFromSurfaceRAS_(tmp);
    MRIfree(&tmp);
  } else {
    // no valid geom - assume it came from the provided volume
    sras2ras = RASFromSurfaceRAS_(vol);
  }

  MATRIX *ras2vox = MRIgetRasToVoxelXform(vol);
  MATRIX *sras2vox = MatrixMultiply(ras2vox, sras2ras, NULL);

  affinematrix pymat = affinematrix({4, 4}, cstrides({4, 4}, sizeof(float)), sras2vox->data);
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
  return IsMRISselfIntersecting(this);
}


/*
  Wrapper for MRISfillInterior. Returns a binary numpy array.
*/
py::array PySurface::fillInterior()
{
  MRI mri = MRI({256, 256, 256}, MRI_INT);
  MRISfillInterior(this, 1, &mri);
  return PyVolume::copyImage(&mri);
}


/*
  Reads an annotation file into a numpy array.
*/
py::array_t<int> readAnnotation(const std::string& filename)
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


py::array_t<float> PySurface::parameterizeBarycentric(const py::array_t<float, py::array::forcecast>& overlay, float scale)
{
  MRISsaveVertexPositions(this, CANONICAL_VERTICES);
  MRI_SP *mrisp = MRIStoParameterizationBarycentric(this, nullptr, overlay.data(0), scale, 0);
  // convert to array
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  float *const buffer = new float[udim * vdim];
  float *ptr = buffer;
  for (int u = 0; u < udim; u++) {
  for (int v = 0; v < vdim; v++) {
      *ptr++ = *IMAGEFseq_pix(mrisp->Ip, u, v, 0);
  }
  }
  MRISPfree(&mrisp);
  return makeArray({udim, vdim}, MemoryOrder::C, buffer);
}


py::array_t<float> PySurface::computeParameterizationMapBarycentric(const std::vector<int>& shape)
{
  float *buffer = new float[nvertices * 4 * 3];
  float *ptr = buffer;

  int udim = shape[0];
  int vdim = shape[1];

  float radius = MRISaverageRadius(this);
  for (int vno = 0; vno < nvertices; vno++) {

    VERTEX *vertex = &vertices[vno];
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
    int u0 = floor(uf);
    int v0 = floor(vf);
    int u1 = ceil(uf);
    int v1 = ceil(vf);
    float du = uf - (float)u0;
    float dv = vf - (float)v0;

    // enforce spherical topology
    if (u0 < 0) u0 = -u0;
    if (u0 >= udim) u0 = udim - (u0 - udim + 1);
    if (u1 < 0) u1 = -u1;
    if (u1 >= udim) u1 = udim - (u1 - udim + 1);
    if (v0 < 0) v0 += vdim;
    if (v0 >= vdim) v0 -= vdim;
    if (v1 < 0) v1 += vdim;
    if (v1 >= vdim) v1 -= vdim;

    // make bilinear interpolation map
    *ptr++ = du * dv                   ; *ptr++ = u1; *ptr++ = v1;
    *ptr++ = (1.0f - du) * dv          ; *ptr++ = u0; *ptr++ = v1;
    *ptr++ = (1.0f - du) * (1.0f - dv) ; *ptr++ = u0; *ptr++ = v0;
    *ptr++ = du * (1.0f - dv)          ; *ptr++ = u1; *ptr++ = v0;
  }

  return makeArray({nvertices, 4, 3}, MemoryOrder::C, buffer);
}
