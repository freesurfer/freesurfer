#include "surface.h"
#include "volume.h"
#include "fio.h"


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
