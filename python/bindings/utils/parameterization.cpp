#include "surface.h"
#include "mrisp.h"

namespace surf {


/*
  Parameterizes an nvertices-length overlay to an image. Default method is barycentric.
*/
py::array parameterize(Bridge surf, const arrayf<float>& overlay, int scale)
{
  // get frames and allocate mrisp
  int nframes = (overlay.ndim() == 2) ? overlay.shape(1) : 1;
  MRI_SP *mrisp = MRISPalloc(scale, nframes);

  // parameterize
  MRIS *mris = surf.mris();
  BarycentricSphericalProjector projector = BarycentricSphericalProjector(mris, mrisp);
  for (int frame = 0; frame < nframes ; frame++) {
    const float * array = overlay.data() + frame * overlay.shape(0);
    projector.projectOverlay(array, frame);
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


/*
  Samples a parameterization into an nvertices-length overlay. Sampling method is barycentric.
*/
py::array sampleParameterization(Bridge surf, const arrayf<float>& image)
{
  // get number of frames
  int nframes = (image.ndim() == 3) ? image.shape(2) : 1;

  // init MRISP
  int scale = int(image.shape(0) / 256);
  MRI_SP *mrisp = MRISPalloc(scale, nframes);
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  if ((image.shape(0) != udim) || (image.shape(1) != vdim)) throw py::value_error("parameterization image does not match scale");

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
  MRIS *mris = surf.mris();
  float *const buffer = new float[mris->nvertices * nframes];
  float *vptr = buffer;
  for (int frame = 0; frame < nframes ; frame++) {
    MRISfromParameterizationBarycentric(mrisp, mris, frame);
    for (int v = 0 ; v < mris->nvertices ; v++) *vptr++ = mris->vertices[v].curv;
  }
  return makeArray({mris->nvertices, nframes}, MemoryOrder::Fortran, buffer);
}


}  // end namespace surf
