#include "surface.h"
#include "mrisp.h"

namespace surf {


/*
  Parameterizes an nvertices-length overlay to an image. Interp methods can be 'nearest' or
  'barycentric'.
*/
py::array parameterize(Bridge surf, const arrayf<float>& overlay, int scale, std::string interp)
{
  // get frames and allocate mrisp
  int nframes = (overlay.ndim() == 2) ? overlay.shape(1) : 1;
  MRI_SP *mrisp = MRISPalloc(scale, nframes);

  // configure projector
  MRIS *mris = surf.mris();
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
  for (int frame = 0; frame < nframes ; frame++) {
    const float * array = overlay.data() + frame * overlay.shape(0);
    projector.parameterizeOverlay(array, frame, interpmethod);
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

  return makeArray({udim, vdim, nframes}, MemoryOrder::Fortran, buffer);
}


/*
  Samples a parameterization into an nvertices-length overlay. Interp methods can be 'nearest' or
  'barycentric'.
*/
py::array sampleParameterization(Bridge surf, const arrayf<float>& image, std::string interp)
{
  // extract number of frames
  int nframes = (image.ndim() == 3) ? image.shape(2) : 1;

  // init MRISP
  int scale = int(image.shape(0) / 256);
  MRI_SP *mrisp = MRISPalloc(scale, nframes);
  int udim = U_DIM(mrisp);
  int vdim = V_DIM(mrisp);
  if ((image.shape(0) != udim) || (image.shape(1) != vdim)) throw py::value_error("parameterization image does not match scale");

  // convert numpy array image to MRISP
  float const *iptr = image.data();
  for (int f = 0; f < nframes ; f++) {
    for (int v = 0; v < vdim; v++) {
      for (int u = 0; u < udim; u++) {
        *IMAGEFseq_pix(mrisp->Ip, u, v, f) = *iptr++;
      }
    }
  }

  // init spherical projector
  MRIS *mris = surf.mris();
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

  return makeArray({mris->nvertices, nframes}, MemoryOrder::Fortran, buffer);
}


}  // end namespace surf
