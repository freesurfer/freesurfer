#include "morph.h"
#include "xform.h"
#include "gcamorph.h"


namespace morph {


void writeGCAM(py::object warp, const std::string& filename)
{
  // get warp map and alloc morph
  arrayf<double> data = warp.attr("data").cast<arrayf<double>>();
  GCAM *gcam = GCAMalloc(data.shape(0), data.shape(1), data.shape(2));

  // copy geometries
  transform::pythonToVolGeom(warp.attr("source"), &gcam->image);
  transform::pythonToVolGeom(warp.attr("target"), &gcam->atlas);

  // make LTA from matrix
  LTA* lta = LTAalloc(1, nullptr);
  LINEAR_TRANSFORM *lt = &lta->xforms[0];
  lta->type = LINEAR_VOX_TO_VOX;
  arrayc<double> casted = warp.attr("affine").cast<arrayc<double>>();
  const double *affine = casted.data();
  for (int i = 0; i < 16; i++) lt->m_L->data[i] = affine[i];
  transform::pythonToVolGeom(warp.attr("source"), &lt->src);
  transform::pythonToVolGeom(warp.attr("target"), &lt->dst);

  // warp LTA in transform
  TRANSFORM transform;
  transform.type = LINEAR_VOX_TO_VOX;
  transform.xform = (void *)lta;

  // make dataless mri from source geometry and init morph
  MRI *image_header = MRIallocFromVolGeom(&gcam->image, MRI_UCHAR, 1, 1);
  GCAMinit(gcam, image_header, nullptr, &transform, 0);

  // transfer morph field
  double const *dptr = data.data();
  int framesize = gcam->depth * gcam->height * gcam->width;
  for (int zp = 0; zp < gcam->depth; zp++) {
    for (int yp = 0; yp < gcam->height; yp++) {
      for (int xp = 0; xp < gcam->width; xp++) {
        if (*dptr > 0) {
          // assume that any coordinate less than zero should be skipped
          gcam->nodes[xp][yp][zp].x = *(dptr);
          gcam->nodes[xp][yp][zp].y = *(dptr + framesize);
          gcam->nodes[xp][yp][zp].z = *(dptr + framesize * 2);
        }
        dptr++;
      }
    }
  }

  // write file
  GCAMwrite(gcam, filename.c_str());

  // free
  LTAfree(&lta);
  MRIfree(&image_header);
  GCAMfree(&gcam);
}


}  // end namespace morph
