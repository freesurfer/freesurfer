
#ifndef _h_misc_utils_h_
#define _h_misc_utils_h_

// STL includes
#include <vector>

#include "coords.h"
#include "fem_3d.h"



#include "mri.h"
#include "mrisurf.h"


enum SurfaceType
{
  white,
  pial,
  other
};


struct SurfacePointer
{
  MRI_SURFACE* mris;
  SurfaceType  type;
  void GetTypeFromName(const std::string& name);
};

typedef std::vector<SurfacePointer> SurfaceVectorType;

/*

Convert surface coords from RAS to IJK

*/

void convert_surf_to_vox(MRI_SURFACE* mris, MRI* mri);
void convert_surf_to_vox(SurfaceVectorType& vmris, MRI* mri);

void convert_vox_to_surf(MRI_SURFACE* mris, MRI* mri);
void convert_vox_to_surf(SurfaceVectorType& vmris, MRI* mri);

/*

Bounding Box

*/

void compute_bbox(const SurfaceVectorType& vmris,
                  tDblCoords& cmin,
                  tDblCoords& cmax);

void compute_bbox(MRI_SURFACE* mris,
                  tDblCoords& cmin,
                  tDblCoords& cmax,
                  bool initVars=true); // if true, will reset the variables

class SurfaceVertexIterator
{
public:
  SurfaceVertexIterator(MRIS* psurf = NULL)
      : m_psurf(psurf)
  {}
  void SetSurface(MRIS* psurf)
  {
    m_psurf = psurf;
  }

  template<class F> void Execute(F& f)
  {
    if (!m_psurf) throw std::logic_error("SurfaceVertexIterator - no surface");

    VERTEX* pvtx = &(m_psurf->vertices[0]);
    unsigned int ui=0, nvertices= m_psurf->nvertices;

    for (; ui < nvertices; ++ui, ++pvtx)
      f(pvtx);
  }

private:
  MRIS* m_psurf;
};

#endif
