

#include "surf_utils.h"

void
SurfacePointer::GetTypeFromName(const std::string& name)
{
  if ( name.find("white") != std::string::npos )
    this->type = white;
  else if (name.find("pial") != std::string::npos )
    this->type = pial;
  else
    this->type = other;

  std::cout << " surface name = " << name << " type = "
            << this->type << std::endl;
}

void
convert_surf_to_vox(MRI_SURFACE* mris,
                    MRI* vol)
{
  int const nvertices = mris->nvertices;

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (int vno=0; vno < nvertices; ++vno )
  {
    VERTEX const * v = &mris->vertices[vno];

    double vx, vy, vz;
    MRIvoxelToSurfaceRAS( vol,
                          v->x, v->y, v->z,
                          &vx, &vy, &vz);

    MRISsetXYZ(mris, vno, vx,vy,vz);
  }
}


void
convert_surf_to_vox(SurfaceVectorType& vmris,
                    MRI* mri)
{
  for ( SurfaceVectorType::iterator it = vmris.begin();
        it != vmris.end();
        ++it )
    convert_surf_to_vox( it->mris, mri );
}


void
convert_vox_to_surf(MRI_SURFACE* mris,
                    MRI* vol)
{
  int const nvertices = mris->nvertices;

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (int vno=0; vno < nvertices; ++vno )
  {
    VERTEX const * v = &mris->vertices[vno];

    double vx, vy, vz;
    MRIvoxelToSurfaceRAS( vol,
                          v->x, v->y, v->z,
                          &vx, &vy, &vz );

    MRISsetXYZ(mris, vno, vx,vy,vz);
  }
}

void
convert_vox_to_surf(SurfaceVectorType& vmris,
                    MRI* mri)
{
  for ( SurfaceVectorType::iterator it = vmris.begin();
        it != vmris.end();
        ++it )
    convert_vox_to_surf( it->mris, mri);
}

void
compute_bbox(const SurfaceVectorType& vmris,
             tDblCoords& cmin,
             tDblCoords& cmax)
{
  SurfaceVectorType::const_iterator cit = vmris.begin();
  compute_bbox( cit->mris, cmin, cmax, true);

  for ( ++cit; cit != vmris.end(); ++cit )
    compute_bbox( cit->mris, cmin, cmax, false);
}

void
compute_bbox(MRI_SURFACE* mris,
             tDblCoords& cmin,
             tDblCoords& cmax,
             bool initVars)
{
  tDblCoords cbuf;

  VERTEX* pvtx = &( mris->vertices[0] );
  unsigned int nvertices = (unsigned int)mris->nvertices;

  cbuf(0) = pvtx->x;
  cbuf(1) = pvtx->y;
  cbuf(2) = pvtx->z;

  if ( initVars )
  {
    cmin = cbuf;
    cmax = cbuf;
  }
  else
  {
    cmin = min(cmin, cbuf);
    cmax = max(cmax, cbuf);
  }
  ++pvtx;

  for ( unsigned int ui = 1;
        ui < nvertices;
        ++ui, ++pvtx )
  {
    cbuf(0) = pvtx->x;
    cbuf(1) = pvtx->y;
    cbuf(2) = pvtx->z;

    cmin = min(cmin, cbuf);
    cmax = max(cmax, cbuf);
  } // next ui, pvtx
}
