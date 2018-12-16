
#include "surf_energy.h"


void
add_to_container(const MRI_SURFACE* mris_x,
                 const MRI_SURFACE* mris_fx,
                 PointsContainerType& container)
{
  const VERTEX* pvtx_x  = NULL;
  const VERTEX* pvtx_fx = NULL;
  Coords3d pt_x, pt_fx;

  pvtx_x  = &( mris_x->vertices[0] );
  pvtx_fx = &( mris_fx->vertices[0] );

  for ( unsigned int ui(0), nvertices(mris_x->nvertices);
        ui < nvertices; ++ui, ++pvtx_x, ++pvtx_fx )
  {
    pt_x(0) = pvtx_x->x;
    pt_x(1) = pvtx_x->y;
    pt_x(2) = pvtx_x->z;

    pt_fx(0) = pvtx_fx->x;
    pt_fx(1) = pvtx_fx->y;
    pt_fx(2) = pvtx_fx->z;

    container.push_back( std::make_pair(pt_x, pt_fx) );
  } // next ui, pvtx_x, pvtx_fx
}

float energy(const PointsContainerType& container,
             float* transform)
{
  Coords3d delta, img;

  float fenergy = 0.0f;

  for (PointsContainerType::const_iterator cit = container.begin();
       cit != container.end();
       ++cit )
  {
    for (int j=0; j<3; ++j)
      img(j) = transform[j]* cit->first(0) +
               transform[j+3] * cit->first(1) +
               transform[j+6] * cit->first(2) +
               transform[j+9];

    delta = img - cit->second;
    fenergy += delta.norm();
  } // next cit

  fenergy /= (double)container.size();

  return fenergy;
}

void
apply_lin_transform(MRI_SURFACE* mris, float* transform)
{
  std::cout << " linear transform = \n";
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<3; j++)
      std::cout << transform[i*3 + j] << "\t";
    std::cout << std::endl;
  }

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!

  for (int index = 0; index < mris->nvertices; index++)
  {
    VERTEX * const vertex = &mris->vertices[index];

    float  im_pt[3];

    for (int j=0; j<3; j++)
      im_pt[j] = transform[j] * vertex->x + transform[j+3]*vertex->y
                 + transform[j+6]* vertex->z + transform[j+9];

    MRISsetXYZ(mris, index,
      im_pt[0],
      im_pt[1],
      im_pt[2]);
  }
}
