
#include "morph_utils.h"

struct FunctorInitOctree
{
  template<class T>
  T operator()(T p) const
  {
    p->performInit();
    return p;
  }
};

void
initOctree( gmp::VolumeMorph& morph)
{
  std::transform( morph.m_transforms.begin(),
                  morph.m_transforms.end(),
                  morph.m_transforms.begin(),
                  FunctorInitOctree()
                );
}
