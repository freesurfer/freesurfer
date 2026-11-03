#ifndef kvlLikelihoodImageFilterBase_hxx
#define kvlLikelihoodImageFilterBase_hxx


#include "kvlLikelihoodImageFilterBase.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"

namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
LikelihoodImageFilterBase< TInputImage >
::LikelihoodImageFilterBase()
{
}


//----------------------------------------------------------------------------
template< typename TInputImage >
LikelihoodImageFilterBase< TInputImage >
::~LikelihoodImageFilterBase()
{
}

} // end namespace kvl

#endif // KVLLIKELIHOODIMAGEFILTERBASE_HXX
