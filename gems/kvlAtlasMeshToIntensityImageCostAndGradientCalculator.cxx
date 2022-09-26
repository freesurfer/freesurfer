#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculator.h"

#include <itkMath.h>
#include "vnl/vnl_matrix_fixed.h"
#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::AtlasMeshToIntensityImageCostAndGradientCalculator()
{
  /* 
   * 1. m_LikelihoodFilter is declared as LikelihoodImageFilterBase
   * 2. m_LikelihoodFilter is holding a GMMLikelihoodImageFilter object
   */ 

  m_LikelihoodFilter = gmmLikelihoodFilterType::New();

}


//
//
//
AtlasMeshToIntensityImageCostAndGradientCalculator
::~AtlasMeshToIntensityImageCostAndGradientCalculator()
{
}


//
//
//
void 
AtlasMeshToIntensityImageCostAndGradientCalculator
::SetParameters( const std::vector< vnl_vector< double > >& means,
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass )
{
   /* 
    * 1. m_LikelihoodFilter is declared as LikelihoodImageFilterBase
    * 2. m_LikelihoodFilter is holding a GMMLikelihoodImageFilter object
    * 3. dynamic_case m_LikelihoodFilter to GMMLikelihoodImageFilter, so it can call derived class method SetParameters()
    */ 
    dynamic_cast<gmmLikelihoodFilterType*>(m_LikelihoodFilter.GetPointer())
            ->SetParameters( means,
                             variances,
                             mixtureWeights,
                             numberOfGaussiansPerClass );
}

} // end namespace kvl
