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

  m_LikelihoodFilter = LikelihoodFilterType::New();
  
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
    dynamic_cast<LikelihoodFilterType*>(m_LikelihoodFilter.GetPointer())
            ->SetParameters( means,
                             variances,
                             mixtureWeights,
                             numberOfGaussiansPerClass );
}

} // end namespace kvl
