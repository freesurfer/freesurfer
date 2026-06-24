#include "kvlAtlasMeshToPointSetCostAndGradientCalculator.h"

#include <itkMath.h>



namespace kvl
{

//
//
//
AtlasMeshToPointSetCostAndGradientCalculator
::AtlasMeshToPointSetCostAndGradientCalculator()
{

  m_TargetPoints = 0;

}


//
//
//
AtlasMeshToPointSetCostAndGradientCalculator
::~AtlasMeshToPointSetCostAndGradientCalculator()
{
}




//
//
//
void 
AtlasMeshToPointSetCostAndGradientCalculator
::PostProcessCostAndGradient( const AtlasMesh* mesh )
{

  // Sanity check
  if ( !m_TargetPoints )
    {
    itkExceptionMacro( << "No target point set!" );  
    }  
      
  // Get access to cost and gradient
  double&  minLogLikelihoodTimesPrior = this->GetMinLogLikelihoodTimesPrior();
  AtlasPositionGradientContainerType::Pointer  positionGradient = this->GetPositionGradient();
  
  // Loop over all points, and add contribution to cost and gradient
  AtlasMesh::PointsContainer::ConstIterator  targetPointIt = m_TargetPoints->Begin();
  AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
  AtlasPositionGradientContainerType::Iterator  gradientIt = positionGradient->Begin();
  for ( ; targetPointIt != m_TargetPoints->End(); ++ targetPointIt, ++pointIt, ++gradientIt )
    {
    for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )
      {
      // Skip sentinel values (considered as missing a target value)
      if ( targetPointIt.Value()[ dimensionNumber ] == 0 )
        {
        continue;  
        }
        
      // Add contribution to cost and gradient
      const double  distance = pointIt.Value()[ dimensionNumber ] - targetPointIt.Value()[ dimensionNumber ];  
      minLogLikelihoodTimesPrior += 0.5 * distance * distance;
      gradientIt.Value()[ dimensionNumber ] += distance;    
    
      } // Loop over all three dimensions
      
      
    } // Loop over all points
  
  
  
}



} // end namespace kvl
