#ifndef __kvlAtlasMeshToPointSetCostAndGradientCalculator_h
#define __kvlAtlasMeshToPointSetCostAndGradientCalculator_h

#include "kvlAtlasMeshPositionCostAndGradientCalculator.h"


namespace kvl
{


class AtlasMeshToPointSetCostAndGradientCalculator: public AtlasMeshPositionCostAndGradientCalculator
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToPointSetCostAndGradientCalculator  Self;
  typedef AtlasMeshPositionCostAndGradientCalculator Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToPointSetCostAndGradientCalculator, AtlasMeshPositionCostAndGradientCalculator );

  // Set target points
  void SetTargetPoints( const AtlasMesh::PointsContainer* targetPoints )
    {
    m_TargetPoints = targetPoints;
    }
    
  
protected:
  AtlasMeshToPointSetCostAndGradientCalculator();
  virtual ~AtlasMeshToPointSetCostAndGradientCalculator();
  
  //
  void PostProcessCostAndGradient( const AtlasMesh* mesh );
  
private:
  AtlasMeshToPointSetCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  AtlasMesh::PointsContainer::ConstPointer  m_TargetPoints;
  
};


} // end namespace kvl

#endif
