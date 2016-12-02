#ifndef __kvlAtlasMeshDeformationConjugateGradientOptimizer_h
#define __kvlAtlasMeshDeformationConjugateGradientOptimizer_h

#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"
#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{



/**
 *
 */
class AtlasMeshDeformationConjugateGradientOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationConjugateGradientOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationConjugateGradientOptimizer, AtlasMeshDeformationOptimizer );

  /** Some typedefs */
  typedef Superclass::ProbabilityImageType  ProbabilityImageType;
  typedef Superclass::SegmentedImageType  SegmentedImageType;
  typedef Superclass::ImageType  ImageType;
  typedef Superclass::TransformType  TransformType;

  //
  void  SetMaximalDeformationStopCriterion( float maximalDeformationStopCriterion )
    { m_MaximalDeformationStopCriterion = maximalDeformationStopCriterion; }

  //
  float  GetMaximalDeformationStopCriterion() const
    { return m_MaximalDeformationStopCriterion; }

  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    return m_Cost;
    }
  
  /** */
  bool Go();

  /** */
  double PerformOneIteration();
  

protected:
  AtlasMeshDeformationConjugateGradientOptimizer();
  virtual ~AtlasMeshDeformationConjugateGradientOptimizer();
  
  void Initialize();

  //
  void  GetCostAndGradient( const AtlasMesh::PointsContainer* position, double& cost, AtlasPositionGradientContainerType::Pointer& gradient );

  //
  double  BackTrack( AtlasMesh::PointsContainer::Pointer&  position, double& stepSize, double&  cost, AtlasPositionGradientContainerType::Pointer& gradient, 
                     const AtlasPositionGradientContainerType* direction, double directionalDerivative, double tolX, double referenceCost );
 
  
private:
  AtlasMeshDeformationConjugateGradientOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  float  m_MaximalDeformationStopCriterion;
  
  double  m_Cost;
  double  m_OldCost;
  AtlasMesh::PointsContainer::Pointer  m_Position;
  AtlasPositionGradientContainerType::Pointer  m_Gradient;
  AtlasPositionGradientContainerType::Pointer  m_OldGradient;
  AtlasPositionGradientContainerType::Pointer  m_Direction;
  //double  m_StepSize;
  
};


} // end namespace kvl

#endif

