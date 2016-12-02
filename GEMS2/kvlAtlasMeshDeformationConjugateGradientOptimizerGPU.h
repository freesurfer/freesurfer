#ifndef __kvlAtlasMeshDeformationConjugateGradientOptimizerGPU_h
#define __kvlAtlasMeshDeformationConjugateGradientOptimizerGPU_h

#include "kvlAtlasMeshToIntensityImageGradientCalculatorGPU.h"
#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{



/**
 *
 */
class AtlasMeshDeformationConjugateGradientOptimizerGPU: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationConjugateGradientOptimizerGPU  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationConjugateGradientOptimizerGPU, AtlasMeshDeformationOptimizer );

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
  AtlasMeshDeformationConjugateGradientOptimizerGPU();
  virtual ~AtlasMeshDeformationConjugateGradientOptimizerGPU();
  
  void Initialize();

  //
  void  GetCostAndGradient( const AtlasMesh::PointsContainer* position, double& cost, AtlasPositionGradientContainerType::Pointer& gradient );

  //
  double  BackTrack( AtlasMesh::PointsContainer::Pointer&  position, double& stepSize, double&  cost, AtlasPositionGradientContainerType::Pointer& gradient, 
                     const AtlasPositionGradientContainerType* direction, double directionalDerivative, double tolX, double referenceCost );
 
  
private:
  AtlasMeshDeformationConjugateGradientOptimizerGPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  float  m_MaximalDeformationStopCriterion;
  
  double  m_Cost;
  double  m_OldCost;
  AtlasMesh::PointsContainer::Pointer  m_Position;
  AtlasPositionGradientContainerType::Pointer  m_Gradient;
  AtlasPositionGradientContainerType::Pointer  m_OldGradient;
  AtlasPositionGradientContainerType::Pointer  m_Direction;
  //double  m_StepSize;
  
  // TODO: in the future, the pointer to the calculator should point to the 
  // base class. Based on what the user wants, the correct derived class is
  // then set up
  AtlasMeshToIntensityImageGradientCalculatorGPU::Pointer  m_Calculator;

  
};


} // end namespace kvl

#endif

