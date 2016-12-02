#ifndef __kvlAtlasMeshDeformationLevenbergMarquardtOptimizer_h
#define __kvlAtlasMeshDeformationLevenbergMarquardtOptimizer_h

#include "kvlAtlasMeshDeformationLevenbergMarquardt.h"
#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{



/**
 *
 */
class AtlasMeshDeformationLevenbergMarquardtOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationLevenbergMarquardtOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationLevenbergMarquardtOptimizer, AtlasMeshDeformationOptimizer );

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

  //
  void  SetLambda( float lambda )
    { m_Lambda = lambda; }

  //
  float  GetLambda() const
    { return m_Lambda; }

 
  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    return m_LevenbergMarquardt->GetMinLogLikelihoodTimesPrior();
    }
  
  /** */
  bool Go();

  /** */
  float PerformOneSuccessfulStep();

protected:
  AtlasMeshDeformationLevenbergMarquardtOptimizer();
  virtual ~AtlasMeshDeformationLevenbergMarquardtOptimizer();
  
  void Initialize();


private:
  AtlasMeshDeformationLevenbergMarquardtOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  AtlasMeshDeformationLevenbergMarquardt::Pointer  m_LevenbergMarquardt;
  AtlasMeshDeformationLevenbergMarquardt::Pointer  m_TrialLevenbergMarquardt;

  float  m_Lambda;
  
  float  m_MaximalDeformationStopCriterion;
  
  
};


} // end namespace kvl

#endif

