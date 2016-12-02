#ifndef __kvlAtlasMeshAveragingConjugateGradientOptimizer_h
#define __kvlAtlasMeshAveragingConjugateGradientOptimizer_h

#include "kvlAtlasMeshToIntensityImageGradientCalculator.h"
#include "kvlAtlasMeshDeformationOptimizer.h"
#include "kvlAtlasMeshCollection.h"


namespace kvl
{



/**
 *
 */
class AtlasMeshAveragingConjugateGradientOptimizer: public AtlasMeshDeformationOptimizer
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshAveragingConjugateGradientOptimizer  Self;
  typedef AtlasMeshDeformationOptimizer  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshAveragingConjugateGradientOptimizer, AtlasMeshDeformationOptimizer );

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

  void  SetInputMeshCollection( kvl::AtlasMeshCollection* meshCollection )
    {
    m_InputMeshCollection = meshCollection;
    }

  const  kvl::AtlasMeshCollection*  GetInputMeshCollection() const
    {
    return m_InputMeshCollection;
    }

  void  SetInitialMeshCollection( kvl::AtlasMeshCollection* meshCollection )
    {
    m_InitialMeshCollection = meshCollection;
    }

  const  kvl::AtlasMeshCollection*  GetInitialMeshCollection() const
    {
    return m_InitialMeshCollection;
    }
  
  void SetKatlas (float Katlas) { m_Katlas = Katlas; }
  void SetKtimePoints (float KtimePoints) { m_KtimePoints = KtimePoints; }
  float GetKatlas () {return m_Katlas;}
  float GetKtimePoints () {return m_KtimePoints;}

  void Initialize();

protected:
  AtlasMeshAveragingConjugateGradientOptimizer();
  virtual ~AtlasMeshAveragingConjugateGradientOptimizer();
  
  // void Initialize();

  //
  void  GetCostAndGradient( const AtlasMesh::PointsContainer* position, double& cost, AtlasPositionGradientContainerType::Pointer& gradient );

  //
  double  BackTrack( AtlasMesh::PointsContainer::Pointer&  position, double& stepSize, double&  cost, AtlasPositionGradientContainerType::Pointer& gradient, 
                     const AtlasPositionGradientContainerType* direction, double directionalDerivative, double tolX, double referenceCost );
 
private:
  AtlasMeshAveragingConjugateGradientOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  float  m_MaximalDeformationStopCriterion;
  
  double  m_Cost;
  double  m_OldCost;
  AtlasMesh::PointsContainer::Pointer  m_Position;
  AtlasPositionGradientContainerType::Pointer  m_Gradient;
  AtlasPositionGradientContainerType::Pointer  m_OldGradient;
  AtlasPositionGradientContainerType::Pointer  m_Direction;
  //double  m_StepSize;
  kvl::AtlasMeshCollection::Pointer  m_InputMeshCollection;
  kvl::AtlasMeshCollection::Pointer  m_InitialMeshCollection;
  float m_Katlas;
  float m_KtimePoints;
};


} // end namespace kvl

#endif

