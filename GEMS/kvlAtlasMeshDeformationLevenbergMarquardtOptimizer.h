#ifndef __kvlAtlasMeshDeformationLevenbergMarquardtOptimizer_h
#define __kvlAtlasMeshDeformationLevenbergMarquardtOptimizer_h

#include "kvlAtlasMeshDeformationLevenbergMarquardt2.h"


namespace kvl
{



/**
 *
 */
class AtlasMeshDeformationLevenbergMarquardtOptimizer: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshDeformationLevenbergMarquardtOptimizer  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshDeformationLevenbergMarquardtOptimizer, itk::Object );

  /** Some typedefs */
  typedef AtlasMeshDeformationLevenbergMarquardt::LabelImageType  LabelImageType;
  typedef AtlasMeshDeformationLevenbergMarquardt::ProbabilityImageType  ProbabilityImageType;
  typedef AtlasMeshDeformationLevenbergMarquardt::ImageType  ImageType;
  typedef AtlasMeshDeformationLevenbergMarquardt::TransformType  TransformType;

  /** */
  void  SetProbabilityImage( ProbabilityImageType* probabilityImage )
    {
    m_ProbabilityImage = probabilityImage;
    m_Initialized = false;
    }

  /** */
  const ProbabilityImageType*  GetProbabilityImage() const
    {
    return m_ProbabilityImage;
    }

  /** */
  void  SetImage( ImageType* image )
    {
    m_Image = image;
    m_Initialized = false;
    }

  /** */
  const ImageType*  GetImage() const
    {
    return m_Image;
    }

  /** */
  void  SetMesh( AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    m_Initialized = false;
    }

  /** */
  const AtlasMesh*  GetMesh() const
    {
    return m_Mesh;
    }


  /** */
  void SetMeans( const itk::Array< float >& means )
    { 
    m_Means = means;
    m_Initialized = false;
    }

  /** */ 
  void SetVariances( const itk::Array< float >& variances )
    {
    m_Variances = variances;
    m_Initialized = false;
    }

  /** */
  void  SetMeshToImageTransform( TransformType* meshToImageTransform )
    {
    m_MeshToImageTransform = meshToImageTransform;
    m_Initialized = false;  
    }

  const TransformType* GetMeshToImageTransform() const
    {
    return m_MeshToImageTransform;
    }

  /** */
  void  SetUseProbabilityImage( bool  useProbabilityImage )
    {
    m_UseProbabilityImage = useProbabilityImage;
    m_Initialized = false;
    }

  /** */
  bool  GetUseProbabilityImage() const
    {
    return m_UseProbabilityImage;
    }


  //
  void  SetMaximalDeformationStopCriterion( float maximalDeformationStopCriterion )
    { m_MaximalDeformationStopCriterion = maximalDeformationStopCriterion; }

  //
  float  GetMaximalDeformationStopCriterion() const
    { return m_MaximalDeformationStopCriterion; }

  // 
  unsigned int  GetPositionUpdatingIterationNumber() const
    { return m_PositionUpdatingIterationNumber; }
    
  //
  unsigned int  GetPositionUpdatingMaximumNumberOfIterations() const
    { return m_PositionUpdatingMaximumNumberOfIterations; }

  void  SetPositionUpdatingMaximumNumberOfIterations( unsigned int positionUpdatingMaximumNumberOfIterations )
    { m_PositionUpdatingMaximumNumberOfIterations = positionUpdatingMaximumNumberOfIterations; }

  //
  unsigned int  GetPositionUpdatingIterationEventResolution() const
    { return m_PositionUpdatingIterationEventResolution; }
    
  //
  void SetPositionUpdatingIterationEventResolution( unsigned int  positionUpdatingIterationEventResolution )
    { m_PositionUpdatingIterationEventResolution = positionUpdatingIterationEventResolution; }

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
  
  ImageType::Pointer  m_Image;
  ProbabilityImageType::Pointer  m_ProbabilityImage;
  AtlasMesh::Pointer  m_Mesh;
  TransformType::Pointer  m_MeshToImageTransform;

  itk::Array< float >  m_Means;
  itk::Array< float >  m_Variances;
  
  bool  m_Initialized;
  bool  m_UseProbabilityImage;

  AtlasMeshDeformationLevenbergMarquardt::Pointer  m_LevenbergMarquardt;
  AtlasMeshDeformationLevenbergMarquardt::Pointer  m_TrialLevenbergMarquardt;

  float  m_Lambda;
  
  
  int  m_PositionUpdatingIterationNumber;
  int  m_PositionUpdatingMaximumNumberOfIterations;
  int  m_PositionUpdatingIterationEventResolution;
  float  m_MaximalDeformationStopCriterion;

  
};


} // end namespace kvl

#endif

