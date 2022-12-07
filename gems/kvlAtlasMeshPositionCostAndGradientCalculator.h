#ifndef __kvlAtlasMeshPositionCostAndGradientCalculator_h
#define __kvlAtlasMeshPositionCostAndGradientCalculator_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkAffineTransform.h"

#define KVL_ENABLE_TIME_PROBE 0

#if KVL_ENABLE_TIME_PROBE
  #include "itkTimeProbe.h"
#endif


namespace kvl
{


class AtlasMeshPositionCostAndGradientCalculator: public AtlasMeshRasterizor
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshPositionCostAndGradientCalculator  Self;
  typedef AtlasMeshRasterizor Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  typedef itk::Vector< ThreadAccumDataType, 3 > AtlasPositionGradientThreadAccumType;
  typedef itk::VectorContainer< AtlasMesh::PointIdentifier, AtlasPositionGradientThreadAccumType> AtlasPositionGradientThreadAccumContainerType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshPositionCostAndGradientCalculator, AtlasMeshRasterizor );

  /** Some typedefs */
  typedef itk::AffineTransform< double, 3 >  TransformType;

  /** */
  void SetMeshToImageTransform( const TransformType* meshToImageTransform );

  /** */
  void SetIgnoreDeformationPrior( bool ignoreDeformationPrior )
    { m_IgnoreDeformationPrior = ignoreDeformationPrior; }

  /** */
  void SetOnlyDeformationPrior( bool onlyDeformationPrior )
    { m_OnlyDeformationPrior = onlyDeformationPrior; }

  /** */
  const AtlasPositionGradientContainerType* GetPositionGradient() const
    {
    return m_PositionGradient;
    }

  /** */
  AtlasPositionGradientContainerType* GetPositionGradient()
    {
    return m_PositionGradient;
    }
  
  /** */
  double GetMinLogLikelihoodTimesPrior() const
    {
    return m_MinLogLikelihoodTimesPrior;
    }
    
  /** */
  double&  GetMinLogLikelihoodTimesPrior()
    {
    return m_MinLogLikelihoodTimesPrior;  
    }

  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  /**  Boundary conditions applied to gradient */
  enum BoundaryConditionType { NONE, SLIDING, AFFINE, TRANSLATION };
  void  SetBoundaryCondition( const BoundaryConditionType&  boundaryCondition )
    {
    m_BoundaryCondition = boundaryCondition;  
    }
    
  const BoundaryConditionType&  GetBoundaryCondition() const
    {
    return m_BoundaryCondition;  
    }  
  
  // for AtlasMeshToWishartGaussMixtureCostAndGradientCalculator
  virtual void SetDiffusionParameters( const int numberOfContrasts,
                                       const std::vector< double >&  wmmMixtureWeights,
                                       const std::vector< int >&  numberOfWishartsPerClass,
                                       const double& voxratio,
                                       const std::vector< double >& degreesOfFreedom,
                                       const std::vector< vnl_matrix< double > >& scaleMatrices )
  {
  }

  // for AtlasMeshToFrobeniusGaussMixtureCostAndGradientCalculator
  virtual void SetDiffusionParameters( const int numberOfContrasts,
                                       const std::vector< double >&  frobMixtureWeights,
                                       const std::vector< int >&  numberOfFrobeniusPerClass,
                                       const double& voxratio,
                                       const std::vector< double >& frobVariance,
                                       const std::vector< vnl_vector< double > >& frobMeans )
  {
  }

  // for AtlasMeshToDSWbetaGaussMixtureCostAndGradientCalculator
  virtual void SetDiffusionParameters( const int numberOfContrasts,
                                       const std::vector< double >&  DSWbetaMixtureWeights,
                                       const std::vector< int >&  numberOfDSWbetaePerClass,
                                       const double& voxratio,
                                       const std::vector< double >& DSWbetaAlpha,
                                       const std::vector< vnl_vector< double > >& DSWbetaMeans,
                                       const std::vector< double >& DSWbetaBeta,
                                       const std::vector< double >& DSWbetaConcentration,
                                       const std::vector< double >& logKummerSamples,
                                       const double& logKummerIncrement )
  {
  }

  
protected:
  AtlasMeshPositionCostAndGradientCalculator();
  virtual ~AtlasMeshPositionCostAndGradientCalculator();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );
  
  virtual void AddDataContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                                 const AtlasMesh::PointType& p1,
                                                 const AtlasMesh::PointType& p2,
                                                 const AtlasMesh::PointType& p3,
                                                 const AtlasAlphasType&  alphasInVertex0,
                                                 const AtlasAlphasType&  alphasInVertex1,
                                                 const AtlasAlphasType&  alphasInVertex2,
                                                 const AtlasAlphasType&  alphasInVertex3,
                                                 ThreadAccumDataType&  priorPlusDataCost,
                                                 AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                                 AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                                 AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                                 AtlasPositionGradientThreadAccumType&  gradientInVertex3 )
   {
   }  

  virtual bool AddPriorContributionOfTetrahedron( const AtlasMesh::PointType& p0,
                                                  const AtlasMesh::PointType& p1,
                                                  const AtlasMesh::PointType& p2,
                                                  const AtlasMesh::PointType& p3,
                                                  const ReferenceTetrahedronInfo& info,
                                                  ThreadAccumDataType&  priorPlusDataCost,
                                                  AtlasPositionGradientThreadAccumType&  gradientInVertex0,
                                                  AtlasPositionGradientThreadAccumType&  gradientInVertex1,
                                                  AtlasPositionGradientThreadAccumType&  gradientInVertex2,
                                                  AtlasPositionGradientThreadAccumType&  gradientInVertex3 );

  // Let's provide a "hook" for adding non-tetrahdron-based cost and gradient contributions
  virtual void PostProcessCostAndGradient( const AtlasMesh* mesh )
    {
    }

  //
  void ImposeSlidingBoundaryConditions( const AtlasMesh* mesh );
  
  //
  void ImposeAffineBoundaryConditions( const AtlasMesh* mesh );
 
  //
  void ImposeTranslationBoundaryConditions( const AtlasMesh* mesh );
  
  // 
  void ImposeBoundaryCondition( const AtlasMesh* mesh );

  //
  AtlasPositionGradientContainerType::Pointer  m_PositionGradient;
  bool  m_IgnoreDeformationPrior;
  bool  m_OnlyDeformationPrior;
  double  m_MinLogLikelihoodTimesPrior;
  bool  m_Abort;
  BoundaryConditionType  m_BoundaryCondition;

  //
  vnl_matrix< double >  m_AffineProjectionMatrix;
  
private:
  AtlasMeshPositionCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  typedef itk::Matrix< double >  SlidingBoundaryCorrectionMatrixType;
  SlidingBoundaryCorrectionMatrixType  m_SlidingBoundaryCorrectionMatrices[ 8 ]; 
  
  //
  std::vector< AtlasPositionGradientThreadAccumContainerType::Pointer >  m_ThreadSpecificPositionGradients;
  std::vector< ThreadAccumDataType >  m_ThreadSpecificMinLogLikelihoodTimesPriors;

#if KVL_ENABLE_TIME_PROBE  
  //
  std::vector< itk::TimeProbe >  m_ThreadSpecificDataTermRasterizationTimers;
  std::vector< itk::TimeProbe >  m_ThreadSpecificPriorTermRasterizationTimers;
  std::vector< itk::TimeProbe >  m_ThreadSpecificOtherRasterizationTimers;

#endif  
  
};


} // end namespace kvl

#endif
