#include "itkSingleValuedCostFunction.h"
#include "kvlAtlasParameterEstimator.h"
#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"


namespace kvl
{



class AtlasMeshCollectionFastReferencePositionCost : public itk::SingleValuedCostFunction
{
  
public:
  /** Standard class typedefs. */
  typedef AtlasMeshCollectionFastReferencePositionCost     Self;
  typedef itk::SingleValuedCostFunction                 Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionFastReferencePositionCost, itk::SingleValuedCostFunction );
  
  /** */
  itkNewMacro(Self);

  // Some typedefs
  typedef Superclass::ParametersType  ParametersType;
  typedef Superclass::MeasureType  MeasureType;
  typedef Superclass::DerivativeType  DerivativeType;

  //
  virtual MeasureType GetValue( const ParametersType & parameters ) const;

  //
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const
    {
    }

  virtual unsigned int GetNumberOfParameters() const
    { return 3; }


  // Some typedefs
  typedef CompressionLookupTable::ImageType  LabelImageType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages,
                       const CompressionLookupTable*  compressionLookupTable );

  // Set initial collection
  void SetInitialMeshCollection( AtlasMeshCollection* meshCollection );

  // 
  void SetPointId( AtlasMesh::PointIdentifier pointId )
    { m_PointId = pointId; m_VertexNeighborhood.clear(); }

  // Get the atlas mesh collection that was used to calculate the cost the last time GetValue()
  // was called
  const AtlasMeshCollection*  GetCostCalculationMeshCollection() const
    { return m_CostCalculationMeshCollection; }

  //
  bool GetValue( const ParametersType & parameters, 
                 double& dataCost, double& alphasCost, double& positionCost ) const;

protected:
  AtlasMeshCollectionFastReferencePositionCost();
  virtual ~AtlasMeshCollectionFastReferencePositionCost();

  void CalculateVertexNeighborhood() const;

private:
  AtlasMeshCollectionFastReferencePositionCost(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasParameterEstimator::Pointer  m_Estimator;
  AtlasMeshCollectionModelLikelihoodCalculator::Pointer  m_DataAndAlphasCostCalculator;

  AtlasMesh::PointIdentifier  m_PointId;
  std::vector< AtlasMesh::PointsContainer::Pointer >   m_InitialPositions;
  AtlasMesh::PointDataContainer::Pointer  m_InitialPointParameters;
  AtlasMesh::CellsContainer::Pointer  m_Cells;
  AtlasMesh::PointsContainer::Pointer  m_ReferencePosition;  // This is actually what we're gonna change
  double  m_K;

  // 
  struct VertexNeighboringTetrahedronInfo
    {
    AtlasMesh::CellIdentifier  m_TetrahedronId;

    AtlasMesh::PointType::ValueType  m_X1;
    AtlasMesh::PointType::ValueType  m_Y1;
    AtlasMesh::PointType::ValueType  m_Z1;
    
    AtlasMesh::PointType::ValueType  m_X2;
    AtlasMesh::PointType::ValueType  m_Y2;
    AtlasMesh::PointType::ValueType  m_Z2;

    AtlasMesh::PointType::ValueType  m_X3;
    AtlasMesh::PointType::ValueType  m_Y3;
    AtlasMesh::PointType::ValueType  m_Z3;

    // Constructor
    VertexNeighboringTetrahedronInfo() : m_TetrahedronId( 0 ),
                                         m_X1( 0 ), m_Y1( 0 ), m_Z1( 0 ),
                                         m_X2( 0 ), m_Y2( 0 ), m_Z2( 0 ),
                                         m_X3( 0 ), m_Y3( 0 ), m_Z3( 0 ) {}
    };

  mutable std::vector< VertexNeighboringTetrahedronInfo >  m_VertexNeighborhood; // Cache
  mutable AtlasMeshCollection::Pointer  m_CostCalculationMeshCollection;
};



} // End namespace kvl
