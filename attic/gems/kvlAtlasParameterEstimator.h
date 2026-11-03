#ifndef __kvlAtlasParameterEstimator_h
#define __kvlAtlasParameterEstimator_h

#include "itkImage.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlCompressionLookupTable.h"


namespace kvl
{

// Events generated
itkEventMacro( AlphasEstimationStartEvent, itk::UserEvent );
itkEventMacro( AlphasEstimationIterationEvent, itk::UserEvent );
itkEventMacro( AlphasEstimationEndEvent, itk::UserEvent );
itkEventMacro( PositionEstimationStartEvent, itk::UserEvent );
itkEventMacro( PositionEstimationIterationEvent, itk::UserEvent );
itkEventMacro( PositionEstimationEndEvent, itk::UserEvent );


class AtlasParameterEstimator: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasParameterEstimator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasParameterEstimator, itk::Object );

  // Some typedefs
  typedef CompressionLookupTable::ImageType  LabelImageType;
  
  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages, 
                       const CompressionLookupTable*  lookupTable );

  //
  const CompressionLookupTable*  GetCompressionLookupTable() const
    { return m_CompressionLookupTable; }
  
  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
    { return m_LabelImages; }
    
  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;
      
  // Initialize
  void SetInitialMeshCollection( AtlasMeshCollection* initialMeshCollection )
    { m_MeshCollection = initialMeshCollection; }
  
  //
  const AtlasMeshCollection*  GetCurrentMeshCollection() const
    { return m_MeshCollection; }
    
  //   
  void  SetMaximumNumberOfIterations( unsigned int  maximumNumberOfIterations )
    { m_MaximumNumberOfIterations = maximumNumberOfIterations; }
  
  //
  void Estimate( bool verbose=false );
  
  // Get information about internal estimation state
  unsigned int  GetIterationNumber() const
    { return m_IterationNumber; }
    
  unsigned int  GetMaximumNumberOfIterations() const
    { return m_MaximumNumberOfIterations; }
       
  unsigned int  GetLabelImageNumber() const
    { return m_LabelImageNumber; }
  unsigned int  GetNumberOfLabelImages() const
    { return m_NumberOfLabelImages; }
  unsigned int  GetNumberOfClasses() const
    { return m_CompressionLookupTable->GetNumberOfClasses(); }  
  
  // 
  unsigned int  GetAlphasEstimationIterationNumber() const
    { return m_AlphasEstimationIterationNumber; }
    
  //
  unsigned int  GetAlphasEstimationMaximumNumberOfIterations() const
    { return m_AlphasEstimationMaximumNumberOfIterations; }

  //
  void  SetAlphasEstimationMaximumNumberOfIterations( unsigned int alphasEstimationMaximumNumberOfIterations )
    { m_AlphasEstimationMaximumNumberOfIterations = alphasEstimationMaximumNumberOfIterations; }

  // 
  unsigned int  GetPositionEstimationIterationNumber() const
    { return m_PositionEstimationIterationNumber; }
    
  //
  unsigned int  GetPositionEstimationMaximumNumberOfIterations() const
    { return m_PositionEstimationMaximumNumberOfIterations; }

  //
  void  SetPositionEstimationIterationEventResolution( unsigned int  positionEstimationIterationEventResolution )
    { m_PositionEstimationIterationEventResolution = positionEstimationIterationEventResolution; }
    
  //
  unsigned int  GetPositionEstimationIterationEventResolution() const
    { return m_PositionEstimationIterationEventResolution; }
    
  // 
  AtlasPositionGradientContainerType::Pointer GetCurrentPositionGradient( unsigned int labelImageNumber ) const;
    
  //
  void SetAlphaEstimationStopCriterion( double alphaEstimationStopCriterion )
    { m_AlphaEstimationStopCriterion = alphaEstimationStopCriterion; }
    
  double GetAlphaEstimationStopCriterion() const
    { return m_AlphaEstimationStopCriterion; }
  
  //
  double  GetCurrentMinLogLikelihoodTimesPrior() const
    { return m_CurrentMinLogLikelihoodTimesPrior; }
    
  //  
  void SetAlphasSmoothingFactor( double alphasSmoothingFactor )
    { m_AlphasSmoothingFactor = alphasSmoothingFactor; }

  double  GetAlphasSmoothingFactor() const
    { return m_AlphasSmoothingFactor; }

  //
  void SetStopCriterion( double stopCriterion )
    { m_StopCriterion = stopCriterion; }
    
  double GetStopCriterion() const
    { return m_StopCriterion; }

  /** */
  void SetNumberOfThreads( int numberOfThreads )
    { m_NumberOfThreads = numberOfThreads; }
    
  /** */
  int  GetNumberOfThreads() const
    { return m_NumberOfThreads; }

  /**  Position optimizer type */
  enum PositionOptimizerType { FIXED_STEP_GRADIENT_DESCENT, GRADIENT_DESCENT, CONJUGATE_GRADIENT, LBFGS };
  void  SetPositionOptimizer( const PositionOptimizerType&  positionOptimizer )
    {
    m_PositionOptimizer = positionOptimizer;  
    }
    
  const PositionOptimizerType&  GetPositionOptimizer() const
    {
    return m_PositionOptimizer;  
    }  


protected :
  // Constructor
  AtlasParameterEstimator();
  
  // Destructor
  virtual ~AtlasParameterEstimator();
  
  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;  
  
  //
  virtual void EstimateAlphas();
    
  //
  virtual void SmoothAlphas();
    
  // Estimate positions  
  double EstimatePositions();

  // 
  double EstimatePosition( unsigned int  labelImageNumber );
  
  // 
  AtlasPositionGradientContainerType::Pointer  
     CalculateCurrentPositionCostAndGradient( unsigned int labelImageNumber, double& minLogLikelihoodTimesPrior ) const;
  
  //
  void HandleOptimizerEvent( itk::Object* object, const itk::EventObject & event );
     
private :
  AtlasParameterEstimator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Some typedefs
  typedef AtlasMesh::PointsContainer  PointsContainerType;
  typedef AtlasMesh::PointIdentifier  PointIdentifierType;
  typedef AtlasMesh::PointType  PointType;
  typedef AtlasMesh::PointDataContainer  PointDataContainerType;
  typedef AtlasMesh::PixelType  AlphasType;
  
  typedef AtlasMesh::CellsContainer  CellsContainerType;
  typedef AtlasMesh::CellIdentifier  CellIdentifierType;
  typedef AtlasMesh::CellType  CellType;
  typedef AtlasMesh::CellDataContainer  CellDataContainerType;
  

    
  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;
  
  unsigned int  m_IterationNumber;
  unsigned int  m_MaximumNumberOfIterations;
  unsigned int  m_LabelImageNumber;
  unsigned int  m_NumberOfLabelImages;
  unsigned int  m_AlphasEstimationIterationNumber;
  unsigned int  m_AlphasEstimationMaximumNumberOfIterations;
  unsigned int  m_PositionEstimationIterationNumber;
  unsigned int  m_PositionEstimationMaximumNumberOfIterations;
  unsigned int  m_PositionEstimationIterationEventResolution; 
  double  m_CurrentMinLogLikelihoodTimesPrior;
  double  m_AlphaEstimationStopCriterion;
  double  m_AlphasSmoothingFactor;
  double  m_StopCriterion;
  
  PositionOptimizerType  m_PositionOptimizer;
  int  m_NumberOfThreads;
  
};



} // end namespace kvl


#endif
