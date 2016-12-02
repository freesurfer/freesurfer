#ifndef __kvlAtlasParameterEstimator_h
#define __kvlAtlasParameterEstimator_h

#include "itkImage.h"
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshDeformationOptimizer.h"
#include "kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"

namespace kvl
{

// Events generated
itkEventMacro( AlphasEstimationStartEvent, itk::UserEvent );
itkEventMacro( AlphasEstimationIterationEvent, itk::UserEvent );
itkEventMacro( AlphasEstimationEndEvent, itk::UserEvent );


//Eugenio
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


  void SetModeLM(){m_GD=false; m_DeformationOptimizer = AtlasMeshDeformationLevenbergMarquardtOptimizer::New();}
  void SetModeCJ(){m_GD=false; m_DeformationOptimizer = AtlasMeshDeformationConjugateGradientOptimizer::New();}
  void SetModeGD(){m_GD=true;}

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;
  
  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );
  
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
  void Estimate(bool verbose);

  void Estimate(){Estimate(false);}
  
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
    { return m_NumberOfClasses; }  
  

  // Set/Get mapping of collapsed labels.
  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
  { 
    int nc=0;
    bool done=false;
    for(int i=0; i<256; i++)
    {
      m_mapCompToComp[i] = mapCompToComp[i];
      if(mapCompToComp[i].size()==0 && done==false) { nc=i; done=true; }
    }
    m_NumberOfClasses = static_cast< unsigned int >( nc );
    // std::cout << "Calculated number of classes to be: " << m_NumberOfClasses << std::endl;
  }

  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }

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
  //unsigned int  GetPositionEstimationIterationNumber() const
  //  { return m_PositionEstimationIterationNumber; }
    
  //
  //unsigned int  GetPositionEstimationMaximumNumberOfIterations() const
  //  { return m_PositionEstimationMaximumNumberOfIterations; }
  
  // 
  AtlasPositionGradientContainerType::Pointer GetCurrentPositionGradient( unsigned int labelImageNumber ) const;
    
  //
  //void SetPositionGradientDescentStepSize( float stepSize ) 
  //  { m_PositionGradientDescentStepSize = stepSize; }
  
  //
  //float  GetPositionGradientDescentStepSize() const
  //  { return m_PositionGradientDescentStepSize; }
  
  //
  //unsigned int  GetPositionEstimationIterationEventResolution() const
  //  { return m_PositionEstimationIterationEventResolution; }
    
  //
  //void SetPositionEstimationIterationEventResolution( unsigned int  positionEstimationIterationEventResolution )
  //  { m_PositionEstimationIterationEventResolution = positionEstimationIterationEventResolution; }
    
  
  //
  void SetAlphaEstimationStopCriterion( float alphaEstimationStopCriterion )
    { m_AlphaEstimationStopCriterion = alphaEstimationStopCriterion; }
    
  float GetAlphaEstimationStopCriterion() const
    { return m_AlphaEstimationStopCriterion; }
  
  //
  //void SetPositionEstimationStopCriterion( float positionEstimationStopCriterion )
  //  { m_PositionEstimationStopCriterion = positionEstimationStopCriterion; }
    
  //float GetPositionEstimationStopCriterion() const
  //  { return m_PositionEstimationStopCriterion; }

  void SetDeformationOptimizer( AtlasMeshDeformationOptimizer* optimizer )
    { 
    m_DeformationOptimizer = optimizer;
    }
  
  const AtlasMeshDeformationOptimizer* GetDeformationOptimizer() const
    {
    return m_DeformationOptimizer;  
    }

  AtlasMeshDeformationOptimizer* GetDeformationOptimizer()
    {
    return m_DeformationOptimizer;  
    }

  //
  float  GetCurrentMinLogLikelihoodTimesPrior() const
    { return m_CurrentMinLogLikelihoodTimesPrior; }
    
  //
  float GetCurrentMinLogLikelihood() const;
    
  // 
  float GetMinLogLikelihood( const AtlasMeshCollection* meshCollection ) const;

  //  
  void SetAlphasSmoothingFactor( float alphasSmoothingFactor )
    { m_AlphasSmoothingFactor = alphasSmoothingFactor; }

  float  GetAlphasSmoothingFactor() const
    { return m_AlphasSmoothingFactor; }

  //
  void SetStopCriterion( float stopCriterion )
    { m_StopCriterion = stopCriterion; }
    
  float GetStopCriterion() const
    { return m_StopCriterion; }

  //
  void  SetUseGaussians( bool  useGaussians )
    { m_UseGaussians = useGaussians; }

  bool  GetUseGaussians() const
    { return m_UseGaussians; }

  void  SetIgnoreLastLabelImage( bool ignoreLastLabelImage )
    { m_IgnoreLastLabelImage = ignoreLastLabelImage; }

  /** */
  void SetNumberOfThreads( int numberOfThreads )
    { m_NumberOfThreads = numberOfThreads; }
    
  /** */
  int  GetNumberOfThreads() const
    { return m_NumberOfThreads; }

  /** */
  void SetGD( bool gd )
    { m_GD = gd; }
    
  /** */
  bool  GetGD() const
    { return m_GD; }




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
  float EstimatePositions();

  // 
  float EstimatePosition( unsigned int  labelImageNumber );
  
  // 
  AtlasPositionGradientContainerType::Pointer  
     CalculateCurrentPositionGradient( unsigned int labelImageNumber, float& minLogLikelihoodTimesPrior ) const;
  
     
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

  unsigned int  m_NumberOfClasses;
  std::vector<unsigned char > m_mapCompToComp[256];
  unsigned int  m_IterationNumber;
  unsigned int  m_MaximumNumberOfIterations;
  unsigned int  m_LabelImageNumber;
  unsigned int  m_NumberOfLabelImages;
  unsigned int  m_AlphasEstimationIterationNumber;
  unsigned int  m_AlphasEstimationMaximumNumberOfIterations;
  unsigned int  m_PositionEstimationIterationNumber;
  unsigned int  m_PositionEstimationMaximumNumberOfIterations;
  float  m_PositionGradientDescentStepSize;
  unsigned int  m_PositionEstimationIterationEventResolution;
  float  m_CurrentMinLogLikelihoodTimesPrior;
  float  m_AlphaEstimationStopCriterion;
  float  m_PositionEstimationStopCriterion;
  float  m_AlphasSmoothingFactor;
  float  m_StopCriterion;

  bool  m_GD;
  bool  m_UseGaussians;
  bool  m_IgnoreLastLabelImage;
  
  AtlasMeshDeformationOptimizer::Pointer  m_DeformationOptimizer;
  int  m_NumberOfThreads;
  
};



} // end namespace kvl


#endif
