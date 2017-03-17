#ifndef __kvlAtlasMeshBuilder_h
#define __kvlAtlasMeshBuilder_h

#include "kvlMultiResolutionAtlasMesher.h"
#include "vnl/vnl_sample.h"
#include "itkTimeProbe.h"


namespace kvl
{


class AtlasMeshBuilderMutexLock: public itk::SimpleFastMutexLock
{
public:
  /** Standard class typedefs.  */
  typedef AtlasMeshBuilderMutexLock       Self;
  typedef itk::SimpleFastMutexLock        Superclass;

  /** Lock access. */
  void DescriptiveLock( const std::string& description )
  {
    Superclass::Lock();
    m_TimeProbe.Start();
    m_Description = description;
  }

  /** Unlock access. */
  void DescriptiveUnlock()
  {
    m_TimeProbe.Stop();
    std::cout << m_Description << ": unlocking mutex after " << m_TimeProbe.GetMean() << " seconds" << std::endl;
    Superclass::Unlock();
  }

protected:

  std::string  m_Description;
  itk::TimeProbe  m_TimeProbe;

};



// Events generated
itkEventMacro( StartEdgeAnalysisEvent, itk::UserEvent );
//itkEventMacro( EdgeAnalysisProgressEvent, itk::UserEvent );
itkEventMacro( EndEdgeAnalysisEvent, itk::UserEvent );



class AtlasMeshBuilder: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshBuilder  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshBuilder, itk::Object );

  // Some typedefs
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
  {
    return m_LabelImages;
  }

// Set/Get mapping of collapsed labels.
  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { for(int i=0; i<256; i++) m_mapCompToComp[i] = mapCompToComp[i]; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }

  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;

  // Set/Get initial size
  void  SetInitialSize( unsigned int* size );
  const unsigned int* GetInitialSize() const
  {
    return m_InitialSize;
  }

  //
  void  SetNumberOfUpsamplingSteps( unsigned int  numberOfUpsamplingSteps )
  {
    m_NumberOfUpsamplingSteps = numberOfUpsamplingSteps;
  }
  const unsigned int  GetNumberOfUpsamplingSteps() const
  {
    return m_NumberOfUpsamplingSteps;
  }

  // Set/Get initial stiffness
  void  SetInitialStiffnesses( const float* stiffnesses );
  float  GetInitialStiffness( unsigned int upsampleStepNumber ) const
  {
    return m_InitialStiffnesses[ upsampleStepNumber ];
  }

  //
  void  SetInitialMesh( unsigned int* size, const float* initialStiffnesses );

  //
  unsigned int GetNumberOfClasses() const
  {
    return m_NumberOfClasses;
  }

  //
  const AtlasMeshCollection*  GetCurrentMeshCollection() const
  {
    return m_Current;
  }

  //
  const AtlasParameterEstimator* GetEstimator() const
  {
    return m_Mesher->GetEstimator();
  }

  //
  void SetPositionEstimationIterationEventResolution( unsigned int resolution )
  {
    m_Mesher->SetPositionEstimationIterationEventResolution( resolution );
  }
  unsigned int GetPositionEstimationIterationEventResolution() const
  {
    return m_Mesher->GetPositionEstimationIterationEventResolution();
  }

  //
  //void SetPositionGradientDescentStepSize( float stepSize )
  //  { m_Mesher->SetPositionGradientDescentStepSize( stepSize ); }
  //float  GetPositionGradientDescentStepSize() const
  //  { return m_Mesher->GetPositionGradientDescentStepSize(); }

  //
  void SetAlphaEstimationStopCriterion( float stopCriterion )
  {
    m_Mesher->SetAlphaEstimationStopCriterion( stopCriterion );
  }
  float GetAlphaEstimationStopCriterion() const
  {
    return m_Mesher->GetAlphaEstimationStopCriterion();
  }

  //
  //void SetPositionEstimationStopCriterion( float stopCriterion )
  //  { m_Mesher->SetPositionEstimationStopCriterion( stopCriterion ); }
  //float GetPositionEstimationStopCriterion() const
  //  { return m_Mesher->GetPositionEstimationStopCriterion(); }

  //
  void SetPowellAbsolutePrecision( float powellAbsolutePrecision )
  {
    m_PowellAbsolutePrecision = powellAbsolutePrecision;
  }
  float GetPowellAbsolutePrecision() const
  {
    return m_PowellAbsolutePrecision;
  }


  //
  unsigned long AddEstimatorObserver( const itk::EventObject& event, itk::Command* command )
  {
    return m_Mesher->AddEstimatorObserver( event, command );
  }


  //
  float GetCurrentDataCost() const;

  //
  float GetCurrentAlphasCost() const;

  //
  void GetCurrentDataAndAlphasCost( float& currentDataCost, float& currentAlphasCost ) const;

  //
  float GetCurrentPositionCost() const;

  //
  float GetCurrentCost() const;


#if 0
  //
  void PrintPreviousGains( std::ostream& os ) const
  {
    m_Gains.Print( os );
  }

#endif


  //
  unsigned int GetIterationNumber() const
  {
    return m_IterationNumber;
  }

  //
  unsigned int GetMaximumNumberOfIterations() const
  {
    return m_MaximumNumberOfIterations;
  }

  //
  float GetProgress() const
  {
    return m_Progress;
  }

  //
  AtlasMesh::CellIdentifier  GetEdgeId() const
  {
    return m_EdgeId;
  }

  //
  const AtlasMeshCollection*  GetRetainedMiniCollection() const
  {
    return m_RetainedMiniCollection;
  }

  //
  const AtlasMeshCollection*  GetCollapsedMiniCollection() const
  {
    return m_CollapsedMiniCollection;
  }

  //
  float  GetRetainedCosts( float& retainedDataCost, float& retainedAlphasCost, float& retainedPositionCost ) const
  {
    retainedDataCost = m_RetainedDataCost;
    retainedAlphasCost = m_RetainedAlphasCost;
    retainedPositionCost = m_RetainedPositionCost;
    return m_RetainedCost;
  }

  //
  float  GetCollapsedCosts( float& collapsedDataCost, float& collapsedAlphasCost, float& collapsedPositionCost ) const
  {
    collapsedDataCost = m_CollapsedDataCost;
    collapsedAlphasCost = m_CollapsedAlphasCost;
    collapsedPositionCost = m_CollapsedPositionCost;
    return m_CollapsedCost;
  }

  //
  void  Build( AtlasMeshCollection* explicitStartCollection = 0 );


protected :
  // Constructor
  AtlasMeshBuilder();

  // Destructor
  virtual ~AtlasMeshBuilder();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

  //
  virtual void SetUp();

  //
  unsigned int CountNumberOfEdges( const AtlasMeshCollection* meshCollection ) const;

  //
  float  GetPositionCost( const AtlasMeshCollection* meshCollection ) const;

  //
  void  GetDataCostAndAlphasCost( const AtlasMeshCollection* meshCollection, float& dataCost, float& alphasCost ) const;

  //
  AtlasMeshCollection::Pointer
  TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                     AtlasMesh::CellIdentifier  edgeId,
                     float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const;


  //
  AtlasMeshCollection::Pointer
  TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                     AtlasMesh::CellIdentifier  edgeId,
                     float& miniDataCost, float& miniAlphasCost, float& miniPositionCost) const
  {
    std::set< AtlasMesh::CellIdentifier >  disappearingCellsDummy;
    return this->TryToCollapseFast( miniCollection,
                                    edgeId,
                                    miniDataCost, miniAlphasCost, miniPositionCost,
                                    disappearingCellsDummy );
  }

  //
  void ExpandCollapse( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
                       const AtlasMeshCollection*  optimizedOuterChild,
                       AtlasMeshCollection::Pointer& result,
                       std::set< AtlasMesh::CellIdentifier >*  disappearingCells,
                       AtlasMesh::CellIdentifier& unifiedVertexId ) const;


  //
  bool  LoadBalancedAnalyzeEdgeFast( std::set< AtlasMesh::CellIdentifier >&  edges,
                                     std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies,
                                     int threadId=0 );

  //
  void
  AnalyzeEdgeFast( const AtlasMeshCollection*  miniCollection,
                   AtlasMesh::CellIdentifier edgeId,
                   AtlasMesh::PointIdentifier  edgePoint0Id, AtlasMesh::PointIdentifier edgePoint1Id,
                   std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >&  disappearingPointsLookupTable,
                   AtlasMesh::PointsContainer*  newReferencePosition,
                   std::vector< AtlasMesh::PointsContainer::Pointer >& newPositions,
                   AtlasMesh::PointDataContainer*  newPointParameters,
                   std::set< AtlasMesh::CellIdentifier >&  disappearingCells );

  //
  std::vector< AtlasMesh::CellIdentifier >
  GetIndependentEdges( int maximumNumberOfIndependentEdges,
                       std::vector< AtlasMesh::CellIdentifier >&  edges ) const;

  //
  std::vector< AtlasMesh::CellIdentifier >
  GetIndependentEdges( int maximumNumberOfIndependentEdges =  itk::NumericTraits< int >::max() ) const
  {
    std::vector< AtlasMesh::CellIdentifier >  edgesDummy;

    return this->GetIndependentEdges( maximumNumberOfIndependentEdges, edgesDummy );
  }

  //
  std::vector< AtlasMesh::CellIdentifier >   GetRandomizedEdges() const;
  std::set< AtlasMesh::CellIdentifier >   GetRandomizedEdgesAsSet() const;



  //
  AtlasMeshCollection::Pointer
  TryToRetainFast( const AtlasMeshCollection* miniCollectionConst,
                   AtlasMesh::CellIdentifier  edgeId,
                   float& miniDataCost, float& miniAlphasCost, float& miniPositionCost );

  //
  void  ExpandRetain( AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
                      const AtlasMeshCollection*  outerMiniCollection,
                      AtlasMeshCollection::Pointer& result );


  //
  bool  OptimizeReferencePositionFast( AtlasMeshCollection* meshCollection,
                                       AtlasMesh::CellIdentifier  vertexId, bool optimize = true ) const;

  //
  std::vector< AtlasMesh::CellIdentifier >   Permutate(  const std::vector< AtlasMesh::CellIdentifier >&  edgeList ) const;

  //
  AtlasMeshCollection::Pointer   GetFakeCopy( const AtlasMeshCollection*  input ) const;

  //
  AtlasMeshCollection::Pointer
  ApplyParallellMeshOperations( const AtlasMeshCollection*  collection,
                                const std::set< AtlasMesh::CellIdentifier >&   disappearingCells,
                                const std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >&
                                disappearingPointsLookupTable,
                                const AtlasMesh::PointsContainer*  newReferencePosition,
                                const std::vector< AtlasMesh::PointsContainer::Pointer >&  newPositions,
                                const AtlasMesh::PointDataContainer*  newPointParameters ) const;


  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE LoadBalancedThreaderCallback( void *arg );

  /** Internal structure used for passing image data into the threading library */
  struct LoadBalancedThreadStruct
  {
    //Pointer  m_Builder;
    Self*  m_Builder;
    std::set< AtlasMesh::CellIdentifier >  m_Edges;
    std::map< AtlasMesh::PointIdentifier, int >  m_PointOccupancies;
  };



private :
  AtlasMeshBuilder(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int m_stuckCount;

  //
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  unsigned int  m_InitialSize[ 3 ];
  unsigned int  m_NumberOfUpsamplingSteps;
  float  m_InitialStiffnesses[ 5 ];
  MultiResolutionAtlasMesher::Pointer  m_Mesher;

  float  m_PowellAbsolutePrecision;

  std::vector<unsigned char > m_mapCompToComp[256];

  unsigned int  m_NumberOfClasses;

  unsigned int  m_EdgeCollapseNumber;
  unsigned int  m_MaximumNumberOfEdgeCollapses;

  AtlasMeshCollection::Pointer  m_Current;

  unsigned int  m_IterationNumber;
  unsigned int  m_MaximumNumberOfIterations;
  float  m_Progress;
  AtlasMesh::CellIdentifier  m_EdgeId;

  AtlasMeshCollection::Pointer  m_RetainedMiniCollection;
  AtlasMeshCollection::Pointer  m_CollapsedMiniCollection;

  float  m_RetainedCost;
  float  m_RetainedDataCost;
  float  m_RetainedAlphasCost;
  float  m_RetainedPositionCost;

  float  m_CollapsedCost;
  float  m_CollapsedDataCost;
  float  m_CollapsedAlphasCost;
  float  m_CollapsedPositionCost;

  // Useful class for permutating edge list
  class EdgeElement
  {
  public:
    //
    EdgeElement()
    {
      m_EdgeId = 0;
      m_EdgePriority = vnl_sample_uniform( 0,  1 );
    }

    //
    EdgeElement( AtlasMesh::CellIdentifier edgeId )
    {
      m_EdgeId = edgeId;
      m_EdgePriority = vnl_sample_uniform( 0,  1 );
    }

    //
    ~EdgeElement() {}

    // Needed for std::sort()
    bool operator<(const EdgeElement& it) const
    {
      return m_EdgePriority < it.m_EdgePriority;
    }

    //
    AtlasMesh::CellIdentifier  GetEdgeId() const
    {
      return m_EdgeId;
    }

  private:
    AtlasMesh::CellIdentifier   m_EdgeId;
    double  m_EdgePriority;
  };


};



} // end namespace kvl


#endif
