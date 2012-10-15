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
    std::cout << m_Description << ": unlocking mutex after " << m_TimeProbe.GetMeanTime() << " seconds" << std::endl;
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

#if 0
  //
  const AtlasMeshCollection*  GetSplittedMiniCollection() const
  {
    return m_SplittedMiniCollection;
  }
  //
  const AtlasMeshCollection*  GetSwappedMiniCollection() const
  {
    return m_SwappedMiniCollection;
  }
#endif

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

#if 0
  //
  float  GetSplittedCosts( float& splittedDataCost, float& splittedAlphasCost, float& splittedPositionCost ) const
  {
    splittedDataCost = m_SplittedDataCost;
    splittedAlphasCost = m_SplittedAlphasCost;
    splittedPositionCost = m_SplittedPositionCost;
    return m_SplittedCost;
  }

  //
  float  GetSwappedCosts( float& swappedDataCost, float& swappedAlphasCost, float& swappedPositionCost ) const
  {
    swappedDataCost = m_SwappedDataCost;
    swappedAlphasCost = m_SwappedAlphasCost;
    swappedPositionCost = m_SwappedPositionCost;
    return m_SwappedCost;
  }

  //
  void  BuildWithPriorityQueue( AtlasMeshCollection* explicitStartCollection = 0 );
#endif

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
  TryToCollapse( const AtlasMeshCollection* innerMiniCollection,
                 const AtlasMeshCollection* outerMiniCollection,
                 AtlasMesh::CellIdentifier  edgeId,
                 float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                 std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const;

  //
  AtlasMeshCollection::Pointer
  TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                     AtlasMesh::CellIdentifier  edgeId,
                     float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const;


  //
  AtlasMeshCollection::Pointer
  TryToCollapse( const AtlasMeshCollection* innerMiniCollection,
                 const AtlasMeshCollection* outerMiniCollection,
                 AtlasMesh::CellIdentifier  edgeId,
                 float& miniDataCost, float& miniAlphasCost, float& miniPositionCost) const
  {
    std::set< AtlasMesh::CellIdentifier >  disappearingCellsDummy;
    return this->TryToCollapse( innerMiniCollection,
                                outerMiniCollection,
                                edgeId,
                                miniDataCost, miniAlphasCost, miniPositionCost,
                                disappearingCellsDummy );
  }

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
  void  AnalyzeEdge( AtlasMesh::CellIdentifier edgeId );

  //
  bool  LoadBalancedAnalyzeEdge( std::set< AtlasMesh::CellIdentifier >&  edges,
                                 std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies,
                                 int threadId=0 );

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



#if 0
  //
  AtlasMeshCollection::Pointer
  TryToSplit( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
              float& dataGain, float& alphasGain, float& positionGain,
              AtlasMeshCollection::Pointer& result ) const;

  //
  AtlasMeshCollection::Pointer
  TryToSwap( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
             float& dataGain, float& alphasGain, float& positionGain,
             AtlasMeshCollection::Pointer& result ) const;
#endif

  //
  AtlasMeshCollection::Pointer
  TryToRetain( const AtlasMeshCollection* innerMiniCollectionConst,
               const AtlasMeshCollection* outerMiniCollectionConst,
               AtlasMesh::CellIdentifier  edgeId,
               float& miniDataCost, float& miniAlphasCost, float& miniPositionCost );
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
  bool  OptimizeReferencePosition( AtlasMeshCollection* meshCollection,
                                   AtlasMesh::CellIdentifier  vertexId, bool optimize = true ) const;

  //
  bool  OptimizeReferencePositionFast( AtlasMeshCollection* meshCollection,
                                       AtlasMesh::CellIdentifier  vertexId, bool optimize = true ) const;

#if 0
  //
  AtlasMeshCollection::ConstPointer  OptimizeForK( const AtlasMeshCollection*  initialMeshCollection ) const;
#endif

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
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

  /** Internal structure used for passing image data into the threading library */
  struct ThreadStruct
  {
    Pointer  m_Builder;
    std::vector< AtlasMesh::CellIdentifier >  m_EdgesToTry;
  };


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


  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE FastThreaderCallback( void *arg );

  /** Internal structure used for passing image data into the threading library */
  struct FastThreadStructItem
  {
    // Input
    std::vector< AtlasMeshCollection::ConstPointer >  m_MiniCollections;
    std::vector< AtlasMesh::CellIdentifier >  m_EdgeIds;
    std::vector< AtlasMesh::PointIdentifier >  m_EdgePoint0Ids;
    std::vector< AtlasMesh::PointIdentifier >  m_EdgePoint1Ids;

    // Output
    std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >   m_DisappearingPointsLookupTable;
    AtlasMesh::PointsContainer::Pointer   m_NewReferencePosition;
    std::vector< AtlasMesh::PointsContainer::Pointer >   m_NewPositions;
    AtlasMesh::PointDataContainer::Pointer  m_NewPointParameters;
    std::set< AtlasMesh::CellIdentifier >   m_DisappearingCells;
  };

  struct FastThreadStruct
  {
    Pointer  m_Builder;
    FastThreadStructItem  m_Items[ ITK_MAX_THREADS ];
  };

private :
  AtlasMeshBuilder(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  //
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  unsigned int  m_InitialSize[ 3 ];
  unsigned int  m_NumberOfUpsamplingSteps;
  float  m_InitialStiffnesses[ 5 ];
  MultiResolutionAtlasMesher::Pointer  m_Mesher;

  float  m_PowellAbsolutePrecision;

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
#if 0
  AtlasMeshCollection::Pointer  m_SplittedMiniCollection;
  AtlasMeshCollection::Pointer  m_SwappedMiniCollection;
#endif

  float  m_RetainedCost;
  float  m_RetainedDataCost;
  float  m_RetainedAlphasCost;
  float  m_RetainedPositionCost;

  float  m_CollapsedCost;
  float  m_CollapsedDataCost;
  float  m_CollapsedAlphasCost;
  float  m_CollapsedPositionCost;

#if 0
  float  m_SplittedCost;
  float  m_SplittedDataCost;
  float  m_SplittedAlphasCost;
  float  m_SplittedPositionCost;

  float  m_SwappedCost;
  float  m_SwappedDataCost;
  float  m_SwappedAlphasCost;
  float  m_SwappedPositionCost;
#endif


#if 0
  // Stuff related to the gains container
  struct GainContainerElement
  {
    AtlasMesh::CellIdentifier  m_EdgeId;
    float  m_DataGain;
    float  m_AlphasGain;
    float  m_PositionGain;

    GainContainerElement( AtlasMesh::CellIdentifier  edgeId,
                          float dataGain,
                          float alphasGain,
                          float positionGain ) :
      m_EdgeId( edgeId ),
      m_DataGain( dataGain ),
      m_AlphasGain( alphasGain ),
      m_PositionGain( positionGain ) {}

    ~GainContainerElement() {}

  };

  //
  struct GainContainerComparer
  {
    bool operator()( const GainContainerElement& element1, const GainContainerElement& element2 ) const
    {
      float  totalGain1 = element1.m_DataGain + element1.m_AlphasGain + element1.m_PositionGain;
      float  totalGain2 = element2.m_DataGain + element2.m_AlphasGain + element2.m_PositionGain;
      if ( totalGain1 == totalGain2 )
      {
        return element1.m_EdgeId > element2.m_EdgeId;
      }
      else
      {
        return totalGain1 > totalGain2;
      }
    }

  };


  //
  class GainContainer : public std::set< GainContainerElement, GainContainerComparer >
  {
  private :
    //
    struct EdgeIdIs
    {
      inline bool operator()( const GainContainerElement& element )
      {
        return element.m_EdgeId == m_EdgeId;
      }

      EdgeIdIs( AtlasMesh::CellIdentifier edgeId ) : m_EdgeId( edgeId )
      {
      }

      AtlasMesh::CellIdentifier  m_EdgeId;

    };


  public :
    //
    bool Erase( AtlasMesh::CellIdentifier edgeId )
    {
      GainContainer::iterator pos = std::find_if( this->begin(), this->end(), EdgeIdIs( edgeId ) );
      if ( pos != this->end() )
      {
        this->erase( pos );
        return true;
      }
      else
      {
        return false;
      }

    }

    //
    void SetGain( AtlasMesh::CellIdentifier edgeId, float dataGain, float alphasGain, float positionGain )
    {
      // Make sure no element exists corresponding to the edge
      this->Erase( edgeId );

      //
      this->insert( GainContainerElement( edgeId, dataGain, alphasGain, positionGain ) );
    }

    //
    bool GetGain( AtlasMesh::CellIdentifier edgeId, float& dataGain, float& alphasGain, float& positionGain )
    {
      GainContainer::iterator pos = std::find_if( this->begin(), this->end(), EdgeIdIs( edgeId ) );
      if ( pos != this->end() )
      {
        dataGain = ( *pos ).m_DataGain;
        alphasGain = ( *pos ).m_AlphasGain;
        positionGain = ( *pos ).m_PositionGain;

        return true;
      }
      else
      {
        return false;
      }

    }

    //
    void Print( std::ostream& os ) const
    {
      os << "Gain Container: " << std::endl;

      GainContainer::const_iterator  gainIt = this->begin();
      while ( gainIt != this->end() )
      {
        os << "    edge id: " << ( *gainIt ).m_EdgeId
           << ", total gain: " << ( *gainIt ).m_DataGain +
           ( *gainIt ).m_AlphasGain +
           ( *gainIt ).m_PositionGain
           << "  (dataGain: " << ( *gainIt ).m_DataGain
           << ", alphasGain: " << ( *gainIt ).m_AlphasGain
           << ", positionGain: " << ( *gainIt ).m_PositionGain
           << ")" << std::endl;

        ++gainIt;
      }
      os << std::endl;

    }

    //
    AtlasMesh::CellIdentifier  GetBestEdge() const
    {
      if ( this->empty() )
      {
        return 0;
      }

      return ( *( this->begin() ) ).m_EdgeId;
    }


  };

  //
  GainContainer  m_Gains;

#endif


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
