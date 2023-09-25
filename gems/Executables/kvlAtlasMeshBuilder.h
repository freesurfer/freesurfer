#ifndef __kvlAtlasMeshBuilder_h
#define __kvlAtlasMeshBuilder_h

#if ITK_VERSION_MAJOR >= 5
#include <mutex>
#endif

#include "kvlMultiResolutionAtlasMesher.h"
#include "vnl/vnl_sample.h"
#include "itkTimeProbe.h"
#include "kvlCompressionLookupTable.h"


namespace kvl
{


#if ITK_VERSION_MAJOR >= 5
class AtlasMeshBuilderMutexLock: public std::mutex
#else
class AtlasMeshBuilderMutexLock: public itk::SimpleFastMutexLock
#endif
{
public:
  /** Standard class typedefs.  */
  typedef AtlasMeshBuilderMutexLock       Self;
#if ITK_VERSION_MAJOR >= 5
  typedef std::mutex                      Superclass;
#else  
  typedef itk::SimpleFastMutexLock        Superclass;
#endif  

  /** Lock access. */
  void DescriptiveLock( const std::string& description )
    {
#if ITK_VERSION_MAJOR >= 5
    Superclass::lock();      
#else      
    Superclass::Lock();
#endif    
    m_TimeProbe.Start();
    m_Description = description;
    }

  /** Unlock access. */
  void DescriptiveUnlock()
    {
    m_TimeProbe.Stop();
    std::cout << m_Description << ": unlocking mutex after " << m_TimeProbe.GetMean() << " seconds" << std::endl;
#if ITK_VERSION_MAJOR >= 5
    Superclass::unlock();
#else
    Superclass::Unlock();
#endif    
    }

protected:

  std::string  m_Description;
  itk::TimeProbe  m_TimeProbe;

};



//
// Let's define a helper class for multithreading that holds a mutex that is (1) automatically locked
// at construction time; (2) automatically unlocked at destruction time so that even when exceptions 
// are thrown things get wound down properly; (3) also allows for manually locking and unlocking. In
// addition, this class can also holds a shared (among threads) resource that is added and subtracted
// to in a similar exception-safe way
//
class AtlasMeshBuilderHelper
{
public:
  typedef AtlasMeshBuilderHelper Self;

  AtlasMeshBuilderHelper( AtlasMeshBuilderMutexLock&  mutex, std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies )
  : m_Mutex( mutex ), m_MutexIsLocked( true ), m_PointOccupancies( pointOccupancies )
    {
#if ITK_VERSION_MAJOR >= 5
      m_Mutex.lock();
#else      
      m_Mutex.Lock();
#endif
    }


  ~AtlasMeshBuilderHelper()
    {
    this->Lock();
    this->UnregisterPointOccupancies();
    this->Unlock();  
    }
    
  void Lock()
    {
    if ( !m_MutexIsLocked )
      {
      m_MutexIsLocked = true;
#if ITK_VERSION_MAJOR >= 5
      m_Mutex.lock();
#else      
      m_Mutex.Lock();
#endif      
      }
      
    }
    
  void Unlock(  )
    {
    if ( m_MutexIsLocked )
      {
      m_MutexIsLocked = false;
#if ITK_VERSION_MAJOR >= 5
      m_Mutex.unlock();
#else      
      m_Mutex.Unlock();
#endif      
      }
    }  
    
    
  void SetAffectedPoints( const std::vector< AtlasMesh::PointIdentifier >&  affectedPoints )  
    {
    m_AffectedPoints = affectedPoints;
    if ( !m_MutexIsLocked )
      {
      this->Lock();
      this->RegisterPointOccupancies();
      this->Unlock();
      }
    else
      {
      this->RegisterPointOccupancies();  
      }
    
    }

protected:
  AtlasMeshBuilderMutexLock&   m_Mutex;
  bool        m_MutexIsLocked;
  
  std::map< AtlasMesh::PointIdentifier, int >&  m_PointOccupancies;
  std::vector< AtlasMesh::PointIdentifier >  m_AffectedPoints;
  
  //
  void RegisterPointOccupancies()
    {
    // Indicate that we are working on them
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = m_AffectedPoints.begin();
          it != m_AffectedPoints.end(); ++it )
      {
      ( m_PointOccupancies[ *it ] )++;
      }
    }
    
  //
  void UnregisterPointOccupancies()
    {
    // Now "un-protect" the points we flagged as being worked on
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = m_AffectedPoints.begin();
          it != m_AffectedPoints.end(); ++it )
      {
      ( m_PointOccupancies[ *it ] )--;
      }
      
    }
  
private:
  AtlasMeshBuilderHelper(const Self &) ITK_DELETE_FUNCTION;
  void operator=(const Self &) ITK_DELETE_FUNCTION;

};




// Events generated
itkEventMacro( EdgeAnalysisProgressEvent, itk::UserEvent );


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
  typedef CompressionLookupTable::ImageType  LabelImageType;

  // Set label images.
  void SetUp( const std::vector< LabelImageType::ConstPointer >& labelImages,
              const CompressionLookupTable*  compressionLookupTable,
              const itk::Size< 3>&  initialSize, 
              const std::vector< double >& initialStiffnesses,
              const unsigned int maximumNumberOfIterations );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
    {
    return m_LabelImages;
    }

  //
  const AtlasMeshCollection*  GetCurrentMeshCollection() const
    {
    return m_Current;
    }


  //
  void GetCurrentDataAndAlphasCost( double& currentDataCost, double& currentAlphasCost ) const;

  //
  double GetCurrentPositionCost() const;

  //
  double GetCurrentCost() const;


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
  double  GetRetainedCosts( double& retainedDataCost, double& retainedAlphasCost, double& retainedPositionCost ) const
    {
    retainedDataCost = m_RetainedDataCost;
    retainedAlphasCost = m_RetainedAlphasCost;
    retainedPositionCost = m_RetainedPositionCost;
    return m_RetainedCost;
    }

  //
  double  GetCollapsedCosts( double& collapsedDataCost, double& collapsedAlphasCost, double& collapsedPositionCost ) const
    {
    collapsedDataCost = m_CollapsedDataCost;
    collapsedAlphasCost = m_CollapsedAlphasCost;
    collapsedPositionCost = m_CollapsedPositionCost;
    return m_CollapsedCost;
    }

  //
  void  Build( AtlasMeshCollection* explicitStartCollection = 0, double edgeCollapseEncouragementFactor = 1.0 );

  //
  double  GetProgress() const
    {
    return m_Progress;
    }

  //
  void  SetVerbose( bool verbose )
    {
    m_Verbose = verbose;  
    }  
  
  //
  bool  GetVerbose() const
    {
    return m_Verbose;
    }
    
protected :
  // Constructor
  AtlasMeshBuilder();

  // Destructor
  virtual ~AtlasMeshBuilder();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;

    
  //
  void  SetProgress( double progress )
    {
    m_Progress = progress;  
    }  
  
  //
  double  GetPositionCost( const AtlasMeshCollection* meshCollection ) const;

  //
  void  GetDataCostAndAlphasCost( const AtlasMeshCollection* meshCollection, double& dataCost, double& alphasCost ) const;

  //
  bool  LoadBalancedAnalyzeEdgeFast( std::set< AtlasMesh::CellIdentifier >&  edges,
                                     std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies,
                                     int threadId=0 );

  //
  AtlasMeshCollection::Pointer
  TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                     AtlasMesh::CellIdentifier  edgeId,
                     double& miniDataCost, double& miniAlphasCost, double& miniPositionCost,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const;

  //
  AtlasMeshCollection::Pointer
  TryToRetainFast( const AtlasMeshCollection* miniCollectionConst,
                   AtlasMesh::CellIdentifier  edgeId,
                   double& miniDataCost, double& miniAlphasCost, double& miniPositionCost );

  //
  std::vector< AtlasMesh::CellIdentifier >   GetRandomizedEdges() const;
  std::set< AtlasMesh::CellIdentifier >   GetRandomizedEdgesAsSet() const;


  //
  bool  OptimizeReferencePositionFast( AtlasMeshCollection* meshCollection,
                                       AtlasMesh::CellIdentifier  vertexId, bool optimize = true ) const;

  //
  std::vector< AtlasMesh::CellIdentifier >   Permutate(  const std::vector< AtlasMesh::CellIdentifier >&  edgeList ) const;

  //
  AtlasMeshCollection::Pointer   GetFakeCopy( const AtlasMeshCollection*  input ) const;

  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
#if ITK_VERSION_MAJOR >= 5
  static itk::ITK_THREAD_RETURN_TYPE LoadBalancedThreaderCallback( void *arg );
#else  
  static ITK_THREAD_RETURN_TYPE LoadBalancedThreaderCallback( void *arg );
#endif

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

  int m_StuckCount;

  //
  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;
  itk::Size< 3 >  m_InitialSize;
  std::vector< double >  m_InitialStiffnesses;
  MultiResolutionAtlasMesher::Pointer  m_Mesher;

  
  double  m_PowellAbsolutePrecision;

  unsigned int  m_EdgeCollapseNumber;
  unsigned int  m_MaximumNumberOfEdgeCollapses;

  AtlasMeshCollection::Pointer  m_Current;
  double  m_Progress;
  bool  m_Verbose;

  unsigned int  m_IterationNumber;
  unsigned int  m_MaximumNumberOfIterations;

  AtlasMeshCollection::Pointer  m_RetainedMiniCollection;
  AtlasMeshCollection::Pointer  m_CollapsedMiniCollection;

  double  m_RetainedCost;
  double  m_RetainedDataCost;
  double  m_RetainedAlphasCost;
  double  m_RetainedPositionCost;

  double  m_CollapsedCost;
  double  m_CollapsedDataCost;
  double  m_CollapsedAlphasCost;
  double  m_CollapsedPositionCost;
  
  double  m_EdgeCollapseEncouragementFactor;

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
