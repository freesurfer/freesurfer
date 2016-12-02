#ifndef __kvlAtlasMeshLabelStatisticsCollector_h
#define __kvlAtlasMeshLabelStatisticsCollector_h

#include "kvlAtlasMeshRasterizor.h"


namespace kvl
{


namespace FragmentProcessor 
{

/**
 *
 */
class CollectLabelStatistics
{
public:
  
  typedef itk::MapContainer< AtlasMesh::PointIdentifier , AtlasAlphasType >   StatisticsContainerType;
  typedef itk::Image< AtlasAlphasType, 3 >  ProbabilityImageType;
        
  CollectLabelStatistics() 
    {
    m_ProbabilityImage = 0;
    m_SourcePointer = 0;
    m_Mesh = 0;
    m_LabelStatistics = 0;
    m_StatisticsInVertex0 = 0;
    m_StatisticsInVertex1 = 0;
    m_StatisticsInVertex2 = 0;
    m_StatisticsInVertex3 = 0;
    m_MinLogLikelihood = 0.0f;
    m_UseProbabilityImage = false;
    m_mapCompToComp=0;
    }

  ~CollectLabelStatistics() {};

  inline void operator()( const float& pi0, const float& pi1, const float& pi2,  const float& pi3 )
    {
    float  epi0 = pi0;
    float  epi1 = pi1;
    float  epi2 = pi2;
    float  epi3 = pi3;

    if ( epi0 < 0 )
      epi0 = 0.0f;
    if ( epi1 < 0 )
      epi1 = 0.0f;
    if ( epi2 < 0 )
      epi2 = 0.0f;
    if ( epi3 < 0 )
      epi3 = 0.0f;
    const float  sumOfEpis = epi0 + epi1 + epi2 + epi3;
    epi0 /= sumOfEpis;
    epi1 /= sumOfEpis;
    epi2 /= sumOfEpis;
    epi3 /= sumOfEpis;


    if ( m_UseProbabilityImage )
      {
      // Loop over all classes
      const AtlasAlphasType&  probabilities = m_ProbabilityImage->GetPixel( m_Index );
      for ( unsigned int classNumber = 0; classNumber < m_AlphasInVertex0.Size(); classNumber++ )
        {
        // Get the weight of this class's contribution
        const float  weight = probabilities[ classNumber ];

        // Collect the data terms for each vertex
        float  W0 = m_AlphasInVertex0[ classNumber ] * epi0;
        float  W1 = m_AlphasInVertex1[ classNumber ] * epi1;
        float  W2 = m_AlphasInVertex2[ classNumber ] * epi2;
        float  W3 = m_AlphasInVertex3[ classNumber ] * epi3;

        // Calculate the sum
        float  denominator = W0 + W1 + W2 + W3 + 2.2204e-16 /* 0.00001 */ /* 0.001 */;
        m_MinLogLikelihood -= weight * log( denominator );

        // Normalize to obtain W0, W1, W2, and W3
        W0 /= denominator;
        W1 /= denominator;
        W2 /= denominator;
        W3 /= denominator;
        
        // Update the histogram entries in the three vertices accordingly
        ( *m_StatisticsInVertex0 )[ classNumber ] += W0 * weight;
        ( *m_StatisticsInVertex1 )[ classNumber ] += W1 * weight;
        ( *m_StatisticsInVertex2 )[ classNumber ] += W2 * weight;
        ( *m_StatisticsInVertex3 )[ classNumber ] += W3 * weight;
      
        } // End loop over all classes

        
      }
    else
      {

      if (m_mapCompToComp==0) // old version without collapsed classes
        {
        // Calculate the numerator of the Labelification probabilities W0, W1, and W2 in this point
        float  W0 = m_AlphasInVertex0[ *m_SourcePointer ] * epi0;
        float  W1 = m_AlphasInVertex1[ *m_SourcePointer ] * epi1;
        float  W2 = m_AlphasInVertex2[ *m_SourcePointer ] * epi2;
        float  W3 = m_AlphasInVertex3[ *m_SourcePointer ] * epi3;
        
        // Calculate the sum
        float  denominator = W0 + W1 + W2 + W3 + 2.2204e-16;
        m_MinLogLikelihood -= log( denominator );
  
        // Normalize to obtain W0, W1, W2, and W3
        W0 /= denominator;
        W1 /= denominator;
        W2 /= denominator;
        W3 /= denominator;
        
        // Update the histogram entries in the three vertices accordingly
        ( *m_StatisticsInVertex0 )[ *m_SourcePointer ] += W0;
        ( *m_StatisticsInVertex1 )[ *m_SourcePointer ] += W1;
        ( *m_StatisticsInVertex2 )[ *m_SourcePointer ] += W2;
        ( *m_StatisticsInVertex3 )[ *m_SourcePointer ] += W3;
        }
      else
        {
        // This is the version with collapsed classes
        float denominator = 2.2204e-16;
        for(int ind=0; ind<m_mapCompToComp[*m_SourcePointer].size(); ind++)
          {
          unsigned char k = m_mapCompToComp[*m_SourcePointer][ind];
  
// // std::cout  << " kamlsc " << ((int) *m_SourcePointer) << " " << ((int) k) << std::endl;

          float  W0 = m_AlphasInVertex0[ k ] * epi0;
          float  W1 = m_AlphasInVertex1[ k ] * epi1;
          float  W2 = m_AlphasInVertex2[ k ] * epi2;
          float  W3 = m_AlphasInVertex3[ k ] * epi3;
  
          denominator += (W0 + W1 + W2 + W3);
          }
        m_MinLogLikelihood -= log( denominator );
  
        for(int ind=0; ind<m_mapCompToComp[*m_SourcePointer].size(); ind++)
          {
          unsigned char k = m_mapCompToComp[*m_SourcePointer][ind];
          
          ( *m_StatisticsInVertex0 )[ k ] += (m_AlphasInVertex0[ k ] * epi0 / denominator);
          ( *m_StatisticsInVertex1 )[ k ] += (m_AlphasInVertex1[ k ] * epi1 / denominator);
          ( *m_StatisticsInVertex2 )[ k ] += (m_AlphasInVertex2[ k ] * epi2 / denominator);
          ( *m_StatisticsInVertex3 )[ k ] += (m_AlphasInVertex3[ k ] * epi3 / denominator);
          }
        }

      } // End test m_UseProbabilityImage
    
    // Move on to the next pixel      
    m_SourcePointer++;
    m_Index[ 0 ]++;
    }
    
  inline void StartNewSpan( int x, int y, int z, const unsigned char* sourcePointer )
    {
    m_SourcePointer = sourcePointer;

    m_Index[ 0 ] = x;
    m_Index[ 1 ] = y;
    m_Index[ 2 ] = z;
    }
    
  inline bool StartNewTetrahedron( AtlasMesh::CellIdentifier cellId )
    {
    // Cache relevant elements of the vertices of this triangle
    AtlasMesh::CellAutoPointer  cell;
    m_Mesh->GetCell( cellId, cell );
          
    AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
    m_AlphasInVertex0 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex0 = &( m_LabelStatistics->ElementAt( *pit ) );
    ++pit;
    m_AlphasInVertex1 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex1 = &( m_LabelStatistics->ElementAt( *pit ) );
    ++pit;
    m_AlphasInVertex2 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex2 = &( m_LabelStatistics->ElementAt( *pit ) );
    ++pit;
    m_AlphasInVertex3 = m_Mesh->GetPointData()->ElementAt( *pit ).m_Alphas;
    m_StatisticsInVertex3 = &( m_LabelStatistics->ElementAt( *pit ) );


    return true;
    }
    
  inline void SetMesh( const AtlasMesh* mesh )
    {
    m_Mesh = mesh;
    
    // Create a container to hold the statistics
    unsigned int  numberOfLabeles = m_Mesh->GetPointData()->Begin().Value().m_Alphas.Size();
    AtlasAlphasType  zeroEntry( numberOfLabeles );
    zeroEntry.Fill( 0.0f );
 
    m_LabelStatistics = StatisticsContainerType::New();
    AtlasMesh::PointDataContainer::ConstIterator  it = m_Mesh->GetPointData()->Begin();
    while ( it != m_Mesh->GetPointData()->End() )
      {
      m_LabelStatistics->InsertElement( it.Index(), zeroEntry );
      
      ++it;
      }
      
    m_MinLogLikelihood = 0.0f;
    }
    
  void SetProbabilityImage( const ProbabilityImageType* probabilityImage )
    {
    m_ProbabilityImage = probabilityImage;
    }

  const ProbabilityImageType* GetProbabilityImage() const
    { return m_ProbabilityImage; }

  void  SetUseProbabilityImage( bool  useProbabilityImage )
    {
    m_UseProbabilityImage = useProbabilityImage;
    }

  bool  GetUseProbabilityImage() const
    {
    return m_UseProbabilityImage;
    }

  const StatisticsContainerType* GetLabelStatistics() const
    {
    return m_LabelStatistics;
    }
    
  float GetMinLogLikelihood() const
    {
    return m_MinLogLikelihood;
    }
    
  // Set/Get mapping of collapsed labels.
  void SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    { m_mapCompToComp = mapCompToComp; }
  std::vector<unsigned char > * GetMapCompToComp()
    { return m_mapCompToComp; }


private:

  const unsigned char*  m_SourcePointer;
  ProbabilityImageType::ConstPointer  m_ProbabilityImage;
  bool m_UseProbabilityImage;
  ProbabilityImageType::IndexType  m_Index;

  AtlasAlphasType  m_AlphasInVertex0;
  AtlasAlphasType  m_AlphasInVertex1;
  AtlasAlphasType  m_AlphasInVertex2;
  AtlasAlphasType  m_AlphasInVertex3;
  
  AtlasAlphasType*  m_StatisticsInVertex0;
  AtlasAlphasType*  m_StatisticsInVertex1;
  AtlasAlphasType*  m_StatisticsInVertex2;
  AtlasAlphasType*  m_StatisticsInVertex3;
    

  std::vector<unsigned char > *m_mapCompToComp;

  AtlasMesh::ConstPointer  m_Mesh;
  StatisticsContainerType::Pointer  m_LabelStatistics;
  
  float  m_MinLogLikelihood;

};


 
 

} // End namespace FragmentProcessor





/**
 *
 */
class AtlasMeshLabelStatisticsCollector: public AtlasMeshRasterizor< FragmentProcessor::CollectLabelStatistics >
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshLabelStatisticsCollector  Self;
  typedef AtlasMeshRasterizor< FragmentProcessor::CollectLabelStatistics >  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshLabelStatisticsCollector, itk::Object );

  /** Some typedefs */
  typedef Superclass::FragmentProcessorType  FragmentProcessorType;
  typedef Superclass::LabelImageType  LabelImageType;
  typedef FragmentProcessorType::StatisticsContainerType  StatisticsContainerType;
  typedef FragmentProcessorType::ProbabilityImageType  ProbabilityImageType;
  


  /** */
  void  SetMapCompToComp( std::vector<unsigned char > *mapCompToComp )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetMapCompToComp( mapCompToComp );  
      }
    }

  /** */
  std::vector<unsigned char > *  GetMapCompToComp()
    {
    return this->GetFragmentProcessor().GetMapCompToComp();
    }


  /** */
  void  SetProbabilityImage( const ProbabilityImageType* probabilityImage )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetProbabilityImage( probabilityImage );  
      }
    }

  /** */
  const ProbabilityImageType*  GetProbabilityImage() const
    {
    return this->GetFragmentProcessor().GetProbabilityImage();
    }
    
  /** */
  void  SetUseProbabilityImage( bool  useProbabilityImage )
    {
    for ( std::vector< FragmentProcessorType >::iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      it->SetUseProbabilityImage( useProbabilityImage );
      }
    }

  /** */
  bool  GetUseProbabilityImage() const
    {
    return this->GetFragmentProcessor().GetUseProbabilityImage();
    }

  /** */
  const StatisticsContainerType*  GetLabelStatistics() const
    { 
    if ( !m_ContributionsOfThreadsAlreadyAdded )
      {
      // Add contributions of all fragment processors to the results of the first one  
      std::vector< FragmentProcessorType >::const_iterator  it = this->GetFragmentProcessors().begin();
      ++it; // Skip first one
      for (; it != this->GetFragmentProcessors().end(); ++it )
        {
        StatisticsContainerType::ConstIterator  sourceIt = it->GetLabelStatistics()->Begin();
        StatisticsContainerType::Iterator  destIt = const_cast< StatisticsContainerType* >( this->GetFragmentProcessor().GetLabelStatistics() )->Begin();
        for (; sourceIt != it->GetLabelStatistics()->End(); ++sourceIt, ++destIt )
          {
          destIt.Value() += sourceIt.Value();  
          }
          
        }
      
      m_ContributionsOfThreadsAlreadyAdded = true;
      }

    // return result
    return this->GetFragmentProcessor().GetLabelStatistics();
    }
    
  /** */
  float GetMinLogLikelihood() const
    {
    float minLogLikelihood = 0.0f;
    for ( std::vector< FragmentProcessorType >::const_iterator  it = this->GetFragmentProcessors().begin();
          it != this->GetFragmentProcessors().end(); ++it )
      {
      minLogLikelihood += it->GetMinLogLikelihood();
      }
      
    return minLogLikelihood;
    }
  
  
    /**  */
  void Rasterize( const AtlasMesh* mesh )
    {
    m_ContributionsOfThreadsAlreadyAdded = false;
      
    // Make sure to use multi-threading
    Superclass::Rasterize( mesh, true );
    }

  
protected:
  AtlasMeshLabelStatisticsCollector() 
    {
    m_ContributionsOfThreadsAlreadyAdded = false;  
    };
  virtual ~AtlasMeshLabelStatisticsCollector() {};

  mutable bool  m_ContributionsOfThreadsAlreadyAdded;
  
private:
  AtlasMeshLabelStatisticsCollector(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
 
  
};


} // end namespace kvl

#endif

