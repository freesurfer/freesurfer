#include "kvlAtlasMeshStatisticsCollector.h"

#include "kvlTetrahedronInteriorConstIterator.h"



namespace kvl
{

//
//
//
AtlasMeshStatisticsCollector
::AtlasMeshStatisticsCollector()
{

  m_MinLogLikelihood = 0;

}


//
//
//
AtlasMeshStatisticsCollector
::~AtlasMeshStatisticsCollector()
{
}




//
//
//
void
AtlasMeshStatisticsCollector
::Rasterize( const AtlasMesh* mesh )
{

  // Initialize from a clean slate
  m_LabelStatistics = 0;
  m_MinLogLikelihood = 0;
  m_ThreadSpecificLabelStatistics.clear();
  m_ThreadSpecificMinLogLikelihoods.clear();

  // For each thread, create an empty gradient and cost so that
  // different threads never interfere with one another
  const int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
  AtlasAlphasType  zeroEntry( numberOfClasses );
  zeroEntry.Fill( 0.0f );
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    // Initialize cost to zero for this thread
    m_ThreadSpecificMinLogLikelihoods.push_back( 0.0 );  
      
    // Create a container to hold the label statistics of this thread, and initialize to zero
    StatisticsContainerType::Pointer  labelStatistics = StatisticsContainerType::New();
    for ( AtlasMesh::PointDataContainer::ConstIterator  it = mesh->GetPointData()->Begin();
          it != mesh->GetPointData()->End(); ++it )
      {
      labelStatistics->InsertElement( it.Index(), zeroEntry );
      }
    m_ThreadSpecificLabelStatistics.push_back( labelStatistics );
    
    } // End loop over threads
    
  // Now rasterize
  Superclass::Rasterize( mesh );

  // Collect the results of all the threads
  for ( std::vector< double >::const_iterator  it = m_ThreadSpecificMinLogLikelihoods.begin();
        it != m_ThreadSpecificMinLogLikelihoods.end(); ++it )
    {
    m_MinLogLikelihood += *it;
    }  

  m_LabelStatistics = StatisticsContainerType::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  it = mesh->GetPointData()->Begin();
        it != mesh->GetPointData()->End(); ++it )
    {
    m_LabelStatistics->InsertElement( it.Index(), zeroEntry );
    }
  for ( std::vector< StatisticsContainerType::Pointer >::const_iterator  
             it = m_ThreadSpecificLabelStatistics.begin();
        it != m_ThreadSpecificLabelStatistics.end(); ++it )
    {
    StatisticsContainerType::Iterator  sourceIt = ( *it )->Begin(); 
    StatisticsContainerType::Iterator  targetIt = m_LabelStatistics->Begin(); 
    for ( ; targetIt != m_LabelStatistics->End(); ++sourceIt, ++targetIt )
      {
      targetIt.Value() += sourceIt.Value();
      } 
      
    } // End loop over all threads  
    
    
}    
  
    
  

//
//
//
bool
AtlasMeshStatisticsCollector
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{
  // Retrieve necessary info about tetrahedron
  ReferenceTetrahedronInfo  info;
  mesh->GetCellData( tetrahedronId, &info );
 
  AtlasMesh::CellAutoPointer  cell;
  mesh->GetCell( tetrahedronId, cell );

  AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  
  AtlasMesh::PointType p0;
  AtlasMesh::PointType p1;
  AtlasMesh::PointType p2;
  AtlasMesh::PointType p3;
  mesh->GetPoint( id0, &p0 );
  mesh->GetPoint( id1, &p1 );
  mesh->GetPoint( id2, &p2 );
  mesh->GetPoint( id3, &p3 );
  
  const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
  const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
  
  
  // Compute actual contribution of tetrahedron
  AtlasAlphasType  statisticsInVertex0;
  AtlasAlphasType  statisticsInVertex1;
  AtlasAlphasType  statisticsInVertex2;
  AtlasAlphasType  statisticsInVertex3;
  double  minLogLikelihood = 0.0;
  this->GetContributionOfTetrahedron( p0, p1, p2, p3,
                                      alphasInVertex0, 
                                      alphasInVertex1,
                                      alphasInVertex2,
                                      alphasInVertex3,
                                      minLogLikelihood,
                                      statisticsInVertex0,
                                      statisticsInVertex1,
                                      statisticsInVertex2,
                                      statisticsInVertex3 );
 
  
  // Add contribution of this tetrahedron to thread's results
  m_ThreadSpecificMinLogLikelihoods[ threadNumber ] += minLogLikelihood;

  m_ThreadSpecificLabelStatistics[ threadNumber ]->ElementAt( id0 ) += statisticsInVertex0;
  m_ThreadSpecificLabelStatistics[ threadNumber ]->ElementAt( id1 ) += statisticsInVertex1;
  m_ThreadSpecificLabelStatistics[ threadNumber ]->ElementAt( id2 ) += statisticsInVertex2;
  m_ThreadSpecificLabelStatistics[ threadNumber ]->ElementAt( id3 ) += statisticsInVertex3;

  return true;
}






} // end namespace kvl
