#include "kvlAtlasMeshRasterizor.h"

#if ITK_VERSION_MAJOR >= 5
#include <itkMultiThreaderBase.h>
#include <mutex>
static std::mutex rasterizorMutex;
#else
static itk::SimpleFastMutexLock rasterizorMutex;
#endif




namespace kvl
{

//
//
//
AtlasMeshRasterizor
::AtlasMeshRasterizor()
{
#if ITK_VERSION_MAJOR >= 5
  m_NumberOfThreads = itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads();
#else  
  m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
#endif  
}



//
//
//
void
AtlasMeshRasterizor
::Rasterize( const AtlasMesh* mesh )
{

  // Fill in the data structure to pass on to the threads
  ThreadStruct  str;
  str.m_Rasterizor = this;
  str.m_Mesh = mesh;
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = mesh->GetCells()->Begin();
        cellIt != mesh->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() == AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      str.m_TetrahedronIds.push_back( cellIt.Index() );
      //str.m_TetrahedronIds.insert( cellIt.Index() );
      }
    }

  // Set up the multithreader
#if ITK_VERSION_MAJOR >= 5
  itk::MultiThreaderBase::Pointer  threader = itk::MultiThreaderBase::New();
  threader->SetNumberOfWorkUnits( this->GetNumberOfThreads() );
#else
  itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
  threader->SetNumberOfThreads( this->GetNumberOfThreads() );
  //threader->SetNumberOfThreads( 1 );
#endif  
  threader->SetSingleMethod( this->ThreaderCallback, &str );

  // Let the beast go
  threader->SingleMethodExecute();

  
}




//
//
//
#if ITK_VERSION_MAJOR >= 5
itk::ITK_THREAD_RETURN_TYPE
AtlasMeshRasterizor
::ThreaderCallback( void *arg )
#else  
ITK_THREAD_RETURN_TYPE
AtlasMeshRasterizor
::ThreaderCallback( void *arg )
#endif  
{

  // Retrieve the input arguments
#if ITK_VERSION_MAJOR >= 5
  const int  threadNumber = ((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->WorkUnitID;
  const int  numberOfThreads = ((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->NumberOfWorkUnits;
  ThreadStruct*  str = (ThreadStruct *)(((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->UserData);
#else  
  const int  threadNumber = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  const int  numberOfThreads = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;
  ThreadStruct*  str = (ThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);
#endif  

  
#if 1  
  // Compute up-front which tetrahedra this thread should be responsible for. This isn't a 
  // particularly good way of load-balancing, but it does allow us to get the exact same
  // round-off errors (by adding many floating-point contributions) every single time we
  // repeat the same computation on the same computer with the same number of threads.
  const int  numberOfTetrahedra = str->m_TetrahedronIds.size();
#if 0  
  const int  numberOfTetrahedraPerThread = 
          itk::Math::Ceil< int >( static_cast< double >( numberOfTetrahedra ) / numberOfThreads );
  const int  thisThreadStartNumber =  threadNumber * numberOfTetrahedraPerThread;
  int  thisThreadEndNumber =  thisThreadStartNumber + numberOfTetrahedraPerThread - 1;
  if ( thisThreadEndNumber > (numberOfTetrahedra-1) )
    {
    thisThreadEndNumber = (numberOfTetrahedra-1);
    }
    
  // Rasterize all tetrahedra assigned to this thread  
  for ( int tetrahedronNumber = thisThreadStartNumber; 
        tetrahedronNumber <= thisThreadEndNumber; 
        tetrahedronNumber++ )
#else
  // Rasterize all tetrahedra assigned to this thread  
  for ( int tetrahedronNumber = threadNumber; 
        tetrahedronNumber < numberOfTetrahedra; 
        tetrahedronNumber += numberOfThreads )
#endif
    {
    if ( !str->m_Rasterizor->RasterizeTetrahedron( str->m_Mesh, 
                                                   str->m_TetrahedronIds[ tetrahedronNumber ],
                                                   threadNumber ) )
      {
      // Something wrong with this tetrahedron; abort at least this thread
      break;
      }  
      
    }
#else

  while ( true )
    { 
    rasterizorMutex.Lock();

    // Check if there is anything left to do
    if ( str->m_TetrahedronIds.size() == 0 )
      {
      rasterizorMutex.Unlock();
#if ITK_VERSION_MAJOR >= 5
      return itk::ITK_THREAD_RETURN_DEFAULT_VALUE;
#else        
      return ITK_THREAD_RETURN_VALUE;
#endif      
      }

    // Let's define how many tetrahedra this thread is going to take on 
    const int  numberOfTetrahedra = 100;

    // Pop the first numberOfTetrahedra tetrahedra on the list
    std::vector< AtlasMesh::CellIdentifier >  thisThreadTetrahedronIds;
    for ( int i = 0; i < numberOfTetrahedra; i++ )
      {
      const AtlasMesh::CellIdentifier  tetrahedronId = *( str->m_TetrahedronIds.begin() );
      str->m_TetrahedronIds.erase( tetrahedronId );
      thisThreadTetrahedronIds.push_back( tetrahedronId );
      if ( str->m_TetrahedronIds.size() == 0 )
        {
        break;
        }
      }

    rasterizorMutex.Unlock();


    
    for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = thisThreadTetrahedronIds.begin();
          it != thisThreadTetrahedronIds.end(); ++it )
      {
      // 
      if ( !str->m_Rasterizor->RasterizeTetrahedron( str->m_Mesh, *it ) )
        {
        // Something wrong with this tetrahedron; abort this thread and 
        // make sure other threads also stop ASAP
        rasterizorMutex.Lock();
        str->m_TetrahedronIds.clear();
        rasterizorMutex.Unlock();
          
#if ITK_VERSION_MAJOR >= 5
        return itk::ITK_THREAD_RETURN_DEFAULT_VALUE;
#else  	
        return ITK_THREAD_RETURN_VALUE;
#endif	
        }
        
      }  

    } // End infinite loop   
      
#endif
    
  
#if ITK_VERSION_MAJOR >= 5
  return itk::ITK_THREAD_RETURN_DEFAULT_VALUE;
#else  
  return ITK_THREAD_RETURN_VALUE;
#endif  
}





} // end namespace kvl

