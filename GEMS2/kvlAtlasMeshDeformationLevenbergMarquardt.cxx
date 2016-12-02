#include "kvlAtlasMeshDeformationLevenbergMarquardt.h"

#include "itkAffineTransform.h"


namespace kvl
{

static itk::SimpleFastMutexLock levenbergMarquardtMutex;

//
//
//
AtlasMeshDeformationLevenbergMarquardt
::AtlasMeshDeformationLevenbergMarquardt()
{
 
  m_MaximumNumberOfCGIterations = 200;
  m_OldLambda = 0.0;
  m_PeerToCopyHessianFrom = 0;
}


//
//
//
AtlasMeshDeformationLevenbergMarquardt
::~AtlasMeshDeformationLevenbergMarquardt()
{
}



//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasMeshDeformationLevenbergMarquardt
::GetStep( float lambda, bool verbose ) const
{

//  itk::TimeProbe  timeProbe;
//  timeProbe.Start();

  if ( verbose )
    {
    // Print some interesting stuff out about gradient and Hessian
//    std::cout << "size of gradient: " << this->GetFragmentProcessor().GetGradient().size() << std::endl;
    //gmm::row_matrix< gmm::wsvector< float > >
//    std::cout << "number of rows in Hessian: " << this->GetFragmentProcessor().GetHessian().nrows() << std::endl;
//    std::cout << "number of columns in Hessian: " << this->GetFragmentProcessor().GetHessian().ncols() << std::endl;

//    std::cout << "number of non-zero entries in Hessian: " << gmm::nnz( this->GetFragmentProcessor().GetHessian() )
//              << " (" <<  static_cast< float >( gmm::nnz( this->GetFragmentProcessor().GetHessian() ) ) /
//                          static_cast< float >( this->GetFragmentProcessor().GetHessian().nrows() ) /
//                          static_cast< float >( this->GetFragmentProcessor().GetHessian().ncols() ) *
//                          100.0f << " % of all elements)" << std::endl;

    }


  // Multiply the diagonal of the Hessian by ( 1 + lambda ). For efficiency purposes, we don't store
  // the original Hessian, but rather the Hessian with diagonal multiplied by ( 1 + oldLambda ), where
  // oldLambda is the lambda of the last time the current function was called. 
  const double  correctionFactor = ( 1 + lambda ) / ( 1 + m_OldLambda );
    for ( unsigned int i = 0; i < this->GetFragmentProcessor().GetHessian().ncols(); i++ )
    {
    const_cast< HessianType& >( this->GetFragmentProcessor().GetHessian() )( i, i ) *= correctionFactor;
    }
  m_OldLambda = lambda;

#if 0
  // Set up the linear system; rmemember to loose the last entry as that was used to put immobile components
  gmm::iteration  iter( 0.005 );  // relative residual
  iter.set_maxiter( m_MaximumNumberOfCGIterations );  // maximum number of iterations
  if ( verbose )
    {
    iter.set_noisy( 1 );  // level of verbose ( 0 - 2 )
    }

  std::vector< double >  x( gmm::vect_size( this->GetFragmentProcessor().GetGradient() ) - 1 ); // allocate space for the solution vector

  gmm::csr_matrix< double >  compressedLhs;
  gmm::copy( gmm::sub_matrix( this->GetFragmentProcessor().GetHessian(), 
                              gmm::sub_interval( 0, gmm::mat_nrows( this->GetFragmentProcessor().GetHessian() ) - 1 ), 
                              gmm::sub_interval( 0, gmm::mat_ncols( this->GetFragmentProcessor().GetHessian() ) - 1 ) ), 
             compressedLhs );

  gmm::diagonal_precond< gmm::csr_matrix< double > >  P( compressedLhs ); // diagonal preconditioner (aka Jacobi)

  // Do the work
  gmm::bicgstab( compressedLhs, x, 
                 gmm::sub_vector( this->GetFragmentProcessor().GetGradient(), 
                                  gmm::sub_interval( 0, gmm::vect_size( this->GetFragmentProcessor().GetGradient() ) - 1 ) ), 
                 P, iter );
#else

  // Set up the linear system; rememember to loose the last entry as that was used to put immobile components. Since
  // resizing and/or copying matrices in GMM seems to be inefficient, let's just fill the last row of the Hessian
  // with zeros (except for the diagonal element) and the last entry of the gradient with zero: the solution is then
  // the same is that of the "shrunk" system
  
  
  // Zero out the last row
  GradientType&  gradient = const_cast< GradientType& >( this->GetFragmentProcessor().GetGradient() );
  HessianType&  hessian = const_cast< HessianType& >( this->GetFragmentProcessor().GetHessian() );
  const gmm::size_type  lastRowNumber = gmm::vect_size( gradient ) - 1;
  for ( gmm::linalg_traits< gmm::wsvector< double > >::iterator  it = gmm::vect_begin( gmm::mat_row( hessian, lastRowNumber ) );
        it != gmm::vect_end( gmm::mat_row( hessian, lastRowNumber ) ); ++it )
    {
    *it = 1e-12; 

#if 1
    // Let's also zero out the corresponding column, to keep everything symmetric (not sure if solver uses this)
    hessian( it.index(), lastRowNumber ) = 1e-12;
#endif    
  
    }
  gradient[ lastRowNumber ] = 0.0;  
    
  
  gmm::iteration  iter( 0.005 );  // relative residual
  iter.set_maxiter( m_MaximumNumberOfCGIterations );  // maximum number of iterations
  if ( verbose )
    {
    iter.set_noisy( 1 );  // level of verbose ( 0 - 2 )
    }

  std::vector< double >  x( gmm::vect_size( this->GetFragmentProcessor().GetGradient() ) ); // allocate space for the solution vector

  gmm::csr_matrix< double >  compressedLhs;
  gmm::copy( this->GetFragmentProcessor().GetHessian(), compressedLhs );

  gmm::diagonal_precond< gmm::csr_matrix< double > >  P( compressedLhs ); // diagonal preconditioner (aka Jacobi)

  // Do the work
  gmm::bicgstab( compressedLhs, x, 
                 this->GetFragmentProcessor().GetGradient(), 
                 P, iter );
  
  // Loose the last component
  gmm::resize( x, gmm::vect_size( x ) - 1 );

#endif

  if ( verbose )
    {
    // Report how the solving went
    std::cout << "Number of iterations performed: "<< iter.get_iteration() << std::endl;
    std::cout << "Solver converged?: " << iter.converged() << std::endl;
    }

//  timeProbe.Stop();
//  std::cout << "Time taken to solve system of equations: " << timeProbe.GetMeanTime() << std::endl;
//  timeProbe = itk::TimeProbe();
//  timeProbe.Start();

  // Convert the solution into the correct format to return
  const std::map< AtlasMesh::PointIdentifier, int >&  componentXLookupTable =
                                                            this->GetFragmentProcessor().GetComponentXLookupTable();
  const std::map< AtlasMesh::PointIdentifier, int >&  componentYLookupTable =
                                                            this->GetFragmentProcessor().GetComponentYLookupTable();
  const std::map< AtlasMesh::PointIdentifier, int >&  componentZLookupTable =
                                                            this->GetFragmentProcessor().GetComponentZLookupTable();
  AtlasPositionGradientContainerType::Pointer  step = AtlasPositionGradientContainerType::New();
  for ( AtlasMesh::PointsContainer::ConstIterator pointIt = this->GetFragmentProcessor().GetMesh()->GetPoints()->Begin();
        pointIt != this->GetFragmentProcessor().GetMesh()->GetPoints()->End(); ++pointIt )
    {
    const int componentXIndex = ( componentXLookupTable.find( pointIt.Index() ) )->second;
    const int componentYIndex = ( componentYLookupTable.find( pointIt.Index() ) )->second;
    const int componentZIndex = ( componentZLookupTable.find( pointIt.Index() ) )->second;

    AtlasPositionGradientType  entry;

    if ( componentXIndex < static_cast< int >( x.size() ) )
      {
      entry[ 0 ] = -x[ componentXIndex ];
      }
    else
      {
      entry[ 0 ] = 0; // Immobile points map outside solution vector
      }

    if ( componentYIndex < static_cast< int >( x.size() ) )
      {
      entry[ 1 ] = -x[ componentYIndex ];
      }
    else
      {
      entry[ 1 ] = 0; // Immobile points map outside solution vector
      }

    if ( componentZIndex < static_cast< int >( x.size() ) )
      {
      entry[ 2 ] = -x[ componentZIndex ];
      }
    else
      {
      entry[ 2 ] = 0; // Immobile points map outside solution vector
      }


    // Remember that we now have the step in the coordinate system parallell to the mesh. In order
    // to go back into image grid coordinate system, we need to multiply by T (from  x = T * u )
    if ( this->GetMeshToImageTransform() )
      {
      entry = ( this->GetMeshToImageTransform()->GetMatrix() ) * entry;
      }


    // std::cout << "Gradient in point with index " << pointIt.Index() << ": [ "
    //          << entry[ 0 ] << "   " << entry[ 1 ] << "   " << entry[ 2 ] << "]" << std::endl;

    step->InsertElement( pointIt.Index(), entry );
    }

//  timeProbe.Stop();
//  std::cout << "Time taken to map results back into mesh position step: " << timeProbe.GetMeanTime() << std::endl;

  return step;

}



//
//
//
void
AtlasMeshDeformationLevenbergMarquardt
::Rasterize( const AtlasMesh* mesh )
{
  // Rasterize
//  itk::TimeProbe  timeProbe;
//  timeProbe.Start();

  Superclass::Rasterize( mesh, true );
//  timeProbe.Stop();
 // std::cout << "Time taken to rasterize Levenberg-Marquardt: " << timeProbe.GetMeanTime() << std::endl;
    
  // Add gradient (non-threaded)
//  timeProbe = itk::TimeProbe();
//  timeProbe.Start();
  std::vector< FragmentProcessorType >::const_iterator it = this->GetFragmentProcessors().begin();
  GradientType&  totalGradient = const_cast< GradientType& >( it->GetGradient() );
  ++it;
  for ( ; it != this->GetFragmentProcessors().end(); ++it )
    {
    gmm::add( it->GetGradient(), totalGradient );
    }
//  timeProbe.Stop();
//  std::cout << "Time taken to add gradient: " << timeProbe.GetMeanTime() << std::endl;

  
  // Add hessian using multi-threading
  if ( !m_PeerToCopyHessianFrom )
    {
//    timeProbe = itk::TimeProbe();
//    timeProbe.Start();
    ThreadStruct  str;
    str.m_FragmentProcessors = &( this->GetFragmentProcessors() );
    for ( gmm::size_type rowNumber = 0; rowNumber < gmm::vect_size( totalGradient ); ++rowNumber )
      {
      str.m_RowNumbersToAdd.insert( rowNumber );  
      }
    
    // Set up the multithreader
    itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
    threader->SetSingleMethod( this->ThreaderCallback, &str );

    // Let the beast go
    threader->SingleMethodExecute();

//    timeProbe.Stop();
 //   std::cout << "Time taken to add hessian: " << timeProbe.GetMeanTime() << std::endl;

    // At this point m_Lhs is really the original Hessian
    m_OldLambda = 0.0;
    }
  else
    {
    // Make a copy of the Hessian from our peer
    const HessianType&  sourceHessian = m_PeerToCopyHessianFrom->GetFragmentProcessor().GetHessian();
    HessianType&  targetHessian = const_cast< HessianType& >( this->GetFragmentProcessor().GetHessian() );
    for ( unsigned int rowNumber = 0; rowNumber < this->GetFragmentProcessor().GetHessian().nrows(); rowNumber++ )
      {
      gmm::linalg_traits< gmm::wsvector< double > >::const_iterator  sourceIt = gmm::vect_const_begin( gmm::mat_const_row( sourceHessian, rowNumber ) );
      gmm::linalg_traits< gmm::wsvector< double > >::iterator  targetIt = gmm::vect_begin( gmm::mat_row( targetHessian, rowNumber ) );
      for ( ; sourceIt != gmm::vect_const_end( gmm::mat_const_row( sourceHessian, rowNumber ) ); ++sourceIt, ++targetIt )
        {
        //std::cout << "              target: " << targetIt.index() << " -> " << *targetIt << std::endl;  
        //std::cout << "              source: " << sourceIt.index() << " -> " << *sourceIt << std::endl;  
        *targetIt = *sourceIt;  
        }
      } // End loop over all rows


    // Also make sure to inform ourselves that we've only a tampered-with Hessian
    m_OldLambda = m_PeerToCopyHessianFrom->m_OldLambda;  
    }
    
  
}



//
//
//
ITK_THREAD_RETURN_TYPE
AtlasMeshDeformationLevenbergMarquardt
::ThreaderCallback( void *arg )
{

  // Retrieve the input arguments
  const int  threadId = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  //const int  threadCount = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  ThreadStruct*  str = (ThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  //std::cout << "threadId: " << threadId << std::endl;

  const int  numberOfRowsAtATime = 10; // Could do something more intelligent here...
  int  numberOfRowsAdded = 0;
  while ( true )
    {
    //std::cout << "I am thread " << threadId << " and I have added "
    //          <<  numberOfRowsAdded << " rows so far" << std::endl;

    if ( !Self::AddHessianRows( str->m_FragmentProcessors, 
                                str->m_RowNumbersToAdd, 
                                numberOfRowsAtATime, threadId ) )
      {
      //std::cout << "Nothing left to do for thread " << threadId << std::endl;
      break;
      }

    numberOfRowsAdded += numberOfRowsAtATime;
    }

  return ITK_THREAD_RETURN_VALUE;
}




//
//
//


bool
AtlasMeshDeformationLevenbergMarquardt
::AddHessianRows( std::vector< FragmentProcessorType >*  fragmentProcessors, 
                  std::set< gmm::size_type >&  rowNumbersToAdd, 
                  int numberOfRows, int threadId )
{
  levenbergMarquardtMutex.Lock();


  // Check if there is anything left to do
  if ( rowNumbersToAdd.size() == 0 )
    {
    levenbergMarquardtMutex.Unlock();
    return false;
    }


  // Pop the first numberOfRows on the list
  std::vector< gmm::size_type >  rowNumbers;
  for ( int i = 0; i < numberOfRows; i++ )
    {
    const gmm::size_type rowNumber = *( rowNumbersToAdd.begin() );
    rowNumbersToAdd.erase( rowNumber );

    rowNumbers.push_back( rowNumber );
    if ( rowNumbersToAdd.size() == 0 )
      {
      break;
      }
    }

  levenbergMarquardtMutex.Unlock();

  // Loop over all the row numbers to be handled
  for ( std::vector< gmm::size_type >::const_iterator  it = rowNumbers.begin();
        it != rowNumbers.end(); ++it )
    {
    const gmm::size_type  rowNumber = *it;

  
    // Loop over all the fragment processors
    std::vector< FragmentProcessorType >::const_iterator itt = fragmentProcessors->begin();
    HessianType&  totalHessian = const_cast< HessianType& >( itt->GetHessian() );
    ++itt;
    for ( ; itt != fragmentProcessors->end(); ++itt )
      {
      // Add this row
      gmm::linalg_traits< gmm::wsvector< double > >::iterator  destIt = gmm::vect_begin( gmm::mat_row( totalHessian, rowNumber ) );
      gmm::linalg_traits< gmm::wsvector< double > >::const_iterator  sourceIt = gmm::vect_const_begin( gmm::mat_const_row( itt->GetHessian(), rowNumber ) );
      for ( ; destIt != gmm::vect_end( gmm::mat_row( totalHessian, rowNumber ) ); ++destIt, ++sourceIt )
        {
        //std::cout << "              dest: " << destIt.index() << " -> " << *destIt << std::endl;  
        //std::cout << "              source: " << sourceIt.index() << " -> " << *sourceIt << std::endl;  
        *destIt += *sourceIt;  
        }
      } // End loop over all fragment processors

    } // End loop over all row numbers

  return true;

}




} // end namespace kvl

