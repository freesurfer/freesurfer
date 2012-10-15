/**
 * @file  kvlAtlasMeshDeformationLevenbergMarquardt.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
#include "kvlAtlasMeshDeformationLevenbergMarquardt.h"

#include "itkAffineTransform.h"


namespace kvl
{

#if 1

static int  staticDebugNumber = 0;

#endif

//
//
//
AtlasMeshDeformationLevenbergMarquardt
::AtlasMeshDeformationLevenbergMarquardt()
{
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

  itk::TimeProbe  timeProbe;
  timeProbe.Start();

  // Retrieve the Jacobian and the residuals from the fragment processor
  const FragmentProcessorType::GradientType*  gradient = this->GetFragmentProcessor().GetGradient();
  const FragmentProcessorType::HessianType*  hessian = this->GetFragmentProcessor().GetHessian();

  if ( verbose )
  {
    // Print some interesting stuff out about gradient and Hessian
    std::cout << "size of gradient: " << gradient->size() << std::endl;
    //gmm::row_matrix< gmm::wsvector< float > >
    std::cout << "number of rows in Hessian: " << hessian->nrows() << std::endl;
    std::cout << "number of columns in Hessian: " << hessian->ncols() << std::endl;

    std::cout << "number of non-zero entries in Hessian: " << gmm::nnz( *hessian )
              << " (" <<  static_cast< float >( gmm::nnz( *hessian ) ) /
              static_cast< float >( hessian->nrows() ) /
              static_cast< float >( hessian->ncols() ) *
              100.0f << " % of all elements)" << std::endl;

    // std::cout << "Hessian: " << std::endl;
    // for ( int i = 0; i < hessian->nrows(); i++ )
    //   {
    //   std::cout << "row " << i << " : [";
    //   for ( int j = 0; j < hessian->ncols(); j++ )
    //     {
    //     std::cout << ( *hessian )( i, j ) << " ";
    //     }
    //   std::cout << "]" << std::endl;
    //   }
    //
    // std::cout << "gradient: " << std::endl;
    // for ( int i = 0; i < hessian->nrows(); i++ )
    //   {
    //   std::cout << "row " << i << " :" << ( *gradient )[ i ] << std::endl;
    //   }

  }


  // Set up the linear system of equations ( Hessian + lambda * diag( Hessian ) ) * x = gradient
  // in a manner that GMM++ likes
  gmm::row_matrix< gmm::wsvector< double > >  lhs( hessian->ncols(), hessian->ncols() );
  gmm::copy( *hessian, lhs );
  gmm::resize( lhs, lhs.nrows()-1, lhs.ncols()-1 );  // Loose last entry that was used to put immobile components
  for ( unsigned int i = 0; i < lhs.ncols(); i++ )
  {
    lhs( i, i ) *= ( 1 + lambda );
  }


  std::vector< double >  rhs( hessian->ncols() );
  gmm::copy( *gradient, rhs );
  gmm::resize( rhs, lhs.nrows() );  // Loose last entry that was used to put immobile components

#if 0
  {
    // Save lhs and rhs to be read later in Matlab. This is useful to test Matlab's backslash
    // implementation for sparse, symmetric, positive definite matrices, which behind the scenes
    // calls the super-fast Cholesky (LU decomposition for symmetric matrices, where L is really U')
    // decomposition implementation of Tim Davids' CHOLMOD. If that Matlab test shows that it's
    // a lot faster (as I expect it will be), then we should in the future implement CHOLMOD instead
    // of gmm++. If we end up doing that, then remember NOT to fill in the lower triangle of the
    // symmetric lhs, as that won't be used anyway
    //
    // The file format (Matrix Market) can be read into Matlab using mmread.m (obtained from
    // http://math.nist.gov/MatrixMarket/mmio/matlab/mmread.m)
    //
    staticDebugNumber++;
    std::cout << "staticDebugNumber: " << staticDebugNumber << std::endl;
    std::ostringstream  lhsFileNameStream;
    lhsFileNameStream << "lhs_" << staticDebugNumber << ".txt";
    const std::string  lhsFileName = lhsFileNameStream.str();
    std::ostringstream  rhsFileNameStream;
    rhsFileNameStream << "rhs_" << staticDebugNumber << ".txt";
    const std::string  rhsFileName = rhsFileNameStream.str();


    std::cout << "Writing to " << lhsFileName << "..." << std::flush;
    //gmm::Harwell_Boeing_save( "lhs.txt", lhs );
    gmm::csc_matrix< double >  lhsCsC;
    gmm::copy( lhs, lhsCsC );
    gmm::MatrixMarket_save( lhsFileName.c_str(), lhsCsC );
    std::cout << "done." << std::endl;

    gmm::row_matrix< gmm::wsvector< double > >  tmp( lhs.nrows(), 1 );
    for ( int i = 0; i < lhs.nrows(); i++ )
    {
      tmp( i, 0 ) = rhs[ i ];
    }

    std::cout << "Writing to " << rhsFileName << "..." << std::flush;
    //gmm::Harwell_Boeing_save( "rhs.txt", tmp );
    gmm::csc_matrix< double >  rhsCsC;
    gmm::copy( tmp, rhsCsC );
    gmm::MatrixMarket_save( rhsFileName.c_str(), rhsCsC );
    std::cout << "done." << std::endl;
  }
#endif


  // Solve the linear system
  gmm::iteration  iter( 0.005 );  // relative residual
  iter.set_maxiter( 200 );  // maximum number of iterations
  if ( verbose )
  {
    iter.set_noisy( 1 );  // level of verbose ( 0 - 2 )
  }

  std::vector< double >  x( lhs.ncols() ); // allocate space for the solution vector

  gmm::csr_matrix< double >  compressedLhs;
  gmm::copy( lhs, compressedLhs );

#if 0
  gmm::identity_matrix P;  // Pre-conditioner
#else
  gmm::diagonal_precond< gmm::csr_matrix< double > >  P( compressedLhs ); // diagonal preconditioner (aka Jacobi)
#endif

  gmm::bicgstab( compressedLhs, x, rhs, P, iter ); // Do the work


  if ( verbose )
  {
    // Report how the solving went
    std::cout << "Number of iterations performed: "<< iter.get_iteration() << std::endl;
    std::cout << "Solver converged?: " << iter.converged() << std::endl;
  }


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

  timeProbe.Stop();
  std::cout << "Time taken to solve the Levenberg-Marquardt system of equations: " << timeProbe.GetMeanTime() << std::endl;

  return step;
}



} // end namespace kvl

