#include "kvlAtlasMeshDeformationPartiallySeparableOptimizer.h"


#include "vnl/vnl_vector.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#if 0
  #include "vnl/algo/vnl_symmetric_eigensystem.h"
#endif
  
#define KVL_ENABLE_TIME_PROBE2 1

#if KVL_ENABLE_TIME_PROBE2
  #include "itkTimeProbe.h"
#endif

  
  
namespace kvl
{


//
//
//
AtlasMeshDeformationPartiallySeparableOptimizer
::AtlasMeshDeformationPartiallySeparableOptimizer()
{

  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  m_AlphaUsedLastTime = 0.0;
  
  m_StartDistance = 1.0; // Measured in voxels

}


//
//
//
AtlasMeshDeformationPartiallySeparableOptimizer
::~AtlasMeshDeformationPartiallySeparableOptimizer()
{
}



//
//
//
void AtlasMeshDeformationPartiallySeparableOptimizer
::Initialize()
{

  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  
  Superclass::Initialize();
  


  
}

  
  
//
//
//
double
AtlasMeshDeformationPartiallySeparableOptimizer
::FindAndOptimizeNewSearchDirection()
{

  // 
  // Part I: Decide on a new search direction
  //

  
  // Part I.a: Loop over all tetrahedra, updating each tethradron's mini (approximate) 12x12 Hessian 
  // using the SR1 update, and adding its contribution to the global sparse Hessian. If this is the
  // first iteration, use Hessian = 1/gamma * I, with gamma computed based on some m_StartDistance
#if KVL_ENABLE_TIME_PROBE2
  itk::TimeProbe clock;
  clock.Start();
#endif  
  
  //  
  if ( this->GetIterationNumber() == 0 )
    { 
    // Make sure that first try of line search will be given 
    // by distance provided by user (first iteration means
    // p = -gradient, and alpha1 of line search is always 1.0 
    // for L-BFGS
    //gamma = initialAlpha1Distance / max( abs( gradient ) ); 
    const double  gamma = m_StartDistance / this->ComputeMaximalDeformation( m_Gradient );

    const int  numberOfPoints = m_Position->Size();
    int  numberOfTetrahedra = 0;  
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = this->GetMesh()->GetCells()->Begin();
        cellIt != this->GetMesh()->GetCells()->End(); ++cellIt )
      {
      if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        continue;
        }
        
      numberOfTetrahedra++;
      } // End loop over tetrahedra  
    const double  averageNumberOfTetsPerNode = static_cast< double >( numberOfTetrahedra ) / numberOfPoints;
    std::cout << "averageNumberOfTetsPerNode: " << averageNumberOfTetsPerNode << std::endl;
    
  
    const miniApproxHessianType  initialMiniApproxHessian 
                                    = ( 1/gamma ) / averageNumberOfTetsPerNode 
                                      * miniApproxHessianType().set_identity();
    std::cout << "initialMiniApproxHessian: " << initialMiniApproxHessian << std::endl;

    // m_MiniApproxHessians
    // cellIt.Index()

    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = this->GetMesh()->GetCells()->Begin();
        cellIt != this->GetMesh()->GetCells()->End(); ++cellIt )
      {
      if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        continue;
        }
        
      m_MiniApproxHessians.push_back( initialMiniApproxHessian );
      } // End loop over tetrahedra  
    
    }
  else
    {
    // Try to update miniApproxHessian using SR1 update rule
  
    // s = x - xOld; % Equivalent to alphaUsed * pOld;
    AtlasPositionGradientContainerType::ConstPointer  s  
                     = this->ScaleDeformation( m_OldSearchDirection, m_AlphaUsedLastTime ).GetPointer();
      
    // y = gradient - gradientOld;
    AtlasPositionGradientContainerType::ConstPointer  y 
                     =  this->LinearlyCombineDeformations( m_Gradient, 1.0, m_OldGradient, -1.0 ).GetPointer();

    // Loop over all tethradra
    int  numberOfTetrahedra = 0;
    int  numberOfUpdatedTetrahedra = 0;
    std::vector< miniApproxHessianType >::iterator  miniIt = m_MiniApproxHessians.begin();
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = this->GetMesh()->GetCells()->Begin();
        cellIt != this->GetMesh()->GetCells()->End(); ++cellIt )
      {
      //  
      if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
        continue;
        }
      miniApproxHessianType&  miniApproxHessian = *miniIt;
  
      // Retrieve mini s and y vectors (12-dimensional vectors, stacking 3 coordinates of each of the 4 
      // vertices) 
      vnl_vector< double >  s_mini( 12 );
      vnl_vector< double >  y_mini( 12 );
      int index = 0;
      AtlasMesh::CellType::PointIdIterator  pit = cellIt.Value()->PointIdsBegin();
      for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
        {
        const AtlasPositionGradientType&  sElement = s->ElementAt( *pit );
        const AtlasPositionGradientType&  yElement = y->ElementAt( *pit );
      
        for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )   
          {
          s_mini[ index ] = sElement[ dimensionNumber ];
          y_mini[ index ] = yElement[ dimensionNumber ];
            
          //
          ++index;  
          }
          
        ++pit;
        }
      
      if ( false )
        {
        // Actual RS1 update stuff
        // y - B * s
        const vnl_vector< double >  tmp_mini = y_mini - miniApproxHessian * s_mini;
        const double  denominator = inner_product( tmp_mini, s_mini ); // tmp' * s;
        const double  r = 1e-8;
        if ( std::abs( denominator ) > r * s_mini.two_norm() * tmp_mini.two_norm() )
          {
          // Perform the update
          miniApproxHessian += outer_product( tmp_mini, tmp_mini ) / denominator;
          // std::cout << "outer_product( tmp_mini, tmp_mini ):" << outer_product( tmp_mini, tmp_mini ) << std::endl;
          //std::cout << "denominator: " << denominator << std::endl;
          //std::cout << "s_mini.two_norm(): " << s_mini.two_norm() << std::endl;
          //std::cout << "tmp_mini.two_norm(): " << tmp_mini.two_norm() << std::endl;
          numberOfUpdatedTetrahedra++;
          }
        // else
        //   {
        //   std::cout << "denominator: " << denominator << std::endl;
        //   std::cout << "s_mini.two_norm(): " << s_mini.two_norm() << std::endl;
        //   std::cout << "tmp_mini.two_norm(): " << tmp_mini.two_norm() << std::endl;
        //   std::cout << "std::abs( denominator ): " << std::abs( denominator ) << std::endl;
        //   std::cout << "r * s_mini.two_norm() * tmp_mini.two_norm(): " << r * s_mini.two_norm() * tmp_mini.two_norm() << std::endl;
        //   }
        }
      else
        {
        // BFGS
        if ( inner_product( y_mini, s_mini ) > 1e-10 )
          {
          const vnl_vector< double >  tmp_mini = miniApproxHessian * s_mini;
          miniApproxHessian -= outer_product( tmp_mini, tmp_mini ) / inner_product( tmp_mini, s_mini );
          miniApproxHessian += outer_product( y_mini, y_mini ) / inner_product( y_mini, s_mini );
          numberOfUpdatedTetrahedra++;
          }
        }
        
      //  
      //vnl_symmetric_eigensystem< double >  xxx( miniApproxHessian.as_matrix() );
      //for ( int i = 0; i < 12; i++ )
      //  {
      //  std::cout << xxx.get_eigenvalue( i ) << std::endl;  
      //  }
      //if ( xxx.determinant() < 0 )
      //  {
      //  numberOfUpdatedTetrahedra++; 
      //  }
      
      ++numberOfTetrahedra;
      ++miniIt;
      } // End loop over tetrahedra  

    // std::cout << "numberOfUpdatedTetrahedra: " << numberOfUpdatedTetrahedra << std::endl;  
    // std::cout << "numberOfTetrahedra: " << numberOfTetrahedra << std::endl;  
    std::cout << "  Updated mini Hessians in " 
              << ( 100.0 * numberOfUpdatedTetrahedra ) / numberOfTetrahedra  
              << "% of tetrahedra" << std::endl;  
      
      
    } // End test if first iteration
  
#if KVL_ENABLE_TIME_PROBE2
  clock.Stop();
  std::cout << "  Time taken to update miniApproxHessians: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif  
  
  
    
  // Loop over all tetrahedra, adding each 12x12 mini Hessian to the global sparse Hessian
  // Use the Eigen C++ template library for this purpose
  // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
  //std::cout << "m_MiniApproxHessians.size(): " << m_MiniApproxHessians.size() << std::endl;
  //std::cout << "m_MiniApproxHessians[ 100 ]: " << m_MiniApproxHessians[ 100 ] << std::endl;
  // std::cout << "m_MiniApproxHessians[ 1000 ]: " << m_MiniApproxHessians[ 1000 ] << std::endl;
  // std::cout << "m_MiniApproxHessians[ 10000 ]: " << m_MiniApproxHessians[ 10000 ] << std::endl;
  
  
  // First make a dense mapping from each pointId to a contiguous pointNumber. This pointNumber
  // will be the contiguous index of the first element (of three) of each point
  const int  numberOfPoints = m_Position->Size();
  std::map< AtlasMesh::PointIdentifier, int >  nodeNumberLookupTable;
  for ( AtlasMesh::PointsContainer::ConstIterator  it = m_Position->Begin();
        it != m_Position->End(); ++it )
    {
    const int  counter = nodeNumberLookupTable.size();
    nodeNumberLookupTable[ it.Index() ] = counter;
    }
  //std::cout << "numberOfPoints: " <<  numberOfPoints << std::endl;
  //std::cout << "nodeNumberLookupTable.size(): " << nodeNumberLookupTable.size() << std::endl;
  //std::cout << "nodeNumberLookupTable.begin()->second: " << nodeNumberLookupTable.begin()->second << std::endl;
  
  
  // Loop over all tetrahedra, each time adding 12x12=144 entries to the triplets
  typedef Eigen::Triplet< double > TripletType;
  std::vector< TripletType >  triplets;
  std::vector< miniApproxHessianType >::const_iterator  miniIt = m_MiniApproxHessians.begin();
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = this->GetMesh()->GetCells()->Begin();
        cellIt != this->GetMesh()->GetCells()->End(); ++cellIt )
    {
    //  
    if ( cellIt.Value()->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
      {
      continue;
      }

    AtlasMesh::CellType::PointIdIterator  pit = cellIt->Value()->PointIdsBegin();
    // const std::vector< int >  nodeNumbers; 
    std::vector< int >  indicesInHessian;
    for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
      {
      // nodeNumbers.push_back( nodeNumberLookupTable[ *pit ] );
      const int  nodeNumber = nodeNumberLookupTable[ *pit ];
      for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )   
        {
        indicesInHessian.push_back( nodeNumber*3 + dimensionNumber );
      
        // if ( nodeNumber*3 + dimensionNumber >= 26770*3 )
        //   {
        //   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        //   std::cout << "nodeNumber: " << nodeNumber << std::endl;
        //   std::cout << "dimensionNumber: " << dimensionNumber << std::endl;
        //   std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;  
        //   }
      
        }  
      ++pit;
      }
    
    // Copy
    const miniApproxHessianType&  miniApproxHessian = *miniIt;
    for ( int rowNumber = 0; rowNumber < 12; rowNumber++ )
      {
      const int  rowNumberInHessian = indicesInHessian[ rowNumber ];
      for ( int columnNumber = 0; columnNumber < 12; columnNumber++ )
        {
        const int  columnNumberInHessian = indicesInHessian[ columnNumber ];

        triplets.push_back( TripletType( rowNumberInHessian, columnNumberInHessian, 
                                          miniApproxHessian( rowNumber, columnNumber ) ) );
        
        // if ( rowNumber != columnNumber )
        //   {
        //   if ( miniApproxHessian( rowNumber, columnNumber ) != 0 )
        //     {
        //     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        //     std::cout << miniApproxHessian << std::endl;  
        //     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        //     }
        //   }
          
          
        // if ( (triplets.end()-1)->row() != (triplets.end()-1)->col() )
        //   {
        //   if ( (triplets.end()-1)->value() != 0 )
        //     {
        //     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        //     std::cout << miniApproxHessian << std::endl;  
        //     std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        //     }
        //   }
          
          
        
        }
      } // End loop over 12x12 elements

    ++miniIt;
    } // End loop over tetrahedra  
  // std::cout << "triplets.size(): " << triplets.size() << std::endl;
  // std::cout << "triplets[ 400 ]: " << triplets[ 400 ].row() << ", "
  //                                   << triplets[ 400 ].col() << ", " 
  //                                   << triplets[ 400 ].value() 
  //                                   << std::endl;
  
  // int maxHessianRowNumber = 0;                                
  // int maxHessianColNumber = 0;
  // for ( std::vector< TripletType >::const_iterator it = triplets.begin();
  //       it != triplets.end(); ++ it )
  //   {
  //   if ( it->row() > maxHessianRowNumber )
  //     maxHessianRowNumber = it->row();
  //   if ( it->col() > maxHessianColNumber )
  //     maxHessianColNumber = it->col();
  //   }
  //std::cout << "maxHessianRowNumber: " << maxHessianRowNumber << std::endl;  
  //std::cout << "maxHessianColNumber: " << maxHessianColNumber << std::endl;  
#if KVL_ENABLE_TIME_PROBE2
  clock.Stop();
  std::cout << "  Time taken to construct triplets: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif  
  
  
  
  
  // Construct the Hessian from the (row,col,value) triplets, adding their contributions
  typedef Eigen::SparseMatrix< double >  HessianType;
  HessianType  Hessian( 3*numberOfPoints, 3*numberOfPoints );
  //std::cout << "Hessian.rows(): " << Hessian.rows() << std::endl;
  //std::cout << "Hessian.cols(): " << Hessian.cols() << std::endl;
  Hessian.setFromTriplets( triplets.begin(), triplets.end() );
  //std::cout << "Hessian.rows(): " << Hessian.rows() << std::endl;
  //std::cout << "Hessian.cols(): " << Hessian.cols() << std::endl;
  //std::cout << "Hessian.nonZeros(): " << Hessian.nonZeros() << std::endl;
  // for ( int k=0; k < Hessian.outerSize(); ++k )
  //   {
  //   for ( HessianType::InnerIterator it( Hessian, k ); it; ++it)
  //     {
  //     if ( it.row() != k )
  //       {
  //       if ( it.value() != 0 )
  //         {
  //         std::cout << "k, row, value: " << k << ", " << it.row() << ", " << it.value() << std::endl;
  //         }
  //       }
  //     }
  //   }
    
  
#if KVL_ENABLE_TIME_PROBE2
  clock.Stop();
  std::cout << "  Time taken to construct sparse Hessian from triplets: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif  
  
      
  // Also get the gradient in the same vectorized format
  Eigen::VectorXd  vectorizedGradient( 3 * numberOfPoints );   
  for ( AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin();
        gradIt != m_Gradient->End(); ++gradIt )
    {     
    const int  nodeNumber = nodeNumberLookupTable[ gradIt.Index() ];
    for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )   
      {
      vectorizedGradient( nodeNumber*3 + dimensionNumber ) = gradIt.Value()[ dimensionNumber ];  
      }
    }

  // std::cout << vectorizedGradient.rows() << std::endl;
  //std::cout << "  ============================" << std::endl;
  // double myMin = 1e12;
  // double myMax = -1e12;
  // for ( int i = 0; i < vectorizedGradient.rows(); i++ )
  //   { 
  //   if ( vectorizedGradient[ i ] >  myMax )
  //     myMax = vectorizedGradient[ i ];
  //   if ( vectorizedGradient[ i ] <  myMin )
  //     myMin = vectorizedGradient[ i ];   
  //   }
  // std::cout << "  myMin: " << myMin << std::endl;  
  // std::cout << "  myMax: " << myMax << std::endl;  
  
      
  // Part I.b: Solve for the search direction using a modified conjugate gradient algorithm
  // Use the Eigen CG gradient method as a starting point
  // [ https://eigen.tuxfamily.org/dox/ConjugateGradient_8h_source.html ]
  // [ https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html ]
  const double  gradientNorm = std::sqrt( vectorizedGradient.squaredNorm() );
  double  tmp = std::sqrt( gradientNorm );
  if ( tmp > 0.5 )
    {
    tmp = 0.5;  
    }
  const double  epsilon = tmp * gradientNorm;
  const double tolerance = epsilon / gradientNorm; // Reverse engineered to do the correct thing
                                                  // (threshold is epsilon^2, and also tol*tol*rhsNorm2)
  
  //Eigen::ConjugateGradientKvle< HessianType, Eigen::Lower|Eigen::Upper >  solver;
  Eigen::ConjugateGradient< HessianType, Eigen::Lower|Eigen::Upper >  solver;
  solver.compute( Hessian );
  solver.setMaxIterations( 1000 );
  solver.setTolerance( tolerance );
  Eigen::VectorXd  vectorizedSolution( 3 * numberOfPoints );
  vectorizedSolution = solver.solve( vectorizedGradient );
  std::cout << "  Number of CG iterations: " << solver.iterations() << std::endl;
  //std::cout << "  estimated error: " << solver.error()      << std::endl;

  // myMin = 1e12;
  // myMax = -1e12;
  // for ( int i = 0; i < vectorizedGradient.rows(); i++ )
  //   { 
  //   if ( vectorizedSolution[ i ] >  myMax )
  //     myMax = vectorizedSolution[ i ];
  //   if ( vectorizedSolution[ i ] <  myMin )
  //     myMin = vectorizedSolution[ i ];   
  //   }
  // std::cout << "  myMin: " << myMin << std::endl;  
  // std::cout << "  myMax: " << myMax << std::endl;  
  //std::cout << "  ============================" << std::endl;
  
  
  // Get solution back to the correct ITK format. Also, the searchDirection is
  // the *negative* of the solution
  AtlasPositionGradientContainerType::Pointer  searchDirection = AtlasPositionGradientContainerType::New();
  for ( AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin(); 
        gradIt != m_Gradient->End(); ++gradIt )
    {
    AtlasPositionGradientType  tmp;
    const int  nodeNumber = nodeNumberLookupTable[ gradIt.Index() ];
    for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )
      {
      tmp[ dimensionNumber ] = -vectorizedSolution( nodeNumber*3 + dimensionNumber );  
      }
  
    searchDirection->InsertElement( gradIt.Index(), tmp );
    }

#if KVL_ENABLE_TIME_PROBE2
  clock.Stop();
  std::cout << "  Time taken to solve linear system for search direction: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif  
  

                                                      
  //
  // PartII: Make an educated guess of the appropriate step size
  //
  const double  startAlpha = 1.0; // BFGS requires this to proof for convergence properties 
  

  
  //
  // Part III: Line Search
  // 
  const double  c1 = 1e-4;
  const double  c2 = 0.9; 
  m_OldCost = m_Cost;
  m_OldGradient = m_Gradient;
  m_OldSearchDirection = searchDirection;
  // [ x, cost, gradient, alphaUsed ] = tryLineSearch( x, cost, gradient, p, alpha1, c1, c2 );
  double  alphaUsed = 0.0;
  this->DoLineSearch( m_Position, 
                      m_Cost,
                      m_Gradient,                    
                      searchDirection,                    
                      startAlpha,
                      c1,
                      c2,
                      m_Position,
                      m_Cost,
                      m_Gradient,
                      alphaUsed );      
  
  //std::cout << "m_Cost: " << m_Cost << std::endl;


  // Some book keeping
  const double  maximalDeformation = alphaUsed * this->ComputeMaximalDeformation( searchDirection );
  m_AlphaUsedLastTime = alphaUsed;
  
#if KVL_ENABLE_TIME_PROBE2
  clock.Stop();
  std::cout << "  Time taken to do line search: " << clock.GetMean() << std::endl;
#endif  
  
  
  return maximalDeformation;
}







} // end namespace kvl

