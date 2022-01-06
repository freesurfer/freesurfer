#include "kvlAtlasMeshDeformationPartiallySeparableOptimizer.h"


#include "vnl/vnl_vector.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"


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
  
  
  //  
  if ( this->GetIterationNumber() == 0 )
    { 
    // Make sure that first try of line search will be given 
    // by distance provided by user (first iteration means
    // p = -gradient, and alpha1 of line search is always 1.0 
    // for L-BFGS
    //gamma = initialAlpha1Distance / max( abs( gradient ) ); 
    const double  gamma = m_StartDistance / this->ComputeMaximalDeformation( m_Gradient );
    const miniApproxHessianType  initialMiniApproxHessian = 1/gamma * miniApproxHessianType().set_identity();
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
      
      // Actual RS1 update stuff
      // y - B * s
      const vnl_vector< double >  tmp_mini = y_mini - miniApproxHessian * s_mini;
      const double  denominator = inner_product( tmp_mini, s_mini ); // tmp' * s;
      const double  r = 1e-8;
      if ( abs( denominator ) > r * s_mini.two_norm() * tmp_mini.two_norm() ) 
        {
        // Perform the update
        miniApproxHessian += outer_product( tmp_mini, tmp_mini ) / denominator;
        }
  
      ++miniIt;
      } // End loop over tetrahedra  

    } // End test if first iteration
    
#if 0      
    
  // Loop over all tetrahedra, adding each 12x12 mini Hessian to the global sparse Hessian
  // Use the Eigen C++ template library for this purpose
  // https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
  
  // First make a dense mapping from each pointId to a contiguous pointNumber. This pointNumber
  // will be the contiguous index of the first element (of three) of each point
  const int  numberOfPoints = m_Position->Size();
  std::map< AtlasMesh::PointIdentifier, int >  nodeNumberLookupTable;
  for ( AtlasMesh::PointsContainer::ConstIterator  it = m_Position->Begin();
        it != m_Position->End(); ++it )
    {
    nodeNumberLookupTable[ it.Index() ] = nodeNumberLookupTable.size();
    }
  
  
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

    AtlasMesh::CellType::PointIdIterator  pit = cellIt->PointIdsBegin();
    // const std::vector< int >  nodeNumbers; 
    std::vector< int >  indicesInHessian;
    for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
      {
      // nodeNumbers.push_back( nodeNumberLookupTable[ *pit ] );
      nodeNumber = nodeNumberLookupTable[ *pit ];
      for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )   
        {
        indicesInHessian.push_back( nodeNumber*3 + dimensionNumber )
        }  
      ++pit;
      }
    
    // Copy
    const miniApproxHessianType&  miniApproxHessian = *miniIt;
    for ( int rowNumber = 0; rowNumber < 12; rowNumber++ )
      {
      const  rowNumberInHessian = indicesInHessian[ rowNumber ];
      for ( int columnNumber = 0; columnNumber < 12; columnNumber++ )
        {
        const  columnNumberInHessian = indicesInHessian[ columnNumber ];

        triplets.push_back( TripletType( rowNumberInHessian, columnNumberInHessian, 
                                          miniApproxHessian[ rowNumber, columnNumber ] ) );
        }
      } // End loop over 12x12 elements

    ++miniIt;
    } // End loop over tetrahedra  

  
  // Construct the Hessian from the (row,col,value) triplets, adding their contributions
  typedef Eigen::SparseMatrix< double >  HessianType;
  HessianType  Hessian( 3*numberOfPoints, 3*numberOfPoints );
  Hessian.setFromTriplets( triplets.begin(), triplets.end() );
  
  
  // Also get the gradient in the same vectorized format
  Eigen::VectorXd  vectorizedGradient( 3 * numberOfPoints );   
  for ( AtlasPositionGradientContainerType::ConstIterator  gradIt = gradient->Begin();
        gradIt != gradient->End(); ++gradIt )
    {     
    const int  nodeNumber = nodeNumberLookupTable[ gradIt.Index() ];
    for ( int dimensionNumber = 0; dimensionNumber < 3; dimensionNumber++ )   
      {
      vectorizedGradient( nodeNumber*3 + dimensionNumber ) = gradIt.Value()[ dimensionNumber ];  
      }
    }
  
  
  
      
  // Part I.b: Solve for the search direction using a modified conjugate gradient algorithm
  // Use the Eigen CG gradient method as a starting point
  // [ https://eigen.tuxfamily.org/dox/ConjugateGradient_8h_source.html ]
  // [ https://eigen.tuxfamily.org/dox/classEigen_1_1ConjugateGradient.html ]
  Eigen::ConjugateGradient< HessianType, Eigen::Lower|Eigen::Upper >  solver;
  solver.compute( Hessian );
  solver.setMaxIterations( 10 );
  const Eigen::VectorXd  vectorizedSolution( 3 * numberOfPoints );
  vectorizedSolution = solver.solve( vectorizedGradient );
  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  
  
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
      tmp[ dimensionNumber ] = vectorizedSolution( nodeNumber*3 + dimensionNumber );  
      }
  
    searchDirection->InsertElement( gradIt.Index(), tmp );
    }

  

                                                      
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
  
  return maximalDeformation;
#endif  
}







} // end namespace kvl

