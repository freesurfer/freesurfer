#include "kvlAtlasMeshDeformationPartiallySeparableOptimizer.h"


#include "vnl/vnl_vector_fixed.h"


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
    gamma = m_StartDistance / this->ComputeMaximalDeformation( m_Gradient );
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
      vnl_vector_fixed< double, 12 >  s_mini;
      vnl_vector_fixed< double, 12 >  y_mini;
      int index = 0;
      AtlasMesh::CellType::PointIdIterator  pit = cellIt->PointIdsBegin();
      for ( int vertexNumber = 0; vertexNumber < 4; vertexNumber++ )
        {
        const AtlasPositionGradientThreadAccumType&  sElement = s->ElementAt( *pit );
        const AtlasPositionGradientThreadAccumType&  yElement = y->ElementAt( *pit );
      
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
      const vnl_vector_fixed< double, 12 >  tmp_mini = y_mini - miniApproxHessian * s_mini;
      const double  denominator = tmp_mini.transpose() * s_mini; // tmp' * s;
      const double  r = 1e-8;
      if ( abs( denominator ) > r * s.two_norm() * tmp.two_norm() ) 
        {
        // Perform the update
        miniApproxHessian += tmp_mini * tmp_mini.transpose() / denominator;
        }
  
      ++miniIt;
      } // End loop over tetrahedra  

    } // End test if first iteration
    
    
  // TODO: loop over all tetrahedra, adding each 12x12 mini Hessian to the global sparse Hessian
  // Use the Eigen C++ template library for this purpose
  
  
  
  
      
  // Part I.b: Solve for the search direction using a modified conjugate gradient algorithm
  // Use the Eigen CG gradient method as a starting point
  // [https://eigen.tuxfamily.org/dox/ConjugateGradient_8h_source.html ]
  
  
  
  
  
  
  
# if 0
  AtlasMesh::CellType::PointIdIterator  pit = mesh->GetCells()->ElementAt( tetrahedronId )->PointIdsBegin();
  const AtlasMesh::PointIdentifier  id0 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id1 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id2 = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  id3 = *pit;
  
  //AtlasMesh::PointType p0;
  //AtlasMesh::PointType p1;
  //AtlasMesh::PointType p2;
  //AtlasMesh::PointType p3;
  //mesh->GetPoint( id0, &p0 );
  //mesh->GetPoint( id1, &p1 );
  //mesh->GetPoint( id2, &p2 );
  //mesh->GetPoint( id3, &p3 );
  // Implements internally mesh->GetPoints()->GetElementIfIndexExists(ptId, point);
  // More efficient is AtlasMesh::PointType&  p0 = mesh->GetPoints()->ElementAt( id0 );
  const AtlasMesh::PointType&  p0 = mesh->GetPoints()->ElementAt( id0 );
  const AtlasMesh::PointType&  p1 = mesh->GetPoints()->ElementAt( id1 );
  const AtlasMesh::PointType&  p2 = mesh->GetPoints()->ElementAt( id2 );
  const AtlasMesh::PointType&  p3 = mesh->GetPoints()->ElementAt( id3 );
  
  ThreadAccumDataType&  priorPlusDataCost = m_ThreadSpecificMinLogLikelihoodTimesPriors[ threadNumber ];

  AtlasPositionGradientThreadAccumType&  gradientInVertex0 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id0 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex1 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id1 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex2 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id2 );
  AtlasPositionGradientThreadAccumType&  gradientInVertex3 = m_ThreadSpecificPositionGradients[ threadNumber ]->ElementAt( id3 );

#endif  
  
  
  
  
  
  
  // Direction is -r: p = -r;
  AtlasPositionGradientContainerType::Pointer  searchDirection 
                                                      = this->ScaleDeformation( r, -1.0 );

                                                      
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
}







} // end namespace kvl

