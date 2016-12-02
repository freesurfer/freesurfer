#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationConjugateGradientOptimizer
::AtlasMeshDeformationConjugateGradientOptimizer()
{
  m_MaximalDeformationStopCriterion = /* 0.001 */ 0.05;
  
  m_Cost = 0;
  m_OldCost = 0;
  m_Position = 0;
  m_Gradient = 0;
  m_OldGradient = 0;
  m_Direction = 0;
  //m_StepSize = 0;
}


//
//
//
AtlasMeshDeformationConjugateGradientOptimizer
::~AtlasMeshDeformationConjugateGradientOptimizer()
{
}




//
//
//
bool
AtlasMeshDeformationConjugateGradientOptimizer
::Go()
{

  this->InvokeEvent( DeformationStartEvent() );

  //// std::cout << "m_MaximumNumberOfIterations: " << m_MaximumNumberOfIterations << std::endl;

  
  bool haveMoved = false; // Keep track if we've ever moved or not
  for ( int itn = 0; 
        itn < m_MaximumNumberOfIterations;
        itn++ )
    {

m_IterationNumber = itn;

    // Try to make one successful move
    const float maximalDeformation = this->PerformOneIteration();
    if ( maximalDeformation > 0 )
      {
      haveMoved = true;
      }
      
    // Test if we need to stop
    if ( maximalDeformation <= m_MaximalDeformationStopCriterion )
      {
      if ( ( itn > 0 ) || ( maximalDeformation == 0 ) )
        {
        if ( m_Verbose )
          {
          // std::cout << "        maximalDeformation is too small; stopping" << std::endl
          }
        break;
        }
      }

    if ( !( itn % m_IterationEventResolution ) )
      {
      this->InvokeEvent( DeformationIterationEvent() );
      }
      

    } // End loop over all positive moves


  this->InvokeEvent( DeformationEndEvent() );

  return haveMoved;
 
}




//
//
//
void
AtlasMeshDeformationConjugateGradientOptimizer
::GetCostAndGradient( const AtlasMesh::PointsContainer* position, double& cost, AtlasPositionGradientContainerType::Pointer& gradient )
{
  
  // Create an empty, dummy template image needed by anything that rasterizes meshes
  AtlasMeshToIntensityImageGradientCalculator::LabelImageType::Pointer  dummyTemplateImage = AtlasMeshToIntensityImageGradientCalculator::LabelImageType::New();
  if ( m_SegmentedImage )
    {
    dummyTemplateImage->SetRegions( m_SegmentedImage->GetBufferedRegion() );  
    }
  else if ( m_ProbabilityImage )
    {
    dummyTemplateImage->SetRegions( m_ProbabilityImage->GetBufferedRegion() );  
    }
  else
    {
    dummyTemplateImage->SetRegions( m_Images[0]->GetBufferedRegion() );
    }

  if(0){
	std::cout<<"The means and precisions the optimizer has are:"<<std::endl;
        //std::cout << "m_Means.size(): " << m_Means.size() << std::endl;
        //std::cout << "m_Precisions.size(): " << m_Precisions.size() << std::endl;

	for(unsigned int i = 0; i<m_Means.size(); i++){
		vnl_vector<float> miini = m_Means[i];
		for(unsigned int k = 0; k<miini.size(); k++){
			std::cout<<"Means["<<i<<"]["<<k<<"]: "<<miini[k]<<std::endl;
		}
	}

	for(unsigned int i = 0; i<m_Precisions.size(); i++){
		vnl_matrix<float> pressat = m_Precisions[i];
		for(unsigned int k = 0; k<2; k++){
			for(unsigned int l = 0; l<2; l++){
				std::cout<<"Precisions["<<k<<"]["<<l<<"]["<<i<<"]: "<<pressat[k][l]<<std::endl;
			}
		}
	}


	
  }
  
  // Set up calculator
  AtlasMeshToIntensityImageGradientCalculator::Pointer  calculator = AtlasMeshToIntensityImageGradientCalculator::New();
  calculator->SetNumberOfThreads( this->GetNumberOfThreads() );
  calculator->SetLabelImage( dummyTemplateImage );
  if ( m_SegmentedImage )
    {
    calculator->SetMapCompToComp(m_mapCompToComp);
    calculator->SetSegmentedImage( m_SegmentedImage );    
    }
  else if ( m_ProbabilityImage )
    {
    calculator->SetProbabilityImage( m_ProbabilityImage );
    }
  else
    {
    calculator->SetImages( m_Images );
    calculator->SetMeans( m_Means );
    calculator->SetPrecisions( m_Precisions );
    }
  calculator->SetMeshToImageTransform( m_MeshToImageTransform );

  
  // Set up a mesh
  AtlasMesh::Pointer  mesh = AtlasMesh::New();
  mesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( position ) );
  mesh->SetCells( m_Mesh->GetCells() );
  mesh->SetPointData( m_Mesh->GetPointData() );
  mesh->SetCellData( m_Mesh->GetCellData() );
  
  
  if(0){
  	// std::cout<<"Problematic tetrahedron number 1652908 stats: "<<std::endl;
  	AtlasMesh::CellAutoPointer  cell;
  	mesh->GetCell( 1652908, cell );
  	AtlasMesh::CellType::PointIdIterator  pit = cell->PointIdsBegin();
  	AtlasMesh::PointType  p;
  	mesh->GetPoint( *pit, &p );
  	// std::cout<<"Point id: "<<*pit<<" : "<<p[0]<<" , "<<p[1]<<" , "<<p[2]<<std::endl;
  
  	++pit;
  	mesh->GetPoint( *pit, &p );
  	// std::cout<<"Point id: "<<*pit<<" : "<<p[0]<<" , "<<p[1]<<" , "<<p[2]<<std::endl;
  	++pit;
  	mesh->GetPoint( *pit, &p );
  	// std::cout<<"Point id: "<<*pit<<" : "<<p[0]<<" , "<<p[1]<<" , "<<p[2]<<std::endl;
  	++pit;
  	mesh->GetPoint( *pit, &p );
  	// std::cout<<"Point id: "<<*pit<<" : "<<p[0]<<" , "<<p[1]<<" , "<<p[2]<<std::endl;

  	ReferenceTetrahedronInfo  info;
  	mesh->GetCellData( 1652908, &info );

  	const double  referenceVolumeTimesK = info.m_ReferenceVolumeTimesK;

  	const double  z11 = info.m_Z11;
  	const double  z21 = info.m_Z21;
  	const double  z31 = info.m_Z31;
  	const double  z41 = info.m_Z41;

  	const double  z12 = info.m_Z12;
  	const double  z22 = info.m_Z22;
  	const double  z32 = info.m_Z32;
  	const double  z42 = info.m_Z42;

  	const double  z13 = info.m_Z13;
  	const double  z23 = info.m_Z23;
  	const double  z33 = info.m_Z33;
  	const double  z43 = info.m_Z43;
  	// std::cout << "Reference triangle info components: "<<std::endl;
	// std::cout << "z11: "<<z11<<std::endl;
	// std::cout << "z21: "<<z21<<std::endl;
	// std::cout << "z31: "<<z31<<std::endl;
	// std::cout << "z41: "<<z41<<std::endl;
	// std::cout << "z12: "<<z12<<std::endl;
	// std::cout << "z22: "<<z22<<std::endl;
	// std::cout << "z32: "<<z32<<std::endl;
	// std::cout << "z42: "<<z42<<std::endl;
	// std::cout << "z13: "<<z13<<std::endl;
	// std::cout << "z23: "<<z23<<std::endl;
	// std::cout << "z33: "<<z33<<std::endl;
	// std::cout << "z43: "<<z43<<std::endl;
  }
  
  // Rasterize mesh
  calculator->Rasterize( mesh );
  
  // Retrieve results
  cost = calculator->GetMinLogLikelihoodTimesPrior();
  gradient = calculator->GetPositionGradient();
  
}
  
  
  
  
//
//
//
void AtlasMeshDeformationConjugateGradientOptimizer
::Initialize()
{
  //
  m_Cost = 0;
  m_OldCost = 0;
  m_Position = 0;
  m_Gradient = 0;
  m_OldGradient = 0;
  m_Direction = 0;
  //m_StepSize = 0;
    
  //
  m_Initialized = true;
}

  
  
//
//
//
double
AtlasMeshDeformationConjugateGradientOptimizer
::PerformOneIteration()
{

  //
  const double  tolX = m_MaximalDeformationStopCriterion /* 1e-9 */;
  
  //
  if ( !m_Initialized )
    {
    this->Initialize();  
    }

  
  //
  if ( !m_Gradient )
    { 
    // Initialize if we need to
    if ( m_Verbose )
      {
      std::cout << "Obtaining first cost and gradient CG" << std::endl;
      }
    m_Position = m_Mesh->GetPoints();
    this->GetCostAndGradient( m_Position, m_Cost, m_Gradient );
    if ( m_Verbose )
      {
      std::cout << "       cost: " << m_Cost << std::endl;
      }
    }
   
   
  // Compute direction
  if ( !m_OldGradient )
    {  
    // Initially use steepest descent direction
    if ( m_Verbose )
      {
      std::cout << "Initialize to steepest descent direction" << std::endl;
      }
    m_Direction = AtlasPositionGradientContainerType::New();
    for ( AtlasPositionGradientContainerType::ConstIterator  it = m_Gradient->Begin();
          it != m_Gradient->End(); ++it )
      {
      m_Direction->InsertElement( it.Index(), -it.Value() );
      }
    }
  else
    {
    // Hestenes-Stiefel:  beta = (g'*(g-g_old))/((g-g_old)'*d); 
    //// std::cout << "Hestenes-Stiefel" << std::endl;
    AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin();
    AtlasPositionGradientContainerType::ConstIterator  oldGradIt = m_OldGradient->Begin();
    AtlasPositionGradientContainerType::ConstIterator  dirIt = m_Direction->Begin();
    double  betaNumerator = 0.0;
    double  betaDenominator = 1e-15;
    for ( ; gradIt != m_Gradient->End(); ++gradIt, ++oldGradIt, ++dirIt )
      {
      for ( int i = 0; i < 3; i++ )
        {
        double  tmp = gradIt.Value()[i] - oldGradIt.Value()[i];
        betaNumerator += gradIt.Value()[i] * tmp;
        betaDenominator += dirIt.Value()[i] * tmp;
        }
      }
    const double  beta = betaNumerator / betaDenominator;

    // d = -g + beta*d;
    gradIt = m_Gradient->Begin();
    for ( AtlasPositionGradientContainerType::Iterator  dirIt = m_Direction->Begin();
          dirIt != m_Direction->End(); ++dirIt, ++gradIt )
      {
      dirIt.Value() = -gradIt.Value() + beta * dirIt.Value();
      }
      
      
    // Restart if not a direction of sufficient descent (g'*d > -tolX)
    double  directionalDerivative = 0.0;
    gradIt = m_Gradient->Begin();
    for ( AtlasPositionGradientContainerType::ConstIterator  dirIt = m_Direction->Begin();
          dirIt != m_Direction->End(); ++dirIt, ++gradIt )
      {
      for ( int i = 0; i < 3; i++ )
        {
        directionalDerivative += dirIt.Value()[i] * gradIt.Value()[i];
        }
      }
    //// std::cout << "directionalDerivative: " << directionalDerivative << std::endl;
      
    if ( directionalDerivative > -tolX )
      {
      if ( m_Verbose )
        {
        std::cout << "Restarting Conjugate Gradient" << std::endl;
        }

      // d = -g;
      AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin();
      AtlasPositionGradientContainerType::Iterator  dirIt = m_Direction->Begin();
      for ( ; dirIt != m_Direction->End(); ++dirIt, ++gradIt )
        {
        dirIt.Value() = -gradIt.Value();  
        }
        
      }

    } // End test if this is the first direction
    
  
  //
  m_OldGradient = m_Gradient;
  
  // Compute directional derivative
  double  directionalDerivative = 0.0;
  AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin();
  AtlasPositionGradientContainerType::ConstIterator  dirIt = m_Direction->Begin();
  for ( ; dirIt != m_Direction->End(); ++dirIt, ++gradIt )
    {
    for ( int i = 0; i < 3; i++ )
      {
      directionalDerivative += dirIt.Value()[i] * gradIt.Value()[i];
      }
    }
  //// std::cout << "directionalDerivative: " << directionalDerivative << std::endl;
  
  // Check that progress can be made along direction
  if ( directionalDerivative > -tolX )
    {
    if ( m_Verbose )
      {
      std::cout << "Directional Derivative below TolX" << std::endl;
      }
    return 0.0;
    }
    
    // std::cout<<"Directional derivative is: "<<directionalDerivative<<std::endl;
  
  // Make an initial guess of the step size
  
  // Make sure it's not something completely insane
  double  maximalDirectionDeformation = 0.0f;
  for ( AtlasPositionGradientContainerType::ConstIterator  dirIt = m_Direction->Begin(); dirIt != m_Direction->End(); ++dirIt )
    {
    AtlasPositionGradientType  trialStep = dirIt.Value();

    if ( trialStep.GetNorm() > maximalDirectionDeformation )
      {
      maximalDirectionDeformation = trialStep.GetNorm();
      }
    }
    
     // std::cout<<"Maximal direction deformation is: "<<maximalDirectionDeformation<<std::endl; 

  double  stepSize = 0;
  if ( !m_OldCost )
    {
    // t = min(1,1/sum(abs(g)));
    double sumAbsGradient = 0.0;
    for ( AtlasPositionGradientContainerType::ConstIterator  gradIt = m_Gradient->Begin();
          gradIt != m_Gradient->End(); ++gradIt )
      {
      for ( int i = 0; i < 3; i++ )
        {
        sumAbsGradient += fabs( gradIt.Value()[i] );  
        }
      }
    stepSize = 1.0 / sumAbsGradient;

    if ( ( stepSize * maximalDirectionDeformation ) > 20 )
      {
      stepSize = 20.0 / maximalDirectionDeformation;
      if ( m_Verbose )
        {
        std::cout << "Clipping initial step size to 20 voxel maximal deformation" << std::endl;
        }
      }
    else if ( ( stepSize * maximalDirectionDeformation ) < 1 )
      {
      stepSize = 1.0 / maximalDirectionDeformation;
      if ( m_Verbose )
        {
        std::cout << "Boosting initial step size to 1 voxel maximal deformation" << std::endl;        
        }
      }
    else
      {
      if ( m_Verbose )
        {
        std::cout << "Accepting initial step size heuristics: " << stepSize * maximalDirectionDeformation << " voxel maximalDeformation" << std::endl;
        }
      }
  
    }
  else
    {
    // t = min(1,2*(f-f_old)/(gtd));
    stepSize = 2 * ( m_Cost - m_OldCost ) / directionalDerivative;

      
    if ( stepSize < 0 )
      {
      stepSize = 5.0 / maximalDirectionDeformation;  
      //// std::cout << "Step size heuristics gave smth negative; trying 5 voxel maximal deformation instead" << std::endl;        
      }
    else if ( ( stepSize * maximalDirectionDeformation ) > 20 )
      {
      stepSize = 20.0 / maximalDirectionDeformation;
      //// std::cout << "Clipping step size to 20 voxel maximal deformation" << std::endl;        
      }
    else
      {
      //// std::cout << "Accepting step size heuristics: " << stepSize * maximalDirectionDeformation << " voxel maximalDeformation" << std::endl;  
      }
    }
    
  
  

  
  // Do a line search
  m_OldCost = m_Cost;
  const double  referenceCost = m_Cost;
  // std::cout<<"StepSize: "<<stepSize<<std::endl;

  double maximalDeformation = this->BackTrack( m_Position, stepSize, m_Cost, m_Gradient, m_Direction, directionalDerivative, tolX, referenceCost );

  if ( m_Verbose )
    {
    std::cout << "maximalDeformation: " << maximalDeformation << std::endl;
    }

  

  // Make sure the current position is set to the mesh
  m_Mesh->SetPoints( m_Position );
  
  return maximalDeformation;
}



//
//
//
double
AtlasMeshDeformationConjugateGradientOptimizer
::BackTrack( AtlasMesh::PointsContainer::Pointer&  position, double& stepSize, double&  cost, AtlasPositionGradientContainerType::Pointer& gradient, 
             const AtlasPositionGradientContainerType* direction, double directionalDerivative, double tolX, double referenceCost )
{
  
  // Add the step to the current position to obtain the trial position
  double  maximalDeformation = 0.0f;
  AtlasMesh::PointsContainer::Pointer  trialPosition = AtlasMesh::PointsContainer::New();
  AtlasMesh::PointsContainer::ConstIterator  posIt = position->Begin();
  AtlasPositionGradientContainerType::ConstIterator  dirIt = direction->Begin();
  for ( ; posIt != position->End(); ++posIt, ++dirIt )
    {
    AtlasPositionGradientType  trialStep = stepSize * dirIt.Value();

    if ( trialStep.GetNorm() > maximalDeformation )
      {
      maximalDeformation = trialStep.GetNorm();
      }

    trialPosition->InsertElement( posIt.Index(), posIt.Value() + trialStep );
    }
  // std::cout << "Using trialPosition with maximalDeformation: " << maximalDeformation << std::endl;  

  // Try out this new position
  AtlasPositionGradientContainerType::Pointer  newGradient = 0;
  double  newCost = 0.0;
  this->GetCostAndGradient( trialPosition, newCost, newGradient );
  // std::cout<<"BackTrack: newCost: "<<newCost<<" oldCost: "<<referenceCost<<std::endl;
  // Step-halving until we're satisfied

    while ( newCost > ( referenceCost + 1e-4 * stepSize * directionalDerivative ))
    {
              
      // std::cout<<"The value added to the cost: "<<1e-4 * stepSize * directionalDerivative<<std::endl;
      // std::cout<<"Step-halving: newCost: "<<newCost<<" oldCost: "<<referenceCost<<std::endl;
    if ( m_Verbose )
      {
      //// std::cout << "Halving step size" << std::endl;
      // std::cout << "." << std::flush;
      }
    stepSize *= 0.5;

    // Add the step to the current position to obtain the trial position
    maximalDeformation = 0.0f;
    trialPosition = AtlasMesh::PointsContainer::New();
    posIt = position->Begin();
    dirIt = direction->Begin();
    for ( ; posIt != position->End(); ++posIt, ++dirIt )
      {
      AtlasPositionGradientType  trialStep = stepSize * dirIt.Value();

      if ( trialStep.GetNorm() > maximalDeformation )
        {
        maximalDeformation = trialStep.GetNorm();
        }

      trialPosition->InsertElement( posIt.Index(), posIt.Value() + trialStep );
      }
  
    // Try out this new position
    this->GetCostAndGradient( trialPosition, newCost, newGradient );

    // Check whether step size has become too small
    if ( maximalDeformation < tolX ) 
      {
      if ( m_Verbose )
        {
        // std::cout << std::endl;
        // std::cout << "Backtracking Line Search Failed" << std::endl;
        }
      stepSize = 0;
      maximalDeformation = 0;
      newCost = cost;
      newGradient = gradient;
      break;
      }

    } // End loop over step halving
  if ( m_Verbose )
    {
    // std::cout << std::endl;
    // std::cout<<"After step-halving the cost is: "<<newCost<<std::endl;
    }


    position = trialPosition;
    cost = newCost;
    gradient = newGradient;

   return maximalDeformation;
}




} // end namespace kvl

