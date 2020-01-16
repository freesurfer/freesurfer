#include "kvlConditionalGaussianEntropyCostAndGradientCalculator.h"

#include <itkMath.h>
#include "kvlTetrahedronInteriorConstIterator.h"
// #include "itkTimeProbe.h"

#define SUBTRACT_MARGINAL_ENTROPY 0

namespace kvl
{

//
//
//
ConditionalGaussianEntropyCostAndGradientCalculator
::ConditionalGaussianEntropyCostAndGradientCalculator()
{

  m_Image = 0;
  
}


//
//
//
ConditionalGaussianEntropyCostAndGradientCalculator
::~ConditionalGaussianEntropyCostAndGradientCalculator()
{
}



//
//
//
void 
ConditionalGaussianEntropyCostAndGradientCalculator
::SetImage( const ImageType*  image )
{

  m_Image = image;
}


//
//
//
void
ConditionalGaussianEntropyCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{
  //
  //itk::TimeProbe clock;
  //clock.Start();

  // Clean up any previous mess
  m_Abort = false;
  // m_ThreadSpecificNs.clear();
  // m_ThreadSpecificLs.clear();
  // m_ThreadSpecificQs.clear();
  // m_ThreadSpecificNGradients.clear();
  // m_ThreadSpecificLGradients.clear();
  // m_ThreadSpecificQGradients.clear();
  // m_ThreadSpecificPriorCosts.clear();
  // m_ThreadSpecificPriorGradients.clear();
  
  // Make sure thread-specific, class-specific quantities we're going to collect are defined and
  // initialized to zero
  const int  numberOfClasses = mesh->GetPointData()->Begin().Value().m_Alphas.Size();
  AtlasPositionGradientType  zeroEntry( 0.0f );

  //std::cout << "Preparing..." << std::flush;
  const bool  memoryAlreadyAllocated = ( m_ThreadSpecificNs.size() > 0 );
  //std::cout << "memoryAlreadyAllocated: " << memoryAlreadyAllocated << std::endl;
    
  for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
    {
    if ( !m_OnlyDeformationPrior )
      {
      if ( !memoryAlreadyAllocated )
        {
        //std::cout << "Allocating memory" << std::endl;
        std::vector< double >  Ns;
        std::vector< double >  Ls;
        std::vector< double >  Qs;
          
        std::vector< AtlasPositionGradientContainerType::Pointer >  NGradients;
        std::vector< AtlasPositionGradientContainerType::Pointer >  LGradients;
        std::vector< AtlasPositionGradientContainerType::Pointer >  QGradients;
      
        // For each class, initialize N, L, Q, and their gradients
        for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
          {
          // Scalars
          Ns.push_back( 0.0 );  
          Ls.push_back( 0.0 );  
          Qs.push_back( 0.0 );  
            
          // Gradients
          AtlasPositionGradientContainerType::Pointer  NGradient = AtlasPositionGradientContainerType::New();
          AtlasPositionGradientContainerType::Pointer  LGradient = AtlasPositionGradientContainerType::New();
          AtlasPositionGradientContainerType::Pointer  QGradient = AtlasPositionGradientContainerType::New();
          for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
                pointIt != mesh->GetPoints()->End(); ++pointIt )
            {
            NGradient->InsertElement( pointIt.Index(), zeroEntry );
            LGradient->InsertElement( pointIt.Index(), zeroEntry );
            QGradient->InsertElement( pointIt.Index(), zeroEntry );
            }
          NGradients.push_back( NGradient );
          LGradients.push_back( LGradient );
          QGradients.push_back( QGradient );
          
          } // End loop over classes

        //
        m_ThreadSpecificNs.push_back( Ns );
        m_ThreadSpecificLs.push_back( Ls );
        m_ThreadSpecificQs.push_back( Qs );
        
        m_ThreadSpecificNGradients.push_back( NGradients );
        m_ThreadSpecificLGradients.push_back( LGradients );
        m_ThreadSpecificQGradients.push_back( QGradients );
        
        }
      else
        {
        //std::cout << "zeroing existing memory" << std::endl;
          
        // Memory already allocated; just zero out the values  
        m_ThreadSpecificNs[ threadNumber ].assign( numberOfClasses, 0.0 );
        m_ThreadSpecificLs[ threadNumber ].assign( numberOfClasses, 0.0 );
        m_ThreadSpecificQs[ threadNumber ].assign( numberOfClasses, 0.0 );

        for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
          {
          AtlasPositionGradientContainerType::Iterator  NGradientIt =
                ( m_ThreadSpecificNGradients[ threadNumber ] )[ classNumber ]->Begin();
          AtlasPositionGradientContainerType::Iterator  LGradientIt =
                ( m_ThreadSpecificLGradients[ threadNumber ] )[ classNumber ]->Begin();
          AtlasPositionGradientContainerType::Iterator  QGradientIt =
                ( m_ThreadSpecificQGradients[ threadNumber ] )[ classNumber ]->Begin();
          for ( ; 
                NGradientIt != ( m_ThreadSpecificNGradients[ threadNumber ] )[ classNumber ]->End(); 
                ++NGradientIt, ++LGradientIt, ++QGradientIt )
            {
            NGradientIt.Value().Fill( 0 );
            LGradientIt.Value().Fill( 0 );
            QGradientIt.Value().Fill( 0 );
            }
          } // End loop over classes  

        } // End test if memory already allocated  
          
      }
    
    if ( !m_IgnoreDeformationPrior )
      {
      if ( !memoryAlreadyAllocated )
        {
        //std::cout << "Allocating memory" << std::endl;
          
        //
        m_ThreadSpecificPriorCosts.push_back( 0.0 );
        AtlasPositionGradientContainerType::Pointer  priorGradient = AtlasPositionGradientContainerType::New();
        for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
              pointIt != mesh->GetPoints()->End(); ++pointIt )
          {
          priorGradient->InsertElement( pointIt.Index(), zeroEntry );
          }
        m_ThreadSpecificPriorGradients.push_back( priorGradient );
        }
      else
        {
        //std::cout << "Zeroing out memory" << std::endl;
          
        m_ThreadSpecificPriorCosts[ threadNumber ] = 0.0;
        for ( AtlasPositionGradientContainerType::Iterator  priorGradientIt 
                      = m_ThreadSpecificPriorGradients[ threadNumber ]->Begin();
              priorGradientIt != m_ThreadSpecificPriorGradients[ threadNumber ]->End();
              ++priorGradientIt )
          {
          priorGradientIt.Value().Fill( 0.0 );
          }
          
        }  
      
      }
      
    } // End loop over threads  
  
  
  //clock.Stop();
  //std::cout << "Time taken by initialization: " << clock.GetMean() << std::endl;
  
  // Now rasterize
  //clock.Reset();
  //clock.Start();
  AtlasMeshRasterizor::Rasterize( mesh );
  //clock.Stop();
  //std::cout << "Time taken by rasterization: " << clock.GetMean() << std::endl;

  //clock.Reset();
  //clock.Start();

  // Initialize gradient-to-return to zero
  m_PositionGradient = AtlasPositionGradientContainerType::New();
  for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
        pointIt != mesh->GetPoints()->End(); ++pointIt )
    {
    m_PositionGradient->InsertElement( pointIt.Index(), zeroEntry );
    }
  
  
  // Make sure everything has gone smoothly
  if ( m_Abort )
    {
    // Something has gone wrong
    m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
    return;
    }

  // Add contributions of prior
  double  priorCost = 0.0;
  if ( !m_IgnoreDeformationPrior )
    {
    for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
      {
      // Cost
      if ( std::isnan( m_ThreadSpecificPriorCosts[ threadNumber ] ) || 
          std::isinf( m_ThreadSpecificPriorCosts[ threadNumber ] ) )
        {
        // Something has gone wrong
        m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
        return;
        }
      priorCost += m_ThreadSpecificPriorCosts[ threadNumber ];
        
      // Gradient
      AtlasPositionGradientContainerType::ConstIterator  priorGradientIt =  
                                      m_ThreadSpecificPriorGradients[ threadNumber ]->Begin();
      AtlasPositionGradientContainerType::Iterator  gradientIt = m_PositionGradient->Begin();
      for ( ; gradientIt != m_PositionGradient->End(); 
            ++gradientIt, ++priorGradientIt )
        {
        gradientIt.Value() += priorGradientIt.Value();
        }

      } // End loop over threads
    }
  //std::cout << "priorCost: " << priorCost << std::endl;
    
  
  // Stitch together the results based on thread-specific, class-specific quantities
  double  dataCost = 0.0;

  if ( !m_OnlyDeformationPrior )
    {
    double  marginalN = 0.0;
#if SUBTRACT_MARGINAL_ENTROPY
    double  marginalL = 0.0;
    double  marginalQ = 0.0;
#endif  
    
    // Add the various contributions to the cost function
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      double  N = 1e-15;
      double  L = 0.0;
      double  Q = 0.0;
      for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
        {
        //
        N += ( m_ThreadSpecificNs[ threadNumber ] )[ classNumber ];
        L += ( m_ThreadSpecificLs[ threadNumber ] )[ classNumber ];
        Q += ( m_ThreadSpecificQs[ threadNumber ] )[ classNumber ];
        } // end loop over threads
    
      //  
      const double  variance = ( Q * N - pow( L, 2 ) ) / pow( N, 2 ) + 1e-15; 
      const double  entropy = log( variance ) + 1;
    
      // 
      dataCost += N * entropy;
      marginalN += N;
#if SUBTRACT_MARGINAL_ENTROPY
      marginalL += L;
      marginalQ += Q;
#endif      
      // //
      // std::cout << "variance[ " << classNumber << " ]: " << variance << std::endl;
      // std::cout << "N[ " << classNumber << " ]: " << N << std::endl;
      // std::cout << "entropy[ " << classNumber << " ]: " << entropy << std::endl;
      
      } // End loop over classes
    dataCost /= marginalN;
#if SUBTRACT_MARGINAL_ENTROPY
    const double  marginalVariance = ( marginalQ * marginalN - pow( marginalL, 2 ) ) / pow( marginalN, 2 ) + 1e-15; 
    const double  marginalEntropy = log( marginalVariance ) + 1;
    
    // std::cout << "marginalVariance: " << marginalVariance << std::endl;
    // std::cout << "marginalN: " << marginalN << std::endl;
    // std::cout << "marginalEntropy: " << marginalEntropy << std::endl;
#endif
    
    //
    //std::cout << "marginalN: " << marginalN << std::endl;
    //std::cout << "dataCost: " << dataCost << std::endl;
    
    
    // Add the various contributions to the gradient
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      double  N = 1e-15;
      double  L = 0.0;
      double  Q = 0.0;
      for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
        {
        //
        N += ( m_ThreadSpecificNs[ threadNumber ] )[ classNumber ];
        L += ( m_ThreadSpecificLs[ threadNumber ] )[ classNumber ];
        Q += ( m_ThreadSpecificQs[ threadNumber ] )[ classNumber ];
        } // end loop over threads

      //  
      const double  variance = ( Q * N - pow( L, 2 ) ) / pow( N, 2 ) + 1e-15; 
      const double  entropy = log( variance ) + 1;

#if SUBTRACT_MARGINAL_ENTROPY 
      const double  Nweight = ( entropy - dataCost + 
                                Q / ( variance * N ) - 
                                marginalQ / ( marginalVariance * marginalN ) ) / marginalN;
      const double  Lweight = ( -2 * L / ( variance * N ) + 
                                 2 * marginalL / ( marginalVariance * marginalN ) ) / marginalN;
      const double  Qweight = ( 1/variance - 1/marginalVariance ) / marginalN;
#else
      const double  Nweight = ( entropy - dataCost + Q / ( variance * N ) - 2 ) / marginalN;
      const double  Lweight = ( -2 * L / ( variance * N ) ) / marginalN;
      const double  Qweight = ( 1/variance ) / marginalN;
#endif      
      
      //std::cout << "N[ " << classNumber << " ]: " << N << std::endl;
      //std::cout << "L[ " << classNumber << " ]: " << L << std::endl;
      //std::cout << "Q[ " << classNumber << " ]: " << Q << std::endl;
      //std::cout << "Nweight[ " << classNumber << " ]: " << Nweight << std::endl;
      //std::cout << "Lweight[ " << classNumber << " ]: " << Lweight << std::endl;
      //std::cout << "Qweight[ " << classNumber << " ]: " << Qweight << std::endl;
      //std::cout << std::endl;
      
      for ( int threadNumber = 0; threadNumber < this->GetNumberOfThreads(); threadNumber++ )
        {
        AtlasPositionGradientContainerType::ConstIterator  NGradientIt = 
                ( m_ThreadSpecificNGradients[ threadNumber ] )[ classNumber ]->Begin();
        AtlasPositionGradientContainerType::ConstIterator  LGradientIt = 
                ( m_ThreadSpecificLGradients[ threadNumber ] )[ classNumber ]->Begin();
        AtlasPositionGradientContainerType::ConstIterator  QGradientIt =
                ( m_ThreadSpecificQGradients[ threadNumber ] )[ classNumber ]->Begin();
        for ( AtlasPositionGradientContainerType::Iterator  gradientIt = m_PositionGradient->Begin(); 
              gradientIt != m_PositionGradient->End(); 
              ++gradientIt, ++NGradientIt, ++LGradientIt, ++QGradientIt )
          {
          gradientIt.Value() += Nweight * NGradientIt.Value() + 
                                Lweight * LGradientIt.Value() + 
                                Qweight * QGradientIt.Value();
          }

        } // End loop over threads
        
      } // End loop over classes
      
    //std::cout << "- and done" << std::endl;
#if SUBTRACT_MARGINAL_ENTROPY
    // std::cout << "dataCost - marginalEntropy: " << dataCost << " - " << marginalEntropy << std::endl;

    dataCost -= marginalEntropy;  
#endif    
    }
    
  //std::cout << "dataCost: " << dataCost << std::endl;
  
  //clock.Stop();
  //std::cout << "Time taken by results gathering: " << clock.GetMean() << std::endl;

  //
  m_MinLogLikelihoodTimesPrior = dataCost + priorCost;

  
  // Take care of the desired boundary conditions
  switch( m_BoundaryCondition ) 
    {
    case SLIDING: 
      {
      //std::cout << "SLIDING" << std::endl;
      this->ImposeSlidingBoundaryConditions( mesh );      
      break;
      } 
    case AFFINE: 
      {
      //std::cout << "AFFINE" << std::endl;
      this->ImposeAffineBoundaryConditions( mesh );  
      break;
      } 
    case TRANSLATION: 
      {
      //std::cout << "TRANSLATION" << std::endl;
      this->ImposeTranslationBoundaryConditions( mesh );  
      break;
      } 
    default:
      {
      //std::cout << "NONE" << std::endl;
      break;
      }
    }
  
}    
  
    
  
  
//
//
//
bool
ConditionalGaussianEntropyCostAndGradientCalculator
::RasterizeTetrahedron( const AtlasMesh* mesh, 
                        AtlasMesh::CellIdentifier tetrahedronId,
                        int threadNumber )
{

  //
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

  
  if ( !m_IgnoreDeformationPrior )
    {
    if ( !this->AddPriorContributionOfTetrahedron( p0, p1, p2, p3, info,
                                m_ThreadSpecificPriorCosts[ threadNumber ],
                                m_ThreadSpecificPriorGradients[ threadNumber ]->ElementAt( id0 ),
                                m_ThreadSpecificPriorGradients[ threadNumber ]->ElementAt( id1 ),
                                m_ThreadSpecificPriorGradients[ threadNumber ]->ElementAt( id2 ),
                                m_ThreadSpecificPriorGradients[ threadNumber ]->ElementAt( id3 ) ) )
      {
      m_Abort = true;
      return false;
      }
 
    } // End test if we need to include prior term

    
  if( !m_OnlyDeformationPrior )
    {
    //
    const AtlasAlphasType&  alphasInVertex0 = mesh->GetPointData()->ElementAt( id0 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex1 = mesh->GetPointData()->ElementAt( id1 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex2 = mesh->GetPointData()->ElementAt( id2 ).m_Alphas;
    const AtlasAlphasType&  alphasInVertex3 = mesh->GetPointData()->ElementAt( id3 ).m_Alphas;
      
    if ( ( alphasInVertex0.sum() + 
          alphasInVertex1.sum() + 
          alphasInVertex2.sum() + 
          alphasInVertex3.sum() ) < 1e-5 )
      {
      // All background -- can just as well skip this tetrahedron
      return true;
      }

    //
    std::vector< double >&  Ns = m_ThreadSpecificNs[ threadNumber ];
    std::vector< double >&  Ls = m_ThreadSpecificLs[ threadNumber ];
    std::vector< double >&  Qs = m_ThreadSpecificQs[ threadNumber ];
    std::vector< AtlasPositionGradientContainerType::Pointer >&  NGradients 
              = m_ThreadSpecificNGradients[ threadNumber ];
    std::vector< AtlasPositionGradientContainerType::Pointer >&  LGradients 
              = m_ThreadSpecificLGradients[ threadNumber ];
    std::vector< AtlasPositionGradientContainerType::Pointer >&  QGradients 
              = m_ThreadSpecificQGradients[ threadNumber ];
      
    // Loop over all voxels within the tetrahedron and do The Right Thing  
    const int  numberOfClasses = alphasInVertex0.Size();
    TetrahedronInteriorConstIterator< ImageType::PixelType >  it( m_Image, p0, p1, p2, p3 );
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      it.AddExtraLoading( alphasInVertex0[ classNumber ], 
                          alphasInVertex1[ classNumber ], 
                          alphasInVertex2[ classNumber ], 
                          alphasInVertex3[ classNumber ] );
      }  
      
    for ( ; !it.IsAtEnd(); ++it )
      {
      // Skip voxels with zero intensity
      if ( it.Value() == 0 )
        {
        //std::cout << "Skipping: " << it.Value() << std::endl;
        continue;
        }
        
        
      // Retrieve intensity     
      const double y = it.Value();
      const double ySquared = pow( y, 2 );
    
      for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
        {

        // Add contribution to N, L, and Q
        double&  N = Ns[ classNumber ];
        double&  L = Ls[ classNumber ];
        double&  Q = Qs[ classNumber ];
        N += it.GetExtraLoadingInterpolatedValue( classNumber );
        L += y * it.GetExtraLoadingInterpolatedValue( classNumber );
        Q += ySquared * it.GetExtraLoadingInterpolatedValue( classNumber );
        
        
        // Add contribution to NGradient
        AtlasPositionGradientType&  NGradientInVertex0 
                    = NGradients[ classNumber ]->ElementAt( id0 );
        AtlasPositionGradientType&  NGradientInVertex1
                    = NGradients[ classNumber ]->ElementAt( id1 );
        AtlasPositionGradientType&  NGradientInVertex2
                    = NGradients[ classNumber ]->ElementAt( id2 );
        AtlasPositionGradientType&  NGradientInVertex3
                    = NGradients[ classNumber ]->ElementAt( id3 );
        
        const double  xGradientBasisN = -it.GetExtraLoadingNextRowAddition( classNumber );
        const double  yGradientBasisN = -it.GetExtraLoadingNextColumnAddition( classNumber );
        const double  zGradientBasisN = -it.GetExtraLoadingNextSliceAddition( classNumber );

        // Add contribution to gradient in vertex 0
        NGradientInVertex0[ 0 ] += xGradientBasisN * it.GetPi0();
        NGradientInVertex0[ 1 ] += yGradientBasisN * it.GetPi0();
        NGradientInVertex0[ 2 ] += zGradientBasisN * it.GetPi0();

        // Add contribution to gradient in vertex 1
        NGradientInVertex1[ 0 ] += xGradientBasisN * it.GetPi1();
        NGradientInVertex1[ 1 ] += yGradientBasisN * it.GetPi1();
        NGradientInVertex1[ 2 ] += zGradientBasisN * it.GetPi1();
        
        // Add contribution to gradient in vertex 2
        NGradientInVertex2[ 0 ] += xGradientBasisN * it.GetPi2();
        NGradientInVertex2[ 1 ] += yGradientBasisN * it.GetPi2();
        NGradientInVertex2[ 2 ] += zGradientBasisN * it.GetPi2();
        
        // Add contribution to gradient in vertex 3
        NGradientInVertex3[ 0 ] += xGradientBasisN * it.GetPi3();
        NGradientInVertex3[ 1 ] += yGradientBasisN * it.GetPi3();
        NGradientInVertex3[ 2 ] += zGradientBasisN * it.GetPi3();
      
      
        // Add contribution to LGradient
        AtlasPositionGradientType&  LGradientInVertex0 
                    = LGradients[ classNumber ]->ElementAt( id0 );
        AtlasPositionGradientType&  LGradientInVertex1
                    = LGradients[ classNumber ]->ElementAt( id1 );
        AtlasPositionGradientType&  LGradientInVertex2
                    = LGradients[ classNumber ]->ElementAt( id2 );
        AtlasPositionGradientType&  LGradientInVertex3
                    = LGradients[ classNumber ]->ElementAt( id3 );
        
        const double  xGradientBasisL = y * xGradientBasisN;
        const double  yGradientBasisL = y * yGradientBasisN;
        const double  zGradientBasisL = y * zGradientBasisN;

        // Add contribution to gradient in vertex 0
        LGradientInVertex0[ 0 ] += xGradientBasisL * it.GetPi0();
        LGradientInVertex0[ 1 ] += yGradientBasisL * it.GetPi0();
        LGradientInVertex0[ 2 ] += zGradientBasisL * it.GetPi0();

        // Add contribution to gradient in vertex 1
        LGradientInVertex1[ 0 ] += xGradientBasisL * it.GetPi1();
        LGradientInVertex1[ 1 ] += yGradientBasisL * it.GetPi1();
        LGradientInVertex1[ 2 ] += zGradientBasisL * it.GetPi1();
        
        // Add contribution to gradient in vertex 2
        LGradientInVertex2[ 0 ] += xGradientBasisL * it.GetPi2();
        LGradientInVertex2[ 1 ] += yGradientBasisL * it.GetPi2();
        LGradientInVertex2[ 2 ] += zGradientBasisL * it.GetPi2();
        
        // Add contribution to gradient in vertex 3
        LGradientInVertex3[ 0 ] += xGradientBasisL * it.GetPi3();
        LGradientInVertex3[ 1 ] += yGradientBasisL * it.GetPi3();
        LGradientInVertex3[ 2 ] += zGradientBasisL * it.GetPi3();
      
        
        // Add contribution to QGradient
        AtlasPositionGradientType&  QGradientInVertex0 
                    = QGradients[ classNumber ]->ElementAt( id0 );
        AtlasPositionGradientType&  QGradientInVertex1
                    = QGradients[ classNumber ]->ElementAt( id1 );
        AtlasPositionGradientType&  QGradientInVertex2
                    = QGradients[ classNumber ]->ElementAt( id2 );
        AtlasPositionGradientType&  QGradientInVertex3
                    = QGradients[ classNumber ]->ElementAt( id3 );
        
        const double  xGradientBasisQ = ySquared * xGradientBasisN;
        const double  yGradientBasisQ = ySquared * yGradientBasisN;
        const double  zGradientBasisQ = ySquared * zGradientBasisN;

        // Add contribution to gradient in vertex 0
        QGradientInVertex0[ 0 ] += xGradientBasisQ * it.GetPi0();
        QGradientInVertex0[ 1 ] += yGradientBasisQ * it.GetPi0();
        QGradientInVertex0[ 2 ] += zGradientBasisQ * it.GetPi0();

        // Add contribution to gradient in vertex 1
        QGradientInVertex1[ 0 ] += xGradientBasisQ * it.GetPi1();
        QGradientInVertex1[ 1 ] += yGradientBasisQ * it.GetPi1();
        QGradientInVertex1[ 2 ] += zGradientBasisQ * it.GetPi1();
        
        // Add contribution to gradient in vertex 2
        QGradientInVertex2[ 0 ] += xGradientBasisQ * it.GetPi2();
        QGradientInVertex2[ 1 ] += yGradientBasisQ * it.GetPi2();
        QGradientInVertex2[ 2 ] += zGradientBasisQ * it.GetPi2();
        
        // Add contribution to gradient in vertex 3
        QGradientInVertex3[ 0 ] += xGradientBasisQ * it.GetPi3();
        QGradientInVertex3[ 1 ] += yGradientBasisQ * it.GetPi3();
        QGradientInVertex3[ 2 ] += zGradientBasisQ * it.GetPi3();
        
        } // End loop over classes
        
      } // End loop over all voxels within the tetrahedron

    } // End test if ignore data term  
      
 // Done
 return true;
}
  


} // end namespace kvl
