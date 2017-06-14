#include "kvlAtlasMeshSmoother.h"
#include "kvlAtlasMeshProbabilityImageStatisticsCollector.h"
#include "kvlAtlasMeshMultiAlphaDrawer.h"
#include "itkImageRegionIterator.h"
#include "itkDiscreteGaussianImageFilter.h"



namespace kvl
{


//
//
//
AtlasMeshSmoother
::AtlasMeshSmoother()
{
  m_MeshCollection = 0;
  m_Sigma0 = 1.0f;
  m_Sigma1 = 1.0f;
  m_Sigma2 = 1.0f;
}



//
//
//
AtlasMeshSmoother
::~AtlasMeshSmoother()
{
}




//
//
//
void 
AtlasMeshSmoother
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}




//
//
//
AtlasMeshCollection::Pointer
AtlasMeshSmoother
::GetSmoothedMeshCollection()
{

  // Sanity check on input
  if ( !m_MeshCollection )
    {
    itkExceptionMacro( << "No mesh collection set!" );
    }

  //  std::cout << "Smoothing mesh collection with sigma = " << m_Sigma << "..." << std::flush;

  // Construct container to hold smoothed alphas
  AtlasMesh::PointDataContainer::Pointer  smoothedParameters = 0;

  
  if ( ( m_Sigma0 == 0 ) && ( m_Sigma1 == 0 ) && ( m_Sigma2 ==0 ) )
    {
    // Simply copy, which is not only faster but also avoids numerical problems
    smoothedParameters = AtlasMesh::PointDataContainer::New();
    for ( AtlasMesh::PointDataContainer::ConstIterator  it = m_MeshCollection->GetPointParameters()->Begin();
          it != m_MeshCollection->GetPointParameters()->End(); ++it )
      {
      smoothedParameters->InsertElement( it.Index(), it.Value() );  
      }  
    }
  else
    {
    // 
    // Part I: Rasterize the mesh
    //
    
    // First determine size of template image for the alpha rasterizor
    typedef AtlasMeshMultiAlphaDrawer::ImageType  MultiAlphasImageType;
    typedef MultiAlphasImageType::SizeType  SizeType;
    SizeType  size;
    size[ 0 ] = 0;
    size[ 1 ] = 0;
    size[ 2 ] = 0;
    for ( AtlasMesh::PointsContainer::ConstIterator  it = m_MeshCollection->GetReferencePosition()->Begin();
          it != m_MeshCollection->GetReferencePosition()->End(); ++it )
      {
      for ( int i = 0; i < 3; i++ )
        {
        if ( it.Value()[ i ] > size[ i ] )
          {
          size[ i ] = static_cast< MultiAlphasImageType::SizeType::SizeValueType >( it.Value()[ i ] );
          }
        }
      }
    size[ 0 ]++;
    size[ 1 ]++;
    size[ 2 ]++;
      

    // Rasterize the mesh
    kvl::AtlasMeshMultiAlphaDrawer::Pointer  drawer = kvl::AtlasMeshMultiAlphaDrawer::New();
    drawer->SetRegions( size );
    drawer->Rasterize( m_MeshCollection->GetReferenceMesh() );
    MultiAlphasImageType::ConstPointer  alphasImage = drawer->GetImage();
    
    

    // 
    // Part II: Spatially smooth the image of probability vectors
    // 
    MultiAlphasImageType::Pointer  smoothedAlphasImage = const_cast< MultiAlphasImageType* >( alphasImage.GetPointer() );

    // Create an image holding one component of the probability vector image
    typedef itk::Image< float, 3 >  ComponentImageType;
    ComponentImageType::Pointer  componentImage = ComponentImageType::New();
    componentImage->SetRegions( size );
    componentImage->Allocate();

    // Smooth by smoothing each component individually
    //const int  numberOfClasses = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
    //for(int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
    //Loop only through the classes we explicitly want to smooth - unless the vector is empty; in that case, we do all of them (Eugenio added this)
    std::vector< int >  classesToSmooth;
    if ( m_ClassesToSmooth.size() > 0 )
      {
      classesToSmooth = m_ClassesToSmooth;
      }
    else
      {
      for ( int classNumber = 0; 
            classNumber <m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size(); 
            classNumber++ )
        {
        classesToSmooth.push_back( classNumber );
        }
      }
    for ( std::vector<int>::iterator classIt = classesToSmooth.begin(); 
          classIt != classesToSmooth.end(); 
          ++classIt )
      {

      std::cout<<"Smoothing class: "<<*(classIt)<<std::endl;
      // Extract the probability image
      typedef itk::Image< float, 3 >  ComponentImageType;
      ComponentImageType::Pointer  componentImage = ComponentImageType::New();
      componentImage->SetRegions( size );
      componentImage->Allocate();

      itk::ImageRegionIterator< MultiAlphasImageType >   vectorIt( smoothedAlphasImage, 
                                                                   smoothedAlphasImage->GetBufferedRegion() );
      itk::ImageRegionIterator< ComponentImageType >   componentIt( componentImage, 
                                                                    componentImage->GetBufferedRegion() );
      for ( ; !vectorIt.IsAtEnd(); ++vectorIt, ++componentIt )
        {
        componentIt.Value() = vectorIt.Value()[*(classIt)];//[ classNumber ];
        }
        
      // Smooth it
      typedef itk::DiscreteGaussianImageFilter< ComponentImageType, ComponentImageType >  SmootherType;
      SmootherType::Pointer smoother = SmootherType::New();
      smoother->SetInput( componentImage );
      smoother->SetMaximumError( 0.1 );
      smoother->SetUseImageSpacingOff();
      float  variances[ 3 ];
      variances[ 0 ] = m_Sigma0 * m_Sigma0;
      variances[ 1 ] = m_Sigma1 * m_Sigma1;
      variances[ 2 ] = m_Sigma2 * m_Sigma2;
      smoother->SetVariance( variances );
      smoother->Update();
      componentImage = smoother->GetOutput();
      
      // Put it back into the vector image
      vectorIt = itk::ImageRegionIterator< MultiAlphasImageType >( smoothedAlphasImage, 
                                                                   smoothedAlphasImage->GetBufferedRegion() );
      componentIt = itk::ImageRegionIterator< ComponentImageType >( componentImage, 
                                                                    componentImage->GetBufferedRegion() );
      for ( ; !vectorIt.IsAtEnd(); ++vectorIt, ++componentIt )
        {
        vectorIt.Value()[*(classIt)] = componentIt.Value();//[ classNumber ] = componentIt.Value();
        }
      
        
      } // End loop over classes
    
    
    //
    // Part III: Do EM estimation of the mesh node alphas that best fit the smoothed probability vectors 
    //
    
    // Get unnormalized counts
    double  cost = 0.0;
    m_MeshCollection->FlattenAlphas(); // Initialization to flat alphas
    for ( int iterationNumber = 0; iterationNumber < 10; iterationNumber++ )
      {
      AtlasMeshProbabilityImageStatisticsCollector::Pointer  statisticsCollector = 
                                              AtlasMeshProbabilityImageStatisticsCollector::New();
      statisticsCollector->SetProbabilityImage( smoothedAlphasImage );
      statisticsCollector->Rasterize( m_MeshCollection->GetReferenceMesh() );
      cost = statisticsCollector->GetMinLogLikelihood();
      //std::cout << "   EM iteration " << iterationNumber << " -> " << cost << std::endl;
      
      
      // Normalize and assign to the meshCollection's alpha vectors
      smoothedParameters = m_MeshCollection->GetPointParameters();
      AtlasMesh::PointDataContainer::Iterator  pointParamIt = smoothedParameters->Begin();
      AtlasMeshProbabilityImageStatisticsCollector::StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
      for ( ; pointParamIt != smoothedParameters->End(); ++pointParamIt, ++statIt )
        {
        if ( pointParamIt.Value().m_CanChangeAlphas )
          {
          pointParamIt.Value().m_Alphas = statIt.Value();
          pointParamIt.Value().m_Alphas /= ( statIt.Value().sum() + 1e-12 );
          }
          
        }
    
      } // End loop over EM iteration numbers
    
    } // End test if sigma = 0
    

  // Construct a mesh collection to return
  AtlasMeshCollection::Pointer  result = AtlasMeshCollection::New();
  result->SetPointParameters( smoothedParameters );
  result->SetCells( m_MeshCollection->GetCells() );
  result->SetReferencePosition(  m_MeshCollection->GetReferencePosition() );
  result->SetK( m_MeshCollection->GetK() );
  result->SetPositions( m_MeshCollection->GetPositions() );

  //  std::cout << "done!" << std::endl;
  //  m_MeshCollection->Write( "/tmp/original.txt" );  Eugenio: actually, original.txt will be smoothed too...
  //  result->Write( "/tmp/smoothed.txt" );

  return result;
}



} // end namespace kvl
