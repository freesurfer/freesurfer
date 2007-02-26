#include <itkImageSeriesReader.h>

#include "AsymmetricTensorVectorReaderStrategy.h"

AsymmetricTensorVectorReaderStrategy::AsymmetricTensorVectorReaderStrategy(){
  
}
  
AsymmetricTensorVectorReaderStrategy::~AsymmetricTensorVectorReaderStrategy(){
  
}

AsymmetricTensorVectorReaderStrategy::TensorImageType::Pointer
AsymmetricTensorVectorReaderStrategy::GetTensors(){

  typedef itk::VectorImage< float, 4 > FullTensorImageType;
  typedef itk::ImageFileReader< FullTensorImageType > FullTensorReaderType;
  FullTensorReaderType::Pointer tensorReader = FullTensorReaderType::New();
      
  tensorReader->SetFileName( m_TensorFileName );
  if( m_Observer != NULL ) {
    m_Observer->PostMessage( "reading tensors...\n" );
  }
  
  try { 
    tensorReader->Update();
  } catch( itk::ExceptionObject & excp ) {
    std::ostringstream output;
    output << "Error reading the series." << std::endl << excp << std::endl;
    if( m_Observer != NULL ) {
      m_Observer->PostErrorMessage( output.str() );
    } else {
      std::cerr << output.str();
    }
  }
  
  FullTensorImageType::Pointer fullTensors = tensorReader->GetOutput();
        
  // convert the full 9 component tensors to 6 component tensors
  //  - create and allocate a new DiffusionTensor3D image
  FullTensorImageType::RegionType dtiRegion = 
    fullTensors->GetLargestPossibleRegion();
  TensorImageType::SizeType size;
  double origin[ TensorImageType::RegionType::GetImageDimension() ];
  TensorImageType::IndexType start;
  TensorImageType::SpacingType spacing;
  
  for( unsigned int cDim=0; cDim<TensorImageType::RegionType::GetImageDimension(); cDim++ ) {    
    size[ cDim ] = dtiRegion.GetSize()[ cDim ];
    origin[ cDim ] = fullTensors->GetOrigin()[ cDim ];
    start[ cDim ] = 0;
    spacing[ cDim ] = fullTensors->GetSpacing()[ cDim ];
  }
      
  TensorImageType::RegionType region;
  region.SetSize( size );    
  region.SetIndex( start );
  
  TensorImageType::Pointer tensors = TensorImageType::New();
  tensors->SetRegions( region );
  tensors->SetOrigin( origin );
  tensors->SetSpacing( spacing );      
//  tensors->SetDirection( direction );  
  tensors->Allocate();  
  tensors->FillBuffer( 0.0 );
    
  //  - copy the data into the right areas
  const int nTensorRows = 3;
  const int nTensorCols = 3;

  if( m_Observer != NULL ) {
    m_Observer->PostMessage( "converting tensors to symmetric format\n" );
  }
  for( unsigned int cImageRow=0; cImageRow<size[ 0 ]; cImageRow++ ) {
    for( unsigned int cImageCol=0; cImageCol<size[ 1 ]; cImageCol++ ) {
      for( unsigned int cImageSlice=0; cImageSlice<size[ 2 ]; cImageSlice++ ) {

        TensorImageType::IndexType symmetricIndex;
        symmetricIndex[ 0 ] = cImageRow;
        symmetricIndex[ 1 ] = cImageCol;
        symmetricIndex[ 2 ] = cImageSlice;

        TensorPixelType symmetricTensor;

        FullTensorImageType::IndexType fullTensorIndex;
        fullTensorIndex[ 0 ] = cImageRow;
        fullTensorIndex[ 1 ] = cImageCol;
        fullTensorIndex[ 2 ] = cImageSlice;
        fullTensorIndex[ 3 ] = 0;
        
        FullTensorImageType::PixelType pixel = 
          fullTensors->GetPixel( fullTensorIndex );
      
        int cContinuousIndex = 0;
        
        // iterate through each component
        // TODO: actually we'll only need to iterate through some of them
        for( int cTensorRow = 0; cTensorRow<nTensorRows; cTensorRow++ ) {
          
          for( int cTensorCol = 0; cTensorCol<nTensorCols; cTensorCol++ ) {            
            symmetricTensor( cTensorRow, cTensorCol ) = pixel[ cContinuousIndex ];
            cContinuousIndex++;
          }
          
        }
        
        tensors->SetPixel( symmetricIndex, symmetricTensor );

      }
    }
  }
  
  // for some reason, I can't do this above, so I'm setting the cosine direction
  // here instead; it will seg fault otherwise.  
  TensorImageType::DirectionType dir;
  const FullTensorImageType::DirectionType& fullDirection = fullTensors->GetDirection();
  for( unsigned int cDim=0; cDim<TensorImageType::RegionType::GetImageDimension(); cDim++ ) {    
    for( int cDirection=0; cDirection<3; cDirection++ ) {
      dir[ cDim][ cDirection ] = fullDirection[ cDim ][ cDirection ];
    }
  }
  tensors->SetDirection( dir );
  
  return tensors;
}
