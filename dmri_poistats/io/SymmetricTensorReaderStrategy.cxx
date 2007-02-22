#include <itkImageSeriesReader.h>

#include "SymmetricTensorReaderStrategy.h"

SymmetricTensorReaderStrategy::SymmetricTensorReaderStrategy(){
  
}
  
SymmetricTensorReaderStrategy::~SymmetricTensorReaderStrategy(){
  
}

SymmetricTensorReaderStrategy::TensorImageType::Pointer
SymmetricTensorReaderStrategy::GetTensors(){
  // data stored in 6 components
  typedef itk::ImageFileReader< TensorImageType > TensorReaderType;
  TensorReaderType::Pointer tensorReader = TensorReaderType::New();
  tensorReader->SetFileName( m_TensorFileName );

  if( m_Observer != NULL ) {
    m_Observer->PostMessage( "reading symmetric tensors...\n" );
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
  TensorImageType::Pointer tensors = tensorReader->GetOutput();
  
  return tensors;
}
