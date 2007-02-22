#include "TensorReaderStrategy.h"

TensorReaderStrategy::TensorReaderStrategy(){
  this->m_Observer = NULL;
}


TensorReaderStrategy::~TensorReaderStrategy(){
  
}
  
void 
TensorReaderStrategy::SetObserver( CommandUpdate* observer ){
  this->m_Observer = observer;
}

void
TensorReaderStrategy::SetFileName( std::string fileName ) {
  m_TensorFileName = fileName;
}
