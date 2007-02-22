#ifndef SYMMETRICTENSORREADERSTRATEGY_H_
#define SYMMETRICTENSORREADERSTRATEGY_H_

#include "TensorReaderStrategy.h"

class SymmetricTensorReaderStrategy : public TensorReaderStrategy
{
public:

  typedef TensorReaderStrategy::TensorPixelType TensorPixelType;
  typedef TensorReaderStrategy::TensorImageType TensorImageType;

  SymmetricTensorReaderStrategy();
  
  ~SymmetricTensorReaderStrategy();  
  
  virtual TensorImageType::Pointer GetTensors();  
  
  
private:

};

#endif /*SYMMETRICTENSORREADERSTRATEGY_H_*/
