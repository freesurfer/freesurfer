#ifndef ASYMMETRICTENSORREADERSTRATEGY_H_
#define ASYMMETRICTENSORREADERSTRATEGY_H_

#include "TensorReaderStrategy.h"

class AsymmetricTensorReaderStrategy : public TensorReaderStrategy
{
public:

  typedef TensorReaderStrategy::TensorPixelType TensorPixelType;
  typedef TensorReaderStrategy::TensorImageType TensorImageType;

  AsymmetricTensorReaderStrategy();
  
  ~AsymmetricTensorReaderStrategy();  
  
  virtual TensorImageType::Pointer GetTensors();  
  
  
private:

};

#endif /*ASYMMETRICTENSORREADERSTRATEGY_H_*/
