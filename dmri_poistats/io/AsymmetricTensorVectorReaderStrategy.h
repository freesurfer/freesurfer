#ifndef ASYMMETRICTENSORVECTORREADER_H_
#define ASYMMETRICTENSORVECTORREADER_H_

#include "TensorReaderStrategy.h"

class AsymmetricTensorVectorReaderStrategy : public TensorReaderStrategy
{
public:

  typedef TensorReaderStrategy::TensorPixelType TensorPixelType;
  typedef TensorReaderStrategy::TensorImageType TensorImageType;

  AsymmetricTensorVectorReaderStrategy();
  
  ~AsymmetricTensorVectorReaderStrategy();  
  
  virtual TensorImageType::Pointer GetTensors();  
  
  
private:

};


#endif /*ASYMMETRICTENSORVECTORREADER_H_*/
