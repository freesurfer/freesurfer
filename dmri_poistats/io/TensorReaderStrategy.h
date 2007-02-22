#ifndef TENSORREADERSTRATEGY_H_
#define TENSORREADERSTRATEGY_H_

#include <string>

#include <itkDiffusionTensor3D.h>
#include <itkImage.h>

#include "dmri_poistats/datamodel/events/CommandUpdate.h"

class TensorReaderStrategy
{
public:

  typedef itk::DiffusionTensor3D< float > TensorPixelType;
  typedef itk::Image< TensorPixelType, 3 > TensorImageType;

  TensorReaderStrategy();

  virtual ~TensorReaderStrategy();
  
  void SetObserver( CommandUpdate* observer );
  
  virtual TensorImageType::Pointer GetTensors() = 0;  
  
  void SetFileName( std::string fileName );

protected:

  CommandUpdate* m_Observer;
  
  std::string m_TensorFileName;
  
};

#endif /*TENSORREADERSTRATEGY_H_*/
