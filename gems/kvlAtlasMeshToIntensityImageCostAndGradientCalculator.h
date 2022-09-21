#ifndef __kvlAtlasMeshToIntensityImageCostAndGradientCalculator_h
#define __kvlAtlasMeshToIntensityImageCostAndGradientCalculator_h

#include "kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase.h"
#include "itkImage.h"
#include "kvlGMMLikelihoodImageFilter.h"


namespace kvl
{


class AtlasMeshToIntensityImageCostAndGradientCalculator:
        public AtlasMeshToIntensityImageCostAndGradientCalculatorBase
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshToIntensityImageCostAndGradientCalculator  Self;
  typedef AtlasMeshToIntensityImageCostAndGradientCalculatorBase Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshToIntensityImageCostAndGradientCalculator, kvlAtlasMeshToIntensityImageCostAndGradientCalculatorBase );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */
  void SetParameters( const std::vector< vnl_vector< double > >& means, 
                      const std::vector< vnl_matrix< double > >& variances,
                      const std::vector< double >&  mixtureWeights,
                      const std::vector< int >&  numberOfGaussiansPerClass );

  
protected:
  AtlasMeshToIntensityImageCostAndGradientCalculator();
  virtual ~AtlasMeshToIntensityImageCostAndGradientCalculator();
  
private:
  AtlasMeshToIntensityImageCostAndGradientCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  typedef GMMLikelihoodImageFilter< ImageType >  LikelihoodFilterType;
  
};


} // end namespace kvl

#endif
