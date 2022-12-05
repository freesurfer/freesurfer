#ifndef kvlFrobMMLikelihoodImageFilter_hxx
#define kvlFrobMMLikelihoodImageFilter_hxx

#include "kvlFrobMMLikelihoodImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"
#include <cmath>
//#include <iomanip>

namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
FrobMMLikelihoodImageFilter< TInputImage >
::FrobMMLikelihoodImageFilter()
{
    //m_likelihoodFilter = GMMLikelihoodFilterType::New();
}



#if 0
//
//
//
template< typename TInputImage >
void
FrobMMLikelihoodImageFilter< TInputImage >
::SetGaussianImages( const std::vector< ImageType::ConstPointer >& images )
{
  // it is operating on a different object. The method can be removed???
  //
  for ( unsigned int contrastNumber = 0; contrastNumber < images.size(); contrastNumber++ )
    {
    m_likelihoodFilter->SetInput( contrastNumber, images[ contrastNumber ] );
    }

}
#endif


//
//
//
template< typename TInputImage >
void
FrobMMLikelihoodImageFilter< TInputImage >
::SetParameters( const int numberOfContrasts,
                 const std::vector< double >& frobVariance,
                 const std::vector< vnl_vector< double > >& frobMeans,
                 const std::vector< double >&  frobMixtureWeights,
                 const std::vector< int >&  numberOfFrobeniusPerClass,
                 const double voxratio)
{


   /* 1. m_likelihoodFilter->SetParameters() is updating private variables of a different object
    *      m_Means, m_Precisions, m_piTermMultiv, m_OneOverSqrtDetCov, m_MixtureWeights, m_NumberOfGaussiansPerClass
    * 2. variances, mixtureWeights, numberOfGaussiansPerClass are not used at all other than passed to m_likelihoodFilter->SetParameters()
    *    means is used to get numberOfContrasts
    * 3. m_likelihoodFilter->SetParameters() is not necessary     
    */
    //m_likelihoodFilter->SetParameters(means,
    //                                  variances,
    //                                  mixtureWeights,
    //                                  numberOfGaussiansPerClass);


  // Sanity check on the input parameters
  const size_t  numberOfFrobenius = frobVariance.size();
  if ( numberOfFrobenius == 0 )
    {
    itkExceptionMacro(<< "Empty variances provided" );
    }

  //const size_t  numberOfContrasts = means[ 0 ].size();
  if ( numberOfContrasts == 0 )
    {
    itkExceptionMacro(<< "Empty means provided" );
    }

  size_t  sum = 0;
  for ( size_t i = 0; i < numberOfFrobeniusPerClass.size(); i++ )
    {
    sum += numberOfFrobeniusPerClass[ i ];
    }

  if ( ( frobMeans.size() != numberOfFrobenius ) ||
       ( frobMixtureWeights.size() != numberOfFrobenius ) ||
       ( frobMeans[ 0 ].size()!= 6 ) ||
       ( sum != numberOfFrobenius ) )
    {

      std::cout<<"frob means size = "<<frobMeans.size()<<std::endl;
      std::cout<<"numberOfFrobenius = "<<numberOfFrobenius<<std::endl;
      std::cout<<"frobMixtureWeights size = "<<frobMixtureWeights.size()<<std::endl;
      std::cout<<"frobMeans dimensions = "<<frobMeans[ 0 ].size()<<std::endl;
      std::cout<<"numberOfFrobeniusPerClass sum = "<<sum<<std::endl;


    itkExceptionMacro(<< "Inconsistent parameters" );
    }


  // some parameters are simply copied
  m_means = frobMeans;
  m_frobMixtureWeights = frobMixtureWeights;
  m_NumberOfFrobeniusPerClass = numberOfFrobeniusPerClass;
  m_voxratio = voxratio;


  // Now the variances -- we compute and store precisions instead of variances.
  m_Precisions.resize( numberOfFrobenius );

  // We also calculate the log normaliser
  m_logNormaliser.resize( numberOfFrobenius );

  for(int frobeniusNumber=0; frobeniusNumber< numberOfFrobenius; frobeniusNumber++)
    {

    double precision =  1/(frobVariance[frobeniusNumber]);

    m_Precisions[frobeniusNumber] = precision;

    m_logNormaliser[frobeniusNumber]= 3 * log(2.0* itk::Math::pi);
    m_logNormaliser[frobeniusNumber]+= 3*log(frobVariance[frobeniusNumber]);


    }


  //
  this->Modified();
}




//----------------------------------------------------------------------------
template< typename TInputImage >
void
FrobMMLikelihoodImageFilter< TInputImage >
::BeforeThreadedGenerateData()
{
  // Check to verify all inputs are specified and have the same metadata,
  // spacing etc...
  const unsigned int numberOfInputs = this->GetNumberOfIndexedInputs();
  RegionType         region;

  for ( unsigned int i = 0; i < numberOfInputs; i++ )
    {
    InputImageType *input = itkDynamicCastInDebugMode< InputImageType * >
      (this->itk::ProcessObject::GetInput(i) );
    if ( !input )
      {
      itkExceptionMacro(<< "Input " << i << " not set!");
      }
    if ( i == 0 )
      {
      region = input->GetLargestPossibleRegion();
      }
    else if ( input->GetLargestPossibleRegion() != region )
      {
      itkExceptionMacro(<< "All Inputs must have the same dimensions.");
      }
    }

  // Also check that we have 6 independent tensor variables for 3x3 symmetric log DTI tensors.
  if ( numberOfInputs != 6 )
    {
    itkExceptionMacro(<< "Number of tensor variables don't match number of input channels" )
    }


}

//----------------------------------------------------------------------------
template< typename TInputImage >
void
FrobMMLikelihoodImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  //std::cout << "Executing GMMLikelihoodImageFilter::ThreadedGenerateData()" << std::endl;

  //
  const int  numberOfFrobenius = m_Precisions.size();
  const int  numberOfClasses = m_NumberOfFrobeniusPerClass.size();
  const int  numberOfContrasts = this->GetNumberOfIndexedInputs();
  // std::cout << "numberOfGaussians: " << numberOfGaussians << std::endl;
  // std::cout << "numberOfClasses: " << numberOfClasses << std::endl;
  // std::cout << "numberOfContrasts: " << numberOfContrasts << std::endl;

  //
  itk::ProgressReporter progress( this, threadId, outputRegionForThread.GetNumberOfPixels() );

  // Retreive the output image
  typename OutputImageType::Pointer outputImage =
    static_cast< OutputImageType * >( this->itk::ProcessObject::GetOutput(0) );

  // Initialize iterator over output
  itk::ImageRegionIterator< OutputImageType > oit(outputImage, outputRegionForThread);
  oit.GoToBegin();

  // Als initialize iterator over each input
  typedef itk::ImageRegionConstIterator< InputImageType > InputIteratorType;
  typedef std::vector< InputIteratorType >           InputIteratorContainerType;
  InputIteratorContainerType inputItContainer;
  for ( int contrastNumber = 0; contrastNumber < numberOfContrasts; contrastNumber++ )
    {
    const InputImageType* inputImage = this->GetInput( contrastNumber );

    InputIteratorType iit( inputImage, outputRegionForThread );
    iit.GoToBegin();
    inputItContainer.push_back(iit);
    }

  // Now loop over all pixels
  while ( !oit.IsAtEnd() )
    {
    // Retrieve the input intensity. At the same time, detect the number and pattern of
    // zeroes (interpreted as missing intensities) in the various input channels
    bool isEmpty=true;
    vnl_vector< InputPixelType > aux_v( numberOfContrasts, 0.0 );
    for ( int contrastNumber = 0; contrastNumber < numberOfContrasts; contrastNumber++ )
      {
      const InputPixelType p = inputItContainer[ contrastNumber ].Get();
      ++( inputItContainer[ contrastNumber ] );

      if( p != 0 )
      {
          if( contrastNumber!=0 )
          {
              isEmpty = false;
          }

          aux_v[ contrastNumber ] = p;
      }

      if(std::isinf(p))
      {
          isEmpty = true;
      }
      } // End loop over all contrasts


    // If none of the contrast has any intensity available, fill in the output with
    // some sentinel value, and move on to the next pixel
    if ( isEmpty )
      {
      // Set output pixel and move on to the next pixel
      //std::cout << "Nothing present: " << oit.Value().Size() << std::endl;
      ++oit;
      progress.CompletedPixel();
      continue;
      }


    // Move on with what we actually have
    OutputPixelType pix( numberOfClasses );

    // intensity vector should be of form
    // [log(det(D)),D00,D01,D02,D11,D12,D22]
    vnl_vector< InputPixelType >  intensity_v = aux_v.extract( numberOfContrasts );
    int  shift = 0;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Evaluate the Gaussian mixture model likelihood of this class at the intensity of this pixel
      double  likelihoodAgregate = 0.0;
      double  maxExponent;
      const int  numberOfComponents = m_NumberOfFrobeniusPerClass[ classNumber ];
      for ( int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++ )
        {
          const int  frobNumber = shift + componentNumber;
          double  frob = 0.0;
          const double precision = m_Precisions[frobNumber];
          vnl_vector<double> mean_class = m_means[frobNumber];

          vnl_vector<double> dataV(intensity_v.size());

          for (int dimensionNumber=0; dimensionNumber < 6; dimensionNumber++)
          {
              dataV[dimensionNumber] = intensity_v[dimensionNumber] - mean_class[dimensionNumber];
          }


          // -0.5*(sigma^-2)*(x-m)'*(x-m)
          frob -= 0.5*precision*dot_product(dataV,dataV);

          // normaliser
          frob -= m_logNormaliser[frobNumber];

          // voxratio
          //frob*=m_voxratio;


          //
          const double  mixtureWeight = m_frobMixtureWeights[ frobNumber ];

          if (componentNumber == 0)
          {
              maxExponent = frob;
              likelihoodAgregate = mixtureWeight;
          }
          else if (frob>maxExponent)
          {
              likelihoodAgregate = likelihoodAgregate*exp(maxExponent-frob) + mixtureWeight;
              maxExponent = frob;
          }
          else
          {
              likelihoodAgregate += mixtureWeight*exp(frob-maxExponent);
          }

        } // End loop over components in mixture model for the current class


      //
      pix[ classNumber ] = m_voxratio*(log(likelihoodAgregate) + maxExponent);
      shift += numberOfComponents;

//      if(std::isinf(likelihoodAgregate * exp(maxExponent)))
//      {
//          std::cout<<"agregate = "<<likelihoodAgregate<<std::endl;
//          std::cout<<"maxexp = "<<maxExponent<<std::endl;
//          std::cout<<"exp = "<<exp(maxExponent)<<std::endl;
//          itkExceptionMacro(<<"pix");
//      }

#if 0
      {
      const double d = likelihood;
      uint64_t u;
      memcpy(&u, &d, sizeof(d));
      std::cout << std::hex << u << std::endl;
      }
#endif

      } // End loop over classes

    // Fill in the output pixel and move on
    //std::cout << "pix: " << pix << std::endl;
    oit.Value() = pix;
    ++oit;
    progress.CompletedPixel();
    } // End loop over all pixels
}


} // end namespace kvl

#endif
