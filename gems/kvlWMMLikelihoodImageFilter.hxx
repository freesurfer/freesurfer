#ifndef kvlWMMLikelihoodImageFilter_hxx
#define kvlWMMLikelihoodImageFilter_hxx

#include "kvlWMMLikelihoodImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"
#include <cmath>
//#include <iomanip>

namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
WMMLikelihoodImageFilter< TInputImage >
::WMMLikelihoodImageFilter()
{
    //m_likelihoodFilter = GMMLikelihoodFilterType::New();
}



#if 0
//
//
//
template< typename TInputImage >
void
WMMLikelihoodImageFilter< TInputImage >
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
WMMLikelihoodImageFilter< TInputImage >
::SetParameters( const int numberOfContrasts,
                 const std::vector< double >& degreesOfFreedom,
                 const std::vector< vnl_matrix< double > >& scaleMatrices,
                 const std::vector< double >&  wmmMixtureWeights,
                 const std::vector< int >&  numberOfWishartsPerClass,
                 const double& voxratio)
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
  const size_t  numberOfWisharts = degreesOfFreedom.size();
  if ( numberOfWisharts == 0 )
    {
    itkExceptionMacro(<< "Empty degreesOfFreedom provided" );
    }

  //const size_t  numberOfContrasts = means[ 0 ].size();
  if ( numberOfContrasts == 0 )
    {
    itkExceptionMacro(<< "Empty means provided" );
    }

  size_t  sum = 0;
  for ( size_t i = 0; i < numberOfWishartsPerClass.size(); i++ )
    {
    sum += numberOfWishartsPerClass[ i ];
    }

  if ( ( scaleMatrices.size() != numberOfWisharts ) ||
       ( wmmMixtureWeights.size() != numberOfWisharts ) ||
       ( scaleMatrices[ 0 ].rows() != 3 ) ||
       ( sum != numberOfWisharts ) )
    {
    itkExceptionMacro(<< "Inconsistent parameters" );
    }


  // some parameters are simply copied
  m_wmmMixtureWeights = wmmMixtureWeights;
  m_NumberOfWishartsPerClass = numberOfWishartsPerClass;
  m_voxratio = voxratio;

  // Now the scaleMatrices -- we compute and store negative half inverse for quick trace calculation.
  m_negInverseScaleOverTwo.resize( numberOfWisharts );

  // We also calculate the log normaliser including gamma function evaluation and determinant.
  m_logNormaliser.resize( numberOfWisharts );

  // Finally we calculate the exponent for the data determinant
  m_degreesOfFreedomExponent.resize( numberOfWisharts);

  for(int wishartNumber=0; wishartNumber< numberOfWisharts; wishartNumber++)
    {
    vnl_matrix<double>  FullScale = scaleMatrices[ wishartNumber ];
    double dof = degreesOfFreedom[wishartNumber];

    m_negInverseScaleOverTwo[wishartNumber]=-0.5*vnl_inverse<double>(FullScale);



    m_logNormaliser[wishartNumber]= 1.5 * dof *log(2.0);
    m_logNormaliser[wishartNumber]+= lgamma(0.5*dof) + lgamma(0.5*(dof-1.0))
            + lgamma(0.5*(dof-2.0)) + ( 1.5 * log( itk::Math::pi)) ;
    m_logNormaliser[wishartNumber]+= 0.5*dof*log(vnl_determinant(FullScale));

    m_degreesOfFreedomExponent[wishartNumber]= 0.5*(dof-4);
#if 0


    std::cout<<"m_negInverseScaleOverTwo[1][1] = "<<m_negInverseScaleOverTwo[wishartNumber][0][0]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[1][2] = "<<m_negInverseScaleOverTwo[wishartNumber][0][1]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[1][3] = "<<m_negInverseScaleOverTwo[wishartNumber][0][2]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[2][1] = "<<m_negInverseScaleOverTwo[wishartNumber][1][0]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[2][2] = "<<m_negInverseScaleOverTwo[wishartNumber][1][1]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[2][3] = "<<m_negInverseScaleOverTwo[wishartNumber][1][2]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[3][1] = "<<m_negInverseScaleOverTwo[wishartNumber][2][0]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[3][2] = "<<m_negInverseScaleOverTwo[wishartNumber][2][1]<<std::endl;
    std::cout<<"m_negInverseScaleOverTwo[3][3] = "<<m_negInverseScaleOverTwo[wishartNumber][2][2]<<std::endl;

    std::cout<<"m_logNormaliser = "<<m_logNormaliser[wishartNumber]<<std::endl;

    std::cout<<"m_degreesOfFreedomExponent = "<<m_degreesOfFreedomExponent[wishartNumber]<<std::endl;
#endif

    }

#if 0
  //std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
  //std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 2);
  //std::cout << std::scientific;
  for ( int gaussianNumber = 0; gaussianNumber < numberOfGaussians; gaussianNumber++ )
    {
    std::cout << "m_Precisions[ " << gaussianNumber << " ]: " << std::endl;
    for ( int row = 0; row < 2; row++ )
      {
      std::cout << "      ";
      for ( int col = 0; col < 2; col++ )
        {
        // std::cout << m_Precisions[ gaussianNumber ][3][row][col] << " ";
        const double d = m_Precisions[ gaussianNumber ][3][row][col];
        uint64_t u;
        memcpy(&u, &d, sizeof(d));
        std::cout << std::hex << u << " ";
        }
      std::cout << std::endl;
      }

    // std::cout << "m_OneOverSqrtDetCov: " << m_OneOverSqrtDetCov[ gaussianNumber ][ 3 ] << std::endl;

    }

#endif

  //
  this->Modified();
}




//----------------------------------------------------------------------------
template< typename TInputImage >
void
WMMLikelihoodImageFilter< TInputImage >
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

  // Also check that we have 6 independent tensor variables for 3x3 symmetric DTI tensors
  // and one expensive log determinant.
  if ( numberOfInputs != 7 )
    {
    itkExceptionMacro(<< "Number of tensor variables don't match number of input channels" )
    }


}

//----------------------------------------------------------------------------
template< typename TInputImage >
void
WMMLikelihoodImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  //std::cout << "Executing GMMLikelihoodImageFilter::ThreadedGenerateData()" << std::endl;

  //
  const int  numberOfWisharts = m_degreesOfFreedomExponent.size();
  const int  numberOfClasses = m_NumberOfWishartsPerClass.size();
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
      const int  numberOfComponents = m_NumberOfWishartsPerClass[ classNumber ];
      for ( int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++ )
        {
          const int  wishartNumber = shift + componentNumber;
          double  wish = 0.0;
          const vnl_matrix< double > negInv = m_negInverseScaleOverTwo[wishartNumber];

          // 0.5*(dof-4)*log(det(D))
          wish += m_degreesOfFreedomExponent[wishartNumber]*intensity_v[0];

          // -0.5*trace(W^(-1)*D)
          wish += negInv[0][0]*intensity_v[1] + ((negInv[0][1]+negInv[1][0])*intensity_v[2])
                  + ((negInv[0][2]+negInv[2][0])*intensity_v[3])
                  + negInv[1][1]*intensity_v[4] + ((negInv[1][2]+negInv[2][1])*intensity_v[5])
                  + negInv[2][2]*intensity_v[6];

          // normaliser
          wish -= m_logNormaliser[wishartNumber];

          //wish*= m_voxratio;


          //
          const double  mixtureWeight = m_wmmMixtureWeights[ wishartNumber ];

          if (componentNumber == 0)
          {
              maxExponent = wish;
              likelihoodAgregate = mixtureWeight;
          }
          else if (wish>maxExponent)
          {
              likelihoodAgregate = likelihoodAgregate*exp(maxExponent-wish) + mixtureWeight;
              maxExponent = wish;
          }
          else
          {
              likelihoodAgregate += mixtureWeight*exp(wish-maxExponent);
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
