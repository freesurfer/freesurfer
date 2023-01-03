#ifndef kvlDSWbetaMMLikelihoodImageFilter_hxx
#define kvlDSWbetaMMLikelihoodImageFilter_hxx

#include "kvlDSWbetaMMLikelihoodImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"
#include <cmath>
//#include <iomanip>

namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
DSWbetaMMLikelihoodImageFilter< TInputImage >
::DSWbetaMMLikelihoodImageFilter()
{
    //m_likelihoodFilter = GMMLikelihoodFilterType::New();
}



#if 0
//
//
//
template< typename TInputImage >
void
DSWbetaMMLikelihoodImageFilter< TInputImage >
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
DSWbetaMMLikelihoodImageFilter< TInputImage >
::SetParameters( const int numberOfContrasts,
                 const std::vector< double >& dswbetaAlpha,
                 const std::vector< double >& dswbetaBeta,
                 const std::vector< double >& dswbetaConcentration,
                 const std::vector< vnl_vector< double > >& dswbetaMeans,
                 const std::vector< double >& logKummerSamples,
                 const double& negLogKummerIncrement,
                 const std::vector< double >&  dswbetaMixtureWeights,
                 const std::vector< int >&  numberOfDSWbetaPerClass,
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
  const size_t  numberOfDSWbeta = dswbetaConcentration.size();
  if ( numberOfDSWbeta == 0 )
    {
    itkExceptionMacro(<< "Empty Concentrations provided" );
    }

  //const size_t  numberOfContrasts = means[ 0 ].size();
  if ( numberOfContrasts == 0 )
    {
    itkExceptionMacro(<< "Empty means provided" );
    }

  size_t  sum = 0;
  for ( size_t i = 0; i < numberOfDSWbetaPerClass.size(); i++ )
    {
    sum += numberOfDSWbetaPerClass[ i ];
    }

  if ( ( dswbetaMeans.size() != numberOfDSWbeta ) ||
       ( dswbetaMixtureWeights.size() != numberOfDSWbeta ) ||
       ( dswbetaAlpha.size()!= numberOfDSWbeta ) ||
       ( dswbetaBeta.size()!= numberOfDSWbeta ) ||
       ( dswbetaMeans[ 0 ].size()!= 3 ) ||
       ( sum != numberOfDSWbeta ) )
    {

      std::cout<<"dswbeta alphas size = "<<dswbetaAlpha.size()<<std::endl;
      std::cout<<"dswbeta betas size = "<<dswbetaBeta.size()<<std::endl;
      std::cout<<"dswbeta means size = "<<dswbetaMeans.size()<<std::endl;
      std::cout<<"numberOfDSWbeta = "<<numberOfDSWbeta<<std::endl;
      std::cout<<"dswbetaMixtureWeights size = "<<dswbetaMixtureWeights.size()<<std::endl;
      std::cout<<"dswbetaMeans dimensions = "<<dswbetaMeans[ 0 ].size()<<std::endl;
      std::cout<<"numberOfDSWbetaPerClass sum = "<<sum<<std::endl;


    itkExceptionMacro(<< "Inconsistent parameters" );
    }


  // some parameters are simply copied
  m_Alphas = dswbetaAlpha;
  m_Betas = dswbetaBeta;
  m_Concentrations = dswbetaConcentration;
  m_means = dswbetaMeans;
  m_logNormaliserSamples = logKummerSamples;
  m_logNormaliserIncrement = negLogKummerIncrement;
  m_DSWbetaMixtureWeights = dswbetaMixtureWeights;
  m_NumberOfDSWbetaPerClass = numberOfDSWbetaPerClass;
  m_voxratio = voxratio;


  // We also calculate the beta function
  m_LogBetaFunctions.resize( numberOfDSWbeta );

  for(int dswbetaNumber=0; dswbetaNumber< numberOfDSWbeta; dswbetaNumber++)
    {
      const double alpha = dswbetaAlpha[dswbetaNumber];
      const double beta = dswbetaBeta[dswbetaNumber];


      double betaFunction =  (std::lgamma(alpha)+std::lgamma(beta))-(std::lgamma(alpha+beta));

      m_LogBetaFunctions[dswbetaNumber] = betaFunction;

      //std::cout << "m_LogBetaFunctions[" << dswbetaNumber <<"] = "<< m_LogBetaFunctions[dswbetaNumber] << std::endl;


    }


  //
  this->Modified();
}




//----------------------------------------------------------------------------
template< typename TInputImage >
void
DSWbetaMMLikelihoodImageFilter< TInputImage >
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

  // Also check that we have 4 channels for FA and V1 eigenvector.
  if ( numberOfInputs != 4 )
    {
    itkExceptionMacro(<< "Number of tensor variables don't match number of input channels" )
    }


}

//----------------------------------------------------------------------------
template< typename TInputImage >
void
DSWbetaMMLikelihoodImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  //std::cout << "Executing GMMLikelihoodImageFilter::ThreadedGenerateData()" << std::endl;

  //
  const int  numberOfDSWbeta = m_Concentrations.size();
  const int  numberOfClasses = m_NumberOfDSWbetaPerClass.size();
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
          isEmpty = false;

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
    // [FA,V1,V2,V3]
    vnl_vector< InputPixelType >  intensity_v = aux_v.extract( numberOfContrasts );
    int  shift = 0;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Evaluate the Gaussian mixture model likelihood of this class at the intensity of this pixel
      double  likelihoodAgregate = 0.0;
      double  maxExponent;
      const int  numberOfComponents = m_NumberOfDSWbetaPerClass[ classNumber ];
      for ( int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++ )
        {
          const int  dswbetaNumber = shift + componentNumber;
          double  dswbeta = 0.0;
          double  logpdf_dsw = 0.0;
          double  logpdf_beta = 0.0;
          double  vectinnerprod = 0.0;
          double  effectiveConcentration = 0.0;
          double  vint = 0.0;
          double  diff = 0.0;
          double  prop = 0.0;
          double  logNormaliser = 0.0;
          int     left = 0;
          int     right = 0;
          const double alpha_prm = m_Alphas[dswbetaNumber];
          const double beta_prm = m_Betas[dswbetaNumber];
          const double logBet_fnc = m_LogBetaFunctions[dswbetaNumber];
          const double concentration = m_Concentrations[dswbetaNumber];
          vnl_vector<double> mean_class = m_means[dswbetaNumber];
          const double  mixtureWeight = m_DSWbetaMixtureWeights[ dswbetaNumber ];


          logpdf_beta = (alpha_prm-1.0)*log(intensity_v[0])
                  + (beta_prm-1.0)*log(1.0-intensity_v[0]);

          logpdf_beta -= logBet_fnc;

          for (int dimensionNumber=0; dimensionNumber < 3; dimensionNumber++)
          {
              vectinnerprod += intensity_v[dimensionNumber+1] * mean_class[dimensionNumber];
          }

          effectiveConcentration = concentration*intensity_v[0];

          logpdf_dsw += effectiveConcentration*(vectinnerprod*vectinnerprod);

          vint = effectiveConcentration/m_logNormaliserIncrement;

          left = floor(vint);
          right = ceil(vint);

          // interpolate the log Kummer function
          if (right>m_logNormaliserSamples.size()-1) {
              // extrapolate in case of very large concentrations
              left = m_logNormaliserSamples.size()-2;
              right = m_logNormaliserSamples.size()-1;

          }

          diff = m_logNormaliserSamples[right] - m_logNormaliserSamples[left];

          prop = vint - double(left);

          logNormaliser = m_logNormaliserSamples[left] + prop*diff;


          dswbeta += logpdf_beta + logpdf_dsw;

          // normaliser
          dswbeta -= logNormaliser;

          // voxratio
          // dswbeta *= m_voxratio;


          // add to class mixture in voxel
          if (componentNumber == 0)
          {
              maxExponent = dswbeta;
              likelihoodAgregate = mixtureWeight;
          }
          else if (dswbeta>maxExponent)
          {
              likelihoodAgregate = likelihoodAgregate*exp(maxExponent-dswbeta) + mixtureWeight;
              maxExponent = dswbeta;
          }
          else
          {
              likelihoodAgregate += mixtureWeight*exp(dswbeta-maxExponent);
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
