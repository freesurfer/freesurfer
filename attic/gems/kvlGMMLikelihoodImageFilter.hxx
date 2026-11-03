#ifndef kvlGMMLikelihoodImageFilter_hxx
#define kvlGMMLikelihoodImageFilter_hxx

#include "kvlGMMLikelihoodImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"
//#include <iomanip>

namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
GMMLikelihoodImageFilter< TInputImage >
::GMMLikelihoodImageFilter()
{
#if ITK_VERSION_MAJOR >= 5
  // use classic void ThreadedGenerateData( const OutputRegionType& threadRegion, ThreadIdType threadId )
  // instead of default new signature void DynamicThreadedGenerateData( const OutputRegionType& threadRegion )
  this->DynamicMultiThreadingOff();
#endif  
}




//
//
//
template< typename TInputImage >
void
GMMLikelihoodImageFilter< TInputImage >
::SetParameters( const std::vector< vnl_vector< double > >& means, 
                 const std::vector< vnl_matrix< double > >& variances,
                 const std::vector< double >&  mixtureWeights,
                 const std::vector< int >&  numberOfGaussiansPerClass )
{
  
  // Sanity check on the input parameters
  const int  numberOfGaussians = means.size();
  if ( numberOfGaussians == 0 )
    {
    itkExceptionMacro(<< "Empty means provided" );  
    }

  const int  numberOfContrasts = means[ 0 ].size();
  if ( numberOfContrasts == 0 )
    {
    itkExceptionMacro(<< "Empty means provided" ); 
    }

  int  sum = 0;
  for ( int i = 0; i < numberOfGaussiansPerClass.size(); i++ )
    {
    sum += numberOfGaussiansPerClass[ i ];  
    }  
    
  if ( ( variances.size() != numberOfGaussians ) ||
       ( mixtureWeights.size() != numberOfGaussians ) ||
       ( variances[ 0 ].rows() != numberOfContrasts ) ||
       ( sum != numberOfGaussians ) )
    {
    itkExceptionMacro(<< "Inconsistent parameters" );  
    }
    
  
  
  
  
  
  // All parameters except for variances are simply copied
  m_Means = means;
  m_MixtureWeights = mixtureWeights;
  m_NumberOfGaussiansPerClass = numberOfGaussiansPerClass;

  // Now the variances -- we compute and store precisions instead of variances.
  // In addition, we allow for certain contrasts to be present and others not in a single pixel --
  // we precompute (and store) everything that's need to efficiently evaluate the GMM likelihood in
  // such cases
  m_Precisions.resize( numberOfGaussians ); 

  // We are going to compute 1/sqrt(det(COV)) for all possible covariances given all possible combinations of available channels
  // We use a binary representation for this. For instance, 6 = [1 1 0] means that we have channel 1 not available, but channels 2 and 3 available.
  m_OneOverSqrtDetCov.resize( numberOfGaussians );
  int nCombos = (int)(pow(2,numberOfContrasts));
  std::vector<bool> presentChannels(numberOfContrasts);
  for(int gaussianNumber=0; gaussianNumber< numberOfGaussians; gaussianNumber++)
    {
    vnl_matrix<double>  FullCov = variances[ gaussianNumber ];
    m_OneOverSqrtDetCov[gaussianNumber].resize(nCombos);
    m_Precisions[gaussianNumber].resize(nCombos);
    m_OneOverSqrtDetCov[gaussianNumber][0]=0;
    for(int n=1; n<nCombos; n++) 
      {
      // decode integer -> binary vector of present channels
      int k = n;
      int nPresent=0;
      for(int c=0; c<numberOfContrasts; c++)
        {
        if (k & 1) {presentChannels[c]=true; nPresent++; }
        else {presentChannels[c]=false;}
        k = (k >> 1);
        }

      // Extract sub-matrix
      vnl_matrix<double> PartialCov(nPresent,nPresent);
      int r=0; 
      for(int i=0; i<numberOfContrasts; i++)
        {
        if(presentChannels[i])
          {
          // copy from row i to row r              
          int c=0;
          for(int j=0; j<numberOfContrasts; j++)
            {
            if(presentChannels[j])
              {
              // copy from i,j to r,c
              PartialCov[r][c]=FullCov[i][j];
              c++;
              }
            }

          r++;
          }
        }
      
      m_Precisions[gaussianNumber][n]=vnl_inverse<double>(PartialCov);
      m_OneOverSqrtDetCov[gaussianNumber][n]=1.0/sqrt(vnl_determinant(PartialCov));
      }
    }  

  // We also compute the constant term for number of channels from 0 to numberOfContrasts
  m_piTermMultiv.resize(numberOfContrasts+1);
  for(int i=0; i<=numberOfContrasts; i++) 
    {
    m_piTermMultiv[i] = pow( 2 * itk::Math::pi, -0.5*i );
    std::cout << m_piTermMultiv[i] << std::endl;
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
GMMLikelihoodImageFilter< TInputImage >
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
    
  // Also check that the parameters match the data
  if ( m_Means[0].size() != numberOfInputs )
    {
    itkExceptionMacro(<< "Parameters don't match number of input channels" )  
    }

  
}


//----------------------------------------------------------------------------
template< typename TInputImage >
void
GMMLikelihoodImageFilter< TInputImage >
::BeforeThreadedGenerateData(const RegionType region)
{
  // Check to verify all inputs are specified and have the same metadata,
  // spacing etc...
  const unsigned int numberOfInputs = this->GetNumberOfIndexedInputs();

  for ( unsigned int i = 0; i < numberOfInputs; i++ )
    {
    InputImageType *input = itkDynamicCastInDebugMode< InputImageType * >
      (this->itk::ProcessObject::GetInput(i) );
    if ( !input )
      {
      itkExceptionMacro(<< "Input " << i << " not set!");
      }
    else if ( input->GetLargestPossibleRegion() != region )
      {
      itkExceptionMacro(<< "All Inputs must have the same dimensions.");
      }
    }

  // Also check that the parameters match the data
  if ( m_Means[0].size() != numberOfInputs )
    {
    itkExceptionMacro(<< "Parameters don't match number of input channels" )
    }


}

//----------------------------------------------------------------------------
template< typename TInputImage >
void
GMMLikelihoodImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  //std::cout << "Executing GMMLikelihoodImageFilter::ThreadedGenerateData()" << std::endl;
  
  //
  const int  numberOfGaussians = m_Means.size();
  const int  numberOfClasses = m_NumberOfGaussiansPerClass.size();
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
    std::vector< bool >  isThere( numberOfContrasts );
    int nPresent=0;
    int index=0;
    int aux=1;
    vnl_vector< InputPixelType > aux_v( numberOfContrasts, 0.0 );
    for ( int contrastNumber = 0; contrastNumber < numberOfContrasts; contrastNumber++ )
      {
      const InputPixelType p = inputItContainer[ contrastNumber ].Get();
      ++( inputItContainer[ contrastNumber ] );

      if( p != 0 )
        {
        isThere[ contrastNumber ] = true;
        aux_v[ nPresent ] = p;
        nPresent++;
        index += aux;
        }
      else
        {
        isThere[ contrastNumber ] = false;
        }
      aux = aux << 1;
      } // End loop over all contrasts
      
      
    // If none of the contrast has any intensity available, fill in the output with
    // some sentinel value, and move on to the next pixel
    if ( nPresent==0 )
      {
      // Set output pixel and move on to the next pixel
      //std::cout << "Nothing present: " << oit.Value().Size() << std::endl;
      ++oit;
      progress.CompletedPixel();
      continue;
      }

      
    // Move on with what we actually have
    OutputPixelType pix( numberOfClasses );
    vnl_vector< InputPixelType >  intensity_v = aux_v.extract( nPresent ); 
    int  shift = 0;
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Evaluate the Gaussian mixture model likelihood of this class at the intensity of this pixel
      double  likelihood = 0.0;
      const int  numberOfComponents = m_NumberOfGaussiansPerClass[ classNumber ];
      for ( int componentNumber = 0; componentNumber < numberOfComponents; componentNumber++ )
        {
        const int  gaussianNumber = shift + componentNumber;
        double  gauss = 0.0;
        if ( numberOfContrasts == 1 )
          {
          gauss = exp( -0.5 * m_Precisions[ gaussianNumber ][1][0][0] * pow( intensity_v[ 0 ] - m_Means[ gaussianNumber ][0] , 2 )) / 
                                      sqrt(  2 * itk::Math::pi / m_Precisions[ gaussianNumber ][1][0][0] );
          }
        else
          {
          
          vnl_vector<double>  dataV(intensity_v.size());
          int c=0;
          for( int contrastNumber=0; contrastNumber < numberOfContrasts; contrastNumber++)
            {
            if( isThere[ contrastNumber ] )
              {
              dataV[c]=intensity_v[c]-m_Means[ gaussianNumber ][ contrastNumber ];
              c++;
              }
            }

          gauss= exp( -0.5 * dot_product(dataV,m_Precisions[ gaussianNumber ][index]*dataV)) * m_OneOverSqrtDetCov[ gaussianNumber ][index] * m_piTermMultiv[nPresent];

          } // End test number of contrasts
     
        //
        const double  mixtureWeight = m_MixtureWeights[ gaussianNumber ];
        likelihood += gauss * mixtureWeight;
     
        } // End loop over components in mixture model for the current class
     
     
      //
      pix[ classNumber ] = likelihood;
      shift += numberOfComponents;

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
