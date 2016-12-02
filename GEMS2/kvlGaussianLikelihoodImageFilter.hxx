#ifndef kvlGaussianLikelihoodImageFilter_hxx
#define kvlGaussianLikelihoodImageFilter_hxx

#include "kvlGaussianLikelihoodImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "vnl/vnl_inverse.h"


namespace kvl
{
//----------------------------------------------------------------------------
template< typename TInputImage >
GaussianLikelihoodImageFilter< TInputImage >
::GaussianLikelihoodImageFilter()
{
}




//
//
//
template< typename TInputImage >
void
GaussianLikelihoodImageFilter< TInputImage >
::SetPrecisions( const std::vector< vnl_matrix< float > >& precisions )
{
  if ( !precisions.size() )
    {
    this->Modified();
    }
  
  const int  numberOfContrasts = precisions[ 0 ].rows();
  m_Precisions.resize(precisions.size()); 

  // We are going to compute 1/sqrt(det(COV)) for all possible covariances given all possible combinations of available channels
  // We use a binary representation for this. For instance, 6 = [1 1 0] means that we have channel 1 not available, but channels 2 and 3 available.
  m_OneOverSqrtDetCov.resize(precisions.size());
  int nCombos = (int)(pow(2,numberOfContrasts));
  std::vector<bool> presentChannels(numberOfContrasts);
  for(int classNumber=0; classNumber<precisions.size(); classNumber++)
    {
    vnl_matrix<float> FullCov=vnl_inverse<float>(precisions[classNumber]);
    m_OneOverSqrtDetCov[classNumber].resize(nCombos);
    m_Precisions[classNumber].resize(nCombos);
    m_OneOverSqrtDetCov[classNumber][0]=0;
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
      vnl_matrix<float> PartialCov(nPresent,nPresent);
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
      
      m_Precisions[classNumber][n]=vnl_inverse<float>(PartialCov);
      m_OneOverSqrtDetCov[classNumber][n]=1.0/sqrt(vnl_determinant(PartialCov));
      }
    }  

  // We also compute the constant term for number of channels from 0 to numberOfContrasts
  m_piTermMultiv.resize(numberOfContrasts+1);
  for(int i=0; i<=numberOfContrasts; i++) 
    {
    m_piTermMultiv[i] = pow( 2 * itk::Math::pi, -0.5*i );
    }
    
  this->Modified();    
}
  



//----------------------------------------------------------------------------
template< typename TInputImage >
void
GaussianLikelihoodImageFilter< TInputImage >
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
    
  // Also check that sizes of means and precisions match the data
  if ( ( m_Means.size() != m_Precisions.size() ) || 
       ( m_Means[0].size() != numberOfInputs ) )
    {
    itkExceptionMacro(<< "Means and/or precision matrices don't match number of input channels" )  
    }

  
}

//----------------------------------------------------------------------------
template< typename TInputImage >
void
GaussianLikelihoodImageFilter< TInputImage >
::ThreadedGenerateData(const RegionType & outputRegionForThread,
                       itk::ThreadIdType threadId)
{
  //std::cout << "Executing GaussianLikelihoodImageFilter::ThreadedGenerateData()" << std::endl;
  
  //
  const int  numberOfClasses = m_Means.size();
  const int  numberOfContrasts = this->GetNumberOfIndexedInputs();
  
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
    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
      // Evaluate the Gaussian of this class at the intensity of this pixel
      double gauss = 0.0;
      if( numberOfContrasts == 1 )
        {
        gauss = exp( -0.5 * m_Precisions[classNumber][1][0][0] * pow( intensity_v[ 0 ] - m_Means[ classNumber ][0] , 2 )) / 
                                    sqrt(  2 * itk::Math::pi / m_Precisions[ classNumber ][1][0][0] );
        }
      else
        {
        
        vnl_vector<float>  dataV(intensity_v.size());
        int c=0;
        for( int contrastNumber=0; contrastNumber < numberOfContrasts; contrastNumber++)
          {
          if(isThere[ contrastNumber ])
            {
            dataV[c]=intensity_v[c]-m_Means[classNumber][ contrastNumber ];
            c++;
            }
          }

        gauss= exp( -0.5 * dot_product(dataV,m_Precisions[classNumber][index]*dataV)) * m_OneOverSqrtDetCov[classNumber][index] * m_piTermMultiv[nPresent];

        } // End test number of contrasts
      
      pix[ classNumber ] = gauss;
      } // End loop over contrasts  
          
    // Fill in the output pixel and move on
    //std::cout << "pix: " << pix << std::endl;
    oit.Value() = pix;
    ++oit;
    progress.CompletedPixel();
    } // End loop over all pixels
}

} // end namespace kvl

#endif
