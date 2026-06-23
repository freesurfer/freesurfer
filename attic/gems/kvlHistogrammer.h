#ifndef __kvlHistogrammer_h
#define __kvlHistogrammer_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkImage.h"

namespace kvl
{


class Histogrammer: public AtlasMeshRasterizor
{
public :
  
  /** Standard class typedefs */
  typedef Histogrammer  Self;
  typedef AtlasMeshRasterizor Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( Histogrammer, itk::Object );

  
  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;
  typedef itk::Image< int, 3 >  BinnedImageType;
  typedef std::vector< std::vector< double > >  HistogramType;
  typedef std::vector< std::vector< ThreadAccumDataType > >  HistogramThreadAccumType;
  typedef std::vector< double >  ConditionalIntensityDistributionType;
  

  /** */  
  void SetImage( const ImageType* image )
    {
    m_Image = image;  
    m_BinnedImage = 0;
    m_NumberOfBins = 0;
    }
      
  /** */
  int  GetNumberOfBins() const
    {
    return m_NumberOfBins;  
    }  
      
  /** */
  void  SetConditionalIntensityDistributions( const std::vector< ConditionalIntensityDistributionType >&
                                              conditionalIntensityDistributions )
    {
    if ( m_BinnedImage )
      {
      // Check if number of bins has changed. If so, forgot cached binned image
      if ( conditionalIntensityDistributions[ 0 ].size() != m_NumberOfBins )
        {
        m_BinnedImage = 0;
        m_NumberOfBins = 0;
        }  
      }  
      
    m_ConditionalIntensityDistributions = conditionalIntensityDistributions;
    
    }
    
  /** */
  const std::vector< ConditionalIntensityDistributionType >&  GetConditionalIntensityDistributions() const
    {
    return m_ConditionalIntensityDistributions;
    }
  
  /** */
  const HistogramType&  GetHistogram() const
    {
    return m_Histogram;
    }

  /** */
  double GetMinLogLikelihood() const
    {
    return m_MinLogLikelihood;
    }

  /** */
  const BinnedImageType* GetBinnedImage() const
    {
    return m_BinnedImage;
    }  
    
  /** */  
  void Rasterize( const AtlasMesh* mesh );
  
  
protected:
  Histogrammer();
  virtual ~Histogrammer();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );
  
private:
  Histogrammer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  void  ComputeRobustRange( const ImageType* image, double& robustMin, double& robustMax );
  
  //
  void  UpdateBinnedImage();
  
  //
  ImageType::ConstPointer  m_Image;
  std::vector< ConditionalIntensityDistributionType >  m_ConditionalIntensityDistributions;
  BinnedImageType::Pointer  m_BinnedImage;
  int  m_NumberOfBins;
  HistogramType  m_Histogram;
  double  m_MinLogLikelihood;

   //
  std::vector< HistogramThreadAccumType >  m_ThreadSpecificHistograms;
  std::vector< ThreadAccumDataType >  m_ThreadSpecificMinLogLikelihoods;

};


} // end namespace kvl

#endif
