#ifndef __kvlAtlasMeshCollectionModelLikelihoodCalculator_h
#define __kvlAtlasMeshCollectionModelLikelihoodCalculator_h

#include "kvlAtlasMeshCollection.h"
#include "kvlCompressionLookupTable.h"


namespace kvl
{


class AtlasMeshCollectionModelLikelihoodCalculator: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshCollectionModelLikelihoodCalculator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionModelLikelihoodCalculator, itk::Object );

  /** Some typedefs */
  typedef CompressionLookupTable::ImageType  LabelImageType;

  // Set label images.
  void SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages,
                       const CompressionLookupTable*  compressionLookupTable );

  // Get label images
  const std::vector< LabelImageType::ConstPointer >&  GetLabelImages() const
    { return m_LabelImages; }

  // Get label image
  const LabelImageType*  GetLabelImage( unsigned int labelImageNumber ) const;

  // Initialize
  void SetMeshCollection( const AtlasMeshCollection* meshCollection )
    { m_MeshCollection = meshCollection; }

  //
  const AtlasMeshCollection*  GetMeshCollection() const
    { return m_MeshCollection; }

  /** */
  void GetDataCostAndAlphasCost( double& dataCost, double& alphasCost, bool verbose=false ) const;
  
protected:
  AtlasMeshCollectionModelLikelihoodCalculator();
  virtual ~AtlasMeshCollectionModelLikelihoodCalculator();
  
private:
  AtlasMeshCollectionModelLikelihoodCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  AtlasMeshCollection::ConstPointer  m_MeshCollection;

  std::vector< LabelImageType::ConstPointer >  m_LabelImages;
  CompressionLookupTable::ConstPointer  m_CompressionLookupTable;
  
  
};




} // end namespace kvl

#endif

