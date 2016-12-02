#ifndef __kvlAtlasMeshAlphaDrawerGPU_h
#define __kvlAtlasMeshAlphaDrawerGPU_h

#include "kvlAtlasMeshRasterizorGPU.h"
#include "itkImage.h"


namespace kvl
{


/**
 *
 */
class AtlasMeshAlphaDrawerGPU: public AtlasMeshRasterizorGPU
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshAlphaDrawerGPU  Self;
  typedef AtlasMeshRasterizorGPU Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshAlphaDrawerGPU, itk::Object );

  /** Some typedefs */
  typedef itk::Image< float, 3 >  ImageType;

  /** */
  void SetLabelNumber( unsigned char labelNumber )
    {
    m_LabelNumber = labelNumber;
    }
    
  /** */
  void SetRegions( const ImageType::RegionType&  region )
    {
    m_Image = ImageType::New();
    m_Image->SetRegions( region );
    m_Image->Allocate();
    m_Image->FillBuffer( 0 );
    }
  
  /** */
  const ImageType*  GetImage() const
    { return m_Image; }
    
  
protected:
  AtlasMeshAlphaDrawerGPU();
  virtual ~AtlasMeshAlphaDrawerGPU();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );

private:
  AtlasMeshAlphaDrawerGPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  unsigned char  m_LabelNumber;
  ImageType::Pointer  m_Image;
  
};


} // end namespace kvl

#endif

