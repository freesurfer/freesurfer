#ifndef __kvlAtlasMeshMultiAlphaDrawer_h
#define __kvlAtlasMeshMultiAlphaDrawer_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkImage.h"


namespace kvl
{


class AtlasMeshMultiAlphaDrawer: public AtlasMeshRasterizor
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshMultiAlphaDrawer  Self;
  typedef AtlasMeshRasterizor Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshMultiAlphaDrawer, itk::Object );

  /** Some typedefs */
  typedef itk::Image< AtlasAlphasType, 3 >  ImageType;

  /** */
  void SetRegions( const ImageType::RegionType&  region )
    {
    m_Image = ImageType::New();
    m_Image->SetRegions( region );
    m_Image->Allocate();
    }
  
  /** */
  const ImageType*  GetImage() const
    { return m_Image; }
  
  //
  void Rasterize( const AtlasMesh* mesh );

  
protected:
  AtlasMeshMultiAlphaDrawer();
  virtual ~AtlasMeshMultiAlphaDrawer();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );

private:
  AtlasMeshMultiAlphaDrawer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  ImageType::Pointer  m_Image;
  
};


} // end namespace kvl

#endif
