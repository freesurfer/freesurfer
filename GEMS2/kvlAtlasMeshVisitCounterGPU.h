#ifndef __kvlAtlasMeshVisitCounterGPU_h
#define __kvlAtlasMeshVisitCounterGPU_h

#include "kvlAtlasMeshRasterizorGPU.h"
#include "itkImage.h"


namespace kvl
{


/**
 *
 */
class AtlasMeshVisitCounterGPU: public AtlasMeshRasterizorGPU
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshVisitCounterGPU  Self;
  typedef AtlasMeshRasterizorGPU Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshVisitCounterGPU, itk::Object );

  /** Some typedefs */
  typedef itk::Image< int, 3 >  ImageType;

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
  AtlasMeshVisitCounterGPU();
  virtual ~AtlasMeshVisitCounterGPU();
  
  //
  bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                             AtlasMesh::CellIdentifier tetrahedronId,
                             int threadNumber );

private:
  AtlasMeshVisitCounterGPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  //
  ImageType::Pointer  m_Image;
  
};



} // end namespace kvl

#endif
