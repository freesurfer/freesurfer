#ifndef __kvlAtlasMeshValueDrawer_h
#define __kvlAtlasMeshValueDrawer_h

#include "kvlAtlasMeshRasterizor.h"
#include "itkImage.h"


namespace kvl {


class AtlasMeshValueDrawer: public AtlasMeshRasterizor
{
public :
  AtlasMeshValueDrawer()  {}
  ~AtlasMeshValueDrawer() {}

  // standard class typedefs
  typedef AtlasMeshValueDrawer Self;
  typedef AtlasMeshRasterizor Superclass;
  typedef itk::Array<double> AtlasValuesType;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  // method for creation through the object factory
  itkNewMacro(Self);

  // run-time type information (and related methods)
  itkTypeMacro(AtlasMeshValueDrawer, itk::Object);

  // image type
  typedef itk::Image<AtlasValuesType, 3>  ImageType;

  // set image region
  void SetRegions(const ImageType::RegionType& region, int nframes) {
    m_Image = ImageType::New();
    m_Image->SetRegions(region);
    m_Image->Allocate();
    AtlasValuesType empties(nframes);
    empties.Fill(0);
    m_Image->FillBuffer(empties);
    m_NumFrames = nframes;
  }
  
  // set point values
  void SetValues(double const * values) { m_Values = values; }

  // return internal image
  const ImageType* GetImage() const { return m_Image; }

protected:

  // AtlasMeshValueDrawer()  {}
  // ~AtlasMeshValueDrawer() {}

  bool RasterizeTetrahedron(const AtlasMesh* mesh, AtlasMesh::CellIdentifier tetrahedronId, int threadNumber);

private:
  AtlasMeshValueDrawer(const Self&);  // purposely not implemented
  void operator=(const Self&);        // purposely not implemented
  
  int m_NumFrames = 0;
  double const * m_Values = 0;
  ImageType::Pointer m_Image = nullptr;
  
};


} // end namespace kvl

#endif
