#ifndef __kvlAtlasMeshRasterizor_h
#define __kvlAtlasMeshRasterizor_h

#include "kvlAtlasMesh.h"
#include "itkImage.h"



namespace kvl
{


template < class TFragmentProcessor >
class AtlasMeshRasterizor : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef AtlasMeshRasterizor  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >   Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshRasterizor, itk::Object );

  /** Some convenient typedefs. */
  typedef TFragmentProcessor   FragmentProcessorType;
  typedef itk::Image< unsigned char, 3 >  LabelImageType;

  /** Get one of the FragmentProcessor objects. */
  FragmentProcessorType& GetFragmentProcessor( int threadNumber=0 ) 
    { 
    return m_FragmentProcessors[ threadNumber ]; 
    }
  
  const FragmentProcessorType& GetFragmentProcessor( int threadNumber=0 ) const 
    { 
    return m_FragmentProcessors[ threadNumber ]; 
    }

  /** Get all the FragmentProcessor objects. */
  std::vector< FragmentProcessorType >& GetFragmentProcessors()
    {
    return m_FragmentProcessors;
    }
  const std::vector< FragmentProcessorType >& GetFragmentProcessors() const
    {
    return m_FragmentProcessors;
    }

  
  /** */  
  virtual void SetLabelImage( const LabelImageType*  labelImage )
    {
    m_LabelImage = labelImage;
    this->Modified();
    }
    
  /** */
  const LabelImageType*  GetLabelImage() const
    {
    return m_LabelImage;
    }
  
  /** */
  virtual void Rasterize( const AtlasMesh* mesh, bool useMultiThreading=false );


  /** */
  void RasterizeTetrahedron( const AtlasMesh::PointType& p0, const AtlasMesh::PointType& p1,  
                             const AtlasMesh::PointType& p2, const AtlasMesh::PointType& p3 )
    {
    Self::RasterizeTetrahedron( p0, p1, p2, p3, m_FragmentProcessors[ 0 ], m_LabelImage );
    }

protected:
  AtlasMeshRasterizor();
  virtual ~AtlasMeshRasterizor() {};

  /** */
  static void RasterizeTetrahedron( const AtlasMesh::PointType& p0, const AtlasMesh::PointType& p1,  
                                    const AtlasMesh::PointType& p2, const AtlasMesh::PointType& p3, 
                                    FragmentProcessorType& fragmentProcessor, const LabelImageType* labelImage );

  /** */
  static void RasterizeTriangle( const float* vertex0, const float* vertex1, const float* vertex2, 
                          const float* pisIn0, const float* pisIn1, const float* pisIn2, 
                          const int zLevel,
                          FragmentProcessorType& fragmentProcessor, const LabelImageType* labelImage );

private:
  AtlasMeshRasterizor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  
  /** Internal structure used for passing information to the threading library */
  struct ThreadStruct
    {
    LabelImageType::ConstPointer  m_LabelImage;
    AtlasMesh::ConstPointer  m_Mesh;
    std::vector< FragmentProcessorType >*  m_FragmentProcessors;
    std::set< AtlasMesh::CellIdentifier >  m_TetrahedronsToRasterize;
    };

  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
  
  
  //
  static bool RasterizedTetrahedra( const LabelImageType* labelImage, const AtlasMesh* mesh,
                                    FragmentProcessorType& fragmentProcessor, 
                                    std::set< AtlasMesh::CellIdentifier >&  tetrahedronsToRasterize, 
                                    int numberOfTetrahedra, int threadId );

  
  std::vector< FragmentProcessorType >  m_FragmentProcessors;
  LabelImageType::ConstPointer  m_LabelImage;

  
};


} // end namespace kvl

#include "kvlAtlasMeshRasterizor.txx"

#endif
