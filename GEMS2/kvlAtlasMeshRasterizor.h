#ifndef __kvlAtlasMeshRasterizor_h
#define __kvlAtlasMeshRasterizor_h

#include "kvlAtlasMesh.h"



namespace kvl
{


class AtlasMeshRasterizor : public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef AtlasMeshRasterizor  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >   Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  //itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshRasterizor, itk::Object );

  /** */
  virtual void Rasterize( const AtlasMesh* mesh );

  /** */
  void SetNumberOfThreads( int numberOfThreads )
    { 
    m_NumberOfThreads = numberOfThreads;
    }
  
  /** */
  int GetNumberOfThreads() const
    {
    return m_NumberOfThreads;
    }

protected:
  AtlasMeshRasterizor();
  virtual ~AtlasMeshRasterizor() {};

  /** */
  //
  virtual bool RasterizeTetrahedron( const AtlasMesh* mesh, 
                                     AtlasMesh::CellIdentifier tetrahedronId,
                                     int threadNumber=0 ) = 0;
                                     
  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );
  
  /** Internal structure used for passing information to the threading library */
  struct ThreadStruct
    {
    Pointer  m_Rasterizor;
    AtlasMesh::ConstPointer  m_Mesh;
    std::vector< AtlasMesh::CellIdentifier >  m_TetrahedronIds;
    //std::set< AtlasMesh::CellIdentifier >  m_TetrahedronIds;
    };

                                     

private:
  AtlasMeshRasterizor(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  int  m_NumberOfThreads;
  
};


} // end namespace kvl

#endif
