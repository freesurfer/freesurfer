#ifndef __kvlAtlasMeshRasterizor_h
#define __kvlAtlasMeshRasterizor_h

#include "kvlAtlasMesh.h"


/*
  If defined, this enables complete reproducibility across
  number of threads. Normally, results are always deterministic
  for a given number of threads, but not across threads, as
  floating-point arithmetic will produce small but cascading
  errors within the per-thread accumulators. However, using a
  high-precision data type (float128) for the thread accumulator
  is enough to keep floating-point rounding errors at bay.
  Enabling this increases average runtime by 20% or more.
*/
#ifdef CROSS_THREAD_REPRODUCIBLE
  #define ThreadAccumDataType __float128
#else
  #define ThreadAccumDataType double
#endif


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
