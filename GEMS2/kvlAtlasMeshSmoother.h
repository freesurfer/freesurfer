#ifndef __kvlAtlasMeshSmoother_h
#define __kvlAtlasMeshSmoother_h

#include "kvlAtlasMeshCollection.h"


namespace kvl
{


class AtlasMeshSmoother: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshSmoother  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshSmoother, itk::Object );
  
  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
    { 
    m_MeshCollection = meshCollection;   
    }

  //
  const AtlasMeshCollection* GetMeshCollection() const
    { return m_MeshCollection; }

  // 
  void  SetSigma( float sigma )
    { m_Sigma = sigma; }

  //
  void SetClassesToSmooth( std::vector<int> classesToSmooth )
    { m_classesToSmooth = classesToSmooth; } 
    
  //
  float  GetSigma() const
    { return m_Sigma; }  
    
  // 
  AtlasMeshCollection::Pointer  GetSmoothedMeshCollection();


protected :
  // Constructor
  AtlasMeshSmoother();
  
  // Destructor
  virtual ~AtlasMeshSmoother();
  
  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;


    
private :
  AtlasMeshSmoother(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  float  m_Sigma;
  std::vector<int> m_classesToSmooth;

};



} // end namespace kvl


#endif
