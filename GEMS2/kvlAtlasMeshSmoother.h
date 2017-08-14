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
    { 
    this->SetSigmas( sigma, sigma, sigma ); 
    }

  // 
  void  SetSigmas( float sigma0, float sigma1, float sigma2 )
    { 
    m_Sigma0 = sigma0; 
    m_Sigma1 = sigma1; 
    m_Sigma2 = sigma2; 
    }

  //
  void  SetClassesToSmooth( const std::vector<int>&  classesToSmooth )
    { m_ClassesToSmooth = classesToSmooth; } 
    
    
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
  float  m_Sigma0;
  float  m_Sigma1;
  float  m_Sigma2;
  std::vector<int> m_ClassesToSmooth;

};



} // end namespace kvl


#endif
