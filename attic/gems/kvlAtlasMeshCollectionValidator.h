#ifndef __kvlAtlasMeshCollectionValidator_h
#define __kvlAtlasMeshCollectionValidator_h

#include "kvlAtlasMeshCollection.h"


namespace kvl
{


class AtlasMeshCollectionValidator: public itk::Object
{
public :
  
  /** Standard class typedefs */
  typedef AtlasMeshCollectionValidator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionValidator, itk::Object );
  
  //
  bool Validate( const AtlasMeshCollection* meshCollection );

protected :
  // Constructor
  AtlasMeshCollectionValidator();
  
  // Destructor
  virtual ~AtlasMeshCollectionValidator();
  
  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;


    
private :
  AtlasMeshCollectionValidator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



} // end namespace kvl


#endif
