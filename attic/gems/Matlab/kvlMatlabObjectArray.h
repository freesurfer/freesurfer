#ifndef __kvlMatlabObjectArray_h
#define __kvlMatlabObjectArray_h

#include "itkObject.h"


namespace kvl
{


/**
  * Object to hold an array of ITK objects in Matlab.
  *
  * This class is a so-called Singleton, i.e. it is intelligent enough to 
  * ensure that only one single instance of it can be created.
  *
  */


class MatlabObjectArray : public itk::Object
{
public:
  /** Smart pointer typedef support. */
  typedef MatlabObjectArray         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MatlabObjectArray, itk::Object);

  /** Intercept calls to create new instances to ensure that we have a Singleton */
  static Pointer New();
  
  /** Return the singleton instance. */
  static Pointer GetInstance()
    { return New(); }

  //
  int AddObject( itk::Object* object );

  //
  int AddObject( const itk::Object* object )
    {
    return this->AddObject( const_cast< itk::Object* >( object ) );
    }
  
  //
  void RemoveObject( int number );

  //
  itk::Object* GetObject( int number );

  //
  void SetObject( int number, itk::Object* object );

  //
  void Clear()
    { m_Array.clear(); }
  
protected:
  MatlabObjectArray() {};
  virtual ~MatlabObjectArray() {};


  MatlabObjectArray(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

  std::map< int, itk::Object::Pointer >    m_Array;
  static Pointer  m_Instance;

};

} // end namespace kvl

#endif


