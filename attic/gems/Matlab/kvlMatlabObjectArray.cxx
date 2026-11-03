#include "kvlMatlabObjectArray.h"


namespace kvl
{

MatlabObjectArray::Pointer MatlabObjectArray::m_Instance = 0;



//
// Return the single instance of the MatlabObjectArray
//
MatlabObjectArray::Pointer
MatlabObjectArray
::New()
{
  
   if ( !MatlabObjectArray::m_Instance )
    {
    MatlabObjectArray::m_Instance = new MatlabObjectArray;

    // Remove extra reference from construction.
    MatlabObjectArray::m_Instance->UnRegister();
    }

  return MatlabObjectArray::m_Instance;
 
}


//
// 
//
int 
MatlabObjectArray
::AddObject( itk::Object* object )
{
  // Look for first avaible place to put this in
  int number = 0;
  while ( m_Array.find( number ) != m_Array.end() )
    { 
    number++;
    }  
  
  //std::cout << "     found available number: " << number << std::endl;
  m_Array[ number ] = object;
  
  return number;
}

  
//
// Remove an image, and evoke an event
//
void 
MatlabObjectArray
::RemoveObject( int number )
{
  
  if ( m_Array.find( number ) == m_Array.end() )
    {
    itkExceptionMacro( << "Trying to remove non-existing object (" << number << ")" );
    }

  m_Array.erase( number );

}

  
//
// 
//
itk::Object*
MatlabObjectArray
::GetObject( int number )
{

  std::map< int, itk::Object::Pointer >::const_iterator  it = m_Array.find( number );
  if ( it == m_Array.end() )
    {
    itkExceptionMacro( << "Trying to access non-existing object (" << number << ")" );
    }

  return it->second;

 
}



//
//
//
void
MatlabObjectArray
::SetObject( int number, itk::Object* object )
{

  m_Array[ number ] = object;

}




} // end namespace kvl


