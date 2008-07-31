project(itkio)

SET (libitkio_a_SOURCES
  itkMGHImageIO.cxx
  itkMGHImageIOFactory.cxx 
 )
 
add_library (itkio ${libitkio_a_SOURCES})

