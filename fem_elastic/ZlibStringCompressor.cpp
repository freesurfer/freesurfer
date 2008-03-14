#include <zlib.h>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
//#include "mrs/message.h"

#include "ZlibStringCompressor.h"
//#include "ZlibException.h"
using std::string;

ZlibStringCompressor::ZlibStringCompressor(unsigned long bufferSize)
    : m_bufferAllocationMultiplier(10), m_bufferSize(0), m_buffer(0), m_debugLevel(0)
{
  if (bufferSize) setBufferSize(bufferSize);
}

unsigned long ZlibStringCompressor::getBufferSize()
{
  return m_bufferSize;
}

Bytef* ZlibStringCompressor::setBufferSize(const unsigned long size)
{
  delete m_buffer;
  if (getDebugLevel()>0) std::cout << "Deleted buffer" << std::endl;
  m_bufferSize=0;
  try
  {
    m_buffer = new Bytef[size];
  }
  catch (std::bad_alloc&)
  {
    std::cerr << " requested buffer size = " << size << std::endl;
    throw " ZlibStringCompressor::setBufferSize - Caught std::bad_alloc";
    //ZlibException("Caught std::bad_alloc", __FILE__, __LINE__);
  }
  if (getDebugLevel()>0) std::cout << "Set buffer size to " << size  << std::endl;
  m_bufferSize=size;
  return m_buffer;
}

void ZlibStringCompressor::checkZlibResult(const int result, const char* file, const int line)
{
  switch (result)
  {
  case  (Z_OK) :
  { /* "zlib returned ok" */ break;
  }
  case (Z_MEM_ERROR) :
  {
    //throw ZlibException("zlib reports memory error", file, line);
    throw " ZlibStringCompressor::checkZlibResult - zlib reports memory error";
  }
  case (Z_BUF_ERROR) :
  {
    //throw ZlibException("zlib buffer insufficient", file, line);
    throw " ZlibStringCompressor::checkZlibResult - zlib buffer insufficient";
  }
  case (Z_STREAM_ERROR) :
  {
    // throw ZlibException("zlib reports stream error", file, line);
    throw "ZlibStringCompressor::checkZlibResult - zlib reports stream error";
  }
  case (Z_DATA_ERROR) :
  {
    //throw ZlibException("zlib reports data corrupted error", file,  line);
    throw " ZlibStringCompressor::checkZlibResult - zlib reports data corrupter error";
  }
  default :
  {
    std::ostringstream oss;
    oss << "ZlibStringCompressor::checkZlibResult  - zlib returned the unknown value of `" << result << "'";
    //throw ZlibException(oss.str(), file, line);
    throw oss.str();
  }
  }
}

const string ZlibStringCompressor::compress(const string& s_in, int level)
{
  // find minimum buffer size for compressed string
  unsigned long required_length=findRequiredBufferSize(findSizeInBytes(s_in));
  // make sure buffer is at least that size - if its too small add a factor of 2 safety margin!
  if ( required_length>getBufferSize() )
  {
    setBufferSize( required_length*2 );
  }

  // reinterpret input string as an array of bytes.
  const char* input_char = s_in.c_str();
  const Bytef* input_bytes = reinterpret_cast<const Bytef*> ( input_char );

  // make sure level is in the range 1->9
  int checked_level = (level<1) ? 1 : (level>9) ? 9 : level;

  unsigned long output_length = getBufferSize();
  int result = compress2 (m_buffer, &output_length, input_bytes,
                          findSizeInBytes(s_in), checked_level);

  checkZlibResult(result, __FILE__, __LINE__);
  const char* result_char=reinterpret_cast<const char*>(m_buffer);
  unsigned long char_length=output_length/sizeof(char);
  const string s_out(result_char, char_length);
  return s_out;
}

const string ZlibStringCompressor::inflate(const string& s)
{
  long unsigned compressed_length=findSizeInBytes(s);
  if (compressed_length*3 > getBufferSize() )
  {
    setBufferSize(compressed_length*m_bufferAllocationMultiplier);
  }

  long unsigned inflated_length=0;
  const char* compressed_c_str=s.c_str();
  const Bytef* compressed_bytes=reinterpret_cast<const Bytef*>(compressed_c_str);

  int result = Z_OK;
  do
  {
    inflated_length=getBufferSize();
    result = uncompress (m_buffer, &inflated_length, compressed_bytes, findSizeInBytes(s));
    if (result==Z_BUF_ERROR) setBufferSize(getBufferSize()*2);
  }
  while (result==Z_BUF_ERROR);

  checkZlibResult(result,__FILE__,__LINE__);

  unsigned long inflated_char_length=inflated_length/sizeof(char);
  const string s_inflated((char*)m_buffer, inflated_char_length);

  return s_inflated;
}

