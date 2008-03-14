#ifndef _h_tag_fio_h
#define _h_tag_fio_h

#include <fstream>
#include <string>
#include <sstream>

#include "streamIoMisc.h"

namespace ftags
{

class TagReader
{
protected:
  std::istream& m_istream;

public:
  char* m_data;
  long  m_len;
  int   m_tag;
  TagReader(std::istream& is) : m_istream(is),
      m_data(NULL), m_len(0), m_tag(-1)
  {}

  ~TagReader()
  {
    if (m_len) delete[] m_data;
  }

  bool Read()
  {
    m_tag = TRead<int>(m_istream);
    if ( m_istream.eof() ) return false;
    m_len = TRead<long>(m_istream);

    if (m_data) delete[] m_data;
    m_data = new char[m_len+1];
    m_istream.read(m_data, m_len);

    return true;
  }
};


std::string
CreateTag(int tag,
          std::string contents)
{
  std::ostringstream oss;
  oss.write( (const char*)(&tag),
             sizeof(int) );

  long len = contents.size();
  oss.write( (const char*)(&len),
             sizeof(long) );

  oss.write( contents.c_str(),
             len );
  return oss.str();
}


} // end namespace fio

#endif
