#ifndef ZLIBSTRINGCOMPRESSOR_H
#define ZLIBSTRINGCOMPRESSOR_H
#include <zlib.h>
#include <string>
using std::string;

class ZlibStringCompressor
{
public:
  ZlibStringCompressor(unsigned long bufferSize=0);
  const string compress(const string& toCompress, int level);
  const string inflate(const string& compressed);
  unsigned long getBufferSize();
  void setDebugLevel(int level);
  int getDebugLevel();
  Bytef* setBufferSize(unsigned long size);
  unsigned int m_bufferAllocationMultiplier;
private:
  static void checkZlibResult(const int result,
                              const char* file,
                              const int line);
  long unsigned m_bufferSize;
  /*
    pointer to start of buffer
  */
  Bytef* m_buffer;
  static unsigned long findSizeInBytes(const string& s);
  static unsigned long findRequiredBufferSize
  (const unsigned long uncompressedSize);
  int m_debugLevel;
};


inline unsigned long ZlibStringCompressor::findSizeInBytes(const string&s )
{
  return sizeof(char)*s.size();
}

// buffer size myst be at least original size + 0.1% + 12 bytes
inline unsigned long ZlibStringCompressor::findRequiredBufferSize
(long unsigned unbufferedSize)
{
  return unbufferedSize + unbufferedSize/1000 + 12;
}

inline void ZlibStringCompressor::setDebugLevel(int level)
{
  m_debugLevel=level;
}

inline int ZlibStringCompressor::getDebugLevel()
{
  return m_debugLevel;
}

#endif //  ZLIBSTRINGCOMPRESSOR_H
