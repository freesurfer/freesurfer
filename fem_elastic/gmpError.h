
#ifndef h_gmpError_h_
#define h_gmpError_h_

#include <string>

class gmpErr
{
public:
  gmpErr(std::string msg)
  {
    m_msg = msg;
  }
  std::string what() const
  {
    return m_msg;
  }
protected:
  std::string m_msg;
};

#endif
