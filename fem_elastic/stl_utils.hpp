/*

Gheorghe Postelnicu, 2007

Mainly intended to hold general usage functors

*/

#ifndef _h_stl_utils_hpp
#define _h_stl_utils_hpp

struct FunctorDeletePointer
{
  template<class T>
  T* operator()(T* p) const
  {
    delete p;
    return 0;
  }
};


#endif
