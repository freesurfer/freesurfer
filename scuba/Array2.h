#ifndef Array2_h
#define Array2_h

#include "Point2.h"

template <typename T>
class Array2 { 
 public:
  Array2 ( int const izX, int const izY );
  Array2 ( int const izX, int const izY, T const iInitValue );
  ~Array2 () { if( this->mArray ) delete[] this->mArray; };

  // Get/Set functions 
  T Get ( int const inX, int const inY ) const
    { return this->mArray[inX + inY * this->mzX]; };
  void Set ( int const inX, int const inY, T const iValue ) 
    { this->mArray[inX + inY * this->mzX] = iValue; };

  void SetAll ( T const iValue );

  // return the array element itself
  T& operator () ( int const inX, int const inY )
    { return this->mArray[inX + inY * this->mzX]; };

  // return a pointer to the array element
  T* Element ( int const inX, int const inY )
    { return( this->mArray + inX + inY * this->mzX ); };


  T Get ( Point2<int>& iLocation ) const
    { return this->mArray[iLocation.x() + iLocation.y() * this->mzX]; };
  void Set ( Point2<int>& iLocation, T const iValue ) 
    { this->mArray[iLocation.x() + iLocation.y() * this->mzX] = iValue; };

  T& operator () ( Point2<int>& iLocation )
    { return this->mArray[ iLocation.x() + iLocation.y() * this->mzX]; };
  T* Element ( Point2<int>& iLocation )
    { return( this->mArray + iLocation.x() + iLocation.y() * this->mzX ); };


 private:
  T* mArray;
  int mzX, mzY; 
};


#endif
