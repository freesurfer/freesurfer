#include "Array2.h"

template <class T>
Array2<T>::Array2 ( int izX, int izY ) {
  this->mzX = izX;
  this->mzY = izY;
  this->mArray = new T[this->mzX * this->mzY];
}

template <class T>
Array2<T>::Array2 ( int izX, int izY, T iInitValue ) {

  this->mzX = izX;
  this->mzY = izY;
  this->mArray = new T[this->mzX * this->mzY];

  for( int n = 0; n < this->mzX * this->mzY; n++){    
    this->mArray[n]= iInitValue;
  }
}

template <class T>
void 
Array2<T>::SetAll ( T const iValue ) {
  for( int n = 0; n < this->mzX * this->mzY; n++){    
    this->mArray[n]= iValue;
  }
}
