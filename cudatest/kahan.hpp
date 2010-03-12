/*! \file
  Kahan summation object
*/


#ifndef KAHAN_CLASS_H
#define KAHAN_CLASS_H

template<typename T>
class KahanSum {
public:
  //! Default constructor
  KahanSum( const T initVal=0 ) : sum(initVal),
				  correction(0) {};
  
  //! Accessor for the sum
  T GetSum( void ) const {
    return( sum );
  }

  //! Assignment of a native type
  inline KahanSum& operator=( const T val ) {
    this->sum = val;
    this->correction = 0;

    return( *this );
  }


  //! Self addition of a native type
  inline KahanSum& operator+=( const T val ) {
    T y = val - this->correction;
    T t = this->sum + y;
    this->correction = ( t - this->sum ) - y;
    this->sum = t;

    return( *this );
  }

  //! Addition of a native type
  inline const KahanSum operator+( const T val ) const {
    KahanSum lhs = *this;
    lhs += val;
    return( lhs );
  }

private:
  //! Internal sum
  volatile T sum;
  //! Internal correction
  volatile T correction;
 
};


#endif
