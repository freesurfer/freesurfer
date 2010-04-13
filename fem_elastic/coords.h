
#ifndef H_COORDS_H
#define H_COORDS_H

#include <assert.h>
#include <cmath>
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stack>
#include <vector>

typedef enum
{
  cValid,
  cInvalid,
  cOutOfBounds,
  cOther
}
CoordsStatusType;

//
// generic Coords class
// although the first type is declared as a parameter,
// this class is mainly intended for usage with scalar types,
// such as int, float or double.
//
// The implementation should be checked before other instantiations are used!
template<class T, int n>
class TCoords
{
public:
  TCoords();
  TCoords(T _x);
  TCoords(const TCoords& _x);

  void set(T _x);

  TCoords operator=(const TCoords& _x);

  double norm() const;
  void   div(const TCoords& _d);

  T   operator()(int i) const
  {
    assert(i>-1&&i<n);
    return m_pdata[i];
  }
  T&  operator()(int i)
  {
    assert(i>-1&&i<n);
    return m_pdata[i];
  }

  void print(std::ostream& os) const;

  CoordsStatusType& status()
  {
    return m_status;
  }
  const CoordsStatusType& status() const
  {
    return m_status;
  }
  void invalidate()
  {
    m_status = cInvalid;
  }
  void validate()
  {
    m_status = cValid;
  }
  bool isValid() const
  {
    return m_status==cValid;
  }
private:
  T  m_pdata[n];
  CoordsStatusType m_status;

  void clone(const TCoords& _x);
};

// typedefs

typedef TCoords<double, 3> Coords3d;

template<class T, int n>
double dist(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> operator-(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> operator-(const TCoords<T,n>& a, const T& b);
template<class T, int n>
TCoords<T,n> operator*(const TCoords<T,n>& a, const T& b);
template<class T, int n>
TCoords<T,n> operator/(const TCoords<T,n>& a, const T& b);
template<class T, int n>
TCoords<T,n> operator+(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> operator+(const TCoords<T,n>& a, const T& b);
template<class T, int n>
TCoords<T,n>& operator+=(TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n>& operator-=(TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> operator*(const T& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> min(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> max(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
bool operator==(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
TCoords<T,n> operator!=(const TCoords<T,n>& a, const TCoords<T,n>& b);
template<class T, int n>
bool operator<=(const TCoords<T,n>& a, const TCoords<T,n>& b);

template<class T, int n>
std::ostream& operator<<( std::ostream& os, const TCoords<T,n>& c);

template<class T, int n>
std::istream& operator>>( std::istream& is, TCoords<T,n>& c);

//----------------------------------
//
// special implementations
//

template<int n>
double dot(const TCoords<double,n>& a,
           const TCoords<double,n>& b); // dot product
template<int n>
TCoords<double,n> operator*(double val, const TCoords<double,n>& c);

template<int n>
std::stack<TCoords<int,n> >*
neighbors(const TCoords<int,n>& c);

//------------------------------------

// lexicographic order -> mainly for the int instanciation
template<class T, int n>
struct cless : std::binary_function< TCoords<T,n>, TCoords<T,n>, bool>
{
  bool operator()(const TCoords<T,n>& a, const TCoords<T,n>& b) const
  {
    bool bEq = false;
    for (int i=0; i<n; ++i)
      if ( a(i)>b(i) ) return false;
      else if ( a(i) == b(i) ) bEq = true;
      else if ( a(i) < b(i) ) return true;

    return (!bEq); // less is a strict relationship -> if here,
    // then all the components are equal.
  }
};

//-------------------------------------
//-
// Definitions
//
//-------------------------------------


template<class T, int n>
TCoords<T,n>::TCoords()
    : m_status(cValid)
{
  assert(n>0);
}

template<class T, int n>
TCoords<T,n>::TCoords(T _x)
    : m_status(cValid)
{
  assert(n>0);
  std::fill_n( m_pdata, n, _x);
}

template<class T, int n>
void
TCoords<T,n>::set(T _x)
{
  std::fill_n( m_pdata, n, _x);
}

template<class T, int n>
void
TCoords<T,n>::clone(const TCoords<T,n>& _x)
{
  for (int i=0; i<n; ++i)
    m_pdata[i] = _x.m_pdata[i];
  this->status() = _x.status();
}

template<class T, int n>
TCoords<T,n>::TCoords(const TCoords<T,n>& _x)
{
  assert(n>0);
  clone(_x);
}

template<class T, int n>
TCoords<T,n>
TCoords<T,n>::operator=(const TCoords<T,n>& _x)
{
  clone(_x);
  return *this;
}

template<class T, int n>
void
TCoords<T,n>::div(const TCoords<T,n>& _d)
{
  T* p = m_pdata;
  const T* pd = _d.m_pdata;
  for ( int i=0; i<n; ++i, ++p, ++pd)
  {
    *p /= *pd;
  }
}

template<class T, int n>
double
TCoords<T,n>::norm() const
{
  double dret = 0.0;
  double dbuf;
  int i=0;
  for ( const T* p = m_pdata; i<n; ++i, ++p)
  {
    dbuf = double(*p);
    dret += dbuf*dbuf;
  }
  return std::sqrt(dret);
}

template<class T, int n>
void
TCoords<T,n>::print(std::ostream& os) const
{
  //os << "( ";
  os << m_pdata[0];
  for (int i=1; i<n; ++i)
    os  << " " << m_pdata[i];
  //os << " )";
}

template<class T, int n>
double
dist(const TCoords<T,n>& a,
     const TCoords<T,n>& b)
{
  return (a-b).norm();
}

template<class T, int n>
TCoords<T,n>
operator-(const TCoords<T,n>& a,
          const TCoords<T,n>& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; i++)
    ret(i) = a(i) - b(i);

  return ret;
}

template<class T, int n>
TCoords<T,n>
operator-(const TCoords<T,n>& a,
          const T& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = a(i)-b;

  return ret;
}

template<class T, int n>
TCoords<T,n>
operator*(const TCoords<T,n>& a,
          const T& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = a(i)*b;

  return ret;
}

template<class T, int n>
TCoords<T,n>
operator/(const TCoords<T,n>& a,
          const T& b)
{
  assert(std::abs(b)>1.0e-10);

  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = a(i)/b;

  return ret;
}

template<class T, int n>
TCoords<T,n>
operator+(const TCoords<T,n>& a,
          const TCoords<T,n>& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; i++)
    ret(i) = a(i) + b(i);

  return ret;
}

template<class T, int n>
TCoords<T,n>
operator+(const TCoords<T,n>& a,
          const T& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = a(i) + b;

  return ret;
}

template<class T, int n>
TCoords<T,n>
min(const TCoords<T,n>& a,
    const TCoords<T,n>& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = std::min( a(i), b(i));

  return ret;
}

template<class T, int n>
TCoords<T,n>
max(const TCoords<T,n>& a,
    const TCoords<T,n>& b)
{
  TCoords<T,n> ret;
  for (int i=0; i<n; ++i)
    ret(i) = std::max( a(i), b(i));

  return ret;
}

template<class T, int n>
bool
operator==(const TCoords<T,n>& a,
           const TCoords<T,n>& b)
{
  for (int i=0; i<n; ++i)
    if ( a(i)!=b(i) )
      return false;

  return true;
}

template<class T, int n>
bool
operator<=(const TCoords<T,n>& a,
           const TCoords<T,n>& b)
{
  for (int i=0; i<n; ++i)
    if ( a(i) > b(i) ) return false;

  return true;
}

template<class T, int n>
bool
operator!=(const TCoords<T,n>& a,
           const TCoords<T,n>& b)
{
  return !(a==b);
}

template<class T, int n>
std::ostream& operator<<( std::ostream& os, const TCoords<T,n>& c)
{
  c.print(os);
  return os;
}

template<class T, int n>
std::istream& operator>>( std::istream& is, TCoords<T,n>& coords)
{
#if 1
  for (int i=0; i<n; ++i)
    is >> coords(i);
#else
  char c = 0;
  is >> c;
  if (c=='(')
  {
    for ( int i=0; i<n-1; ++i)
    {
      is >> coords(i) >> c;
      if ( c!=',')
      {
        std::cerr << "error reading TCoords\n";
        exit(1);
      }
    }
    is >> coords(n-1) >> c;
    if ( c!=')' )
    {
      std::cerr << " error reading TCoords -> did not find terminal )\n";
      exit(1);
    }
  }
#endif
  return is;
}

//------------------------

template<int n>
double
dot(const TCoords<double,n>& a,
    const TCoords<double,n>& b)
{
  double dret = 0.0;
  for (int i=0; i<n; ++i)
    dret += a(i)*b(i);

  return dret;
}

template<int n>
TCoords<double,n>
operator*(double val,
          const TCoords<double,n>& c)
{
  TCoords<double, n> _c(c);
  for (int i=0; i<n; ++i)
    _c(i) *= val;

  return _c;
}


template<int n>
std::stack<TCoords<int,n> >*
neighbors(const TCoords<int,n>& c)
{
  std::stack<TCoords<int,n> >* pstack =
    new std::stack<TCoords<int,n> >;

  // for each coordinate, 3 possibilities -> -1, 0, 1
  int p3[n];
  int count = 1;
  for (int i=0; i<n; ++i)
  {
    p3[i] = count;
    count*=3;
  }
  TCoords<int,n> tic_buf;
  for ( int index=0; index<count; ++index)
  {
    for (int alpha=0; alpha<n; ++alpha)
      tic_buf(alpha) = (index/p3[alpha])%3 -1;
    pstack->push(tic_buf+c);
  }

  return pstack;
}

template<class T, int n>
TCoords<T,n>& operator+=(TCoords<T,n>& a, const TCoords<T,n>& b)
{
  for (int i=0; i<n; ++i)
    a(i) += b(i);
  return a;
}

template<class T, int n>
TCoords<T,n>& operator-=(TCoords<T,n>& a, const TCoords<T,n>& b)
{
  for (int i=0; i<n; ++i)
    a(i) -= b(i);
  return a;
}

template<class T, int n>
TCoords<T,n> operator*(const T& a, const TCoords<T,n>& b)
{
  TCoords<T,n> retVal;
  for (int i=0; i<n; ++i)
    retVal(i) = a* b(i);

  return retVal;
}

#endif // H_COORDS_H
