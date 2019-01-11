/*------------------------------------------------------------------------------*/
/** 
*  \file   fm_config.h
*  \brief  Main configuration file.
*  \author Gabriel Peyre
*  \date   2004
*/ 
/*------------------------------------------------------------------------------*/

#ifndef _MY_CONFIG_H_
#define _MY_CONFIG_H_

#define INLINE __inline

//-------------------------------------------------------------------------
/** \name debug & inline directive */
//-------------------------------------------------------------------------
//@{ 
#ifdef _DEBUG
  #ifndef GW_DEBUG
    #define GW_DEBUG
  #endif // GW_DEBUG
#endif // _DEBUG
// #undef GW_DEBUG
//@}

//-------------------------------------------------------------------------
/** \name numerical macros */
//-------------------------------------------------------------------------
//@{
#undef    MIN
#undef    MAX
#define   MIN(a, b)       ((a) < (b) ? (a) : (b))     //!<  Returns the min value between a and b
#define   MAX(a, b)       ((a) > (b) ? (a) : (b))     //!<  Returns the max value between a and b
#define   MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c))
#define   GW_MIN(a, b)       MIN(a,b)           //!<  Returns the min value between a and b
#undef    GW_MAX  // already defined by Windows.h ...
#define   GW_MAX(a, b)       MAX(a,b)           //!<  Returns the max value between a and b
#define   GW_MAXMAX(a,b,c)   MAXMAX(a,b,c)

#define GW_SCALE_01(x,rMin,rMax) ((x-rMin)/(rMax-rMin))

#define   GW_ABS(a)       ((a) > 0 ? (a) : -(a))      //!<  Returns the absolute value a
#define   GW_SIGN(a)       ((a) > 0 ? 1 : -1)       //!<  Returns the sign of a
#define   SQR(x)      ((x)*(x))           //!<  Returns x square
#define   CUBE(x)     ((x)*(x)*(x))         //!<  Returns x cube
#define   GW_SQR(x)   SQR(x)              //!<  Returns x square
#define   GW_CUBE(x)    CUBE(x)             //!<  Returns x cube

#define GW_CLAMP_01(x)  if( (x)<0 ) x=0; if( (x)>1 ) x=1
#define GW_CLAMP(x, a,b)  if( (x)<a ) x=a; if( (x)>b ) x=b

#define GW_SWAP(x,y) x^=y; y^=x; x^=y
#define GW_ORDER(x,y) if(x>y){ GW_SWAP(x,y); }
//@}

//-------------------------------------------------------------------------
/** \name generic macros */
//-------------------------------------------------------------------------
//@{
/** a random number in [0-1] */
#define GW_RAND ((double) (rand()%10000))/10000
/** a random number in [a,b] */
#define GW_RAND_RANGE(a,b) (a)+((b)-(a))*((GW_Float) (rand()%10000))/10000
/** delete a single pointer */
#define GW_DELETE(p) {if (p!=NULL) delete p; p=NULL;}
/** delete an array pointer */
#define GW_DELETEARRAY(p) {if (p!=NULL) delete [] p; p=NULL;}
//@}

//-------------------------------------------------------------------------
/** \name some constants */
//-------------------------------------------------------------------------
//@{
#define GW_True  true
#define GW_False false
/** to make aproximate computations (derivation, GW_Float comparaisons ...) */
#define GW_EPSILON 1e-9
/** very big number */
#define GW_INFINITE 1e9
//@}

//-------------------------------------------------------------------------
/** \name numerical  constants */
//-------------------------------------------------------------------------
//@{
/** pi */
#define GW_PI   3.1415926535897932384626433832795028841971693993751f
/** pi/2 */
#define GW_HALFPI 1.57079632679489661923f
/** 2*pi */
#define GW_TWOPI  6.28318530717958647692f
/** 1/pi */
#define GW_INVPI  0.31830988618379067154f
/** 180/pi */
#define GW_RADTODEG(x)  (x)*57.2957795130823208768f
/** pi/180 */
#define GW_DEGTORAD(x)  (x)*0.01745329251994329577f
/** e */
#define GW_EXP    2.71828182845904523536f
/** 1/log10(2) */
#define GW_ILOG2  3.32192809488736234787f
/** 1/3 */
#define GW_INV3   0.33333333333333333333f
/** 1/6 */
#define GW_INV6   0.16666666666666666666f
/** 1/9 */
#define GW_INV7   0.14285714285714285714f
/** 1/9 */
#define GW_INV9   0.11111111111111111111f
/** 1/255 */
#define GW_INV255 0.00392156862745098039f
/** sqrt(2) */
#define GW_SQRT2    1.41421356237f
//@}

//-------------------------------------------------------------------------
/** \name assertion macros */
//-------------------------------------------------------------------------
//@{
#ifdef GW_DEBUG
  #define GW_ASSERT(expr) _ASSERT(expr)
  #define GW_DEBUG_ONLY(expr) expr
#else
  #define GW_ASSERT(expr) // if(!(expr)) cerr << "Error in file " << __FILE__ << " line " << __LINE__ << "." << endl 
  #define GW_DEBUG_ONLY(expr)
#endif // GW_DEBUG
//@}

INLINE
int log2(int n)
{
  int x = 0;
  while (n > 1) 
  {
    x++;
    n /= 2;
  }
  return x;
}

INLINE
int cmp(double a, double b)
{ 
  if (a < b)
    return -1;
  if (a == b)
    return 0;
  return 1;
}

#endif // _MY_CONFIG_H_