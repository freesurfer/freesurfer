#include <boost/test/unit_test.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/mpl/list.hpp>

#include "randomsupply.hpp"
#include "testiosupport.hpp"

#include "dimensioncuda.hpp"

BOOST_AUTO_TEST_SUITE( DimensionCuda )

typedef boost::mpl::list<unsigned char,unsigned short,unsigned int,unsigned long> TestLengthTypes;

RandomSupply gPRNG;

BOOST_AUTO_TEST_CASE_TEMPLATE( testSetAndGet1d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<1,LengthType> testObject;
    
    const LengthType anyLength = gPRNG.GetInteger<LengthType>(12);
  
    testObject[0] = anyLength;

    BOOST_CHECK_EQUAL( testObject[0], anyLength );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testEqualityOperators1d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<1,LengthType> d1, d2, d3;

    const LengthType aLength = gPRNG.GetInteger<LengthType>(32);
    const LengthType bLength = aLength + 1;
  
    d1[0] = aLength;
    d2[0] = aLength;
    d3[0] = bLength;
    
    BOOST_TEST( d1 == d2 );
    BOOST_TEST( d1 != d3 );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testOutOfRange1d, LengthType, TestLengthTypes )
{
  kvl::cuda::Dimension<1,LengthType> testObject;

  BOOST_CHECK_THROW( testObject[10], std::range_error ); 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testSetAndGet2d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<2,LengthType> d;

    const LengthType l0 = gPRNG.GetInteger<LengthType>(34);
    const LengthType l1 = l1 + 2;
    
    d[0] = l0;
    d[1] = l1;
  
    BOOST_CHECK_EQUAL( l0,  d[0] );
    BOOST_CHECK_EQUAL( l1,  d[1] );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testGetLinearIndexOutOfRange, LengthType, TestLengthTypes )
{
  const LengthType n0 = 4;
  const LengthType n1 = 5;
  const LengthType n2 = 6;
  
  const LengthType zero = 0;
  
  kvl::cuda::Dimension<3,LengthType> d;

  d[0] = n0;
  d[1] = n1;
  d[2] = n2;

  BOOST_CHECK_THROW( d.GetLinearIndex((LengthType)(n0+1),zero,zero), std::range_error );
  BOOST_CHECK_THROW( d.GetLinearIndex(zero,(LengthType)(n1+1),zero), std::range_error );
  BOOST_CHECK_THROW( d.GetLinearIndex(zero,zero,(LengthType)(n2+1)), std::range_error );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testGetLinearIndex1d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    const LengthType n0 = 1+gPRNG.GetInteger<LengthType>(101);
    const LengthType i0 = gPRNG.GetInteger<LengthType>(n0);

    kvl::cuda::Dimension<1,LengthType> d;

    d[0] = n0;
    
    BOOST_CHECK_EQUAL( static_cast<size_t>(i0), d.GetLinearIndex(i0) );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testGetLinearIndex2d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    const LengthType n0 = 1+gPRNG.GetInteger<LengthType>(101);
    const LengthType n1 = 1+gPRNG.GetInteger<LengthType>(82);
    const LengthType i0 = gPRNG.GetInteger<LengthType>(n0);
    const LengthType i1 = gPRNG.GetInteger<LengthType>(1);
    
    kvl::cuda::Dimension<2,LengthType> d;
    
    d[0] = n0;
    d[1] = n1;
    
    size_t expected = (static_cast<size_t>(i0)*n1) + i1;
    
    BOOST_CHECK_EQUAL( expected, d.GetLinearIndex(i0,i1) );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testGetLinearIndex3d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    const LengthType n0 = 1+gPRNG.GetInteger<LengthType>(101);
    const LengthType n1 = 1+gPRNG.GetInteger<LengthType>(105);
    const LengthType n2 = 1+gPRNG.GetInteger<LengthType>(107);
    const LengthType i0 = gPRNG.GetInteger<LengthType>(n0);
    const LengthType i1 = gPRNG.GetInteger<LengthType>(n1);
    const LengthType i2 = gPRNG.GetInteger<LengthType>(n2);
    
    kvl::cuda::Dimension<3,LengthType> d;
    
    d[0] = n0;
    d[1] = n1;
    d[2] = n2;
  
    size_t expected =  i2 + (n2 * (i1 + (n1*static_cast<size_t>(i0))));
    
    BOOST_CHECK_EQUAL( expected, d.GetLinearIndex(i0,i1,i2) );
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testPointInRangePositive, LengthType, TestLengthTypes )
{
  const LengthType n0 = 4;
  const LengthType n1 = 13;
  const LengthType n2 = 23;

  const LengthType zero = 0;
  
  kvl::cuda::Dimension<3,LengthType> d;
  d[0] = n0;
  d[1] = n1;
  d[2] = n2;
  
  BOOST_TEST( d.PointInRange(zero,zero,zero) );
  BOOST_TEST( d.PointInRange((LengthType)(n0-1),(LengthType)(n1-1),(LengthType)(n2-1)) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testPointInRangeNegative, LengthType, TestLengthTypes )
{
  const LengthType n0 = 4;
  const LengthType n1 = 7;
  const LengthType n2 = 31;

  const LengthType inRange = 1;
  
  
  kvl::cuda::Dimension<3,LengthType> d;
  d[0] = n0;
  d[1] = n1;
  d[2] = n2;

  BOOST_CHECK_EQUAL( false, d.PointInRange(n0,inRange,inRange) );
  BOOST_CHECK_EQUAL( false, d.PointInRange(inRange,n1,inRange) );
  BOOST_CHECK_EQUAL( false, d.PointInRange(inRange,inRange,n2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testElementCount, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    kvl::cuda::Dimension<1,LengthType> d1;
    kvl::cuda::Dimension<3,LengthType> d3;
    
    const LengthType n0 = gPRNG.GetInteger<LengthType>(std::numeric_limits<LengthType>::max());
    const LengthType n1 = gPRNG.GetInteger<LengthType>(std::numeric_limits<LengthType>::max());
    const LengthType n2 = gPRNG.GetInteger<LengthType>(std::numeric_limits<LengthType>::max());
    
    const size_t total1 = n0;
    const size_t total3 = static_cast<size_t>(n0) *
      static_cast<size_t>(n1) *
      static_cast<size_t>(n2);
  
    d1[0] = d3[0] = n0;
    d3[1] = n1;
    d3[2] = n2;
    
    BOOST_CHECK_EQUAL( total1, d1.ElementCount() );
    BOOST_CHECK_EQUAL( total3, d3.ElementCount() );
  }
}

template<typename LengthType,int nDims>
void checkLinearIndexToLocation() {
  kvl::cuda::Dimension<nDims,LengthType> d;
  
  LengthType n[nDims];

  long double exponent = 1.0 / nDims;
  long double base = std::numeric_limits<size_t>::max();
  long double maxLength = pow(base,exponent);
  
  LengthType maxVal = std::numeric_limits<LengthType>::max()-1;
  if( maxLength < maxVal ) {
    maxVal = maxLength;
  }
  
  for( unsigned char i=0; i<nDims; i++ ) {
    n[i] = 1+gPRNG.GetInteger<LengthType>(maxVal);
    d[i] = n[i];
  }
  
  // Basic checks
  LengthType result[nDims];
  size_t elementCount = 1;
  for( unsigned char i=0; i<nDims; i++ ) {
    elementCount *= n[i];
  }
  
  // Check in range
  BOOST_CHECK_THROW( d.LinearIndexToLocation(elementCount, result), std::range_error );

  // Check first element
  d.LinearIndexToLocation(0, result);
  for( unsigned char i=0; i<nDims; i++ ) {
    BOOST_CHECK_EQUAL( result[i], 0 );
  }

  // Check last element
  d.LinearIndexToLocation(elementCount-1, result);
  for( unsigned char i=0; i<nDims; i++ ) {
    BOOST_CHECK_EQUAL( result[i], n[i]-1 );
  }

  // Check a random element
  size_t location = gPRNG.GetInteger<size_t>(elementCount);

  LengthType loc[nDims];
  d.LinearIndexToLocation(location, loc);
  BOOST_CHECK_EQUAL( location, d.GetLinearIndexFromArray(loc) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testLinearIndexToLocation1d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    checkLinearIndexToLocation<LengthType,1>();
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testLinearIndexToLocation2d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    checkLinearIndexToLocation<LengthType,2>();
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testLinearIndexToLocation3d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    checkLinearIndexToLocation<LengthType,3>();
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testLinearIndexToLocation4d, LengthType, TestLengthTypes )
{
  BOOST_TEST_CONTEXT("Seed: " << gPRNG.Reseed()) {
    checkLinearIndexToLocation<LengthType,4>();
  }
}


BOOST_AUTO_TEST_SUITE_END();
