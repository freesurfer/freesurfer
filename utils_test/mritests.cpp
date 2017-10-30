#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include "error.h"
#include "mri.h"

// #############################################################

const unsigned char mriTypes[] = { MRI_UCHAR, MRI_SHORT, MRI_INT, MRI_FLOAT };

// #############################################################

BOOST_AUTO_TEST_SUITE( MRIstruct )

BOOST_AUTO_TEST_SUITE( Runtime )

BOOST_AUTO_TEST_CASE( SimpleAllocate )
{
  MRI* target = NULL;

  const int size = 128;
  const int nFrames = 1;

  target = MRIallocSequence( size, size+1, size+2, MRI_UCHAR, nFrames );

  BOOST_REQUIRE( target != NULL );
  BOOST_CHECK_EQUAL( target->width, size );
  BOOST_CHECK_EQUAL( target->height, size+1 );
  BOOST_CHECK_EQUAL( target->depth, size+2 );
  BOOST_CHECK_EQUAL( target->nframes, nFrames );
  BOOST_CHECK_EQUAL( target->type, MRI_UCHAR );
  
  // Deallocate
  BOOST_CHECK( NO_ERROR == MRIfree( &target ) );
  BOOST_CHECK( target == NULL );
}

BOOST_DATA_TEST_CASE( ValueFill, boost::unit_test::data::make( mriTypes ), mriType )
{
  MRI* target = NULL;

  const int size = 4;
  const int nFrames = 1;

  const float value = 10.1;

  // Allocate the target structure
  target = MRIallocSequence( size, size, size, mriType, nFrames );
  BOOST_REQUIRE( target != NULL );

  // Call fill routine
  BOOST_REQUIRE_EQUAL( MRIvalueFill( target, value ), 0 );

  // Check the results
  for( int iy=0; iy<size; iy++ ) {
    for( int iz=0; iz<size; iz++ ) {
      for( int ix=0; ix<size; ix++ ) {
	switch (target->type) {
	case MRI_UCHAR:
	  BOOST_CHECK_EQUAL( MRIvox(target,ix,iy,iz), (unsigned char)(nint(value)) );
	  break;
	case MRI_SHORT:
	  BOOST_CHECK_EQUAL( MRISvox(target,ix,iy,iz), (short)(nint(value)) );
	  break;
	case MRI_INT:
	  BOOST_CHECK_EQUAL( MRIIvox(target,ix,iy,iz), (int)(nint(value)) );
	  break;
	case MRI_FLOAT:
	  BOOST_CHECK_EQUAL( MRIFvox(target,ix,iy,iz), value);
	  break;
	default:
	  BOOST_FAIL( "Unsupported target->type" );
	}
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END() // Runtime

BOOST_AUTO_TEST_SUITE_END() // MRIstruct
