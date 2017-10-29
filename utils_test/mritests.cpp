#include <boost/test/unit_test.hpp>

#include "error.h"
#include "mri.h"

BOOST_AUTO_TEST_SUITE( MRIstruct )

BOOST_AUTO_TEST_SUITE( Runtime )

BOOST_AUTO_TEST_CASE( SimpleAllocate )
{
  MRI* target = NULL;

  const int size = 128;
  const int nFrames = 1;

  target = MRIallocSequence( size, size, size, MRI_UCHAR, nFrames );

  BOOST_CHECK( target != NULL );
  
  // Deallocate
  BOOST_CHECK( NO_ERROR == MRIfree( &target ) );
  BOOST_CHECK( target == NULL );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
