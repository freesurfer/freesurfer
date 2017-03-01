#define BOOST_TEST_MODULE GEMS2 Tests
#include <boost/test/included/unit_test.hpp>

// -------------------------

#ifdef CUDA_FOUND
#include "cudaglobalfixture.hpp"
BOOST_GLOBAL_FIXTURE(CUDAGlobalFixture);
#endif
