#ifndef _SBL_UNIT_TEST_H_
#define _SBL_UNIT_TEST_H_
namespace sbl {


/*! \file UnitTest.h
    \brief The UnitTest module provides a very simple mechanism for creating 
    unit tests.  The "unittest" command can be run to execute all currently registered
    unit tests.
*/


// register commands, etc. defined in this module
void initUnitTest();


/// register a unit test; the test name is obtained from the callback function name
#define registerUnitTest( callback ) registerUnitTestInternal( #callback, callback, __FILE__ )
void registerUnitTestInternal( const char *name, bool (*callback)(), const char *fileName );


/// assert that a condition is true in a unit test
#define unitAssert( condition ) if (!(condition)) { disp( 3, "unit test failed: %s", #condition ); return false; }


} // end namespace sbl
#endif // _SBL_UNIT_TEST_H_

