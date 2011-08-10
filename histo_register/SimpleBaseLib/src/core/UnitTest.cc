// Licensed under MIT license; see license.txt.

#include <sbl/core/UnitTest.h>
#include <sbl/core/Command.h>
#include <sbl/system/Timer.h>
namespace sbl {


//-------------------------------------------
// UNIT TEST CLASS
//-------------------------------------------


/// The UnitTest class stores a single unit test callback.
class UnitTest {
public:

    // basic constructor
    UnitTest( const String &name, bool (*callback)(), const String &fileName );

    /// run the unit test
    bool run();

private:

    // the callback and associated data
    bool (*m_callback)();
    String m_name;
    String m_fileName;
};


// basic constructor
UnitTest::UnitTest( const String &name, bool (*callback)(), const String &fileName ) { 
    m_name = name; 
    m_callback = callback; 
    m_fileName = fileName; 
}


/// run the unit test
bool UnitTest::run() {
    Timer timer;
    timer.start();
    bool result = m_callback(); 
    timer.stop();
    disp( 2, "%s (%s) time: %1.1fms, result: %s", m_name.c_str(), m_fileName.c_str(), timer.timeSum(), result ? "passed" : "FAILED" );
    return result;
}


//-------------------------------------------
// UNIT TEST MANAGEMENT
//-------------------------------------------


// the set of currently registered unit tests
Array<UnitTest> &unitTestSet() {
    static Array<UnitTest> s_unitTests;
    return s_unitTests;
}


// register a unit test
void registerUnitTestInternal( const char *name, bool (*callback)(), const char *fileName ) {
    unitTestSet().append( new UnitTest( name, callback, fileName ));
}


// run all unit tests that are currently registered
void runUnitTests( Config &conf ) {
    disp( 1, "unit tests: %d", unitTestSet().count() );
    int passCount = 0;
    for (int i = 0; i < unitTestSet().count(); i++) {
        if (unitTestSet()[ i ].run())
            passCount++;
    }
    int failCount = unitTestSet().count() - passCount;
    disp( 1, "%d passed, %d failed", passCount, failCount );
}


// register commands, etc. defined in this module
void initUnitTest() {
    registerCommand( "unittest", runUnitTests );
}


} // end namespace sbl

