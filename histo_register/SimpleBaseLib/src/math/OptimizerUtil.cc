// Licensed under MIT license; see license.txt.

#include <sbl/math/OptimizerUtil.h>
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h> 
#include <sbl/math/Optimizer.h>
namespace sbl {


//-------------------------------------------
// OPTIMIZER TESTING
//-------------------------------------------


// http://en.wikipedia.org/wiki/Rosenbrock_function
// global min at (1, 1)
class Objective1 : public Objective {
    double eval( const VectorD &point ) {
        double x = point[ 0 ];
        double y = point[ 1 ];
        double a = 1 - x;
        double b = y - x * x;
        return a * a + 100 * b * b;
    }
};


// http://en.wikipedia.org/wiki/Himmelblau's_function
// four local min (x and y in [-6, 6])
class Objective2 : public Objective {
    double eval( const VectorD &point ) {
        double x = point[ 0 ];
        double y = point[ 1 ];
        double a = x * x + y - 11;
        double b = x + y * y - 7;
        return a * a + b * b;
    }
};


// test the optimizer on a 2D sample problem
void testOptimizer( const String &outputPrefix, Objective &objective, double xStart, double yStart, double min, double max ) {
    SimplexOptimizer optimizer( objective );
    optimizer.storeHistory();

    // set starting point
    VectorD start( 2 );
    start[ 0 ] = xStart;
    start[ 1 ] = yStart;
    optimizer.setStart( start );

    // set bounds
    VectorD lBound( 2 ), uBound( 2 );
    lBound.clear( min );
    uBound.clear( max );
    optimizer.setBounds( lBound, uBound );

    // run the optimization
    VectorD result = optimizer.run();

    // display diagnostics
    disp( 1, "result: %f, %f", result[ 0 ], result[ 1 ] );
    optimizer.dispHistoryStats( 2 );
    optimizer.plotBestHistory( addDataPath( outputPrefix + "_best.svg" ) );
    optimizer.plotEvalHistory( addDataPath( outputPrefix + "_eval.svg" ) );
    optimizer.plotObjectiveSlices( result, 100, addDataPath( outputPrefix ) );
    optimizer.plotObjectiveMap( result, 0, 1, 50, 50, 4, addDataPath( outputPrefix + ".png" ) );
}


// test the optimizer on some 2D sample problems
void testOptimizer( Config &conf ) {
    Objective1 obj1;
    testOptimizer( "obj1", obj1, 0.1, 0.1, 0, 2 );
    Objective2 obj2;
    testOptimizer( "obj2", obj2, 0, 0, -6, 6 );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initOptimizerUtil() {
#ifdef REGISTER_TEST_COMMANDS
    registerCommand( "opttest", testOptimizer );
#endif
}


} // end namespace sbl

