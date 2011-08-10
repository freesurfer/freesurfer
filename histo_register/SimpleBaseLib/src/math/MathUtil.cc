// Licensed under MIT license; see license.txt.

#include <sbl/math/MathUtil.h>
#include <sbl/core/Display.h>
#include <stdlib.h>
namespace sbl {


/// returns a random float between 0.0 and 1.0;
/// note: this isn't a very good random number
float randomFloat() {
    double r = (double) randomInt( 0, 10000 );
    return (float) (r / 10000.0);
}


/// returns a random float between the specified bounds;
/// note: this isn't a very good random number
float randomFloat( float min, float max ) {
    return min + randomFloat() * (max - min);
}


/// returns a random integer within specified bounds;
/// note: this isn't a very good random number
int randomInt( int min, int max ) {
    assertDebug( max > min );
    return min + rand() % (max - min + 1);
}


/// set the random number generator seed;
/// note: this isn't a very good random number
void randomSeed( int seed ) {
    srand( seed );
}


} // end namespace sbl

