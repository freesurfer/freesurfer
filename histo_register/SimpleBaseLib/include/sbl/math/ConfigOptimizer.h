#ifndef _SBL_CONFIG_OPTIMIZER_H_
#define _SBL_CONFIG_OPTIMIZER_H_
#include <sbl/core/Config.h>
namespace sbl {


/// The ConfigOptimizerObjective class provides an interface for objective functions used by the runConfigOptimizer function.
class ConfigOptimizerObjective {
public:

    // avoid non-virtual distructor warnings
    virtual ~ConfigOptimizerObjective() {}

    /// evaluate a particular set of Config parameters
    virtual double eval( Config &conf, const String &paramUpdate ) = 0;
};


/// optimize an objective function by varying parameters in a Config object
void runConfigOptimizer( ConfigOptimizerObjective &objective, Config &conf, 
                         double factor, 
                         const String &paramFileName, 
                         const String &outputPath );


} // end namespace sbl
#endif // _SBL_CONFIG_OPTIMIZER_H_

