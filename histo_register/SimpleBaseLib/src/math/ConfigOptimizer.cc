// Licensed under MIT license; see license.txt.

#include <sbl/math/ConfigOptimizer.h>
#include <sbl/core/StringUtil.h>
#include <sbl/core/Command.h>
namespace sbl {


/// evaluate a particular set of Config parameters
double evalConfig( ConfigOptimizerObjective &objective, Config &conf, 
                   const String &paramUpdate, const String &outputPath ) {

    // compute the score
    double score = objective.eval( conf, paramUpdate );

    // add to log file
    File logFile( outputPath + "/confOpt.txt", FILE_APPEND, FILE_TEXT );
    if (logFile.openSuccess()) {
        logFile.writeF( "%s, %f\n", paramUpdate.c_str(), score );
    }
    return score;
}


/// optimize an objective function by varying parameters in a Config object
void runConfigOptimizer( ConfigOptimizerObjective &objective, Config &conf, 
                         double factor,
                         const String &paramFileName, 
                         const String &outputPath ) {

    // get baseline
    double baseline = evalConfig( objective, conf, "baseline", outputPath );

    // update each parameters
    double bestScore = baseline;
    String bestParamUpdate = "baseline";
    for (int entryIndex = 0; entryIndex < conf.entryCount(); entryIndex++) {

        // get info about this config entry
        const String &name = conf.entry( entryIndex ).name;
        float origValue = conf.entry( entryIndex ).value.toFloat();

        // try decreasing and increasing value
        for (int dir = -1; dir <= 1; dir += 2) {

            // evaluate new parameter value
            double newValue = dir > 0 ? origValue * factor : origValue / factor;
            String paramUpdate = sprintF( "%s=%f", name.c_str(), newValue );
            conf.update( paramUpdate );
            double score = evalConfig( objective, conf, paramUpdate, outputPath );

            // check for user cancel
            if (checkCommandEvents())
                break;

            // store best score
            if (score < bestScore) {
                bestScore = score;
                bestParamUpdate = paramUpdate;
            }
        }

        // check for user cancel
        if (checkCommandEvents())
            break;

        // store original value and return best
        conf.update( sprintF( "%s=%f", name.c_str(), origValue ) );
    }

    // pick best
    File logFile( outputPath + "/confOpt.txt", FILE_APPEND, FILE_TEXT );
    if (logFile.openSuccess()) {
        logFile.writeF( "\nbest: %s, %f\n\n", bestParamUpdate.c_str(), bestScore );
    }
}


} // end namespace sbl

