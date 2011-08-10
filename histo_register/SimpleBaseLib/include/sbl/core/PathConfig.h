#ifndef _SBL_PATH_CONFIG_H_
#define _SBL_PATH_CONFIG_H_
#include <sbl/core/Config.h>
namespace sbl {


/*! \file PathConfig.h
    \brief The PathConfig module handles a path.conf object that specifies
    the location(s) of application-specific data and other global configuration
    parameters.
*/


/// the program's main configuration file
/// (with data paths and other application-specific parameters);
/// loaded from the application's directory
Config &pathConfig();


/// save the main config file
void savePathConfig();


/// the application's log path, as defined in the main config
String logPath();


/// the application's data path, as defined in the main config
String dataPath();


/// if relative path, adds data path; otherwise returns unmodified
String addDataPath( const String &fileName );


/// set the application's data path, so that we don't need the main config
/// (if this isn't called before the first call to dataPath(), the earlier dataPath() call will use the main config)
void setDataPath( const String &dataPath );


} // end namespace sbl
#endif // _SBL_PATH_CONFIG_H_

