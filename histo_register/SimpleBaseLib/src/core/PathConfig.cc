// Licensed under MIT license; see license.txt.

#include <sbl/core/PathConfig.h>
#include <sbl/core/File.h>
#include <sbl/system/FileSystem.h>
namespace sbl {


/// a global config instance returned by pathConfig()
Config g_pathConfig;


/// the program's main configuration file
/// (with data paths and other application-specific parameters);
/// loaded from the application's directory
Config &pathConfig() {
    if (g_pathConfig.entryCount() == 0) {
        g_pathConfig.load( "path.conf" );
        if (g_pathConfig.entryCount() == 0)
            fatalError( "unable to open path.conf" );
    }
    return g_pathConfig;
}


/// save the main config file (to the application's directory)
void savePathConfig() {
    g_pathConfig.save( "path.conf" );
}



// a global instance of the log path, loaded from the path config
String g_logPath;


/// the application's log path, as defined in the path config
String logPath() {
    if (g_logPath.length() == 0) {
        g_logPath = pathConfig().readString( "logPath" );

        // make sure path ends with a slash (or backslash)
        if (g_logPath.endsWith( "/" ) == false || g_logPath.endsWith( "\\" ) == false)
            g_logPath += "/";
    }
    return g_logPath;
}


// a global instance of the data path, loaded from the main config
String g_dataPath;


/// the application's data path, as defined in the path config
String dataPath() {
    if (g_dataPath.length() == 0) {
        g_dataPath = pathConfig().readString( "dataPath" );

        // make sure path ends with a slash (or backslash)
        if (g_dataPath.endsWith( "/" ) == false || g_dataPath.endsWith( "\\" ) == false)
            g_dataPath += "/";
    }
    return g_dataPath;
}


/// if relative path, adds data path; otherwise returns unmodified
String addDataPath( const String &fileName ) {
    String result;
    if (isAbsoluteFileName( fileName ))
        result = fileName;
    else 
        result = dataPath() + fileName;
    return result;
}


/// set the application's data path, so that we don't need the path config
/// (if this isn't called before the first call to dataPath(), the earlier dataPath() call will use the main config)
void setDataPath( const String &dataPath ) {
    g_dataPath = dataPath;
}


} // end namespace sbl

