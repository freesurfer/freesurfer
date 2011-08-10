// Licensed under MIT license; see license.txt.

#include <sbl/core/Config.h>
#include <sbl/core/Display.h>
#include <sbl/core/StringUtil.h>
namespace sbl {


//-------------------------------------------
// CONFIG CLASS
//-------------------------------------------


// basic constructor
Config::Config() {
    m_allowDefault = false;
    m_allowNewEntries = true;
    m_commandArgPosition = 0;
    m_dualPass = false;
    m_initialPass = false;
    m_checkedInitialPass = false;
    m_missingValue = false;
}


/// read value from parameter with given name
bool Config::readBool( const String &name, int defaultValue ) {
    const ConfigEntry *configEntry = findEntry( name );
    if (configEntry == NULL) {
        if (m_commandArgPosition < m_commandArgs.count())
            return m_commandArgs[ m_commandArgPosition++ ].toBool();
        m_missingValue = true;
        if (m_allowDefault) {
            if (defaultValue == -7777) {
                if (m_initialPass == false)
                    warning( "Config::readBool: entry not found: %s", name.c_str() );
                if (m_dualPass)
                    writeBool( name, false );
                return false;
            } else {
                if (m_dualPass) 
                    writeBool( name, defaultValue ? true : false );
                return defaultValue ? true : false;
            }
        } else {
            fatalError( "Config::readBool: entry not found: %s", name.c_str() );
        }
    }
    return configEntry->value.toBool();
}


/// read value from parameter with given name
int Config::readInt( const String &name, int defaultValue ) {
    const ConfigEntry *configEntry = findEntry( name );
    if (configEntry == NULL) {
        if (m_commandArgPosition < m_commandArgs.count())
            return m_commandArgs[ m_commandArgPosition++ ].toInt();
        m_missingValue = true;
        if (m_allowDefault) {
            if (defaultValue == -7777) {
                if (m_initialPass == false)
                    warning( "Config::readInt: entry not found: %s; no default specified", name.c_str() );
                if (m_dualPass) 
                    writeInt( name, 0 );
                return 0;
            } else {
                if (m_dualPass) 
                    writeInt( name, defaultValue );
                return defaultValue;
            }
        } else {
            fatalError( "Config::readInt: entry not found: %s", name.c_str() );
        }
    }
    return configEntry->value.toInt();
}


/// read value from parameter with given name
float Config::readFloat( const String &name, float defaultValue ) {
    return (float) readDouble( name, defaultValue );
}


/// read value from parameter with given name
double Config::readDouble( const String &name, double defaultValue ) {
    const ConfigEntry *configEntry = findEntry( name );
    if (configEntry == NULL) {
        if (m_commandArgPosition < m_commandArgs.count())
            return m_commandArgs[ m_commandArgPosition++ ].toDouble();
        m_missingValue = true;
        if (m_allowDefault) {
            if (defaultValue > -7777.78 && defaultValue < -7777.76) {
                if (m_initialPass == false)
                    warning( "Config::readDouble: entry not found: %s; no default specified", name.c_str() );
                if (m_dualPass) 
                    writeDouble( name, 0.0 );
                return 0.0f;
            } else {
                if (m_dualPass) 
                    writeDouble( name, defaultValue );
                return defaultValue;
            }
        } else {
            fatalError( "Config::readDouble: entry not found: %s", name.c_str() );
        }
    }
    return configEntry->value.toDouble();
}


/// read value from parameter with given name
String Config::readString( const String &name, const String &defaultValue ) {
    const ConfigEntry *configEntry = findEntry( name );
    if (configEntry == NULL) {
        if (m_commandArgPosition < m_commandArgs.count())
            return m_commandArgs[ m_commandArgPosition++ ];
        m_missingValue = true;
        if (m_allowDefault) {
            if (defaultValue == "[none]") {
                if (m_initialPass == false)
                    warning( "Config::readString: entry not found: %s", name.c_str() );
                if (m_dualPass) 
                    writeString( name, "" );
                return "";
            } else {
                if (m_dualPass) 
                    writeString( name, defaultValue );
                return defaultValue;
            }
        } else {
            fatalError( "Config::readString: entry not found: %s", name.c_str() );
        }
    }
    return configEntry->value;
}


/// write value to parameter with given name
void Config::writeString( const String &name, const String &val, ConfigEntryType type ) {
    if (name.length() == 0) {
        warning( "config entry name: empty string" );
        return;
    }
    ConfigEntry *configEntry = findEntry( name );
    if (configEntry == NULL) {
        if (m_allowNewEntries) {
            configEntry = new ConfigEntry;
            configEntry->name = name;
            if (type == CONFIG_ENTRY_TEXT) {
                if (name.endsWith( "fileName" ) || name.endsWith( "FileName" )) {
                    type = CONFIG_ENTRY_FILE;
                } else if (name.endsWith( "path" ) || name.endsWith( "Path" )) {
                    type = CONFIG_ENTRY_PATH;
                }
            }
            configEntry->type = type;
            m_configEntries.append( configEntry );
        } else {
            warning( "failed attempt to add new config entry: %s=%s", name.c_str(), val.c_str() );
            return;
        }
    }
    configEntry->value = val;
}


/// load config from text file
// fix(later): more error checking
bool Config::load( const String &fileName ) {
    File file( fileName, FILE_READ, FILE_TEXT );
    while (file.endOfFile() == false) {
        String line = file.readLine().strip();

        // if blank line, store blank config entry (so we can save the config with blank lines included)
        if (line.length() == 0) {
            if (file.endOfFile() == false) {
                ConfigEntry *configEntry = new ConfigEntry;
                configEntry->type = CONFIG_ENTRY_BLANK;
                m_configEntries.append( configEntry );
            }

        // if not blank, try to parse the line
        } else {
            ConfigEntry *configEntry = new ConfigEntry;

            // split line into body and comment
            String body, comment;
            int hashPos = line.firstCharPos( '#' );
            if (hashPos >= 0 && hashPos + 1 < line.length()) {
                body = line.leftOf( hashPos ).strip();
                comment = line.rightOf( hashPos ).strip();
            } else {
                body = line.strip();
            }

            // if section line (no name or value, just comment)
            if (body.length() == 0 && comment.length() > 0) {
                configEntry->name = comment;
                configEntry->type = CONFIG_ENTRY_SECTION;

            // if normal name/value line
            } else {
                int spacePos = body.firstCharPos( ' ' );
                if (spacePos > 0) {

                    // parse meta-data tags (if any)
                    int leftBracketPos = body.firstCharPos( '[' );
                    int rightBracketPos = body.lastCharPos( ']' );
                    if (leftBracketPos >= 0 && rightBracketPos >= 0) {
                        String meta = body.leftOf( rightBracketPos ).rightOf( leftBracketPos );
                        if (meta == "bool")
                            configEntry->type = CONFIG_ENTRY_BOOL;
                        else if (meta == "path")
                            configEntry->type = CONFIG_ENTRY_PATH;
                        else if (meta == "file")
                            configEntry->type = CONFIG_ENTRY_FILE;
                        body = body.leftOf( leftBracketPos );
                    }

                    // store main entry properties
                    configEntry->name = body.leftOf( spacePos ).strip();
                    configEntry->value = body.rightOf( spacePos ).strip();
                    configEntry->description = comment;
                }
            }

            // if entry is good, store it
            if (configEntry->name.length() && (configEntry->value.length() || configEntry->type == CONFIG_ENTRY_SECTION)) {
                m_configEntries.append( configEntry );
            } else {
                delete configEntry;
                warning( "invalid config entry" );
                return false;
            }
        }
    }
    return true;
}


/// save config to text file
void Config::save( const String &fileName ) const {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {
        for (int i = 0; i < m_configEntries.count(); i++) {
            const ConfigEntry &configEntry = m_configEntries[ i ];
            if (configEntry.name.length() && configEntry.value.length()) {
                String line = configEntry.name + " " + configEntry.value;
                if (configEntry.type == CONFIG_ENTRY_BOOL)
                    line += " [bool]";
                else if (configEntry.type == CONFIG_ENTRY_PATH)
                    line += " [path]";
                else if (configEntry.type == CONFIG_ENTRY_FILE)
                    line += " [file]";
                if (configEntry.description.length()) {
                    line += " # ";
                    line += configEntry.description;
                }
                line += "\n";
                file.writeRawString( line );
            } else if (configEntry.type == CONFIG_ENTRY_SECTION) {
                String line = "# ";
                line += configEntry.name;
                line += "\n";
                file.writeRawString( line );
            } else if (configEntry.type == CONFIG_ENTRY_BLANK) {
                file.writeRawString( "\n" );        
            } else {
                warning( "invalid config entry" );
            }
        }
    }
}


/// updates the config using this syntax: param1=value1+param2=value2+param3=value3
void Config::update( const String &update ) {
    Array<String> split = update.split( "+" );
    for (int i = 0; i < split.count(); i++) {
        const String &s = split[ i ];
        Array<String> subSplit = s.split( "=" );
        if (subSplit.count() == 2) {
            disp( 1, "update param: %s, value: %s", subSplit[ 0 ].c_str(), subSplit[ 1 ].c_str() );
            String val = subSplit[ 1 ].strip( '"' ); // strip quotes from value
            writeString( subSplit[ 0 ].c_str(), val );
        } else {
            warning( "invalid update string: %s", s.c_str() );
        }
    }
}


/// copy values for entries that are specified in other config
// fix(later): doesn't write to local if super
void Config::updateFrom( const Config &conf ) {
    for (int i = 0; i < m_configEntries.count(); i++) {
        ConfigEntry &thisConfigEntry = m_configEntries[ i ];
        const ConfigEntry *otherConfigEntry = conf.findEntry( thisConfigEntry.name ); 
        if (otherConfigEntry)
            thisConfigEntry.value = otherConfigEntry->value;
    }
}


/// set an array of command arguments, to be assigned to subsequent calls to read() functions
void Config::setCommandArgs( const Array<String> &commandArgs ) {
    for (int i = 0; i < commandArgs.count(); i++) {
        const String &arg = commandArgs[ i ];
        if (arg.contains( '=' )) {
            String name = arg.leftOfFirst( '=' );
            String value = arg.rightOfFirst( '=' );
            writeString( name, value );
        } else {
            m_commandArgs.appendCopy( commandArgs[ i ] );
        }
    }
}


/// find an entry by name
const ConfigEntry *Config::findEntry( const String &name ) const {
    assertDebug( name.length() );
    for (int i = 0; i < m_configEntries.count(); i++) {
        const ConfigEntry *configEntry = &m_configEntries[ i ];
        if (configEntry->name == name)
            return configEntry;
    }

    // not found
    return NULL;
}


/// find an entry by name
ConfigEntry *Config::findEntry( const String &name ) {
    assertDebug( name.length() );
    for (int i = 0; i < m_configEntries.count(); i++) {
        ConfigEntry *configEntry = &m_configEntries[ i ];
        if (configEntry->name == name)
            return configEntry;
    }

    // not found
    return NULL;
}


} // end namespace sbl

