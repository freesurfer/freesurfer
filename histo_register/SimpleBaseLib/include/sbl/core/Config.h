#ifndef _SBL_CONFIG_H_
#define _SBL_CONFIG_H_
#include <sbl/core/Array.h>
#include <sbl/core/File.h>
#include <sbl/core/StringUtil.h>
#include <sbl/math/VectorUtil.h>
namespace sbl {


// config entry types are currently used mostly for display/editing
enum ConfigEntryType {
    CONFIG_ENTRY_TEXT,
    CONFIG_ENTRY_BOOL,
    CONFIG_ENTRY_SECTION,
    CONFIG_ENTRY_BLANK,
    CONFIG_ENTRY_PATH,
    CONFIG_ENTRY_FILE
};


//-------------------------------------------
// CONFIG ENTRY CLASS
//-------------------------------------------


/// A ConfigEntry holds a named parameter in a Config object.
class ConfigEntry {
public:

    /// basic constructor
    ConfigEntry() { type = CONFIG_ENTRY_TEXT; }

    // public members
    // fix(clean): make private
    String name;
    String value;
    String description;
    ConfigEntryType type;
};


//-------------------------------------------
// CONFIG CLASS
//-------------------------------------------


/// The Config class represents a configuration file with a set of named parameters.
class Config {
public:

    // basic constructor
    Config();

    //-------------------------------------------
    // ACCESS ENTRIES
    //-------------------------------------------

    /// read value from parameter with given name
    bool readBool( const String &name, int defaultValue = -7777 );
    int readInt( const String &name, int defaultValue = -7777 );
    float readFloat( const String &name, float defaultValue = -7777.77f );
    double readDouble( const String &name, double defaultValue = -7777.77 );
    String readString( const String &name, const String &defaultValue = "[none]" );
    inline Array<String> readStrings( const String &name ) { return readString( name ).split( "," ); }
    inline VectorI readVectorI( const String &name ) { return splitNumbersI( readString( name, "" ) ); }
    inline VectorF readVectorF( const String &name ) { return splitNumbersF( readString( name, "" ) ); }

    /// write value to parameter with given name
    inline void writeBool( const String &name, bool val ) { writeString( name, sprintF( "%d", (int) val ), CONFIG_ENTRY_BOOL ); }
    inline void writeInt( const String &name, int val ) { writeString( name, sprintF( "%d", val ) ); }
    inline void writeFloat( const String &name, float val ) { writeString( name, sprintF( "%f", val ) ); }
    inline void writeDouble( const String &name, double val ) { writeString( name, sprintF( "%f", val ) ); }
    void writeString( const String &name, const String &val, ConfigEntryType type = CONFIG_ENTRY_TEXT );
    inline void writeStrings( const String &name, const Array<String> &strArr ) { writeString( name, join( strArr, "," ) ); }
    inline void writeVectorI( const String &name, const VectorI &vector ) { writeString( name, toString( vector ) ); }
    inline void writeVectorF( const String &name, const VectorF &vector ) { writeString( name, toString( vector, 0, 1000, "%f" ) ); }
    inline void writeSection( const String &name ) { writeString( name, "", CONFIG_ENTRY_SECTION ); }

    /// checks whether entry exists in this config
    inline bool entryExists( const String &name ) const { return findEntry( name ) ? true : false; }

    /// the number of entries
    inline int entryCount() const { return m_configEntries.count(); }

    /// access entry by index
    inline const ConfigEntry &entry( int index ) const { return m_configEntries[ index ]; }
    inline ConfigEntry &entry( int index ) { return m_configEntries[ index ]; }

    /// remove all entries (note: doesn't clear command args)
    inline void reset() { m_configEntries.reset(); }

    //-------------------------------------------
    // LOAD/SAVE
    //-------------------------------------------

    /// load config from text file
    bool load( const String &fileName );

    /// save config to text file
    void save( const String &fileName ) const;

    //-------------------------------------------
    // OTHER METHODS
    //-------------------------------------------

    /// allow write() calls to create new entries
    inline void allowNewEntries( bool allow ) { m_allowNewEntries = allow; }

    /// allow read() calls to return default values if entry not found
    inline void allowDefault( bool allow ) { m_allowDefault = allow; }

    /// enabled dual-pass mode (first pass lets command specify args for GUI)
    inline void enableDualPass() { m_dualPass = true; }

    /// returns true if first pass (in dual-pass mode), meaning command should terminate after reading params (it will be run again in for the second pass)
    inline bool initialPass() { m_checkedInitialPass = true; return m_initialPass; }

    /// set first/second pass (in dual-pass mode)
    // fix(clean): better to write args+names first time around rather than reset arg pos?
    inline void setInitialPass( bool initialPass ) { m_initialPass = initialPass; m_commandArgPosition = 0; }

    /// true if checked whether initial pass
    inline bool checkedInitialPass() { return m_checkedInitialPass; }

    /// true if encountered unspecified parameters
    inline bool missingValue() const { return m_missingValue; }
    
    /// updates the config using this syntax: param1=value1+param2=value2+param3=value3
    void update( const String &update );

    /// copy values for entries that are specified in other config
    void updateFrom( const Config &conf );

    /// set an array of command arguments, to be assigned to subsequent calls to read() functions
    void setCommandArgs( const Array<String> &commandArgs );

    /// an array of command arguments, to be assigned to subsequent calls to read() functions
    inline const Array<String> &commandArgs() const { return m_commandArgs; }

private:

    /// find an entry by name
    const ConfigEntry *findEntry( const String &name ) const;
    ConfigEntry *findEntry( const String &name );

    // the entries
    Array<ConfigEntry> m_configEntries;

    // if specified, fill in config entries from command args
    Array<String> m_commandArgs;
    mutable int m_commandArgPosition;

    // allow write() calls to create new entries
    bool m_allowNewEntries;

    // allow read() calls to return default values if entry not found
    bool m_allowDefault;

    // dual-pass mode: first pass lets command specify args for GUI, second pass actually runs the command
    bool m_dualPass;

    // true if first pass (in dual-pass mode), meaning command should terminate after reading params (it will be run again in for the second pass)
    bool m_initialPass;

    // true if checked whether initial pass
    bool m_checkedInitialPass;

    // true if encountered unspecified parameters (even if defaults specified)
    bool m_missingValue;

    // disable copy constructor and assignment operator
    Config( const Config &x );
    Config &operator=( const Config &x );
};


} // end namespace sbl
#endif // _SBL_CONFIG_H_

