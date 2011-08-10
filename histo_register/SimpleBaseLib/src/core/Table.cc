// Licensed under MIT license; see license.txt.

#include <sbl/core/Table.h>
#include <sbl/core/Command.h>
#include <sbl/system/FileSystem.h>
#include <sbl/system/TimeUtil.h>
#include <sbl/math/VectorUtil.h>
namespace sbl {


// used to automatically identify timestamp columns
#define MIN_TIMESTAMP_YEAR 2000
#define MAX_TIMESTAMP_YEAR 2030


//-------------------------------------------
// TABLE COLUMN CLASS
//-------------------------------------------


// basic constructor
TableColumn::TableColumn( const String &name ) {

    // init meta-data
    m_name = name;
    m_nameHash = strHash( name );
    m_visible = true;
    m_isTimestamp = false;
    m_highlightMinMax = false;

    // init stats
    m_nonZeroCount = 0;
    m_minInt = 0;
    m_maxInt = 0;
    m_meanInt = 0;
    m_minDouble = 0.0;
    m_maxDouble = 0.0;
    m_meanDouble = 0.0;
}


/// load from tagged file
TableColumn::TableColumn( TaggedFile &file ) {
    m_visible = true;
    m_isTimestamp = false;
    m_highlightMinMax = false;
    int tag = 0;
    do {
        tag = file.readTag();
        switch (tag) {
        case TAG_NAME: m_name = file.readString(); break;
        case TAG_VISIBLE: m_visible = file.readBool(); break;
        case TAG_TOOL_TIP: m_toolTip = file.readString(); break;
        case TAG_FORMAT: m_format = file.readString(); break;
        case TAG_IS_TIMESTAMP: m_isTimestamp = file.readBool(); break;
        case TAG_HIGHLIGHT_MIN_MAX: m_highlightMinMax = file.readBool(); break;
        case TAG_DATA_INT: m_dataInt = file.readVector<int>(); break;
        case TAG_DATA_DOUBLE: m_dataDouble = file.readVector<double>(); break;
        case TAG_DATA_STRING: file.readStrings( m_dataString ); break;
        default: file.skipTag( tag );
        }
    } while (tag != TAG_END_SECTION);
    m_visible = true; // fix(later): remove this

    // init other properties
    m_nameHash = strHash( m_name );
    m_nonZeroCount = 0;
    m_minInt = 0;
    m_maxInt = 0;
    m_meanInt = 0;
    m_minDouble = 0.0;
    m_maxDouble = 0.0;
    m_meanDouble = 0.0;
}


/// save to tagged file
void TableColumn::save( TaggedFile &file ) const {

    // save meta-data
    file.tag( TAG_NAME ).writeString( m_name );
    file.tag( TAG_VISIBLE ).writeBool( m_visible );
    file.tag( TAG_TOOL_TIP ).writeString( m_toolTip );
    file.tag( TAG_FORMAT ).writeString( m_format );
    file.tag( TAG_IS_TIMESTAMP ).writeBool( m_isTimestamp );
    file.tag( TAG_HIGHLIGHT_MIN_MAX ).writeBool( m_highlightMinMax );

    // save data
    file.tag( TAG_DATA_INT ).writeVector( m_dataInt );
    file.tag( TAG_DATA_DOUBLE ).writeVector( m_dataDouble );
    file.tag( TAG_DATA_STRING ).writeStrings( m_dataString );
}


/// the number of items in the column
int TableColumn::pointCount() const {
    int pointCount = 0;
    if (m_dataInt.length())
        pointCount = m_dataInt.length();
    if (m_dataDouble.length())
        pointCount = m_dataDouble.length();
    if (m_dataString.count())
        pointCount = m_dataString.count();
    return pointCount;
}


/// compute value statistics for display
void TableColumn::computeStats() {
    if (m_dataInt.length()) {
        m_nonZeroCount = nonZeroCount( m_dataInt );
        m_minInt = m_dataInt.min();
        m_maxInt = m_dataInt.max();
        m_meanInt = m_dataInt.mean();

        // fix(later): must be better way to do this; create timestamp type?
        if (m_minInt > dateUTC( 2000, 1, 1 ) && m_maxInt < dateUTC( 2030, 1, 1 )) 
            m_isTimestamp = true;
    }
    if (m_dataDouble.length()) {
        m_nonZeroCount = nonZeroCount( m_dataDouble );
        m_minDouble = m_dataDouble.min();
        m_maxDouble = m_dataDouble.max();
        m_meanDouble = m_dataDouble.mean();
    }
}


/// web-formatted stats
String TableColumn::webStatText() const {
    String result;
    if (m_dataInt.length()) {
        if (m_isTimestamp) {
            result = sprintF( "min: %s<br>mean: %s<br>max: %s<br>frac: %1.2f<br>count: %d",
                tsToStrLocal( m_minInt, true ).c_str(), tsToStrLocal( m_meanInt, true ).c_str(), tsToStrLocal( m_maxInt, true ).c_str(), (float) m_nonZeroCount / (float) m_dataInt.length(), m_dataInt.length() );
        } else {
            result = sprintF( "min: %d<br>mean: %d<br>max: %d<br>frac: %1.2f<br>count: %d",
                m_minInt, m_meanInt, m_maxInt, (float) m_nonZeroCount / (float) m_dataInt.length(), m_dataInt.length() );
        }
    } 
    if (m_dataDouble.length()) {
        result = sprintF( "min: %1.4f<br>mean: %1.4f<br>max: %1.4f<br>frac: %1.2f<br>count: %d",
            m_minDouble, m_meanDouble, m_maxDouble, (float) m_nonZeroCount / (float) m_dataDouble.length(), m_dataDouble.length() );
    }
    return result;
}


/// value formatted for display 
String TableColumn::displayValue( int index, String &target, bool alignRight, String &color ) const {
    String value;
    if (index >= 0) {

        // handle int data
        if (index < m_dataInt.length()) {
            if (m_isTimestamp)
                value = tsToStrLocal( m_dataInt[ index ], true );
            else
                value = sprintF( "%d", m_dataInt[ index ] );
            if (m_highlightMinMax && m_minInt != m_maxInt) {
                if (m_dataInt[ index ] == m_minInt) 
                    color = "FFCCCC";
                else if (m_dataInt[ index ] == m_maxInt) 
                    color = "CCFFCC";
            }
        }

        // handle double data
        if (index < m_dataDouble.length()) {
            if (m_format.length())
                value = sprintF( m_format.c_str(), m_dataDouble[ index ] );
            else
                value = sprintF( "%1.4f", m_dataDouble[ index ] );
            if (m_highlightMinMax && m_minDouble != m_maxDouble) {
                if (m_dataDouble[ index ] == m_minDouble) 
                    color = "FFCCCC";
                else if (m_dataDouble[ index ] == m_maxDouble) 
                    color = "CCFFCC";
            }
        }

        // handle string data
        if (index < m_dataString.count()) {
            value = String( m_dataString[ index ] );
            if (value.contains( "||" )) {
                
                // fix(later): assuming no extra pipes
                int pipePos = value.firstCharPos( '|' );
                target = value.rightOf( pipePos + 1 );
                value = value.leftOf( pipePos );
            }
        }
    }

    // determine left/right alignment
    alignRight = m_dataInt.length() || m_dataDouble.length();

    // return value to display
    return value;
}


/// value formatted for export
String TableColumn::simpleDisplayValue( int index ) const {
    String value;
    if (index < m_dataInt.length()) {
        if (m_isTimestamp)
            value = tsToStrLocal( m_dataInt[ index ], true );
        else
            value = sprintF( "%d", m_dataInt[ index ] );
    }
    if (index < m_dataDouble.length()) {
        if (m_format.length())
            value = sprintF( m_format.c_str(), m_dataDouble[ index ] );
        else
            value = sprintF( "%4.2f", m_dataDouble[ index ] );
    }
    if (index < m_dataString.count()) {
        value = String( m_dataString[ index ] );
        if (value.contains( "||" )) 
            value = value.leftOfFirst( '|' );
    }
    return value;
}


/// value formatted for export
String TableColumn::exportValue( int index ) const {
    String value;
    if (index < m_dataInt.length())
        value = sprintF( "%d", m_dataInt[ index ] );
    if (index < m_dataDouble.length()) {
           value = sprintF( "%g", m_dataDouble[ index ] );
        if (value.contains( '.' ) == false) // fix(later): this is a hack to make sure double fields always have decimal point, so that we can load the csv correctly
            value += ".0";
    }
    if (index < m_dataString.count())
        value = String( m_dataString[ index ] );
    return value;
}


//-------------------------------------------
// TABLE CLASS
//-------------------------------------------


// basic constructor
Table::Table( const String &title ) {
    m_title = title;
    m_startIndex = 0;
    m_dispLength = 50;
    m_timestampColumn = -1;
}


/// access column by name
const TableColumn &Table::column( const String &colName ) const {
    int nameHash = strHash( colName );
    for (int i = 0; i < m_columns.count(); i++) 
        if (m_columns[ i ].nameHash() == nameHash) 
            return m_columns[ i ];
    fatalError( "column not found: %s", colName.c_str() );
    return m_columns[ 0 ];
}


/// access column by name
TableColumn &Table::column( const String &colName ) {
    int nameHash = strHash( colName );
    for (int i = 0; i < m_columns.count(); i++) 
        if (m_columns[ i ].nameHash() == nameHash) 
            return m_columns[ i ];
    fatalError( "column not found: %s", colName.c_str() );
    return m_columns[ 0 ];
}


/// true if column is already defined
bool Table::columnDefined( const String &colName ) const {
    int nameHash = strHash( colName );
    for (int i = 0; i < m_columns.count(); i++) 
        if (m_columns[ i ].nameHash() == nameHash) 
            return true;
    return false;
}


/// read a value from a table entry
int Table::readInt( const char *colName, int rowIndex ) const { 
    const TableColumn &col = column( colName );
    assertAlways( rowIndex >= 0 && rowIndex < col.dataI().length() );
    return col.dataI()[ rowIndex ]; 
}


/// read a value from a table entry
double Table::readDouble( const char *colName, int rowIndex ) const { 
    const TableColumn &col = column( colName );
    assertAlways( rowIndex >= 0 && rowIndex < col.dataD().length() );
    return col.dataD()[ rowIndex ]; 
}


/// read a value from a table entry
String Table::readString( const char *colName, int rowIndex ) const { 
    const TableColumn &col = column( colName );
    assertAlways( rowIndex >= 0 && rowIndex < col.dataStr().count() );
    return col.dataStr()[ rowIndex ]; 
}


/// find a column that contains timestamps
void Table::findTimestampColumn() {
    for (int i = 0; i < m_columns.count(); i++) {
        TableColumn &col = m_columns[ i ];
        if (col.isTimestamp()) {
            m_timestampColumn = i;
            break;
        }
    }
}


/// compute stats for all columns
void Table::computeStats() {
    for (int i = 0; i < m_columns.count(); i++) 
        m_columns[ i ].computeStats();
}


/// retrieve column with given name; create if doesn't exist
TableColumn &Table::findAddColumn( const char *colName, bool &newCol ) {
    newCol = false;
    int nameHash = strHash( colName );
    for (int i = 0; i < m_columns.count(); i++) 
        if (m_columns[ i ].nameHash() == nameHash) 
            return m_columns[ i ];
    
    // if not found, add
    TableColumn *column = new TableColumn( colName );
    m_columns.append( column );
    newCol = true;
    return *column;
}


/// add a double value
void Table::add( const String &colName, double val, bool highlightMinMax, const String &format ) {
    bool newCol = false;
    TableColumn &column = findAddColumn( colName.c_str(), newCol );
    column.add( (float) val ); // fix(later): store double data?

    // if new column, set display properties
    if (newCol) {
        column.setHighlightMinMax( highlightMinMax );
        if (format.length())
            column.setFormat( format );
        else
            column.setFormat( m_defaultFormat );
    }
}


/// add a float value
void Table::add( const String &colName, float val, bool highlightMinMax, const String &format ) {
    bool newCol = false;
    TableColumn &column = findAddColumn( colName.c_str(), newCol );
    column.add( val );

    // if new column, set display properties
    if (newCol) {
        column.setHighlightMinMax( highlightMinMax );
        if (format.length())
            column.setFormat( format );
        else
            column.setFormat( m_defaultFormat );
    }
}


/// add an int value
void Table::add( const String &colName, int val, bool highlightMinMax, const String &format ) {
    bool newCol = false;
    TableColumn &column = findAddColumn( colName.c_str(), newCol );
    column.add( val );

    // if new column, set display properties
    if (newCol) {
        column.setHighlightMinMax( highlightMinMax );
        if (format.length())
            column.setFormat( format );
        else
            column.setFormat( m_defaultFormat );
    }
}


/// add a string value
void Table::add( const String &colName, const String &val ) {
    bool newCol = false;
    TableColumn &column = findAddColumn( colName.c_str(), newCol );
    column.add( val );
}


/// add a string value with a command
void Table::add( const String &colName, const String &val, const String &cmd ) {
    add( colName, val + "||" + cmd );
}


/// add a blank column between columns
void Table::addBlank( int index ) {
    add( sprintF( "[blank:%d]", index ).c_str(), "" );
}


/// the number of rows (assumes all columns have the same number of rows)
int Table::pointCount() const {
    if (m_columns.count())
        return m_columns[ 0 ].pointCount();
    return 0;
}


/// adjust start index by specified amount
void Table::step( int step ) {
    int pointCount = m_columns[ 0 ].pointCount();
    switch (step) {
    case -2:
        m_startIndex = 0;
        break;
    case -1:
        m_startIndex -= m_dispLength;
        break;
    case 1:
        m_startIndex += m_dispLength;
        break;
    case 2:
        m_startIndex = pointCount - m_dispLength;
        break;
    }
    if (m_startIndex > pointCount - m_dispLength)
        m_startIndex = pointCount - m_dispLength;
    if (m_startIndex < 0) 
        m_startIndex = 0;
}


/// compute sort index using the specified column 
void Table::sort( const TableColumn &col, int sortDir ) {
    if (col.dataI().length()) {
        if (sortDir > 0)
            m_sortIndex = sbl::sortIndex( col.dataI() );
        else
            m_sortIndex = reverseSortIndex( col.dataI() );
    }
    if (col.dataD().length()) {
        if (sortDir > 0)
            m_sortIndex = sbl::sortIndex( col.dataD() );
        else
            m_sortIndex = reverseSortIndex( col.dataD() );
    }
}


/// load from tagged file
void Table::load( const String &fileName ) {
    TaggedFile file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess()) {
        int tag = 0;
        do {
            tag = file.readTag();
            switch (tag) {
            case TAG_TITLE: m_title = file.readString(); break;
            case TAG_COLUMNS: file.readArray( m_columns ); break;
            default: file.skipTag( tag );
            }
        } while (tag != TAG_END_SECTION);        
    }
}


/// save to tagged file
void Table::save( const String &fileName ) const {
    TaggedFile file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess()) {
        file.tag( TAG_TITLE ).writeString( m_title );
        file.tag( TAG_COLUMNS ).writeArray( m_columns );
    }
}


/// true if string appears to be numeric
// fix(later): do better, move to util
bool isNumeric( const String &str ) {
    bool foundNumeric = false;
    for (int i = 0; i < str.length(); i++) {
        int c = str.get( i );
        if (c >= 'A' && c != 'E' && c != 'e')
            return false;
        if (c >= '0' && c <= '9')
            foundNumeric = true;
    }
    return foundNumeric;
}


// column types used for CSV parsing
#define COL_TYPE_STRING 0
#define COL_TYPE_INT 1
#define COL_TYPE_DOUBLE 2


/// load from comma-separated file
void Table::loadCSV( const String &fileName ) {
    File file( fileName, FILE_READ, FILE_TEXT );
    if (file.openSuccess()) {
    
        // get labels
        Array<String> labels = file.readLine().split( "," );
        int colCount = labels.count();
        VectorI colType( colCount );
        colType.clear( COL_TYPE_STRING );
        bool firstLine = true;

        // write data
        while (file.endOfFile() == false) {
            Array<String> split = file.readLine().split( "," );
            if (split.count() == colCount) {

                // if first line, determine type of each column
                if (firstLine) {
                    for (int i = 0; i < split.count(); i++) {
                        if (isNumeric( split[ i ] )) {
                            if (split[ i ].firstCharPos( '.' ) >= 0) {
                                colType[ i ] = COL_TYPE_DOUBLE;
                            } else {
                                colType[ i ] = COL_TYPE_INT;
                            }
                        }
                    }
                    firstLine = false;
                }

                // add data to table
                for (int i = 0; i < split.count(); i++) {
                    if (colType[ i ] == COL_TYPE_DOUBLE) {
                        add( labels[ i ].c_str(), split[ i ].toDouble() );
                    } else if (colType[ i ] == COL_TYPE_INT) {
                        add( labels[ i ].c_str(), split[ i ].toInt() );
                    } else {
                        add( labels[ i ].c_str(), split[ i ] );
                    }
                }
            } else if (split.count() > 1) {
                warning( "%s: discarding line", fileName.c_str() );
            }
        }
    } else {
        warning( "unable to open: %s", fileName.c_str() );
    }
}


/// save to comma-separated file
void Table::saveCSV( const String &fileName, bool append ) const {
    bool skipHeader = append && fileExists( fileName );
    File file( fileName, append ? FILE_APPEND : FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {
        
        // write headers
        if (skipHeader == false) {
            for (int colIndex = 0; colIndex < columnCount(); colIndex++) {
                const TableColumn &col = column( colIndex );
                file.writeF( "%s", col.name().c_str() );
                if (colIndex + 1 < columnCount())
                    file.writeF( "," );
            }
            file.writeF( "\n" );
        }

        // write data
        int count = pointCount();
        for (int i = 0; i < count; i++) {
            for (int colIndex = 0; colIndex < columnCount(); colIndex++) {
                const TableColumn &col = column( colIndex );
                file.writeF( col.exportValue( i ).c_str() );
                if (colIndex + 1 < columnCount())
                    file.writeF( "," );
            }
            file.writeF( "\n" );
        }
    }
}


/// display a table using the currently register callback (takes ownership)
void Table::disp( aptr<Table> table ) { 
    if (s_tableDisplay) 
        s_tableDisplay->display( table ); 
}


// init table disp callback object (used by all table instances)
Display< aptr<Table> > *Table::s_tableDisplay = NULL;


} // end namespace sbl

