#ifndef _SBL_TABLE_H_
#define _SBL_TABLE_H_
#include <sbl/core/String.h>
#include <sbl/core/Array.h>
#include <sbl/core/Pointer.h>
#include <sbl/other/TaggedFile.h>
#include <sbl/math/Vector.h>
namespace sbl {


//-------------------------------------------
// TABLE COLUMN CLASS 
//-------------------------------------------


/// A TableColumn stores a named sequences of values in a Table.
class TableColumn {
public:

    // basic constructor 
    explicit TableColumn( const String &name );

    /// file load constructor
    explicit TableColumn( TaggedFile &file );

    /// the column's name name
    inline const String &name() const { return m_name; }

    /// a hash of the name 
    // fix(later): make unique? (we are assuming it is unique for now)
    inline int nameHash() const { return m_nameHash; }

    /// save to file
    void save( TaggedFile &file ) const;

    //-------------------------------------------
    // ACCESS COLUMN DATA 
    //-------------------------------------------

    /// direct data access
    inline const VectorI &dataI() const { return m_dataInt; } // fix(clean): do we need these vector versions?
    inline const VectorD &dataD() const { return m_dataDouble; }
    inline const Array<String> &dataStr() const { return m_dataString; }
    inline Array<String> &dataStr() { return m_dataString; }

    /// the number of items in the column
    int pointCount() const;

    /// add data
    inline void add( int val ) { m_dataInt.append( val ); }
    inline void add( double val ) { m_dataDouble.append( val ); }
    inline void add( const String &val ) { m_dataString.append( new String( val )); }

    /// true if the column contains timestamps (assuming timestamps are integer valued)
    inline bool isTimestamp() const { return m_isTimestamp; }

    /// string representations of the data
    String displayValue( int index, String &target, bool alignRight, String &color ) const;
    String simpleDisplayValue( int index ) const; // fix(later): condense these
    String exportValue( int index ) const;

    //-------------------------------------------
    // DATA SUMMARY STATS
    //-------------------------------------------

    /// get stats about data
    inline int minInt() const { return m_minInt; }
    inline int maxInt() const { return m_maxInt; }
    inline int meanInt() const { return m_meanInt; }
    inline double minDouble() const { return m_minDouble; }
    inline double maxDouble() const { return m_maxDouble; }
    inline double meanDouble() const { return m_meanDouble; }

    /// compute value statistics for display
    void computeStats();

    //-------------------------------------------
    // DISPLAY PROPERTIES
    //-------------------------------------------

    /// access display info
    inline bool visible() const { return m_visible; }
    inline void setVisible( bool visible ) { m_visible = visible; }
    String webStatText() const;    

    /// if true, highlight the min/max values
    inline bool highlightMinMax() const { return m_highlightMinMax; }
    inline void setHighlightMinMax( bool v ) { m_highlightMinMax = v; }

    /// get/set the printf format for this column's data (e.g. "%4.2f")
    inline const String &format() const { return m_format; }
    inline void setFormat( const String &format ) { m_format = format; }

    /// get/set a text description 
    inline const String &description() const { return m_description; }
    inline void setDescription( const String &description ) { m_description = description; }

private:

    // column's data values
    VectorI m_dataInt;
    VectorD m_dataDouble;
    Array<String> m_dataString;
    bool m_isTimestamp;

    // stats about data
    int m_nonZeroCount;
    int m_minInt;
    int m_maxInt;
    int m_meanInt;
    double m_minDouble;
    double m_maxDouble;
    double m_meanDouble;

    /// column's display info
    String m_name;
    int m_nameHash;
    bool m_visible;
    String m_toolTip;
    bool m_highlightMinMax;
    String m_format;
    String m_description;

    // file format tags
    enum {
        TAG_NAME = 101,
        TAG_VISIBLE = 200,
        TAG_TOOL_TIP = 201,
        TAG_FORMAT = 202,
        TAG_IS_TIMESTAMP = 203,
        TAG_HIGHLIGHT_MIN_MAX = 204,
        TAG_DATA_INT = 301,
        TAG_DATA_DOUBLE = 302,
        TAG_DATA_STRING = 303
    };

    // disable copy constructor and assignment operator
    TableColumn( const TableColumn &x );
    TableColumn &operator=( const TableColumn &x );
};


//-------------------------------------------
// TABLE CLASS 
//-------------------------------------------


/// The Table class consists of a set of named columns (like a CSV file), 
/// with mechanisms for viewing and interacting with the data.
// fix(clean): move paging and other web stuff into WebTableView class? (or just remove?)
class Table {
public:

    /// create an empty table
    explicit Table( const String &title = "[none]" );

    /// get/set the table title
    inline const String &title() const { return m_title; }
    inline void setTitle( const String &title ) { m_title = title; }

    //-------------------------------------------
    // READ/WRITE DATA
    //-------------------------------------------

    /// access columns
    inline const TableColumn &column( int index ) const { return m_columns[ index ]; }
    inline TableColumn &column( int index ) { return m_columns[ index ]; }
    inline int columnCount() const { return m_columns.count(); }
    inline int timestampColumn() const { return m_timestampColumn; }
    const TableColumn &column( const String &colName ) const;
    TableColumn &column( const String &colName );
    bool columnDefined( const String &colName ) const;

    /// read a value
    int readInt( const char *colName, int rowIndex ) const;
    double readDouble( const char *colName, int rowIndex ) const;
    String readString( const char *colName, int rowIndex ) const;

    /// add a value
    void add( const String &colName, int val, bool highlightMinMax = false, const String &format = "" );
    void add( const String &colName, float val, bool highlightMinMax = false, const String &format = "" );
    void add( const String &colName, double val, bool highlightMinMax = false, const String &format = "" );
    void add( const String &colName, const String &val );
    void add( const String &colName, const String &val, const String &cmd );

    /// the number of rows (assumes all columns have the same number of rows)
    int pointCount() const; 

    /// load/save 
    void load( const String &fileName );
    void save( const String &fileName ) const;

    /// load/save as CSV file
    void loadCSV( const String &fileName );
    void saveCSV( const String &fileName, bool append = false ) const;

    //-------------------------------------------
    // DISPLAY / USER-INTERFACE
    //-------------------------------------------

    /// find a column that contains timestamps
    void findTimestampColumn();

    /// set default format used by new columns
    void setDefaultFormat( const String &defaultFormat ) { m_defaultFormat = defaultFormat; } 

    /// add a blank column between columns
    void addBlank( int index );

    /// the current sorting (use these indices in order)
    inline const VectorI &sortIndex() const { return m_sortIndex; }

    /// compute sort index using the specified column 
    void sort( const TableColumn &col, int sortDir );

    /// compute stats for all columns
    void computeStats();

    /// display a table using the currently register callback (takes ownership)
    static void disp( aptr<Table> table );

    /// set disp callback object (used by all table instances)
    static void setDisplay( Display< aptr<Table> > *tableDisplay ) { s_tableDisplay = tableDisplay; }

    //-------------------------------------------
    // TABLE PAGING
    //-------------------------------------------

    /// the current start index when viewing a page of table rows
    inline int startIndex() const { return m_startIndex; }

    /// the number of rows currently displayed
    inline int dispLength() const { return m_dispLength; }

    /// adjust start index by specified amount
    void step( int step );

private:

    /// retrieve column with given name; create if doesn't exist
    TableColumn &findAddColumn( const char *colName, bool &newCol );

    // table title
    String m_title;

    // the columns
    Array<TableColumn> m_columns;

    // display info
    VectorI m_sortIndex;
    int m_startIndex;
    int m_dispLength;
    int m_timestampColumn;
    String m_defaultFormat;

    // display callback (static so used by all table instances)
    static Display< aptr<Table> > *s_tableDisplay; 

    // file format tags
    enum {
        TAG_TITLE = 101,
        TAG_COLUMNS = 201 + TAG_OFFSET_ARRAY
    };

    // disable copy constructor and assignment operator
    Table( const Table &x );
    Table &operator=( const Table &x );
};


} // end namespace sbl
#endif // _SBL_TABLE_H_

