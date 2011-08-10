#ifndef _SBL_TIME_SERIES_H_
#define _SBL_TIME_SERIES_H_
#include <sbl/math/Vector.h>
#include <sbl/core/Array.h>
#include <sbl/core/File.h>
namespace sbl {


//-------------------------------------------
// TIME SERIES CLASS
//-------------------------------------------


/// The TimeSeries class represents a sequence of values, each associated with a timestamp.
/// The timestamps may be irregularly spaced.  We assume appended items have monotonically increasing timestamps.
class TimeSeries {
public:

    // basic constructors
    TimeSeries() {}
    explicit TimeSeries( const String &name ) { m_name = name; }
    explicit TimeSeries( File &file ) { load( file ); }

    /// the name of the time series
    inline const String &name() const { return m_name; }

    /// set the name of the time series
    inline void setName( const String &name ) { m_name = name; }

    /// get a time series value by item index
    inline double value( int index ) const { return m_values[ index ]; }

    /// get a timestamp by item index
    inline double timestamp( int index ) const { return m_timestamps[ index ]; }

    /// the number of items in the series
    inline int size() const { return m_values.length(); }

    /// interpolate the time series to compute a value at the given timestamp (we assume timestamps increase monotonically)
    double interpolate( double timestamp ) const;

    /// add an item to the series (we assume timestamps increase monotonically)
    inline void append( double timestamp, double value ) { m_timestamps.append( timestamp ); m_values.append( value ); }

    /// load from raw data file
    void load( File &file );

    /// save to raw data file
    void save( File &file ) const; 

    /// load from CSV file
    void loadCSV( const String &fileName );

    /// save to CSV file
    void saveCSV( const String &fileName ) const;

private:

    // the name of the time series
    String m_name;

    // the time series data
    VectorD m_values;
    VectorD m_timestamps;

    // disable copy constructor and assignment operator
    TimeSeries( const TimeSeries &x );
    TimeSeries &operator=( const TimeSeries &x );
};


//-------------------------------------------
// TIME SERIES UTILITIES
//-------------------------------------------


/// load a set of TimeSeries objects from a raw data file
void loadTimeSeriesSet( const String &fileName, Array<TimeSeries> &timeSeriesSet );


/// save a set of TimeSeries objects to a raw data file
void saveTimeSeriesSet( const String &fileName, const Array<TimeSeries> &timeSeriesSet );


} // end namespace sbl
#endif // _SBL_TIME_SERIES_H_

