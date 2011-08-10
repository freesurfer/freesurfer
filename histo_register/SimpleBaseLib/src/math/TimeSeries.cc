// Licensed under MIT license; see license.txt.

#include <sbl/math/TimeSeries.h> 
#include <sbl/math/VectorUtil.h>
namespace sbl {


//-------------------------------------------
// TIME SERIES CLASS
//-------------------------------------------


/// interpolate the time series to compute a value at the given timestamp (we assume timestamps increase monotonically)
double TimeSeries::interpolate( double timestamp ) const {
    if (m_values.length() == 0)
        return 0;
    if (timestamp <= m_timestamps[ 0 ])
        return m_values[ 0 ];
    if (timestamp >= m_timestamps.endValue())
        return m_values.endValue();
    int len = m_timestamps.size();
    double value = 0;
    for (int i = 0; i < len - 1; i++) {
        double tPrev = m_timestamps[ i ];
        double tNext = m_timestamps[ i + 1 ];
        if (timestamp >= tPrev && timestamp <= tNext) {
            double vPrev = m_values[ i ];
            double vNext = m_values[ i + 1 ];
            double tDiff = tNext - tPrev;
            if (tDiff > 1e-9) {
                double frac = (timestamp - tPrev) / tDiff;
                value = (1.0 - frac) * vPrev + frac * vNext;
            } else {
                value = vPrev;
            }
            break;
        }
    }
    return value;
}


/// load from raw data file
void TimeSeries::load( File &file ) {
    m_name = file.readString();
    m_values = file.readVector<double>();
    m_timestamps = file.readVector<double>();
}


/// save to raw data file
void TimeSeries::save( File &file ) const {
    file.writeString( m_name );
    file.writeVector( m_values );
    file.writeVector( m_timestamps );
}


/// load from CSV file
void TimeSeries::loadCSV( const String &fileName ) {
    File file( fileName, FILE_READ, FILE_TEXT );
    if (file.openSuccess()) {
        while (file.endOfFile() == false) {
            Array<String> split = file.readLine().split( "," );
            if (split.size() == 2) {
                double timestamp = split[ 0 ].toDouble();
                double value = split[ 1 ].toDouble();
                append( timestamp, value );
            }
        }
    }
}


/// save to CSV file
void TimeSeries::saveCSV( const String &fileName ) const {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess()) {
        for (int i = 0; i < m_values.size(); i++) {
            file.writeF( "%f, %f\n", m_timestamps[ i ], m_values[ i ] );
        }
    }
}


//-------------------------------------------
// TIME SERIES UTILITIES
//-------------------------------------------


/// load a set of TimeSeries objects from a raw data file
void loadTimeSeriesSet( const String &fileName, Array<TimeSeries> &timeSeriesSet ) {
    File file( fileName, FILE_READ, FILE_TEXT );
    if (file.openSuccess())
        file.readArray( timeSeriesSet );
}


/// save a set of TimeSeries objects to a raw data file
void saveTimeSeriesSet( const String &fileName, const Array<TimeSeries> &timeSeriesSet ) {
    File file( fileName, FILE_WRITE, FILE_TEXT );
    if (file.openSuccess())
        file.writeArray( timeSeriesSet );
}


} // end namespace sbl

