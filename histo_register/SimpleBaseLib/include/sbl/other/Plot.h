#ifndef _SBL_PLOT_H_
#define _SBL_PLOT_H_
#include <sbl/core/Array.h>
#include <sbl/core/String.h>
#include <sbl/core/Pointer.h>
#include <sbl/math/Vector.h>
#include <sbl/math/Matrix.h>
#include <sbl/other/TaggedFile.h>
#include <sbl/other/SVG.h>
namespace sbl {


/*! \file Plot.h
    \brief The Plot module is used for generating plots of data, including line plots, scatter plots, and histograms.
    The plots can be saved as SVG files.
*/


/// drawing style of plot items
enum PlotStyle {
    PLOT_LINES = 1,
    PLOT_DOTS = 2,
    PLOT_CIRCLES = 4
};


/// The PlotItem class represents a set of data with a uniform appearance (a line or set of dots of a given color).
class PlotItem {
public:

    // basic constructor
    PlotItem( const VectorD &x, const VectorD &y, const Array<String> *labels, 
              const String &caption, int r, int g, int b, PlotStyle plotStyle );

    /// load plot data
    PlotItem( TaggedFile &file );

    /// save plot data
    void save( TaggedFile &file ) const;

    /// the lot data
    inline const VectorD &x() const { return m_x; } 
    inline const VectorD &y() const { return m_y; } 

    /// the labels assigned to the data points (if any)
    inline const Array<String> &labels() const { return m_labels; }

    /// display color
    inline int r() const { return m_r; }
    inline int g() const { return m_g; }
    inline int b() const { return m_b; }

    /// display style (lines, dots, etc.)
    inline PlotStyle style() const { return m_plotStyle; }

private:

    // the data to plot
    VectorD m_x;
    VectorD m_y;
    Array<String> m_labels; 

    // appearance: color/style
    String m_caption;
    int m_r;
    int m_g;
    int m_b;
    PlotStyle m_plotStyle;

    // file format tags
    enum {
        TAG_X = 1101,
        TAG_Y = 1102,
        TAG_LABELS = 1103,
        TAG_CAPTION = 1201, 
        TAG_R = 1301,
        TAG_G = 1302,
        TAG_B = 1303,
        TAG_PLOT_STYLE = 1304
    };

    // disable copy constructor and assignment operator
    PlotItem( const PlotItem &x );
    PlotItem &operator=( const PlotItem &x );
};


/// The Plot class can be used to generate a 2D graphical plot from a set of data.
class Plot {
public:

    // basic constructor 
    explicit Plot( const String &title = "" );

    /// plot data
    inline int itemCount() const { return m_plotItems.count(); }
    inline const PlotItem &item( int index ) const { return m_plotItems[ index ]; }

    /// the visible data bounds
    inline double xMin() const { return m_xMin; }
    inline double xMax() const { return m_xMax; }
    inline double yMin() const { return m_yMin; }
    inline double yMax() const { return m_yMax; }

    /// add data to plot
    void add( const VectorD &x, const VectorD &y, const Array<String> *labels = NULL, const String &caption = "", bool updateBounds = true );
    void add( const VectorD &y, const Array<String> *labels = NULL, const String &caption = "", bool updateBounds = true );
    void add( double x, double y, bool updateBounds = true );
    void add( double x1, double y1, double x2, double y2, bool updateBounds = true );
    void addColumns( const MatrixF &m );
    void addRows( const MatrixF &m );
    void addVertLine( double x ) { add( x, m_yMin, x, m_yMax, false ); }
    void addHorizLine( double y ) { add( m_xMin, y, m_xMax, y, false ); }

    /// add caption to key (using current color)
    void addKeyCaption( const String &caption );

    /// set item appearance (used for next item added)
    inline void setColor( int r, int g, int b ) { m_r = r; m_g = g; m_b = b; }
    inline void setStyle( int plotStyle ) { m_plotStyle = (PlotStyle) plotStyle; }

    /// set plot appearance
    inline void setTitle( const String &title ) { m_title = title; }
    inline void setAxisLabels( const String &xLabel, const String &yLabel ) { m_xLabel = xLabel; m_yLabel = yLabel; }
    inline void setXBoundLabels( const String &xMinLabel, const String &xMaxLabel ) { m_xMinLabel = xMinLabel; m_xMaxLabel = xMaxLabel; }
    inline void setYBoundLabels( const String &yMinLabel, const String &yMaxLabel ) { m_yMinLabel = yMinLabel; m_yMaxLabel = yMaxLabel; }
    inline void setTimeAxis( bool timeAxis ) { m_timeAxis = timeAxis; }
    inline void setZeroLines( bool zeroLines ) { m_zeroLines = zeroLines; }
    inline void setXBounds( double xMin, double xMax ) { m_xMin = xMin; m_xMax = xMax; }
    inline void setYBounds( double yMin, double yMax ) { m_yMin = yMin; m_yMax = yMax; }

    /// save to SVG file
    void save( const String &fileName, int width = 950, int height = 600 );

    /// load/save plot data, allowing the plot object to be reconstructed
    void loadData( const String &fileName );
    void saveData( const String &fileName );

    /// display the plot (takes ownership of the plot)
    static void disp( aptr<Plot> plot );

    /// set disp callback object (used by all plot instances)
    static void setDisplay( Display< aptr<Plot> > *plotDisplay ) { s_plotDisplay = plotDisplay; }

private:

    // plot data
    Array<PlotItem> m_plotItems;

    // plot formatting info
    double m_xMin;
    double m_xMax;
    double m_yMin;
    double m_yMax;
    bool m_boundSet;
    float m_xLabelWidth;
    float m_xLabelHeight;
    float m_yLabelWidth;
    float m_yLabelHeight;
    float m_titleWidth;
    float m_titleHeight;
    bool m_timeAxis;
    bool m_zeroLines;

    // labels and titles
    String m_title;
    String m_xLabel;
    String m_yLabel;
    String m_xMinLabel;
    String m_xMaxLabel;
    String m_yMinLabel;
    String m_yMaxLabel;

    // item appearance (used for next item added)
    PlotStyle m_plotStyle;
    int m_r;
    int m_g;
    int m_b;

    // convert data coord to image coord
    void toImageCoord( double x, double y, float &xImg, float &yImg, int width, int height ) const;

    // display callback (static so used by all plot instances)
    static Display< aptr<Plot> > *s_plotDisplay; 

    // file format tags
    enum {
        TAG_TITLE = 101,
        TAG_X_LABEL = 102,
        TAG_Y_LABEL = 103,
        TAG_X_MIN_LABEL = 104,
        TAG_X_MAX_LABEL = 105,
        TAG_Y_MIN_LABEL = 106,
        TAG_Y_MAX_LABEL = 107,
        TAG_X_MIN = 301,
        TAG_X_MAX = 302,
        TAG_Y_MIN = 303,
        TAG_Y_MAX = 304,
        TAG_BOUND_SET = 305,
        TAG_TIME_AXIS = 401,
        TAG_ZERO_LINES = 402,
        TAG_ITEMS = 1000
    };

    // disable copy constructor and assignment operator
    Plot( const Plot &x );
    Plot &operator=( const Plot &x );
};


/// create a simple line plot of the data
aptr<Plot> simplePlot( const VectorD &v );


/// create a scatter plot of the data
aptr<Plot> scatterPlot( const VectorD &v1, const VectorD &v2, const String &name1, const String &name2 );


/// create a histogram of the data
aptr<Plot> histogramPlot( const VectorD &v, int bucketCount = -1 );


} // end namespace sbl
#endif // _SBL_PLOT_H_

