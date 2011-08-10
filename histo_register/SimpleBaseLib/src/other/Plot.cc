// Licensed under MIT license; see license.txt.

#include <sbl/other/Plot.h>
#include <sbl/core/StringUtil.h>
#include <sbl/core/Command.h>
#include <sbl/system/TimeUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/other/SVG.h>
namespace sbl {


//-------------------------------------------
// PLOT ITEM CLASS
//-------------------------------------------


// basic constructor
PlotItem::PlotItem( const VectorD &x, const VectorD &y, const Array<String> *labels, const String &caption, int r, int g, int b, PlotStyle plotStyle ) {
    m_x = x;
    m_y = y;
    m_r = r;
    m_g = g;
    m_b = b;
    m_plotStyle = plotStyle;
    m_caption = caption;
}


/// load plot data
PlotItem::PlotItem( TaggedFile &file ) {
    // fix(later): should init everything
    int tag = 0;
    do {
        tag = file.readTag();
        switch (tag) {
        case TAG_X: m_x = file.readVector<double>(); break;
        case TAG_Y: m_y = file.readVector<double>(); break;
        case TAG_LABELS: file.readStrings( m_labels ); break;
        case TAG_CAPTION: m_caption = file.readString(); break;
        case TAG_R: m_r = file.readInt(); break;
        case TAG_G: m_g = file.readInt(); break;
        case TAG_B: m_b = file.readInt(); break;
        case TAG_PLOT_STYLE: m_plotStyle = (PlotStyle) file.readInt(); break;
        default: file.skipTag( tag );
        }
    } while (tag != TAG_END_SECTION);
}


/// save plot data
void PlotItem::save( TaggedFile &file ) const {

    // save data
    file.tag( TAG_X ).writeVector( m_x );
    file.tag( TAG_Y ).writeVector( m_y );
    file.tag( TAG_LABELS ).writeStrings( m_labels );

    // save appearance
    file.tag( TAG_CAPTION ).writeString( m_caption );
    file.tag( TAG_R ).writeInt( m_r );
    file.tag( TAG_G ).writeInt( m_g );
    file.tag( TAG_B ).writeInt( m_b );
    file.tag( TAG_PLOT_STYLE ).writeInt( m_plotStyle );
}


//-------------------------------------------
// PLOT CLASS
//-------------------------------------------


// basic constructor
Plot::Plot( const String &title ) {
    m_title = title;
    m_xMin = 0;
    m_yMin = 0;
    m_xMax = 1;
    m_yMax = 1;
    m_boundSet = false;
    m_xLabelWidth = 0;
    m_xLabelHeight = 0;
    m_yLabelWidth = 0;
    m_yLabelHeight = 0;
    m_titleWidth = 0;
    m_titleHeight = 0;
    m_timeAxis = false;
    m_plotStyle = PLOT_LINES;
    m_zeroLines = true;
    m_r = 0;
    m_g = 200;
    m_b = 0;
}


/// add data to plot
void Plot::add( const VectorD &x, const VectorD &y, const Array<String> *labels, const String &caption, bool updateBounds ) {
    PlotItem *plotItem = new PlotItem( x, y, labels, caption, m_r, m_g, m_b, m_plotStyle );
    m_plotItems.append( plotItem );

    // update data bounds if requested
    if (updateBounds) {
        double xMin = x.min();
        double yMin = y.min();
        double xMax = x.max();
        double yMax = y.max();
        if ((m_boundSet == false) || xMin < m_xMin)
            m_xMin = xMin;
        if ((m_boundSet == false) || xMax > m_xMax)
            m_xMax = xMax;
        if ((m_boundSet == false) || yMin < m_yMin)
            m_yMin = yMin;
        if ((m_boundSet == false) || yMax > m_yMax)
            m_yMax = yMax;
        m_boundSet = true;
    }
}


/// add data to plot
void Plot::add( const VectorD &y, const Array<String> *labels, const String &caption, bool updateBounds ) {
    int len = y.length();
    VectorD x( len );
    for (int i = 0; i < len; i++)
        x[ i ] = i + 1.0f;
    add( x, y, labels, caption, updateBounds );
}


/// add data to plot
void Plot::add( double x, double y, bool updateBounds ) {
    VectorD xVect( 1 ), yVect( 1 );
    xVect[ 0 ] = x;
    yVect[ 0 ] = y;
    add( xVect, yVect, NULL, "", updateBounds );
}


/// add line to plot
void Plot::add( double x1, double y1, double x2, double y2, bool updateBounds ) {
    VectorD xVect( 2 ), yVect( 2 );
    xVect[ 0 ] = x1;
    yVect[ 0 ] = y1;
    xVect[ 1 ] = x2;
    yVect[ 1 ] = y2;
    add( xVect, yVect, NULL, "", updateBounds );
}


/// add columns of matrix to plot
void Plot::addColumns( const MatrixF &m ) {
    int cols = m.cols();
    for (int i = 0; i < cols; i++) 
        add( toDouble( m.col( i ) ) );
}


/// add rows of matrix to plot
void Plot::addRows( const MatrixF &m ) {
    int rows = m.rows();
    for (int i = 0; i < rows; i++) 
        add( toDouble( m.row( i ) ) );
}


/// convert data coord to image coord
void Plot::toImageCoord( double x, double y, float &xImg, float &yImg, int width, int height ) const {
    xImg = (float) ((width - 7.0 - m_yLabelWidth) * (x - m_xMin) / (m_xMax - m_xMin) + m_yLabelWidth + 1.0);
    yImg = (float) ((height - 8.0 - m_xLabelHeight - m_titleHeight) * (1.0 - (y - m_yMin) / (m_yMax - m_yMin)) + m_titleHeight + 6.0);
}


/// save to SVG file
void Plot::save( const String &fileName, int width, int height ) {

    // create image
    SVG svg( width, height, fileName );

    // check bounds
    if (m_xMin == m_xMax)
        m_xMax = m_xMin + 1;
    if (m_yMin == m_yMax)
        m_yMax = m_yMin + 1;

    // create axis labels
    String xMinLabel, xMaxLabel, yMinLabel, yMaxLabel;
    if (m_timeAxis) {
        xMinLabel = tsToStrLocal( m_xMin, true );
        xMaxLabel = tsToStrLocal( m_xMax, true );
    } else {
        xMinLabel = sprintF( "%1.4f", m_xMin );
        xMaxLabel = sprintF( "%1.4f", m_xMax );
    }
    yMinLabel = sprintF( "%1.4f", m_yMin );
    yMaxLabel = sprintF( "%1.4f", m_yMax );
    if (m_xMinLabel.length()) 
        xMinLabel = m_xMinLabel;
    if (m_xMaxLabel.length()) 
        xMaxLabel = m_xMaxLabel;
    if (m_yMinLabel.length()) 
        yMinLabel = m_yMinLabel;
    if (m_yMaxLabel.length()) 
        yMaxLabel = m_yMaxLabel;

    // fix(clean): revise or replace this
    // get axis label size
    m_xLabelWidth = 120;
    m_xLabelHeight = 20;
    m_yLabelWidth = 120;
    m_yLabelHeight = 20;
    m_titleWidth = 120;
    m_titleHeight = 20;

    // draw border
    svg.startGroup( 0, 0, 0, 0, 0, 0, 1.0f );
    svg.addRectangle( m_yLabelWidth, 5.0f + m_titleHeight, width - 5.0f, height - 1.0f - m_xLabelHeight );

    // draw zero lines
    if (m_zeroLines) {
        svg.startGroup( 220, 220, 220, 0, 0, 0, 1.0f );
        if (m_xMin < 0 && m_xMax > 0) {
            float xImg1 = 0, yImg1 = 0, xImg2 = 0, yImg2 = 0;
            toImageCoord( 0, m_yMin, xImg1, yImg1, width, height );
            toImageCoord( 0, m_yMax, xImg2, yImg2, width, height );
            svg.addLine( xImg1, yImg1, xImg2, yImg2 );
        }
        if (m_yMin < 0 && m_yMax > 0) {
            float xImg1 = 0, yImg1 = 0, xImg2 = 0, yImg2 = 0;
            toImageCoord( m_xMin, 0, xImg1, yImg1, width, height );
            toImageCoord( m_xMax, 0, xImg2, yImg2, width, height );
            svg.addLine( xImg1, yImg1, xImg2, yImg2 );
        }
    }

    // draw the data
    for (int i = 0; i < m_plotItems.count(); i++) {
        const PlotItem &plotItem = m_plotItems[ i ];

        // set visual properties for this plot item 
        int rFill = 255, gFill = 255, bFill = 255;
        if ((plotItem.style() & PLOT_CIRCLES) == 0) {
            rFill = plotItem.r(); 
            gFill = plotItem.g();
            bFill = plotItem.b();
        }
        svg.startGroup( plotItem.r(), plotItem.g(), plotItem.b(), rFill, gFill, bFill, 2.0f );

        // draw plot item
        const VectorD &x = plotItem.x();
        const VectorD &y = plotItem.y();
        assertAlways( x.length() == y.length() );
        if (plotItem.style() & PLOT_LINES) {
            for (int j = 0; j < x.length() - 1; j++) {
                float xImg1 = 0, yImg1 = 0, xImg2 = 0, yImg2 = 0;
                toImageCoord( x[ j ], y[ j ], xImg1, yImg1, width, height );
                toImageCoord( x[ j + 1 ], y[ j + 1 ], xImg2, yImg2, width, height );
                svg.addLine( xImg1, yImg1, xImg2, yImg2 );
            }
        } 
        if (plotItem.style() & (PLOT_CIRCLES | PLOT_DOTS )) {
            for (int j = 0; j < x.length(); j++) {
                float xImg = 0, yImg = 0;
                toImageCoord( x[ j ], y[ j ], xImg, yImg, width, height );
                String label;
                if (plotItem.labels().count())
                    label = plotItem.labels()[ j ];
                if (plotItem.style() & PLOT_CIRCLES)
                    svg.addCircle( xImg, yImg, 4, label );
                else
                    svg.addCircle( xImg, yImg, 1, label );
            }
        }
    }

    // draw bound captions
    svg.startGroup( 0, 0, 0, 0, 0, 0, 0.01f );
    float yLabelStart = m_yLabelHeight + 3 + m_titleHeight; 
    float yLabelEnd = height - m_xLabelHeight - 2;
    float yLabelX = m_yLabelWidth - 4;
    svg.addText( yLabelX, yLabelStart, yMaxLabel, ANCHOR_END );
    svg.addText( yLabelX, yLabelEnd, yMinLabel, ANCHOR_END );
    float xLabelStart = m_yLabelWidth + 1.0f;
    float xLabelEnd = width - 6.0f;
    float xLabelY = height - 6.0f;
    svg.addText( xLabelStart, xLabelY, xMinLabel );
    svg.addText( xLabelEnd, xLabelY, xMaxLabel, ANCHOR_END );

    // add axis labels
    if (m_yLabel.length())
        svg.addText( yLabelX, (yLabelStart + yLabelEnd) * 0.5f, m_yLabel, ANCHOR_END );
    if (m_xLabel.length())
        svg.addText( (xLabelStart + xLabelEnd) * 0.5f, xLabelY, m_xLabel, ANCHOR_MIDDLE );
    if (m_title.length())
        svg.addText( (xLabelStart + xLabelEnd) * 0.5f, m_titleHeight - 2, m_title, ANCHOR_MIDDLE );
}


/// load plot data
void Plot::loadData( const String &fileName ) {
    TaggedFile file( fileName, FILE_READ, FILE_BINARY );
    if (file.openSuccess()) {
        int tag = 0;
        do {
            tag = file.readTag();
            switch (tag) {
            case TAG_TITLE: m_title = file.readString(); break;
            case TAG_X_LABEL: m_xLabel = file.readString(); break;
            case TAG_Y_LABEL: m_yLabel = file.readString(); break;
            case TAG_X_MIN_LABEL: m_xMinLabel = file.readString(); break;
            case TAG_X_MAX_LABEL: m_xMaxLabel = file.readString(); break;
            case TAG_Y_MIN_LABEL: m_yMinLabel = file.readString(); break;
            case TAG_Y_MAX_LABEL: m_yMaxLabel = file.readString(); break;
            case TAG_X_MIN: m_xMin = file.readFloat(); break;
            case TAG_X_MAX: m_xMax = file.readFloat(); break;
            case TAG_Y_MIN: m_yMin = file.readFloat(); break;
            case TAG_Y_MAX: m_yMax = file.readFloat(); break;
            case TAG_BOUND_SET: m_boundSet = file.readBool(); break;
            case TAG_TIME_AXIS: m_timeAxis = file.readBool(); break;
            case TAG_ZERO_LINES: m_zeroLines = file.readBool(); break;
            case TAG_ITEMS: file.readArray( m_plotItems ); break;
            default: file.skipTag( tag );
            }
        } while (tag != TAG_END_SECTION);        
    }
}


/// save plot data, allowing the plot object to be reconstructed
void Plot::saveData( const String &fileName ) {
    TaggedFile file( fileName, FILE_WRITE, FILE_BINARY );
    if (file.openSuccess()) {
        file.tag( TAG_TITLE ).writeString( m_title );
        file.tag( TAG_X_LABEL ).writeString( m_xLabel );
        file.tag( TAG_Y_LABEL ).writeString( m_yLabel );
        file.tag( TAG_X_MIN_LABEL ).writeString( m_xMinLabel );
        file.tag( TAG_X_MAX_LABEL ).writeString( m_xMaxLabel );
        file.tag( TAG_Y_MIN_LABEL ).writeString( m_yMinLabel );
        file.tag( TAG_Y_MAX_LABEL ).writeString( m_yMaxLabel );
        file.tag( TAG_X_MIN ).writeDouble( m_xMin );
        file.tag( TAG_X_MAX ).writeDouble( m_xMax );
        file.tag( TAG_Y_MIN ).writeDouble( m_yMin );
        file.tag( TAG_Y_MAX ).writeDouble( m_yMax );
        file.tag( TAG_BOUND_SET ).writeBool( m_boundSet );
        file.tag( TAG_TIME_AXIS ).writeBool( m_timeAxis );
        file.tag( TAG_ZERO_LINES ).writeBool( m_zeroLines );
        file.tag( TAG_ITEMS ).writeArray( m_plotItems );
    }
}


/// display the plot (takes ownership of the plot)
void Plot::disp( aptr<Plot> plot ) { 
    if (s_plotDisplay) 
        s_plotDisplay->display( plot ); 
}


//-------------------------------------------
// HIGH-LEVEL PLOTTING FUNCTIONS
//-------------------------------------------


/// create a simple line plot of the data
aptr<Plot> simplePlot( const VectorD &v ) {
    aptr<Plot> plot( new Plot );
    plot->setStyle( PLOT_LINES );
    plot->add( v );
    return plot;
}


/// create a scatter plot of the data
aptr<Plot> scatterPlot( const VectorD &v1, const VectorD &v2, const String &name1, const String &name2 ) {
    aptr<Plot> plot( new Plot );
    plot->setStyle( PLOT_DOTS );
    plot->add( v1, v2 );
    plot->setAxisLabels( name1, name2 );
    return plot;
}


/// create a histogram of the data
aptr<Plot> histogramPlot( const VectorD &v, int bucketCount ) {
    double min = v.min();
    double max = v.max();
    int len = v.length();

    // if not specified, compute bucket count
    if (bucketCount < 1) {
        bucketCount = len / 10;
        if (bucketCount < 5)
            bucketCount = 5;
        if (bucketCount > 1000)
            bucketCount = 1000;
    }

    // create histogram
    VectorD hist( bucketCount );
    hist.clear( 0 );
    for (int i = 0; i < len; i++) {
        double val = v[ i ];
        val = (val - min) / (max - min);
        int index = (int) (val * (double) (bucketCount - 1));
        if (index < 0) index = 0;
        if (index > bucketCount - 1) index = bucketCount - 1;
        hist[ index ]++;
    }

    // create range vector
    VectorD range( bucketCount );
    for (int i = 0; i < bucketCount; i++) {
        double frac = (double) i / (double) (bucketCount - 1);
        range[ i ] = min + (max - min) * frac;
    }

    // plot the histogram
    aptr<Plot> plot( new Plot() );
    plot->setStyle( PLOT_LINES );
    plot->add( range, hist );
    return plot;
}


// init plot disp callback object (used by all plot instances)
Display< aptr<Plot> > *Plot::s_plotDisplay = NULL;


} // end namespace sbl

