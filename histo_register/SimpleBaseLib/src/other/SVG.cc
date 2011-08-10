// Licensed under MIT license; see license.txt.

#include <sbl/other/SVG.h>
#include <sbl/core/File.h>
#include <sbl/core/StringUtil.h>
#include <sbl/math/MathUtil.h>
#include <math.h>
namespace sbl {


/// create SVG file (we open file before writing shapes into the file)
SVG::SVG( int width, int height, const String &fileName ) {
    m_file = new File( fileName, FILE_WRITE, FILE_TEXT );
    m_scale = 1;
    if (m_file->openSuccess()) {
        m_file->writeF( "<?xml version=\"1.0\"?>\n" );
        m_file->writeF( "<svg xmlns=\"http://www.w3.org/2000/svg\"\n" );
        m_file->writeF( "    xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n" );
        m_file->writeF( "    xmlns:xhtml=\"http://www.w3.org/1999/xhtml\"\n" );
        m_file->writeF( "    width=\"%d\" height=\"%d\">\n", width, height );
        m_file->writeF( "  <g>\n" );
    } else {
        delete m_file;
        m_file = NULL;
    }
}


// basic destructor
SVG::~SVG() {
    close();
}


/// close file (can no longer add shapes)
void SVG::close() {
    if (m_file) {
        m_file->writeF( "  </g>\n" );
        m_file->writeF( "</svg>\n" );
        delete m_file;
        m_file = NULL;
    }
}


/// create a group of elements with similar visual properties
void SVG::startGroup( int r, int g, int b, int rFill, int gFill, int bFill, float lineWidth ) {
    if (m_file) {
        m_file->writeF( "  </g>\n" );
        m_file->writeF( "  <g style=\"fill-opacity:1.0; stroke:#%02x%02x%02x; fill:#%02x%02x%02x; stroke-width:%4.2f;\">\n",
                            r, g, b, rFill, gFill, bFill, lineWidth );
    }
}


/// add a line to the current group
void SVG::addLine( float x1, float y1, float x2, float y2 ) {
    x1 *= m_scale;
    y1 *= m_scale;
    x2 *= m_scale;
    y2 *= m_scale;
    addObjectStr( sprintF( "    <line x1='%5.3f' y1='%5.3f' x2='%5.3f' y2='%5.3f'/>\n", x1, y1, x2, y2 ));
}


/// add a circle to the current group
void SVG::addCircle( float x, float y, float radius, const String &tooltip ) {
    x *= m_scale;
    y *= m_scale;
    radius *= m_scale;
    String tooltipText;
    if (tooltip.length())
        tooltipText = sprintF( " title='%s'", tooltip.c_str() );
    addObjectStr( sprintF( "    <circle cx='%5.3f' cy='%5.3f' r='%5.3f' %s/>\n", 
                    x, y, radius, tooltipText.c_str() ));
}


/// add text to the current group
void SVG::addText( float x, float y, const String &text, TextAnchor anchor, float size, const String &font, float rotate, float xRotate, float yRotate, float baselineShift ) {
    x *= m_scale;
    y *= m_scale;
    xRotate *= m_scale;
    yRotate *= m_scale;
    size *= m_scale;
    String extraText;
    if (anchor == ANCHOR_END)
        extraText = " text-anchor='end'";
    if (anchor == ANCHOR_MIDDLE)
        extraText = " text-anchor='middle'";
    if (rotate) 
        extraText += sprintF( " transform='rotate(%5.3f,%5.3f,%5.3f)'", rotate, xRotate, yRotate );
    if (baselineShift)
        extraText += sprintF( " dominant-baseline='middle'", baselineShift * m_scale );
    if (font.length())
        extraText += sprintF( " font-family='%s'", font.c_str() );
    addObjectStr( sprintF( "    <text x='%5.3f' y='%5.3f' font-size='%5.3f'%s>%s</text>\n",
                        x, y, size, extraText.c_str(), text.c_str() ));
}


/// add a rectangle to the current group
// fix(later): use rectangle shape
void SVG::addRectangle( float x1, float y1, float x2, float y2 ) {
    x1 *= m_scale;
    y1 *= m_scale;
    x2 *= m_scale;
    y2 *= m_scale;
    addLine( x1, y1, x1, y2 );
    addLine( x2, y1, x2, y2 );
    addLine( x1, y1, x2, y1 );
    addLine( x1, y2, x2, y2 );
}


/// add an arc to the current group
void SVG::addArc( float x, float y, float radius, float startDegrees, float endDegrees ) {
    x *= m_scale;
    y *= m_scale;
    radius *= m_scale;
    float startRadians = degToRad( startDegrees );
    float endRadians = degToRad( endDegrees );
    float xStart = x + cosf( startRadians ) * radius;
    float yStart = y + sinf( startRadians ) * radius;
    float xEnd = x + cosf( endRadians ) * radius;
    float yEnd = y + sinf( endRadians ) * radius;
    addObjectStr( sprintF( "    <path d='M%5.3f,%5.3f A%5.3f,%5.3f 0 0,1 %5.3f,%5.3f'/>", xStart, yStart, radius, radius, xEnd, yEnd ) );
}


/// utility function used to add data to the file
void SVG::addObjectStr( const String &s ) {
    if (m_file) {
        String sCopy( s );
        sCopy.replace( '\'', '"' );
        m_file->writeF( sCopy.c_str() );
    }
}


} // end namespace sbl

