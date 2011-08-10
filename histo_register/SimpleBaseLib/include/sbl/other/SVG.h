#ifndef _SBL_SVG_H_
#define _SBL_SVG_H_
#include <sbl/core/File.h>
#include <sbl/core/String.h>
namespace sbl {


/// horizontal text alignment
enum TextAnchor {
    ANCHOR_START,
    ANCHOR_MIDDLE,
    ANCHOR_END
};


/// The SVG class is used to create SVG (scalable vector graphics) files.
class SVG {
public:

    /// create SVG file (we open file before writing shapes into the file)
    SVG( int width, int height, const String &fileName );

    // basic destructor
    ~SVG();

    /// close file (can no longer add shapes)
    void close();

    /// a scale factor applied to all subsequently added shapes
    inline void setScale( float scale ) { m_scale = scale; }

    /// create a group of elements with similar visual properties
    void startGroup( int r, int g, int b, int rFill, int gFill, int bFill, float lineWidth );

    /// draw shapes
    void addRectangle( float x1, float y1, float x2, float y2 );    
    void addLine( float x1, float y1, float x2, float y2 );
    void addCircle( float x, float y, float radius, const String &tooltip = "" );
    void addText( float x, float y, const String &text, TextAnchor anchor = ANCHOR_START, float size = 16, const String &font = "", float rotate = 0, float xRotate = 0, float yRotate = 0, float baselineShift = 0 );
    void addArc( float x, float y, float radius, float startDegrees, float endDegrees );
    
private:

    /// utility function used to add data to the file
    void addObjectStr( const String &s );

    // a scale factor applied to all subsequently added shapes
    float m_scale;

    // the currently open file; appends to this file as shapes/etc. added to SVG object
    File *m_file;

    // disable copy constructor and assignment operator
    SVG( const SVG &x );
    SVG &operator=( const SVG &x );
};


} // end namespace sbl
#endif // _SBL_SVG_H_

