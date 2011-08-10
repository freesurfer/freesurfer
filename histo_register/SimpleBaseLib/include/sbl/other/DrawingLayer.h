#ifndef _SBL_DRAWING_LAYER_H_
#define _SBL_DRAWING_LAYER_H_
#include <sbl/core/String.h>
namespace sbl {


// The HorizAlignment specifies the horizontal alignment of text draw with a DrawingLayer.
enum HorizAlignment {
    HORIZ_ALIGN_LEFT,
    HORIZ_ALIGN_CENTER,
    HORIZ_ALIGN_RIGHT
};


/// The DrawingLayer class provides a wrapper for different drawing mediums (SVG, GUI, GL, images, etc.)
class DrawingLayer {
public:

    // avoid virtual distructor warning
    virtual ~DrawingLayer() {}

    /// draw a line between the two points
    virtual void drawLine( float x1, float y1, float x2, float y2, float r, float g, float b, float lineWidth ) = 0;

    /// draw a rectangle
    virtual void drawRect( float x, float y, float width, float height, float r, float g, float b, float lineWidth ) = 0;

    /// draw a solid-filled rectangle
    virtual void drawRectFilled( float x, float y, float width, float height, float r, float g, float b ) = 0;

    /// draw a circle of given radius
    virtual void drawCircle( float x, float y, float radius, float r, float g, float b, float lineWidth ) = 0;

    /// draw a solid-filled circle of given radius
    virtual void drawCircleFilled( float x, float y, float radius, float r, float g, float b ) = 0;

    /// draw text 
    virtual void drawText( const String &text, float x, float y, float r, float g, float b, HorizAlignment alignment = HORIZ_ALIGN_LEFT ) = 0;
};


} // end namespace sbl
#endif // _SBL_DRAWING_LAYER_H_

