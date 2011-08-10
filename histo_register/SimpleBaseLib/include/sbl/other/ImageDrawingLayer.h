#ifndef _SBL_IMAGE_DRAWING_LAYER_H_
#define _SBL_IMAGE_DRAWING_LAYER_H_
#include <sbl/image/Image.h>
#include <sbl/other/DrawingLayer.h>
namespace sbl {


/// The ImageDrawingLayer allows drawing on a color image.
class ImageDrawingLayer : public DrawingLayer {
public:

    // basic constructor (image is owned outside this class)
    ImageDrawingLayer( ImageColorU &image ) : m_image( image ) {}
    virtual ~ImageDrawingLayer() {}

    /// draw a line between two points
    void drawLine( float x1, float y1, float x2, float y2, float r, float g, float b, float lineWidth );

    /// draw a rectangle
    void drawRect( float x, float y, float width, float height, float r, float g, float b, float lineWidth );

    /// draw a solid-filled rectangle
    void drawRectFilled( float x, float y, float width, float height, float r, float g, float b );

    /// draw a circle of given radius
    void drawCircle( float x, float y, float radius, float r, float g, float b, float lineWidth );

    /// draw a solid-filled circle of given radius
    void drawCircleFilled( float x, float y, float radius, float r, float g, float b );

    /// draw text
    void drawText( const String &text, float x, float y, float r, float g, float b, HorizAlignment alignment = HORIZ_ALIGN_LEFT );

private:

    // the image on which we are drawing
    ImageColorU &m_image;
};

    
} // end namespace sbl
#endif // _SBL_IMAGE_DRAWING_LAYER_H_

