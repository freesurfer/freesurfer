// Licensed under MIT license; see license.txt.

#include <sbl/other/ImageDrawingLayer.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageDraw.h> // fix(clean): eliminate other interface and only use this?
namespace sbl {


/// draw a line between two points
void ImageDrawingLayer::drawLine( float x1, float y1, float x2, float y2, float r, float g, float b, float lineWidth ) {
    sbl::drawLine( m_image, round( x1 ), round( y1 ), round( x2 ), round( y2 ), round( r ), round( g ), round( b ), lineWidth > 1 );
}


/// draw a rectangle
void ImageDrawingLayer::drawRect( float x, float y, float width, float height, float r, float g, float b, float lineWidth ) {
    int xMin = round( x );
    int xMax = round( x + width );
    int yMin = round( y );
    int yMax = round( y + height );
    sbl::drawRect( m_image, xMin, xMax, yMin, yMax, round( r ), round( g ), round( b ) );
}


    /// draw a solid-filled rectangle
void ImageDrawingLayer::drawRectFilled( float x, float y, float width, float height, float r, float g, float b ) {
    int xMin = round( x );
    int xMax = round( x + width );
    int yMin = round( y );
    int yMax = round( y + height );
    sbl::drawRectFilled( m_image, xMin, xMax, yMin, yMax, round( r ), round( g ), round( b ) );
}


/// draw a circle of given radius
void ImageDrawingLayer::drawCircle( float x, float y, float radius, float r, float g, float b, float lineWidth ) {
    fatalError( "not implemented" );
    // fix(later): add this
}


/// draw a solid-filled circle of given radius
void ImageDrawingLayer::drawCircleFilled( float x, float y, float radius, float r, float g, float b ) {
    sbl::drawCircleFilled( m_image, round( x ), round( y ), round( radius ), round( r ), round( g ), round( b ) );
}


/// draw text
void ImageDrawingLayer::drawText( const String &text, float x, float y, float r, float g, float b, HorizAlignment alignment ) {
    sbl::drawText( m_image, text, round( x ), round( y ), round( r ), round( g ), round( b ), alignment == HORIZ_ALIGN_RIGHT );
}


} // end namespace sbl

