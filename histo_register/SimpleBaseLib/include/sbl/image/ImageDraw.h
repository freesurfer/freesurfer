#ifndef _SBL_IMAGE_DRAW_H_
#define _SBL_IMAGE_DRAW_H_
#include <sbl/image/Image.h>
namespace sbl {


/*! \file ImageDraw.h
    \brief The ImageDraw module includes functions for drawing shapes, etc. on images.
*/


/// draws a filled rectangle on the image
// fix(clean): use templates to reduce number of func versions
void drawRectFilled( ImageColorU &img, int xMin, int xMax, int yMin, int yMax, int r, int g, int b );
void drawRectFilled( ImageGrayU &img, int xMin, int xMax, int yMin, int yMax, int v );
void drawRectFilled( ImageGrayI &img, int xMin, int xMax, int yMin, int yMax, int v );
void drawRectFilled( ImageGrayF &img, int xMin, int xMax, int yMin, int yMax, float v );


/// draws a non-filled rectangle on the image
void drawRect( ImageColorU &img, int xMin, int xMax, int yMin, int yMax, int r, int g, int b );


/// draws a line between the two points
void drawLine( ImageGrayU &img, int xStart, int yStart, int xEnd, int yEnd, int v, bool thick );
void drawLine( ImageColorU &img, int xStart, int yStart, int xEnd, int yEnd, int r, int g, int b, bool thick );


/// draws a pair of crossing lines centered on the given point
void drawCross( ImageColorU &img, int x, int y, int radius, int r, int g, int b, bool thick );


/// draws a line with an arrow head at the end
void drawArrow( ImageColorU &img, int xStart, int yStart, int xEnd, int yEnd, int r, int g, int b, bool thick );


/// draws a filled circle in the image
void drawCircleFilled( ImageColorU &img, int x, int y, int radius, int r, int g, int b );
void drawCircleFilled( ImageGrayU &img, int x, int y, int radius, int v );


/// draws text onto the image
void drawText( ImageColorU &img, const String &text, int x, int y, int r, int g, int b, bool alignRight = false, int scale = 1 );


/// maps a value in the range to an RGB color using a color table
void colorize( float v, float min, float max, int &r, int &g, int &b );


/// maps a value in [0, 1000] to an RGB color using a color table
void colorize( int v, int &r, int &g, int &b );


/// maps an integer to one of a set of colors
void colorizeDiscrete( int v, int &r, int &g, int &b );


} // end namespace sbl
#endif // _SBL_IMAGE_DRAW_H_

