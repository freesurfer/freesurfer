// Licensed under MIT license; see license.txt.

#include <sbl/image/ImageDraw.h>
#include <sbl/image/ImageUtil.h> // for drawText
#include <math.h>
namespace sbl {


/// draws a filled rectangle on the image
void drawRectFilled( ImageColorU &img, int xMin, int xMax, int yMin, int yMax, int r, int g, int b ) {

    // check bounds
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;

    // perform draw
    for (int y = yMin; y <= yMax; y++)
        for (int x = xMin; x <= xMax; x++)
            img.setRGB( x, y, r, g, b );
}


/// draws a filled rectangle on the image
void drawRectFilled( ImageGrayU &img, int xMin, int xMax, int yMin, int yMax, int v ) {

    // check bounds
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;

    // perform draw
    for (int y = yMin; y <= yMax; y++)
        for (int x = xMin; x <= xMax; x++)
            img.data( x, y ) = v;
}


/// draws a filled rectangle on the image
void drawRectFilled( ImageGrayI &img, int xMin, int xMax, int yMin, int yMax, int v ) {

    // check bounds
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;

    // perform draw
    for (int y = yMin; y <= yMax; y++)
        for (int x = xMin; x <= xMax; x++)
            img.data( x, y ) = v;
}


/// draws a filled rectangle on the image
void drawRectFilled( ImageGrayF &img, int xMin, int xMax, int yMin, int yMax, float v ) {

    // check bounds
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;

    // perform draw
    for (int y = yMin; y <= yMax; y++)
        for (int x = xMin; x <= xMax; x++)
            img.data( x, y ) = v;
}


/// draws a non-filled rectangle on the image
void drawRect( ImageColorU &img, int xMin, int xMax, int yMin, int yMax, int r, int g, int b ) {

    // check bounds
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMax < 0) xMax = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMax < 0) yMax = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;

    // perform draw
    for (int t = xMin; t <= xMax; t++) {
        img.setRGB( t, yMin, r, g, b );
        img.setRGB( t, yMax, r, g, b );
    }
    for (int t = yMin + 1; t < yMax; t++) {
        img.setRGB( xMin, t, r, g, b );
        img.setRGB( xMax, t, r, g, b );
    }
}


/// draws a line between the two points
// fix(faster): make faster
void drawLine( ImageGrayU &img, int xStart, int yStart, int xEnd, int yEnd, int v, bool thick ) {
    float x = (float) xStart;
    float y = (float) yStart;
    float xDiff = (float) (xEnd - xStart);
    float yDiff = (float) (yEnd - yStart);
    float dist = sqrtf( xDiff * xDiff + yDiff * yDiff );
    xDiff /= dist;
    yDiff /= dist;
    int count = (int) dist;
    int width = img.width(), height = img.height();

    // loop over line
    for (int i = 0; i < count; i++) {
        int xInt = (int) x;
        int yInt = (int) y;
        if (thick) {
            if (xInt > 0 && xInt < width - 1 && yInt > 0 && yInt < height - 1) {
                img.data( xInt, yInt ) = v;
                img.data( xInt - 1, yInt ) = v;
                img.data( xInt + 1, yInt ) = v;
                img.data( xInt, yInt - 1 ) = v;
                img.data( xInt, yInt + 1 ) = v;
                img.data( xInt - 1, yInt + 1 ) = v;
                img.data( xInt + 1, yInt + 1 ) = v;
                img.data( xInt + 1, yInt - 1 ) = v;
                img.data( xInt - 1, yInt - 1 ) = v;
            }
        } else {
            if (xInt >= 0 && xInt < width && yInt >= 0 && yInt < height) 
                img.data( xInt, yInt ) = v;
        }


        // move along line
        x += xDiff;
        y += yDiff;
    }
}


/// draws a line between the two points
// fix(faster): make faster
void drawLine( ImageColorU &img, int xStart, int yStart, int xEnd, int yEnd, int r, int g, int b, bool thick ) {
    float x = (float) xStart;
    float y = (float) yStart;
    float xDiff = (float) (xEnd - xStart);
    float yDiff = (float) (yEnd - yStart);
    float dist = sqrtf( xDiff * xDiff + yDiff * yDiff );
    xDiff /= dist;
    yDiff /= dist;
    int count = (int) dist;
    int width = img.width(), height = img.height();

    // loop over line
    for (int i = 0; i < count; i++) {
        int xInt = (int) x;
        int yInt = (int) y;
        if (thick) {
            if (xInt > 0 && xInt < width - 1 && yInt > 0 && yInt < height - 1) {
                img.setRGB( xInt, yInt, r, g, b );
                img.setRGB( xInt - 1, yInt, r, g, b );
                img.setRGB( xInt + 1, yInt, r, g, b );
                img.setRGB( xInt, yInt - 1, r, g, b );
                img.setRGB( xInt, yInt + 1, r, g, b );
                img.setRGB( xInt - 1, yInt + 1, r, g, b );
                img.setRGB( xInt + 1, yInt + 1, r, g, b );
                img.setRGB( xInt + 1, yInt - 1, r, g, b );
                img.setRGB( xInt - 1, yInt - 1, r, g, b );
            }
        } else {
            if (xInt >= 0 && xInt < width && yInt >= 0 && yInt < height) 
                img.setRGB( xInt, yInt, r, g, b );
        }


        // move along line
        x += xDiff;
        y += yDiff;
    }
}


/// draws a pair of crossing lines centered on the given point
void drawCross( ImageColorU &img, int x, int y, int radius, int r, int g, int b, bool thick ) {
    drawLine( img, x - radius, y, x + radius, y, r, g, b, thick );
    drawLine( img, x, y - radius, x, y + radius, r, g, b, thick );
}


/// draws a line with an arrow head at the end
void drawArrow( ImageColorU &img, int xStart, int yStart, int xEnd, int yEnd, int r, int g, int b, bool thick ) {
    drawLine( img, xStart, yStart, xEnd, yEnd, r, g, b, thick );
    int xDiff = xEnd - xStart;
    int yDiff = yEnd - yStart;
    drawLine( img, xStart + xDiff * 4 / 5 - yDiff / 5, yStart + yDiff * 4 / 5 + xDiff / 5, xEnd, yEnd, r, g, b, thick );
    drawLine( img, xStart + xDiff * 4 / 5 + yDiff / 5, yStart + yDiff * 4 / 5 - xDiff / 5, xEnd, yEnd, r, g, b, thick );
}


/// draws a filled circle in the image
void drawCircleFilled( ImageColorU &img, int x, int y, int radius, int r, int g, int b ) {

    // check bounds
    int xMin = x - radius, xMax = x + radius;
    int yMin = y - radius, yMax = y + radius;
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;
    int rsq = radius * radius;

    // draw the circle
    for (int yc = yMin; yc <= yMax; yc++) {
        for (int xc = xMin; xc <= xMax; xc++) {
            int xd = xc - x;
            int yd = yc - y;
            if (xd * xd + yd * yd < rsq) {
                img.setRGB( xc, yc, r, g, b );
            }
        }
    }
}


/// draws a filled circle in the image
void drawCircleFilled( ImageGrayU &img, int x, int y, int radius, int v ) {

    // check bounds
    int xMin = x - radius, xMax = x + radius;
    int yMin = y - radius, yMax = y + radius;
    int width = img.width(), height = img.height();
    if (xMin < 0) xMin = 0;
    if (xMin >= width) xMin = width - 1;
    if (xMax >= width) xMax = width - 1;
    if (yMin < 0) yMin = 0;
    if (yMin >= height) yMin = height - 1;
    if (yMax >= height) yMax = height - 1;
    int rsq = radius * radius;

    // draw the circle
    for (int yc = yMin; yc <= yMax; yc++) {
        for (int xc = xMin; xc <= xMax; xc++) {
            int xd = xc - x;
            int yd = yc - y;
            if (xd * xd + yd * yd < rsq) {
                img.data( xc, yc ) = v;
            }
        }
    }
}


//-------------------------------------------
// DRAW TEXT
//-------------------------------------------


// data used by drawText
aptr<ImageGrayU> g_textImg;
const int g_letterWidth = 9, g_letterHeight = 14;


// returns the size of the specified text if it were drawn by DrawText
void textSize( const String &text, int *width, int *height ) {
    *width = text.length() * g_letterWidth;
    *height = g_letterHeight;
}


/// draws text onto the image
void drawText( ImageColorU &img, const String &text, int x, int y, int r, int g, int b, bool alignRight, int scale ) {

    // load text image if not done already
    if (g_textImg.get() == NULL) {
        g_textImg = load<ImageGrayU>( "text.bmp" );
    }

    // skip out if out of image vertically
    int width = img.width(), height = img.height();
    if (y < 0 || y + g_letterHeight * scale > height) 
        return;

    // draw the text
    if (g_textImg.get()) {
        int len = text.length();
        if (alignRight)
            x -= len * g_letterWidth * scale;
        for (int i = 0; i < len; i++) {
            int asc = text.get( i );
            int xDest = i * g_letterWidth * scale;
            if ((asc > 32 || asc < 128) && x + xDest >= 0 && x + xDest + g_letterWidth * scale <= width) {
                int xSrc = (asc - 32) * g_letterWidth;
                for (int dy = 0; dy < g_letterHeight * scale; dy++)
                    for (int dx = 0; dx < g_letterWidth * scale; dx++) {
                        if (g_textImg->data( dx / scale + xSrc, dy / scale ) == 0)
                            img.setRGB( x + dx + xDest, y + dy, r, g, b );
                    }
            }
        }
    }
}


//-------------------------------------------
// COMPUTE VISUALIZATION COLORS
//-------------------------------------------


/// maps a value in the range to an RGB color using a color table
void colorize( float v, float min, float max, int &r, int &g, int &b ) {
    colorize( (int) ((v - min) * 1000 / (max - min)), r, g, b );
}


/// maps a value in [0, 1000] to an RGB color using a color table
void colorize( int v, int &r, int &g, int &b ) {
    if (v < 0) v = 0;
    if (v > 1000) v = 1000;

    // color table
    const int count = 4;
    const int color[count][3] = {{200, 200, 200}, {0, 0, 150}, {150, 150, 0}, {200, 0, 0}};
    const int value[count] = {0, 400, 600, 1000};

    // find position in table
    int i = 0;
    for (; i < count - 1; i++) {
        if (v <= value[i + 1])
            break;
    }

    // compute color
    int lowValue = value[i];
    int highValue = value[i + 1];
    r = ((highValue - v) * (color[i][0]) + (v - lowValue) * (color[i + 1][0])) / (highValue - lowValue);
    g = ((highValue - v) * (color[i][1]) + (v - lowValue) * (color[i + 1][1])) / (highValue - lowValue);
    b = ((highValue - v) * (color[i][2]) + (v - lowValue) * (color[i + 1][2])) / (highValue - lowValue);
}


/// maps an integer to one of a set of colors
void colorizeDiscrete( int v, int &r, int &g, int &b ) {
    switch (v & 15) {
    case 0:        r = 255;    g = 0;        b = 0;        break;
    case 1:        r = 0;        g = 255;    b = 0;        break;
    case 2:        r = 0;        g = 0;        b = 255;    break;
    case 3:        r = 255;    g = 255;    b = 0;        break;
    case 4:        r = 255;    g = 0;        b = 255;    break;
    case 5:        r = 0;        g = 255;    b = 255;    break;
    case 6:        r = 255;    g = 127;    b = 0;        break;
    case 7:        r = 255;    g = 0;        b = 127;    break;
    case 8:        r = 127;    g = 255;    b = 0;        break;
    case 9:        r = 0;        g = 255;    b = 127;    break;
    case 10:    r = 127;    g = 0;        b = 255;    break;
    case 11:    r = 0;        g = 127;    b = 255;    break;
    case 12:    r = 127;    g = 0;        b = 0;        break;
    case 13:    r = 0;        g = 127;    b = 0;        break;
    case 14:    r = 0;        g = 0;        b = 127;    break;
    case 15:    r = 127;    g = 127;    b = 127;    break;
    };
}


} // end namespace sbl

