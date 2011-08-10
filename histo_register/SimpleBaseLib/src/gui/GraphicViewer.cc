// Licensed under MIT license; see license.txt.

#include <sbl/gui/GraphicViewer.h>
namespace sbl {


// This code uses two coordinate systems: data coords and window coords.
// Unless specified, all coordinates are data coordinates.


// basic constructor
GraphicViewer::GraphicViewer( wxWindow *parent, wxWindowID id ) : wxPanel( parent, id ) {
    Connect( wxEVT_PAINT, wxPaintEventHandler( GraphicViewer::onPaint ) );
    Connect( wxEVT_LEFT_DOWN, wxMouseEventHandler( GraphicViewer::onLeftClick ) );
    Connect( wxEVT_LEFT_DCLICK, wxMouseEventHandler( GraphicViewer::onLeftClick ) );
    Connect( wxEVT_RIGHT_DOWN, wxMouseEventHandler( GraphicViewer::onRightClick ) );
    Connect( wxEVT_RIGHT_DCLICK, wxMouseEventHandler( GraphicViewer::onRightClick ) );
    Connect( wxEVT_MOTION, wxMouseEventHandler( GraphicViewer::onMove ) );
    Connect( wxEVT_MOUSEWHEEL,  wxMouseEventHandler( GraphicViewer::onMouseWheel ) );

    // set initial viewing area
    m_xCenter = 0.5f;
    m_yCenter = 0.5f;
    m_viewWidth = 1.0f;
    m_viewHeight = 1.0f;

    // init everything else
    m_xOverlayWin = -1;
    m_yOverlayWin = -1;
    m_xDragStart = -1;
    m_yDragStart = -1;
    m_handler = NULL;

    // set widget style
    SetWindowStyle( wxBORDER_DOUBLE );
    SetBackgroundStyle( wxBG_STYLE_COLOUR );
    SetBackgroundColour( wxColor( 255, 255, 255 ) );
}


/// fit the current viewing area to the given data bounds
void GraphicViewer::viewBox( float xMin, float xMax, float yMin, float yMax ) {
    m_xCenter = (xMin + xMax) * 0.5f;
    m_yCenter = (yMin + yMax) * 0.5f;
    m_viewWidth = xMax - xMin;
    m_viewHeight = yMax - yMin;
    if (m_viewWidth == 0)
        m_viewWidth = 1;
    if (m_viewHeight == 0)
        m_viewHeight = 1;
}


// call sub-class draw function
void GraphicViewer::onPaint( wxPaintEvent &event ) {
    draw();
}


// left click: zoom in/out
void GraphicViewer::onLeftClick( wxMouseEvent &event ) {

    // get coordinate of click before scaling
    Point2 winCoord = Point2( event.GetX(), event.GetY() );
    Point2 origDataCoord = windowToDataCoord( winCoord );

    // zoom in
    if (event.ShiftDown()) {
        m_viewWidth *= 0.5f;
        m_viewHeight *= 0.5f;

    // zoom out
    } else if (event.ControlDown()) {
        m_viewWidth *= 2.0f;
        m_viewHeight *= 2.0f;
    }

    // either case: center on click
    if (event.ShiftDown() || event.ControlDown()) {
        Point2 newDataCoord = windowToDataCoord( winCoord );
        m_xCenter += origDataCoord.x - newDataCoord.x;
        m_yCenter += origDataCoord.y - newDataCoord.y;
    }

    // redraw
    Refresh();
}


// right click: pass to click handler
void GraphicViewer::onRightClick( wxMouseEvent &event ) {
    if (m_handler) {
        Point2 winCoord = Point2( event.GetX(), event.GetY() );
        Point2 dataCoord = windowToDataCoord( winCoord );
        int keyModifier = 0;
        if (event.ShiftDown())
            keyModifier |= 1;
        if (event.ControlDown())
            keyModifier |= 2;
        if (event.AltDown())
            keyModifier |= 4;
        m_handler->onClick( dataCoord.x, dataCoord.y, keyModifier );
    }
}


// left drag: reposition view based on mouse drag
void GraphicViewer::onMove( wxMouseEvent &event ) {
    if (event.LeftIsDown() && event.ShiftDown() == false && event.ControlDown() == false) {

        // convert mouse window coordinates to data coordinates
        Point2 winCoord = Point2( event.GetX(), event.GetY() );
        Point2 dataCoord = windowToDataCoord( winCoord );

        // if drag start, store data coords
        if (m_xDragStart == -1 && m_yDragStart == -1) { // fix(clean): use flag instead of special values
            m_xDragStart = dataCoord.x;
            m_yDragStart = dataCoord.y;
        }

        // adjust center so that data coords match initial drag point
        m_xCenter += m_xDragStart - dataCoord.x;
        m_yCenter += m_yDragStart - dataCoord.y;

        // redraw
        Refresh();
    } else {
        m_xDragStart = -1;
        m_yDragStart = -1;
    }
}


void GraphicViewer::onMouseWheel( wxMouseEvent &event ) {

    // get coordinate of click before scaling
    Point2 winCoord = Point2( event.GetX(), event.GetY() );
    Point2 origDataCoord = windowToDataCoord( winCoord );

    // zoom in/out
    float zoom = (float) event.GetWheelRotation() / (float) event.GetWheelDelta();
    float scale = powf( 0.8f, zoom );
    m_viewWidth *= scale;
    m_viewHeight *= scale;

    // shift center so that new click coord matches old click coord
    Point2 newDataCoord = windowToDataCoord( winCoord );
    m_xCenter += origDataCoord.x - newDataCoord.x;
    m_yCenter += origDataCoord.y - newDataCoord.y;

    // redraw
    Refresh();
}


/// convert data coordinate to window coordinate (note: we take care of y flip here)
Point2 GraphicViewer::dataToWindowCoord( Point2 dataPt ) {
    return Point2( (dataPt.x - m_xCenter + m_viewWidth * 0.5) * (double) windowWidth() / m_viewWidth,  
                   windowHeight() - (dataPt.y - m_yCenter + m_viewHeight * 0.5) * (double) windowHeight() / m_viewHeight );
}


/// convert window coordinate to data coordinate
Point2 GraphicViewer::windowToDataCoord( Point2 winPt ) const {
    return Point2( m_xCenter - m_viewWidth * 0.5 + m_viewWidth * winPt.x / (double) windowWidth(), 
                   m_yCenter + m_viewHeight * 0.5 - m_viewHeight * winPt.y / (double) windowHeight() );
}


} // end namespace sbl

