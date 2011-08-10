#ifndef _SBL_GRAPHIC_VIEWER_H_
#define _SBL_GRAPHIC_VIEWER_H_
#include <wx/wx.h>
#include <sbl/core/String.h>
#include <sbl/math/Geometry.h>
namespace sbl {


/// The GraphicViewerHandler class defines an interface for objects that wish to be notified of clicks within a GraphicViewer.
class GraphicViewerHandler {
public:

    // called on a user mouse click in a GraphicViewer
    virtual void onClick( double x, double y, int keyModifier ) = 0;
};


/// The GraphicViewer widget is a sub-class for other widgets that display graphics, allowing zooming, panning, etc.
class GraphicViewer : public wxPanel {
public:

    // basic constructor
    GraphicViewer( wxWindow *parent, wxWindowID id );

    /// the current viewing center
    inline float xCenter() const { return m_xCenter; }
    inline float yCenter() const { return m_yCenter; }

    /// the current viewing area size
    inline float viewWidth() const { return m_viewWidth; }
    inline float viewHeight() const { return m_viewHeight; }

    /// fit the current viewing area to the given data bounds
    void viewBox( float xMin, float xMax, float yMin, float yMax );

    /// the size of this viewer widget
    inline int windowWidth() const { return GetSize().GetWidth() - GetWindowBorderSize().GetX(); }
    inline int windowHeight() const { return GetSize().GetHeight() - GetWindowBorderSize().GetY(); }

    /// set the object that will receive mouse event notifications
    inline void setHandler( GraphicViewerHandler *handler ) { m_handler = handler; }

    /// draw/redraw the graphic; implemented by subclass
    virtual void draw() = 0;

    // wxWidgets event handlers
    void onPaint( wxPaintEvent &event );
    void onLeftClick( wxMouseEvent &event );
    void onRightClick( wxMouseEvent &event );
    void onMove( wxMouseEvent &event );
    void onMouseWheel( wxMouseEvent &event );

protected:

    /// convert data coordinate to window coordinate (note: we take care of y flip here)
    Point2 dataToWindowCoord( Point2 dataPt );

    /// convert window coordinate to image coordinate
    Point2 windowToDataCoord( Point2 winPt ) const;

private:

    // current viewing area
    float m_xCenter;
    float m_yCenter;
    float m_viewWidth;
    float m_viewHeight;

    // the center position when dragging started
    float m_xDragStart;
    float m_yDragStart;

    // the window coordinates of the mouse for displaying the overlay
    int m_xOverlayWin;
    int m_yOverlayWin;

    // the object currently registered to handle events from this viewer
    GraphicViewerHandler *m_handler;
};


} // end namespace sbl
#endif // _SBL_GRAPHIC_VIEWER_H_

