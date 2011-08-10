// Licensed under MIT license; see license.txt.

#include <sbl/gui/ImageViewer.h>
#include <sbl/math/MathUtil.h>
#include <wx/graphics.h>
namespace sbl {


// basic constructor
ImageViewer::ImageViewer( wxWindow *parent, wxWindowID id ) : GraphicViewer( parent, id ) {
    m_bitmap = NULL;
    SetDoubleBuffered( true );
}


// does not delete or take ownership of bitmap (must delete externally)
void ImageViewer::setBitmap( wxBitmap *bitmap ) { 

    // first time through, set view to full image (expanding to preserve aspect ratio)
    if (m_bitmap == NULL) {
        double windowAspect = (double) windowWidth() / (double) windowHeight();
        double imageAspect = (double) bitmap->GetWidth() / (double) bitmap->GetHeight();
        if (windowAspect > imageAspect) {
            int viewWidth = round( bitmap->GetHeight() * windowAspect );
            viewBox( 0, viewWidth - 1, 0, bitmap->GetHeight() - 1 );
        } else if (windowAspect) {
            int viewHeight = round( bitmap->GetHeight() / windowAspect );
            viewBox( 0, bitmap->GetWidth() - 1, 0, viewHeight - 1 );
        } else {
            viewBox( 0, bitmap->GetWidth() - 1, 0, bitmap->GetHeight() - 1 );
        }
    }

    // store the bitmap pointer (we'll assume it remains valid until the next setBitmap call)
    m_bitmap = bitmap;
    Refresh();
}


/// draw/redraw the image
void ImageViewer::draw() {
    wxPaintDC dc( this );
    wxGraphicsContext *gc = wxGraphicsContext::Create( dc );
    if (gc && m_bitmap) {
        double xScale = (double) windowWidth() / (double) viewWidth();
        double yScale = (double) windowHeight() / (double) viewHeight();
        double xOffset = -xCenter() + (double) viewWidth() * 0.5;
        double yOffset = -yCenter() - (double) viewHeight() * 0.5f;
        gc->Scale( xScale, -yScale );
        gc->Translate( xOffset, yOffset );
        gc->DrawBitmap( *m_bitmap, 0, 0, m_bitmap->GetWidth(), m_bitmap->GetHeight() );
    }
}


/// create a wxBitmap object from the image data (copies pixels)
wxBitmap *createBitmap( const ImageColorU &img ) {
    int width = img.width(), height = img.height();
    wxImage image( width, height );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int r = img.r( x, y );
            int g = img.g( x, y );
            int b = img.b( x, y );
            image.SetRGB( x, y, r, g, b );
        }
    }
    return new wxBitmap( image );
}


} // end namespace sbl

