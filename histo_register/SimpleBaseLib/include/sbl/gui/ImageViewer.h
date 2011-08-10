#ifndef _SBL_IMAGE_VIEWER_H_
#define _SBL_IMAGE_VIEWER_H_
#include <sbl/gui/GraphicViewer.h>
#include <sbl/image/Image.h>
namespace sbl {


/// The ImageViewer widget displays a color image, allowing zooming, panning, etc.
class ImageViewer : public GraphicViewer {
public:

    // basic constructor
    ImageViewer( wxWindow *parent, wxWindowID id );

    /// does not delete / take ownership of bitmap (must delete externally)
    void setBitmap( wxBitmap *bitmap );
    
private:

    /// draw/redraw the image
    void draw();

    // the bitmap being displayed (we own this object)
    wxBitmap *m_bitmap;
};


/// create a wxBitmap object from the image data (copies pixels)
wxBitmap *createBitmap( const ImageColorU &img );


} // end namespace sbl
#endif // _SBL_IMAGE_VIEWER_H_

