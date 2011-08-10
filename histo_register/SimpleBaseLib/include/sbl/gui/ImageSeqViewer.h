#ifndef _SBL_IMAGE_SEQ_VIEWER_H_
#define _SBL_IMAGE_SEQ_VIEWER_H_
#include <sbl/gui/ImageViewer.h>
#include <sbl/image/Image.h>
namespace sbl {


// register commands, etc. defined in this module
void initImageSeqViewer();


/// The ImageSeqViewer provides a user interface for viewing a sequence of images.
/// The viewer doesn't actually contain the data being displayed; that is stored and controlled externally.
class ImageSeqViewer : public wxPanel, public GraphicViewerHandler {
public:

    // basic constructor
    ImageSeqViewer( wxWindow *parent, wxWindowID id, bool storeInstance );
    ~ImageSeqViewer();

    /// display an image with optional label
    void dispImage( const ImageColorU &image, const String &label );

    /// the index of the image currently being displayed
    inline int currentImageIndex() const { return m_currentImageIndex; }

    /// set the index of the image currently being displayed
    void setImageIndex( int imageIndex );

    /// set the number of images to display
    void setImageCount( int imageCount );

    /// set a callback to call when we need to draw an image
    inline void setDrawCallback( void (*drawCallback)( int imageIndex ) ) { m_drawCallback = drawCallback; }

    /// a callback to call when the user clicks on the image
    inline void setClickCallback( void (*clickCallback)( int imageIndex, double x, double y, int keyModifier ) ) { m_clickCallback = clickCallback; }

    // handle graphic viewer events
    void onClick( double x, double y, int keyModifier );

    // handlers for wxWidgets events
    void onSlider( wxCommandEvent &event );
    void onMouseWheel( wxMouseEvent &event );

    /// the current instance (if requested set in constructor)
    inline static ImageSeqViewer *instance() { return s_instance; }

private:

    // the currently displayed image index
    int m_currentImageIndex;

    // the total number of images to display
    int m_imageCount;

    // internal gui elements
    ImageViewer *m_imageViewer;
    wxSlider *m_slider;
    wxStaticText *m_statusText;
    wxBitmap *m_bitmap;

    // the current instance (if requested set in constructor)
    static ImageSeqViewer *s_instance;

    // a callback to call when we need to draw an image
    void (*m_drawCallback)( int imageIndex );

    // a callback to call when the user clicks on the image
    void (*m_clickCallback)( int imageIndex, double x, double y, int keyModifier );
};


} // end namespace sbl
#endif // _SBL_IMAGE_SEQ_VIEWER_H_

