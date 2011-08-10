// Licensed under MIT license; see license.txt.

#include <sbl/gui/ImageSeqViewer.h>
#include <sbl/core/StringUtil.h>
#include <sbl/core/Command.h> // for test command
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/gui/MainWindow.h>
namespace sbl {


//-------------------------------------------
// IMAGE SEQUENCE VIEWER CLASS
//-------------------------------------------


// the current instance (if requested set in constructor)
ImageSeqViewer *ImageSeqViewer::s_instance = NULL;


// basic constructor
ImageSeqViewer::ImageSeqViewer( wxWindow *parent, wxWindowID id, bool storeInstance ) : wxPanel( parent, id ) {
    if (storeInstance)
        s_instance = this;
    m_currentImageIndex = 0;
    m_bitmap = NULL;
    m_drawCallback = NULL;
    m_clickCallback = NULL;
    m_imageCount = 10;

    // add vertical sizer
    wxBoxSizer *sizer = new wxBoxSizer( wxVERTICAL );
    SetSizer( sizer );

    // add image viewer
    m_imageViewer = new ImageViewer( this, -1 );
    sizer->Add( m_imageViewer, 1, wxEXPAND );
    m_imageViewer->setHandler( this );

    // add slider
    m_slider = new wxSlider( this, 1, 0, 0, 99 );
    sizer->Add( m_slider, 0, wxEXPAND );
    Connect( 1, wxEVT_COMMAND_SLIDER_UPDATED, wxCommandEventHandler( ImageSeqViewer::onSlider ) );

    // add status text area
    m_statusText = new wxStaticText( this, -1, "" );
    sizer->Add( m_statusText, 0, wxEXPAND );
}


// basic destructor
ImageSeqViewer::~ImageSeqViewer() {
    if (m_bitmap)
        delete m_bitmap;
}


/// display an image with optional label
void ImageSeqViewer::dispImage( const ImageColorU &image, const String &label ) {
    String status = sprintF( "%d / %d", m_currentImageIndex, m_imageCount );
    if (label.length()) {
        status += ": ";
        status += label;
    }
    m_statusText->SetLabel( status.c_str() );
    if (m_bitmap) delete m_bitmap;
    m_bitmap = createBitmap( image );
    m_imageViewer->setBitmap( m_bitmap );
}


/// set the index of the image currently being displayed
void ImageSeqViewer::setImageIndex( int imageIndex ) {
    imageIndex = bound( imageIndex, 0, m_imageCount - 1 );
    m_currentImageIndex = imageIndex;
    m_slider->SetValue( imageIndex );
}


/// set the number of images to display
void ImageSeqViewer::setImageCount( int imageCount ) { 
    m_imageCount = imageCount;
    m_slider->SetMax( imageCount );
}


// handle click on image
void ImageSeqViewer::onClick( double x, double y, int keyModifier ) { 
    if (m_clickCallback) {
        m_clickCallback( m_currentImageIndex, x, y, keyModifier ); 
        if (m_drawCallback) {
            m_drawCallback( m_currentImageIndex );
        }
    }
}


// handle slider value change
void ImageSeqViewer::onSlider( wxCommandEvent &event ) {
    int imageIndex = event.GetInt();
    m_currentImageIndex = bound( imageIndex, 0, m_imageCount - 1 );
    if (m_drawCallback)
        m_drawCallback( m_currentImageIndex );
}


// handle mouse wheel movement; redirect to image viewer
void ImageSeqViewer::onMouseWheel( wxMouseEvent &event ) {
    m_imageViewer->onMouseWheel( event );
}


//-------------------------------------------
// TEST IMAGE SEQUENCE VIEWER
//-------------------------------------------


// a collection of images used to test the image sequence viewer
Array<ImageColorU> g_testImages;


// draw one of the test images
void testImageDrawCallback( int imageIndex ) {
    String label = sprintF( "label-%d", imageIndex );
    dispImage( g_testImages[ imageIndex ], label );
}


// make a white dot on one of the test images
void testImageClickCallback( int imageIndex, double x, double y, int keyModifier ) {
    int xInt = round( x );
    int yInt = round( y );
    ImageColorU &image = g_testImages[ imageIndex ];
    if (image.inBounds( xInt, yInt ))
        image.setRGB( xInt, yInt, 255, 255, 255 );
}


// command to set up the image seq viewer test
void testImageSeqViewer( Config &conf ) {
    int count = 50;
    int width = 200, height = 200;
    for (int i = 0; i < count; i++) {
        ImageColorU *image = new ImageColorU( width, height );
        image->clear( 0 );
        for (int y = 1; y < height - 1; y++) {
            for (int x = 1; x < width - 1; x++) {
                image->setRGB( x, y, x, y, i * 4 );
            }
        }
        g_testImages.append( image );
    }
    ImageSeqViewer::instance()->setClickCallback( testImageClickCallback );
    ImageSeqViewer::instance()->setDrawCallback( testImageDrawCallback );
    ImageSeqViewer::instance()->setImageCount( count );
    testImageDrawCallback( 0 );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initImageSeqViewer() {
//#ifdef REGISTER_TEST_COMMANDS
    registerCommand( "testiseq", testImageSeqViewer );
//#endif
}


} // end namespace sbl

