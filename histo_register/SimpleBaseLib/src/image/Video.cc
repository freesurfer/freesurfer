// Licensed under MIT license; see license.txt.

#include <sbl/image/Video.h>
#include <sbl/core/StringUtil.h>
#include <sbl/math/MathUtil.h>
#include <sbl/system/FileSystem.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
#ifdef USE_OPENCV
    #include <opencv/highgui.h>
#endif
namespace sbl {


//-------------------------------------------
// INPUT VIDEO CLASS
//-------------------------------------------


/// open a video file
InputVideo::InputVideo( const String &fileName ) {
    m_cvCapture = NULL;
    m_width = 0;
    m_height = 0;
    m_length = 0;
    if (fileName.length())
        open( fileName );
}


/// open a video file
void InputVideo::open( const String &fileName ) {
    if (fileName.contains( '*' )) { 
        m_imagePath = fileName.leftOfLast( '/' ) + "/"; 
        String filePattern = fileName.rightOfLast( '/' );
        String prefix = filePattern.leftOfFirst( '*' );
        String suffix = filePattern.rightOfLast( '*' );
        m_imageList = dirFileList( m_imagePath, prefix, suffix );
        disp( 1, "image list path: %s, prefix: %s, suffix: %s, count: %d", m_imagePath.c_str(), prefix.c_str(), suffix.c_str(), m_imageList.count() );
    } else if (fileName.endsWith( ".txt" )) { 
        m_imagePath = fileName.leftOfLast( '/' ) + "/"; // assume images are in same folder as the list
        m_imageList = loadStrings( fileName, false );
        disp( 1, "image list file count: %d, path: %s", m_imageList.count(), m_imagePath.c_str() );
    } else {
#ifdef USE_OPENCV
        m_cvCapture = cvCreateFileCapture( fileName.c_str() );
#endif
    }
}


/// close the video file
InputVideo::~InputVideo() {
#ifdef USE_OPENCV
    if (m_cvCapture)
        cvReleaseCapture( &m_cvCapture );
#endif
}


/// retrieve frame by frame index
aptr<ImageColorU> InputVideo::frame( int frameIndex ) {
    if (frameIndex < 0 || frameIndex >= length())
        fatalError( "invalid input frame: %d, video length: %d", frameIndex, length() );
    aptr<ImageColorU> img;

    // get frame from OpenCV
    if (m_cvCapture) {
#ifdef USE_OPENCV
        cvSetCaptureProperty( m_cvCapture, CV_CAP_PROP_POS_FRAMES, frameIndex );
        IplImage *iplImg = cvQueryFrame( m_cvCapture );
        if (iplImg) {
            ImageColorU imgMustCopy( iplImg, false ); // need to copy this image, because it contains reference to iplImage which we do not own
            img = flipVert( imgMustCopy );
        }
#endif

    // image file
    } else if (m_imageList.count()) {
        img = load<ImageColorU>( m_imagePath + m_imageList[ frameIndex ] );
    }
    return img;
}


/// retrieve frame by frame index, as grayscale image
aptr<ImageGrayU> InputVideo::frameGray( int frameIndex ) {
    return toGray( *frame( frameIndex ) );
}


/// the number of frames
int InputVideo::length() { 
    if (m_length == 0) {
        if (m_cvCapture) {
#ifdef USE_OPENCV
            m_length = round( cvGetCaptureProperty( m_cvCapture, CV_CAP_PROP_FRAME_COUNT ) );
#endif
        } else if (m_imageList.count()) {
            m_length = m_imageList.count();
        }
    }
    return m_length;
}


/// the width of the video frames
int InputVideo::width() { 
    if (m_width == 0) {
        aptr<ImageColorU> img = frame( 0 ); // don't want to use CV_CAP_PROP_FRAME_WIDTH/HEIGHT because may not return valid value if not already loaded a frame
        if (img.get()) {
            m_width = img->width();
            m_height = img->height();
        }
    }
    return m_width;
}


/// the height of the video frames
int InputVideo::height() { 
    if (m_height == 0) {
        aptr<ImageColorU> img = frame( 0 ); // don't want to use CV_CAP_PROP_FRAME_WIDTH/HEIGHT because may not return valid value if not already loaded a frame
        if (img.get()) {
            m_width = img->width();
            m_height = img->height();
        }
    }
    return m_height;
}


//-------------------------------------------
// OUTPUT VIDEO CLASS
//-------------------------------------------


/// create unopened output video object
OutputVideo::OutputVideo() {
    m_videoWriter = NULL;
    m_width = 0;
    m_height = 0;
}


/// open the output video file
OutputVideo::OutputVideo( const String &fileName, int outputWidth, int outputHeight, double fps ) {
    m_videoWriter = NULL;
    open( fileName, outputWidth, outputHeight, fps );
}


/// open the output video file
void OutputVideo::open( const String &fileName, int outputWidth, int outputHeight, double fps ) {
#ifdef USE_OPENCV
    int fourcc = CV_FOURCC_DEFAULT;
    if (fileName.endsWith( ".avi" ))
        fourcc = CV_FOURCC( 'H', 'F', 'Y', 'U' ); // HuffYUV (lossless)
    m_videoWriter = cvCreateVideoWriter( fileName.c_str(), fourcc, fps, cvSize( outputWidth, outputHeight ) );
#endif
    m_width = outputWidth;
    m_height = outputHeight;
}


/// close the video file
OutputVideo::~OutputVideo() {
#ifdef USE_OPENCV
    if (m_videoWriter)
        cvReleaseVideoWriter( &m_videoWriter );
#endif
}


/// add a frame to the video
void OutputVideo::append( const ImageColorU &img ) {
    assertAlways( img.width() == m_width && img.height() == m_height );
    if (m_videoWriter) {
        aptr<ImageColorU> flip = flipVert( img ); // fix(asap): this shouldn't be needed!
#ifdef USE_OPENCV
        cvWriteFrame( m_videoWriter, flip->iplImage() );
#endif
    }
}


/// add a frame to the video
void OutputVideo::append( const ImageGrayU &img ) {
    aptr<ImageColorU> color = toColor( img );
    append( *color );
}


/// returns true if video file successfully created
bool OutputVideo::openSuccess() { 
    return m_videoWriter ? true : false;
}


} // end namespace sbl

