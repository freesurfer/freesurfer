#ifndef _SBL_VIDEO_H_
#define _SBL_VIDEO_H_
#include <sbl/core/Pointer.h>
#include <sbl/image/Image.h>
struct CvCapture;
struct CvVideoWriter;
namespace sbl {


//-------------------------------------------
// INPUT VIDEO CLASS
//-------------------------------------------


/// The InputVideo class provides an interface to an input video file (or image sequence).
class InputVideo {
public:

    /// open a video file
    explicit InputVideo( const String &fileName = "" );
    void open( const String &fileName );

    /// close the video file
    ~InputVideo();

    /// retrieve frame by frame index
    aptr<ImageColorU> frame( int frameIndex );
    aptr<ImageGrayU> frameGray( int frameIndex );

    /// the number of frames
    int length();

    /// the size of the video frames
    int width();
    int height();

    /// return true if successfully open input video
    bool openSuccess() { return m_cvCapture || m_imageList.count(); }

private:

    // internal data
    CvCapture *m_cvCapture;
    String m_imagePath;
    Array<String> m_imageList;
    int m_width;
    int m_height;
    int m_length;

    // disable copy constructor and assignment operator
    InputVideo( const InputVideo &x );
    InputVideo &operator=( const InputVideo &x );
};


//-------------------------------------------
// OUTPUT VIDEO CLASS
//-------------------------------------------


/// The OutputVideo class provides an interface to an output video file.
class OutputVideo {
public:

    /// open the output video file
    OutputVideo();
    OutputVideo( const String &fileName, int outputWidth, int outputHeight, double fps = 30 );
    void open( const String &fileName, int outputWidth, int outputHeight, double fps = 30 );

    /// close the output video file
    ~OutputVideo();

    /// add a frame to the video
    void append( const ImageColorU &img );
    void append( const ImageGrayU &img );

    /// returns true if video file successfully created
    bool openSuccess();

    /// video image size
    inline int width() const { return m_width; }
    inline int height() const { return m_height; }

private:

    // an OpenCV video object
    CvVideoWriter *m_videoWriter;

    // video image size
    int m_width;
    int m_height;

    // disable copy constructor and assignment operator
    OutputVideo( const OutputVideo &x );
    OutputVideo &operator=( const OutputVideo &x );
};


} // end namespace sbl
#endif // _SBL_INPUT_VIDEO_H_

