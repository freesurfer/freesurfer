#ifndef _SBL_IMAGE_H_
#define _SBL_IMAGE_H_
#include <sbl/core/String.h>
#ifdef USE_OPENCV
    #include <opencv2/core/types_c.h>
    CVAPI(void) cvReleaseImage( IplImage** image ); // we need this, but cxcore seems to be incompatible with something else (windows?)
#else
    #include <string.h> // for memcpy
#endif
#include <typeinfo>
namespace sbl {


// static assert; used here to make sure we don't use color methods on gray images or vice versa
#define assertStatic(condition) switch(0){case 0:case condition:;}


// debug mode image bound checking
#define IMGCHK assertDebug( inBounds( x, y ) )
#define IMGCHK_INTERP assertDebug( inBounds( x, y ) )


// note: we use BGR order for RGB images
#define R_CHANNEL 2
#define G_CHANNEL 1
#define B_CHANNEL 0


//-------------------------------------------
// IMAGE CLASS 
//-------------------------------------------


/// The ImageColor class represents a multi-channel (e.g. RGB) image.
template <typename T, int CHANNEL_COUNT > class Image {
public:

    /// basic constructor; does not initialize image
    Image( int width, int height ) { alloc( width, height ); }

    // basic copy constructor
    Image( const Image &img ) {    alloc( img.width(), img.height() );    memcpy( m_raw, img.rawConst(), m_rowBytes * m_height ); }

    // basic destructor
    ~Image();

    /// get/set pixel values for gray images (assumes CHANNEL_COUNT == 1)
    inline T &data( int x, int y ) { assertStatic( CHANNEL_COUNT == 1 ); IMGCHK return m_ptr[ y ][ x ]; }
    inline const T &data( int x, int y ) const { assertStatic( CHANNEL_COUNT == 1 ); IMGCHK return m_ptr[ y ][ x ]; }
    float interp( float x, float y ) const;

    /// get/set pixel values for multi-channel images
    inline T &data( int x, int y, int c ) { IMGCHK return m_ptr[ y ][ x * CHANNEL_COUNT + c ]; }
    inline const T &data( int x, int y, int c ) const { IMGCHK return m_ptr[ y ][ x * CHANNEL_COUNT + c ]; }
    float interp( float x, float y, int c ) const;

    /// get/set pixel values for RGB (actually BGR) color images (assumes CHANNEL_COUNT >= 3)
    inline T b( int x, int y ) const { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK return m_ptr[ y ][ x * CHANNEL_COUNT + B_CHANNEL ]; }
    inline T g( int x, int y ) const { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK return m_ptr[ y ][ x * CHANNEL_COUNT + G_CHANNEL ]; }
    inline T r( int x, int y ) const { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK return m_ptr[ y ][ x * CHANNEL_COUNT + R_CHANNEL ]; }
    inline void setRGB( int x, int y, T r, T g, T b ) { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK m_ptr[ y ][ x * CHANNEL_COUNT + B_CHANNEL ] = b; m_ptr[ y ][ x * CHANNEL_COUNT + G_CHANNEL ] = g; m_ptr[ y ][ x * CHANNEL_COUNT + R_CHANNEL ] = r; }
    inline void setB( int x, int y, T v ) { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK m_ptr[ y ][ x * CHANNEL_COUNT + B_CHANNEL ] = v; }
    inline void setG( int x, int y, T v ) { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK m_ptr[ y ][ x * CHANNEL_COUNT + G_CHANNEL ] = v; }
    inline void setR( int x, int y, T v ) { assertStatic( CHANNEL_COUNT >= 3 ); IMGCHK m_ptr[ y ][ x * CHANNEL_COUNT + R_CHANNEL ] = v; }

    // get image info
    inline int width() const { return m_width; }
    inline int height() const { return m_height; }
    inline int rowBytes() const { return m_rowBytes; }
    inline int memUsed() const { return sizeof( Image<T,CHANNEL_COUNT> ) + m_rowBytes * m_height; }
    inline int depth() const { return sizeof(T) * 8; }
    inline int channelCount() const { return CHANNEL_COUNT; }
    inline bool isFloat() const { return typeid( T ) == typeid( float ); }

    /// check whether position is in image bounds
    inline bool inBounds( int x, int y ) const { return x >= 0 && x < m_width && y >= 0 && y < m_height; }
    inline bool inBounds( float x, float y ) const { return x >= 0 && x < m_width - 1 && y >= 0 && y < m_height - 1; }

    /// access image data via pointers
    inline T *raw() { return m_raw; }
    inline const T *rawConst() const { return m_raw; }
    inline T **pointers() { return m_ptr; }

    /// clear the image to the specified color
    void clear( T r, T g, T b );
    void clear( T v );

private:

    // common constructor code
    void alloc( int width, int height );

    // image data
    T *m_raw;
    T **m_ptr;
    int m_width;
    int m_height;
    int m_rowBytes;

    // indicates whether to delete m_raw
    bool m_deleteRaw;

    // disable assignment operator
    Image &operator=( const Image &img );

// the rest of the class is only used for OpenCV interaction 
#ifdef USE_OPENCV
public:

    /// create an Image object wrapping an IplImage object
    Image( IplImage *iplImage, bool deallocIplImage );

    /// return IplImage object; create object if needed
    inline IplImage *iplImage() { if (m_iplImage == NULL) createIplImage(); return m_iplImage; } 
    inline const IplImage *iplImage() const { if (m_iplImage == NULL) createIplImage(); return m_iplImage; }

private:

    /// create object if needed (m_iplImage is mutable so this func can be const)
    void createIplImage() const;

    /// the stored IplImage, if any
    mutable IplImage *m_iplImage;
    mutable bool m_deallocIplImage;

#endif // USE_OPENCV
};


//-------------------------------------------
// IMAGE CLASS IMPLEMENTATION 
//-------------------------------------------


// deallocate image data
template<typename T, int CHANNEL_COUNT> Image<T, CHANNEL_COUNT>::~Image() {
    if (m_deleteRaw)
        delete [] m_raw;
    delete [] m_ptr;
#ifdef USE_OPENCV
    if (m_deallocIplImage) {
        if (m_deleteRaw)
            delete m_iplImage;
        else
            cvReleaseImage( &m_iplImage );
    }
#endif
}


// common constructor code
template<typename T, int CHANNEL_COUNT> void Image<T, CHANNEL_COUNT>::alloc( int width, int height ) {
    m_width = width;
    m_height = height;

    // require quad-word alignment for rows 
    // (round up to nearest mult of 8)
    m_rowBytes = (m_width * sizeof( T ) * CHANNEL_COUNT + 7) & (0xffff - 7);
    int rowWidth = m_rowBytes / sizeof( T );
    assertDebug( rowWidth >= m_width );
    int size = rowWidth * m_height; // total number of elements

    // alloc pixel data
    m_raw = new T[ size ];
    if (m_raw == NULL) fatalError( "error allocating ImageColor data" );
    m_deleteRaw = true;

    // check alignment
    if (((long) m_raw) & 7)
        fatalError("not quad-word aligned");

    // create row pointers
    m_ptr = new T*[ m_height ];
    if (m_ptr == NULL) fatalError( "error allocating ImageColor pointers" );
    for (int i = 0; i < m_height; i++)
        m_ptr[ i ] = m_raw + i * rowWidth;

    // init IPL pointer, if defined
#ifdef USE_OPENCV
    m_iplImage = NULL;
    m_deallocIplImage = false;
#endif
}


/// clear the image to the specified color
template<typename T, int CHANNEL_COUNT> void Image<T, CHANNEL_COUNT>::clear( T r, T g, T b ) {
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            setRGB( x, y, r, g, b );
        }
    }
}


/// clear the image to the specified color
template<typename T, int CHANNEL_COUNT> void Image<T, CHANNEL_COUNT>::clear( T v ) {
    for (int y = 0; y < m_height; y++) {
        for (int x = 0; x < m_width; x++) {
            m_ptr[ y ][ x ] = v;
        }
    }
}


/// performs bilinear interpolation at specified point
template<typename T, int CHANNEL_COUNT> float Image<T, CHANNEL_COUNT>::interp( float x, float y ) const {
    IMGCHK_INTERP
    int xInt = (int) x; 
    int yInt = (int) y; 
    float val = 0;
    assertDebug( xInt >= 0 && yInt >= 0 && xInt < m_width && yInt < m_height );
    if (xInt == m_width - 1) { // fix(faster): make faster?
        if (yInt == m_height - 1) {
            val = (float) m_ptr[ yInt ][ xInt ];
        } else {
            float yFrac = y - yInt;
            val = (1.0f - yFrac) * m_ptr[ yInt     ][ xInt ] + 
                          yFrac  * m_ptr[ yInt + 1 ][ xInt ];
        }
    } else if (yInt == m_height - 1) {
        float xFrac = x - xInt;
        val = (1.0f - xFrac) * m_ptr[ yInt ][ xInt     ] + 
                      xFrac  * m_ptr[ yInt ][ xInt + 1 ];
    } else {
        float xFrac = x - xInt;
        float yFrac = y - yInt;
        val = (1.0f - xFrac) * (1.0f - yFrac) * m_ptr[ yInt     ][ xInt     ] + 
              (1.0f - xFrac) *         yFrac  * m_ptr[ yInt + 1 ][ xInt     ] + 
                      xFrac  * (1.0f - yFrac) * m_ptr[ yInt     ][ xInt + 1 ] + 
                      xFrac  *         yFrac  * m_ptr[ yInt + 1 ][ xInt + 1 ];
    }
    return val;
}


/// performs bilinear interpolation at specified point
template<typename T, int CHANNEL_COUNT> float Image<T, CHANNEL_COUNT>::interp( float x, float y, int c ) const {
    IMGCHK_INTERP
    int xInt = (int) x; 
    int yInt = (int) y; 
    float xFrac = x - xInt;
    float yFrac = y - yInt;
    float val = (1.0f - xFrac) * (1.0f - yFrac) * m_ptr[ yInt     ][ (xInt    ) * CHANNEL_COUNT + c ] + 
                  (1.0f - xFrac) *         yFrac  * m_ptr[ yInt + 1 ][ (xInt    ) * CHANNEL_COUNT + c ] + 
                        xFrac  * (1.0f - yFrac) * m_ptr[ yInt     ][ (xInt + 1) * CHANNEL_COUNT + c ] + 
                        xFrac  *         yFrac  * m_ptr[ yInt + 1 ][ (xInt + 1) * CHANNEL_COUNT + c ];
    return val;
}


//-------------------------------------------
// IMAGE CLASS IMPLEMENTATION - WITH OPENCV
//-------------------------------------------
#ifdef USE_OPENCV


/// create an ImageColor object wrapping an IplImage object
template<typename T, int CHANNEL_COUNT> Image<T, CHANNEL_COUNT>::Image( IplImage *iplImage, bool deallocIplImage ) {
    assertAlways( iplImage->nChannels == CHANNEL_COUNT ); 
//    assertAlways( iplImage->origin == IPL_ORIGIN_BL );
//    assertAlways( sizeof( T ) == 1 ? iplImage->depth == 8 : iplImage->depth == 4 ); // assuming U or F type
    m_width = iplImage->width;
    m_height = iplImage->height;
    int rowWidth = iplImage->widthStep;
    m_rowBytes = rowWidth * sizeof( T );
    m_raw = (T *) iplImage->imageData;
    m_iplImage = iplImage;
    m_deleteRaw = false;
    m_deallocIplImage = deallocIplImage;

    // create row pointers
    m_ptr = new T*[ m_height ];
    if (m_ptr == NULL) fatalError( "error allocating ImageColor pointers" );
    for (int i = 0; i < m_height; i++)
        m_ptr[ i ] = m_raw + i * rowWidth;
}


/// create object if needed (m_iplImage is mutable so this func can be const)
template<typename T, int CHANNEL_COUNT> void Image<T, CHANNEL_COUNT>::createIplImage() const {
    assert( m_iplImage == NULL );
    m_iplImage = new IplImage;
    assert( m_iplImage );
    m_deallocIplImage = true;
    // fix(clean): make sure every field is initialized; use cvInitImageHeader?
    m_iplImage->nSize = sizeof( IplImage );
    m_iplImage->ID = 0;
    m_iplImage->nChannels = CHANNEL_COUNT;
    m_iplImage->width = m_width;
    m_iplImage->height = m_height;
    m_iplImage->dataOrder = IPL_DATA_ORDER_PIXEL;
    m_iplImage->origin = IPL_ORIGIN_BL;
    if (isFloat()) {
        m_iplImage->depth = IPL_DEPTH_32F;
    } else if (depth() == 8) {
        m_iplImage->depth = IPL_DEPTH_8U;
    } else if (depth() == 16) {
        m_iplImage->depth = IPL_DEPTH_16U;
    } else {
        m_iplImage->depth = IPL_DEPTH_32S;
    }
    m_iplImage->align = IPL_ALIGN_QWORD;
    m_iplImage->imageSize = m_rowBytes * m_height;
    m_iplImage->imageData = (char *) m_raw;
    m_iplImage->widthStep = m_rowBytes;
    m_iplImage->imageDataOrigin = (char *) m_raw;
    m_iplImage->imageData = (char *) m_raw;
    m_iplImage->roi = NULL;
    m_iplImage->maskROI = NULL;
    m_iplImage->imageId = NULL;
    m_iplImage->tileInfo = NULL;
}


#endif // USE_OPENCV


//-------------------------------------------
// IMAGE TYPEDEFS 
//-------------------------------------------


/// common image types
typedef Image<unsigned char, 3> ImageColorU;
typedef Image<unsigned short, 3> ImageColorS;
typedef Image<int, 3> ImageColorI;
typedef Image<float, 3> ImageColorF;
typedef Image<unsigned char, 1> ImageGrayU;
typedef Image<unsigned short, 1> ImageGrayS;
typedef Image<int, 1> ImageGrayI;
typedef Image<float, 1> ImageGrayF;


} // end namespace sbl
#endif // _SBL_IMAGE_H_

