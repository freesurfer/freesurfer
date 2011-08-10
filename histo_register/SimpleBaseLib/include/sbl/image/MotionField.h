#ifndef _SBL_MOTION_FIELD_H_
#define _SBL_MOTION_FIELD_H_
#include <sbl/core/Pointer.h>
#include <sbl/image/Image.h>
namespace sbl {


/// The MotionField class represents a 2D correspondence field between a pair of images (e.g., optical flow).
class MotionField {
public:

    /// create an uninitialized motion field
    MotionField( int width, int height );

    // basic copy constructor
    MotionField( const MotionField &mf );

    // deallocate motion field
    ~MotionField();

    // access motion data
    inline int width() const { return m_width; }
    inline int height() const { return m_height; }
    inline float u( int x, int y ) const { return m_u->data( x, y ); }
    inline float v( int x, int y ) const { return m_v->data( x, y ); }
    inline float uInterp( float x, float y ) const { return m_u->interp( x, y ); }
    inline float vInterp( float x, float y ) const { return m_v->interp( x, y ); }
    inline void setU( int x, int y, float val ) { m_u->data( x, y ) = val; }
    inline void setV( int x, int y, float val ) { m_v->data( x, y ) = val; }
    inline ImageGrayF &uRef() { return *m_u; }
    inline ImageGrayF &vRef() { return *m_v; }
    inline const ImageGrayF &uRef() const { return *m_u; }
    inline const ImageGrayF &vRef() const { return *m_v; }
    inline bool inBounds( int x, int y ) const { return m_u->inBounds( x, y ); }
    inline bool inBounds( float x, float y ) const { return m_u->inBounds( x, y ); }

    /// enable occlusion mask
    void enableOcc();

    // access occlusion mask
    inline bool occDefined() const { return m_occMask != NULL; }
    inline int occ( int x, int y ) const { return m_occMask->data( x, y ); }
    inline int occSafe( int x, int y ) const { return m_occMask ? m_occMask->data( x, y ) :  0; }
    inline void setOcc( int x, int y, int val ) { m_occMask->data( x, y ) = val; }
    inline ImageGrayU &occRef() { return *m_occMask; }
    inline const ImageGrayU &occRef() const { return *m_occMask; }

    /// compute number of occluded pixels (according to occlusion mask)
    int occCount() const;

    /// increases/decreases the size of the motion field
    void resize( int newWidth, int newHeight, bool rescale );

    /// multiplies each vector by given factor
    void rescale( float scale );

    /// reset motion vectors to zero
    void clear();

    // get/set meta-data
    inline void setBuildTime( float buildTime ) { m_buildTime = buildTime; }
    inline float buildTime() const { return m_buildTime; }
    inline const String &algorithmName() const { return m_algorithmName; }
    inline const String &videoName() const { return m_videoName; }
    inline int srcFrameIndex() const { return m_srcFrameIndex; }
    inline int destFrameIndex() const { return m_destFrameIndex; }
    inline int originalWidth() const { return m_originalWidth; }
    inline int originalHeight() const { return m_originalHeight; }

    /// returns number of bytes used by this object
    int memUsed() const;

    /// map the given destination image back to the source coordinates
    aptr<ImageGrayU> mapBackward( const ImageGrayU &img, int fillColor, float frac = 1.0f ) const;
    aptr<ImageColorU> mapBackward( const ImageColorU &img, int fillColor, float frac = 1.0f ) const;

private:

    // motion data
    int m_width;
    int m_height;
    ImageGrayF *m_u;
    ImageGrayF *m_v;
    ImageGrayU *m_occMask;

    // meta-data
    float m_buildTime;
    String m_algorithmName;
    String m_videoName;
    int m_srcFrameIndex;
    int m_destFrameIndex;
    int m_originalWidth;
    int m_originalHeight;

    // disable assignment operator
    MotionField &operator=( const MotionField &x );
};


} // end namespace sbl
#endif // _SBL_MOTION_FIELD_H_

