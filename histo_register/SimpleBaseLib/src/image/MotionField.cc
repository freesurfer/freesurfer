// Licensed under MIT license; see license.txt.

#include <sbl/image/MotionField.h>
#include <sbl/image/ImageUtil.h>
#include <sbl/image/ImageTransform.h>
namespace sbl {


//-------------------------------------------
// MOTION FIELD CLASS
//-------------------------------------------


/// create an uninitialized motion field
MotionField::MotionField( int width, int height ) { 
    m_width = width;
    m_height = height;
    m_u = new ImageGrayF( width, height );
    m_v = new ImageGrayF( width, height );
    m_occMask = NULL;

    // init meta-data
    m_buildTime = 0;
    m_srcFrameIndex = -1;
    m_destFrameIndex = -1;
    m_originalWidth = -1;
    m_originalHeight = -1;
}


// basic copy constructor
MotionField::MotionField( const MotionField &mf ) { 
    m_width = mf.width();
    m_height = mf.height();
    m_u = new ImageGrayF( mf.uRef() );
    m_v = new ImageGrayF( mf.vRef() );
    m_occMask = mf.occDefined() ? new ImageGrayU( mf.occRef() ) : NULL;

    // copy meta-data
    m_buildTime = mf.buildTime();
    m_algorithmName = mf.algorithmName();
    m_videoName = mf.videoName();
    m_srcFrameIndex = mf.srcFrameIndex();
    m_destFrameIndex = mf.destFrameIndex();
    m_originalWidth = mf.originalWidth();
    m_originalHeight = mf.originalHeight();
}


// deallocate motion field
MotionField::~MotionField() {
    delete m_u;
    delete m_v;
    if (m_occMask) delete m_occMask;
}


/// enable occlusion mask
void MotionField::enableOcc() {
    if (m_occMask == NULL) {
        m_occMask = new ImageGrayU( width(), height() );
        m_occMask->clear( 0 );
    }
}


/// compute number of occluded pixels (according to occlusion mask)
int MotionField::occCount() const {
    int count = 0;
    if (occDefined()) {
        int width = m_occMask->width(), height = m_occMask->height();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                if (m_occMask->data( x, y ) > 128)
                    count++;
            }
        }
    }
    return count;
}


/// increases the size of the motion field
void MotionField::resize( int newWidth, int newHeight, bool rescale ) {

    // if no size change, do nothing
    if (newWidth == m_width && newHeight == m_height) 
        return;

    // allocate new flow field
    ImageGrayF *u = sbl::resize( *m_u, newWidth, newHeight, true ).release();
    ImageGrayF *v = sbl::resize( *m_v, newWidth, newHeight, true ).release();
    ImageGrayU *occMask = NULL;
    if (m_occMask)
        occMask = sbl::resize( *m_occMask, newWidth, newHeight, true ).release();

    // rescale motion vectors 
    if (rescale) {
        multiply( *u, (float) newWidth / (float) m_width, *u );
        multiply( *v, (float) newHeight / (float) m_height, *v );
    }

    // replace existing data with new data
    delete m_u;
    delete m_v;
    if (m_occMask) delete m_occMask;
    m_u = u;
    m_v = v;
    m_occMask = occMask;
    m_width = newWidth;
    m_height = newHeight;
}


/// multiplies each vector by given factor
void MotionField::rescale(float scale) {
    multiply( *m_u, scale, *m_u );
    multiply( *m_v, scale, *m_v );
}


/// reset motion vectors to zero
void MotionField::clear() {
    m_u->clear( 0 );
    m_v->clear( 0 );
}


/// map the given destination image back to the source coordinates
// fix(faster): use opencv?
aptr<ImageGrayU> MotionField::mapBackward( const ImageGrayU &img, int fillColor, float frac ) const {
    int width = m_u->width(), height = m_u->height();
    aptr<ImageGrayU> mapped( new ImageGrayU( width, height ) );
    mapped->clear( fillColor );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float xDest = x + m_u->data( x, y ) * frac;
            float yDest = y + m_v->data( x, y ) * frac;
            if (img.inBounds( xDest, yDest )) 
                mapped->data( x, y ) = (int) img.interp( xDest, yDest );
        }
    }
    return mapped;
}


/// map the given destination image back to the source coordinates
// fix(faster): use opencv?
aptr<ImageColorU> MotionField::mapBackward( const ImageColorU &img, int fillColor, float frac ) const {
    int width = m_u->width(), height = m_u->height();
    aptr<ImageColorU> mapped( new ImageColorU( width, height ) );
    mapped->clear( fillColor );
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            float xDest = x + m_u->data( x, y ) * frac;
            float yDest = y + m_v->data( x, y ) * frac;
            if (img.inBounds( xDest, yDest )) {
                for (int c = 0; c < 3; c++)
                    mapped->data( x, y, c ) = (int) img.interp( xDest, yDest, c );
            }
        }
    }
    return mapped;
}


} // end namespace sbl

