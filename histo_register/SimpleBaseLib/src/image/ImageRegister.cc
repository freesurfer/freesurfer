// Licensed under MIT license; see license.txt.

#include <sbl/image/ImageRegister.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/math/Optimizer.h>
#include <sbl/math/MathUtil.h>
#include <sbl/image/ImageUtil.h>
namespace sbl {


/// computes the mean-abs image difference given an image transformation;
/// a step value greater than one allows faster (less accurate) optimization by ignoring some pixels;
/// the xBorder and yBorder specify areas to be ignored for the objective function;
/// the offsetBound value specifies the largest allowed translation between the images
double evalImageTransform( const ImageTransform &transform, const ImageGrayU &src, const ImageGrayU &dest, int step, int xBorder, int yBorder, const ImageGrayU *srcMask, const ImageGrayU *destMask, bool interp, bool verbose ) {
    int xDestMin = xBorder; // note that we're using border on both source and dest
    int xDestMax = dest.width() - 1 - xBorder;
    int yDestMin = yBorder;
    int yDestMax = dest.height() - 1 - yBorder;
    float sumDiff = 0;
    int count = 0;
    int width = src.width(), height = src.height();
    for (int y = yBorder; y < height - yBorder; y += step) {
        for (int x = xBorder; x < width - xBorder; x += step) {
            if (srcMask == NULL || srcMask->data( x, y )) {
                float xSrc = (float) x;
                float ySrc = (float) y;
                float xDest = transform.xTransform( xSrc, ySrc );
                float yDest = transform.yTransform( xSrc, ySrc );
                int xDestInt = round( xDest );
                int yDestInt = round( yDest );
                if (xDestInt >= xDestMin && xDestInt <= xDestMax && yDestInt >= yDestMin && yDestInt <= yDestMax && (destMask == NULL || destMask->data( xDestInt, yDestInt ) )) {
                    int vDest = interp ? (int) dest.interp( xDest, yDest ) : dest.data( xDestInt, yDestInt );
                    int diff = src.data( x, y ) - vDest;
                    if (diff < 0)
                        diff = -diff;
                    sumDiff += (float) diff;
                    count++;
                }
            }
        }
    }
    if (verbose)
        disp( 1, "count: %d, sumDiff: %f", count, sumDiff );
    return count ? sumDiff / count : 1000;
}


/// The RegistrationObjective class represents an image registration objective function.
class RegistrationObjective : public Objective {
public:

    // basic constructor
    RegistrationObjective( const ImageGrayU &src, const ImageGrayU &dest, int step, int xBorder, int yBorder, const ImageGrayU *srcMask, const ImageGrayU *destMask, bool interp ) 
                : m_src( src ), m_dest( dest ) {
        m_step = step;
        m_xBorder = xBorder;
        m_yBorder = yBorder;
        m_srcMask = srcMask;
        m_destMask = destMask;
        m_interp = interp;
        m_verbose = false;
    }

    /// evaluate objective function at given point
    double eval( const VectorD &point ) {

        // convert parameter point to transform
        VectorF transformParams = toFloat( point );
        ImageTransform transform( transformParams );

        // evalutate the match given this transform
        return evalImageTransform( transform, m_src, m_dest, m_step, m_xBorder, m_yBorder, m_srcMask, m_destMask, m_interp, m_verbose );
    }

    // enable diagnostic display
    inline void setVerbose() { m_verbose = true; }

private:

    // internal data
    const ImageGrayU &m_src;
    const ImageGrayU &m_dest;
    int m_step;
    int m_xBorder;
    int m_yBorder;
    const ImageGrayU *m_srcMask;
    const ImageGrayU *m_destMask;
    bool m_interp;
    bool m_verbose;
};


/// registers a pair of images using to minimize mean-abs difference;
/// a step value greater than one allows faster (less accurate) optimization by ignoring some pixels;
/// the xBorder and yBorder specify areas to be ignored for the objective function;
/// the offsetBound value specifies the largest allowed translation between the images
aptr<ImageTransform> registerUsingImageTransform( const ImageGrayU &src, const ImageGrayU &dest, int transformParamCount, int step, int xBorder, int yBorder, float offsetBound, const ImageTransform *initTransform, const ImageGrayU *srcMask, const ImageGrayU *destMask, bool interp ) {
    assertAlways( transformParamCount == 2 || transformParamCount == 6 );

    // prepare objective function
    RegistrationObjective obj( src, dest, step, xBorder, yBorder, srcMask, destMask, interp );

    // determine starting point
    VectorD start( transformParamCount );
    start.clear( 0 );
    if (transformParamCount == 6) {
        start[ 2 ] = 1;
        start[ 5 ] = 1;
    }
    if (initTransform) {
        assertAlways( initTransform->paramCount() == transformParamCount );
        for (int i = 0; i < transformParamCount; i++)
            start[ i ] = initTransform->param( i );
    }

    // determine bounds
    VectorD lBound( transformParamCount ), uBound( transformParamCount );
    lBound[ 0 ] = start[ 0 ] - offsetBound;
    lBound[ 1 ] = start[ 1 ] - offsetBound;
    uBound[ 0 ] = start[ 0 ] + offsetBound;
    uBound[ 1 ] = start[ 1 ] + offsetBound;
    if (transformParamCount == 6) {
        lBound[ 2 ] = start[ 2 ] - 0.2;
        lBound[ 5 ] = start[ 2 ] - 0.2;
        lBound[ 3 ] = start[ 2 ] - 0.1;
        lBound[ 4 ] = start[ 2 ] - 0.1;
        uBound[ 2 ] = start[ 2 ] + 0.2;
        uBound[ 5 ] = start[ 2 ] + 0.2;
        uBound[ 3 ] = start[ 2 ] + 0.1;
        uBound[ 4 ] = start[ 2 ] + 0.1;
    }

    // create and configure the optimizer
    SimplexOptimizer opt( obj );
    opt.setStart( start );
    opt.setBounds( lBound, uBound );

    // run optimizer
    VectorD result = opt.run();
//    obj.setVerbose();
//    obj.eval( result );

    // compute resulting transform 
    aptr<ImageTransform> transform( new ImageTransform( toFloat( result ) ) );
    return transform;
}


} // end namespace sbl

