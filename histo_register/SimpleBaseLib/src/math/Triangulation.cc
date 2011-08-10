// Licensed under MIT license; see license.txt.

#ifdef USE_CDT
#include <sbl/math/Triangulation.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/math/MathUtil.h>
#include <../external/CDT/CDT.h>
namespace sbl {


/// returns matrix of edges; each row is (x1, y1, x2, y2)
aptr<MatrixF> triangulate( const MatrixF &points ) {

    // get bounds
    int pointCount = points.rows();
    float xMin = points.data( 0, 0 ), xMax = points.data( 0, 0 );
    float yMin = points.data( 0, 1 ), yMax = points.data( 0, 1 );
    for (int i = 1; i < pointCount; i++) {
        float x = points.data( i, 0 );
        float y = points.data( i, 1 );
        if (x < xMin) xMin = x;
        if (x > xMax) xMax = x;
        if (y < yMin) yMin = y;
        if (y > yMax) yMax = y;
    }

    // create mesh using bounds
    Mesh *mesh = new Mesh( Point2d( xMin, yMin ),
                           Point2d( xMax, yMin ),
                           Point2d( xMax, yMax ),
                           Point2d( xMin, yMax ) );

    // add points
    for (int i = 0; i < pointCount; i++) 
        mesh->InsertSite( Point2d( points.data( i, 0 ), points.data( i, 1 ) ));

    // count edges so we can allocate a matrix of the correct size
    int edgeCount = 0;
    for (LlistPos p = mesh->edges.first(); !mesh->edges.isEnd(p); p = mesh->edges.next(p)) 
        edgeCount++;

    // add edges to matrix
    aptr<MatrixF> edges( new MatrixF( edgeCount, 4 ) );
    int edgeIndex = 0;
    for (LlistPos p = mesh->edges.first(); !mesh->edges.isEnd( p ); p = mesh->edges.next( p )) {
        QuadEdge *qp = (QuadEdge *) mesh->edges.retrieve( p );
        float x1 = (float) qp->edges()[ 0 ].Org2d()[ 0 ];
        float y1 = (float) qp->edges()[ 0 ].Org2d()[ 1 ];
        float x2 = (float) qp->edges()[ 2 ].Org2d()[ 0 ];
        float y2 = (float) qp->edges()[ 2 ].Org2d()[ 1 ];
        edges->data( edgeIndex, 0 ) = x1;
        edges->data( edgeIndex, 1 ) = y1;
        edges->data( edgeIndex, 2 ) = x2;
        edges->data( edgeIndex, 3 ) = y2;
        edgeIndex++;
    }

    // clean up
    delete mesh;
    return edges;
}


/// returns matrix of edges: each row is pair (i1, i2) of point indices
aptr<MatrixI> triangulateIndex( const MatrixF &points ) {

    // get bounds
    int pointCount = points.rows();
    float xMin = points.data( 0, 0 ), xMax = points.data( 0, 0 );
    float yMin = points.data( 0, 1 ), yMax = points.data( 0, 1 );
    for (int i = 1; i < pointCount; i++) {
        float x = points.data( i, 0 );
        float y = points.data( i, 1 );
        if (x < xMin) xMin = x;
        if (x > xMax) xMax = x;
        if (y < yMin) yMin = y;
        if (y > yMax) yMax = y;
    }

    // extend bounds a bit
    xMin -= 0.1f * (xMax - xMin);
    xMax += 0.1f * (xMax - xMin);
    yMin -= 0.1f * (yMax - yMin);
    yMax += 0.1f * (yMax - yMin);

    // store indices in mesh
    VectorI seq = sequenceI( 0, pointCount - 1 );

    // create mesh using bounds (bound point have no tags)
    Mesh *mesh = new Mesh( Point2d( xMin, yMin ),
                           Point2d( xMax, yMin ),
                           Point2d( xMax, yMax ),
                           Point2d( xMin, yMax ) );

    // add points
    for (int i = 0; i < pointCount; i++) 
        mesh->InsertSite( Point2d( points.data( i, 0 ), points.data( i, 1 ), &(seq.data( i )) ));

    // count edges so we can allocate a matrix of the correct size
    int edgeCount = 0;
    for (LlistPos p = mesh->edges.first(); !mesh->edges.isEnd(p); p = mesh->edges.next(p)) {
        QuadEdge *qp = (QuadEdge *)mesh->edges.retrieve(p);
        int *t1 = (int *) qp->edges()[0].Org2d().tag;
        int *t2 = (int *) qp->edges()[2].Org2d().tag;
        if (t1 && t2) {
            edgeCount++;
        }
    }

    // add edges to matrix
    aptr<MatrixI> edges( new MatrixI( edgeCount, 2 ) );
    int edgeIndex = 0;
    for (LlistPos p = mesh->edges.first(); !mesh->edges.isEnd( p ); p = mesh->edges.next( p )) {
        QuadEdge *qp = (QuadEdge *) mesh->edges.retrieve( p );
        int *t1 = (int *) qp->edges()[ 0 ].Org2d().tag;
        int *t2 = (int *) qp->edges()[ 2 ].Org2d().tag;

        // if tags defined, add to edge set
        if (t1 && t2) {
            edges->data( edgeIndex, 0 ) = *t1;
            edges->data( edgeIndex, 1 ) = *t2;
            edgeIndex++;
        }
    }

    // clean up
    delete mesh;
    return edges;
}


// stores results computed by PerTriangle
int g_triCount = 0;
ImageGrayF *g_triResult = NULL;


/// callback used for TriangulationInterpolation
void perTriangle( Point2d p1, Point2d p2, Point2d p3 ) {

    // get values at points
    float *pval1 = (float *) p1.tag;
    float *pval2 = (float *) p2.tag;
    float *pval3 = (float *) p3.tag;
    float val1 = pval1 ? *pval1 : 0;
    float val2 = pval2 ? *pval2 : 0;
    float val3 = pval3 ? *pval3 : 0;

    // get locations of points
    float x1 = (float) p1[ 0 ];
    float y1 = (float) p1[ 1 ];
    float x2 = (float) p2[ 0 ];
    float y2 = (float) p2[ 1 ];
    float x3 = (float) p3[ 0 ];
    float y3 = (float) p3[ 1 ];

    // prepare for interpolation
    float ux = x2 - x1;
    float uy = y2 - y1;
    float vx = x3 - x1;
    float vy = y3 - y1;
    float uVal = val2 - val1;
    float vVal = val3 - val1;
    float uStep = 0.5f / ((float) fabs( ux ) + (float) fabs( uy ) + 1);
    float vStep = 0.5f / ((float) fabs( vx ) + (float) fabs( vy ) + 1);

    // do interpolation
    for (float uFrac = 0; uFrac <= 1; uFrac += uStep) {
        for (float vFrac = 0; vFrac <= 1 - uFrac; vFrac += vStep) {
            float x = x1 + ux * uFrac + vx * vFrac;
            float y = y1 + uy * uFrac + vy * vFrac;
            int xInt = round( x ), yInt = round( y );
            float val = val1 + uVal * uFrac + vVal * vFrac;
            g_triResult->data( xInt, yInt ) = val;
        }
    }

    // keep diagnostic info
    g_triCount++;
}


/// get value of point nearest to query point
float *nearestValue( float x, float y, const MatrixF &points, VectorF &values ) {
    int count = points.rows(); 
    float *value = NULL;
    float bestDistSqd = 0;
    for (int i = 0; i < count; i++) {
        float xDiff = x - (float) points.data( i, 0 );
        float yDiff = y - (float) points.data( i, 1 );
        float distSqd = xDiff * xDiff + yDiff * yDiff;
        if (distSqd < bestDistSqd || i == 0) {
            bestDistSqd = distSqd;
            value = &(values.data( i ));
        }
    }
    return value;
}


/// use triangulation to interpolate values across an image
aptr<ImageGrayF> triangulationInterpolation( const MatrixF &points, VectorF &values, int width, int height ) {
    ImageGrayF *result = new ImageGrayF( width, height );
    result->clear( -777 );

    // create mesh using image bounds
    float xMin = 0, xMax = (float) (width - 1), yMin = 0, yMax = (float) (height - 1);
    Mesh *mesh = new Mesh( Point2d( xMin, yMin, nearestValue( xMin, yMin, points, values ) ),
                           Point2d( xMax, yMin, nearestValue( xMax, yMin, points, values ) ),
                           Point2d( xMax, yMax, nearestValue( xMax, yMax, points, values ) ),
                           Point2d( xMin, yMax, nearestValue( xMin, yMax, points, values ) ) );

    // add points
    int pointCount = points.rows();
    for (int i = 0; i < pointCount; i++) 
        mesh->InsertSite( Point2d( points.data( i, 0 ), points.data( i, 1 ), &(values.data( i )) ));

    // store data for per-triangle callback
    g_triCount = 0;
    g_triResult = result;

    // execute per-triangle callback for each triangle
    mesh->ApplyTriangles( perTriangle );

    // fill holes 
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            if (result->data( x, y ) == -777) {
                if (x > 0 && result->data( x - 1, y ) != -777) {
                    result->data( x, y ) = result->data( x - 1, y );
                } else if (y > 0 && result->data( x, y - 1 ) != -777) {
                    result->data( x, y ) = result->data( x, y - 1 );
                } else if (x < width - 1 && result->data( x + 1, y ) != -777) {
                    result->data( x, y ) = result->data( x + 1, y );
                } else if (y < height - 1 && result->data( x, y + 1 ) != -777) {
                    result->data( x, y ) = result->data( x, y + 1 );
                } else {
                    result->data( x, y ) = 0;
                }
            }
        }
    }

    // done; clean up
    delete mesh;
    return aptr<ImageGrayF>( result );
}


} // end namespace sbl
#endif // USE_CDT

