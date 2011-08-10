// Licensed under MIT license; see license.txt.

#include <sbl/math/KMeans.h>
#include <sbl/core/Command.h>
#include <sbl/math/MathUtil.h>
#include <sbl/math/VectorUtil.h>
#include <sbl/math/MatrixUtil.h>
#include <sbl/other/Plot.h> // for test command
#include <sbl/image/ImageDraw.h> // for test command
namespace sbl {


/// sum of squared differences
float distSqd( float *p1, float *p2, int len ) {
    float distSqd = 0;
    for (int i = 0; i < len; i++) {
        float diff = p1[ i ] - p2[ i ];
        distSqd += diff * diff;
    }
    return distSqd;
}


/// pick initial cluster centers, such that each is furthest from the others
void chooseInitialMeans( const MatrixF &points, MatrixF &means ) {
    int pointCount = points.rows();
    int clusterCount = means.rows();
    int dimCount = means.cols();
    
    // pick first point at random
    int index = randomInt( 0, pointCount - 1 );
    for (int d = 0; d < dimCount; d++) 
        means.data( 0, d ) = points.data( index, d );
    
    // pick each other mean
    for (int j = 1; j < clusterCount; j++) {

        // find furthest point
        float maxDist = 0;
        int maxDistIndex = 0;
        for (int i = 0; i < pointCount; i++) {

            // find nearest mean
            float minDist = 1e20f;
            int minDistIndex = 0;
            for (int k = 0; k < j; k++) {
                float dist = distSqd( points.dataRow( i ), means.dataRow( k ), dimCount );
                if (dist < minDist) {
                    minDist = dist;
                    minDistIndex = k;
                }
            }
            if (minDist > maxDist) {
                maxDist = minDist;
                maxDistIndex = i;
            }
        }

        // copy this point
        for (int d = 0; d < dimCount; d++) 
            means.data( j, d ) = points.data( maxDistIndex, d );
    }
}


/// run the K-means cluster algorithm; returns the cluster means and assignment of points to clusters
void kMeans( const MatrixF &points, int clusterCount, aptr<MatrixF> &means, aptr<VectorI> &assign ) {
    int maxIter = 100;
    int pointCount = points.rows();
    int dimCount = points.cols();
    float **pointsData = points.dataPtr();

    // allocate cluster assignments
    assign.reset( new VectorI( pointCount ) );
    int *assignData = assign->dataPtr();

    // allocate means
    means.reset( new MatrixF( clusterCount, dimCount ) );
    float **meansData = means->dataPtr();

    // choose well-distributed initial means
    chooseInitialMeans( points, *means );

    // loop until done
    for (int iter = 0; iter < maxIter; iter++) {
        
        // assign each point to nearest mean
        for (int i = 0; i < pointCount; i++) {

            // find closest mean
            float bestDist = 1e20f;
            int bestCluster = 0;
            for (int j = 0; j < clusterCount; j++) {
                float dist = distSqd( meansData[ j ], pointsData[ i ], dimCount );
                if (dist < bestDist) {
                    bestDist = dist;
                    bestCluster = j;
                }
            }
            assignData[ i ] = bestCluster;
        }

        // update means
        for (int j = 0; j < clusterCount; j++) {
            float *mean = meansData[ j ];
            for (int d = 0; d < dimCount; d++) {
                mean[ d ] = 0;
            }
            int clusterMemberCount = 0;
            for (int i = 0; i < pointCount; i++) {
                if (assignData[ i ] == j) {
                    clusterMemberCount++;
                    for (int d = 0; d < dimCount; d++) 
                        mean[ d ] += pointsData[ i ][ d ];
                }
            }
            for (int d = 0; d < dimCount; d++) {
                mean[ d ] /= clusterMemberCount;
            }
        }
    }
}


/// transforms data into space suitable for clustering as described in "On Spectral Cluster: Analysis and an Algorithm" by Ng, Jordan, and Weiss;
/// assumes affinity is symmetric and 0 on diagonal
aptr<MatrixF> spectralTransform( MatrixF &affinity, int outputDimCount, bool verbose ) {
    int pointCount = affinity.rows();

    // normalize affinity matrix; multiply each entry by sqrt( column sum ) * sqrt( row sum )
    VectorF cSum = colSum( affinity );
    for (int i = 0; i < pointCount; i++)
        cSum[ i ] = sqrtf( cSum[ i ] );
    for (int i = 0; i < pointCount; i++)
        for (int j = 0; j < pointCount; j++) 
            affinity.data( i, j ) *= cSum[ i ] * cSum[ j ];

    // compute eigenvectors 
    VectorF eigenVals( pointCount );
    aptr<MatrixF> eigenVects = eigenSymmetric( affinity, eigenVals );
    
    // sort eigenvalue magnitude
    for (int i = 0; i < pointCount; i++)
        if (eigenVals[ i ] < 0)
            eigenVals[ i ] = -eigenVals[ i ];
    VectorI sortInd = reverseSortIndex( eigenVals );
    if (verbose)
        disp( 2, "abs eigs min/max/mean: %f, %f, %f", eigenVals.min(), eigenVals.mean(), eigenVals.max() );

    // create output space by stacking eigenvectors as columns
    aptr<MatrixF> output( new MatrixF( pointCount, outputDimCount ) );
    for (int k = 0; k < outputDimCount; k++) {
        int eigenIndex = sortInd[ k ];
        if (verbose)
            disp( 2, "dim %d: eigen value: %f", k, eigenVals[ eigenIndex ] );
        for (int i = 0; i < pointCount; i++) 
            output->data( i, k ) = eigenVects->data( i, eigenIndex );
    }

    // normalize each output point (row) to have unit length
    for (int i = 0; i < pointCount; i++) {
        float lenSqd = 0;
        for (int j = 0; j < outputDimCount; j++) {
            float val = output->data( i, j );
            lenSqd += val * val;
        }
        if (lenSqd > 1e-5f) {
            float factor = 1.0f / sqrtf( lenSqd );
            for (int j = 0; j < outputDimCount; j++) 
                output->data( i, j ) *= factor;
        }
    }
    return output;
}


//-------------------------------------------
// DIAGNOSTIC COMMANDS
//-------------------------------------------


// test KMeans function on synthetic data
void testKMeans( Config &conf ) {
    int dimCount = 2;
    int pointCount = 1000;
    int clusterCount = 5;

    // generate data
    MatrixF points( pointCount, dimCount );
    for (int i = 0; i < pointCount; i++) {
        points( i, 0 ) = (float) randomFloat();
        points( i, 1 ) = (float) randomFloat();
    }

    // run k-means
    aptr<MatrixF> means;
    aptr<VectorI> assign;
    kMeans( points, clusterCount, means, assign );

    // plot data
    aptr<Plot> plot( new Plot() );
    plot->setXBounds( 0, 1 );
    plot->setYBounds( 0, 1 );
    plot->setStyle( PLOT_DOTS );
    for (int i = 0; i < pointCount; i++) {
        float x = points( i, 0 );
        float y = points( i, 1 );
        int r = 0, g = 0, b = 0;
        int clusterIndex = assign->data( i );
        colorizeDiscrete( clusterIndex, r, g, b );
        plot->setColor( r, g, b );
        plot->add( x, y );
    }
    plot->setStyle( PLOT_CIRCLES );
    for (int i = 0; i < clusterCount; i++) {
        float x = means->data( i, 0 );
        float y = means->data( i, 1 );
        int r = 0, g = 0, b = 0;
        colorizeDiscrete( i, r, g, b );
        plot->setColor( r, g, b );
        plot->add( x, y );
    }
    Plot::disp( plot );
}


//-------------------------------------------
// INIT / CLEANUP
//-------------------------------------------


// register commands, etc. defined in this module
void initKMeans() {
    registerCommand( "testkmeans", testKMeans );
}


} // end namespace sbl

