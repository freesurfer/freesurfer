#ifndef _SBL_K_MEANS_H_
#define _SBL_K_MEANS_H_
#include <sbl/core/Pointer.h>
#include <sbl/math/Matrix.h>
namespace sbl {


/*! \file KMeans.h
    \brief The KMeans module includes in an implementation of the K-means clustering
    algorithm.
*/


// register commands, etc. defined in this module
void initKMeans();


/// run the K-means cluster algorithm; returns the cluster means and assignment of points to clusters
void kMeans( const MatrixF &points, int clusterCount, aptr<MatrixF> &means, aptr<VectorI> &assign );


/// transforms data into space suitable for clustering as described in "On Spectral Cluster: Analysis and an Algorithm" by Ng, Jordan, and Weiss;
/// assumes affinity is symmetric and 0 on diagonal
aptr<MatrixF> spectralTransform( MatrixF &affinity, int outputDimCount, bool verbose );


} // end namespace sbl
#endif // _SBL_K_MEANS_H_

