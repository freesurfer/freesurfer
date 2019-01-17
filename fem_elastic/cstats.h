
#ifndef _h_coords_stats_h
#define _h_coords_stats_h

#include "coords.h"
#include "small_matrix.h"

template<int n>
double
coords_statistics( const std::vector<TCoords<double, n> >& vcoords,
                   TCoords<double,n>& mean)
{
  typedef std::vector<TCoords<double,n> > VectorType;

  mean.set(0.0);

  for ( typename VectorType::const_iterator cit = vcoords.begin();
        cit != vcoords.end(); ++cit )
  {
    mean += *cit;
  } // next cit
  mean = mean / (double)vcoords.size();

  SmallMatrix mat(n,1);
  mat.set(.0);

  SmallMatrix varcova(n,n);
  varcova.set(.0);

  for ( typename VectorType::const_iterator cit = vcoords.begin();
        cit != vcoords.end(); ++cit )
  {
    for ( int i=0; i<n; ++i)
      mat(i,0) = (*cit)(i) - mean(i);

    SmallMatrix cov_local = mat * mat.transpose();
    varcova += cov_local;
  }

  // compute the norm of the matrix
  varcova *= 1.0/(double)vcoords.size();

  return varcova.norm();
}

#endif
