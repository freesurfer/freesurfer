/**
 * @file  gcanode.hpp
 * @brief Class to hold a volume of GCA nodes in linear memory on the GPU
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:23 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef GCA_NODE_GPU_HPP
#define GCA_NODE_GPU_HPP

#include <iostream>

#include "gca.h"
#include "gcalinearnode.hpp"


namespace GPU
{

namespace Classes
{

// Forward declarations
class const_GCAnode;
class const_GCAnode_GC1D;


//! Device class for GCA node data
class GCAnodeGPUarg
{
public:

  //! Constructor from contents
  GCAnodeGPUarg( const long _xDim, const long _yDim, const long _zDim,
                 const size_t _n4D, const size_t _n5D, const size_t _n6D,
                 size_t *const _offsets4D,
                 size_t *const _offsets5D,
                 size_t *const _offsets6D,
                 int *const _nodeMaxLabels,
                 int *const _nodeTotalTraining,
                 unsigned short *const _nodeLabels,
                 float *const _means,
                 float *const _variances,
                 short *const _nJustPriors,
                 int *const _nTraining,
                 char *const _regularised,
                 float *const _gc1dDirecLabelPriors,
                 unsigned short *const _gc1dDirecLabels ) :
    xDim(_xDim), yDim(_yDim), zDim(_zDim),
    n4D(_n4D), n5D(_n5D), n6D(_n6D),
    gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
    offsets4D(_offsets4D), offsets5D(_offsets5D), offsets6D(_offsets6D),
    nodeMaxLabels(_nodeMaxLabels), nodeTotalTraining(_nodeTotalTraining),
    nodeLabels(_nodeLabels),
    means(_means), variances(_variances),
    nJustPriors(_nJustPriors),
    nTraining(_nTraining), regularised(_regularised),
    gc1dDirecLabelPriors(_gc1dDirecLabelPriors),
    gc1dDirecLabels(_gc1dDirecLabels) {};

  // ---------------------------------------------
  // 3D data access

  //! Accessor for nodeMaxLabels
  __device__
  int maxLabels( const int ix, const int iy, const int iz ) const
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeMaxLabels[idx] );
  }


  //! Accessor for nodeTotalTraining
  __device__
  int totalTraining( const int ix,
                     const int iy,
                     const int iz ) const
  {
    const size_t idx = this->index3D(ix,iy,iz);

    return( this->nodeTotalTraining[idx] );
  }


  //! Dynamic computation of gc1dCount
  __device__
  int gc1dCount( const int ix,
                 const int iy,
                 const int iz ) const
  {
    const size_t idx3D = this->index3D(ix,iy,iz);
    const size_t currOffset = this->offsets4D[idx3D];
    const size_t nextOffset = this->offsets4D[idx3D+1];

    return( nextOffset - currOffset );
  }

  // ---------------------------------------------
  // 4D data access

  //! Access for nodeLabels
  __device__
  unsigned short labelsAtNode( const int ix,
                               const int iy,
                               const int iz,
                               const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nodeLabels[idx] );
  }


  //! Accessor for means
  __device__
  float meansAtNodeGC1D( const int ix,
                         const int iy,
                         const int iz,
                         const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->means[idx] );
  }


  //! Accessor for variances
  __device__
  float variancesAtNodeGC1D( const int ix,
                             const int iy,
                             const int iz,
                             const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->variances[idx] );
  }


  //! Accessor for nJustPriors
  __device__
  short nJustPriorsAtNodeGC1D( const int ix,
                               const int iy,
                               const int iz,
                               const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nJustPriors[idx] );
  }


  //! Accessor for nTraining
  __device__
  int nTrainingAtNodeGC1D( const int ix,
                           const int iy,
                           const int iz,
                           const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->nTraining[idx] );
  }


  //! Accessor for regularised
  __device__
  char regularisedAtNodeGC1D( const int ix,
                              const int iy,
                              const int iz,
                              const int iGC1D ) const
  {
    const size_t idx = this->index4D(ix,iy,iz,iGC1D);

    return( this->regularised[idx] );
  }

  // ---------------------------------------------
  // 5D data access

  //! Dynamic computation of nLabelsAtNodeGC1Ddirection
  __device__
  int nLabelsAtNodeGC1Ddirection( const int ix,
                                  const int iy,
                                  const int iz,
                                  const int iGC1D,
                                  const int iDirec ) const
  {
    /*!
      This is computed as the difference between
      consecutive entries of the offsets6D array
    */
    const size_t idx5D = this->index5D(ix,iy,iz,iGC1D,iDirec);

    const size_t currOffset = this->offsets6D[idx5D];
    const size_t nextOffset = this->offsets6D[idx5D+1];

    return( nextOffset - currOffset );
  }

  // ---------------------------------------------
  // 6D data access

  //! Accesor for gc1dDirecLabelPriors
  __device__
  float labelPriorsAtNodeGC1Ddirection( const int ix,
                                        const int iy,
                                        const int iz,
                                        const int iGC1D,
                                        const int iDirec,
                                        const int iLabel ) const
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabelPriors[idx] );
  }


  //! Accessor for gc1dDirecLabels
  __device__
  unsigned short labelsAtNodeGC1Ddirection( const int ix,
      const int iy,
      const int iz,
      const int iGC1D,
      const int iDirec,
      const int iLabel ) const
  {
    const size_t idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

    return( this->gc1dDirecLabels[idx] );
  }

#ifdef __CUDA_ARCH__
  __device__
  const_GCAnode GetConstNode( const int ix,
                              const int iy,
                              const int iz ) const;
#else
#warning "Not declaring GetConstNode method on host"
#endif

private:
  // Dimensions
  const long xDim, yDim, zDim;
  const size_t n4D, n5D, n6D;
  const unsigned int gc1dNeighbourDim;

  // Offset arrays
  size_t *const offsets4D;
  size_t *const offsets5D;
  size_t *const offsets6D;

  // 3D arrays
  int *const nodeMaxLabels;
  int *const nodeTotalTraining;

  // 4D arrays
  unsigned short *const nodeLabels;
  float *const means;
  float *const variances;
  short *const nJustPriors;
  int *const nTraining;
  char *const regularised;

  // 6D arrays
  float *const gc1dDirecLabelPriors;
  unsigned short *const gc1dDirecLabels;

  // -------------------------------------------

  //! Index computation for 3D arrays
  __device__
  size_t index3D( const int ix,
                  const int iy,
                  const int iz ) const
  {
    size_t index;
    index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );

    return( index );
  }

  //! Index computation for 4D arrays
  __device__
  size_t index4D( const int ix,
                  const int iy,
                  const int iz,
                  const int iGC1D ) const
  {
    const size_t idx3D = this->index3D(ix,iy,iz);
    const size_t currOffset = this->offsets4D[ idx3D ];

    return( currOffset + iGC1D );
  }

  //! Index computation for 5D arrays
  __device__
  size_t index5D( const int ix,
                  const int iy,
                  const int iz,
                  const int iGC1D,
                  const int iDirec ) const
  {
    const size_t idx3D = this->index3D(ix,iy,iz);
    const size_t currGC1dOffset = this->offsets4D[ idx3D ];
    const size_t currOffset = this->offsets5D[ currGC1dOffset+iGC1D ];

    return( currOffset+iDirec );
  }

  //! Index computation for 6D arrays
  __device__
  size_t index6D( const int ix,
                  const int iy,
                  const int iz,
                  const int iGC1D,
                  const int iDirec,
                  const int iLabel ) const
  {
    const size_t idx3D = this->index3D(ix,iy,iz);
    const size_t currGC1dOffset = this->offsets4D[ idx3D ];
    const size_t currDirecOffset = this->offsets5D[ currGC1dOffset+iGC1D ];
    const size_t currOffset = this->offsets6D[ currDirecOffset+iDirec ];

    return( currOffset+iLabel );
  }

  friend const_GCAnode;
  friend const_GCAnode_GC1D;
};


// ===============

#ifdef __CUDA_ARCH__
//! Mirror of the CPU const_GCAnode
class const_GCAnode
{
public:
  //! Constructor from location
  __device__
  const_GCAnode( const int _ix,
                 const int _iy,
                 const int _iz,
                 const GCAnodeGPUarg& src ) : gnga(src),
    idx3d(src.index3D(_ix,_iy,_iz)),
    offset4d(src.offsets4D[idx3d]),
    myGC1Dcount(src.offsets4D[idx3d+1]-offset4d),
    ix(_ix),
    iy(_iy),
    iz(_iz) {}


  //! Accessor for max_labels
  __device__
  int maxLabels( void ) const
  {
    return( this->gnga.nodeMaxLabels[this->idx3d] );
  }

  //! Accessor for totalTraining
  __device__
  int totalTraining( void ) const
  {
    return( this->gnga.nodeTotalTraining[this->idx3d] );
  }

  //! Accessor for gc1dCount
  __device__
  int gc1dCount( void ) const
  {
    return( this->myGC1Dcount );
  }

  //! Access for nodeLabels
  __device__
  unsigned short labels( const int iGC1D ) const
  {
    return( this->gnga.nodeLabels[this->offset4d+iGC1D] );
  }

  //! Create the GC1D
  const_GCAnode_GC1D GetConstGC1D( const int iGC1D ) const;

private:
  const GCAnodeGPUarg& gnga;
  const size_t idx3d;
  const size_t offset4d;
  const int myGC1Dcount;

  const int ix, iy, iz;

};


__device__
const_GCAnode GCAnodeGPUarg::GetConstNode( const int ix,
    const int iy,
    const int iz ) const
{
  return( const_GCAnode( ix, iy, iz, *this ) );
}
#else
#warning "const_GCAnode does not exist in host code"
#endif


// ===============

#ifdef __CUDA_ARCH__
class const_GCAnode_GC1D
{
public:
  //! Construction from location
  __device__
  const_GCAnode_GC1D( const int _ix,
                      const int _iy,
                      const int _iz,
                      const int _iGC1D,
                      const GCAnodeGPUarg& src ) : gnga(src),
    offset4d(src.index4D(_ix,_iy,_iz,_iGC1D)) {}

  // -----------------------------
  // 4D data

  //! Accessor for the mean
  __device__
  float mean( void ) const
  {
    return( this->gnga.means[this->offset4d] );
  }

  //! Accessor for the variance
  __device__
  float variance( void ) const
  {
    return( this->gnga.variances[this->offset4d] );
  }

  //! Accessor for nJustPriors
  __device__
  short nJustPriors( void ) const
  {
    return( this->gnga.nJustPriors[this->offset4d] );
  }

  //! Accessor for nTraining
  __device__
  int nTraining( void ) const
  {
    return( this->gnga.nTraining[this->offset4d] );
  }

  //! Accessor for regularised
  __device__
  int regularised( void ) const
  {
    return( this->gnga.regularised[this->offset4d] );
  }

  // -----------------------------
  // 5D data
  __device__
  short nLabels( const int iDirec ) const
  {
    /*!
      This is basically an alternative for
      GCAlinearNode::nLabelsAtNodeGC1Ddirection
    */
    const size_t idx5D = this->gnga.offsets5D[this->offset4d] + iDirec;

    const size_t currOffset = this->gnga.offsets6D[idx5D];
    const size_t nextOffset = this->gnga.offsets6D[idx5D+1];

    return( nextOffset-currOffset );
  }


  // -----------------------------
  // 6D data

  //! Accessor for labels
  __device__
  unsigned short labels( const int iDirec, const int iLabel ) const
  {
    const size_t idx6D = this->index6D(iDirec,iLabel);
    return( this->gnga.gc1dDirecLabels[idx6D] );
  }

  //! Accessor for labelPriors
  __device__
  float labelPriors( const int iDirec, const int iLabel ) const
  {
    const size_t idx6D = this->index6D(iDirec,iLabel);
    return( this->gnga.gc1dDirecLabelPriors[idx6D] );
  }

private:
  const GCAnodeGPUarg& gnga;
  const size_t offset4d;

  //! Index computation for 6D data
  __device__
  size_t index6D( const int iDirec,
                  const int iLabel ) const
  {
    const size_t idx5D = this->gnga.offsets5D[this->offset4d] + iDirec;
    const size_t currOffset = this->gnga.offsets6D[idx5D];
    return( currOffset + iLabel );
  }
};



__host__ __device__
const_GCAnode_GC1D const_GCAnode::GetConstGC1D( const int iGC1D ) const
{
  const_GCAnode_GC1D gc1d( this->ix,
                           this->iy,
                           this->iz,
                           iGC1D,
                           this->gnga );
  return( gc1d );
}
#else
#warning "const_GCAnode_GC1D does not exsist in host code"
#endif




// =============================================================

//! Management class for GCA node data on the GPU
class GCAnodeGPU
{
public:
  //! Default constructor
  GCAnodeGPU( void ) : xDim(0), yDim(0), zDim(0),
    n4D(0), n5D(0), n6D(0),
    gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
    d_offsets4D(NULL), d_offsets5D(NULL), d_offsets6D(NULL),
    d_nodeMaxLabels(NULL), d_nodeTotalTraining(NULL),
    d_nodeLabels(NULL),
    d_means(NULL), d_variances(NULL),
    d_nJustPriors(NULL), d_nTraining(NULL),
    d_regularised(NULL),
    d_gc1dDirecLabelPriors(NULL),
    d_gc1dDirecLabels(NULL) {};

  //! Destructor
  ~GCAnodeGPU( void )
  {
    this->Release();
  }

  //! Conversion operator
  operator GCAnodeGPUarg( void ) const
  {
    GCAnodeGPUarg gnga( this->xDim, this->yDim, this->zDim,
                        this->n4D, this->n5D, this->n6D,
                        this->d_offsets4D,
                        this->d_offsets5D,
                        this->d_offsets6D,
                        this->d_nodeMaxLabels,
                        this->d_nodeTotalTraining,
                        this->d_nodeLabels,
                        this->d_means,
                        this->d_variances,
                        this->d_nJustPriors,
                        this->d_nTraining,
                        this->d_regularised,
                        this->d_gc1dDirecLabelPriors,
                        this->d_gc1dDirecLabels );

    return( gnga );
  }


  // -------------------
  // Memory management

  //! Allocate all arrays
  void Allocate( const int nxDim, const int nyDim, const int nzDim,
                 const size_t num4D, const size_t num5D, const size_t num6D );

  //! Release all arrays
  void Release( void );

  // -------------------
  // Transfers

  //! Send data to the GPU
  void Send( const Freesurfer::GCAlinearNode& src );


  // -------------------------------------------
private:
  // Dimensions
  long long xDim;
  long long yDim;
  long long zDim;
  size_t n4D;
  size_t n5D;
  size_t n6D;
  const unsigned int gc1dNeighbourDim;

  // Offsets
  size_t *d_offsets4D;
  size_t *d_offsets5D;
  size_t *d_offsets6D;

  // 3D arrays
  int *d_nodeMaxLabels;
  int *d_nodeTotalTraining;

  // 4D arrays
  unsigned short *d_nodeLabels;
  float *d_means;
  float *d_variances;
  short *d_nJustPriors;
  int *d_nTraining;
  char *d_regularised;

  // 6D arrays
  float *d_gc1dDirecLabelPriors;
  unsigned short *d_gc1dDirecLabels;


  // -----------------------------------

  //! Inhibit copy constructor
  GCAnodeGPU( const GCAnodeGPU& src ): xDim(0), yDim(0), zDim(0),
    n4D(0), n5D(0), n6D(0),
    gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
    d_offsets4D(NULL), d_offsets5D(NULL), d_offsets6D(NULL),
    d_nodeMaxLabels(NULL), d_nodeTotalTraining(NULL),
    d_nodeLabels(NULL),
    d_means(NULL), d_variances(NULL),
    d_nJustPriors(NULL), d_nTraining(NULL),
    d_regularised(NULL),
    d_gc1dDirecLabelPriors(NULL),
    d_gc1dDirecLabels(NULL)
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy"
              << std::endl;
    abort();
  }

  //! Inhibit assignment
  GCAnodeGPU& operator=( const GCAnodeGPU& src )
  {
    std::cerr << __FUNCTION__
              << ": Please don't copy"
              << std::endl;
    abort();
  }

};

}
}




#endif
