/**
 * @brief Class to hold a volume of GCA priors in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

#ifndef GCA_LINEAR_PRIOR
#define GCA_LINEAR_PRIOR

#include <vector>

#include "timer.h"
#include "gca.h"

// Forward declaration
namespace GPU
{
namespace Classes
{
class GCApriorGPU;
}
}


namespace Freesurfer
{


// Forward declaration
class const_GCAprior;



//! Class to hold a 3D volume of GCA priors in linear memory
class GCAlinearPrior
{
public:

  GCAlinearPrior( void ) : xDim(0), yDim(0), zDim(0), n4D(0),
    bytes(),
    offsets4D(),
    maxLabels(),
    labels(),
    priors(),
    totTraining() {};

  // -----------------------------------------------------

  //! Accessor for totTraining array
  inline int& totalTraining( const int ix,
                             const int iy,
                             const int iz )
  {
    /*!
    The totTraining array holds per voxel data
    */
    size_t idx = index3D(ix,iy,iz);

    return( this->totTraining.at(idx) );
  }

  //! Const accessor for totTraining array
  inline int totalTraining( const int ix,
                            const int iy,
                            const int iz ) const
  {
    /*!
    The totTraining array holds per voxel data
    */
    size_t idx = index3D(ix,iy,iz);

    return( this->totTraining.at(idx) );
  }


  // --


  //! Accessor for maxLabels array
  inline short& maxVoxelLabel( const int ix,
                               const int iy,
                               const int iz )
  {
    /*!
    The maxLabels array holds per voxel data
    */
    size_t idx = index3D(ix,iy,iz);

    return( this->maxLabels.at(idx) );
  }

  //! Const accessor for maxLabels array
  inline short maxVoxelLabel( const int ix,
                              const int iy,
                              const int iz ) const
  {
    /*!
    The maxLabels array holds per voxel data
    */
    size_t idx = index3D(ix,iy,iz);

    return( this->maxLabels.at(idx) );
  }


  // --

  //! Dynamically compute nlabels for each prior voxel
  inline short voxelLabelCount( const int ix,
                                const int iy,
                                const int iz ) const
  {
    /*!
    This is computed as a difference between consecutive
    entries on the offsets4D array;
    */
    long long currOffset = this->offsets4D.at( this->index3D(ix,iy,iz) );
    long long nextOffset = this->offsets4D.at( this->index3D(ix,iy,iz) + 1 );

    return( nextOffset - currOffset );
  }

  // --

  //! Accessor for the labels array
  inline unsigned short& voxelLabel( const int ix,
                                     const int iy,
                                     const int iz,
                                     const int iLabel )
  {
    /*!
    The labels array is 4D.
    Each voxel has a 1D array of labels hanging from it,
    which are indexed according to the offsets4D array.
    */
    size_t idx = this->index4D(ix,iy,iz,iLabel);

    return( this->labels.at(idx) );
  }

  //! Const accessor for the labels array
  inline unsigned short voxelLabel( const int ix,
                                    const int iy,
                                    const int iz,
                                    const int iLabel ) const
  {
    /*!
    The labels array is 4D.
    Each voxel has a 1D array of labels hanging from it,
    which are indexed according to the offsets4D array.

    */
    size_t idx = this->index4D(ix,iy,iz,iLabel);

    return( this->labels.at(idx) );
  }


  // --

  //! Accessor for the priors array
  inline float& voxelPrior(  const int ix, const int iy, const int iz,
                             const int iLabel )
  {
    /*!
    The priors array is 4D.
    Each voxel has a 1D array of priors hanging from it,
    which are indexed according to the offsets4D array.
    */
    size_t idx = this->index4D(ix,iy,iz,iLabel);

    return( this->priors.at(idx) );
  }

  //! Const accessor for the priors array
  inline float voxelPrior(  const int ix, const int iy, const int iz,
                            const int iLabel ) const
  {
    /*!
    The priors array is 4D.
    Each voxel has a 1D array of priors hanging from it,
    which are indexed according to the offsets4D array.
    */
    size_t idx = this->index4D(ix,iy,iz,iLabel);

    return( this->priors.at(idx) );
  }

  // --



  // -------------------------

  void Exhume( const GCA* const src );
  void Inhume( GCA* targ ) const;

  // -------------------------
  void PrintStats( std::ostream& os = std::cout ) const;

  const_GCAprior GetConstPrior( const int ix,
                                const int iy,
                                const int iz ) const;


  // -----------------------------------------------------
private:

  //! Method to find required array sizes
  void ExtractDims( const GCA* const src );

  //! Method to allocate memory for the priors
  void Allocate( void );

  //! Method to remove priors from GCA
  void ScorchPriors( GCA* targ ) const;

  // -----------------------------------------------------

  //! Index computation for 3D indices
  inline size_t index3D( const int ix,
                         const int iy,
                         const int iz ) const
  {
    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    size_t index;
    index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );

    return( index );
  }

  //! Index computation for 4D indices
  inline size_t index4D( const int ix,
                         const int iy,
                         const int iz,
                         const int iLabel ) const
  {
    if ( (ix<0) || (ix>=this->xDim) ||
         (iy<0) || (iy>=this->yDim) ||
         (iz<0) || (iz>=this->zDim ) )
    {
      std::cerr << __FUNCTION__
                << ": Index out of range" << std::endl;
      abort();
    }

    const size_t idx3D = this->index3D(ix,iy,iz);

    const long long currOffset = this->offsets4D.at( idx3D );
    const long long nextOffset = this->offsets4D.at( idx3D + 1 );

    if ( (iLabel<0) || (iLabel>=(nextOffset-currOffset) ) )
    {
      std::cerr << __FUNCTION__
                << ": iLabel out of range" << std::endl;
      abort();
    }

    return( currOffset+iLabel );
  }



  // =======================

  // Dimensions
  long long xDim;
  long long yDim;
  long long zDim;
  size_t n4D;

  //! Count of bytes allocated
  size_t bytes;

  //! Stores offsets of the (variable length) 4th dimensions
  std::vector<size_t> offsets4D;

  //! Stores max_labels of GCA_PRIOR
  std::vector<short> maxLabels;
  //! Stores labels of GCA_PRIOR (4D, variable length 4th dimension)
  std::vector<unsigned short> labels;
  //! Stores priors field of GCA_PRIOR (4D, variable length 4th dimension)
  std::vector<float> priors;
  //! Stores total_training field of GCA_PRIOR
  std::vector<int> totTraining;

  //! Time for exhumation of data
  mutable long exhumeTime;
  //! Time for inhumation of data
  mutable long inhumeTime;

  friend class const_GCAprior;
  friend class GPU::Classes::GCApriorGPU;
};




//! Equivalent of 'const GCA_PRIOR'
/*!
  This class provides quick access to a particular 3D point
  of a GCAlinearPrior.
  The indexing must match that used in GCAlinearPrior, or
  everything will fall apart horribly.
  A lot of the precomputations are performed in the constructor
*/
class const_GCAprior
{
public:
  const_GCAprior( const int _ix,
                  const int _iy,
                  const int _iz,
                  const GCAlinearPrior& src ) : gcalp(src),
    idx3d(src.index3D(_ix,_iy,_iz)),
    currOffset(src.offsets4D.at(idx3d)),
    myLabelCount(src.offsets4D.at(idx3d+1)-currOffset) {}


  //! Accessor for totalTraining
  inline int totalTraining( void ) const
  {
    return( this->gcalp.totTraining.at(this->idx3d) );
  }

  //! Accessor for maxLabels
  inline short maxLabel( void ) const
  {
    return( this->gcalp.maxLabels.at(this->idx3d) );
  }

  //! Accessor for nlabels
  inline short labelCount( void ) const
  {
    return( this->myLabelCount );
  }

  //! Accessor for labels
  inline unsigned short labels( const int iLabel ) const
  {
    if ( (iLabel<0) || (iLabel>=this->myLabelCount) )
    {
      std::cerr << __FUNCTION__
                << ": iLabel out of range " << iLabel
                << std::endl;
      abort();
    }

    return( this->gcalp.labels.at( this->currOffset+iLabel ) );
  }

  //! Accessor for priors
  inline float priors( const int iLabel ) const
  {
    if ( (iLabel<0) || (iLabel>=this->myLabelCount) )
    {
      std::cerr << __FUNCTION__
                << ": iLabel out of range " << iLabel
                << std::endl;
      abort();
    }

    return( this->gcalp.priors.at( this->currOffset+iLabel ) );
  }


private:

  //! The GCAlinearPrior we're part of
  const GCAlinearPrior& gcalp;
  //! Precomputed linear index for 3D data
  const size_t idx3d;
  //! Precomputed linear start index for 4D data
  const size_t currOffset;
  //! Length of the 4th dimension
  const int myLabelCount;
};



}

#endif
