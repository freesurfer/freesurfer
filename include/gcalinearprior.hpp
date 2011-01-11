/**
 * @file  gcalinearprior.hpp
 * @brief Class to hold a volume of GCA priors in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/01/11 18:24:21 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2010,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 *
 */

#ifndef GCA_LINEAR_PRIOR
#define GCA_LINEAR_PRIOR

#include <vector>

#include "chronometer.hpp"

#include "gca.h"



namespace Freesurfer {

  //! Class to hold a 3D volume of GCA priors in linear memory
  class GCAlinearPrior {
  public:

    GCAlinearPrior( void ) : xDim(0), yDim(0), zDim(0), n4D(0),
			     bytes(),
			     offsets4D(),
			     maxLabels(),
			     labels(),
			     priors(),
			     totTraining(),
			     tExhume() {};

    // -----------------------------------------------------

    //! Accessor for totTraining array
    inline int& totalTraining( const int ix,
			       const int iy,
			       const int iz ) {
      /*!
	The totTraining array holds per voxel data
      */
      unsigned int idx = index3D(ix,iy,iz);
      
      return( this->totTraining.at(idx) );
    }
    
    //! Const accessor for totTraining array
    inline int totalTraining( const int ix,
			      const int iy,
			      const int iz ) const {
      /*!
	The totTraining array holds per voxel data
      */
      unsigned int idx = index3D(ix,iy,iz);
      
      return( this->totTraining.at(idx) );
    }


    // --
    

    //! Accessor for maxLabels array
    inline short& maxVoxelLabel( const int ix,
				 const int iy,
				 const int iz ) {
      /*!
	The maxLabels array holds per voxel data
      */
      unsigned int idx = index3D(ix,iy,iz);

      return( this->maxLabels.at(idx) );
    }

    //! Const accessor for maxLabels array
    inline short maxVoxelLabel( const int ix,
				const int iy,
				const int iz ) const {
      /*!
	The maxLabels array holds per voxel data
      */
      unsigned int idx = index3D(ix,iy,iz);

      return( this->maxLabels.at(idx) );
    }

    
    // --
    
    //! Dynamically compute nlabels for each prior voxel
    inline short voxelLabelCount( const int ix,
				  const int iy,
				  const int iz ) const {
      /*!
	This is computed as a difference between consecutive
	entries on the offsets4D array;
      */
      short currOffset = this->offsets4D.at( this->index3D(ix,iy,iz) );
      short nextOffset = this->offsets4D.at( this->index3D(ix,iy,iz) + 1 );

      return( nextOffset - currOffset );
    }

    // --

    //! Accessor for the labels array
    inline unsigned short& voxelLabel( const int ix,
				       const int iy,
				       const int iz,
				       const int iLabel ) {
      /*!
	The labels array is 4D.
	Each voxel has a 1D array of labels hanging from it,
	which are indexed according to the offsets4D array.
      */
      unsigned int idx = this->index4D(ix,iy,iz,iLabel);

      return( this->labels.at(idx) );
    }

    //! Const accessor for the labels array
    inline unsigned short voxelLabel( const int ix,
				      const int iy,
				      const int iz,
				      const int iLabel ) const {
      /*!
	The labels array is 4D.
	Each voxel has a 1D array of labels hanging from it,
	which are indexed according to the offsets4D array.
	
      */
      unsigned int idx = this->index4D(ix,iy,iz,iLabel);
      
      return( this->labels.at(idx) );
    }
    

    // --

    //! Accessor for the priors array
    inline float& voxelPrior(  const int ix, const int iy, const int iz,
			       const int iLabel ) {
      /*!
	The priors array is 4D.
	Each voxel has a 1D array of priors hanging from it,
	which are indexed according to the offsets4D array.
      */
      unsigned int idx = this->index4D(ix,iy,iz,iLabel);

      return( this->priors.at(idx) );
    }

    //! Const accessor for the priors array
    inline float voxelPrior(  const int ix, const int iy, const int iz,
			      const int iLabel ) const {
      /*!
	The priors array is 4D.
	Each voxel has a 1D array of priors hanging from it,
	which are indexed according to the offsets4D array.
      */
      unsigned int idx = this->index4D(ix,iy,iz,iLabel);

      return( this->priors.at(idx) );
    }

    // --

    

    // -------------------------

    void Exhume( const GCA* const src );
    void Inhume( GCA* targ ) const;

    // -------------------------
    void PrintStats( std::ostream& os = std::cout ) const;

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
    inline unsigned int index3D( const int ix,
				 const int iy,
				 const int iz ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	abort();
      }

      unsigned int index;
      index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );

      return( index );
    }

    //! Index computation for 4D indices
    inline unsigned int index4D( const int ix,
				 const int iy,
				 const int iz,
				 const int iLabel ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	abort();
      }

      const unsigned int idx3D = this->index3D(ix,iy,iz);

      const int currOffset = this->offsets4D.at( idx3D );
      const int nextOffset = this->offsets4D.at( idx3D + 1 );

      if( (iLabel<0) || (iLabel>=(nextOffset-currOffset) ) ) {
	cerr << __FUNCTION__
	     << ": iLabel out of range" << endl;
	abort();
      }
      
      return( currOffset+iLabel );
    }



    // =======================

    // Dimensions
    int xDim;
    int yDim;
    int zDim;
    unsigned int n4D;

    //! Count of bytes allocated
    size_t bytes;

    //! Stores offsets of the (variable length) 4th dimensions
    std::vector<unsigned int> offsets4D;

    //! Stores max_labels of GCA_PRIOR
    std::vector<short> maxLabels;
    //! Stores labels of GCA_PRIOR (4D, variable length 4th dimension)
    std::vector<unsigned short> labels;
    //! Stores priors field of GCA_PRIOR (4D, variable length 4th dimension)
    std::vector<float> priors;
    //! Stores total_training field of GCA_PRIOR
    std::vector<int> totTraining;

    //! Time required for exhumation of data
    mutable SciGPU::Utilities::Chronometer tExhume;
    //! Time require for inhumation of data
    mutable SciGPU::Utilities::Chronometer tInhume;
  };

}

#endif
