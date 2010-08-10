/**
 * @file  gcalinearprior.hpp
 * @brief Class to hold a volume of GCA priors in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/10 18:55:39 $
 *    $Revision: 1.1 $
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

    GCAlinearPrior( void ) : xDim(0), yDim(0), zDim(0), labelDim(0),
			     nLabels(),
			     maxLabels(),
			     labels(),
			     priors(),
			     totTraining(),
			     tExhume() {};

    // -----------------------------------------------------

    //! Accessor for nLabels array
    short& voxelLabelCount( const int ix, const int iy, const int iz ) {
      /*!
	The nLabels array holds per voxel data
      */
      unsigned int index = nLabelsIndex(ix,iy,iz);

      return( this->nLabels.at(index) );
    }
    
    //! Const accessor for nLabels array
    short voxelLabelCount( const int ix, const int iy, const int iz ) const {
      /*!
	The nLabels array holds per voxel data
      */
      unsigned int index = nLabelsIndex(ix,iy,iz);

      return( this->nLabels.at(index) );
    }

    // --
    

    //! Accessor for maxLabels array
    short& maxVoxelLabel( const int ix, const int iy, const int iz ) {
      /*!
	The maxLabels array holds per voxel data
      */
      unsigned int idx = maxLabelsIndex(ix,iy,iz);

      return( this->maxLabels.at(idx) );
    }

    //! Const accessor for maxLabels array
    short maxVoxelLabel( const int ix, const int iy, const int iz ) const {
      /*!
	The maxLabels array holds per voxel data
      */
      unsigned int idx = maxLabelsIndex(ix,iy,iz);

      return( this->maxLabels.at(idx) );
    }

    
    // --

    //! Accessor for the labels array
    unsigned short& voxelLabel( const int ix, const int iy, const int iz,
				const int iLabel ) {
      /*!
	The labels array is 4D.
	Each voxel has a 1D array of labels hanging from it,
	the length of which is given by the corresponding entry
	in nLabels (max. held in this->labelDim).
      */
      unsigned int idx = labelsIndex(ix,iy,iz,iLabel);

      return( this->labels.at(idx) );
    }

    //! Const accessor for the labels array
    unsigned short voxelLabel( const int ix, const int iy, const int iz,
			       const int iLabel ) const {
      /*!
	The labels array is 4D.
	Each voxel has a 1D array of labels hanging from it,
	the length of which is given by the corresponding entry
	in nLabels.
      */
      unsigned int idx = labelsIndex(ix,iy,iz,iLabel);
      
      return( this->labels.at(idx) );
    }
    

    // --

    //! Accessor for the priors array
    float& voxelPrior(  const int ix, const int iy, const int iz,
			const int iLabel ) {
      /*!
	The priors array is 4D.
	Each voxel has a 1D array of priors hanging from it,
	the length of which is given by the corresponding entry
	in nLabels (max. held in this->labelDim).
      */
      unsigned int idx = priorsIndex(ix,iy,iz,iLabel);

      return( this->priors.at(idx) );
    }

    //! Const accessor for the priors array
    float voxelPrior(  const int ix, const int iy, const int iz,
		       const int iLabel ) const {
      /*!
	The priors array is 4D.
	Each voxel has a 1D array of priors hanging from it,
	the length of which is given by the corresponding entry
	in nLabels.
      */
      unsigned int idx = priorsIndex(ix,iy,iz,iLabel);

      return( this->priors.at(idx) );
    }

    // --

    //! Accessor for totTraining array
    int& totalTraining( const int ix, const int iy, const int iz ) {
      /*!
	The totTraining array holds per voxel data
      */
      unsigned int idx = totTrainingIndex(ix,iy,iz);

      return( this->totTraining.at(idx) );
    }

    //! Const accessor for totTraining array
    int totalTraining( const int ix, const int iy, const int iz ) const {
      /*!
	The totTraining array holds per voxel data
      */
      unsigned int idx = totTrainingIndex(ix,iy,iz);

      return( this->totTraining.at(idx) );
    }

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

    //! Index computation for nLabels
    unsigned int nLabelsIndex( const int ix,
			       const int iy,
			       const int iz ) const {

      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      unsigned int index;
      index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );

      return( index );
    }



    //! Index computation for maxLabels
    unsigned int maxLabelsIndex( const int ix,
				 const int iy,
				 const int iz ) const {
      
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int index;
      index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );
      
      return( index );
    }



    //! Index computation for labels
    unsigned int labelsIndex( const int ix,
			      const int iy,
			      const int iz,
			      const int iLabel ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iLabel<0) || (iLabel>=this->labelDim) ||
	  (iLabel>=this->voxelLabelCount(ix,iy,iz) ) ) {
	cerr << __FUNCTION__
	     << ": Label index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      unsigned int index;
      index = ix + ( this->xDim *
		     ( iy + ( this->yDim *
			      ( iz + ( this->zDim *
				       iLabel ) )
			      ) )
		     );
      return( index );
    }



    //! Index computation for priors
    unsigned int priorsIndex( const int ix,
			      const int iy,
			      const int iz,
			      const int iLabel ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iLabel<0) || (iLabel>=this->labelDim) ||
	  (iLabel>=this->voxelLabelCount(ix,iy,iz) ) ) {
	cerr << __FUNCTION__
	     << ": Label index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int index;
      index = ix + ( this->xDim *
		     ( iy + ( this->yDim *
			      ( iz + ( this->zDim *
				       iLabel ) )
			      ) )
		     );
      return( index );
    }
    

    //! Index computation for totTraining
    unsigned int totTrainingIndex( const int ix,
				   const int iy,
				   const int iz ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      unsigned int index;
      index = ix + ( this->xDim * ( iy + ( this->yDim * iz ) ) );
      
      return( index );
    }

    // =======================

    // Dimensions
    int xDim;
    int yDim;
    int zDim;
    short labelDim;
    

    //! Count of bytes allocated
    size_t bytes;

    //! Stores nlabels field of GCA_PRIOR
    std::vector<short> nLabels;
    //! Stores max_labels of GCA_PRIOR
    std::vector<short> maxLabels;
    //! Stores labels of GCA_PRIOR (4D, bounded by labelDim and nLabels)
    std::vector<unsigned short> labels;
    //! Stores priors field of GCA_PRIOR (4D, bounded by labelDim and nLabels)
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
