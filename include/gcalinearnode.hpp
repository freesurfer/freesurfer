/**
 * @file  gcalinearnode.hpp
 * @brief Class to hold a volume of GCA nodes in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/01/11 20:56:18 $
 *    $Revision: 1.5 $
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

#ifndef GCA_LINEAR_NODE_HPP
#define GCA_LINEAR_NODE_HPP

#include <vector>

#include "chronometer.hpp"

#include "gca.h"


//! Catch-all namespace
namespace Freesurfer {

  //! Class to hold a 3D volume of GCA nodes in linear memory (ninputs=1)
  class GCAlinearNode {
  public:
    GCAlinearNode( void ) : xDim(0), yDim(0), zDim(0),
			    n4D(0), n5D(0), n6D(0),
			    gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
			    offsets4D(), offsets5D(), offsets6D(),
			    nodeMaxLabels(),
			    nodeTotalTraining(),
			    nodeLabels(),
			    means(), variances(),
			    nJustPriors(), nTraining(), regularised(),
			    gc1dDirecLabelPriors(),
			    gc1dDirecLabels(),
			    tExhume(), tInhume() {};

   
   

    // -------------------------------------------------
    // 3D data access

    //! Accessor for nodeMaxLabels
    inline int maxLabels( const int ix, const int iy, const int iz ) const {
      unsigned int idx = this->index3D(ix,iy,iz);

      return( this->nodeMaxLabels.at(idx) );
    }

    //! Mutator for nodeMaxLabels
    inline int& maxLabels( const int ix, const int iy, const int iz ) {
      unsigned int idx = this->index3D(ix,iy,iz);

      return( this->nodeMaxLabels.at(idx) );
    }

    // --------------

    //! Accessor for nodeTotalTraining
    inline int totalTraining( const int ix,
			      const int iy,
			      const int iz ) const {
      unsigned int idx = this->index3D(ix,iy,iz);

      return( this->nodeTotalTraining.at(idx) );
    }

    //! Mutator for nodeTotalTraining
    inline int& totalTraining( const int ix,
			       const int iy,
			       const int iz ) {
      unsigned int idx = this->index3D(ix,iy,iz);
      
      return( this->nodeTotalTraining.at(idx) );
    }

    // -------------------------------------------------


    //! Dynamic computation of gc1dCount
    inline int gc1dCount( const int ix,
			  const int iy,
			  const int iz ) const {
      /*!
	This is computed as the difference between
	consecutive entries of the offsets4D array
      */
      const int idx3D = this->index3D(ix,iy,iz);

      int currOffset = this->offsets4D.at( idx3D );
      int nextOffset = this->offsets4D.at( idx3D+1 );

      return( nextOffset - currOffset );
    }


    // -------------------------------------------------
    // 4D data access
    
    //! Accessor for nodeLabels
    inline unsigned short labelsAtNode( const int ix,
					const int iy,
					const int iz,
					const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->nodeLabels.at(idx) );
    }

    //! Mutator for nodeLabels
    inline unsigned short& labelsAtNode( const int ix,
					 const int iy,
					 const int iz,
					 const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->nodeLabels.at(idx) );
    }

    // --

    //! Accessor for means
    inline float meansAtNodeGC1D( const int ix,
				  const int iy,
				  const int iz,
				  const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->means.at(idx) );
    }

    //! Mutator for means
    inline float& meansAtNodeGC1D( const int ix,
				   const int iy,
				   const int iz,
				   const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->means.at(idx) );
    }

    // --

    //! Accessor for variances
    inline float variancesAtNodeGC1D( const int ix,
				      const int iy,
				      const int iz,
				      const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->variances.at(idx) );
    }

    //! Mutator for variances
    inline float& variancesAtNodeGC1D( const int ix,
				       const int iy,
				       const int iz,
				       const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->variances.at(idx) );
    }

    // --

    //! Accessor for nJustPriors
    inline short nJustPriorsAtNodeGC1D( const int ix,
					const int iy,
					const int iz,
					const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->nJustPriors.at(idx) );
    }

    //! Mutator for nJustPriors
    inline short& nJustPriorsAtNodeGC1D( const int ix,
					 const int iy,
					 const int iz,
					 const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->nJustPriors.at(idx) );
    }

    // --

    //! Accessor for nTraining
    inline int nTrainingAtNodeGC1D( const int ix,
				    const int iy,
				    const int iz,
				    const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->nTraining.at(idx) );
    }

    //! Mutator for nTraining
    inline int& nTrainingAtNodeGC1D( const int ix,
				     const int iy,
				     const int iz,
				     const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);
      
      return( this->nTraining.at(idx) );
    }

    // --

    //! Accessor for regularised
    inline char regularisedAtNodeGC1D( const int ix,
				       const int iy,
				       const int iz,
				       const int iGC1D ) const {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);

      return( this->regularised.at(idx) );
    }

    //! Mutator for regularised
    inline char& regularisedAtNodeGC1D( const int ix,
					const int iy,
					const int iz,
					const int iGC1D ) {
      const unsigned int idx = this->index4D(ix,iy,iz,iGC1D);
      
      return( this->regularised.at(idx) );
    }

    // -------------------------------------------------


    //! Dynamic computation of nLabelsAtNodeGC1Ddirection
    inline int nLabelsAtNodeGC1Ddirection( const int ix,
					   const int iy,
					   const int iz,
					   const int iGC1D,
					   const int iDirec ) const {
      /*!
	This is computed as the difference between
	consecutive entries of the offsets6D array
      */
      const unsigned int idx5D = this->index5D(ix,iy,iz,iGC1D,iDirec);

      int currOffset = this->offsets6D.at( idx5D );
      int nextOffset = this->offsets6D.at( idx5D+1 );

      return( nextOffset - currOffset );
    }

    // -------------------------------------------------
    // 6D data access


    //! Accesor for gc1dDirecLabelPriors
    inline float labelPriorsAtNodeGC1Ddirection( const int ix,
						 const int iy,
						 const int iz,
						 const int iGC1D,
						 const int iDirec,
						 const int iLabel ) const {
      const unsigned int idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

      return( this->gc1dDirecLabelPriors.at(idx) );
    }

    //! Mutator for gc1dDirecLabelPriors
    inline float& labelPriorsAtNodeGC1Ddirection( const int ix,
						  const int iy,
						  const int iz,
						  const int iGC1D,
						  const int iDirec,
						  const int iLabel ) {
      const unsigned int idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);

      return( this->gc1dDirecLabelPriors.at(idx) );
    }

    // --

    //! Accessor for gc1dDirecLabels
    inline unsigned short labelsAtNodeGC1Ddirection( const int ix,
						     const int iy,
						     const int iz,
						     const int iGC1D,
						     const int iDirec,
						     const int iLabel ) const {
      const unsigned int idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);
      
      return( this->gc1dDirecLabels.at(idx) );
    }
    
    //! Mutator for gc1dDirecLabels
    inline unsigned short& labelsAtNodeGC1Ddirection( const int ix,
						      const int iy,
						      const int iz,
						      const int iGC1D,
						      const int iDirec,
						      const int iLabel ) {
      const unsigned int idx = this->index6D(ix,iy,iz,iGC1D,iDirec,iLabel);
      
      return( this->gc1dDirecLabels.at(idx) );
    }

    // -------------------------------------------------

    //! Method to extract node data from a GCA
    void Exhume( const GCA* const src );

    //! Method to place node data into a GCA
    void Inhume( GCA* dst ) const;

    //! Method to remove nodes from a GCA prior to inhumation of new data
    void ScorchNodes( GCA* targ ) const;


    // -------------------------------------------------
  private:

    // Dimensions
    int xDim;
    int yDim;
    int zDim;
    unsigned int n4D;
    unsigned int n5D;
    unsigned int n6D;
    const int gc1dNeighbourDim;
    

    //! Stores offsets of (variable length) 4th dimension
    std::vector<unsigned int> offsets4D;

    //! Stores offsets of 5th dimension (5th dim itself always of length gc1dNeighbourDim)
    std::vector<unsigned int> offsets5D;

    //! Stores offsets of (variable length) 6th dimension
    std::vector<unsigned int> offsets6D;


    //! Stores the value of max_labels for each node (3D)
    std::vector<int> nodeMaxLabels;
    //! Stores value of total_training for each node (3D)
    std::vector<int> nodeTotalTraining;

    
    //! Stores a list of labels for each node (4D)
    std::vector<unsigned short> nodeLabels;

    //! Stores the mean for each GC1D of each node (4D)
    std::vector<float> means;
    //! Stores the variance for each GC1D of each node (4D)
    std::vector<float> variances;
    //! Stores n_just_priors for each GC1D of each node (4D)
    std::vector<short> nJustPriors;
    //! Stores ntraining for each GC1D of each node (4D)
    std::vector<int> nTraining;
    //! Stores regularized for each GC1D of each node (4D)
    std::vector<char> regularised;
    

    //! Stores label_priors for each direction of each GC1D of each node (6D)
    std::vector<float> gc1dDirecLabelPriors;
    //! Stores labels for each direction of each GC1D of each node (6D)
    std::vector<unsigned short> gc1dDirecLabels;


    // -------------------------------------------------

    //! Method to extract the max. dimensions from a GCA
    void ExtractDims( const GCA* const src );

    //! Method to allocate memory according to internals
    void Allocate( void );

    // -------------------------------------------------

    //! Timer for exhumation
    mutable SciGPU::Utilities::Chronometer tExhume;
    //! Inhumation timer
    mutable SciGPU::Utilities::Chronometer tInhume;
    
    // -------------------------------------------------

    //! Index computation for 3D matrices
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

    //! Index computation for 4D arrays
    inline unsigned int index4D( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D ) const {
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

      if( (iGC1D<0) || (iGC1D>=(nextOffset-currOffset) ) ) {
	std::cerr << __FUNCTION__
		  << ": iGC1D out of range" << std::endl;
	abort();
      }

      return( currOffset + iGC1D );
    }
    

    //! Index computation for 5D arrays
    inline unsigned int index5D( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D,
				 const int iDirec ) const {
      /*!
	No actual 5D arrays exist in the class, since nlabels in GC1D
	can be dynamically computed.
	However, the Exhume method does have a 5D temporary
	to generate the 6D offset array
      */
      
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	abort();
      }

      const unsigned int idx3D = this->index3D(ix,iy,iz);

      const int currGC1dOffset = this->offsets4D.at( idx3D );
      const int nextGC1dOffset = this->offsets4D.at( idx3D + 1 );

      if( (iGC1D<0) || (iGC1D>=(nextGC1dOffset-currGC1dOffset) ) ) {
	std::cerr << __FUNCTION__
		  << ": iGC1D out of range" << std::endl;
	abort();
      }
      
      const int currOffset = this->offsets5D.at( currGC1dOffset+iGC1D );
      const int nextOffset = this->offsets5D.at( currGC1dOffset+iGC1D+1 );

      if( (iDirec<0) || (iDirec>=(nextOffset-currOffset) ) ) {
	std::cerr << __FUNCTION__
		  << ": iDirec out of range" << std::endl;
	abort();
      }

      return( currOffset + iDirec );
    }
      
  
      
    //! Index computation for 6D arrays
    inline unsigned int index6D( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D,
				 const int iDirec,
				 const int iLabel ) const {

      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Index out of range" << endl;
	abort();
      }

      const unsigned int idx3D = this->index3D(ix,iy,iz);

      const int currgc1DOffset = this->offsets4D.at( idx3D );
      const int nextgc1DOffset = this->offsets4D.at( idx3D + 1 );

      if( (iGC1D<0) || (iGC1D>=(nextgc1DOffset-currgc1DOffset) ) ) {
	std::cerr << __FUNCTION__
		  << ": iGC1D outof range" << std::endl;
	abort();
      }

      const int currDirecOffset = this->offsets5D.at( currgc1DOffset+iGC1D );
      const int nextDirecOffset = this->offsets5D.at( currgc1DOffset+iGC1D+1 );
      
      if( (iDirec<0) || (iDirec>=(nextDirecOffset-currDirecOffset) ) ) {
	std::cerr << __FUNCTION__
		  << ": iDirec out of range" << std::endl;
	abort();
      }

      const int currOffset = this->offsets6D.at( currDirecOffset+iDirec );
      const int nextOffset = this->offsets6D.at( currDirecOffset+iDirec+1 );

      if( (iLabel<0) || (iLabel>=(nextOffset-currOffset)) ) {
	std::cerr << __FUNCTION__
		  << ": iLabel out of range" << std::endl;
	abort();
      }
      
      return( currOffset + iLabel );
    }

  };

}

#endif
