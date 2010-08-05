/**
 * @file  gcalinearnode.hpp
 * @brief Class to hold a volume of GCA nodes in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/05 18:25:30 $
 *    $Revision: 1.4 $
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
			    gc1dDim(0), gc1dNeighbourDim(GIBBS_NEIGHBORHOOD),
			    gc1dLabelDim(0),
			    nGC1Dnode(), maxLabelsNode(), totTrainNode(),
			    nodeLabels(),
			    means(), variances(),
			    nJustPriors(),
			    nTraining(),
			    regularised(),
			    nLabelsGC1D(),
			    gc1dDirecLabels(), gc1dDirecLabelPriors(),
			    bytes(0), tExhume(), tInhume() {};

   
   

    // -------------------------------------------------


    // -------------------------------------------------


    //! Accessor for nGC1Dnode array
    int& gc1dCount( const int ix, const int iy, const int iz ) {
      /*!
	The nGC1Dnode array holds per voxel data - the number
	of GC1D's (and also labels) hanging off each voxel.
       */

      unsigned int index = this->gc1dCountIndex(ix,iy,iz);

      return( this->nGC1Dnode.at(index) );
    }


    //! Const accessor for nGC1Dnode array
    int gc1dCount( const int ix, const int iy, const int iz ) const {
      /*!
	The nGC1Dnode array holds per voxel data - the number
	of GC1D's (and also labels) hanging off each voxel.
       */
      unsigned int index = this->gc1dCountIndex(ix,iy,iz);

      return( this->nGC1Dnode.at(index) );
    }

    // ---

    //! Accessor for maxLabelsNode array
    int& maxLabels( const int ix, const int iy, const int iz ) {
      /*!
	The maxLabels array holds per voxel data 
      */

      unsigned int index = this->maxLabelsIndex(ix,iy,iz);

      return( this->maxLabelsNode.at(index) );
    }

    //! Const accessor for maxLabelsNode array
    int maxLabels( const int ix, const int iy, const int iz ) const {
      /*!
	The maxLabels array holds per voxel data 
      */

      unsigned int index = this->maxLabelsIndex(ix,iy,iz);

      return( this->maxLabelsNode.at(index) );
    }
    

    // ---

    //! Accessor for totTrainNode array
    int& totalTraining( const int ix, const int iy, const int iz ) {
      /*!
	The totTrainNode array holds per voxel data 
      */

      unsigned int index = this->totalTrainingIndex(ix,iy,iz);

      return( this->totTrainNode.at(index) );
    }

    //! Const accessor for totTrainNode array
    int totalTraining( const int ix, const int iy, const int iz ) const {
      /*!
	The totTrainNode array holds per voxel data 
      */

      unsigned int index = this->totalTrainingIndex(ix,iy,iz);

      return( this->totTrainNode.at(index) );
    }


    // ---

    
    //! Accessor for nodeLabels array
    unsigned short& labelsAtNode( const int ix,
				  const int iy,
				  const int iz,
				  const int iGC1D ) {
      /*!
	The nodeLabels array is 4D.
	From each voxel there is a linear array, with a
	number of entries given by the corresponding
	voxel in the nGC1Dnode array.
      */
      
      unsigned int idx;
      idx = this->labelsAtNodeIndex(ix,iy,iz,iGC1D);

      return( this->nodeLabels.at(idx) );
    }
    
    //! Const accessor for nodeLabels array
    unsigned short labelsAtNode( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D ) const {
      /*!
	The nodeLabels array is 4D.
	From each voxel there is a linear array, with a
	number of entries given by the corresponding
	voxel in the nGC1Dnode array.
      */
      unsigned int idx;
      idx = this->labelsAtNodeIndex(ix,iy,iz,iGC1D);

      return( this->nodeLabels.at(idx) );
    }

    // ---


    //! Accessor for the means array
    float& meansAtNodeGC1D( const int ix,
			    const int iy,
			    const int iz,
			    const int iGC1D ) {
      /*!
	The means array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a mean associated with
	it.
	Remember that we are assuming that ninputs==1.
      */
      unsigned int idx;
      
      idx = this->meansAtNodeGC1Dindex(ix,iy,iz,iGC1D);
      
      return( this->means.at(idx) );
    }


    //! Const accessor for the means array
    float meansAtNodeGC1D( const int ix,
			   const int iy,
			   const int iz,
			   const int iGC1D ) const {
      /*!
	The means array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a mean associated with
	it.
	Remember that we are assuming that ninputs==1.
      */
      unsigned int idx;
      
      idx = this->meansAtNodeGC1Dindex(ix,iy,iz,iGC1D);
      
      return( this->means.at(idx) );
    }


    
    // ---



    //! Accessor for the variances array
    float& variancesAtNodeGC1D( const int ix,
				const int iy,
				const int iz,
				const int iGC1D ) {
      /*!
	The variances array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a mean and variance
	associated with it.
	Remember that we are assuming that ninputs==1,
	so the covariance matrix is just a single
	variance.
      */
      
      unsigned int idx;
      
      idx = this->variancesAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->variances.at(idx) );
    }


    //! Accessor for the variances array
    float variancesAtNodeGC1D( const int ix,
			       const int iy,
			       const int iz,
			       const int iGC1D ) const {
      /*!
	The variances array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a mean and variance
	associated with it.
	Remember that we are assuming that ninputs==1,
	so the covariance matrix is just a single
	variance.
      */
      
      unsigned int idx;
      
      idx = this->variancesAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->variances.at(idx) );
    }


    // ---


    //! Accessor for the nJustPriors array
    short& nJustPriorsAtNodeGC1D( const int ix,
				  const int iy,
				  const int iz,
				  const int iGC1D ) {
      /*!
	The nJustPriors array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an n_just_priors value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->nJustPriorsAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->nJustPriors.at(idx) );
    }

    //! Const accessor for the nJustPriors array
    short nJustPriorsAtNodeGC1D( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D ) const {
      /*!
	The nJustPriors array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an n_just_priors value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->nJustPriorsAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->nJustPriors.at(idx) );
    }

    // ---


    //! Accessor for the nTraining array
    int& nTrainingAtNodeGC1D( const int ix,
			      const int iy,
			      const int iz,
			      const int iGC1D ) {
      /*!
	The nTraining array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an ntraining value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->nTrainingAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->nTraining.at(idx) );
    }

     //! Const accessor for the nTraining array
    int nTrainingAtNodeGC1D( const int ix,
			     const int iy,
			     const int iz,
			     const int iGC1D ) const {
      /*!
	The nTraining array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an ntraining value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->nTrainingAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->nTraining.at(idx) );
    }


    // ---

    //! Accessor for the regularised array
    char& regularisedAtNodeGC1D( const int ix,
				 const int iy,
				 const int iz,
				 const int iGC1D ) {
      /*!
	The regularised array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an regularized value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->regularisedAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->regularised.at(idx) );
    }

    //! Const accessor for the regularised array
    char regularisedAtNodeGC1D( const int ix,
				const int iy,
				const int iz,
				const int iGC1D ) const {
      /*!
	The regularised array is 4D.
	For each voxel, there is an array of GC1Ds,
	and each of these has an regularized value
	associated with it.
      */
      
      unsigned int idx;
      
      idx = this->regularisedAtNodeGC1Dindex(ix,iy,iz,iGC1D);

      return( this->regularised.at(idx) );
    }



    // ---


    //! Accessor for the nLabelsGC1D array
    short& nLabelsAtNodeGC1Ddirection( const int ix,
				       const int iy,
				       const int iz,
				       const int iGC1D,
				       const int iDirec ) {
      /*!
	The nLabelsGC1D array is 5-D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	The nLabelsGC1D array holds the length of these
	arrays in each direction
      */

      unsigned int idx;
      
      idx = this->nLabelsAtNodeGC1DdirectionIndex(ix,iy,iz,iGC1D,iDirec);
      return( this->nLabelsGC1D.at(idx) );
    }



    //! Const accessor for the nLabelsGC1D array
    short nLabelsAtNodeGC1Ddirection( const int ix,
				      const int iy,
				      const int iz,
				      const int iGC1D,
				      const int iDirec ) const {
      /*!
	The nLabelsGC1D array is 5-D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	The nLabelsGC1D array holds the length of these
	arrays in each direction
      */

      unsigned int idx;
      
      idx = this->nLabelsAtNodeGC1DdirectionIndex(ix,iy,iz,iGC1D,iDirec);
      return( this->nLabelsGC1D.at(idx) );
    }


    // ---



    //! Accessor for gc1dDirecLabels array
    unsigned short& labelsAtNodeGC1Ddirection( const int ix,
					       const int iy,
					       const int iz,
					       const int iGC1D,
					       const int iDirec,
					       const int iLabel ) {
      /*!
	The gc1dDirecLabels array is 6D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	This is the accessor for the labels themselves
      */

      unsigned int idx;
      idx = this->labelsAtNodeGC1DdirectionIndex(ix,iy,iz,iGC1D,iDirec,iLabel);
      return( this->gc1dDirecLabels.at(idx) );
    }


    //! Const accessor for gc1dDirecLabels array
    unsigned short labelsAtNodeGC1Ddirection( const int ix,
					      const int iy,
					      const int iz,
					      const int iGC1D,
					      const int iDirec,
					      const int iLabel ) const {
      /*!
	The gc1dDirecLabels array is 6D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	This is the accessor for the labels themselves
      */
      unsigned int idx;
      idx = this->labelsAtNodeGC1DdirectionIndex(ix,iy,iz,iGC1D,iDirec,iLabel);
      return( this->gc1dDirecLabels.at(idx) );
    }


    // ---
    

    //! Accessor for gc1dDirecLabelPriors array
    float& labelPriorsAtNodeGC1Ddirection( const int ix,
					   const int iy,
					   const int iz,
					   const int iGC1D,
					   const int iDirec,
					   const int iLabel ) {
      /*!
	The gc1dDirecLabels array is 6D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	This is the accessor for the labels themselves
      */

      unsigned int idx;
      idx = this->labelPriorsAtNodeGC1DdirectionIndex( ix, iy, iz,
						       iGC1D,
						       iDirec,
						       iLabel );
      
      return( this->gc1dDirecLabelPriors.at(idx) );
    }


    //! Const accessor for gc1dDirecLabelPriors array
    float labelPriorsAtNodeGC1Ddirection( const int ix,
					  const int iy,
					  const int iz,
					  const int iGC1D,
					  const int iDirec,
					  const int iLabel ) const {
      /*!
	The gc1dDirecLabels array is 6D.
	For each voxel, there is an array of GC1Ds,
	and each of these has a label and label prior
	array hanging off it in each of six directions.
	This is the accessor for the labels themselves
      */
      unsigned int idx;
      idx = this->labelPriorsAtNodeGC1DdirectionIndex( ix, iy, iz,
						       iGC1D,
						       iDirec,
						       iLabel );

      return( this->gc1dDirecLabelPriors.at(idx) );
    }


    // -------------------------------------------------

    //! Method to extract node data from a GCA
    void Exhume( const GCA* const src );

    //! Method to place node data into a GCA
    void Inhume( GCA* dst ) const;

    //! Method to remove nodes from a GCA prior to inhumation of new data
    void ScorchNodes( GCA* targ ) const;

    // -------------------------------------------------

    //! Method to print out statistics & timers
    void PrintStats( std::ostream& os = std::cout ) const ;

    // -------------------------------------------------
  private:

    // Dimensions
    int xDim;
    int yDim;
    int zDim;
    int gc1dDim;
    const int gc1dNeighbourDim;
    int gc1dLabelDim;

    //! Stores number of GC1D hanging off each node (nlabels in GCA_NODE)
    std::vector<int> nGC1Dnode;

    //! Stores max_labels in each node
    std::vector<int> maxLabelsNode;

    //! Stores total_training in each node
    std::vector<int> totTrainNode;

    //! Stores the labels for the node (4D, controlled via gca1dDim and gc1dCount)
    std::vector<unsigned short> nodeLabels;

    //! Stores the means for each GC1D of the node (4-D)
    std::vector<float> means;
    //! Stores the variances for each GC1D of the node
    std::vector<float> variances;

    //! Stores n_just_priors for each GC1D of the node (4D)
    std::vector<short> nJustPriors;

    //! Stores ntraining for each GC1D of the node (4D)
    std::vector<int> nTraining;

    //! Stores the 'regularizied' member for each GC1D of the node (4D)
    std::vector<char> regularised;

    //! Stores the number of labels for each direction for each GC1D (5D)
    std::vector<short> nLabelsGC1D;

    //! Stores the labels for each direction of each GC1D (6D)
    std::vector<unsigned short> gc1dDirecLabels;

    //! Stores the label priors for each direction of each GC1D (6D)
    std::vector<float> gc1dDirecLabelPriors;

    // -------------------------------------------------

    //! Method to extract the max. dimensions from a GCA
    void ExtractDims( const GCA* const src );

    //! Method to allocate memory according to internals
    void Allocate( void );

    // -------------------------------------------------
    
    //! Total size of the allocated arrays
    size_t bytes;

    //! Timer for exhumation
    mutable SciGPU::Utilities::Chronometer tExhume;
    //! Inhumation timer
    mutable SciGPU::Utilities::Chronometer tInhume;
    
    // -------------------------------------------------

    //! Index computation for gc1dCount
    unsigned int gc1dCountIndex( const int ix,
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

    // ---

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

    // ---

    //! Index computation for total training
    unsigned int totalTrainingIndex( const int ix,
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
    

    // ---

    //! Index computation for labelsAtNode
    unsigned int labelsAtNodeIndex( const int ix,
				    const int iy,
				    const int iz,
				    const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      unsigned int idx;

      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );
      
      return(idx);
    }

    
    // ---

    //! Index computation for meansAtNodeGC1D
    unsigned int meansAtNodeGC1Dindex( const int ix,
				       const int iy,
				       const int iz,
				       const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );

      return( idx );
    }

    // ---

    //! Index computation for variancesAtNodeGC1D
    unsigned int variancesAtNodeGC1Dindex( const int ix,
					   const int iy,
					   const int iz,
					   const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );
      
      return( idx );
    }

    
    // ---

    //! Index computation for nJustPriorsAtNodeGC1D
    unsigned int nJustPriorsAtNodeGC1Dindex( const int ix,
					     const int iy,
					     const int iz,
					     const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );
      
      return( idx );
    }
      


    // ---

    //! Index computation for nTrainingAtNodeGC1D
    unsigned int nTrainingAtNodeGC1Dindex( const int ix,
					   const int iy,
					   const int iz,
					   const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );
      
      return( idx );
    }

    // ---

    //! Index computation for regularisedAtNodeGC1D
    unsigned int regularisedAtNodeGC1Dindex( const int ix,
					     const int iy,
					     const int iz,
					     const int iGC1D ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim * iGC1D ) ) ) ) );
      
      return( idx );
    }


    
    // ---

    //! Index computation for nLabelsAtNodeGC1Ddirection
    unsigned int nLabelsAtNodeGC1DdirectionIndex( const int ix,
						  const int iy,
						  const int iz,
						  const int iGC1D,
						  const int iDirec ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iDirec<0) || (iDirec>=this->gc1dNeighbourDim) ) {
	cerr << __FUNCTION__
	     << ": iDirec index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim *
				     ( iGC1D + ( this->gc1dDim *
						 iDirec ) )
				     ) )
			    ) )
		   );
      
      return( idx );
    }

    // ---

    //! Index computation for labelsAtNodeGC1Ddirection
    unsigned int labelsAtNodeGC1DdirectionIndex( const int ix,
						 const int iy,
						 const int iz,
						 const int iGC1D,
						 const int iDirec,
						 const int iLabel ) const {
      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      if( (iDirec<0) || (iDirec>=this->gc1dNeighbourDim) ) {
	cerr << __FUNCTION__
	     << ": iDirec index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      if( (iLabel<0) || (iLabel>=this->gc1dLabelDim) ||
	  (iLabel>=this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec) ) ) {
	cerr << __FUNCTION__
	     << ": iLabel index out of range" << endl;
	exit( EXIT_FAILURE );
      }

      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim *
				     ( iGC1D + ( this->gc1dDim *
						 ( iDirec + ( this->gc1dNeighbourDim *
							      iLabel ) )
						 ) )
				     ) )
			    ) )
		   );


      return( idx );
    }

    // ---

    //! Index computation for labelPriorsAtNodeGC1Ddirection
    unsigned int labelPriorsAtNodeGC1DdirectionIndex( const int ix,
						      const int iy,
						      const int iz,
						      const int iGC1D,
						      const int iDirec,
						      const int iLabel ) const {

      if( (ix<0) || (ix>=this->xDim) ||
	  (iy<0) || (iy>=this->yDim) ||
	  (iz<0) || (iz>=this->zDim ) ) {
	cerr << __FUNCTION__
	     << ": Voxel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iGC1D<0) || (iGC1D>=this->gc1dDim) ||
	  (iGC1D >= this->gc1dCount(ix,iy,iz)) ) {
	cerr << __FUNCTION__
	     << ": GC1D index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iDirec<0) || (iDirec>=this->gc1dNeighbourDim) ) {
	cerr << __FUNCTION__
	     << ": iDirec index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      if( (iLabel<0) || (iLabel>=this->gc1dLabelDim) ||
	  (iLabel>=this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec) ) ) {
	cerr << __FUNCTION__
	     << ": iLabel index out of range" << endl;
	exit( EXIT_FAILURE );
      }
      
      unsigned int idx;
      
      idx = ix + ( this->xDim *
		   ( iy + ( this->yDim *
			    ( iz + ( this->zDim *
				     ( iGC1D + ( this->gc1dDim *
						 ( iDirec + ( this->gc1dNeighbourDim *
							      iLabel ) )
						 ) )
				     ) )
			    ) )
		   );

      return( idx );
    }

  };

}

#endif
