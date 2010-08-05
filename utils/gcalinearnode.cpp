/**
 * @file  gcalinearnode.cpp
 * @brief Class to hold a volume of GCA nodes in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/05 16:05:24 $
 *    $Revision: 1.3 $
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


#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
using namespace std;


#include "minmax.hpp"

#include "gcalinearnode.hpp"

namespace Freesurfer {

  // ====================================================

  void GCAlinearNode::Exhume( const GCA* const src ) {
    /*!
      This method is responsible for extracting data from
      the source GCA, and packing into the linear arrays
    */

    this->tExhume.Start();

    this->ExtractDims( src );
    this->Allocate();

    for( int ix=0; ix<this->xDim; ix++ ) {
      for( int iy=0; iy<this->yDim; iy++ ) {
	for( int iz=0; iz<this->zDim; iz++ ) {
	  const GCA_NODE* const gcan = &(src->nodes[ix][iy][iz]);

	  this->gc1dCount(ix,iy,iz) = gcan->nlabels;

	  
	  this->maxLabels(ix,iy,iz) = gcan->max_labels;
	  this->totalTraining(ix,iy,iz) = gcan->total_training;

	  for( int iGC1D=0; iGC1D<this->gc1dCount(ix,iy,iz); iGC1D++ ) {
	    
	    // Deal with the label held inside each node
	    this->labelsAtNode(ix,iy,iz,iGC1D) = gcan->labels[iGC1D];
	    

	    // Grab hold of the current GC1D
	    const GC1D* const gc1d = &(gcan->gcs[iGC1D]);
	    
	    // Deal with the means and variances
	    // Recall that ninputs==1
	    this->meansAtNodeGC1D(ix,iy,iz,iGC1D) = gc1d->means[0];
	    this->variancesAtNodeGC1D(ix,iy,iz,iGC1D) = gc1d->covars[0];

	    this->nJustPriorsAtNodeGC1D(ix,iy,iz,iGC1D) = gc1d->n_just_priors;
	    this->nTrainingAtNodeGC1D(ix,iy,iz,iGC1D) = gc1d->ntraining;
	    this->regularisedAtNodeGC1D(ix,iy,iz,iGC1D) = gc1d->regularized;

	    for( int iDirec=0;
		 iDirec<this->gc1dNeighbourDim;
		 iDirec++ ) {
	      // Unpack the number of labels for each direction
	      this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec) =
		gc1d->nlabels[iDirec];

	      for( int iLabel=0;
		   iLabel<this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec);
		   iLabel++ ) {
		this->labelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec,iLabel) =
		  gc1d->labels[iDirec][iLabel];
		this->labelPriorsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec,iLabel) =
		  gc1d->label_priors[iDirec][iLabel];
	      }
		
	      
	    }
	    
	  }

	}
      }
    }

    this->tExhume.Stop();

  }

  // ====================================================

  void GCAlinearNode::Inhume( GCA* dst ) const {
    /*!
      Stores data back from the linear arrays into a GCA.
      The GCA must already have been allocated to be
      of appropriate dimensions
    */

    this->tInhume.Start();

    if( dst->ninputs != 1 ) {
      cerr << __FUNCTION__
	   << ": Must have ninputs==1!" << endl;
      exit( EXIT_FAILURE );
    }

    if( (this->xDim!=dst->node_width) ||
	(this->yDim!=dst->node_height) ||
	(this->zDim!=dst->node_depth) ) {
      cerr << __FUNCTION__
	   << ": Volume dimension mismatch!"
	   << endl;
      exit( EXIT_FAILURE );
    }

    

    for( int ix=0; ix<this->xDim; ix++ ) {
      for( int iy=0; iy<this->yDim; iy++ ) {
	for( int iz=0; iz<this->zDim; iz++ ) {
	  GCA_NODE* const gcan = &(dst->nodes[ix][iy][iz]);

	  // Check the number of labels at the node
	  if( gcan->nlabels != this->gc1dCount(ix,iy,iz) ) {
	    cerr << __FUNCTION__
		 << ": nlabels mismatch at "
		 << "( " << ix << ", " << iy << ", " << iz << " )"
		 << endl;
	    exit( EXIT_FAILURE );
	  }

	  gcan->max_labels = this->maxLabels(ix,iy,iz);
	  gcan->total_training = this->totalTraining(ix,iy,iz);

	  // Loop over the GC1Ds
	  for( int iGC1D=0; iGC1D<this->gc1dCount(ix,iy,iz); iGC1D++ ) {
	    // Do the labels on the side
	    gcan->labels[iGC1D] = this->labelsAtNode(ix,iy,iz,iGC1D);

	    GC1D* const gc1d = &(gcan->gcs[iGC1D]);

	    // Do the mean and variance (recall ninputs==1)
	    gc1d->means[0] = this->meansAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->covars[0] = this->variancesAtNodeGC1D(ix,iy,iz,iGC1D);

	    gc1d->n_just_priors = this->nJustPriorsAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->ntraining = this->nTrainingAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->regularized = this->regularisedAtNodeGC1D(ix,iy,iz,iGC1D);

	    for( int iDirec=0;
		 iDirec<this->gc1dNeighbourDim;
		 iDirec++ ) {
	      // Check the number of labels in each direction
	      if( this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec) !=
		  gc1d->nlabels[iDirec] ) {
		cerr << __FUNCTION__
		     << ": Mismatch in labels at node "
		     << "( " << ix << ", " << iy << ", " << iz << " )"
		     << " iGC1D = " << iGC1D
		     << " iDirec = " << iDirec
		     << endl;
		exit( EXIT_FAILURE );
	      }

	      for( int iLabel=0;
		   iLabel<this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec);
		   iLabel++ ) {
		gc1d->labels[iDirec][iLabel] =
		  this->labelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec,iLabel);
		gc1d->label_priors[iDirec][iLabel] =
		  this->labelPriorsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec,iLabel);
	      }

	    }

	  }

	}
      }
    }

    this->tInhume.Stop();

  }

    
  // ====================================================

  void GCAlinearNode::Allocate( void ) {
    /*!
      Allocates the main arrays according to
      the values stored in the dimension member
      variables.
      All vectors are cleared first, in order to
      fill with defaults
    */

    this->bytes = 0;

    const size_t nVoxels = this->xDim * this->yDim * this->zDim;

    // Indicate the number of GC1D and labels for this node
    this->nGC1Dnode.clear();
    this->nGC1Dnode.resize( nVoxels, 0 );
    this->bytes += this->nGC1Dnode.size() * sizeof(int);

    // Allocate the max_labels held in the GCA_NODE structure
    this->maxLabelsNode.clear();
    this->maxLabelsNode.resize( nVoxels,
				numeric_limits<int>::max() );
    this->bytes += this->maxLabelsNode.size() * sizeof(int);
    
    // Allocate the total_training held in the GCA_NODE structure
    this->totTrainNode.clear();
    this->totTrainNode.resize( nVoxels,
			       numeric_limits<int>::max() );
    this->bytes += this->totTrainNode.size() * sizeof(int);

    // Allocate the labels held in the GCA_NODE structure
    this->nodeLabels.clear();
    this->nodeLabels.resize( nVoxels * this->gc1dDim,
			     numeric_limits<unsigned short>::max() );
    this->bytes += this->nodeLabels.size() * sizeof(unsigned short);


    // Allocate the means held in the GC1D structure (ninputs=1)
    this->means.clear();
    this->means.resize( nVoxels * this->gc1dDim,
			numeric_limits<float>::quiet_NaN() );
    this->bytes += this->means.size() * sizeof(float);

    // Allocate the variances held in the GC1D structure (ninputs=1)
    this->variances.clear();
    this->variances.resize( nVoxels * this->gc1dDim,
			    numeric_limits<float>::quiet_NaN() );
    this->bytes += this->variances.size() * sizeof(float);

    //! Allocate nJustPriors for the GC1D structures
    this->nJustPriors.clear();
    this->nJustPriors.resize( nVoxels * this->gc1dDim,
			      numeric_limits<short>::max() );
    this->bytes += this->nJustPriors.size() * sizeof(short);

    //! Allocate nTraining for the GC1D structures
    this->nTraining.clear();
    this->nTraining.resize( nVoxels * this->gc1dDim,
			    numeric_limits<int>::max() );
    this->bytes += this->nTraining.size() * sizeof(int);

    //! Allocate regularised for the GC1D structures
    this->regularised.clear();
    this->regularised.resize( nVoxels * this->gc1dDim,
			      numeric_limits<char>::max() );
    this->bytes += this->regularised.size() * sizeof(char);

    // Allocate the number of labels for each direction for each GC1D
    this->nLabelsGC1D.clear();
    this->nLabelsGC1D.resize( nVoxels * this->gc1dDim * this->gc1dNeighbourDim,
			      numeric_limits<short>::max() );
    this->bytes += this->nLabelsGC1D.size() * sizeof(short);

    // Allocate the labels for each direction of each GC1D
    this->gc1dDirecLabels.clear();
    this->gc1dDirecLabels.resize( nVoxels * this->gc1dDim * this->gc1dNeighbourDim * this->gc1dLabelDim,
				  numeric_limits<unsigned short>::max() );
    this->bytes += this->gc1dDirecLabels.size() * sizeof(unsigned short);

    // Allocate the label priors for each direction of each GC1D
    this->gc1dDirecLabelPriors.clear();
    this->gc1dDirecLabelPriors.resize(  nVoxels * this->gc1dDim * this->gc1dNeighbourDim * this->gc1dLabelDim,
					numeric_limits<float>::quiet_NaN() );
    this->bytes += this->gc1dDirecLabelPriors.size() * sizeof(float);
  }




  // ====================================================
  
  void GCAlinearNode::ExtractDims( const GCA* const src ) {
    /*!
      Fills in the dimensions required from the given GCA.
      Does this by looping over all voxels, and finding the
      maximum number of elements allocated for each.
    */
    if( src->ninputs != 1 ) {
      cerr << __FUNCTION__
	   << ": Must have ninputs==1!" << endl;
      exit( EXIT_FAILURE );
    }

    this->xDim = src->node_width;
    this->yDim = src->node_height;
    this->zDim = src->node_depth;

    MinMax<int> nLabelsNode;
    MinMax<short> nLabelsGC1D;

    for( int ix=0; ix<this->xDim; ix++ ) {
      for( int iy=0; iy<this->yDim; iy++ ) {
	for( int iz=0; iz<this->zDim; iz++ ) {
	  const GCA_NODE* const gcan = &(src->nodes[ix][iy][iz]);

	  nLabelsNode.Accumulate( gcan->nlabels );
	  for( int iGC1D=0; iGC1D<gcan->nlabels; iGC1D++ ) {
	    const GC1D* const gc1d = &(gcan->gcs[iGC1D]);

	    for( unsigned int iNeighbour=0;
		 iNeighbour<GIBBS_NEIGHBORHOOD;
		 iNeighbour++ ) {
	      nLabelsGC1D.Accumulate( gc1d->nlabels[iNeighbour] );
	    }
	  }

	}
      }
    }

    this->gc1dDim = nLabelsNode.GetMax();
    this->gc1dLabelDim = nLabelsGC1D.GetMax();

    double mean = nLabelsGC1D.GetDoubleTotal() / nLabelsGC1D.GetN();
    cout << "nLabels GC1D: "
	 << nLabelsGC1D
	 << " (avg. "
	 << setprecision(3) << mean
	 << " )" << endl;
      
  }


  // ====================================================

  void GCAlinearNode::PrintStats( ostream& os ) const {

    os << "Stats for GCAlinearNode" << endl;
    os << "  Bytes allocated = " << this->bytes << endl;
    os << "  Exhumation time = " << this->tExhume << endl;
    os << "  Inhumation time = " << this->tInhume << endl;

  }
}
