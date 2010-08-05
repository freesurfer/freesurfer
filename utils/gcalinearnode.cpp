/**
 * @file  gcalinearnode.cpp
 * @brief Class to hold a volume of GCA nodes in linear memory
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/08/05 18:25:33 $
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
      This removes the existing node structure, and reallocates
      it iself.
      One might prefer this not to happen, and to use the
      routines from gca.c.
      However, that would involve writing some of them first.
    */

    this->tInhume.Start();

    if( (this->xDim!=dst->node_width) ||
	(this->yDim!=dst->node_height) ||
	(this->zDim!=dst->node_depth) ) {
      cerr << __FUNCTION__
	   << ": Volume dimension mismatch!"
	   << endl;
      exit( EXIT_FAILURE );
    }

    
    // Dispose of the old node data
    this->ScorchNodes( dst );


    // Set dimensions
    dst->ninputs = 1;
    dst->node_width = this->xDim;
    dst->node_height = this->yDim;
    dst->node_depth = this->zDim;

    // Dispose of old node data
    dst->nodes = (GCA_NODE***)calloc( this->xDim, sizeof(GCA_NODE**) );
    if( !(dst->nodes) ) {
      cerr << __FUNCTION__
	   << ": dst->nodes allocation failed" << endl;
      exit( EXIT_FAILURE );
    }

    for( int ix=0; ix<this->xDim; ix++ ) {
      // Allocate pointer block
      dst->nodes[ix] = (GCA_NODE**)calloc( this->yDim, sizeof(GCA_NODE*) );
      if( !(dst->nodes[ix]) ) {
	cerr << __FUNCTION__
	     << ": dst->nodes[" << ix << "] allocation failed" << endl;
	exit( EXIT_FAILURE );
      }

      for( int iy=0; iy<this->yDim; iy++ ) {
	// Allocate pointer block
	dst->nodes[ix][iy] = (GCA_NODE*)calloc( this->zDim, sizeof(GCA_NODE) );
	if( !(dst->nodes[ix][iy]) ) {
	  cerr << __FUNCTION__
	       << ": dst->nodes"
	       << "[" << ix << "]"
	       << "[" << iy << "]"
	       << " allocation failed" << endl;
	  exit( EXIT_FAILURE );
	}
	
	for( int iz=0; iz<this->zDim; iz++ ) {
	  // Allocate pointer block
	  GCA_NODE* const gcan = &(dst->nodes[ix][iy][iz]);

	  gcan->nlabels = this->gc1dCount(ix,iy,iz);
	  gcan->max_labels = this->maxLabels(ix,iy,iz);
	  gcan->total_training = this->totalTraining(ix,iy,iz);

	  // Allocate labels array
	  gcan->labels = (unsigned short*)calloc( this->gc1dCount(ix,iy,iz),
						  sizeof(unsigned short) );
	  if( !(gcan->labels) ) {
	    cerr << __FUNCTION__
		 << ": dst->nodes"
		 << "[" << ix << "]"
		 << "[" << iy << "]"
		 << "[" << iz << "].labels"
		 << " allocation failed" << endl;
	    exit( EXIT_FAILURE );
	  }

	  // Allocate GC1D array
	  gcan->gcs = (GC1D*)calloc( this->gc1dCount(ix,iy,iz),
				     sizeof(GC1D) );
	  if( !(gcan->gcs) ) {
	    cerr << __FUNCTION__
		 << ": dst->nodes"
		 << "[" << ix << "]"
		 << "[" << iy << "]"
		 << "[" << iz << "].gcs"
		 << " allocation failed" << endl;
	    exit( EXIT_FAILURE );
	  }
				     

	  // Loop over the GC1Ds
	  for( int iGC1D=0; iGC1D<this->gc1dCount(ix,iy,iz); iGC1D++ ) {
	    // Do the labels on the side
	    gcan->labels[iGC1D] = this->labelsAtNode(ix,iy,iz,iGC1D);

	    GC1D* const gc1d = &(gcan->gcs[iGC1D]);

	    gc1d->means = (float*)calloc( dst->ninputs, // Always 1
					  sizeof(float) );
	    if( !(gc1d->means) ) {
	      cerr << __FUNCTION__
		   << ": Allocation failure of means"
		   << endl;
	      exit( EXIT_FAILURE );
	    }

	    gc1d->covars = (float*)calloc( dst->ninputs, // Always 1
					   sizeof(float) );
	    if( !(gc1d->covars) ) {
	      cerr << __FUNCTION__
		   << ": Allocation failure of covars"
		   << endl;
	      exit( EXIT_FAILURE );
	    }

	    // Do the mean and variance (recall ninputs==1)
	    gc1d->means[0] = this->meansAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->covars[0] = this->variancesAtNodeGC1D(ix,iy,iz,iGC1D);

	    gc1d->n_just_priors = this->nJustPriorsAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->ntraining = this->nTrainingAtNodeGC1D(ix,iy,iz,iGC1D);
	    gc1d->regularized = this->regularisedAtNodeGC1D(ix,iy,iz,iGC1D);

	    // Allocate the nlabels array
	    gc1d->nlabels = (short*)calloc( this->gc1dNeighbourDim, // Always 6/GIBBS_NEIGHBORHOOD
					    sizeof(short) );
	    if( !(gc1d->nlabels) ) {
	      cerr << __FUNCTION__
		   << ": Allocation failure of nlabels"
		   << endl;
	      exit( EXIT_FAILURE );
	    }

	    // Allocate pointers for label_priors
	    gc1d->label_priors = (float**)calloc( this->gc1dNeighbourDim,
						  sizeof(float*) );
	    if( !(gc1d->label_priors) ) {
	      cerr << __FUNCTION__
		   << ": Allocation failure of label_priors"
		   << endl;
	      exit( EXIT_FAILURE );
	    }

	    // Allocate pointers for labels
	    gc1d->labels = (unsigned short**)calloc( this->gc1dNeighbourDim,
						     sizeof(unsigned short*) );
	    if( !(gc1d->labels) ) {
	      cerr << __FUNCTION__
		   << ": Allocation failure of labels"
		   << endl;
	      exit( EXIT_FAILURE );
	    }

	    for( int iDirec=0;
		 iDirec<this->gc1dNeighbourDim;
		 iDirec++ ) {

	      // Set the number
	      gc1d->nlabels[iDirec] = 
		this->nLabelsAtNodeGC1Ddirection(ix,iy,iz,iGC1D,iDirec);

	      // Allocate the memory
	      gc1d->label_priors[iDirec] = (float*)calloc( gc1d->nlabels[iDirec],
							   sizeof(float) );
	      if( !(gc1d->label_priors[iDirec]) ) {
		cerr << __FUNCTION__
		     << ": Allocation failure of label_priors"
		     << endl;
		exit( EXIT_FAILURE );
	      }

	      gc1d->labels[iDirec] = (unsigned short*)calloc( gc1d->nlabels[iDirec],
							      sizeof(unsigned short) );
	      if( !(gc1d->labels[iDirec]) ) {
		cerr << __FUNCTION__
		     << ": Allocation failure of labels"
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

  void GCAlinearNode::ScorchNodes( GCA* targ ) const {
    /*!
      Deletes all of the node related things from
      a GCA, prior to inhumation of new data
    */

    for( int ix=0; ix<targ->node_width; ix++ ) {
      for( int iy=0; iy<targ->node_height; iy++ ) {
	for( int iz=0; iz<targ->node_depth; iz++ ) {
	  GCANfree( &(targ->nodes[ix][iy][iz]), targ->ninputs );
	}
	free( targ->nodes[ix][iy] );
      }
      free( targ->nodes[ix] );
    }
    free( targ->nodes );
  }


  // ====================================================

  void GCAlinearNode::PrintStats( ostream& os ) const {

    os << "Stats for GCAlinearNode" << endl;
    os << "  Bytes allocated = " << this->bytes << endl;
    os << "  Exhumation time = " << this->tExhume << endl;
    os << "  Inhumation time = " << this->tInhume << endl;

  }
}
