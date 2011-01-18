/**
 * @file  gcanode.hpp
 * @brief Class to hold a volume of GCA nodes in linear memory on the GPU
 *
 */
/*
 * Original Authors: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2011/01/18 17:53:43 $
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

#ifndef GCA_NODE_GPU_HPP
#define GCA_NODE_GPU_HPP

#include <iostream>

#include "gca.h"
#include "gcalinearnode.hpp"


namespace GPU {

  namespace Classes {




    class GCAnodeGPU {
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
      ~GCAnodeGPU( void ) {
	this->Release();
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
					   d_gc1dDirecLabels(NULL) {
	std::cerr << __FUNCTION__
		  << ": Please don't copy"
		  << std::endl;
	abort();
      }

      //! Inhibit assignment
      GCAnodeGPU& operator=( const GCAnodeGPU& src ) {
	std::cerr << __FUNCTION__
		  << ": Please don't copy"
		  << std::endl;
	abort();
      }
      
    };

  }
}




#endif
