/**
 * @file  gcamorphgpu.cu
 * @brief Holds GCA morph data on the GPU
 *
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/02/24 15:23:35 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2008,
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

#include "chronometer.hpp"


#include "gcamorphgpu.hpp"

// ==============================================================

namespace GPU {
  namespace Classes {

    // --------------------------------------------

    void GCAmorphGPU::CheckIntegrity( void ) const {
      /*!
	Checks that all the allocated member arrays have
	the same dimensions.
	Aborts the program if the check fails
      */

      const dim3 myDims = this->d_r.GetDims();

      bool good = ( myDims == this->d_invalid.GetDims() );
      good = ( good && ( myDims == this->d_area.GetDims() ) );
      good = ( good && ( myDims == this->d_area1.GetDims() ) );
      good = ( good && ( myDims == this->d_area2.GetDims() ) );

      if( !good ) {
	std::cerr << __FUNCTION__
		  << ": Dimension mismatch"
		  << std::endl;
	exit( EXIT_FAILURE );
      }
    }

    // --------------------------------------------

    void GCAmorphGPU::AllocateAll( const dim3& dims ) {
      /*!
	Allocates GPU memory to hold a volume
	of the given size.
	If possible, it keeps the current allocation.
      */

      // Start by seeing if the current allocation is consistent
      this->CheckIntegrity();

      // See if we can re-use existing allocation
      if( dims == this->d_r.GetDims() ) {
	return;
      }

      // Release existing memory
      this->ReleaseAll();

      // Allocate anew
      this->d_r.Allocate( dims );
      this->d_invalid.Allocate( dims );
      this->d_area.Allocate( dims );
      this->d_area1.Allocate( dims );
      this->d_area2.Allocate( dims );
    }


    void GCAmorphGPU::ReleaseAll( void ) {
      /*!
	Releases each of the members.
	Recall that the VolumeGPU::Release method
	will also release any CUDA arrays.
      */
      this->d_r.Release();
      this->d_invalid.Release();
      this->d_area.Release();
      this->d_area1.Release();
      this->d_area2.Release();
    }

    // --------------------------------------------

    void GCAmorphGPU::SendAll( const GCAM* src ) {
      /*!
	Sends all supported data in the given GCAM
	to the GPU.
	This involves a lot of packing data, and hence
	is going to be painfully slow
      */
      
      SciGPU::Utilities::Chronometer t_tot;
      SciGPU::Utilities::Chronometer t_mem, t_pack, t_send;

      t_tot.Start();

      // Extract the dimensions
      const dim3 dims = make_uint3( src->width,
				    src->height,
				    src->depth );

      t_mem.Start();
      // Allocate device memory
      this->AllocateAll( dims );

      // Allocate some page-locked host buffers
      float3* h_r = this->d_r.AllocateHostBuffer();
      unsigned char* h_invalid = this->d_invalid.AllocateHostBuffer();
      float* h_area = this->d_area.AllocateHostBuffer();
      float* h_area1 = this->d_area1.AllocateHostBuffer();
      float* h_area2 = this->d_area2.AllocateHostBuffer();
      t_mem.Stop();

      t_pack.Start();
      for( unsigned int i=0; i<dims.x; i++ ) {
	for( unsigned int j=0; j<dims.y; j++ ) {
	  for( unsigned int k=0; k<dims.z; k++ ) {

	    // Get the 1d index (same for all arrays)
	    const unsigned int i1d = this->d_r.Index1D( i, j, k );
	    // Get the current node
	    const GCA_MORPH_NODE gcamn = src->nodes[i][j][k];
	    
	    // Pack the data
	    h_r[i1d] = make_float3( gcamn.x,
				    gcamn.y,
				    gcamn.z );

	    h_invalid[i1d] = gcamn.invalid;
	    h_area[i1d] = gcamn.area;
	    h_area1[i1d] = gcamn.area1;
	    h_area2[i1d] = gcamn.area2;
	  }
	}
      }
      t_pack.Stop();

      t_send.Start();
      // Send the data
      this->d_r.SendBuffer( h_r );
      this->d_invalid.SendBuffer( h_invalid );
      this->d_area.SendBuffer( h_area );
      this->d_area1.SendBuffer( h_area1 );
      this->d_area2.SendBuffer( h_area2 );

      // Wait for the copies to complete
      CUDA_SAFE_CALL( cudaThreadSynchronize() );
      t_send.Stop();

      // Release page-locked host memory
      t_mem.Start();
      CUDA_SAFE_CALL( cudaFreeHost( h_r ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_invalid ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area1 ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area2 ) );
      t_mem.Stop();

      t_tot.Stop();

      std::cout << __FUNCTION__ << std::endl;
      std::cout << "t_mem " << t_mem << std::endl;
      std::cout << "t_pack " << t_pack << std::endl;
      std::cout << "t_send " << t_send << std::endl;
      std::cout << "Total : " << t_tot << std:: endl;
      std::cout << "-------------" << std::endl;
    }

    // --------------------------------------------

    void GCAmorphGPU::RecvAll( GCAM* dst ) const {
      /*!
	Retrieves all supported data in the given GCAM
	from the GPU.
	This involves a lot of packing data, and hence
	is going to be painfully slow
      */
      
      SciGPU::Utilities::Chronometer t_tot;
      SciGPU::Utilities::Chronometer t_mem, t_pack, t_recv;

      t_tot.Start();

      // Extract the dimensions
      const dim3 dims = this->d_r.GetDims();

      t_mem.Start();
     
      // Allocate some page-locked host buffers
      float3* h_r = this->d_r.AllocateHostBuffer();
      unsigned char* h_invalid = this->d_invalid.AllocateHostBuffer();
      float* h_area = this->d_area.AllocateHostBuffer();
      float* h_area1 = this->d_area1.AllocateHostBuffer();
      float* h_area2 = this->d_area2.AllocateHostBuffer();
      t_mem.Stop();

      // Fetch the data
      t_recv.Start();
      this->d_r.RecvBuffer( h_r );
      this->d_invalid.RecvBuffer( h_invalid );
      this->d_area.RecvBuffer( h_area );
      this->d_area1.RecvBuffer( h_area1 );
      this->d_area2.RecvBuffer( h_area2 );
      CUDA_SAFE_CALL( cudaThreadSynchronize() );
      t_recv.Stop();

      t_pack.Start();
      for( unsigned int i=0; i<dims.x; i++ ) {
	for( unsigned int j=0; j<dims.y; j++ ) {
	  for( unsigned int k=0; k<dims.z; k++ ) {

	    // Get the 1d index (same for all arrays)
	    const unsigned int i1d = this->d_r.Index1D( i, j, k );
	    // Get the current node
	    GCA_MORPH_NODE gcamn = dst->nodes[i][j][k];

	    gcamn.x = h_r[i1d].x;
	    gcamn.y = h_r[i1d].y;
	    gcamn.z = h_r[i1d].z;

	    gcamn.invalid = h_invalid[i1d];
	    gcamn.area = h_area[i1d];
	    gcamn.area1 = h_area1[i1d];
	    gcamn.area2 = h_area2[i1d];
	  }
	}
      }
      t_pack.Stop();


      // Release page-locked host memory
      t_mem.Start();
      CUDA_SAFE_CALL( cudaFreeHost( h_r ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_invalid ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area1 ) );
      CUDA_SAFE_CALL( cudaFreeHost( h_area2 ) );
      t_mem.Stop();

      t_tot.Stop();

      std::cout << __FUNCTION__ << std::endl;
      std::cout << "t_mem " << t_mem << std::endl;
      std::cout << "t_pack " << t_pack << std::endl;
      std::cout << "t_recv " << t_recv << std::endl;
      std::cout << "Total : " << t_tot << std:: endl;
      std::cout << "-------------" << std::endl;
    }
  }
}

#include "testgpu.h"

void TestGCAMorphGPU( GCAM* src ) {

  GPU::Classes::GCAmorphGPU myMorph;

  myMorph.SendAll( src );

  myMorph.RecvAll( src );

}
