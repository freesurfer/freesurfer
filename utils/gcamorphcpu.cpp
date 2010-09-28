/**
 * @file  gcamorphcpu.cpp
 * @brief Holds GCA morph data on the CPU
 *
 * 
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/09/28 19:40:25 $
 *    $Revision: 1.2 $
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

#ifdef GCAMORPH_ON_GPU


#include "gcamorphcpu.hpp"



namespace Freesurfer {


  void GCAmorphCPU::CheckIntegrity( void ) const {
    bool good;

    good = this->rx.MatchDims( this->ry );
    good = good && ( this->rx.MatchDims( this->rz ) );
    
    good = good && ( this->rx.MatchDims( this->origx ) );
    good = good && ( this->rx.MatchDims( this->origy ) );
    good = good && ( this->rx.MatchDims( this->origz ) );

    good = good && ( this->rx.MatchDims( this->dx ) );
    good = good && ( this->rx.MatchDims( this->dy ) );
    good = good && ( this->rx.MatchDims( this->dz ) );

    good = good && ( this->rx.MatchDims( this->odx ) );
    good = good && ( this->rx.MatchDims( this->ody ) );
    good = good && ( this->rx.MatchDims( this->odz ) );

    good = good && ( this->rx.MatchDims( this->origArea ) );
    good = good && ( this->rx.MatchDims( this->origArea1 ) );
    good = good && ( this->rx.MatchDims( this->origArea2 ) );

    good = good && ( this->rx.MatchDims( this->area ) );
    good = good && ( this->rx.MatchDims( this->area1 ) );
    good = good && ( this->rx.MatchDims( this->area2 ) );

    good = good && ( this->rx.MatchDims( this->invalid ) );
    good = good && ( this->rx.MatchDims( this->label ) );
    good = good && ( this->rx.MatchDims( this->status ) );
    good = good && ( this->rx.MatchDims( this->labelDist ) );

    good = good && ( this->rx.MatchDims( this->mean ) );
    good = good && ( this->rx.MatchDims( this->variance ) );

    if( !good ) {
      std::cerr << __FUNCTION__
		<< ": Dimension mismatch"
		<< std::endl;
      exit( EXIT_FAILURE );
    }
  }

  // -----------------------------

  void GCAmorphCPU::AllocateAll( const unsigned int nx,
				 const unsigned int ny,
				 const unsigned int nz ) {
    this->rx.Allocate( nx, ny, nz );
    this->ry.Allocate( nx, ny, nz );
    this->rz.Allocate( nx, ny, nz );

    this->origx.Allocate( nx, ny, nz );
    this->origy.Allocate( nx, ny, nz );
    this->origz.Allocate( nx, ny, nz );

    this->dx.Allocate( nx, ny, nz );
    this->dy.Allocate( nx, ny, nz );
    this->dz.Allocate( nx, ny, nz );

    this->odx.Allocate( nx, ny, nz );
    this->ody.Allocate( nx, ny, nz );
    this->odz.Allocate( nx, ny, nz );

    this->origArea.Allocate( nx, ny, nz );
    this->origArea1.Allocate( nx, ny, nz );
    this->origArea2.Allocate( nx, ny, nz );

    this->area.Allocate( nx, ny, nz );
    this->area1.Allocate( nx, ny, nz );
    this->area2.Allocate( nx, ny, nz );

    this->invalid.Allocate( nx, ny, nz );
    this->label.Allocate( nx, ny, nz );
    this->status.Allocate( nx, ny, nz );
    this->labelDist.Allocate( nx, ny, nz );

    this->mean.Allocate( nx, ny, nz );
    this->variance.Allocate( nx, ny, nz );
  }



  // --------------------------------

  void GCAmorphCPU::GetFromGPU( const GPU::Classes::GCAmorphGPU& src ) {
    /*!
      Retrieves linearly packed data from the GPU.
      Assumes that the GCAmorphGPU has already allocated its page-locked
      memory
    */

    // Get the data to page-locked memory
    src.d_rx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_rx );
    src.d_ry.RecvBuffer( GPU::Classes::GCAmorphGPU::h_ry );
    src.d_rz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_rz );

    src.d_origx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origx );
    src.d_origy.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origy );
    src.d_origz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origz );
    
    src.d_dx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dx );
    src.d_dy.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dy );
    src.d_dz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dz );
    
    src.d_odx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_odx );
    src.d_ody.RecvBuffer( GPU::Classes::GCAmorphGPU::h_ody );
    src.d_odz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_odz );
    
    src.d_origArea.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea );
    src.d_origArea1.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea1 );
    src.d_origArea2.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea2 );
    
    src.d_area.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area );
    src.d_area1.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area1 );
    src.d_area2.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area2 );

    src.d_invalid.RecvBuffer( GPU::Classes::GCAmorphGPU::h_invalid );
    src.d_status.RecvBuffer( GPU::Classes::GCAmorphGPU::h_status );
    src.d_label.RecvBuffer( GPU::Classes::GCAmorphGPU::h_label );
    src.d_labelDist.RecvBuffer( GPU::Classes::GCAmorphGPU::h_labelDist );
    
    src.d_mean.RecvBuffer( GPU::Classes::GCAmorphGPU::h_mean );
    src.d_variance.RecvBuffer( GPU::Classes::GCAmorphGPU:: h_variance );
    
    CUDA_SAFE_CALL( cudaThreadSynchronize() );
    
    // Transfer to internal buffers
    this->rx.SendBuffer( GPU::Classes::GCAmorphGPU::h_rx );
    this->ry.SendBuffer( GPU::Classes::GCAmorphGPU::h_ry );
    this->rz.SendBuffer( GPU::Classes::GCAmorphGPU::h_rz );

    this->origx.SendBuffer( GPU::Classes::GCAmorphGPU::h_origx );
    this->origy.SendBuffer( GPU::Classes::GCAmorphGPU::h_origy );
    this->origz.SendBuffer( GPU::Classes::GCAmorphGPU::h_origz );

    this->dx.SendBuffer( GPU::Classes::GCAmorphGPU::h_dx );
    this->dy.SendBuffer( GPU::Classes::GCAmorphGPU::h_dy );
    this->dz.SendBuffer( GPU::Classes::GCAmorphGPU::h_dz );
    
    this->odx.SendBuffer( GPU::Classes::GCAmorphGPU::h_odx );
    this->ody.SendBuffer( GPU::Classes::GCAmorphGPU::h_ody );
    this->odz.SendBuffer( GPU::Classes::GCAmorphGPU::h_odz );
    
    this->origArea.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea );
    this->origArea1.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea1 );
    this->origArea2.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea2 );
    
    this->area.SendBuffer( GPU::Classes::GCAmorphGPU::h_area );
    this->area1.SendBuffer( GPU::Classes::GCAmorphGPU::h_area1 );
    this->area2.SendBuffer( GPU::Classes::GCAmorphGPU::h_area2 );

    this->invalid.SendBuffer( GPU::Classes::GCAmorphGPU::h_invalid );
    this->status.SendBuffer( GPU::Classes::GCAmorphGPU::h_status );
    this->label.SendBuffer( GPU::Classes::GCAmorphGPU::h_label );
    this->labelDist.SendBuffer( GPU::Classes::GCAmorphGPU::h_labelDist );
    
    this->mean.SendBuffer( GPU::Classes::GCAmorphGPU::h_mean );
    this->variance.SendBuffer( GPU::Classes::GCAmorphGPU:: h_variance );

  }




  void GCAmorphCPU::PutOnGPU( GPU::Classes::GCAmorphGPU& dst ) const {
    /*!
      Puts the linearly packed data back on the GPU.
      Assumes that the GCAmorphGPU already has its page-locked
      memory allocated.
    */

    // Copy data to page-locked memory
    this->rx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_rx );
    this->ry.RecvBuffer( GPU::Classes::GCAmorphGPU::h_ry );
    this->rz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_rz );

    this->origx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origx );
    this->origy.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origy );
    this->origz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origz );

    this->dx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dx );
    this->dy.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dy );
    this->dz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_dz );
    
    this->odx.RecvBuffer( GPU::Classes::GCAmorphGPU::h_odx );
    this->ody.RecvBuffer( GPU::Classes::GCAmorphGPU::h_ody );
    this->odz.RecvBuffer( GPU::Classes::GCAmorphGPU::h_odz );
    
    this->origArea.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea );
    this->origArea1.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea1 );
    this->origArea2.RecvBuffer( GPU::Classes::GCAmorphGPU::h_origArea2 );
    
    this->area.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area );
    this->area1.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area1 );
    this->area2.RecvBuffer( GPU::Classes::GCAmorphGPU::h_area2 );

    this->invalid.RecvBuffer( GPU::Classes::GCAmorphGPU::h_invalid );
    this->status.RecvBuffer( GPU::Classes::GCAmorphGPU::h_status );
    this->label.RecvBuffer( GPU::Classes::GCAmorphGPU::h_label );
    this->labelDist.RecvBuffer( GPU::Classes::GCAmorphGPU::h_labelDist );
    
    this->mean.RecvBuffer( GPU::Classes::GCAmorphGPU::h_mean );
    this->variance.RecvBuffer( GPU::Classes::GCAmorphGPU:: h_variance );


    // Now, do the GPU transfers
    
    dst.d_rx.SendBuffer( GPU::Classes::GCAmorphGPU::h_rx );
    dst.d_ry.SendBuffer( GPU::Classes::GCAmorphGPU::h_ry );
    dst.d_rz.SendBuffer( GPU::Classes::GCAmorphGPU::h_rz );

    dst.d_origx.SendBuffer( GPU::Classes::GCAmorphGPU::h_origx );
    dst.d_origy.SendBuffer( GPU::Classes::GCAmorphGPU::h_origy );
    dst.d_origz.SendBuffer( GPU::Classes::GCAmorphGPU::h_origz );
    
    dst.d_dx.SendBuffer( GPU::Classes::GCAmorphGPU::h_dx );
    dst.d_dy.SendBuffer( GPU::Classes::GCAmorphGPU::h_dy );
    dst.d_dz.SendBuffer( GPU::Classes::GCAmorphGPU::h_dz );
    
    dst.d_odx.SendBuffer( GPU::Classes::GCAmorphGPU::h_odx );
    dst.d_ody.SendBuffer( GPU::Classes::GCAmorphGPU::h_ody );
    dst.d_odz.SendBuffer( GPU::Classes::GCAmorphGPU::h_odz );
    
    dst.d_origArea.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea );
    dst.d_origArea1.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea1 );
    dst.d_origArea2.SendBuffer( GPU::Classes::GCAmorphGPU::h_origArea2 );
    
    dst.d_area.SendBuffer( GPU::Classes::GCAmorphGPU::h_area );
    dst.d_area1.SendBuffer( GPU::Classes::GCAmorphGPU::h_area1 );
    dst.d_area2.SendBuffer( GPU::Classes::GCAmorphGPU::h_area2 );

    dst.d_invalid.SendBuffer( GPU::Classes::GCAmorphGPU::h_invalid );
    dst.d_status.SendBuffer( GPU::Classes::GCAmorphGPU::h_status );
    dst.d_label.SendBuffer( GPU::Classes::GCAmorphGPU::h_label );
    dst.d_labelDist.SendBuffer( GPU::Classes::GCAmorphGPU::h_labelDist );
    
    dst.d_mean.SendBuffer( GPU::Classes::GCAmorphGPU::h_mean );
    dst.d_variance.SendBuffer( GPU::Classes::GCAmorphGPU:: h_variance );

    CUDA_SAFE_CALL( cudaThreadSynchronize() );
  }


}





#endif
