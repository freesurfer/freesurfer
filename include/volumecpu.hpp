/**
 * @file  volumecpu.hpp
 * @brief Holds templated datatype for volume data on the CPU
 *
 * Holds a templated datatype for volume data on the CPU.
 * It is suspiciously similar to VolumeGPU - indeed it exists
 * because I can't come up with a simple method of implementing
 * a portion of mri_ca_register on the GPU (gcamLabelTerm, to be exact)
 * and I don't want a full repacking penalty.
 * You may well complain that this is horribly inelegant....
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/09/28 16:44:28 $
 *    $Revision: 1.1 $
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


#include <cstdlib>

#include <iostream>

namespace Freesurfer {

  //! Templated "mutator" class for the CPU
  template<typename T>
  class VolumeArgCPU {
  public:

    //! Volume sizes
    const unsigned int nx, ny, nz;

    // ------------------------------------------

    //! Default constructor
    VolumeArgCPU( void ) : nx(0), ny(0), nz(0),
			   data(NULL) {};

    //! Constructor from inputs
    VolumeArgCPU( const unsigned int _nx,
		  const unsigned int _ny,
		  const unsigned int _nz,
		  T* const _data ) : nx(_nx), ny(_ny), nz(_nz),
				     data(_data) {};

    // ---------------------------------------------------
    // Subscripting operators (unsafe)
    
    //! RHS subscripting operator
    T operator()( const unsigned int ix,
		  const unsigned int iy,
		  const unsigned int iz ) const {
      return( data[ix + ( nx * ( iy + (ny*iz) ) )] );
    }

    //! LHS subscripting operator
    T& operator()( const unsigned int ix,
		   const unsigned int iy,
		   const unsigned int iz ) {
      return( data[ix + ( nx * ( iy + (ny*iz) ) )] );
    }

  private:
    //! Data pointer
    T* const data;

  };



  //! Templated management class for the CPU
  template<typename T>
  class VolumeCPU {
  public:

    //! Default constructor
    VolumeCPU( void ) : nx(0), ny(0), nz(0),
			data(NULL) {};

    //! Destructor
    ~VolumeCPU( void ) {
      this->Release();
    }

    //! Conversion operator
    operator VolumeArgCPU<T>( void ) const {
      VolumeArgCPU<T> vac( this->nx, this->ny, this->nz,
			   this->data );

      return( vac );
    }

    // --------------------------------

    //! Allocates storage of given dimensions
    void Allocate( const unsigned int xn,
		   const unsigned int yn,
		   const unsigned int zn ) {
      
      // Check for prior allocation
      if( this->data != NULL ) {
	if( (xn == this->nx) &&
	    (yn == this->ny) &&
	    (zn == this->nz) ) {
	  // Same size
	  return;
	} else {
	  this->Release();
	}
      }

      this->nx = xn;
      this->ny = yn;
      this->nz = zn;

      this->data = new T[xn*yn*zn];
    }

    //! Releases memory
    void Relase( void ) {
      if( this->data != NULL ) {
	this->nx = 0;
	this->ny = 0;
	this->nz = 0;
	delete[] this->data;
	this->data = NULL;
      }
    }



  private:
    unsigned int nx, ny, nz;
    T* data;
  };


}
