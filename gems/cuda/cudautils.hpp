#pragma once

namespace kvl {
  namespace cuda {
    static unsigned int GetBlockSize( const size_t nVoxels, const size_t nTetrahedra ) {
      unsigned int blockSize = 8;
      
      // Compute the average number of voxels in each tetrahedron
      const float avgTetVol = nVoxels / nTetrahedra;
      // Compute the equivalent average cube (factor of 6 to get volume of bounding box)
      const float avgCubeSize = powf( 6 * avgTetVol, 0.3333333f );
      
      // Check to see if we need to increase blockSize
      if( avgCubeSize > 32 ) {
	blockSize = 16;
	if( avgCubeSize > 64 ) {
	  blockSize = 32;
	}
      }
      return blockSize;
    }
  }
}
