#pragma once

#include "cudaimage.hpp"
#include "cudatetrahedralmesh.hpp"

namespace kvl {
  namespace cuda {
    void RunAtlasMeshAlphaDrawerCUDA( CudaImage<float,3,unsigned short>& d_output,
				      const CudaTetrahedralMesh<double,unsigned long,float>& ctm,
				      const unsigned int iAlpha );
  }
}
