#include "atlasmeshalphadrawercudaimpl.hpp"

#include "alphadraweraction.hpp"
#include "simplesharedtetrahedroninterior.hpp"

namespace kvl {
  namespace cuda {
    void RunAtlasMeshAlphaDrawerCUDA( CudaImage<float,3,unsigned short>& d_output,
				      const CudaTetrahedralMesh<double,unsigned long,float>& ctm,
				      const unsigned int iAlpha ) {
      AlphaDrawerAction<float> ada(d_output.getArg(), iAlpha );
      auto domain = d_output.GetDimensions();

      RunSimpleSharedTetrahedron( domain[0], domain[1], domain[2], ctm, ada );
    }
  }
}
