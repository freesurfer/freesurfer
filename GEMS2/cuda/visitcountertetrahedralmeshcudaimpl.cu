#include "visitcountertetrahedralmeshcudaimpl.hpp"

#include "visitcounteraction.hpp"
#include "simplesharedtetrahedroninterior.hpp"


namespace kvl {
  namespace cuda {
    void RunVisitCounterTetrahedralMeshCUDA( CudaImage<int,3,unsigned short>& d_output,
					     const CudaTetrahedralMesh<double,unsigned long,float>& ctm ) {
      
      VisitCounterAction vca(d_output.getArg());
      auto domain = d_output.GetDimensions();
      
      // Run the kernel
      
      RunSimpleSharedTetrahedron( domain[0], domain[1], domain[2],
				  ctm, vca );
    } 
  }
}
