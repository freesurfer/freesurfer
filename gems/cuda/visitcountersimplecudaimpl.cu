#include <stdexcept>

#include "visitcountersimplecudaimpl.hpp"

#include "visitcounteraction.hpp"
#include "simplesharedtetrahedroninterior.hpp"

const unsigned int nDims = 3;
const unsigned int nVertices = 4;

// ---------

template<typename ArgType>
class SimpleMesh_GPU {
public:
  typedef size_t MeshIdxType; 

  __device__ __host__
  SimpleMesh_GPU( const kvl::cuda::Image_GPU<ArgType,3,size_t>& tetrahedra ) : tetInfo(tetrahedra) {}

  __device__
  ArgType GetVertexCoordinate( size_t iTet, size_t iVert, size_t iDim ) const {
    return this->tetInfo(iTet,iVert,threadIdx.x);
  }

  __device__
  size_t GetTetrahedraCount() const {
    return this->tetInfo.dims[0];
  }

private:
  const kvl::cuda::Image_GPU<ArgType,3,size_t> tetInfo;
};

namespace kvl {
  namespace cuda {
    
    template<typename CoordinateType>
    class SimpleMeshSupply {
    public:
      typedef SimpleMesh_GPU<CoordinateType> GPUType;
      typedef CoordinateType CoordType;
      
      SimpleMeshSupply( const kvl::cuda::CudaImage<CoordinateType,3,size_t>& tetrahedra ) : d_tetInfo(tetrahedra) {}
      
      SimpleMesh_GPU<CoordinateType> getArg() const {
	return SimpleMesh_GPU<CoordinateType>(this->d_tetInfo.getArg());
      }
      
      size_t GetTetrahedraCount() const {
	return this->d_tetInfo.GetDimensions()[0];
       }
      
    private:
      const kvl::cuda::CudaImage<CoordinateType,3,size_t>& d_tetInfo;
    };

    template<>
    void RunVisitCounterSimpleCUDA<float,float>( CudaImage<int,3,unsigned short>& d_output,
						 const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      VisitCounterAction vca(d_output.getArg());
      auto domain = d_output.GetDimensions();
      auto mesh = SimpleMeshSupply<float>(d_tetrahedra);

      RunSimpleSharedTetrahedron<SimpleMeshSupply<float>,VisitCounterAction,unsigned short,float>(domain[0], domain[1], domain[2], mesh, vca );
    }

    template<>
    void RunVisitCounterSimpleCUDA<double,double>( CudaImage<int,3,unsigned short>& d_output,
						   const CudaImage<double,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      VisitCounterAction vca(d_output.getArg());
      auto domain = d_output.GetDimensions();
      auto mesh = SimpleMeshSupply<double>(d_tetrahedra);

      RunSimpleSharedTetrahedron<SimpleMeshSupply<double>,VisitCounterAction,unsigned short, double>(domain[0], domain[1], domain[2], mesh, vca );
    }

    template<>
    void RunVisitCounterSimpleCUDA<float,double>( CudaImage<int,3,unsigned short>& d_output,
						  const CudaImage<float,3,size_t>& d_tetrahedra ) {
      if( d_tetrahedra.GetDimensions()[1] != nVertices ) {
	throw std::runtime_error("Must have four vertices per tetrahedron!");
      }
      if( d_tetrahedra.GetDimensions()[2] != nDims ) {
	throw std::runtime_error("Only implemented for 3D space");
      }

      VisitCounterAction vca(d_output.getArg());
      auto domain = d_output.GetDimensions();
      auto mesh = SimpleMeshSupply<float>(d_tetrahedra);

      RunSimpleSharedTetrahedron<SimpleMeshSupply<float>,VisitCounterAction,unsigned short, double>(domain[0], domain[1], domain[2], mesh, vca );
    }
  }
}
